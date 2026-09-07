/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2026 Danny van Dyk
 * Copyright (c) 2011      Frederik Beaujean
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/utils/expression-cacher.hh>
#include <eos/utils/expression-observable.hh>
#include <eos/utils/lock.hh>
#include <eos/utils/log.hh>
#include <eos/utils/mutex.hh>
#include <eos/utils/observable_cache.hh>
#include <eos/utils/observable_set.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/thread_pool.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <format>
#include <iterator>
#include <limits>
#include <map>
#include <tuple>
#include <typeindex>
#include <unordered_set>
#include <utility>
#include <vector>

namespace eos
{
    template <> struct WrappedForwardIteratorTraits<ObservableCache::IteratorTag>
    {
            using UnderlyingIterator = std::vector<ObservablePtr>::iterator;
    };
    template class WrappedForwardIterator<ObservableCache::IteratorTag, ObservablePtr>;

    template <> struct Implementation<ObservableCache>
    {
            // Parameters which are common to all observables in the cache.
            Parameters parameters;

            // Contains each observable that needs to be calculated exactly once
            std::vector<ObservablePtr> observables;

            // Contains each regular observable and its associated index
            std::vector<std::tuple<ObservablePtr, ObservableCache::ObservableId>> regular_observables;

            // Contains each cacheable observable and its associated index
            std::multimap<std::type_index, std::tuple<CacheableObservable *, ObservableCache::ObservableId>> cacheable_observables;

            // Contains each cached observable and its associated index
            std::vector<std::tuple<ObservablePtr, ObservableCache::ObservableId>> cached_observables;

            // Contains each expression observable and its associated index
            std::vector<std::tuple<ObservablePtr, ObservableCache::ObservableId>> expression_observables;

            // Contains values of all observables
            std::vector<double> predictions;

            // The observables and predictions of all batches, stored back to back so that update()
            // can partition them without building a work list.
            std::vector<ObservablePtr> batch_observables;
            std::vector<double>        batch_predictions;

            // A batch is a contiguous block of the above, addressed by its BatchId (i.e. its index).
            struct Batch
            {
                    std::size_t offset;
                    std::size_t size;
            };

            std::vector<Batch> batches;

            // Identifies this cache; the counter advances once per successful update(), never on a failed one.
            uint64_t instance, counter;

            Implementation(const Parameters & parameters) :
                parameters(parameters),
                instance(_next_instance()),
                counter(0u)
            {
            }

            static uint64_t
            _next_instance()
            {
                static std::atomic<uint64_t> next(1u);

                return next.fetch_add(1u, std::memory_order_relaxed);
            }

            ~Implementation() {}

            static bool
            identical_observables(const ObservablePtr & lhs, const ObservablePtr & rhs)
            {
                const KinematicUser &                     kinematic_user_lhs = static_cast<const KinematicUser &>(*lhs);
                const KinematicUser &                     kinematic_user_rhs = static_cast<const KinematicUser &>(*rhs);
                std::unordered_set<KinematicVariable::Id> kinematic_ids_lhs(kinematic_user_lhs.begin_kinematics(), kinematic_user_lhs.end_kinematics());
                std::unordered_set<KinematicVariable::Id> kinematic_ids_rhs(kinematic_user_rhs.begin_kinematics(), kinematic_user_rhs.end_kinematics());

                // compare used_kinematics
                if (kinematic_ids_lhs != kinematic_ids_rhs)
                {
                    return false;
                }

                // compare name
                if (lhs->name() != rhs->name())
                {
                    return false;
                }

                // compare kinematics
                if (lhs->kinematics() != rhs->kinematics())
                {
                    return false;
                }

                // compare options
                if (lhs->options() != rhs->options())
                {
                    return false;
                }

                return true;
            }

            ObservableCache::ObservableId
            add(const ObservablePtr & observable, const ObservableCache & cache)
            {
                if (observable->parameters() != parameters)
                {
                    throw InternalError("ObservableCache::add(): Mismatch of Parameters between different observables detected.");
                }

                // compare each observable for options, kinematics and name
                unsigned index = 0;
                for (auto i = observables.begin(), i_end = observables.end(); i != i_end; ++i, ++index)
                {
                    if (identical_observables(*i, observable))
                    {
                        return ObservableCache::ObservableId(index);
                    }
                }

                CacheableObservable *  cacheable_observable  = dynamic_cast<CacheableObservable *>(observable.get());
                ExpressionObservable * expression_observable = dynamic_cast<ExpressionObservable *>(observable.get());

                if (nullptr != expression_observable) // is the new observable an expression?
                {
                    ObservablePtr cached_expression_observable(new ExpressionObservable(expression_observable->name(),
                                                                                        cache,
                                                                                        expression_observable->kinematics(),
                                                                                        expression_observable->options(),
                                                                                        expression_observable->expression()));

                    // ensure that the new index is correct, since the ExpressionCacher is capable to modify our cache
                    index = observables.size();

                    observables.push_back(cached_expression_observable);
                    predictions.push_back(std::numeric_limits<double>::quiet_NaN());
                    expression_observables.push_back(std::make_tuple(cached_expression_observable, ObservableCache::ObservableId(index)));

                    return ObservableCache::ObservableId(index);
                }
                else if (nullptr != cacheable_observable) // is the new observable cacheable?
                {
                    std::type_index type_index(cacheable_observable->prepare_type_index());

                    // have we encountered a cacheable observable with a compatible intermediate result before?
                    auto range = cacheable_observables.equal_range(type_index);
                    for (auto c = range.first, c_end = range.second; c != c_end; ++c)
                    {
                        // attempt to adopt its intermediate result...
                        ObservablePtr cached_observable = cacheable_observable->make_cached_observable(std::get<0>(c->second));
                        if (! cached_observable)
                        {
                            continue;
                        }

                        // add the newly created cached observable
                        observables.push_back(cached_observable);
                        predictions.push_back(std::numeric_limits<double>::quiet_NaN());
                        cached_observables.push_back(std::make_tuple(cached_observable, ObservableCache::ObservableId(index)));

                        return ObservableCache::ObservableId(index);
                    }

                    // else add this new cacheable observable
                    observables.push_back(observable);
                    predictions.push_back(std::numeric_limits<double>::quiet_NaN());
                    cacheable_observables.insert(std::make_pair(type_index, std::make_tuple(cacheable_observable, ObservableCache::ObservableId(index))));

                    return ObservableCache::ObservableId(index);
                }
                else
                {
                    // add this new regular observable
                    observables.push_back(observable);
                    predictions.push_back(std::numeric_limits<double>::quiet_NaN());
                    regular_observables.push_back(std::make_tuple(observable, ObservableCache::ObservableId(index)));

                    return ObservableCache::ObservableId(index);
                }

                throw InternalError("should not be reached");
            }

            ObservableCache::BatchId
            add_batch(std::vector<ObservablePtr> && computations)
            {
                for (const auto & observable : computations)
                {
                    if (observable->parameters() != parameters)
                    {
                        throw InternalError("ObservableCache::add_batch(): Mismatch of Parameters between different observables detected.");
                    }
                }

                const unsigned index = batches.size();

                const Batch batch{ batch_observables.size(), computations.size() };
                batch_observables.insert(batch_observables.end(), std::make_move_iterator(computations.begin()), std::make_move_iterator(computations.end()));
                batch_predictions.resize(batch_observables.size(), std::numeric_limits<double>::quiet_NaN());
                batches.push_back(batch);

                return ObservableCache::BatchId(index);
            }
    };

    ObservableCache::ObservableCache(const Parameters & parameters) :
        PrivateImplementationPattern<ObservableCache>(new Implementation<ObservableCache>(parameters))
    {
    }

    ObservableCache::~ObservableCache() {}

    bool
    ObservableCache::operator== (const ObservableCache & rhs) const
    {
        return _imp.get() == rhs._imp.get();
    }

    ObservableCache::Generation
    ObservableCache::generation() const
    {
        return Generation(_imp->instance, _imp->counter);
    }

    ObservableCache::ObservableId
    ObservableCache::add(const ObservablePtr & observable)
    {
        return _imp->add(observable, *this);
    }

    ObservableCache::BatchId
    ObservableCache::add_batch(std::vector<ObservablePtr> && computations)
    {
        return _imp->add_batch(std::move(computations));
    }

    void
    ObservableCache::update()
    {
        // an observable may call back into the embedding runtime from one of the pool's threads
        const auto guard = ThreadPool::instance()->wait_guard();

        // The exception raised first by any evaluation path below, kept across threads so that
        // update() can rethrow it once every ticket has been waited on.
        Mutex              exception_mutex;
        std::exception_ptr exception;

        // Catches whatever evaluate() throws -- not only eos::Exception, since anything else would
        // otherwise escape into ThreadPool::thread_function and terminate the process. Logs the
        // observable and the error, then stashes the first exception seen for update() to rethrow.
        auto record_exception = [&exception_mutex, &exception](const char * kind, Observable & o)
        {
            const std::exception_ptr current = std::current_exception();

            std::string what = "unknown exception";
            try
            {
                std::rethrow_exception(current);
            }
            catch (const std::exception & e)
            {
                what = e.what();
            }
            catch (...)
            {
            }

            Log::instance()->message("ObservableCache::update", ll_error) << "Exception encountered when evaluating " << kind << " observable '" << o.name() << "["
                                                                          << o.kinematics().as_string() << "];" << o.options().as_string() << "': " << what;

            Lock lock(exception_mutex);
            if (! exception)
            {
                exception = current;
            }
        };

        // parallelize the evaluation of the observables
        std::vector<Ticket> cacheable_tickets;
        cacheable_tickets.reserve(_imp->cacheable_observables.size());

        // evaluate all cacheable observables in parallel
        for (auto co : _imp->cacheable_observables)
        {
            auto f = [=, this, &record_exception]()
            {
                auto & o  = std::get<0>(co.second);
                auto & id = std::get<1>(co.second);
                try
                {
                    _imp->predictions[id.value()] = o->evaluate();
                }
                catch (...)
                {
                    record_exception("cacheable", *o);
                    _imp->predictions[id.value()] = std::numeric_limits<double>::quiet_NaN();
                }
            };
            cacheable_tickets.push_back(ThreadPool::instance()->enqueue(std::function<void(void)>(f)));
        }

        std::vector<Ticket> regular_tickets;
        regular_tickets.reserve(_imp->regular_observables.size());

        // evaluate all regular observables in parallel
        for (auto ro : _imp->regular_observables)
        {
            auto f = [=, this, &record_exception]()
            {
                auto & o  = std::get<0>(ro);
                auto & id = std::get<1>(ro);
                try
                {
                    _imp->predictions[id.value()] = o->evaluate();
                }
                catch (...)
                {
                    record_exception("regular", *o);
                    _imp->predictions[id.value()] = std::numeric_limits<double>::quiet_NaN();
                }
            };
            regular_tickets.push_back(ThreadPool::instance()->enqueue(std::function<void(void)>(f)));
        }

        // Evaluate all batched computations in roughly one contiguous slice per thread; distinct
        // slices write to disjoint prediction slots, so no locking is required.
        struct BatchSlicing
        {
                Implementation<ObservableCache> * imp;
                decltype(record_exception) *      record;
                std::size_t                       total;
                std::size_t                       slice;
        };

        const std::size_t   number_of_slices = std::max<std::size_t>(1u, ThreadPool::instance()->number_of_threads());
        const std::size_t   total            = _imp->batch_observables.size();
        const BatchSlicing  slicing{ _imp.get(), &record_exception, total, (total + number_of_slices - 1) / number_of_slices };
        std::vector<Ticket> batch_tickets;
        batch_tickets.reserve(std::min(total, number_of_slices));

        for (std::size_t start = 0; start < total; start += slicing.slice)
        {
            // two words of captures fit into std::function's inline storage, avoiding a heap allocation
            auto f = [s = &slicing, start]()
            {
                const std::size_t end = std::min(s->total, start + s->slice);
                for (std::size_t k = start; k < end; ++k)
                {
                    const auto & o = s->imp->batch_observables[k];
                    try
                    {
                        s->imp->batch_predictions[k] = o->evaluate();
                    }
                    catch (...)
                    {
                        (*s->record)("batched", *o);
                        s->imp->batch_predictions[k] = std::numeric_limits<double>::quiet_NaN();
                    }
                }
            };
            batch_tickets.push_back(ThreadPool::instance()->enqueue(std::function<void(void)>(f)));
        }

        // await completion of the cacheable observables
        for (auto ticket : cacheable_tickets)
        {
            ticket.wait();
        }

        std::vector<Ticket> cached_tickets;
        cached_tickets.reserve(_imp->cached_observables.size());

        // evaluate all cached observables in parallel
        for (auto co : _imp->cached_observables)
        {
            auto f = [=, this, &record_exception]()
            {
                auto & o  = std::get<0>(co);
                auto & id = std::get<1>(co);
                try
                {
                    _imp->predictions[id.value()] = o->evaluate();
                }
                catch (...)
                {
                    record_exception("cached", *o);
                    _imp->predictions[id.value()] = std::numeric_limits<double>::quiet_NaN();
                }
            };
            cached_tickets.push_back(ThreadPool::instance()->enqueue(std::function<void(void)>(f)));
        }

        // await completion of the regular observables
        for (auto ticket : regular_tickets)
        {
            ticket.wait();
        }

        // await completion of the cached observables
        for (auto ticket : cached_tickets)
        {
            ticket.wait();
        }

        // await completion of the batched computations
        for (auto ticket : batch_tickets)
        {
            ticket.wait();
        }

        // evaluate all expression observables in a serial fashion
        //
        // This is necessary, since an expression observable can rely on
        // another expression observable, which would be located earlier in
        // the sequence.
        // Serial evaluation ensures that no race conditions arise.
        // There is not reason to optimize this, since expression observables
        // are evaluated very quickly.
        for (auto eo : _imp->expression_observables)
        {
            auto & o  = std::get<0>(eo);
            auto & id = std::get<1>(eo);
            try
            {
                _imp->predictions[id.value()] = o->evaluate();
            }
            catch (...)
            {
                record_exception("expression", *o);
                _imp->predictions[id.value()] = std::numeric_limits<double>::quiet_NaN();
            }
        }

        // every ticket has been waited on; surface the first failure on the calling thread rather
        // than letting a NaN prediction propagate into a log-likelihood
        if (exception)
        {
            std::rethrow_exception(exception);
        }

        ++_imp->counter;
    }

    Parameters
    ObservableCache::parameters() const
    {
        return _imp->parameters;
    }

    double
    ObservableCache::operator[] (const ObservableCache::ObservableId & id) const
    {
        return _imp->predictions[id.value()];
    }

    std::span<const double>
    ObservableCache::operator[] (const ObservableCache::BatchId & id) const
    {
        if (id.value() >= _imp->batches.size())
        {
            throw InternalError(std::format("ObservableCache: BatchId {} is not valid for this cache", id.value()));
        }

        const auto & batch = _imp->batches[id.value()];

        return std::span<const double>(_imp->batch_predictions.data() + batch.offset, batch.size);
    }

    ObservablePtr
    ObservableCache::observable(const ObservableCache::ObservableId & id) const
    {
        return _imp->observables[id.value()];
    }

    unsigned
    ObservableCache::size() const
    {
        return _imp->observables.size();
    }

    ObservableCache::Iterator
    ObservableCache::begin() const
    {
        return _imp->observables.begin();
    }

    ObservableCache::Iterator
    ObservableCache::end() const
    {
        return _imp->observables.end();
    }

    ObservableCache
    ObservableCache::clone(const Parameters & parameters) const
    {
        ObservableCache result(parameters);

        // batches are not cloned; their owners (e.g. LogLikelihoodBlocks) register them again
        for (auto o = _imp->observables.begin(), o_end = _imp->observables.end(); o != o_end; ++o)
        {
            // cloning cached observables creates independent *cacheable* observables
            // adding them back creates new and independent cached observables
            result._imp->add((*o)->clone(parameters), *this);
        }

        result.update();

        return result;
    }
} // namespace eos
