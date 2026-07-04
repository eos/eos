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
#include <eos/utils/log.hh>
#include <eos/utils/observable_cache.hh>
#include <eos/utils/observable_set.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/thread_pool.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <algorithm>
#include <cstddef>
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

            // A batch bundles a group of independent computations whose predictions are stored in a
            // single contiguous block of memory. Batches are addressed by their own BatchId and are
            // kept separate from the single-observable id space above.
            struct Batch
            {
                    std::vector<ObservablePtr> observables;
                    std::vector<double>        predictions;
            };

            // Contains each batch of computations, addressed by its BatchId (i.e. its index).
            std::vector<Batch> batches;

            Implementation(const Parameters & parameters) :
                parameters(parameters)
            {
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
                    std::type_index type_index(typeid(*cacheable_observable));

                    // have we encountered this type of cacheable observable before?
                    auto range = cacheable_observables.equal_range(type_index);
                    for (auto c = range.first, c_end = range.second; c != c_end; ++c)
                    {
                        // have we encountered this cacheable observable with the same properties before?
                        if (std::get<0>(c->second)->kinematics() != cacheable_observable->kinematics())
                        {
                            continue;
                        }

                        if (std::get<0>(c->second)->options() != cacheable_observable->options())
                        {
                            continue;
                        }
                        const KinematicUser &                     kinematic_user_cacheable_obs = static_cast<const KinematicUser &>(*cacheable_observable);
                        const KinematicUser &                     kinematic_user_cached_obs    = static_cast<const KinematicUser &>(*std::get<0>(c->second));
                        std::unordered_set<KinematicVariable::Id> kinematic_ids_cacheable_obs(kinematic_user_cacheable_obs.begin_kinematics(),
                                                                                              kinematic_user_cacheable_obs.end_kinematics());
                        std::unordered_set<KinematicVariable::Id> kinematic_ids_cached_obs(kinematic_user_cached_obs.begin_kinematics(),
                                                                                           kinematic_user_cached_obs.end_kinematics());
                        if (kinematic_ids_cacheable_obs != kinematic_ids_cached_obs)
                        {
                            continue;
                        }

                        // yes! cache it...
                        ObservablePtr cached_observable = cacheable_observable->make_cached_observable(std::get<0>(c->second));
                        if (! cached_observable)
                        {
                            throw InternalError("make_cached_observable() failed");
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

                Batch batch;
                batch.predictions.assign(computations.size(), std::numeric_limits<double>::quiet_NaN());
                batch.observables = std::move(computations);
                batches.push_back(std::move(batch));

                return ObservableCache::BatchId(index);
            }
    };

    ObservableCache::ObservableCache(const Parameters & parameters) :
        PrivateImplementationPattern<ObservableCache>(new Implementation<ObservableCache>(parameters))
    {
    }

    ObservableCache::~ObservableCache() {}

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
        // parallelize the evaluation of the observables
        std::vector<Ticket> cacheable_tickets;
        cacheable_tickets.reserve(_imp->cacheable_observables.size());

        // evaluate all cacheable observables in parallel
        for (auto co : _imp->cacheable_observables)
        {
            auto f = [=, this]()
            {
                auto & o  = std::get<0>(co.second);
                auto & id = std::get<1>(co.second);
                try
                {
                    _imp->predictions[id.value()] = o->evaluate();
                }
                catch (eos::Exception & e)
                {
                    Log::instance()->message("ObservableCache::update", ll_error) << "Exception encountered when evaluating cacheable observable '" << o->name() << "["
                                                                                  << o->kinematics().as_string() << "];" << o->options().as_string() << "': " << e.what();
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
            auto f = [=, this]()
            {
                auto & o  = std::get<0>(ro);
                auto & id = std::get<1>(ro);
                try
                {
                    _imp->predictions[id.value()] = o->evaluate();
                }
                catch (eos::Exception & e)
                {
                    Log::instance()->message("ObservableCache::update", ll_error) << "Exception encountered when evaluating regular observable '" << o->name() << "["
                                                                                  << o->kinematics().as_string() << "];" << o->options().as_string() << "': " << e.what();
                    _imp->predictions[id.value()] = std::numeric_limits<double>::quiet_NaN();
                }
            };
            regular_tickets.push_back(ThreadPool::instance()->enqueue(std::function<void(void)>(f)));
        }

        // Flatten the computations of all batches into a single work list of (batch index, index
        // within batch) pairs. This lets us partition the total amount of batched work across the
        // available threads irrespective of how it is distributed over the individual batches.
        std::vector<std::pair<std::size_t, std::size_t>> batch_work;
        for (std::size_t bi = 0; bi < _imp->batches.size(); ++bi)
        {
            for (std::size_t oi = 0; oi < _imp->batches[bi].observables.size(); ++oi)
            {
                batch_work.emplace_back(bi, oi);
            }
        }

        std::vector<Ticket> batch_tickets;

        // Evaluate all batched computations in parallel by handing each thread a contiguous slice of
        // the work list; the computations within a slice are evaluated in series. Aiming for roughly
        // one slice per thread keeps every thread busy while amortising the per-ticket overhead over
        // many computations. Distinct slices write to disjoint prediction slots, so no locking is
        // required.
        if (! batch_work.empty())
        {
            const unsigned    number_of_threads = ThreadPool::instance()->number_of_threads();
            const std::size_t target_slices     = std::max<std::size_t>(1u, number_of_threads);
            const std::size_t slice             = (batch_work.size() + target_slices - 1) / target_slices;

            for (std::size_t start = 0; start < batch_work.size(); start += slice)
            {
                const std::size_t end = std::min(batch_work.size(), start + slice);

                // 'batch_work' outlives every ticket, since update() waits for all batch tickets
                // before returning; capturing it by reference therefore avoids copying the slice.
                auto f = [start, end, &batch_work, this]()
                {
                    for (std::size_t k = start; k < end; ++k)
                    {
                        auto &       batch = _imp->batches[batch_work[k].first];
                        const auto & o     = batch.observables[batch_work[k].second];
                        try
                        {
                            batch.predictions[batch_work[k].second] = o->evaluate();
                        }
                        catch (eos::Exception & e)
                        {
                            Log::instance()->message("ObservableCache::update", ll_error) << "Exception encountered when evaluating batched observable '" << o->name() << "["
                                                                                          << o->kinematics().as_string() << "];" << o->options().as_string() << "': " << e.what();
                            batch.predictions[batch_work[k].second] = std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                };
                batch_tickets.push_back(ThreadPool::instance()->enqueue(std::function<void(void)>(f)));
            }
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
            auto f = [=, this]()
            {
                auto & o  = std::get<0>(co);
                auto & id = std::get<1>(co);
                try
                {
                    _imp->predictions[id.value()] = o->evaluate();
                }
                catch (eos::Exception & e)
                {
                    Log::instance()->message("ObservableCache::update", ll_error) << "Exception encountered when evaluating cached observable '" << o->name() << "["
                                                                                  << o->kinematics().as_string() << "];" << o->options().as_string() << "': " << e.what();
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
            catch (eos::Exception & e)
            {
                Log::instance()->message("ObservableCache::update", ll_error) << "Exception encountered when evaluating expression observable '" << o->name() << "["
                                                                              << o->kinematics().as_string() << "];" << o->options().as_string() << "': " << e.what();
                _imp->predictions[id.value()] = std::numeric_limits<double>::quiet_NaN();
            }
        }
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
        const auto & batch = _imp->batches[id.value()];

        return std::span<const double>(batch.predictions.data(), batch.predictions.size());
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

        // Note: batches are intentionally not cloned here. Their owners (e.g. LogLikelihoodBlocks)
        // re-register them via add_batch() on the cloned cache, which reproduces the original BatchId
        // ordering. This mirrors how single observables are re-added by their owners after cloning.
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
