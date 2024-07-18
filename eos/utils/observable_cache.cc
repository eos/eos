/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2024 Danny van Dyk
 * Copyright (c) 2011 Frederik Beaujean
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
#include <limits>
#include <map>
#include <tuple>
#include <typeindex>
#include <vector>

namespace eos
{
    template <>
    struct WrappedForwardIteratorTraits<ObservableCache::IteratorTag>
    {
        using UnderlyingIterator = std::vector<ObservablePtr>::iterator;
    };
    template class WrappedForwardIterator<ObservableCache::IteratorTag, ObservablePtr>;

    template <> struct
    Implementation<ObservableCache>
    {
        // Parameters which are common to all observables in the cache.
        Parameters parameters;

        // Contains each observable that needs to be calculated exactly once
        std::vector<ObservablePtr> observables;

        // Contains each regular observable and its associated index
        std::vector<std::tuple<ObservablePtr, ObservableCache::Id>> regular_observables;

        // Contains each cacheable observable and its associated index
        std::multimap<std::type_index, std::tuple<CacheableObservable *, ObservableCache::Id>> cacheable_observables;

        // Contains each cached observable and its associated index
        std::vector<std::tuple<ObservablePtr, ObservableCache::Id>> cached_observables;

        // Contains each expression observable and its associated index
        std::vector<std::tuple<ObservablePtr, ObservableCache::Id>> expression_observables;

        // Contains values of all observables
        std::vector<double> predictions;

        Implementation(const Parameters & parameters) :
            parameters(parameters)
        {
        }

        ~Implementation()
        {
        }

        static bool identical_observables(const ObservablePtr & lhs, const ObservablePtr & rhs)
        {
            // compare name
            if (lhs->name() != rhs->name())
            {
                return false;
            }

            // compare kinematics
            if( lhs->kinematics() != rhs->kinematics())
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

        ObservableCache::Id add(const ObservablePtr & observable, const ObservableCache & cache)
        {
            if (observable->parameters() != parameters)
                throw InternalError("ObservableCache::add(): Mismatch of Parameters between different observables detected.");

            // compare each observable for options, kinematics and name
            unsigned index = 0;
            for (auto i = observables.begin(), i_end = observables.end() ; i != i_end ; ++i, ++index)
            {
                if (identical_observables(*i, observable))
                    return index;
            }

            CacheableObservable * cacheable_observable = dynamic_cast<CacheableObservable *>(observable.get());
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
                expression_observables.push_back(std::make_tuple(cached_expression_observable, index));

                return index;
            }
            else if (nullptr != cacheable_observable) // is the new observable cacheable?
            {
                std::type_index type_index(typeid(*cacheable_observable));

                // have we encountered this type of cacheable observable before?
                auto range = cacheable_observables.equal_range(type_index);
                for (auto c = range.first, c_end = range.second ; c != c_end ; ++c)
                {
                    // have we encountered this cacheable observable with the same properties before?
                    if (std::get<0>(c->second)->kinematics() != cacheable_observable->kinematics())
                        continue;

                    if (std::get<0>(c->second)->options() != cacheable_observable->options())
                        continue;

                    // yes! cache it...
                    ObservablePtr cached_observable = cacheable_observable->make_cached_observable(std::get<0>(c->second));
                    if (! cached_observable)
                        throw InternalError("make_cached_observable() failed");

                    // add the newly created cached observable
                    observables.push_back(cached_observable);
                    predictions.push_back(std::numeric_limits<double>::quiet_NaN());
                    cached_observables.push_back(std::make_tuple(cached_observable, index));

                    return index;
                }

                // else add this new cacheable observable
                observables.push_back(observable);
                predictions.push_back(std::numeric_limits<double>::quiet_NaN());
                cacheable_observables.insert(std::make_pair(type_index, std::make_tuple(cacheable_observable, index)));

                return index;
            }
            else
            {
                // add this new regular observable
                observables.push_back(observable);
                predictions.push_back(std::numeric_limits<double>::quiet_NaN());
                regular_observables.push_back(std::make_tuple(observable, index));

                return index;
            }

            throw InternalError("should not be reached");
        }
    };

    ObservableCache::ObservableCache(const Parameters & parameters) :
        PrivateImplementationPattern<ObservableCache>(new Implementation<ObservableCache>(parameters))
    {
    }

    ObservableCache::~ObservableCache()
    {
    }

    ObservableCache::Id
    ObservableCache::add(const ObservablePtr & observable)
    {
        return _imp->add(observable, *this);
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
            auto f = [=, this]() {
                auto & o   = std::get<0>(co.second);
                auto & idx = std::get<1>(co.second);
                try
                {
                    _imp->predictions[idx] = o->evaluate();
                }
                catch (eos::Exception & e)
                {
                    Log::instance()->message("ObservableCache::update", ll_error)
                        << "Exception encountered when evaluating cacheable observable '" << o->name() << "[" << o->kinematics().as_string() << "];" << o->options().as_string() << "': "
                        << e.what();
                    _imp->predictions[idx] = std::numeric_limits<double>::quiet_NaN();

                }
            };
            cacheable_tickets.push_back(ThreadPool::instance()->enqueue(std::function<void (void)>(f)));
        }

        std::vector<Ticket> regular_tickets;
        regular_tickets.reserve(_imp->regular_observables.size());

        // evaluate all regular observables in parallel
        for (auto ro : _imp->regular_observables)
        {
            auto f = [=, this]() {
                auto & o   = std::get<0>(ro);
                auto & idx = std::get<1>(ro);
                try
                {
                    _imp->predictions[idx] = o->evaluate();
                }
                catch (eos::Exception & e)
                {
                    Log::instance()->message("ObservableCache::update", ll_error)
                        << "Exception encountered when evaluating regular observable '" << o->name() << "[" << o->kinematics().as_string() << "];" << o->options().as_string() << "': "
                        << e.what();
                    _imp->predictions[idx] = std::numeric_limits<double>::quiet_NaN();

                }
            };
            regular_tickets.push_back(ThreadPool::instance()->enqueue(std::function<void (void)>(f)));
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
            auto f = [=, this]() {
                auto & o   = std::get<0>(co);
                auto & idx = std::get<1>(co);
                try
                {
                    _imp->predictions[idx] = o->evaluate();
                }
                catch (eos::Exception & e)
                {
                    Log::instance()->message("ObservableCache::update", ll_error)
                        << "Exception encountered when evaluating cached observable '" << o->name() << "[" << o->kinematics().as_string() << "];" << o->options().as_string() << "': "
                        << e.what();
                    _imp->predictions[idx] = std::numeric_limits<double>::quiet_NaN();

                }
            };
            cached_tickets.push_back(ThreadPool::instance()->enqueue(std::function<void (void)>(f)));
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
            auto & o   = std::get<0>(eo);
            auto & idx = std::get<1>(eo);
            try
            {
                _imp->predictions[idx] = o->evaluate();
            }
            catch (eos::Exception & e)
            {
                Log::instance()->message("ObservableCache::update", ll_error)
                    << "Exception encountered when evaluating expression observable '" << o->name() << "[" << o->kinematics().as_string() << "];" << o->options().as_string() << "': "
                    << e.what();
                _imp->predictions[idx] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }

    Parameters
    ObservableCache::parameters() const
    {
        return _imp->parameters;
    }

    double
    ObservableCache::operator[] (const ObservableCache::Id & id) const
    {
        return _imp->predictions[id];
    }

    ObservablePtr
    ObservableCache::observable(const ObservableCache::Id & id) const
    {
        return _imp->observables[id];
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

        for (auto o = _imp->observables.begin(), o_end = _imp->observables.end() ; o != o_end ; ++o)
        {
            // cloning cached observables creates independent *cacheable* observables
            // adding them back creates new and independent cached observables
            result._imp->add((*o)->clone(parameters), *this);
        }

        result.update();

        return result;
    }
}
