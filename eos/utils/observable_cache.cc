/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2016 Danny van Dyk
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

#include <eos/utils/observable_cache.hh>
#include <eos/utils/observable_set.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <algorithm>
#include <limits>
#include <map>
#include <tuple>
#include <vector>

namespace eos
{
    template <> struct
    Implementation<ObservableCache>
    {
        // Parameters which are common to all observables in the cache.
        Parameters parameters;

        // Store each observable that needs to be calculated exactly once
        ObservableSet observables;

        // Store values of observables
        std::vector<double> predictions;

        Implementation(const Parameters & parameters) :
            parameters(parameters)
        {
        }

        ~Implementation()
        {
        }

        ObservableCache::Id add(const ObservablePtr & observable)
        {
            auto result = observables.add(observable);

            // new observable
            if (result.second)
            {
                predictions.push_back(std::numeric_limits<double>::quiet_NaN());
            }

            return result.first;
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
        return _imp->add(observable);
    }

    void
    ObservableCache::update()
    {
        // evaluate all observables
        auto p = _imp->predictions.begin();

        for (auto o = _imp->observables.begin(), o_end = _imp->observables.end() ; o != o_end ; ++o, ++p)
        {
            *p = (*o)->evaluate();
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

    ObservableCache
    ObservableCache::clone(const Parameters & parameters) const
    {
        ObservableCache result(parameters);

        for (auto o = _imp->observables.begin(), o_end = _imp->observables.end() ; o != o_end ; ++o)
        {
            result._imp->add((*o)->clone(parameters));
        }

        result.update();

        return result;
    }
}
