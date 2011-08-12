/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

        // Same order as in predictions
        // <last result of evaluation>
        std::vector<double> predictions_dirty;

        // Store which observables(index list for observations) are needed for a parameter id
        std::map<unsigned, std::vector<unsigned>> index_lists;

        Implementation(const Parameters & parameters) :
            parameters(parameters)
        {
        }

        ~Implementation()
        {
        }

        ObservableCache::Id add(const ObservablePtr & observable)
        {
            if (! index_lists.empty())
                throw InternalError("ObservableCache::add(): Add all observables before evaluating the cache for the first time!");

            auto result = observables.add(observable);

            // new observable
            if (result.second)
            {
                predictions.push_back(std::numeric_limits<double>::quiet_NaN());
                predictions_dirty.push_back(std::numeric_limits<double>::quiet_NaN());
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

    void
    ObservableCache::update(const Parameter::Id & id)
    {
        // get list of observables for this id
        auto index_list = _imp->index_lists.find(id);

        // if parameter changes for the first time, determine the  list of observables for this id
        if (index_list == _imp->index_lists.end())
        {
            std::vector<unsigned> indices;

            // loop over all observables
            unsigned index = 0;
            for (auto o = _imp->observables.begin(), o_end = _imp->observables.end() ; o != o_end ; ++o, ++index)
            {
                // search through all parameters that this observable uses
                auto result = std::find((*o)->begin(), (*o)->end(), id);

                // no dependence on parameter
                if (result == (*o)->end())
                   continue;

                indices.push_back(index);
            }

            // add the index list
            index_list = _imp->index_lists.insert(std::make_pair(id, indices)).first;
        }

        // evaluate observables one by one
        for (auto o = index_list->second.begin(), o_end = index_list->second.end(); o != o_end; ++o)
        {
            // store last value of observable
            _imp->predictions_dirty[*o] = _imp->predictions[*o];

            // the new value
            _imp->predictions[*o] = _imp->observables[*o]->evaluate();
        }
    }

    Parameters
    ObservableCache::parameters() const
    {
        return _imp->parameters;
    }

    void
    ObservableCache::reset()
    {
        for (unsigned i = 0 ; i < _imp->predictions.size() ; ++i)
        {
            _imp->predictions[i] = _imp->predictions_dirty[i];
        }
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
