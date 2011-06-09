/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#include <src/utils/unique_observable_vector.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/wrapped_forward_iterator-impl.hh>

#include <map>
#include <vector>

namespace eos
{
    template <>
    struct Implementation<UniqueObservableVector>
    {
        // The list of observables
        std::vector<ObservablePtr> observables;

        // <observable name, index>
        std::map<std::string, unsigned> observable_names;

        Implementation()
        {
        }

        std::pair<unsigned, bool> add(const ObservablePtr & observable)
        {
            if (! observables.empty() && (observable->parameters() != observables.front()->parameters()))
                throw InternalError("ObservableEvaluator::add(): Mismatch of Parameters between different observables detected.");

            // compare each observable for options, kinematics and name
            unsigned index = 0;
            for (auto i = observables.begin(), i_end = observables.end() ; i != i_end ; ++i, ++index)
            {
                if (identical_observables(*i, observable))
                    return std::make_pair(index, false);
            }

            // new observable
            observables.push_back(observable);

            return std::make_pair(index, true);
        }

        static bool identical_observables(const ObservablePtr & lhs, const ObservablePtr & rhs)
        {
            // compare name
            if( lhs->name() != rhs->name())
                return false;

            // compare kinematics
            if( lhs->kinematics() != rhs->kinematics())
                return false;

            // compare options
            if (lhs->options() != rhs->options())
                return false;

            return true;
        }
    };

    UniqueObservableVector::UniqueObservableVector() :
    PrivateImplementationPattern<UniqueObservableVector>(new Implementation<UniqueObservableVector>())
    {
    }

    UniqueObservableVector::~UniqueObservableVector()
    {
    }

    std::pair<unsigned, bool>
    UniqueObservableVector::add(const ObservablePtr & observable)
    {
        return _imp->add(observable);
    }

    template class WrappedForwardIterator<UniqueObservableVector::IteratorTag, ObservablePtr>;

    UniqueObservableVector::Iterator
    UniqueObservableVector::begin() const
    {
        return Iterator(_imp->observables.begin());
    }

    UniqueObservableVector::Iterator
    UniqueObservableVector::end() const
    {
        return Iterator(_imp->observables.end());
    }

    ObservablePtr &
    UniqueObservableVector::operator[] (const unsigned & index) const
    {
        return _imp->observables[index];
    }

    Parameters
    UniqueObservableVector::parameters()
    {
        return _imp->observables.front()->parameters();
    }

    unsigned
    UniqueObservableVector::size() const
    {
        return _imp->observables.size();
    }
}
