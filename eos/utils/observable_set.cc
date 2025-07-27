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

#include <eos/utils/observable_set.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <map>
#include <vector>

namespace eos
{
    template <> struct Implementation<ObservableSet>
    {
            // The list of observables
            std::vector<ObservablePtr> observables;

            // <observable name, index>
            std::map<QualifiedName, unsigned> observable_names;

            Implementation() {}

            std::pair<unsigned, bool>
            add(const ObservablePtr & observable)
            {
                if (! observables.empty() && (observable->parameters() != observables.front()->parameters()))
                {
                    throw InternalError("ObservableSet::add(): Mismatch of Parameters between different observables detected.");
                }

                // compare each observable for options, kinematics and name
                unsigned index = 0;
                for (auto i = observables.begin(), i_end = observables.end(); i != i_end; ++i, ++index)
                {
                    if (identical_observables(*i, observable))
                    {
                        return std::make_pair(index, false);
                    }
                }

                // new observable
                observables.push_back(observable);

                return std::make_pair(index, true);
            }

            static bool
            identical_observables(const ObservablePtr & lhs, const ObservablePtr & rhs)
            {
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
    };

    ObservableSet::ObservableSet() :
        PrivateImplementationPattern<ObservableSet>(new Implementation<ObservableSet>())
    {
    }

    ObservableSet::~ObservableSet() {}

    std::pair<unsigned, bool>
    ObservableSet::add(const ObservablePtr & observable)
    {
        return _imp->add(observable);
    }

    template <> struct WrappedForwardIteratorTraits<ObservableSet::IteratorTag>
    {
            using UnderlyingIterator = std::vector<ObservablePtr>::iterator;
    };
    template class WrappedForwardIterator<ObservableSet::IteratorTag, ObservablePtr>;

    ObservableSet::Iterator
    ObservableSet::begin() const
    {
        return Iterator(_imp->observables.begin());
    }

    ObservableSet::Iterator
    ObservableSet::end() const
    {
        return Iterator(_imp->observables.end());
    }

    ObservablePtr &
    ObservableSet::operator[] (const unsigned & index) const
    {
        return _imp->observables[index];
    }

    Parameters
    ObservableSet::parameters()
    {
        return _imp->observables.front()->parameters();
    }

    unsigned
    ObservableSet::size() const
    {
        return _imp->observables.size();
    }
} // namespace eos
