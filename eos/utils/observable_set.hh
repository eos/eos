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

#ifndef EOS_GUARD_EOS_UTILS_OBSERVABLE_SET_HH
#define EOS_GUARD_EOS_UTILS_OBSERVABLE_SET_HH 1

#include <eos/observable.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

namespace eos
{
    /*!
     * A container class to hold only unique observables,
     * i.e. they differ by at least one of the following:
     * name, kinematics, or options.
     */
    class ObservableSet : public PrivateImplementationPattern<ObservableSet>
    {
        public:
            ///@name Basic Functions
            ///@{
            /// Constructor.
            ObservableSet();

            /// Destructor.
            ~ObservableSet();
            ///@}

            ///@name Iteration and Access
            ///@{
            struct IteratorTag;
            using Iterator = WrappedForwardIterator<IteratorTag, ObservablePtr>;

            /// Iterator to the first observable.
            Iterator begin() const;

            /// Iterator pointing past the last observable.
            Iterator end() const;

            /*!
             * Random access to an observable.
             *
             * @param index The position in the vector.
             */
            ObservablePtr & operator[] (const unsigned & index) const;

            /*!
             * Add an observable to the vector.
             *
             * @param observable The observable to be added.
             *
             * @return If observable is found to be an existing observable,
             * the second value will be false. If the observable is added successfully,
             * the second value is true, and the index of the new element in the vector
             * is returned for use with operator[].
             */
            std::pair<unsigned, bool> add(const ObservablePtr & observable);

            /// Access to the underlying Parameters object.
            Parameters parameters();

            /// The total number of elements.
            unsigned size() const;
            ///@}
    };

    extern template class WrappedForwardIterator<ObservableSet::IteratorTag, ObservablePtr>;
} // namespace eos

#endif
