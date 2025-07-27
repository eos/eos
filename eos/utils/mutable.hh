/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_MUTABLE_HH
#define EOS_GUARD_EOS_UTILS_MUTABLE_HH 1

#include <eos/utils/mutable-fwd.hh>
#include <eos/utils/private_implementation_pattern.hh>

#include <memory>

namespace eos
{
    using MutablePtr = std::shared_ptr<Mutable>;

    /// Base class for all mutable entities, e.g. a Parameter.
    class Mutable
    {
        public:
            ///@name Basic Operatios
            ///@{
            /// Destructor.
            virtual ~Mutable() = 0;

            /// Make a copy of this Mutable.
            virtual MutablePtr clone() const = 0;

            ///@name Access & Modification of the numeric Value
            ///@{
            /// Cast a Mutable to a double.
            virtual operator double () const = 0;

            /// Retrieve a Mutable's numeric value.
            virtual double operator() () const = 0;

            /// Retrieve a Mutable's numeric value.
            // TODO: evaluate -> get
            virtual double evaluate() const = 0;

            /// Set a Mutable's numeric value.
            virtual const Mutable & operator= (const double &) = 0;

            /// Set a Mutable's numeric value.
            virtual void set(const double &) = 0;
            ///@}

            ///@name Access to Meta Data
            ///@{
            /// Retrieve the Parameter's name.
            virtual const std::string & name() const = 0;
            ///@}
    };
} // namespace eos

#endif
