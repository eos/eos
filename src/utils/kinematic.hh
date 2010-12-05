/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_KINEMATIC_HH
#define EOS_GUARD_SRC_UTILS_KINEMATIC_HH 1

#include <src/utils/exception.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace eos
{
    struct UnknownKinematicVariableError :
        public Exception
    {
        UnknownKinematicVariableError(const std::string & variable) throw ();
    };

    class Kinematics :
        public PrivateImplementationPattern<Kinematics>
    {
        public:
            Kinematics();

            ~Kinematics();

            double operator[] (const std::string & variable) const;

            std::string as_string() const;

            void declare(const std::string & variable);

            void set(const std::string & variable, const double & value);
    };

    inline double lambda(const double & a, const double & b, const double & c)
    {
        return a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c);
    }
}

#endif
