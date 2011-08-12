/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_OBSERVABLE_HH
#define EOS_GUARD_SRC_UTILS_OBSERVABLE_HH 1

#include <eos/utils/exception.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>

#include <string>
#include <memory>

namespace eos
{
    class Observable;

    class Options;

    typedef std::shared_ptr<Observable> ObservablePtr;

    class Observable :
        public ParameterUser
    {
        public:
            virtual const std::string & name() const = 0;

            virtual double evaluate() const = 0;

            virtual Kinematics kinematics() = 0;

            virtual Parameters parameters() = 0;

            virtual Options options() = 0;

            virtual ObservablePtr clone() const = 0;

            virtual ObservablePtr clone(const Parameters & parameters) const = 0;

            static ObservablePtr make(const std::string & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options);
    };

    typedef std::shared_ptr<Observable> ObservablePtr;

    class ObservableFactory
    {
        public:
            ObservableFactory();

            virtual ~ObservableFactory();

            virtual ObservablePtr make(const Parameters &, const Kinematics &, const Options &) const = 0;
    };

    /*!
     * ObservableNameError is thrown when Observable::make encounters a malformed observable name.
     */
    struct ObservableNameError :
        public Exception
    {
        ///@name Basic Functions
        ///@{
        /*!
         * Constructor.
         *
         * @param name The offending malformed observable name.
         */
        ObservableNameError(const std::string & name);
        ///@}
    };
}

#endif
