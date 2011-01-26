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

#ifndef EOS_GUARD_SRC_UTILS_OBSERVABLE_HH
#define EOS_GUARD_SRC_UTILS_OBSERVABLE_HH 1

#include <src/utils/exception.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/parameters.hh>

#include <string>
#include <tr1/memory>

namespace eos
{
    class Observable;

    class ObservableOptions;

    typedef std::shared_ptr<Observable> ObservablePtr;

    class Observable
    {
        public:
            Observable();

            virtual ~Observable();

            virtual const std::string & name() const = 0;

            virtual double evaluate(const Kinematics &) const = 0;

            virtual Parameters parameters() = 0;

            virtual ObservableOptions options() = 0;

            virtual ObservablePtr clone() const = 0;
    };

    typedef std::shared_ptr<Observable> ObservablePtr;

    struct UnknownOptionError :
        public Exception
    {
        UnknownOptionError(const std::string & key) throw ();
    };

    class ObservableOptions :
        public PrivateImplementationPattern<ObservableOptions>
    {
        public:
            ObservableOptions();

            ~ObservableOptions();

            const std::string & operator[] (const std::string & key) const;

            bool has(const std::string & key) const;

            void set(const std::string & key, const std::string & value = "");

            std::string get(const std::string & key, const std::string & default_value = "") const;

            std::string as_string() const;
    };

    class ObservableFactory
    {
        public:
            ObservableFactory();

            virtual ~ObservableFactory();

            virtual ObservablePtr make(const Parameters &, const ObservableOptions &) const = 0;
    };
}

#endif
