/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016-2025 Danny van Dyk
 * Copyright (c) 2016      Rafael Silva Coutinho
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

#ifndef EOS_GUARD_EOS_UTILS_QUALIFIED_NAME_PARTS_HH
#define EOS_GUARD_EOS_UTILS_QUALIFIED_NAME_PARTS_HH 1

#include <string>

namespace eos
{
    namespace qnp
    {
        class Prefix
        {
            private:
                std::string _prefix;

            public:
                Prefix(const std::string &);

                const std::string & str() const { return _prefix; };

                inline bool operator<  (const Prefix & rhs) const { return this->_prefix <  rhs._prefix; };
                inline bool operator== (const Prefix & rhs) const { return this->_prefix == rhs._prefix; };
        };

        class Name
        {
            private:
                std::string _name;

            public:
                Name(const std::string &);

                const std::string & str() const { return _name; };

                inline bool operator<  (const Name & rhs) const { return this->_name <  rhs._name; };
                inline bool operator== (const Name & rhs) const { return this->_name == rhs._name; };
        };

        class Suffix
        {
            private:
                std::string _suffix;

            public:
                Suffix();
                Suffix(const std::string &);

                const std::string & str() const { return _suffix; };
                bool empty() const { return _suffix.empty(); };

                inline bool operator<  (const Suffix & rhs) const { return this->_suffix <  rhs._suffix; };
        };

        class OptionKey
        {
            private:
                std::string _key;

            public:
                OptionKey(const std::string &);

                const std::string & str() const { return _key; };
        };

        class OptionValue
        {
            private:
                std::string _value;

            public:
                OptionValue(const std::string &);

                const std::string & str() const { return _value; };
        };
    }
}

#endif
