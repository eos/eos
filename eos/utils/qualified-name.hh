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

#ifndef EOS_GUARD_EOS_UTILS_QUALIFIED_NAME_HH
#define EOS_GUARD_EOS_UTILS_QUALIFIED_NAME_HH 1

#include <eos/utils/exception.hh>
#include <eos/utils/options.hh>
#include <eos/utils/qualified-name-parts.hh>

#include <string>
#include <vector>

namespace eos
{
    /*
     * Holds a syntactically-correct qualified name.
     *
     * The syntax of a qualified name follows:
     *
     *   PREFIX '::' NAME
     *   PREFIX '::' NAME '@' SUFFIX
     *   PREFIX '::' NAME ';' OPTIONSLIST
     *   PREFIX '::' NAME '@' SUFFIX ';' OPTIONSLIST
     *
     * with
     *
     *   PREFIX      := ['a'-'z', 'A'-'Z', '0'-'9', '<', '>', '^', '_', '*', '+', '-', '(', ')']
     *                  but "<<" and ">>" are prohibited (they delimit observables for the expression parser)
     *   NAME        := ['a'-'z', 'A'-'Z', '0'-'9', '(', ')', '[', ']', '{', '}', '|',
     *                   '\'', '.', ',', '/', '^', '_', '*', '+', '-']
     *   SUFFIX      := ['a'-'z', 'A'-'Z', '0'-'9', '.', ':', '-', '(', ')']
     *   OPTIONSLIST := OPTION | OPTIONSLIST ',' OPTION
     *   OPTION      := KEY '=' VALUE
     *   KEY         := ['a'-'z', 'A'-'Z', '0'-'9', '-']
     *   VALUE       := ['a'-'z', 'A'-'Z', '0'-'9', '+', '-', '/', '.', '^', '_']
     */
    class QualifiedName
    {
            friend std::ostream & operator<< (std::ostream &, const QualifiedName &);

        private:
            std::string _str;  // short hand name, excluding possible options
            std::string _full; // full name, including all given options
            qnp::Prefix _prefix;
            qnp::Name   _name;
            qnp::Suffix _suffix;
            Options     _options;

        public:
            QualifiedName(const std::string & name);
            QualifiedName(const char * name);
            QualifiedName(const QualifiedName & other);
            QualifiedName(const qnp::Prefix & prefix, const qnp::Name & name, const qnp::Suffix & suffix = qnp::Suffix());
            ~QualifiedName();

            inline const std::string &
            str() const
            {
                return _str;
            }

            inline const std::string &
            full() const
            {
                return _full;
            }

            inline const qnp::Prefix &
            prefix_part() const
            {
                return _prefix;
            }

            inline const qnp::Name &
            name_part() const
            {
                return _name;
            }

            inline const qnp::Suffix &
            suffix_part() const
            {
                return _suffix;
            }

            inline const Options &
            options() const
            {
                return _options;
            }

            /*
             * Two qualified names are compared based on their short names only.
             * As a consequence, two qualified names can be identical, even if their
             * full names aren't.
             */
            inline bool
            operator< (const QualifiedName & rhs) const
            {
                return this->_str < rhs._str;
            }

            inline bool
            operator== (const QualifiedName & rhs) const
            {
                return this->_str == rhs._str;
            }

            inline bool
            operator!= (const QualifiedName & rhs) const
            {
                return this->_str != rhs._str;
            }
    };

    class QualifiedNameSyntaxError : public Exception
    {
        public:
            QualifiedNameSyntaxError(const std::string & msg);
    };

    inline std::ostream &
    operator<< (std::ostream & lhs, const QualifiedName & rhs)
    {
        lhs << rhs._str;

        return lhs;
    }
} // namespace eos

#endif
