/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016, 2017 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2006 Stephen Bennett
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

#include "cli_dumper.hh"

#include "cli_option.hh"
#include <sstream>

namespace eos
{
    namespace cli
    {
        Dumper::Dumper(std::ostream & os) :
            _os(os)
        {
        }

        void
        Dumper::generic_visit(const Option & a)
        {
            std::stringstream p;
            p << "  --" << a.long_name();
            if (a.short_name())
            {
                p << ", -" << a.short_name();
            }
            if (p.str().length() < 24)
            {
                p << std::string(24 - p.str().length(), ' ');
            }
            else
            {
                p << std::endl << std::string(24, ' ');
            }
            _os << p.str();
            _os << " " << a.description() << std::endl;
        }

        void
        Dumper::visit(const AliasArg & a)
        {
            if (! a.hidden())
            {
                generic_visit(a);
            }
        }

        void
        Dumper::visit(const EnumArg & a)
        {
            generic_visit(a);

            for (EnumArg::AllowedArgConstIterator it = a.begin_allowed_args(), it_end = a.end_allowed_args(); it != it_end; ++it)
            {
                std::stringstream p;
                p << "      " << it->long_name();
                if (it->short_name())
                {
                    p << " (" << std::string(1, it->short_name()) << ")";
                }
                if (p.str().length() < 26)
                {
                    p << std::string(26 - p.str().length(), ' ');
                }
                _os << p.str();
                _os << " " << it->description();
                if (it->long_name() == a.default_arg())
                {
                    _os << " (default)";
                }
                _os << std::endl;
            }
        }

        void
        Dumper::visit(const IntegerArg & a)
        {
            generic_visit(a);
        }

        void
        Dumper::visit(const KeyValueArg & a)
        {
            generic_visit(a);
        }

        void
        Dumper::visit(const StringArg & a)
        {
            generic_visit(a);
        }

        void
        Dumper::visit(const StringListArg & a)
        {
            generic_visit(a);
        }

        void
        Dumper::visit(const SwitchArg & a)
        {
            generic_visit(a);
        }
    } // namespace cli
} // namespace eos
