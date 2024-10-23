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

#ifndef EOS_GUARD_SRC_CLIENTS_CLI_DUMPER_HH
#define EOS_GUARD_SRC_CLIENTS_CLI_DUMPER_HH 1

#include "cli_visitor.hh"
#include <iosfwd>

namespace eos
{
    namespace cli
    {
        class Option;

        class AliasArg;
        class EnumArg;
        class IntegerArg;
        class StringArg;
        class SwitchArg;

        /**
         * Prints help text appropriate to each command line option.
         */
        class Dumper
        {
            private:
                std::ostream & _os;

                void generic_visit(const Option &);

            public:
                /**
                 * Constructor.
                 */
                Dumper(std::ostream & os);

                /// Visit an AliasArg.
                void visit(const AliasArg &);

                /// Visit an EnumArg.
                void visit(const EnumArg &);

                /// Visit an IntegerArg.
                void visit(const IntegerArg &);

                /// Visit a KeyValueArg.
                void visit(const KeyValueArg &);

                /// Visit a StringArg.
                void visit(const StringArg &);

                /// Visit a StringListArg.
                void visit(const StringListArg &);

                /// Visit a SwitchArg.
                void visit(const SwitchArg &);
        };
    } // namespace cli
} // namespace eos
#endif
