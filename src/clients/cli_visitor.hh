/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016, 2017 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2005-2011 Ciaran McCreesh
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

#ifndef EOS_GUARD_SRC_CLIENTS_CLI_VISITOR_HH
#define EOS_GUARD_SRC_CLIENTS_CLI_VISITOR_HH 1

#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/visitor.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include "cli_handler.hh"
#include <string>

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

        class Visitor : public PrivateImplementationPattern<cli::Visitor>
        {
            private:
                const std::string & get_param(const Option &);

            public:
                Visitor(Handler::ArgsIterator *, Handler::ArgsIterator, std::string & remaining_chars, bool no);

                ~Visitor();

                /// Visit an AliasArg.
                void visit(AliasArg &);

                /// Visit an EnumArg.
                void visit(EnumArg &);

                /// Visit an IntegerArg.
                void visit(IntegerArg &);

                /// Visit a KeyValueArg.
                void visit(KeyValueArg &);

                /// Visit a StringArg.
                void visit(StringArg &);

                /// Visit a StringListArg.
                void visit(StringListArg &);

                /// Visit a SwitchArg.
                void visit(SwitchArg &);
        };
    } // namespace cli
} // namespace eos

#endif
