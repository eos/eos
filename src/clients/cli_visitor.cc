/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016, 2017 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2005-2011 Ciaran McCreesh
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

#include "cli_visitor.hh"

#include <eos/utils/destringify.hh>
#include <eos/utils/log.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include "cli_error.hh"
#include "cli_handler.hh"
#include "cli_option.hh"
#include <algorithm>
#include <sstream>
#include <stdlib.h>

using namespace eos;
using namespace eos::cli;

namespace eos
{
    template <> struct Implementation<cli::Visitor>
    {
            Handler::ArgsIterator * args_index;
            Handler::ArgsIterator   args_end;
            std::string &           remaining_chars;
            bool                    no;

            Implementation(Handler::ArgsIterator * i, Handler::ArgsIterator e, std::string & s, bool n) :
                args_index(i),
                args_end(e),
                remaining_chars(s),
                no(n)
            {
            }
    };

    namespace cli
    {
        Visitor::Visitor(Handler::ArgsIterator * ai, Handler::ArgsIterator ae, std::string & s, bool n) :
            PrivateImplementationPattern<cli::Visitor>(new Implementation<cli::Visitor>(ai, ae, s, n))
        {
        }

        Visitor::~Visitor() = default;

        const std::string &
        Visitor::get_param(const Option & arg)
        {
            if (++(*_imp->args_index) == _imp->args_end)
            {
                throw MissingValue("--" + arg.long_name());
            }

            return **_imp->args_index;
        }

        void
        Visitor::visit(EnumArg & arg)
        {
            if (arg.specified())
            {
                Log::instance()->message("args.specified_twice", ll_warning)
                        << "Option '--" << arg.long_name() << "' was specified more than once, but it does not take multiple values";
            }

            arg.set_specified(true);

            std::string p;
            if (_imp->remaining_chars.length() == 1)
            {
                p = _imp->remaining_chars;
                _imp->remaining_chars.clear();
            }
            else
            {
                p = get_param(arg);
            }

            arg.set_argument(p);
        }

        void
        Visitor::visit(IntegerArg & arg)
        {
            if (arg.specified())
            {
                Log::instance()->message("args.specified_twice", ll_warning)
                        << "Option '--" << arg.long_name() << "' was specified more than once, but it does not take multiple values";
            }

            arg.set_specified(true);
            std::string param;

            if ((! _imp->remaining_chars.empty()) && (std::string::npos == _imp->remaining_chars.find_first_not_of("0123456789")))
            {
                param = _imp->remaining_chars;
                _imp->remaining_chars.clear();
            }
            else
            {
                param = get_param(arg);
            }

            try
            {
                int a(destringify<int>(param));
                arg.set_argument(a);
            }
            catch (const DestringifyError &)
            {
                throw BadValue("--" + arg.long_name(), param);
            }
        }

        void
        Visitor::visit(KeyValueArg & arg)
        {
            std::string key(get_param(arg));
            std::string value(get_param(arg));
            arg.validate_and_set_arguments(key, value);
        }

        void
        Visitor::visit(StringArg & arg)
        {
            if (arg.specified())
            {
                Log::instance()->message("args.specified_twice", ll_warning)
                        << "Option '--" << arg.long_name() << "' was specified more than once, but it does not take multiple values";
            }

            std::string p(get_param(arg));
            arg.set_specified(true);
            arg.set_argument(p);
        }

        void
        Visitor::visit(StringListArg & arg)
        {
            std::string p(get_param(arg));
            arg.validate_and_add_argument(p);
        }

        void
        Visitor::visit(AliasArg & arg)
        {
            arg.other()->accept(*this);
        }

        void
        Visitor::visit(SwitchArg & arg)
        {
            arg.set_specified(true);
        }
    } // namespace cli
} // namespace eos
