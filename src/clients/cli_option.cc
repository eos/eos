/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016, 2017 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2005-2010 Ciaran McCreesh
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

#include "cli_error.hh"
#include "cli_group.hh"
#include "cli_handler.hh"
#include "cli_option.hh"
#include "cli_section.hh"

#include <eos/utils/destringify.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <set>
#include <vector>
#include <algorithm>
#include <list>

namespace
{
    struct ArgIs
    {
        const std::string arg;

        ArgIs(const std::string & a) :
            arg(a)
        {
        }

        bool operator() (const std::pair<std::string, std::string> & p) const
        {
            return p.first == arg;
        }

        bool operator() (const eos::cli::AllowedEnumArg & p) const
        {
            return p.long_name() == arg || (p.short_name() && std::string(1, p.short_name()) == arg);
        }
    };
}

namespace eos
{
    template <>
    struct Implementation<cli::EnumArg>
    {
        std::vector<cli::AllowedEnumArg> allowed_args;
    };

    template <>
    struct Implementation<cli::EnumArg::EnumArgOptions>
    {
        std::vector<cli::AllowedEnumArg> options;
    };

    template <>
    struct WrappedForwardIteratorTraits<cli::EnumArg::AllowedArgConstIteratorTag>
    {
        typedef std::vector<cli::AllowedEnumArg>::const_iterator UnderlyingIterator;
    };

    template class WrappedForwardIterator<cli::EnumArg::AllowedArgConstIteratorTag, const cli::AllowedEnumArg>;

    template <>
    struct Implementation<cli::StringListArg>
    {
        std::vector<std::string> args;
    };

    template <>
    struct WrappedForwardIteratorTraits<cli::StringListArg::ConstIteratorTag>
    {
        typedef std::vector<std::string>::const_iterator UnderlyingIterator;
    };

    template class WrappedForwardIterator<cli::StringListArg::ConstIteratorTag, const std::string>;

    namespace cli
    {
        // Option
        Option::Option(Group * const g, const std::string & our_long_name,
                const char our_short_name, const std::string & our_description) :
            _group(g),
            _long_name(our_long_name),
            _short_name(our_short_name),
            _description(our_description),
            _specified(false)
        {
            g->add(this);
            g->section()->handler()->add_option(this, our_long_name, our_short_name);
        }

        Option::~Option() = default;

        void
        Option::remove()
        {
            _group->remove(this);
            _group->section()->handler()->remove_option(_long_name, _short_name);
        }

        // AliasArg
        AliasArg::AliasArg(Option * const o, const std::string & our_long_name, bool is_hidden) :
            Option(o->group(), our_long_name, '\0', "Alias for --" + o->long_name()),
            _other(o), _hidden(is_hidden)
        {
        }

        AliasArg::~AliasArg() = default;

        bool
        AliasArg::can_be_negated() const
        {
            return _other->can_be_negated();
        }

        // EnumArg
        EnumArg::EnumArg(Group * const our_group, const std::string & our_long_name,
                const char our_short_name, const std::string & our_description,
                const EnumArgOptions & opts, const std::string & our_default_arg) :
            cli::Option(our_group, our_long_name, our_short_name, our_description),
            PrivateImplementationPattern<cli::EnumArg>(new Implementation<cli::EnumArg>()),
            _argument(our_default_arg),
            _default_arg(our_default_arg)
        {
            _imp->allowed_args = opts._imp->options;
        }

        EnumArg::~EnumArg() = default;

        bool
        EnumArg::can_be_negated() const
        {
            return false;
        }

        void
        EnumArg::set_argument(const std::string & arg)
        {
//            Context context("When handling argument '" + arg + "' for '--" + long_name() + "':");

            /* if we're given the short arg, turn it magically into the long one */
            AllowedArgConstIterator i(std::find_if(_imp->allowed_args.begin(), _imp->allowed_args.end(), ArgIs(arg)));
            if (i == _imp->allowed_args.end())
                throw BadValue("--" + long_name(), arg);

            _argument = i->long_name();
        }

        void
        EnumArg::set_default_arg(const std::string & arg)
        {
            _argument = arg;
            _default_arg = arg;
        }

        EnumArg::AllowedArgConstIterator
        EnumArg::begin_allowed_args() const
        {
            return AllowedArgConstIterator(_imp->allowed_args.begin());
        }

        EnumArg::AllowedArgConstIterator
        EnumArg::end_allowed_args() const
        {
            return AllowedArgConstIterator(_imp->allowed_args.end());
        }

        // EnumArg::EnumArgOptions
        EnumArg::EnumArgOptions::EnumArgOptions(const std::string & opt, const std::string & desc) :
            PrivateImplementationPattern<cli::EnumArg::EnumArgOptions>(new Implementation<cli::EnumArg::EnumArgOptions>())
        {
            _imp->options.push_back(make_named_values<AllowedEnumArg>(
                        n::description() = desc,
                        n::long_name() = opt,
                        n::short_name() = '\0'
                        ));
        }

        EnumArg::EnumArgOptions::EnumArgOptions(const std::string & opt, const char s, const std::string & desc) :
            PrivateImplementationPattern<cli::EnumArg::EnumArgOptions>(new Implementation<cli::EnumArg::EnumArgOptions>())
        {
            _imp->options.push_back(make_named_values<AllowedEnumArg>(
                        n::description() = desc,
                        n::long_name() = opt,
                        n::short_name() = s
                        ));
        }

        EnumArg::EnumArgOptions::~EnumArgOptions() = default;

        EnumArg::EnumArgOptions &
        EnumArg::EnumArgOptions::operator() (const std::string & opt, const std::string & desc)
        {
            _imp->options.push_back(make_named_values<AllowedEnumArg>(
                        n::description() = desc,
                        n::long_name() = opt,
                        n::short_name() = '\0'
                        ));
            return *this;
        }

        EnumArg::EnumArgOptions &
        EnumArg::EnumArgOptions::operator() (const std::string & opt, const char s, const std::string & desc)
        {
            _imp->options.push_back(make_named_values<AllowedEnumArg>(
                        n::description() = desc,
                        n::long_name() = opt,
                        n::short_name() = s
                        ));
            return *this;
        }

        // IntegerArg
        IntegerArg::IntegerArg(Group * const our_group, const std::string & our_long_name,
                        char our_short_name, const std::string & our_description) :
            Option(our_group, our_long_name, our_short_name, our_description)
        {
        }

        IntegerArg::~IntegerArg() = default;

        bool
        IntegerArg::can_be_negated() const
        {
            return false;
        }

        // KeyValueArg
        KeyValueArg::KeyValueArg(Group * const our_group, const std::string & our_long_name,
                        char our_short_name, const std::string & our_description) :
            Option(our_group, our_long_name, our_short_name, our_description)
        {
        }

        KeyValueArg::~KeyValueArg() = default;

        // StringArg
        StringArg::StringArg(Group * const g, const std::string & our_long_name,
                const char our_short_name, const std::string & our_description,
                const bool neg) :
            Option(g, our_long_name, our_short_name, our_description),
            _can_be_negated(neg),
            _validator(nullptr)
        {
        }

        StringArg::StringArg(Group * const g, const std::string & our_long_name,
                const char our_short_name, const std::string & our_description,
                void (* v) (const std::string &), const bool neg) :
            Option(g, our_long_name, our_short_name, our_description),
            _can_be_negated(neg),
            _validator(v)
        {
        }

        StringArg::~StringArg() = default;

        bool
        StringArg::can_be_negated() const
        {
            return _can_be_negated;
        }

        void
        StringArg::set_argument(const std::string & arg)
        {
//            Context context("When handling argument '" + arg + "' for '--" + long_name() + "':");

            if (_validator)
                (*_validator)(arg);

            _argument = arg;
        }

        // StringListArg
        StringListArg::StringListArg(Group * const g, const std::string & our_long_name,
                const char our_short_name, const std::string & our_description) :
            Option(g, our_long_name, our_short_name, our_description),
            PrivateImplementationPattern<cli::StringListArg>(new Implementation<cli::StringListArg>()),
            _validator(nullptr)
        {
        }

        StringListArg::StringListArg(Group * const g, const std::string & our_long_name,
                const char our_short_name, const std::string & our_description,
                void (* v) (const std::string &)) :
            Option(g, our_long_name, our_short_name, our_description),
            PrivateImplementationPattern<cli::StringListArg>(new Implementation<cli::StringListArg>()),
            _validator(v)
        {
        }

        StringListArg::~StringListArg() = default;

        bool
        StringListArg::can_be_negated() const
        {
            return false;
        }

        void
        StringListArg::validate_and_add_argument(const std::string & arg)
        {
//            Context context("When handling argument '" + arg + "' for '--" + long_name() + "':");
//
            if (_validator)
                (*_validator)(arg);

            _imp->args.push_back(arg);
        }

        StringListArg::ConstIterator
        StringListArg::begin_args() const
        {
            return _imp->args.cbegin();
        }

        StringListArg::ConstIterator
        StringListArg::end_args() const
        {
            return _imp->args.cend();
        }

        // SwitchArg
        SwitchArg::SwitchArg(Group * const our_group, const std::string & our_long_name, char our_short_name,
                const std::string & our_description, const bool c) :
            Option(our_group, our_long_name, our_short_name, our_description),
            _can_be_negated(c)
        {
        }

        SwitchArg::~SwitchArg() = default;

        bool
        SwitchArg::can_be_negated() const
        {
            return _can_be_negated;
        }
    }

    template <>
    struct WrappedForwardIteratorTraits<cli::ParameterBudgetArg::ParameterBudgetArgConstIteratorTag>
    {
        typedef std::vector<cli::ParameterBudgetArg::ParameterBudget>::const_iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<cli::ParameterBudgetArg::ParameterBudgetArgConstIteratorTag, const cli::ParameterBudgetArg::ParameterBudget>;


    namespace cli
    {
        // LogLevelArg
        LogLevelArg::LogLevelArg(Group * const our_group, const std::string & our_long_name, char our_short_name) :
            EnumArg(our_group, our_long_name, our_short_name, "Specify the log level",
            EnumArgOptions
            ("debug",   'd', "Show debug output (noisy)")
            ("info",    'i', "Show informations and and warnings only")
            ("warning", 'w', "Show warnings only")
            ("error",   'e', "Show errors only")
            ("silent",  's', "Suppress all log messages (UNSAFE)"),
            "info")
        {
        }

        LogLevelArg::~LogLevelArg() = default;

        LogLevel
        LogLevelArg::option() const
        {
            if ("debug" == argument())
                return ll_debug;
            if ("info" == argument())
                return ll_informational;
            if ("warning" == argument())
                return ll_warning;
            if ("error" == argument())
                return ll_error;
            if ("silent" == argument())
                return ll_silent;

            throw DoHelp("Bad value for --" + long_name());
        }

        // KinematicVariableArg
        KinematicVariableArg::KinematicVariableArg(Group * const our_group, const std::string & our_long_name,
                char our_short_name, const Kinematics & k) :
            KeyValueArg(our_group, our_long_name, our_short_name, "Set the value of a kinematic variable"),
            _kinematics(k)
        {
        }

        KinematicVariableArg::~KinematicVariableArg() = default;

        void
        KinematicVariableArg::validate_and_set_arguments(const std::string & key, const std::string & value)
        {
            try
            {
                _kinematics.declare(key, destringify<double>(value));
            }
            catch (DestringifyError & e)
            {
                throw DoHelp("Bad value for --" + long_name());
            }
        }

        Kinematics
        KinematicVariableArg::kinematics() const
        {
            return _kinematics;
        }

        // ParameterBudgetArg
        ParameterBudgetArg::ParameterBudgetArg(Group * const our_group, const std::string & our_long_name,
                char our_short_name, const std::string & our_description, const Parameters & p) :
            StringArg(our_group, our_long_name, our_short_name, our_description),
            _parameters(p)
        {
        }

        ParameterBudgetArg::~ParameterBudgetArg() = default;

        ParameterBudgetArg::ParameterBudgetArgConstIterator
        ParameterBudgetArg::begin_budgets() const
        {
            return _budgets.begin();
        }

        ParameterBudgetArg::ParameterBudgetArgConstIterator
        ParameterBudgetArg::end_budgets() const
        {
            return _budgets.end();
        }
    }
}
