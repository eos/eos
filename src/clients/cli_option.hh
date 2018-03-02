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

#ifndef EOS_GUARD_SRC_CLIENTS_CLI_OPTION_HH
#define EOS_GUARD_SRC_CLIENTS_CLI_OPTION_HH 1

#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/named-value.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/type-list.hh>
#include <eos/utils/visitor.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <memory>
#include <vector>

namespace eos
{
    namespace n
    {
        typedef Name<struct name_description> description;
        typedef Name<struct name_long_name> long_name;
        typedef Name<struct name_short_name> short_name;
    }

    namespace cli
    {
        class Group;

        // Generic args
        class AliasArg;
        class EnumArg;
        class IntegerArg;
        class KeyValueArg;
        class LogLevelArg;
        class StringArg;
        class StringListArg;
        class SwitchArg;

        // EOS specific args
        class KinematicsSingleVariableArg;
        class KinematicsVariableRangeArg;
        class ObservableArg;
        class ParameterBudgetArg;
        class ParameterVariationArg;

        class Option :
            public virtual DeclareAbstractAcceptMethods<Option, MakeTypeList<
                    AliasArg, EnumArg, IntegerArg, KeyValueArg, StringArg, StringListArg, SwitchArg
                    >::Type>
        {
            friend class Handler;

            private:
                Group * const _group;

                const std::string _long_name;
                const char _short_name;
                const std::string _description;

                bool _specified;

                Option(const Option &);
                void operator= (const Option &);

            protected:
                /**
                 * Constructor.
                 */
                Option(Group * const, const std::string & long_name,
                        const char short_name, const std::string & description);

                /**
                 * Destructor.
                 */
                virtual ~Option();

            public:
                /**
                 * Remove this option.  Removes our group from its
                 * section if the group would be left empty.
                 */
                void remove();

                /**
                 * Fetch our long name.
                 */
                const std::string & long_name() const
                {
                    return _long_name;
                }

                /**
                 * Fetch our short name (may be 0).
                 */
                char short_name() const
                {
                    return _short_name;
                }

                /**
                 * Fetch our description.
                 */
                const std::string & description() const
                {
                    return _description;
                }

                /**
                 * Fetch whether or not we were specified on the
                 * command line (or as an env var).
                 */
                virtual bool specified() const
                {
                    return _specified;
                }

                /**
                 * Set the value returned by specified().
                 */
                virtual void set_specified(bool value)
                {
                    _specified = value;
                }

                /**
                 * Fetch our group.
                 */
                Group * group()
                {
                    return _group;
                }

                /**
                 * Can we be negated?
                 *
                 * Needs to match up with Visitor logic.
                 */
                virtual bool can_be_negated() const = 0;
        };

        /**
         * An AliasArg is an alias for another argument.
         */
        class AliasArg :
            public Option,
            public ImplementAcceptMethods<Option, AliasArg>
        {
            private:
                Option * const _other;
                bool _hidden;

            public:
                /**
                 * Constructor.
                 */
                AliasArg(Option * const other, const std::string & new_long_name, bool is_hidden = false);

                /**
                 * Destructor.
                 */
                virtual ~AliasArg();

                virtual bool specified() const
                {
                    return _other->specified();
                }

                virtual void set_specified(const bool value)
                {
                    _other->set_specified(value);
                }

                virtual bool hidden() const
                {
                    return _hidden;
                }

                virtual void set_hidden(const bool value)
                {
                    _hidden = value;
                }

                /**
                 * Fetch our associated option.
                 */
                Option * other() const
                {
                    return _other;
                }

                virtual bool can_be_negated() const;
        };

        /**
         * A SwitchArg is an option that can either be specified or not
         * specified, and that takes no value (for example, --help).
         */
        class SwitchArg :
            public Option,
            public ImplementAcceptMethods<Option, SwitchArg>
        {
            private:
                bool _can_be_negated;

            public:
                /**
                 * Constructor.
                 */
                SwitchArg(Group * const group, const std::string & long_name, char short_name,
                        const std::string & description, const bool can_be_negated);

                /**
                 * Destructor.
                 */
                virtual ~SwitchArg();

                virtual bool can_be_negated() const;
        };

        /**
         * An option that takes a string argument.
         */
        class StringArg :
            public Option,
            public ImplementAcceptMethods<Option, StringArg>
        {
            private:
                std::string _argument;
                bool _can_be_negated;
                void (* _validator) (const std::string &);

            public:
                /**
                 * Constructor.
                 */
                StringArg(Group * const, const std::string & long_name,
                       const char short_name, const std::string & description,
                       const bool can_be_negated = false);

                /**
                 * Constructor with validator.
                 */
                StringArg(Group * const, const std::string & long_name,
                       const char short_name, const std::string & description,
                       void (* validator) (const std::string &),
                       const bool can_be_negated = false);

                /**
                 * Destructor.
                 */
                virtual ~StringArg();

                /**
                 * Fetch the argument that was given to this option.
                 */
                const std::string & argument() const { return _argument; }

                /**
                 * Set the argument returned by argument().
                 */
                void set_argument(const std::string & arg);

                virtual bool can_be_negated() const;
        };

        /**
         * An option that takes a list of strings.
         */
        class StringListArg :
            public Option,
            public ImplementAcceptMethods<Option, StringListArg>,
            public PrivateImplementationPattern<StringListArg>
        {
            private:
                void (* _validator) (const std::string &);

            public:
                /**
                 * Constructor
                 */
                StringListArg(Group * const, const std::string & long_name,
                        const char short_name, const std::string & description);

                /**
                 * Constructor with validator.
                 */
                StringListArg(Group * const, const std::string & long_name,
                        const char short_name, const std::string & description,
                        void (* validator) (const std::string &));

                /**
                 * Destructor.
                 */
                virtual ~StringListArg();

                /**
                 * Iterate over our args.
                 */
                struct ConstIteratorTag;
                typedef WrappedForwardIterator<ConstIteratorTag, const std::string> ConstIterator;

                ConstIterator begin_args() const;
                ConstIterator end_args() const;

                /**
                 * Add an argument to the list.
                 */
                void validate_and_add_argument(const std::string & arg);

                virtual bool can_be_negated() const;
        };

        /**
         * An option that takes an integer argument.
         */
        class IntegerArg :
            public Option,
            public ImplementAcceptMethods<Option, IntegerArg>
        {
            private:
                int _argument;

            public:
                /**
                 * Constructor
                 */
                IntegerArg(Group * const, const std::string & long_name,
                        const char short_name, const std::string & description);

                /**
                 * Destructor.
                 */
                virtual ~IntegerArg();

                /**
                 * Fetch the argument that was given to this option.
                 */
                int argument() const { return _argument; }

                /**
                 * Set the argument returned by argument().
                 */
                void set_argument(const int arg) { _argument = arg; }

                virtual bool can_be_negated() const;
        };

        /**
         * An option that takes a key and a value.
         */
        class KeyValueArg :
            public Option,
            public ImplementAcceptMethods<Option, KeyValueArg>
        {
            public:
                /**
                 * Constructor
                 */
                KeyValueArg(Group * const, const std::string & long_name,
                        const char short_name, const std::string & description);

                /**
                 * Destructor.
                 */
                virtual ~KeyValueArg();

                /**
                 * Validate the correctness of key and value, and set the arguments if
                 * validated.
                 */
                virtual void validate_and_set_arguments(const std::string & key, const std::string & value) = 0;

                virtual bool can_be_negated() const { return false; }
        };

        /**
         * An allowed argument for an EnumArg.
         */
        struct AllowedEnumArg
        {
            NamedValue<n::description, std::string> description;
            NamedValue<n::long_name, std::string> long_name;

            /// Might be '@0', for none.
            NamedValue<n::short_name, char> short_name;
        };

        /**
         * An option that takes one of a predefined set of string arguments.
         */
        class EnumArg :
            public Option,
            public ImplementAcceptMethods<Option, EnumArg>,
            public PrivateImplementationPattern<EnumArg>
        {
            private:
                std::string _argument;
                std::string _default_arg;

            public:
                /**
                 * Helper class for passing available options and associated descriptions
                 * to the EnumArg constructor.
                 */
                class EnumArgOptions :
                    public PrivateImplementationPattern<EnumArgOptions>
                {
                    friend class EnumArg;

                    public:
                        /**
                         * Constructor
                         */
                        EnumArgOptions(const std::string &, const std::string &);

                        /**
                         * Constructor, with short arg.
                         *
                         * @since 0.40
                         */
                        EnumArgOptions(const std::string &, const char, const std::string &);

                        /**
                         * Destructor.
                         */
                        ~EnumArgOptions();

                        /**
                         * Adds another (option, description).
                         */
                        EnumArgOptions & operator() (const std::string &, const std::string &);

                        /**
                         * Adds another (option, short-option, description).
                         */
                        EnumArgOptions & operator() (const std::string &, const char, const std::string &);
                };

                /**
                 * Constructor.
                 */
                EnumArg(Group * const group, const std::string & long_name,
                        const char short_name, const std::string & description,
                        const EnumArgOptions & opts, const std::string & default_arg);

                /**
                 * Destructor.
                 */
                virtual ~EnumArg();

                /**
                 * Fetch the argument that was given to this option.
                 */
                const std::string & argument() const
                {
                    return _argument;
                }

                /**
                 * Set the argument returned by argument(), having verified that
                 * it is one of the arguments allowed for this option.
                 */
                void set_argument(const std::string & arg);

                /**
                 * Change the default option (should be called before
                 * set_argument()).
                 */
                void set_default_arg(const std::string & arg);

                /**
                 * Fetch the default option, as specified to the
                 * constructor or set_default_arg().
                 */
                const std::string & default_arg() const
                {
                    return _default_arg;
                }

                ///@name Iterate over our allowed arguments and associated descriptions
                ///@{

                struct AllowedArgConstIteratorTag;
                typedef WrappedForwardIterator<AllowedArgConstIteratorTag,
                        const AllowedEnumArg> AllowedArgConstIterator;

                AllowedArgConstIterator begin_allowed_args() const;

                AllowedArgConstIterator end_allowed_args() const;

                ///@}

                virtual bool can_be_negated() const;
        };

        /**
         * The '--log-level' standard command line argument.
         */
        class LogLevelArg :
            public EnumArg
        {
            public:
                ///@name Basic operations
                ///@{

                LogLevelArg(Group * const, const std::string &, char);
                virtual ~LogLevelArg();

                ///@}

                /**
                 * Our selected value, as a LogLevel.
                 */
                LogLevel option() const;
        };

        /**
         * The '--kinematic-variable' EOS-specific command line argument.
         */
        class KinematicVariableArg :
            public KeyValueArg
        {
            private:
                /// Our set of kinematic variables.
                Kinematics _kinematics;

            public:
                ///@name Basic operations
                ///@{

                KinematicVariableArg(Group * const, const std::string &, char, const Kinematics &);
                virtual ~KinematicVariableArg();

                ///@}

                /**
                 * Validate the correctness of the key (KinematicVariable name) and value (double),
                 * and set the arguments if validated.
                 */
                virtual void validate_and_set_arguments(const std::string & key, const std::string & value);

                /**
                 * Our current set of kinematic variables, as Kinematics.
                 */
                Kinematics kinematics() const;
        };

        /**
         * The '--parameter-budget' EOS-specific command line argument.
         */
        class ParameterBudgetArg :
            public StringArg
        {
            public:
                struct ParameterBudget
                {
                    std::string name;
                    std::vector<Parameter> parameters;
                };

            private:
                /// Our set of parameters.
                Parameters _parameters;

                /// Our list of parameter budgets.
                std::vector<ParameterBudget> _budgets;

            public:
                ///@name Basic operations
                ///@{

                ParameterBudgetArg(Group * const, const std::string &, char, const std::string &, const Parameters &);
                virtual ~ParameterBudgetArg();

                ///@}

                ///@name Iterate over our list of parameter budgets.
                ///@{

                struct ParameterBudgetArgConstIteratorTag;
                typedef WrappedForwardIterator<ParameterBudgetArgConstIteratorTag,
                        const ParameterBudget> ParameterBudgetArgConstIterator;

                ParameterBudgetArgConstIterator begin_budgets() const;

                ParameterBudgetArgConstIterator end_budgets() const;


                ///@}
        };
    }

    extern template class WrappedForwardIterator<cli::EnumArg::AllowedArgConstIteratorTag, const cli::AllowedEnumArg>;

    extern template class WrappedForwardIterator<cli::ParameterBudgetArg::ParameterBudgetArgConstIteratorTag, const cli::ParameterBudgetArg::ParameterBudget>;
}

#endif
