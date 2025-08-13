/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013-2025 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_FORM_FACTOR_ADAPTER_HH
#define EOS_GUARD_EOS_FORM_FACTORS_FORM_FACTOR_ADAPTER_HH 1

#include <eos/observable.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/utils/tuple-maker.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <tuple>

namespace eos
{
    /* Form factor adapter class for interfacing Observable */
    template <typename Transition_, typename ... Args_>
    class FormFactorAdapter :
        public Observable
    {
        private:
            QualifiedName _name;

            qnp::Prefix _process;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            SpecifiedOption _opt_form_factors;

            std::shared_ptr<FormFactors<Transition_>> _form_factors;

            std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> _form_factor_function;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::tuple<const FormFactors<Transition_> *, typename impl::ConvertTo<Args_, KinematicVariable>::Type ...> _argument_tuple;

        public:
            FormFactorAdapter(const QualifiedName & name,
                    const qnp::Prefix & process,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> & form_factor_function,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _process(process),
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _opt_form_factors(options, FormFactorFactory<Transition_>::option_specification(process)),
                _form_factors(FormFactorFactory<Transition_>::create(process.str() + "::" + _opt_form_factors.value(), _parameters, _options)),
                _form_factor_function(form_factor_function),
                _kinematics_names(kinematics_names),
                _argument_tuple(impl::TupleMaker<sizeof...(Args_)>::make(_kinematics, _kinematics_names, _form_factors.get()))
            {
                if (! _form_factors)
                    throw NoSuchFormFactorError(process.str(), options["form-factors"_ok]);

                uses(*_form_factors);
                auto _register_kinematics = [this](const FormFactors<Transition_> *, typename impl::ConvertTo<Args_, KinematicVariable>::Type... args)
                {
                    std::array<const KinematicVariable, sizeof...(Args_)> kinematics_array = { args... };
                    for (const auto & kinematic_variable : kinematics_array)
                    {
                        this->uses_kinematic(kinematic_variable.id());
                    }
                };
                std::apply(_register_kinematics, _argument_tuple);
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                std::tuple<const FormFactors<Transition_> *, typename impl::ConvertTo<Args_, double>::Type ...> values = _argument_tuple;

                return std::apply(_form_factor_function, values);
            };

            virtual Parameters parameters()
            {
                return _parameters;
            };

            virtual Kinematics kinematics()
            {
                return _kinematics;
            };

            virtual Options options()
            {
                return _options;
            }

            virtual ObservablePtr clone() const
            {
                return ObservablePtr(new FormFactorAdapter(_name, _process, _parameters.clone(), _kinematics.clone(), _options, _form_factor_function, _kinematics_names));
            }

            virtual ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new FormFactorAdapter(_name, _process, parameters, _kinematics.clone(), _options, _form_factor_function, _kinematics_names));
            }
    };

    template <typename Transition_, typename ... Args_>
    class FormFactorAdapterEntry :
        public ObservableEntry
    {
        private:
            QualifiedName _name;

            std::string _latex;

            Unit _unit;

            qnp::Prefix _process;

            std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> _form_factor_function;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::array<const std::string, sizeof...(Args_)> _kinematics_names_array;

            const std::vector<OptionSpecification> _options;

        public:
            FormFactorAdapterEntry(const QualifiedName & name,
                    const std::string & latex,
                    const Unit & unit,
                    const qnp::Prefix & process,
                    const std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> & form_factor_function,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _latex(latex),
                _unit(unit),
                _process(process),
                _form_factor_function(form_factor_function),
                _kinematics_names(kinematics_names),
                _kinematics_names_array(impl::make_array<const std::string>(kinematics_names)),
                _options{ FormFactorFactory<Transition_>::option_specification(process) }
            {
            }

            ~FormFactorAdapterEntry()
            {
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual const std::string & latex() const
            {
                return _latex;
            }

            virtual const Unit & unit() const
            {
                return _unit;
            }

            virtual ObservableEntry::KinematicVariableIterator begin_kinematic_variables() const
            {
                return _kinematics_names_array.begin();
            }

            virtual ObservableEntry::KinematicVariableIterator end_kinematic_variables() const
            {
                return _kinematics_names_array.end();
            }

            virtual ObservableEntry::OptionIterator begin_options() const
            {
                return _options.begin();
            }

            virtual ObservableEntry::OptionIterator end_options() const
            {
                return _options.end();
            }

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                return ObservablePtr(new FormFactorAdapter<Transition_, Args_ ...>(_name, _process, parameters, kinematics, options, _form_factor_function, _kinematics_names));
            }
    };
}

#endif
