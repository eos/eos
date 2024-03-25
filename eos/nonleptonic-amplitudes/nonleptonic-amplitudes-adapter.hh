/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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

#ifndef EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_NONLEPTONIC_AMPLITUDES_ADAPTER_HH
#define EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_NONLEPTONIC_AMPLITUDES_ADAPTER_HH 1

#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes.hh>
#include <eos/observable.hh>
#include <eos/utils/tuple-maker.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <tuple>

namespace eos
{
    /* Amplitude adapter class for interfacing Observable */
    template <typename Transition_, typename ... Args_>
    class NonleptonicAmplitudesAdapter :
        public Observable
    {
        private:
            QualifiedName _name;

            qnp::Prefix _process;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            std::shared_ptr<NonleptonicAmplitudes<Transition_>> _nonleptonic_amplitudes;

            std::function<double (const NonleptonicAmplitudes<Transition_> *, const Args_ & ...)> _nonleptonic_amplitudes_function;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::tuple<const NonleptonicAmplitudes<Transition_> *, typename impl::ConvertTo<Args_, KinematicVariable>::Type ...> _argument_tuple;

        public:
            NonleptonicAmplitudesAdapter(const QualifiedName & name,
                    const qnp::Prefix & process,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<double (const NonleptonicAmplitudes<Transition_> *, const Args_ & ...)> & nonleptonic_amplitudes_function,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _process(process),
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _nonleptonic_amplitudes(NonleptonicAmplitudeFactory<Transition_>::create(process.str() + "::" + options["representation"], _parameters, _options)),
                _nonleptonic_amplitudes_function(nonleptonic_amplitudes_function),
                _kinematics_names(kinematics_names),
                _argument_tuple(impl::TupleMaker<sizeof...(Args_)>::make(_kinematics, _kinematics_names, _nonleptonic_amplitudes.get()))
            {
                if (! _options.has("representation"))
                    throw UnknownOptionError("representation");

                if (! _nonleptonic_amplitudes)
                    throw NoSuchNonleptonicAmplitudeError(process.str(), options["representation"]);

                uses(*_nonleptonic_amplitudes);
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                std::tuple<const NonleptonicAmplitudes<Transition_> *, typename impl::ConvertTo<Args_, double>::Type ...> values = _argument_tuple;

                return std::apply(_nonleptonic_amplitudes_function, values);
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
                return ObservablePtr(new NonleptonicAmplitudesAdapter(_name, _process, _parameters.clone(), _kinematics.clone(), _options, _nonleptonic_amplitudes_function, _kinematics_names));
            }

            virtual ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new NonleptonicAmplitudesAdapter(_name, _process, parameters, _kinematics.clone(), _options, _nonleptonic_amplitudes_function, _kinematics_names));
            }
    };

    template <typename Transition_, typename ... Args_>
    class NonleptonicAmplitudesAdapterEntry :
        public ObservableEntry
    {
        private:
            QualifiedName _name;

            std::string _latex;

            Unit _unit;

            qnp::Prefix _process;

            std::function<double (const NonleptonicAmplitudes<Transition_> *, const Args_ & ...)> _nonleptonic_amplitudes_function;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::array<const std::string, sizeof...(Args_)> _kinematics_names_array;

            const std::vector<OptionSpecification> _options;

        public:
            NonleptonicAmplitudesAdapterEntry(const QualifiedName & name,
                    const std::string & latex,
                    const Unit & unit,
                    const qnp::Prefix & process,
                    const std::function<double (const NonleptonicAmplitudes<Transition_> *, const Args_ & ...)> & nonleptonic_amplitudes_function,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _latex(latex),
                _unit(unit),
                _process(process),
                _nonleptonic_amplitudes_function(nonleptonic_amplitudes_function),
                _kinematics_names(kinematics_names),
                _kinematics_names_array(impl::make_array<const std::string>(kinematics_names)),
                _options{ NonleptonicAmplitudeFactory<Transition_>::option_specification(process) }
            {
            }

            ~NonleptonicAmplitudesAdapterEntry()
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
                return ObservablePtr(new NonleptonicAmplitudesAdapter<Transition_, Args_ ...>(_name, _process, parameters, kinematics, options, _nonleptonic_amplitudes_function, _kinematics_names));
            }

            virtual std::ostream & insert(std::ostream & os) const
            {
                os << "    type: adapter for one nonloeptonic amplitude" << std::endl;

                return os;
            }
    };
}

#endif
