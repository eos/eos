/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_CONCRETE_SIGNAL_PDF_HH
#define EOS_GUARD_EOS_UTILS_CONCRETE_SIGNAL_PDF_HH 1

#include <eos/signal-pdf.hh>
#include <eos/utils/apply.hh>
#include <eos/utils/density-impl.hh>
#include <eos/utils/tuple-maker.hh>

#include <functional>
#include <limits>
#include <string>
#include <vector>

#include <iostream>

namespace eos
{
    namespace impl
    {
        template <unsigned i_, typename Tuple_> struct DescriptionFiller
        {
            static std::vector<ParameterDescription> fill(Kinematics & k, const Tuple_ & kinematic_ranges)
            {
                auto result = DescriptionFiller<i_ - 1, Tuple_>::fill(k, kinematic_ranges);

                auto & kinematic_range = std::get<i_ - 1>(kinematic_ranges);
                std::string kname = std::string(kinematic_range.name);
                double kvalue = (kinematic_range.max - kinematic_range.min) / 2.0;
                MutablePtr kvar(new KinematicVariable(k.declare(kname, kvalue)));
                result.push_back(ParameterDescription{ kvar, kinematic_range.min, kinematic_range.max, false });

                return result;
            }
        };

        template <typename Tuple_> struct DescriptionFiller<1u, Tuple_>
        {
            static std::vector<ParameterDescription> fill(Kinematics & k, const Tuple_ & kinematic_ranges)
            {
                std::vector<ParameterDescription> result;

                auto & kinematic_range = std::get<0u>(kinematic_ranges);
                std::string kname = std::string(kinematic_range.name);
                double kvalue = (kinematic_range.max - kinematic_range.min) / 2.0;
                MutablePtr kvar(new KinematicVariable(k.declare(kname, kvalue)));
                result.push_back(ParameterDescription{ kvar, kinematic_range.min, kinematic_range.max, false });

                return std::move(result);
            }
        };

        template <unsigned i_, typename Tuple_> struct KinematicRangePrinter
        {
            static void print(std::ostream & os, const Tuple_ & kinematic_ranges)
            {
                KinematicRangePrinter<i_ - 1, Tuple_>::print(os, kinematic_ranges);
                auto & kinematic_range = std::get<i_ - 1>(kinematic_ranges);

                os << "    " << kinematic_range.name << "\t" << kinematic_range.description << std::endl;
            }
        };

        template <typename Tuple_> struct KinematicRangePrinter<1u, Tuple_>
        {
            static void print(std::ostream & os, const Tuple_ & kinematic_ranges)
            {
                auto & kinematic_range = std::get<0u>(kinematic_ranges);

                os << "    " << kinematic_range.name << "\t" << kinematic_range.description << std::endl;
            }
        };
    }

    template <typename Decay_, typename ... Args_>
    class ConcreteSignalPDF :
        public SignalPDF
    {
        public:

        private:
            std::string _name;

            Parameters _parameters;

            Kinematics _kinematics;

            std::vector<ParameterDescription> _descriptions;

            Options _options;

            Decay_ _decay;

            std::function<double (const Decay_ *, const Args_ & ...)> _function;

            std::tuple<typename impl::ConvertTo<Args_, KinematicRange>::Type ...> _kinematic_ranges;

            std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, KinematicVariable>::Type ...> _argument_tuple;

        public:
            ConcreteSignalPDF(const std::string & name,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & function,
                    const std::tuple<typename impl::ConvertTo<Args_, KinematicRange>::Type ...> & kinematic_ranges) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _descriptions(impl::DescriptionFiller<sizeof...(Args_), std::tuple<typename impl::ConvertTo<Args_, KinematicRange>::Type ...>>::fill(_kinematics, kinematic_ranges)),
                _options(options),
                _decay(parameters, options),
                _function(function),
                _kinematic_ranges(kinematic_ranges),
                _argument_tuple(impl::TupleMaker<sizeof...(Args_)>::make(_kinematics, _kinematic_ranges, &_decay))
            {
            }

            virtual const std::string & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values = _argument_tuple;

                double result = apply(_function, values);

                return (result > 0 ? std::log(result) : -std::numeric_limits<double>::max());
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

            virtual DensityPtr clone() const
            {
                return DensityPtr(new ConcreteSignalPDF(_name, _parameters.clone(), _kinematics.clone(), _options, _function, _kinematic_ranges));
            }

            virtual DensityPtr clone(const Parameters & parameters) const
            {
                return DensityPtr(new ConcreteSignalPDF(_name, parameters, _kinematics.clone(), _options, _function, _kinematic_ranges));
            }

            virtual Density::Iterator begin() const
            {
                return Density::Iterator(_descriptions.cbegin());
            }

            virtual Density::Iterator end() const
            {
                return Density::Iterator(_descriptions.cend());
            }
    };

    template <typename Decay_, typename ... Args_>
    class ConcreteSignalPDFEntry :
        public SignalPDFEntry
    {
        private:
            std::string _name;

            std::function<double (const Decay_ *, const Args_ & ...)> _function;

            std::tuple<typename impl::ConvertTo<Args_, KinematicRange>::Type ...> _kinematic_ranges;

        public:
            ConcreteSignalPDFEntry(const std::string & name,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & function,
                    const std::tuple<typename impl::ConvertTo<Args_, KinematicRange>::Type ...> & kinematic_ranges) :
                _name(name),
                _function(function),
                _kinematic_ranges(kinematic_ranges)
            {
            }

            ~ConcreteSignalPDFEntry()
            {
            }

            virtual SignalPDFPtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {

                return SignalPDFPtr(new ConcreteSignalPDF<Decay_, Args_ ...>(_name, parameters, kinematics, options, _function, _kinematic_ranges));
            }

            virtual std::ostream & insert(std::ostream & os) const
            {
                os << "    " << Decay_::description << std::endl;

                impl::KinematicRangePrinter<sizeof...(Args_), std::tuple<typename impl::ConvertTo<Args_, KinematicRange>::Type ...>>::print(os, _kinematic_ranges);

                return os;
            }
    };

    template <typename Decay_, typename Tuple_, typename ... Args_>
    SignalPDFEntry * make_concrete_signal_pdf_entry(const std::string & name, double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematic_ranges)
    {
        static_assert(sizeof...(Args_) == impl::TupleSize<Tuple_>::size, "Need as many function arguments as kinematics names!");
        return new ConcreteSignalPDFEntry<Decay_, Args_ ...>(name,
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(function)),
                kinematic_ranges);
    }
}

#endif
