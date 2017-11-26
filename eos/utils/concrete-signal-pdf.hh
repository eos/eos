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

        template <typename Result_, unsigned n_> struct ArgumentMaker
        {
            template <typename Ranges_, typename ... TempArgs_>
            static Result_ make(Kinematics & k, const Ranges_ & r, TempArgs_ ... args)
            {
                return ArgumentMaker<Result_, n_ - 1>::make(k, r, k[r[n_ - 1].name], args ...);
            }
        };

        template <typename Result_> struct ArgumentMaker<Result_, 0>
        {
            template <typename Ranges_, typename ... TempArgs_>
            static Result_ make(Kinematics &, const Ranges_ &, TempArgs_ ... args)
            {
                return Result_{ args ... };
            }
        };

        template <unsigned long n_> auto make_arguments(Kinematics & kinematics, const std::array<KinematicRange, n_> & kinematic_ranges)
             -> std::array<KinematicVariable, n_>
        {
            return ArgumentMaker<std::array<KinematicVariable, n_>, n_>::make(kinematics, kinematic_ranges);
        }

        template <unsigned long n_> std::vector<ParameterDescription> make_descriptions(Kinematics & kinematics, const std::array<KinematicRange, n_> & kinematic_ranges)
        {
            std::vector<ParameterDescription> result;

            for (auto krange : kinematic_ranges)
            {
                auto kname = std::string(krange.name);
                auto kvalue = (krange.max - krange.min) / 2.0;

                MutablePtr kvar(new KinematicVariable(kinematics.declare(kname, kvalue)));

                result.push_back(ParameterDescription{ kvar, krange.min, krange.max, false });
            }

            return result;
        }

        template <typename Result_, unsigned n_> struct Evaluator
        {
            template <typename Variables_, typename ... TempArgs_>
            static Result_ evaluate(const Variables_ & v, TempArgs_ ... args)
            {
                return Evaluator<Result_, n_ - 1>::evaluate(v, v[n_ - 1].evaluate(), args ...);
            }
        };

        template <typename Result_> struct Evaluator<Result_, 0>
        {
            template <typename Variables_, typename ... TempArgs_>
            static Result_ evaluate(const Variables_ &, TempArgs_ ... args)
            {
                return Result_{ args ... };
            }
        };

        template <unsigned long n_> auto evaluate(const std::array<KinematicVariable, n_> & kinematic_variables)
             -> std::array<double, n_>
        {
            return Evaluator<std::array<double, n_>, n_>::evaluate(kinematic_variables);
        }

    }

    template <typename Decay_, typename ... Args_>
    class ConcreteSignalPDF :
        public SignalPDF
    {
        public:

        private:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            std::vector<ParameterDescription> _descriptions;

            Options _options;

            Decay_ _decay;

            std::function<double (const Decay_ *, const Args_ & ...)> _function;

            std::array<KinematicRange, sizeof...(Args_)> _kinematic_ranges;

            std::array<KinematicVariable, sizeof...(Args_)> _arguments;

        public:
            ConcreteSignalPDF(const QualifiedName & name,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & function,
                    const std::array<KinematicRange, sizeof...(Args_)> & kinematic_ranges) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _descriptions(impl::make_descriptions(_kinematics, kinematic_ranges)),
                _options(options),
                _decay(parameters, options),
                _function(function),
                _kinematic_ranges(kinematic_ranges),
                _arguments(impl::make_arguments(_kinematics, kinematic_ranges))
            {
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                std::array<double, sizeof...(Args_)> arguments = impl::evaluate(_arguments);

                double result = apply(_function, &_decay, arguments);

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

    template <typename Decay_, typename ... FunctionArgs_>
    class ConcreteSignalPDFEntry :
        public SignalPDFEntry
    {
        private:
            QualifiedName _name;

            std::function<double (const Decay_ *, const FunctionArgs_ & ...)> _function;

            Options _default_options;

            std::array<KinematicRange, sizeof...(FunctionArgs_)> _kinematic_ranges;

        public:
            ConcreteSignalPDFEntry(const QualifiedName & name,
                    const std::function<double (const Decay_ *, const FunctionArgs_ & ...)> & function,
                    const Options & default_options,
                    const std::array<KinematicRange, sizeof...(FunctionArgs_)> & kinematic_ranges) :
                _name(name),
                _function(function),
                _default_options(default_options),
                _kinematic_ranges(kinematic_ranges)
            {
            }

            ~ConcreteSignalPDFEntry()
            {
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual const std::string & description() const
            {
                return Decay_::description;
            }

            virtual SignalPDFEntry::KinematicRangeIterator begin_kinematic_ranges() const
            {
                return _kinematic_ranges.begin();
            }

            virtual SignalPDFEntry::KinematicRangeIterator end_kinematic_ranges() const
            {
                return _kinematic_ranges.end();
            }

            virtual SignalPDFPtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {

                return SignalPDFPtr(new ConcreteSignalPDF<Decay_, FunctionArgs_ ...>(_name, parameters, kinematics, _default_options + options, _function, _kinematic_ranges));
            }

            virtual std::ostream & insert(std::ostream & os) const
            {
                return os;
            }
    };

    template <typename Decay_, typename ... FunctionArgs_, typename ... KinematicRanges_>
    SignalPDFEntry * make_concrete_signal_pdf_entry(const QualifiedName & name,
            double (Decay_::* function)(const FunctionArgs_ & ...) const,
            const Options & default_options,
            const KinematicRanges_ & ... kinematic_ranges)
    {
        static_assert(sizeof...(FunctionArgs_) == sizeof...(KinematicRanges_), "Need as many function arguments as kinematics ranges!");

        return new ConcreteSignalPDFEntry<Decay_, FunctionArgs_...>(name,
                std::function<double (const Decay_ *, const FunctionArgs_ & ...)>(std::mem_fn(function)),
                default_options,
                std::array<KinematicRange, sizeof...(FunctionArgs_)>{ kinematic_ranges... });
    }
}

#endif
