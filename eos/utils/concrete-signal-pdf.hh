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
#include <eos/utils/density-impl.hh>
#include <eos/utils/tuple-maker.hh>

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace eos
{
    namespace impl
    {
        template <unsigned i_, typename Tuple_> struct DescriptionFiller
        {
                static std::vector<ParameterDescription>
                fill(Kinematics & k, const Tuple_ & kinematic_ranges)
                {
                    auto result = DescriptionFiller<i_ - 1, Tuple_>::fill(k, kinematic_ranges);

                    auto &      kinematic_range = std::get<i_ - 1>(kinematic_ranges);
                    std::string kname           = std::string(kinematic_range.name);
                    double      kvalue          = (kinematic_range.max - kinematic_range.min) / 2.0;
                    MutablePtr  kvar(new KinematicVariable(k.declare(kname, kvalue)));
                    result.push_back(ParameterDescription{ kvar, kinematic_range.min, kinematic_range.max, false });

                    return result;
                }
        };

        template <typename Tuple_> struct DescriptionFiller<1u, Tuple_>
        {
                static std::vector<ParameterDescription>
                fill(Kinematics & k, const Tuple_ & kinematic_ranges)
                {
                    std::vector<ParameterDescription> result;

                    auto &      kinematic_range = std::get<0u>(kinematic_ranges);
                    std::string kname           = std::string(kinematic_range.name);
                    double      kvalue          = (kinematic_range.max - kinematic_range.min) / 2.0;
                    MutablePtr  kvar(new KinematicVariable(k.declare(kname, kvalue)));
                    result.push_back(ParameterDescription{ kvar, kinematic_range.min, kinematic_range.max, false });

                    return std::move(result);
                }
        };

        template <unsigned i_, typename Tuple_> struct KinematicRangePrinter
        {
                static void
                print(std::ostream & os, const Tuple_ & kinematic_ranges)
                {
                    KinematicRangePrinter<i_ - 1, Tuple_>::print(os, kinematic_ranges);
                    auto & kinematic_range = std::get<i_ - 1>(kinematic_ranges);

                    os << "    " << kinematic_range.name << "\t" << kinematic_range.description << std::endl;
                }
        };

        template <typename Tuple_> struct KinematicRangePrinter<1u, Tuple_>
        {
                static void
                print(std::ostream & os, const Tuple_ & kinematic_ranges)
                {
                    auto & kinematic_range = std::get<0u>(kinematic_ranges);

                    os << "    " << kinematic_range.name << "\t" << kinematic_range.description << std::endl;
                }
        };

        template <typename Result_, unsigned n_> struct ArgumentFromRange
        {
                template <typename Ranges_, typename... TempArgs_>
                static Result_
                make(Kinematics & k, const Ranges_ & r, TempArgs_... args)
                {
                    return ArgumentFromRange<Result_, n_ - 1>::make(k, r, k[r[n_ - 1].name], args...);
                }
        };

        template <typename Result_> struct ArgumentFromRange<Result_, 0>
        {
                template <typename Ranges_, typename... TempArgs_>
                static Result_
                make(Kinematics &, const Ranges_ &, TempArgs_... args)
                {
                    return Result_{ args... };
                }
        };

        template <unsigned long n_>
        auto
        make_arguments(Kinematics & kinematics, const std::array<KinematicRange, n_> & kinematic_ranges) -> std::array<KinematicVariable, n_>
        {
            return ArgumentFromRange<std::array<KinematicVariable, n_>, n_>::make(kinematics, kinematic_ranges);
        }

        template <typename Result_, unsigned n_> struct ArgumentFromName
        {
                template <typename Names_, typename... TempArgs_>
                static Result_
                make(Kinematics & k, const Names_ & n, TempArgs_... args)
                {
                    return ArgumentFromName<Result_, n_ - 1>::make(k, n, k[n[n_ - 1]], args...);
                }
        };

        template <typename Result_> struct ArgumentFromName<Result_, 0>
        {
                template <typename Names_, typename... TempArgs_>
                static Result_
                make(Kinematics &, const Names_ &, TempArgs_... args)
                {
                    return Result_{ args... };
                }
        };

        template <unsigned long n_>
        auto
        make_arguments(Kinematics & kinematics, const std::array<std::string, n_> & kinematic_names) -> std::array<KinematicVariable, n_>
        {
            return ArgumentFromName<std::array<KinematicVariable, n_>, n_>::make(kinematics, kinematic_names);
        }

        template <unsigned long n_>
        std::vector<ParameterDescription>
        make_descriptions(Kinematics & kinematics, const std::array<KinematicRange, n_> & kinematic_ranges)
        {
            std::vector<ParameterDescription> result;

            for (auto krange : kinematic_ranges)
            {
                auto kname  = std::string(krange.name);
                auto kvalue = (krange.max - krange.min) / 2.0;

                MutablePtr kvar(new KinematicVariable(kinematics.declare(kname, kvalue)));

                result.push_back(ParameterDescription{ kvar, krange.min, krange.max, false });
            }

            return result;
        }

        template <typename Result_, unsigned n_> struct Evaluator
        {
                template <typename Variables_, typename... TempArgs_>
                static Result_
                evaluate(const Variables_ & v, TempArgs_... args)
                {
                    return Evaluator<Result_, n_ - 1>::evaluate(v, v[n_ - 1].evaluate(), args...);
                }
        };

        template <typename Result_> struct Evaluator<Result_, 0>
        {
                template <typename Variables_, typename... TempArgs_>
                static Result_
                evaluate(const Variables_ &, TempArgs_... args)
                {
                    return Result_{ args... };
                }
        };

        template <unsigned long n_>
        auto
        evaluate(const std::array<KinematicVariable, n_> & kinematic_variables) -> std::array<double, n_>
        {
            return Evaluator<std::array<double, n_>, n_>::evaluate(kinematic_variables);
        }

        // convert std::tuple to std::array
        template <int... Indices_> struct indices
        {
                using next = indices<Indices_..., sizeof...(Indices_)>;
        };

        template <int Size_> struct build_indices
        {
                using type = typename build_indices<Size_ - 1>::type::next;
        };

        template <> struct build_indices<0>
        {
                using type = indices<>;
        };

        template <typename T_> using Bare = typename std::remove_cv<typename std::remove_reference<T_>::type>::type;

        template <typename Tuple_>
        constexpr typename build_indices<std::tuple_size<Bare<Tuple_>>::value>::type
        make_indices()
        {
            return {};
        }

        template <typename Tuple_, int... Indices_>
        std::array<typename std::tuple_element<0, Bare<Tuple_>>::type, std::tuple_size<Bare<Tuple_>>::value>
        to_array(Tuple_ && tuple, indices<Indices_...>)
        {
            using std::get;
            return { { get<Indices_>(std::forward<Tuple_>(tuple))... } };
        }

        template <typename Tuple_>
        auto
        to_array(Tuple_ && tuple) -> decltype(to_array(std::declval<Tuple_>(), make_indices<Tuple_>()))
        {
            return to_array(std::forward<Tuple_>(tuple), make_indices<Tuple_>());
        }

        // convert Decay_ + std::array to std::tuple such that std::apply can be used
        template <typename T_, std::size_t N_, std::size_t... indices_>
        auto
        make_tuple(const T_ * t, const std::array<double, N_> & a, std::index_sequence<indices_...>)
        {
            return std::make_tuple(t, std::get<indices_>(a)...);
        }

        template <typename T_, std::size_t N_>
        auto
        convert_to_tuple(const T_ * t, const std::array<double, N_> & args)
        {
            return make_tuple(t, args, std::make_index_sequence<N_>());
        }
    } // namespace impl

    template <typename Decay_, typename PDFSignature_, unsigned pdf_args_, typename NormSignature_, unsigned norm_args_> class ConcreteSignalPDF : public SignalPDF
    {
        public:
        private:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            std::vector<ParameterDescription> _descriptions;

            Options _options;

            Decay_ _decay;

            std::function<PDFSignature_> _pdf;

            std::array<KinematicRange, pdf_args_> _pdf_kinematic_ranges;

            std::array<KinematicVariable, pdf_args_> _pdf_arguments;

            std::function<NormSignature_> _norm;

            std::array<std::string, norm_args_> _norm_kinematic_names;

            std::array<KinematicVariable, norm_args_> _norm_arguments;

        public:
            ConcreteSignalPDF(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options,
                              const std::function<PDFSignature_> & pdf, const std::array<KinematicRange, pdf_args_> & pdf_kinematic_ranges,
                              const std::function<NormSignature_> & norm, const std::array<std::string, norm_args_> & norm_kinematic_names) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _descriptions(impl::make_descriptions(_kinematics, pdf_kinematic_ranges)),
                _options(options),
                _decay(parameters, options),
                _pdf(pdf),
                _pdf_kinematic_ranges(pdf_kinematic_ranges),
                _pdf_arguments(impl::make_arguments(_kinematics, pdf_kinematic_ranges)),
                _norm(norm),
                _norm_kinematic_names(norm_kinematic_names),
                _norm_arguments(impl::make_arguments(_kinematics, norm_kinematic_names))
            {
            }

            virtual const QualifiedName &
            name() const
            {
                return _name;
            }

            virtual double
            evaluate() const
            {
                std::array<double, pdf_args_> pdf_arguments = impl::evaluate(_pdf_arguments);

                double result = std::apply(_pdf, impl::convert_to_tuple(&_decay, pdf_arguments));

                return (result > 0 ? std::log(result) : -std::numeric_limits<double>::max());
            }

            virtual double
            normalization() const
            {
                std::array<double, norm_args_> norm_arguments = impl::evaluate(_norm_arguments);

                double result = std::apply(_norm, impl::convert_to_tuple(&_decay, norm_arguments));

                return (result > 0 ? std::log(result) : -std::numeric_limits<double>::max());
            }

            virtual Parameters
            parameters()
            {
                return _parameters;
            }

            virtual Kinematics
            kinematics()
            {
                return _kinematics;
            }

            virtual Options
            options()
            {
                return _options;
            }

            virtual DensityPtr
            clone() const
            {
                return DensityPtr(new ConcreteSignalPDF(_name, _parameters.clone(), _kinematics.clone(), _options, _pdf, _pdf_kinematic_ranges, _norm, _norm_kinematic_names));
            }

            virtual DensityPtr
            clone(const Parameters & parameters) const
            {
                return DensityPtr(new ConcreteSignalPDF(_name, parameters, _kinematics.clone(), _options, _pdf, _pdf_kinematic_ranges, _norm, _norm_kinematic_names));
            }

            virtual Density::Iterator
            begin() const
            {
                return Density::Iterator(_descriptions.cbegin());
            }

            virtual Density::Iterator
            end() const
            {
                return Density::Iterator(_descriptions.cend());
            }
    };

    template <typename Decay_, typename PDFSignature_, unsigned pdf_args_, typename NormSignature_, unsigned norm_args_> class ConcreteSignalPDFEntry : public SignalPDFEntry
    {
        private:
            QualifiedName _name;

            Options _default_options;

            std::function<PDFSignature_> _pdf;

            std::array<KinematicRange, pdf_args_> _pdf_kinematic_ranges;

            std::function<NormSignature_> _norm;

            std::array<std::string, norm_args_> _norm_kinematic_names;

        public:
            ConcreteSignalPDFEntry(const QualifiedName & name, const Options & default_options, const std::function<PDFSignature_> & pdf,
                                   const std::array<KinematicRange, pdf_args_> & pdf_kinematic_ranges, const std::function<NormSignature_> & norm,
                                   const std::array<std::string, norm_args_> & norm_kinematic_names) :
                _name(name),
                _default_options(default_options),
                _pdf(pdf),
                _pdf_kinematic_ranges(pdf_kinematic_ranges),
                _norm(norm),
                _norm_kinematic_names(norm_kinematic_names)
            {
            }

            ~ConcreteSignalPDFEntry() {}

            virtual const QualifiedName &
            name() const
            {
                return _name;
            }

            virtual const std::string &
            description() const
            {
                return Decay_::description;
            }

            virtual SignalPDFEntry::KinematicRangeIterator
            begin_kinematic_ranges() const
            {
                return _pdf_kinematic_ranges.begin();
            }

            virtual SignalPDFEntry::KinematicRangeIterator
            end_kinematic_ranges() const
            {
                return _pdf_kinematic_ranges.end();
            }

            virtual SignalPDFPtr
            make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                return SignalPDFPtr(new ConcreteSignalPDF<Decay_, PDFSignature_, pdf_args_, NormSignature_, norm_args_>(_name,
                                                                                                                        parameters,
                                                                                                                        kinematics,
                                                                                                                        _default_options + options,
                                                                                                                        _pdf,
                                                                                                                        _pdf_kinematic_ranges,
                                                                                                                        _norm,
                                                                                                                        _norm_kinematic_names));
            }

            virtual std::ostream &
            insert(std::ostream & os) const
            {
                return os;
            }
    };

    template <typename Decay_, typename... PDFArgs_, typename... PDFKinematicRanges_, typename... NormArgs_, typename... NormKinematicNames_>
    SignalPDFEntry *
    make_concrete_signal_pdf_entry(const QualifiedName & name, const Options & default_options, const std::function<double(const Decay_ *, const PDFArgs_ &...)> & pdf,
                                   const std::tuple<PDFKinematicRanges_...> & _pdf_kinematic_ranges, const std::function<double(const Decay_ *, const NormArgs_ &...)> & norm,
                                   const std::tuple<NormKinematicNames_...> & _norm_kinematic_names)
    {
        static_assert(sizeof...(PDFArgs_) == sizeof...(PDFKinematicRanges_), "Need as many function arguments as kinematics ranges!");
        static_assert(sizeof...(NormArgs_) == sizeof...(NormKinematicNames_), "Need as many function arguments for the normalization as kinematics names!");

        using PDFSignature  = double(const Decay_ *, const PDFArgs_ &...);
        using NormSignature = double(const Decay_ *, const NormArgs_ &...);
        std::array<KinematicRange, sizeof...(PDFArgs_)> pdf_kinematic_ranges(impl::to_array(_pdf_kinematic_ranges));
        std::array<std::string, sizeof...(NormArgs_)>   norm_kinematic_names;
        std::array<const char *, sizeof...(NormArgs_)>  __norm_kinematic_names(impl::to_array(_norm_kinematic_names));
        for (auto i = 0u; i < sizeof...(NormArgs_); ++i)
        {
            norm_kinematic_names[i] = std::string(__norm_kinematic_names[i]);
        }

        return new ConcreteSignalPDFEntry<Decay_, PDFSignature, sizeof...(PDFArgs_), NormSignature, sizeof...(NormArgs_)>(name,
                                                                                                                          default_options,
                                                                                                                          pdf,
                                                                                                                          pdf_kinematic_ranges,
                                                                                                                          norm,
                                                                                                                          norm_kinematic_names);
    }
} // namespace eos

#endif
