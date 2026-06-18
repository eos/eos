/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2026 Danny van Dyk
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

#include <eos/observable.hh>
#include <eos/signal-pdf.hh>
#include <eos/utils/density-impl.hh>
#include <eos/utils/tuple-maker.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

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
        std::vector<typename std::tuple_element<0, Bare<Tuple_>>::type>
        to_vector(Tuple_ && tuple, indices<Indices_...>)
        {
            using std::get;
            return { { get<Indices_>(std::forward<Tuple_>(tuple))... } };
        }

        template <typename Tuple_>
        auto
        to_vector(Tuple_ && tuple) -> decltype(to_vector(std::declval<Tuple_>(), make_indices<Tuple_>()))
        {
            return to_vector(std::forward<Tuple_>(tuple), make_indices<Tuple_>());
        }

#if 0
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
#endif
    } // namespace impl

    class ConcreteSignalPDF : public SignalPDF
    {
        public:
        private:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            ObservablePtr _unnormalized_pdf;

            ObservablePtr _normalization;

            // TODO: remove
            std::vector<ParameterDescription> _descriptions;

        public:
            ConcreteSignalPDF(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options,
                              const QualifiedName & unnormalized_pdf, const QualifiedName & normalization);

            virtual const QualifiedName & name() const;

            virtual double evaluate() const;

            virtual double evaluate_linear() const;

            virtual double normalization() const;

            virtual Parameters parameters();

            virtual Kinematics kinematics();

            virtual Options options();

            virtual DensityPtr clone() const;

            virtual DensityPtr clone(const Parameters & parameters) const;

            virtual Density::Iterator begin() const;

            virtual Density::Iterator end() const;
    };

    class ConcreteSignalPDFEntry : public SignalPDFEntry
    {
        private:
            QualifiedName _name;

            std::string _description;

            Options _default_options;

            QualifiedName _numerator;

            QualifiedName _normalization;

            std::vector<std::string> _numerator_kinematic_names;

            std::vector<std::string> _normalization_kinematic_names;

        public:
            ConcreteSignalPDFEntry(const QualifiedName & name, const std::string & description, const Options & default_options, const QualifiedName & numerator,
                                   const QualifiedName & normalization, const std::vector<std::string> & numerator_kinematic_names,
                                   const std::vector<std::string> & normalization_kinematic_names);
            ~ConcreteSignalPDFEntry();

            virtual const QualifiedName &                              name() const;
            virtual const std::string &                                description() const;
            virtual SignalPDFEntry::NumeratorKinematicVariableIterator begin_numerator_kinematic_variables() const;

            virtual SignalPDFEntry::NumeratorKinematicVariableIterator end_numerator_kinematic_variables() const;

            virtual SignalPDFEntry::DenominatorKinematicVariableIterator begin_denominator_kinematic_variables() const;

            virtual SignalPDFEntry::DenominatorKinematicVariableIterator end_denominator_kinematic_variables() const;

            virtual SignalPDFPtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const;

            virtual std::ostream & insert(std::ostream & os) const;
    };

    template <typename... NumeratorKinematicNames_, typename... NormalizationKinematicNames_>
    SignalPDFEntry *
    make_concrete_signal_pdf_entry(const QualifiedName & name, const std::string & description, const Options & default_options, const QualifiedName & numerator,
                                   const std::tuple<NumeratorKinematicNames_...> & _numerator_kinematic_names, const QualifiedName & normalization,
                                   const std::tuple<NormalizationKinematicNames_...> & _normalization_kinematic_names)
    {
        std::vector<std::string> __numerator_kinematic_names{ impl::to_vector(_numerator_kinematic_names) };
        std::vector<std::string> __normalization_kinematic_names{ impl::to_vector(_normalization_kinematic_names) };

        auto result = new ConcreteSignalPDFEntry(name, description, default_options, numerator, normalization, __numerator_kinematic_names, __normalization_kinematic_names);

        return result;
    }
} // namespace eos

#endif
