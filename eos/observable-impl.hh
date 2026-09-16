/* vim: set sw=4 sts=4 et tw=150 foldmethod=syntax : */

/*
 * Copyright (c) 2019-2026 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_OBSERVABLE_IMPL_HH
#define EOS_GUARD_EOS_OBSERVABLE_IMPL_HH 1

#include <eos/observable.hh>
#include <eos/utils/expression-observable.hh>
#include <eos/utils/expression-parser.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/tuple-maker.hh>
#include <eos/utils/units.hh>

#include <array>
#include <map>
#include <tuple>

namespace eos
{
    namespace impl
    {
        extern std::map<QualifiedName, ObservableEntryPtr> observable_entries;
    }

    template <> struct Implementation<ObservableGroup>
    {
            std::string name;

            std::string description;

            std::map<QualifiedName, ObservableEntryPtr> entries;

            Implementation(const std::string & name, const std::string & description, std::initializer_list<std::pair<const QualifiedName, ObservableEntryPtr>> && entries) :
                name(name),
                description(description),
                entries(entries)
            {
            }
    };

    template <> struct Implementation<ObservableSection>
    {
            std::string name;

            std::string description;

            std::vector<ObservableGroup> groups;

            Implementation(const std::string & name, const std::string & description, std::initializer_list<ObservableGroup> && groups) :
                name(name),
                description(description),
                groups(groups)
            {
            }
    };

    /* Helper functions to create ObservableEntry for a regular observable */
    template <typename Decay_, typename... Args_>
    std::pair<QualifiedName, ObservableEntryPtr>
    make_observable(const char * name, const Unit & unit, double (Decay_::*function)(const Args_ &...) const, const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        auto result = std::make_pair(qn, make_concrete_observable_entry(qn, "", unit, function, std::make_tuple(), forced_options));

        impl::observable_entries.insert(result);

        return result;
    }

    template <typename Decay_, typename... Args_>
    std::pair<QualifiedName, ObservableEntryPtr>
    make_observable(const char * name, const char * latex, const Unit & unit, double (Decay_::*function)(const Args_ &...) const, const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        auto result = std::make_pair(qn, make_concrete_observable_entry(qn, latex, unit, function, std::make_tuple(), forced_options));

        impl::observable_entries.insert(result);

        return result;
    }

    template <typename Decay_, typename Tuple_, typename... Args_>
    std::pair<QualifiedName, ObservableEntryPtr>
    make_observable(const char * name, const Unit & unit, double (Decay_::*function)(const Args_ &...) const, const Tuple_ & kinematics_names,
                    const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        auto result = std::make_pair(qn, make_concrete_observable_entry(qn, "", unit, function, kinematics_names, forced_options));

        impl::observable_entries.insert(result);

        return result;
    }

    template <typename Decay_, typename Tuple_, typename... Args_>
    std::pair<QualifiedName, ObservableEntryPtr>
    make_observable(const char * name, const char * latex, const Unit & unit, double (Decay_::*function)(const Args_ &...) const, const Tuple_ & kinematics_names,
                    const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        auto result = std::make_pair(qn, make_concrete_observable_entry(qn, latex, unit, function, kinematics_names, forced_options));

        impl::observable_entries.insert(result);

        return result;
    }

    namespace impl
    {
        /* A provider's prepare function, together with the kinematic variables it consumes */
        template <typename Decay_, typename IntermediateResult_, typename... Args_> struct Preparer
        {
                const IntermediateResult_ * (Decay_::*function)(const Args_ &...) const;

                std::tuple<typename ConvertTo<Args_, const char *>::Type...> kinematics_names;
        };

        /* A provider's evaluation function, together with the kinematic variables it consumes */
        template <typename Decay_, typename IntermediateResult_, typename... Args_> struct Evaluator
        {
                double (Decay_::*function)(const IntermediateResult_ *, const Args_ &...) const;

                std::tuple<typename ConvertTo<Args_, const char *>::Type...> kinematics_names;
        };
    } // namespace impl

    /*!
     * Declare the kinematic variables that enter the intermediate result of a cacheable observable.
     */
    template <typename Decay_, typename IntermediateResult_, typename... Args_, typename... Names_>
    impl::Preparer<Decay_, IntermediateResult_, Args_...>
    cache(const IntermediateResult_ * (Decay_::*prepare_fn)(const Args_ &...) const, const Names_ &... names)
    {
        static_assert(sizeof...(Names_) == sizeof...(Args_), "Need as many kinematic variable names as arguments of the prepare function!");

        return { prepare_fn, std::make_tuple(names...) };
    }

    /*!
     * Declare the kinematic variables that enter a cacheable observable beyond its intermediate result.
     */
    template <typename Decay_, typename IntermediateResult_, typename... Args_, typename... Names_>
    impl::Evaluator<Decay_, IntermediateResult_, Args_...>
    evaluate(double (Decay_::*evaluate_fn)(const IntermediateResult_ *, const Args_ &...) const, const Names_ &... names)
    {
        static_assert(sizeof...(Names_) == sizeof...(Args_), "Need as many kinematic variable names as arguments of the evaluation function!");

        return { evaluate_fn, std::make_tuple(names...) };
    }

    /* Helper functions to create ObservableEntry for a cacheable observable */
    template <typename Decay_, typename IntermediateResult_, typename... PrepareArgs_, typename... EvaluateArgs_>
    std::pair<QualifiedName, ObservableEntryPtr>
    make_cacheable_observable(const char * name, const Unit & unit, const impl::Preparer<Decay_, IntermediateResult_, PrepareArgs_...> & preparer,
                              const impl::Evaluator<Decay_, IntermediateResult_, EvaluateArgs_...> & evaluator, const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        auto result = std::make_pair(qn, make_concrete_cacheable_observable_entry(qn, "", unit, preparer, evaluator, forced_options));

        impl::observable_entries.insert(result);

        return result;
    }

    template <typename Decay_, typename IntermediateResult_, typename... PrepareArgs_, typename... EvaluateArgs_>
    std::pair<QualifiedName, ObservableEntryPtr>
    make_cacheable_observable(const char * name, const char * latex, const Unit & unit, const impl::Preparer<Decay_, IntermediateResult_, PrepareArgs_...> & preparer,
                              const impl::Evaluator<Decay_, IntermediateResult_, EvaluateArgs_...> & evaluator, const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        auto result = std::make_pair(qn, make_concrete_cacheable_observable_entry(qn, latex, unit, preparer, evaluator, forced_options));

        impl::observable_entries.insert(result);

        return result;
    }

    /* Deprecated: every kinematic variable enters the intermediate result */
    template <typename Decay_, typename IntermediateResult_, typename Tuple_, typename... Args_>
    std::pair<QualifiedName, ObservableEntryPtr>
    make_cacheable_observable(const char * name, const Unit & unit, const IntermediateResult_ * (Decay_::*prepare_fn)(const Args_ &...) const,
                              double (Decay_::*evaluate_fn)(const IntermediateResult_ *) const, const Tuple_ & kinematics_names, const Options & forced_options = Options{})
    {
        return std::apply([&](auto... names) { return make_cacheable_observable(name, unit, cache(prepare_fn, names...), evaluate(evaluate_fn), forced_options); },
                          kinematics_names);
    }

    template <typename Decay_, typename IntermediateResult_, typename Tuple_, typename... Args_>
    std::pair<QualifiedName, ObservableEntryPtr>
    make_cacheable_observable(const char * name, const char * latex, const Unit & unit, const IntermediateResult_ * (Decay_::*prepare_fn)(const Args_ &...) const,
                              double (Decay_::*evaluate_fn)(const IntermediateResult_ *) const, const Tuple_ & kinematics_names, const Options & forced_options = Options{})
    {
        return std::apply([&](auto... names) { return make_cacheable_observable(name, latex, unit, cache(prepare_fn, names...), evaluate(evaluate_fn), forced_options); },
                          kinematics_names);
    }

    /* expressions involving observables */

    std::pair<QualifiedName, ObservableEntryPtr> make_expression_observable(const char * name, const char * latex, const Unit & unit, const char * _expression);

    template <> struct WrappedForwardIteratorTraits<ObservableEntry::KinematicVariableIteratorTag>
    {
            using UnderlyingIterator = std::array<const std::string, 1u>::iterator;
    };

    template <> struct WrappedForwardIteratorTraits<ObservableEntry::OptionIteratorTag>
    {
            using UnderlyingIterator = std::vector<OptionSpecification>::const_iterator;
    };
} // namespace eos

#endif
