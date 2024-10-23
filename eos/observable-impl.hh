/* vim: set sw=4 sts=4 et tw=150 foldmethod=syntax : */

/*
 * Copyright (c) 2019, 2021, 2022 Danny van Dyk
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
#include <eos/utils/units.hh>

#include <array>
#include <map>

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

    /* Helper functions to create ObservableEntry for a cacheable observable */
    template <typename Decay_, typename Tuple_, typename... Args_>
    std::pair<QualifiedName, ObservableEntryPtr>
    make_cacheable_observable(const char * name, const char * latex, const Unit & unit, const typename Decay_::IntermediateResult * (Decay_::*prepare_fn)(const Args_ &...) const,
                              double (Decay_::*evaluate_fn)(const typename Decay_::IntermediateResult *) const, const Tuple_ & kinematics_names,
                              const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        auto result = std::make_pair(qn, make_concrete_cacheable_observable_entry(qn, latex, unit, prepare_fn, evaluate_fn, kinematics_names, forced_options));

        impl::observable_entries.insert(result);

        return result;
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
