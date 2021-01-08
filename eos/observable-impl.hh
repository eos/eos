/* vim: set sw=4 sts=4 et tw=150 foldmethod=syntax : */

/*
 * Copyright (c) 2019 Danny van Dyk
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
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>

#include <array>
#include <map>

namespace eos
{
    template <>
    struct Implementation<ObservableGroup>
    {
        std::string name;

        std::string description;

        std::map<QualifiedName, ObservableEntryPtr> entries;

        Implementation(const std::string & name, const std::string & description,
                std::initializer_list<std::pair<const QualifiedName, ObservableEntryPtr>> && entries) :
            name(name),
            description(description),
            entries(entries)
        {
        }
    };

    template <>
    struct Implementation<ObservableSection>
    {
        std::string name;

        std::string description;

        std::vector<ObservableGroup> groups;

        Implementation(const std::string & name, const std::string & description,
                std::initializer_list<ObservableGroup> && groups) :
            name(name),
            description(description),
            groups(groups)
        {
        }
    };

    /* Helper functions to create ObservableEntry */
    template <typename Decay_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_observable(const char * name,
            double (Decay_::* function)(const Args_ & ...) const,
            const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        return std::make_pair(qn, make_concrete_observable_entry(qn, "", function, std::make_tuple(), forced_options));
    }

    template <typename Decay_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_observable(const char * name,
            const char * latex,
            double (Decay_::* function)(const Args_ & ...) const,
            const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        return std::make_pair(qn, make_concrete_observable_entry(qn, latex, function, std::make_tuple(), forced_options));
    }

    template <typename Decay_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_observable(const char * name,
            double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names,
            const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        return std::make_pair(qn, make_concrete_observable_entry(qn, "", function, kinematics_names, forced_options));
    }

    template <typename Decay_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_observable(const char * name,
            const char * latex,
            double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names,
            const Options & forced_options = Options{})
    {
        QualifiedName qn(name);

        return std::make_pair(qn, make_concrete_observable_entry(qn, latex, function, kinematics_names, forced_options));
    }

    /* ratios of regular observables */
    template <typename Decay_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_observable_ratio(const char * name,
            const char * latex,
            double (Decay_::* numerator)(const Args_ & ...) const,
            const Options & forced_options_numerator,
            double (Decay_::* denominator)(const Args_ & ...) const,
            const Options & forced_options_denominator
            )
    {
        QualifiedName qn(name);

        return std::make_pair(qn,
                make_concrete_observable_ratio_entry(
                        qn,
                        latex,
                        numerator,   std::make_tuple(), forced_options_numerator,
                        denominator, std::make_tuple(), forced_options_denominator
                       )
                );
    }

    template <typename Decay_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_observable_ratio(const char * name,
            const char * latex,
            double (Decay_::* numerator)(const Args_ & ...) const,
            const Tuple_ & kinematics_names_numerator,
            const Options & forced_options_numerator,
            double (Decay_::* denominator)(const Args_ & ...) const,
            const Tuple_ & kinematics_names_denominator,
            const Options & forced_options_denominator
            )
    {
        QualifiedName qn(name);

        return std::make_pair(qn,
                make_concrete_observable_ratio_entry(
                        qn,
                        latex,
                        numerator,   kinematics_names_numerator,   forced_options_numerator,
                        denominator, kinematics_names_denominator, forced_options_denominator
                        )
                );
    }

    /* sums of regular observables */

    template <typename Decay_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_observable_sum(const char * name,
            const char * latex,
            double (Decay_::* numerator)(const Args_ & ...) const,
            const Tuple_ & kinematics_names_numerator,
            const Options & forced_options_numerator,
            const double & weight_numerator,
            double (Decay_::* denominator)(const Args_ & ...) const,
            const Tuple_ & kinematics_names_denominator,
            const Options & forced_options_denominator,
            const double & weight_denominator
            )
    {
        QualifiedName qn(name);

        return std::make_pair(qn,
                make_concrete_observable_sum_entry(
                        qn,
                        latex,
                        numerator,   kinematics_names_numerator,   forced_options_numerator, weight_numerator,
                        denominator, kinematics_names_denominator, forced_options_denominator, weight_denominator
                        )
                );
    }

    template <>
    struct WrappedForwardIteratorTraits<ObservableEntry::KinematicVariableIteratorTag>
    {
        using UnderlyingIterator = std::array<const std::string, 1u>::iterator;
    };
}

#endif
