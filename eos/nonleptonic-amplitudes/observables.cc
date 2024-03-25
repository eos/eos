/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

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

#include <eos/observable-impl.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes-adapter.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    /* nonleptonic amplitude as observables */
    template <typename Transition_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_nonleptonic_amplitudes_adapter(const char * name,
            const char * latex,
            double (NonleptonicAmplitudes<Transition_>::* _function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        QualifiedName qn(name);
        qnp::Prefix pp = qn.prefix_part();
        std::function<double (const NonleptonicAmplitudes<Transition_> *, const Args_ & ...)> function(_function);

        auto result = std::make_pair(qn, std::make_shared<NonleptonicAmplitudesAdapterEntry<Transition_, Args_ ...>>(qn, latex, Unit::None(), pp, function, kinematics_names));

        impl::observable_entries.insert(result);

        return result;
    }

    template <typename Transition_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_nonleptonic_amplitudes_adapter(const char * name,
            double (NonleptonicAmplitudes<Transition_>::* _function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        QualifiedName qn(name);
        qnp::Prefix pp = qn.prefix_part();
        std::function<double (const NonleptonicAmplitudes<Transition_> *, const Args_ & ...)> function(_function);

        auto result = std::make_pair(qn, std::make_shared<NonleptonicAmplitudesAdapterEntry<Transition_, Args_ ...>>(qn, "", Unit::None(), pp, function, kinematics_names));

        impl::observable_entries.insert(result);

        return result;
    }

    // Pseudo-observables related to the P(seudoscalar) -> P(seudoscalar) P(seudoscalar) amplitudes
    // {{{
    ObservableGroup
    make_p_to_p_p_amplitudes_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Pseudo-observables related to the $P\to PP$ amplitudes)",
            R"()",
            {
                make_nonleptonic_amplitudes_adapter("B^0->pi^+pi^-::Re{amplitude}", R"(\mathrm{Re}\,\mathcal{A}^{B^0\to\pi^+\pi^-})",
                        &NonleptonicAmplitudes<PToPP>::re_amplitude, std::make_tuple()),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_nonleptonic_amplitudes_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Pseudo-observables in nonleptonic amplitudes",
            "",
            {
                // P -> PP amplitudes
                make_p_to_p_p_amplitudes_group(),
            }
        );

        return ObservableSection(imp);
    }
}
