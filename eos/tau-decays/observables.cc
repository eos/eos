/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2024 Danny van Dyk
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
#include <eos/tau-decays/tau-to-k-pi-nu.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Hadronic tau decays
    // {{{

    // tau^- -> [K pi]^- nu
    // {{{
    ObservableGroup
    make_tau_to_k_pi_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(R"(Observables in $\tau \to K \pi^- \bar{\nu}_\ell$ decays)",
                                                       R"(The option "l" selects the neutrino flavor.)",
                                                       { make_observable("tau->K_Spinu::dBR/dq2",
                                                                         R"(d\mathcal{B}(tau^- \to K \pi^- \bar{\nu}_\ell)/dq^2)",
                                                                         Unit::None(),
                                                                         &TauToKPiNeutrino::differential_branching_ratio,
                                                                         std::make_tuple("q2")),
                                                         make_observable("tau->K_Spinu::dGamma/dq2",
                                                                         R"(d\Gamma(tau^- \to K \pi^- \bar{\nu}_\ell)/dq^2)",
                                                                         Unit::InverseGeV2(),
                                                                         &TauToKPiNeutrino::differential_decay_width,
                                                                         std::make_tuple("q2")) });

        return ObservableGroup(imp);
    }

    // }}}

    // }}}

    ObservableSection
    make_tau_decays_section()
    {
        auto imp = new Implementation<ObservableSection>("Observables in $\tau$ decays",
                                                         "",
                                                         { // tau^- -> [K pi]^- nu
                                                           make_tau_to_k_pi_nu_group() });

        return ObservableSection(imp);
    }
} // namespace eos
