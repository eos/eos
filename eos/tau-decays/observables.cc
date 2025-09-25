/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2024 Danny van Dyk
 * Copyright (c) 2025 Matthew Kirk
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
#include <eos/tau-decays/tau-to-k-nu.hh>
#include <eos/tau-decays/tau-to-k-pi-nu.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Hadronic tau decays
    // {{{

    // tau^- -> K^- nu
    // {{{
    ObservableGroup
    make_tau_to_k_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
                R"(Observables in $\tau \to K^- \nu_\tau$ decays)",
                R"()",
                {
                    make_observable("tau->Knu::BR", R"(\mathcal{B}(\tau^- \to K^- \nu))", Unit::None(), &TauToKNeutrino::branching_ratio, std::make_tuple(), Options{}),
                });
        return ObservableGroup(imp);
    }

    // }}}

    // tau^- -> [K pi]^- nu
    // {{{
    ObservableGroup
    make_tau_to_k_pi_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(R"(Observables in $\tau \to K \pi^- \nu_\tau$ decays)",
                                                       R"()",
                                                       {
                                                           make_observable("tau->K_Spinu::dBR/dq2",
                                                                           R"(d\mathcal{B}(tau^- \to K_S \pi^- \nu_\tau)/dq^2)",
                                                                           Unit::InverseGeV2(),
                                                                           &TauToKPiNeutrino::differential_branching_ratio,
                                                                           std::make_tuple("q2"),
                                                                           Options{ { "K"_ok, "K_S" } }
                                                                           ),
                                                           make_observable("tau->K_Spinu::dGamma/dq2",
                                                                           R"(d\Gamma(tau^- \to K_S \pi^- \nu_\tau)/dq^2)",
                                                                           Unit::InverseGeV(),
                                                                           &TauToKPiNeutrino::differential_decay_width,
                                                                           std::make_tuple("q2"),
                                                                           Options{ { "K"_ok, "K_S" } }
                                                                           ),
                                                           make_observable("tau->K_Spinu::BR",
                                                                           R"(\mathcal{B}(tau^- \to K_S \pi^- \nu_\tau))",
                                                                           Unit::None(),
                                                                           &TauToKPiNeutrino::total_branching_ratio,
                                                                           std::make_tuple(),
                                                                           Options{ { "K"_ok, "K_S" } }
                                                                           ),
                                                           make_observable("tau->K_Spinu::P(q2_min,q2_max)",
                                                                           R"(P(\tau^-\to K_S \pi^- \nu_\tau))",
                                                                           Unit::None(),
                                                                           &TauToKPiNeutrino::integrated_pdf_q2,
                                                                           std::make_tuple("q2_min", "q2_max"),
                                                                           Options{ { "K"_ok, "K_S" } }
                                                                           ),
                                                           make_observable("tau->K_Spinu::P(q2)",
                                                                           R"(dP(\tau^- \to K_S \pi^- \nu_\tau)/dq^2)",
                                                                           Unit::InverseGeV2(),
                                                                           &TauToKPiNeutrino::differential_pdf_q2,
                                                                           std::make_tuple("q2"),
                                                                           Options{ { "K"_ok, "K_S" } }
                                                                           ),
                                                           make_observable("tau->K^-pinu::BR",
                                                                           R"(\mathcal{B}(tau^- \to K^- \pi^0 \nu_\tau))",
                                                                           Unit::None(),
                                                                           &TauToKPiNeutrino::total_branching_ratio,
                                                                           std::make_tuple(),
                                                                           Options{ { "K"_ok, "K_u" } }
                                                                           ),
        });
        return ObservableGroup(imp);
    }

    // }}}

    // }}}

    ObservableSection
    make_tau_decays_section()
    {
        auto imp = new Implementation<ObservableSection>("Observables in $\tau$ decays",
                                                         "",
                                                         { // tau^- -> K^- nu
                                                           make_tau_to_k_nu_group(),
                                                           // tau^- -> [K pi]^- nu
                                                           make_tau_to_k_pi_nu_group() });

        return ObservableSection(imp);
    }
} // namespace eos
