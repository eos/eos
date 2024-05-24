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
#include <eos/s-decays/k-to-pi-ll.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Semileptonic K -> pi decays
    // {{{

    // K -> pi l nu
    // {{{
    ObservableGroup
    make_k_to_pi_ll_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $K \to pi \ell^+\ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor.)",
            {
                // B -> K ll, Large Recoil
                make_observable("K->pill::d^2Gamma", R"(d^2\mathcal{\Gamma(\bar{K} \to \pi\ell^+\ell^-)}/(dq^2\, d\cos\theta_\ell))",
                        Unit::InverseGeV2(),
                        &KToPiDilepton::two_differential_decay_width,
                        std::make_tuple("q2", "cos(theta_l)")),

                make_observable("K->pill::dBR/dq2", R"(d\mathcal{B}(\bar{K} \to \pi\ell^+\ell^-)/dq^2)",
                        Unit::InverseGeV2(),
                        &KToPiDilepton::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("K->pill::F_H(q2)", R"(F_\mathrm{H}(\bar{K} \to \pi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &KToPiDilepton::differential_flat_term,
                        std::make_tuple("q2")),

                make_observable("K->pill::A_FB(q2)", R"(A_\mathrm{FB}(\bar{K} \to \pi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &KToPiDilepton::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_expression_observable("K->pill::R_K(q2)", R"(R_K(q^2))",
                        Unit::None(),
                        R"(
                        <<K->pill::dBR/dq2;l=mu>>
                        /
                        <<K->pill::dBR/dq2;l=e>>
                        )"),

                make_observable("K->pill::BR_CP_specific", R"(\mathcal{B}(\bar{K} \to \pi\ell^+\ell^-))",
                        Unit::None(),
                        &KToPiDilepton::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("K->pill::BR", R"(\bar{\mathcal{B}}(\bar{K} \to \pi\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<K->pill::BR_CP_specific;cp-conjugate=false>>
                               +
                               <<K->pill::BR_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_expression_observable("K->pill::A_CP", R"(A_\mathrm{CP}(\bar{K} \to \pi\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        (<<K->pill::BR_CP_specific;cp-conjugate=false>> - <<K->pill::BR_CP_specific;cp-conjugate=true>>)
                        /
                        (<<K->pill::BR_CP_specific;cp-conjugate=false>> + <<K->pill::BR_CP_specific;cp-conjugate=true>>)
                        )"),

                make_observable("K->pill::Gamma", R"(\Gamma(\bar{K} \to \pi\ell^+\ell^-))",
                        Unit::GeV(),
                        &KToPiDilepton::integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("K->pill::F_H_CP_specific", R"(F_\mathrm{H}(\bar{K} \to \pi\ell^+\ell^-))",
                        Unit::None(),
                        &KToPiDilepton::integrated_flat_term,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("K->pill::F_H", R"(\bar F_\mathrm{H}(\bar{K} \to \pi\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<K->pill::F_H_CP_specific;cp-conjugate=false>>
                               +
                               <<K->pill::F_H_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_observable("K->pill::A_FB_CP_specific", R"(A_\mathrm{FB}(\bar{K} \to \pi\ell^+\ell^-))",
                        Unit::None(),
                        &KToPiDilepton::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("K->pill::A_FB", R"(\bar A_\mathrm{FB}(\bar{K} \to \pi\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<K->pill::A_FB_CP_specific;cp-conjugate=false>>
                               +
                               <<K->pill::A_FB_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_expression_observable("K->pill::R_pi", R"(R_\pi)",
                        Unit::None(),
                        R"(
                        <<K->pill::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        /
                        <<K->pill::BR;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                        )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}
    // }}}

    ObservableSection
    make_s_decays_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in (semi)leptonic $s$-hadron decays",
            "",
            {
                // K -> pi l^+ l^-
                make_k_to_pi_ll_group()
            }
        );

        return ObservableSection(imp);
    }
}
