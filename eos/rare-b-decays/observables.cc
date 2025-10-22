/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2019-2025 Danny van Dyk
 * Copyright (c) 2021      Méril Reboud
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
#include <eos/rare-b-decays/decays.hh>
#include <eos/rare-b-decays/b-to-ll.hh>
#include <eos/rare-b-decays/b-to-k-charmonium.hh>
#include <eos/rare-b-decays/b-to-k-ll.hh>
#include <eos/rare-b-decays/b-to-kstar-charmonium.hh>
#include <eos/rare-b-decays/b-to-kstar-gamma.hh>
#include <eos/rare-b-decays/b-to-kstar-ll.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-impl.hh>
#include <eos/rare-b-decays/b-to-psd-nu-nu.hh>
#include <eos/rare-b-decays/b-to-vec-nu-nu.hh>
#include <eos/rare-b-decays/bs-to-phi-charmonium.hh>
#include <eos/rare-b-decays/bs-to-phi-ll.hh>
#include <eos/rare-b-decays/inclusive-b-to-s-dilepton.hh>
#include <eos/rare-b-decays/inclusive-b-to-s-gamma.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-dilepton.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-nu-nu.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-nu-nu-impl.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda1520-ll.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda1520-gamma.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // B_q -> l^+ l^-
    // {{{
    ObservableGroup
    make_b_to_ll_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q \to \ell^+\ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor.)",
            {
                make_observable("B_q->ll::BR", R"(\mathcal{B}(B_q \to \ell^+\ell^-))",
                        Unit::None(),
                        &BToDilepton::branching_ratio_time_zero),

                make_observable("B_q->ll::BR@Untagged", R"(\left\langle\mathcal{B}(B_q \to \ell^+\ell^-)\right\rangle)",
                        Unit::None(),
                        &BToDilepton::branching_ratio_untagged_integrated),

                make_observable("B_q->ll::A_DeltaGamma", R"(\mathcal{A}_{\Delta\Gamma}(B_q \to \ell^+\ell^-))",
                        Unit::None(),
                        &BToDilepton::cp_asymmetry_del_gamma),

                make_observable("B_q->ll::S", R"(\mathcal{S}(B_q \to \ell^+\ell^-))",
                        Unit::None(),
                        &BToDilepton::cp_asymmetry_mixing_S),

                make_observable("B_q->ll::eff_lifetime", R"(\langle\tau\rangle(B_q \to \ell^+\ell^-))",
                        Unit::None(),
                        &BToDilepton::effective_lifetime),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}


    // B_q -> P charmonium
    // {{{
    ObservableGroup
    make_b_to_p_charmonium_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q \to P \psi$ decays)",
            R"(The option "q" selects the spectator quark flavor.)",
            {
                /// Branching ratio of B -> K psi
                make_observable("B->Kpsi::BR", R"(\mathcal{B}(\bar{B} \to \bar{K}\psi))",
                        Unit::None(),
                        &BToKCharmonium::branching_ratio),
                make_observable("B->Kpsi::plus_phase",
                        Unit::None(),
                        &BToKCharmonium::plus_phase)
            }

        );

        return ObservableGroup(imp);
    }
    // }}}


    // B_q -> V charmonium
    // {{{
    ObservableGroup
    make_b_to_v_charmonium_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q \to V \psi$ decays)",
            R"(The option "q" selects the spectator quark flavor.)",
            {
                // Angular observables as detected in the decay B -> K^* psi (-> l^+ l^-)
                make_observable("B->K^*psi::perp_polarization",
                        Unit::None(),
                        &BToKstarCharmonium::perp_polarization),
                make_observable("B->K^*psi::para_polarization",
                        Unit::None(),
                        &BToKstarCharmonium::para_polarization),
                make_observable("B->K^*psi::long_polarization",
                        Unit::None(),
                        &BToKstarCharmonium::long_polarization),
                make_observable("B->K^*psi::long_phase",
                        Unit::None(),
                        &BToKstarCharmonium::long_phase),
                make_observable("B->K^*psi::delta_perp_long",
                        Unit::None(),
                        &BToKstarCharmonium::delta_perp_long),
                make_observable("B->K^*psi::delta_para_long",
                        Unit::None(),
                        &BToKstarCharmonium::delta_para_long),

                // Angular observables as detected in the decay B -> K^* psi (-> l^+ l^-)
                make_observable("B->K^*psi::S_1s@LHCb", R"(S_{1s}(\bar{B} \to \bar{K}^*\psi))",
                        Unit::None(),
                        &BToKstarCharmonium::S_1s_LHCb),
                make_observable("B->K^*psi::S_1c@LHCb", R"(S_{1c}(\bar{B} \to \bar{K}^*\psi))",
                        Unit::None(),
                        &BToKstarCharmonium::S_1c_LHCb),
                make_observable("B->K^*psi::S_3@LHCb", R"(S_{3}(\bar{B} \to \bar{K}^*\psi))",
                        Unit::None(),
                        &BToKstarCharmonium::S_3_LHCb),
                make_observable("B->K^*psi::S_4@LHCb", R"(S_{4}(\bar{B} \to \bar{K}^*\psi))",
                        Unit::None(),
                        &BToKstarCharmonium::S_4_LHCb),
                make_observable("B->K^*psi::S_8@LHCb", R"(S_{8}(\bar{B} \to \bar{K}^*\psi))",
                        Unit::None(),
                        &BToKstarCharmonium::S_8_LHCb),
                make_observable("B->K^*psi::S_9@LHCb", R"(S_{9}(\bar{B} \to \bar{K}^*\psi))",
                        Unit::None(),
                        &BToKstarCharmonium::S_9_LHCb),

                // Branching ratio of B -> K^* psi
                make_observable("B->K^*psi::BR", R"(\mathcal{B}(\bar{B} \to \bar{K}^*\psi))",
                        Unit::None(),
                        &BToKstarCharmonium::branching_ratio),


                // Angular observables as detected in the decay B_s -> phi psi (-> l^+ l^-)
                make_observable("B_s->phipsi::perp_polarization",
                        Unit::None(),
                        &BsToPhiCharmonium::perp_polarization,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phipsi::para_polarization",
                        Unit::None(),
                        &BsToPhiCharmonium::para_polarization,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phipsi::long_polarization",
                        Unit::None(),
                        &BsToPhiCharmonium::long_polarization,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phipsi::long_phase",
                        Unit::None(),
                        &BsToPhiCharmonium::long_phase,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phipsi::delta_perp_long",
                        Unit::None(),
                        &BsToPhiCharmonium::delta_perp_long,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phipsi::delta_para_long",
                        Unit::None(),
                        &BsToPhiCharmonium::delta_para_long,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } }),

                // Branching ratio of B -> phi psi
                make_observable("B_s->phipsi::BR", R"(\mathcal{B}(\bar{B}_s \to \phi\psi))",
                        Unit::None(),
                        &BsToPhiCharmonium::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } })

            }

        );

        return ObservableGroup(imp);
    }
    // }}}


    // B_q -> V gamma
    // {{{
    ObservableGroup
    make_b_to_v_gamma_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q \to V \gamma$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor.)",
            {
                // B -> K^* gamma
                make_observable("B->K^*gamma::BR_CP_specific", R"(\mathcal{B}(\bar{B}\to \bar{K}^*\gamma))",
                        Unit::None(),
                        &BToKstarGamma::branching_ratio),

                make_expression_observable("B->K^*gamma::BR", R"(\bar{\mathcal{B}}(\bar{B}\to \bar{K}^*\gamma))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<B->K^*gamma::BR_CP_specific;cp-conjugate=false>>
                               +
                               <<B->K^*gamma::BR_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_expression_observable("B->K^*gamma::A_CP", R"(A_\mathrm{CP}(\bar{B}\to \bar{K}^*\gamma))",
                        Unit::None(),
                        R"(
                        (<<B->K^*gamma::BR_CP_specific;cp-conjugate=false>> - <<B->K^*gamma::BR_CP_specific;cp-conjugate=true>>)
                        /
                        (<<B->K^*gamma::BR_CP_specific;cp-conjugate=false>> + <<B->K^*gamma::BR_CP_specific;cp-conjugate=true>>)
                        )"),

                // {S,C}_K^*gamma are calculated for B as the first state, Bbar as the second.
                // This is the opposite order than in B->K^*ll.
                make_observable("B->K^*gamma::Gamma_CP_specific",
                        Unit::None(),
                        &BToKstarGamma::decay_rate),
                make_observable("B->K^*gamma::Re{q_over_p}",
                        Unit::None(),
                        &BToKstarGamma::real_q_over_p),
                make_observable("B->K^*gamma::Im{q_over_p}",
                        Unit::None(),
                        &BToKstarGamma::imag_q_over_p),
                make_observable("B->K^*gamma::Re{a_left}",
                        Unit::None(),
                        &BToKstarGamma::real_a_left),
                make_observable("B->K^*gamma::Im{a_left}",
                        Unit::None(),
                        &BToKstarGamma::imag_a_left),
                make_observable("B->K^*gamma::Re{a_right}",
                        Unit::None(),
                        &BToKstarGamma::real_a_right),
                make_observable("B->K^*gamma::Im{a_right}",
                        Unit::None(),
                        &BToKstarGamma::imag_a_right),

                make_expression_observable("B->K^*gamma::S_K^*gamma", R"(S_{K^*\gamma})",
                        Unit::None(),
                        R"(
                        -2.0 * (
                            <<B->K^*gamma::Re{q_over_p}>> * (
                                    <<B->K^*gamma::Re{a_left};cp-conjugate=true>>  * <<B->K^*gamma::Im{a_right};cp-conjugate=false>>
                                  - <<B->K^*gamma::Im{a_left};cp-conjugate=true>>  * <<B->K^*gamma::Re{a_right};cp-conjugate=false>>
                                  + <<B->K^*gamma::Re{a_right};cp-conjugate=true>> * <<B->K^*gamma::Im{a_left};cp-conjugate=false>>
                                  - <<B->K^*gamma::Im{a_right};cp-conjugate=true>> * <<B->K^*gamma::Re{a_left};cp-conjugate=false>>
                            )
                            +
                            <<B->K^*gamma::Im{q_over_p}>> * (
                                    <<B->K^*gamma::Re{a_left};cp-conjugate=true>>  * <<B->K^*gamma::Re{a_right};cp-conjugate=false>>
                                  + <<B->K^*gamma::Re{a_right};cp-conjugate=true>> * <<B->K^*gamma::Re{a_left};cp-conjugate=false>>
                                  + <<B->K^*gamma::Im{a_left};cp-conjugate=true>>  * <<B->K^*gamma::Im{a_right};cp-conjugate=false>>
                                  + <<B->K^*gamma::Im{a_right};cp-conjugate=true>> * <<B->K^*gamma::Im{a_left};cp-conjugate=false>>
                            )
                        )
                        /
                        (<<B->K^*gamma::Gamma_CP_specific;cp-conjugate=false>> + <<B->K^*gamma::Gamma_CP_specific;cp-conjugate=true>>)
                        )"),

                make_expression_observable("B->K^*gamma::C_K^*gamma", R"(C_{K^*\gamma})",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*gamma::A_CP>> )"),

                make_expression_observable("B->K^*gamma::A_I", R"(A_\mathrm{I}(\bar{B}\to \bar{K}^*\gamma))",
                        Unit::None(),
                        R"(
                        (<<B->K^*gamma::BR_CP_specific;q=d>> - <<B->K^*gamma::BR_CP_specific;q=u>>)
                        /
                        (<<B->K^*gamma::BR_CP_specific;q=d>> + <<B->K^*gamma::BR_CP_specific;q=u>>)
                        )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}


    // B_q -> P l^+l^-
    // {{{
    ObservableGroup
    make_b_to_p_ll_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q \to P \ell^+\ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor.)",
            {
                // B -> K ll, Large Recoil
                make_observable("B->Kll::d^2Gamma", R"(d^2\mathcal{\Gamma(\bar{B}\to \bar{K}\ell^+\ell^-)}/(dq^2\, d\cos\theta_\ell))",
                        Unit::InverseGeV2(),
                        &BToKDilepton::two_differential_decay_width,
                        std::make_tuple("q2", "cos(theta_l)")),

                make_observable("B->Kll::dBR/ds", R"(d\mathcal{B}(\bar{B}\to \bar{K}\ell^+\ell^-)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToKDilepton::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B->Kll::F_H(q2)", R"(F_\mathrm{H}(\bar{B}\to \bar{K}\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKDilepton::differential_flat_term,
                        std::make_tuple("q2")),

                make_observable("B->Kll::A_FB(q2)", R"(A_\mathrm{FB}(\bar{B}\to \bar{K}\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKDilepton::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_expression_observable("B->Kll::R_K(q2)", R"(R_K(q^2))",
                        Unit::None(),
                        R"(
                        <<B->Kll::dBR/ds;l=mu>>
                        /
                        <<B->Kll::dBR/ds;l=e>>
                        )"),

                make_observable("B->Kll::BR_CP_specific", R"(\mathcal{B}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        Unit::None(),
                        &BToKDilepton::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("B->Kll::BR", R"(\bar{\mathcal{B}}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<B->Kll::BR_CP_specific;cp-conjugate=false>>
                               +
                               <<B->Kll::BR_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_expression_observable("B->Kll::NormalizedBR", R"(\mathcal{B}(\bar{B}\to \bar{K}\ell^+\ell^-)/\mathcal{B}(\bar{B}\to \bar{K}J/\psi))",
                        Unit::None(),
                        R"(
                        <<B->Kll::BR>> / <<B->Kpsi::BR;psi=J/psi>>
                        )"),

                make_expression_observable("B->Kll::A_CP", R"(A_\mathrm{CP}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        (<<B->Kll::BR_CP_specific;cp-conjugate=false>> - <<B->Kll::BR_CP_specific;cp-conjugate=true>>)
                        /
                        (<<B->Kll::BR_CP_specific;cp-conjugate=false>> + <<B->Kll::BR_CP_specific;cp-conjugate=true>>)
                        )"),

                make_observable("B->Kll::Gamma", R"(\Gamma(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        Unit::GeV(),
                        &BToKDilepton::integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::F_H_CP_specific", R"(F_\mathrm{H}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        Unit::None(),
                        &BToKDilepton::integrated_flat_term,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("B->Kll::F_H", R"(\bar F_\mathrm{H}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<B->Kll::F_H_CP_specific;cp-conjugate=false>>
                               +
                               <<B->Kll::F_H_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_observable("B->Kll::A_FB_CP_specific", R"(A_\mathrm{FB}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        Unit::None(),
                        &BToKDilepton::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("B->Kll::A_FB", R"(\bar A_\mathrm{FB}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<B->Kll::A_FB_CP_specific;cp-conjugate=false>>
                               +
                               <<B->Kll::A_FB_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_expression_observable("B->Kll::R_K", R"(R_K)",
                        Unit::None(),
                        R"(
                        <<B->Kll::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        /
                        <<B->Kll::BR;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                        )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}


    // B_q -> V l^+l^-
    // {{{
    ObservableGroup
    make_b_to_v_ll_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q \to V \ell^+\ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor.)",
            {
                make_observable("B->K^*ll::d^4Gamma",
                        Unit::GeV(),
                        &BToKstarDilepton::decay_width,
                        std::make_tuple("q2", "cos(theta_l)", "cos(theta_k)", "phi")),

                make_observable("B->K^*ll::dBR/ds", R"(d\mathcal{B}/dq^2(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::InverseGeV2(),
                        &BToKstarDilepton::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_FB(q2)", R"(A_\mathrm{FB}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^2(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_transverse_asymmetry_2,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^3(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_transverse_asymmetry_3,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^4(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_transverse_asymmetry_4,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^5(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_transverse_asymmetry_5,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^re(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_transverse_asymmetry_re,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^im(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_transverse_asymmetry_im,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::F_L(q2)", R"(F_L(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_longitudinal_polarisation,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::F_T(q2)", R"(F_T(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_transversal_polarisation,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_1s(q2)", R"(J_{1s}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_1s,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_1c(q2)", R"(J_{1c}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_1c,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_2s(q2)", R"(J_{2s}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_2s,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_2c(q2)", R"(J_{2c}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_2c,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_3(q2)", R"(J_3(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_3,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_4(q2)", R"(J_4(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_4,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_5(q2)", R"(J_5(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_5,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_6s(q2)", R"(J_{6s}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_6s,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_6c(q2)", R"(J_{6c}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_6c,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_7(q2)", R"(J_7(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_7,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_8(q2)", R"(J_8(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_8,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_9(q2)", R"(J_9(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BToKstarDilepton::differential_j_9,
                        std::make_tuple("q2")),

                make_expression_observable("B->K^*ll::P'_4(q2)", R"(P'_4(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q2))",
                        Unit::None(),
                        R"(
                        (<<B->K^*ll::J_4(q2);cp-conjugate=false>> + <<B->K^*ll::J_4(q2);cp-conjugate=true>>)
                        /
                        ( -1.0 *
                          (<<B->K^*ll::J_2c(q2);cp-conjugate=false>> + <<B->K^*ll::J_2c(q2);cp-conjugate=true>>) *
                          (<<B->K^*ll::J_2s(q2);cp-conjugate=false>> + <<B->K^*ll::J_2s(q2);cp-conjugate=true>>)
                         ) ^ 0.5
                        )"),

                make_expression_observable("B->K^*ll::P'_5(q2)", R"(P'_5(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q2))",
                        Unit::None(),
                        R"(
                        0.5 * (<<B->K^*ll::J_5(q2);cp-conjugate=false>> + <<B->K^*ll::J_5(q2);cp-conjugate=true>>)
                        /
                        ( -1.0 *
                          (<<B->K^*ll::J_2c(q2);cp-conjugate=false>> + <<B->K^*ll::J_2c(q2);cp-conjugate=true>>) *
                          (<<B->K^*ll::J_2s(q2);cp-conjugate=false>> + <<B->K^*ll::J_2s(q2);cp-conjugate=true>>)
                         ) ^ 0.5
                        )"),

                make_expression_observable("B->K^*ll::P'_6(q2)", R"(P'_6(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q2))",
                        Unit::None(),
                        R"(
                        -0.5 * (<<B->K^*ll::J_7(q2);cp-conjugate=false>> + <<B->K^*ll::J_7(q2);cp-conjugate=true>>)
                        /
                        ( -1.0 *
                          (<<B->K^*ll::J_2c(q2);cp-conjugate=false>> + <<B->K^*ll::J_2c(q2);cp-conjugate=true>>) *
                          (<<B->K^*ll::J_2s(q2);cp-conjugate=false>> + <<B->K^*ll::J_2s(q2);cp-conjugate=true>>)
                         ) ^ 0.5
                        )"),

                make_expression_observable("B->K^*ll::R_K^*(q2)", R"(R_{K^*}(q^2))",
                        Unit::None(),
                        R"(
                        <<B->K^*ll::dBR/ds;l=mu>>
                        /
                        <<B->K^*ll::dBR/ds;l=e>>
                        )"),

                make_cacheable_observable("B->K^*ll::A_FB_CP_specific", R"(A_\mathrm{FB}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::Abar_FB", R"()",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_unnormalized_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("B->K^*ll::A_FB", R"(\bar{A}_\mathrm{FB}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<B->K^*ll::A_FB_CP_specific;cp-conjugate=false>>
                               +
                               <<B->K^*ll::A_FB_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_cacheable_observable("B->K^*ll::BR_CP_specific", R"(\mathcal{B}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("B->K^*ll::BR", R"(\bar{\mathcal{B}}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<B->K^*ll::BR_CP_specific;cp-conjugate=false>>
                               +
                               <<B->K^*ll::BR_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_expression_observable("B->K^*ll::A_CP", R"(\bar{A}_\mathrm{CP}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        (<<B->K^*ll::BR_CP_specific;cp-conjugate=false>> - <<B->K^*ll::BR_CP_specific;cp-conjugate=true>>)
                        /
                        (<<B->K^*ll::BR_CP_specific;cp-conjugate=false>> + <<B->K^*ll::BR_CP_specific;cp-conjugate=true>>)
                        )"),

                make_cacheable_observable("B->K^*ll::F_L_CP_specific", R"(F_L(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_longitudinal_polarisation,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("B->K^*ll::F_L", R"(\bar{F}_L(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<B->K^*ll::F_L_CP_specific;cp-conjugate=false>>
                               +
                               <<B->K^*ll::F_L_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_cacheable_observable("B->K^*ll::F_T_CP_specific", R"(F_T(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_transversal_polarisation,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("B->K^*ll::F_T", R"(\bar{T}_L(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<B->K^*ll::F_T_CP_specific;cp-conjugate=false>>
                               +
                               <<B->K^*ll::F_T_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_cacheable_observable("B->K^*ll::A_T^2_CP_specific",  R"()",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_transverse_asymmetry_2,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("B->K^*ll::A_T^2", R"(\bar{A}_T^2(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<B->K^*ll::A_T^2_CP_specific;cp-conjugate=false>>
                               +
                               <<B->K^*ll::A_T^2_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_cacheable_observable("B->K^*ll::A_T^3", R"(A_T^3(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_transverse_asymmetry_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::A_T^4", R"(A_T^4(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_transverse_asymmetry_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::A_T^5", R"(A_T^5(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_transverse_asymmetry_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::A_T^re", R"(\mathrm{Re}A_T(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_transverse_asymmetry_re,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::A_T^im", R"(\mathrm{Im}A_T(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_transverse_asymmetry_im,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::H_T^1(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_h_1,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::H_T^2(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_h_2,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::H_T^3(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_h_3,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::H_T^4(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_h_4,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::H_T^5(q2)",
                        Unit::None(),
                        &BToKstarDilepton::differential_h_5,
                        std::make_tuple("q2")),

                make_cacheable_observable("B->K^*ll::H_T^1", R"(H_T^1(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_h_1,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::H_T^2", R"(H_T^2(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_h_2,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::H_T^3", R"(H_T^3(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_h_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::H_T^4", R"(H_T^4(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_h_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::H_T^5", R"(H_T^5(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_h_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::s_0^A_FB",
                        Unit::GeV2(),
                        &BToKstarDilepton::a_fb_zero_crossing),

                make_cacheable_observable("B->K^*ll::Gamma_CP_specific", R"()",
                        Unit::GeV(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::Gamma_CP_specific(q2)", R"()",
                        Unit::GeV(),
                        &BToKstarDilepton::differential_decay_width,
                        std::make_tuple("q2")),

                make_expression_observable("B->K^*ll::Gamma", R"(\Gamma(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (<<B->K^*ll::Gamma_CP_specific;cp-conjugate=false>> + <<B->K^*ll::Gamma_CP_specific;cp-conjugate=true>>)
                        )"),

                make_expression_observable("B->K^*ll::Gamma(q2)", R"(\Gamma^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        0.5 * (<<B->K^*ll::Gamma_CP_specific(q2);cp-conjugate=false>> + <<B->K^*ll::Gamma_CP_specific(q2);cp-conjugate=true>>)
                        )"),

                make_cacheable_observable("B->K^*ll::J_1s", R"(J_{1s}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_1s,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::J_1c", R"(J_{1c}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_1c,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::J_2s", R"(J_{2s}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_2s,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::J_2c", R"(J_{2c}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_2c,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::J_3", R"(J_3(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("B->K^*ll::J_3norm_CP_specific", R"(J_3/\Gamma(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        <<B->K^*ll::J_3>> / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::J_3norm",  R"(\bar{J}/\bar{\Gamma}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        (<<B->K^*ll::J_3;cp-conjugate=false>> + <<B->K^*ll::J_3;cp-conjugate=true>>)
                        /
                        (<<B->K^*ll::Gamma;cp-conjugate=false>> + <<B->K^*ll::Gamma;cp-conjugate=true>>)
                        )"),

                make_cacheable_observable("B->K^*ll::J_4", R"(J_4(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::J_5", R"(J_5(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::J_6s", R"(J_{6s}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_6s,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::J_6c", R"(J_{6c}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_6c,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::J_7", R"(J_7(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_7,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::J_8", R"(J_8(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_8,
                        std::make_tuple("q2_min", "q2_max")),

                make_cacheable_observable("B->K^*ll::J_9", R"(J_9(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        &BToKstarDilepton::prepare,
                        &BToKstarDilepton::integrated_j_9,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("B->K^*ll::J_9norm_CP_specific", R"(J_9/\Gamma(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        <<B->K^*ll::J_9>> / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::J_9norm", R"(\bar{J}/\bar{\Gamma}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        (<<B->K^*ll::J_9;cp-conjugate=false>> + <<B->K^*ll::J_9;cp-conjugate=true>>)
                        /
                        (<<B->K^*ll::Gamma;cp-conjugate=false>> + <<B->K^*ll::Gamma;cp-conjugate=true>>)
                        )"),

                make_expression_observable("B->K^*ll::S_1s(q2)", R"(S_{1s}^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_1s(q2);cp-conjugate=false>> + <<B->K^*ll::J_1s(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_1s", R"(S_{1s}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_1s;cp-conjugate=false>> + <<B->K^*ll::J_1s;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_1c(q2)", R"(S_{1c}^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_1c(q2);cp-conjugate=false>> + <<B->K^*ll::J_1c(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_1c", R"(S_{1c}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_1c;cp-conjugate=false>> + <<B->K^*ll::J_1c;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_2s(q2)", R"(S_{2s}^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_2s(q2);cp-conjugate=false>> + <<B->K^*ll::J_2s(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_2s", R"(S_{2s}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_2s;cp-conjugate=false>> + <<B->K^*ll::J_2s;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_2c(q2)", R"(S_{2c}^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_2c(q2);cp-conjugate=false>> + <<B->K^*ll::J_2c(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_2c", R"(S_{2c}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_2c;cp-conjugate=false>> + <<B->K^*ll::J_2c;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_3(q2)", R"(S_3^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_3(q2);cp-conjugate=false>> + <<B->K^*ll::J_3(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_3", R"(S_3(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_3;cp-conjugate=false>> + <<B->K^*ll::J_3;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_4(q2)", R"(S_4^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_4(q2);cp-conjugate=false>> + <<B->K^*ll::J_4(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_4", R"(S_4(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_4;cp-conjugate=false>> + <<B->K^*ll::J_4;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_5(q2)", R"(S_5^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_5(q2);cp-conjugate=false>> + <<B->K^*ll::J_5(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_5", R"(S_5(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_5;cp-conjugate=false>> + <<B->K^*ll::J_5;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_6s(q2)", R"(S_{6s}^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_6s(q2);cp-conjugate=false>> + <<B->K^*ll::J_6s(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_6s", R"(S_{6s}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_6s;cp-conjugate=false>> + <<B->K^*ll::J_6s;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_6c(q2)", R"(S_{6c}^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_6c(q2);cp-conjugate=false>> + <<B->K^*ll::J_6c(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_6c", R"(S_{6c}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_6c;cp-conjugate=false>> + <<B->K^*ll::J_6c;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_7(q2)", R"(S_7^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_7(q2);cp-conjugate=false>> + <<B->K^*ll::J_7(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_7", R"(S_7(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_7;cp-conjugate=false>> + <<B->K^*ll::J_7;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_8(q2)", R"(S_8^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_8(q2);cp-conjugate=false>> + <<B->K^*ll::J_8(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_8", R"(S_8(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_8;cp-conjugate=false>> + <<B->K^*ll::J_8;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::S_9(q2)", R"(S_9^{\bar{B}\to \bar{K}^*\ell^+\ell^-}(q^2))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_9(q2);cp-conjugate=false>> + <<B->K^*ll::J_9(q2);cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma(q2)>>
                        )"),

                make_expression_observable("B->K^*ll::S_9", R"(S_9(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_9;cp-conjugate=false>> + <<B->K^*ll::J_9;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_1s", R"(A_{1s}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_1s;cp-conjugate=false>> - <<B->K^*ll::J_1s;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_1c", R"(A_{1c}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_1c;cp-conjugate=false>> - <<B->K^*ll::J_1c;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_2s", R"(A_{2s}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_2s;cp-conjugate=false>> - <<B->K^*ll::J_2s;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_2c", R"(A_{2c}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_2c;cp-conjugate=false>> - <<B->K^*ll::J_2c;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_3", R"(A_3(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_3;cp-conjugate=false>> - <<B->K^*ll::J_3;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_4", R"(A_4(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_4;cp-conjugate=false>> - <<B->K^*ll::J_4;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_5", R"(A_5(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_5;cp-conjugate=false>> - <<B->K^*ll::J_5;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_6s", R"(A_{6s}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_6s;cp-conjugate=false>> - <<B->K^*ll::J_6s;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_6c", R"(A_{6c}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_6c;cp-conjugate=false>> - <<B->K^*ll::J_6c;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_7", R"(A_7(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_7;cp-conjugate=false>> - <<B->K^*ll::J_7;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_8", R"(A_8(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_8;cp-conjugate=false>> - <<B->K^*ll::J_8;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::A_9", R"(A_9(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        2.0 / 3.0 * (<<B->K^*ll::J_9;cp-conjugate=false>> - <<B->K^*ll::J_9;cp-conjugate=true>>)
                                  / <<B->K^*ll::Gamma>>
                        )"),

                make_expression_observable("B->K^*ll::N'_bin", R"(\mathcal{N}'_\mathrm{bin}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        ( -1.0 *
                          (<<B->K^*ll::J_2c;cp-conjugate=false>> + <<B->K^*ll::J_2c;cp-conjugate=true>>) *
                          (<<B->K^*ll::J_2s;cp-conjugate=false>> + <<B->K^*ll::J_2s;cp-conjugate=true>>)
                         ) ^ 0.5
                        )"),

                make_expression_observable("B->K^*ll::P_1", R"(P_1(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * <<B->K^*ll::S_3>> / <<B->K^*ll::S_2s>>
                        )"),

                make_expression_observable("B->K^*ll::P_2", R"(P_2(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        1.0 / 8.0 * <<B->K^*ll::S_6s>> / <<B->K^*ll::S_2s>>
                        )"),

                make_expression_observable("B->K^*ll::P_3", R"(P_3(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        -0.25 * <<B->K^*ll::S_9>> / <<B->K^*ll::S_2s>>
                        )"),

                make_expression_observable("B->K^*ll::P'_4", R"(P'_4(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        (<<B->K^*ll::J_4;cp-conjugate=false>> + <<B->K^*ll::J_4;cp-conjugate=true>>) / <<B->K^*ll::N'_bin>>
                        )"),

                make_expression_observable("B->K^*ll::P'_5", R"(P'_5(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (<<B->K^*ll::J_5;cp-conjugate=false>> + <<B->K^*ll::J_5;cp-conjugate=true>>) / <<B->K^*ll::N'_bin>>
                        )"),

                make_expression_observable("B->K^*ll::P'_6", R"(P'_6(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        -0.5 * (<<B->K^*ll::J_7;cp-conjugate=false>> + <<B->K^*ll::J_7;cp-conjugate=true>>) / <<B->K^*ll::N'_bin>>
                        )"),

                make_expression_observable("B->K^*ll::P'_8", R"(P'_8(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        -1.0 * (<<B->K^*ll::J_8;cp-conjugate=false>> + <<B->K^*ll::J_8;cp-conjugate=true>>) / <<B->K^*ll::N'_bin>>
                        )"),


                // Observables in the LHCb angular convention: cf. DHMV:2015A p. 9
                make_expression_observable("B->K^*ll::S_1s(q2)@LHCb", R"(S_{1s}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_1s(q2)>> )"),

                make_expression_observable("B->K^*ll::S_1c(q2)@LHCb", R"(S_{1c}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_1c(q2)>> )"),

                make_expression_observable("B->K^*ll::S_2s(q2)@LHCb", R"(S_{2s}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_2s(q2)>> )"),

                make_expression_observable("B->K^*ll::S_2c(q2)@LHCb", R"(S_{2c}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_2c(q2)>> )"),

                make_expression_observable("B->K^*ll::S_3(q2)@LHCb", R"(S_3^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_3(q2)>> )"),

                make_expression_observable("B->K^*ll::S_4(q2)@LHCb", R"(S_4^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::S_4(q2)>> )"),

                make_expression_observable("B->K^*ll::S_5(q2)@LHCb", R"(S_5^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_5(q2)>> )"),

                make_expression_observable("B->K^*ll::S_6s(q2)@LHCb", R"(S_{6s}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::S_6s(q2)>> )"),

                make_expression_observable("B->K^*ll::S_6c(q2)@LHCb", R"(S_{6c}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::S_6c(q2)>> )"),

                make_expression_observable("B->K^*ll::S_7(q2)@LHCb", R"(S_7^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::S_7(q2)>> )"),

                make_expression_observable("B->K^*ll::S_8(q2)@LHCb", R"(S_8^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_8(q2)>> )"),

                make_expression_observable("B->K^*ll::S_9(q2)@LHCb", R"(S_9^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::S_9(q2)>> )"),

                make_expression_observable("B->K^*ll::A_FB(q2)@LHCb", R"(A_\mathrm{FB}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::A_FB(q2)>> )"),

                make_expression_observable("B->K^*ll::S_1s@LHCb", R"(S_{1s}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_1s>> )"),

                make_expression_observable("B->K^*ll::S_1c@LHCb", R"(S_{1c}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_1c>> )"),

                make_expression_observable("B->K^*ll::S_2s@LHCb", R"(S_{2s}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_2s>> )"),

                make_expression_observable("B->K^*ll::S_2c@LHCb", R"(S_{2c}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_2c>> )"),

                make_expression_observable("B->K^*ll::S_3@LHCb", R"(S_3^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_3>> )"),

                make_expression_observable("B->K^*ll::S_4@LHCb", R"(S_4^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::S_4>> )"),

                make_expression_observable("B->K^*ll::S_5@LHCb", R"(S_5^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_5>> )"),

                make_expression_observable("B->K^*ll::S_6s@LHCb", R"(S_{6s}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::S_6s>> )"),

                make_expression_observable("B->K^*ll::S_6c@LHCb", R"(S_{6c}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::S_6c>> )"),

                make_expression_observable("B->K^*ll::S_7@LHCb", R"(S_7^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::S_7>> )"),

                make_expression_observable("B->K^*ll::S_8@LHCb", R"(S_8^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( <<B->K^*ll::S_8>> )"),

                make_expression_observable("B->K^*ll::S_9@LHCb", R"(S_9^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::S_9>> )"),

                make_expression_observable("B->K^*ll::A_FB@LHCb", R"(A_\mathrm{FB}^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::A_FB>> )"),

                make_expression_observable("B->K^*ll::P_1@LHCb", R"(P_1^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( <<B->K^*ll::P_1>> )"),

                make_expression_observable("B->K^*ll::P_2@LHCb", R"(P_2^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::P_2>> )"),

                make_expression_observable("B->K^*ll::P_3@LHCb", R"(P_3^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( -1.0 * <<B->K^*ll::P_3>> )"),

                make_expression_observable("B->K^*ll::P'_4@LHCb", R"(P'_4^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( -0.5 * <<B->K^*ll::P'_4>> )"),

                make_expression_observable("B->K^*ll::P'_5@LHCb", R"(P'_5^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( <<B->K^*ll::P'_5>> )"),

                make_expression_observable("B->K^*ll::P'_6@LHCb", R"(P'_6^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( <<B->K^*ll::P'_6>> )"),

                make_expression_observable("B->K^*ll::P'_8@LHCb", R"(P'_8^\mathrm{LHCb}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"( -0.5 * <<B->K^*ll::P'_8>> )"),


                make_expression_observable("B->K^*ll::R_K^*", R"(R_{K^*})",
                        Unit::None(),
                        R"(
                        <<B->K^*ll::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        /
                        <<B->K^*ll::BR;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                        )"),

                make_expression_observable("B->K^*ll::NormalizedBR", R"(\mathcal{B}(\bar{B}\to \bar{K}^*\ell^+\ell^-)/\mathcal{B}(\bar{B}\to \bar{K}^*J/\psi))",
                        Unit::None(),
                        R"(
                        <<B->K^*ll::BR>> / <<B->K^*psi::BR;psi=J/psi>>
                        )"),

                make_observable("B->K^*ll::Re{C9_perp}(q2)",
                        Unit::None(),
                        &BToKstarDilepton::real_C9_perp,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::Re{C9_para}(q2)",
                        Unit::None(),
                        &BToKstarDilepton::real_C9_para,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::Im{C9_perp}(q2)",
                        Unit::None(),
                        &BToKstarDilepton::imag_C9_perp,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::Im{C9_para}(q2)",
                        Unit::None(),
                        &BToKstarDilepton::imag_C9_para,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::H_perp_corrections(q2)",
                        Unit::None(),
                        &BToKstarDilepton::H_perp_corrections,
                        std::make_tuple("q2")),
                make_observable("B->K^*ll::H_para_corrections(q2)",
                        Unit::None(),
                        &BToKstarDilepton::H_para_corrections,
                        std::make_tuple("q2")),
                make_observable("B->K^*ll::H_long_corrections(q2)",
                        Unit::None(),
                        &BToKstarDilepton::H_long_corrections,
                        std::make_tuple("q2")),






                // B_s^0 -> \phi \ell^+ \ell^-
                make_observable("B_s->phill::d^4Gamma",
                        Unit::InverseGeV2(),
                        &BsToPhiDilepton::decay_width,
                        std::make_tuple("q2", "cos(theta_l)", "cos(theta_k)", "phi"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::dBR/ds", R"(d\mathcal{B}/dq^2(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::InverseGeV2(),
                        &BsToPhiDilepton::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_FB(q2)", R"(A_\mathrm{FB}(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_forward_backward_asymmetry,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::F_L(q2)", R"(F_L(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_longitudinal_polarisation,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_expression_observable("B_s->phill::R_phi(q2)", R"(R_{\phi}(q^2))",
                        Unit::None(),
                        R"(
                        <<B_s->phill::dBR/ds;l=mu>>
                        /
                        <<B_s->phill::dBR/ds;l=e>>
                        )"),

                make_observable("B_s->phill::A_FB", R"(A_\mathrm{FB}(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::BR_CP_specific", R"(\mathcal{B}(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_expression_observable("B_s->phill::BR", R"(\bar{\mathcal{B}}(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (
                               <<B_s->phill::BR_CP_specific;cp-conjugate=false>>
                               +
                               <<B_s->phill::BR_CP_specific;cp-conjugate=true>>
                               )
                        )"),

                make_observable("B_s->phill::F_L", R"(F_L(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_longitudinal_polarisation,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::Gamma_CP_specific",
                        Unit::GeV(),
                        &BsToPhiDilepton::integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::Gamma_CP_specific(q2)",
                        Unit::GeV(),
                        &BsToPhiDilepton::differential_decay_width,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::Gamma", R"(\Gamma(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::GeV(),
                        &BsToPhiDileptonAndConjugate::integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::Gamma(q2)", R"(\Gamma(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::GeV(),
                        &BsToPhiDileptonAndConjugate::differential_decay_width,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_1s(q2)", R"(J_{1s}(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_1s,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_1c(q2)", R"(J_{1c}(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_1c,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_2s(q2)", R"(J_{2s}(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_2s,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_2c(q2)", R"(J_{2c}(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_2c,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_3(q2)", R"(J_3(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_3,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_4(q2)", R"(J_4(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_4,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_5(q2)", R"(J_5(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_5,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_6s(q2)", R"(J_{6s}(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_6s,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_6c(q2)", R"(J_{6c}(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_6c,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_7(q2)", R"(J_7(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_7,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_8(q2)", R"(J_8(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_8,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_9(q2)", R"(J_9(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &BsToPhiDilepton::differential_j_9,
                        std::make_tuple("q2")),

                make_observable("B_s->phill::J_1s", R"(J_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_1s,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_1c", R"(J_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_1c,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_2s", R"(J_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_2s,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_2c", R"(J_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_2c,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_3", R"(J_3(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_3,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_4", R"(J_4(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_4,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_5", R"(J_5(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_5,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_6s", R"(J_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_6s,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_6c", R"(J_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_6c,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_7", R"(J_7(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_7,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_8", R"(J_8(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_8,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::J_9", R"(J_9(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDilepton::integrated_j_9,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" } }),


                make_observable("B_s->phill::H_1s(q2)", R"(H_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_1s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_1s(q2)@LHCb", R"(H_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_1s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_1s", R"(H_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_1s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_1s@LHCb", R"(H_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::H_1s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_1c(q2)", R"(H_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_1c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_1c(q2)@LHCb", R"(H_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_1c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_1c", R"(H_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_1c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_1c@LHCb", R"(H_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::H_1c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_2s(q2)", R"(H_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_2s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_2s(q2)@LHCb", R"(H_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_2s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_2s", R"(H_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_2s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_2s@LHCb", R"(H_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::H_2s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_2c(q2)", R"(H_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_2c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_2c(q2)@LHCb", R"(H_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_2c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_2c", R"(H_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_2c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_2c@LHCb", R"(H_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::H_2c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_3(q2)", R"(H_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_3,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_3(q2)@LHCb", R"(H_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_3(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_3", R"(H_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_3,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_3@LHCb", R"(H_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::H_3>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_4(q2)", R"(H_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_4,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_4(q2)@LHCb", R"(H_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::H_4(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_4", R"(H_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_4,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_4@LHCb", R"(H_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::H_4>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_5(q2)", R"(H_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_5,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_5(q2)@LHCb", R"(H_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_5(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_5", R"(H_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_5,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_5@LHCb", R"(H_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::H_5>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_6s(q2)", R"(H_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_6s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_6s(q2)@LHCb", R"(H_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::H_6s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_6s", R"(H_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_6s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_6s@LHCb", R"(H_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::H_6s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_6c(q2)", R"(H_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_6c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_6c(q2)@LHCb", R"(H_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::H_6c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_6c", R"(H_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_6c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_6c@LHCb", R"(H_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::H_6c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_7(q2)", R"(H_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_7,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_7(q2)@LHCb", R"(H_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::H_7(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_7", R"(H_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_7,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_7@LHCb", R"(H_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::H_7>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_8(q2)", R"(H_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_8,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_8(q2)@LHCb", R"(H_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_8(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_8", R"(H_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_8,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_8@LHCb", R"(H_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::H_8>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::H_9(q2)", R"(H_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_H_9,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_9(q2)@LHCb", R"(H_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::H_9(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::H_9", R"(H_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_H_9,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::H_9@LHCb", R"(H_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::H_9>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_1s(q2)", R"(Z_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_1s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_1s(q2)@LHCb", R"(Z_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_1s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_1s", R"(Z_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_1s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_1s@LHCb", R"(Z_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::Z_1s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_1c(q2)", R"(Z_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_1c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_1c(q2)@LHCb", R"(Z_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_1c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_1c", R"(Z_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_1c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_1c@LHCb", R"(Z_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::Z_1c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_2s(q2)", R"(Z_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_2s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_2s(q2)@LHCb", R"(Z_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_2s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_2s", R"(Z_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_2s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_2s@LHCb", R"(Z_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::Z_2s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_2c(q2)", R"(Z_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_2c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_2c(q2)@LHCb", R"(Z_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_2c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_2c", R"(Z_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_2c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_2c@LHCb", R"(Z_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::Z_2c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_3(q2)", R"(Z_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_3,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_3(q2)@LHCb", R"(Z_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_3(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_3", R"(Z_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_3,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_3@LHCb", R"(Z_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::Z_3>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_4(q2)", R"(Z_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_4,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_4(q2)@LHCb", R"(Z_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::Z_4(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_4", R"(Z_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_4,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_4@LHCb", R"(Z_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::Z_4>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_5(q2)", R"(Z_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_5,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_5(q2)@LHCb", R"(Z_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_5(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_5", R"(Z_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_5,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_5@LHCb", R"(Z_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::Z_5>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_6s(q2)", R"(Z_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_6s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_6s(q2)@LHCb", R"(Z_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::Z_6s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_6s", R"(Z_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_6s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_6s@LHCb", R"(Z_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::Z_6s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_6c(q2)", R"(Z_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_6c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_6c(q2)@LHCb", R"(Z_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::Z_6c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_6c", R"(Z_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_6c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_6c@LHCb", R"(Z_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::Z_6c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_7(q2)", R"(Z_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_7,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_7(q2)@LHCb", R"(Z_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::Z_7(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_7", R"(Z_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_7,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_7@LHCb", R"(Z_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::Z_7>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_8(q2)", R"(Z_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_8,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_8(q2)@LHCb", R"(Z_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_8(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_8", R"(Z_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_8,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_8@LHCb", R"(Z_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::Z_8>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::Z_9(q2)", R"(Z_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_Z_9,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_9(q2)@LHCb", R"(Z_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::Z_9(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::Z_9", R"(Z_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_Z_9,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::Z_9@LHCb", R"(Z_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::Z_9>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_1s(q2)", R"(A_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_1s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_1s(q2)@LHCb", R"(A_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_1s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_1s", R"(A_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_1s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_1s@LHCb", R"(A_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::A_1s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_1c(q2)", R"(A_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_1c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_1c(q2)@LHCb", R"(A_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_1c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_1c", R"(A_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_1c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_1c@LHCb", R"(A_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::A_1c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_2s(q2)", R"(A_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_2s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_2s(q2)@LHCb", R"(A_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_2s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_2s", R"(A_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_2s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_2s@LHCb", R"(A_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::A_2s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_2c(q2)", R"(A_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_2c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_2c(q2)@LHCb", R"(A_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_2c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_2c", R"(A_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_2c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_2c@LHCb", R"(A_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::A_2c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_3(q2)", R"(A_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_3,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_3(q2)@LHCb", R"(A_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_3(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_3", R"(A_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_3,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_3@LHCb", R"(A_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::A_3>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_4(q2)", R"(A_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_4,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_4(q2)@LHCb", R"(A_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::A_4(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_4", R"(A_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_4,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_4@LHCb", R"(A_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::A_4>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_5(q2)", R"(A_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_5,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_5(q2)@LHCb", R"(A_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_5(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_5", R"(A_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_5,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_5@LHCb", R"(A_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::A_5>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_6s(q2)", R"(A_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_6s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_6s(q2)@LHCb", R"(A_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::A_6s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_6s", R"(A_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_6s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_6s@LHCb", R"(A_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::A_6s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_6c(q2)", R"(A_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_6c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_6c(q2)@LHCb", R"(A_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::A_6c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_6c", R"(A_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_6c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_6c@LHCb", R"(A_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::A_6c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_7(q2)", R"(A_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_7,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_7(q2)@LHCb", R"(A_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::A_7(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_7", R"(A_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_7,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_7@LHCb", R"(A_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::A_7>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_8(q2)", R"(A_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_8,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_8(q2)@LHCb", R"(A_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_8(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_8", R"(A_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_8,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_8@LHCb", R"(A_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::A_8>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::A_9(q2)", R"(A_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_A_9,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_9(q2)@LHCb", R"(A_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::A_9(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::A_9", R"(A_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_A_9,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::A_9@LHCb", R"(A_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::A_9>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_1s(q2)", R"(S_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_1s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_1s(q2)@LHCb", R"(S_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_1s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_1s", R"(S_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_1s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_1s@LHCb", R"(S_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::S_1s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_1c(q2)", R"(S_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_1c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_1c(q2)@LHCb", R"(S_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_1c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_1c", R"(S_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_1c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_1c@LHCb", R"(S_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::S_1c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_2s(q2)", R"(S_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_2s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_2s(q2)@LHCb", R"(S_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_2s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_2s", R"(S_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_2s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_2s@LHCb", R"(S_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::S_2s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_2c(q2)", R"(S_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_2c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_2c(q2)@LHCb", R"(S_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_2c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_2c", R"(S_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_2c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_2c@LHCb", R"(S_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::S_2c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_3(q2)", R"(S_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_3,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_3(q2)@LHCb", R"(S_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_3(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_3", R"(S_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_3,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_3@LHCb", R"(S_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::S_3>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_4(q2)", R"(S_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_4,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_4(q2)@LHCb", R"(S_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::S_4(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_4", R"(S_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_4,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_4@LHCb", R"(S_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::S_4>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_5(q2)", R"(S_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_5,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_5(q2)@LHCb", R"(S_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_5(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_5", R"(S_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_5,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_5@LHCb", R"(S_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::S_5>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_6s(q2)", R"(S_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_6s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_6s(q2)@LHCb", R"(S_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::S_6s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_6s", R"(S_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_6s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_6s@LHCb", R"(S_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::S_6s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_6c(q2)", R"(S_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_6c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_6c(q2)@LHCb", R"(S_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::S_6c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_6c", R"(S_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_6c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_6c@LHCb", R"(S_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::S_6c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_7(q2)", R"(S_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_7,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_7(q2)@LHCb", R"(S_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::S_7(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_7", R"(S_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_7,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_7@LHCb", R"(S_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::S_7>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_8(q2)", R"(S_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_8,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_8(q2)@LHCb", R"(S_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_8(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_8", R"(S_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_8,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_8@LHCb", R"(S_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::S_8>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::S_9(q2)", R"(S_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_S_9,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_9(q2)@LHCb", R"(S_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::S_9(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::S_9", R"(S_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_S_9,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::S_9@LHCb", R"(S_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::S_9>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_1s(q2)", R"(K_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_1s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_1s(q2)@LHCb", R"(K_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_1s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_1s", R"(K_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_1s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_1s@LHCb", R"(K_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::K_1s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_1c(q2)", R"(K_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_1c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_1c(q2)@LHCb", R"(K_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_1c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_1c", R"(K_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_1c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_1c@LHCb", R"(K_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::K_1c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_2s(q2)", R"(K_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_2s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_2s(q2)@LHCb", R"(K_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_2s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_2s", R"(K_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_2s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_2s@LHCb", R"(K_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::K_2s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_2c(q2)", R"(K_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_2c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_2c(q2)@LHCb", R"(K_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_2c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_2c", R"(K_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_2c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_2c@LHCb", R"(K_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::K_2c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_3(q2)", R"(K_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_3,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_3(q2)@LHCb", R"(K_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_3(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_3", R"(K_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_3,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_3@LHCb", R"(K_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::K_3>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_4(q2)", R"(K_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_4,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_4(q2)@LHCb", R"(K_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::K_4(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_4", R"(K_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_4,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_4@LHCb", R"(K_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::K_4>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_5(q2)", R"(K_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_5,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_5(q2)@LHCb", R"(K_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_5(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_5", R"(K_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_5,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_5@LHCb", R"(K_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::K_5>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_6s(q2)", R"(K_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_6s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_6s(q2)@LHCb", R"(K_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::K_6s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_6s", R"(K_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_6s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_6s@LHCb", R"(K_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::K_6s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_6c(q2)", R"(K_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_6c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_6c(q2)@LHCb", R"(K_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::K_6c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_6c", R"(K_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_6c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_6c@LHCb", R"(K_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::K_6c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_7(q2)", R"(K_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_7,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_7(q2)@LHCb", R"(K_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::K_7(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_7", R"(K_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_7,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_7@LHCb", R"(K_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::K_7>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_8(q2)", R"(K_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_8,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_8(q2)@LHCb", R"(K_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_8(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_8", R"(K_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_8,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_8@LHCb", R"(K_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::K_8>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::K_9(q2)", R"(K_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_K_9,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_9(q2)@LHCb", R"(K_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::K_9(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::K_9", R"(K_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_K_9,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::K_9@LHCb", R"(K_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::K_9>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_1s(q2)", R"(W_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_1s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_1s(q2)@LHCb", R"(W_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_1s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_1s", R"(W_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_1s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_1s@LHCb", R"(W_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::W_1s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_1c(q2)", R"(W_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_1c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_1c(q2)@LHCb", R"(W_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_1c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_1c", R"(W_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_1c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_1c@LHCb", R"(W_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::W_1c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_2s(q2)", R"(W_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_2s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_2s(q2)@LHCb", R"(W_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_2s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_2s", R"(W_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_2s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_2s@LHCb", R"(W_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::W_2s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_2c(q2)", R"(W_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_2c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_2c(q2)@LHCb", R"(W_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_2c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_2c", R"(W_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_2c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_2c@LHCb", R"(W_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::W_2c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_3(q2)", R"(W_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_3,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_3(q2)@LHCb", R"(W_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_3(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_3", R"(W_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_3,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_3@LHCb", R"(W_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::W_3>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_4(q2)", R"(W_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_4,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_4(q2)@LHCb", R"(W_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::W_4(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_4", R"(W_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_4,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_4@LHCb", R"(W_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::W_4>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_5(q2)", R"(W_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_5,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_5(q2)@LHCb", R"(W_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_5(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_5", R"(W_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_5,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_5@LHCb", R"(W_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::W_5>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_6s(q2)", R"(W_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_6s,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_6s(q2)@LHCb", R"(W_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::W_6s(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_6s", R"(W_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_6s,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_6s@LHCb", R"(W_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::W_6s>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_6c(q2)", R"(W_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_6c,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_6c(q2)@LHCb", R"(W_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::W_6c(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_6c", R"(W_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_6c,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_6c@LHCb", R"(W_{6c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::W_6c>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_7(q2)", R"(W_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_7,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_7(q2)@LHCb", R"(W_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::W_7(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_7", R"(W_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_7,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_7@LHCb", R"(W_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::W_7>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_8(q2)", R"(W_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_8,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_8(q2)@LHCb", R"(W_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_8(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_8", R"(W_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_8,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_8@LHCb", R"(W_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "(<<B_s->phill::W_8>> / <<B_s->phill::Gamma>>)"),


                make_observable("B_s->phill::W_9(q2)", R"(W_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::differential_W_9,
                                std::make_tuple("q2"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_9(q2)@LHCb", R"(W_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "-1.0 * (<<B_s->phill::W_9(q2)>> / <<B_s->phill::Gamma(q2)>>)"),
                make_observable("B_s->phill::W_9", R"(W_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                Unit::None(),
                                &BsToPhiDileptonAndConjugate::integrated_W_9,
                                std::make_tuple("q2_min", "q2_max"), Options{ { "q"_ok, "s" } }),
                make_expression_observable("B_s->phill::W_9@LHCb", R"(W_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(), 
                                        "-1.0 * (<<B_s->phill::W_9>> / <<B_s->phill::Gamma>>)"),


                make_expression_observable("B_s->phill::M_1s(q2)@LHCb", R"(M_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_1s(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::M_1s@LHCb", R"(M_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_1s@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::M_1c(q2)@LHCb", R"(M_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::H_1c(q2)@LHCb>> / (<<B_s->phill::K_2c(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::M_1c@LHCb", R"(M_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::H_1c@LHCb>> / (<<B_s->phill::K_2c@LHCb>>) )"),


                make_expression_observable("B_s->phill::M_2s(q2)@LHCb", R"(M_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_2s(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::M_2s@LHCb", R"(M_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_2s@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::M_2c(q2)@LHCb", R"(M_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::H_2c(q2)@LHCb>> / (<<B_s->phill::K_2c(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::M_2c@LHCb", R"(M_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::H_2c@LHCb>> / (<<B_s->phill::K_2c@LHCb>>) )"),


                make_expression_observable("B_s->phill::M_3(q2)@LHCb", R"(M_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::H_3(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::M_3@LHCb", R"(M_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::H_3@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::M_4(q2)@LHCb", R"(M_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::H_4(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * (2.0 * <<B_s->phill::K_2s(q2)@LHCb>> - <<B_s->phill::K_3(q2)@LHCb>>) ) ^ 0.5) )"),
                make_expression_observable("B_s->phill::M_4@LHCb", R"(M_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::H_4@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * (2.0 * <<B_s->phill::K_2s@LHCb>> - <<B_s->phill::K_3@LHCb>>) ) ^ 0.5) )"),


                make_expression_observable("B_s->phill::M_5(q2)@LHCb", R"(M_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_5(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * (2.0 * <<B_s->phill::K_2s(q2)@LHCb>> - <<B_s->phill::K_3(q2)@LHCb>>) ) ^ 0.5) )"),
                make_expression_observable("B_s->phill::M_5@LHCb", R"(M_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_5@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * (2.0 * <<B_s->phill::K_2s@LHCb>> - <<B_s->phill::K_3@LHCb>>) ) ^ 0.5) )"),


                make_expression_observable("B_s->phill::M_6s(q2)@LHCb", R"(M_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::H_6s(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::M_6s@LHCb", R"(M_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::H_6s@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::M_7(q2)@LHCb", R"(M_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::H_7(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * (2.0 * <<B_s->phill::K_2s(q2)@LHCb>> - <<B_s->phill::K_3(q2)@LHCb>>) ) ^ 0.5) )"),
                make_expression_observable("B_s->phill::M_7@LHCb", R"(M_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::H_7@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * (2.0 * <<B_s->phill::K_2s@LHCb>> - <<B_s->phill::K_3@LHCb>>) ) ^ 0.5) )"),


                make_expression_observable("B_s->phill::M_8(q2)@LHCb", R"(M_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(1.0/(2.0^(0.5)) * <<B_s->phill::H_8(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * (2.0 * <<B_s->phill::K_2s(q2)@LHCb>> - <<B_s->phill::K_3(q2)@LHCb>>) ) ^ 0.5) )"),
                make_expression_observable("B_s->phill::M_8@LHCb", R"(M_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(1.0/(2.0^(0.5)) * <<B_s->phill::H_8@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * (2.0 * <<B_s->phill::K_2s@LHCb>> - <<B_s->phill::K_3@LHCb>>) ) ^ 0.5) )"),


                make_expression_observable("B_s->phill::M_9(q2)@LHCb", R"(M_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::H_9(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::M_9@LHCb", R"(M_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::H_9@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::Q_1s(q2)@LHCb", R"(Q_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_1s(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::Q_1s@LHCb", R"(Q_{1s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_1s@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::Q_1c(q2)@LHCb", R"(Q_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::Z_1c(q2)@LHCb>> / (<<B_s->phill::K_2c(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::Q_1c@LHCb", R"(Q_{1c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::Z_1c@LHCb>> / (<<B_s->phill::K_2c@LHCb>>) )"),


                make_expression_observable("B_s->phill::Q_2s(q2)@LHCb", R"(Q_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_2s(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::Q_2s@LHCb", R"(Q_{2s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_2s@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::Q_2c(q2)@LHCb", R"(Q_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::Z_2c(q2)@LHCb>> / (<<B_s->phill::K_2c(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::Q_2c@LHCb", R"(Q_{2c}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::Z_2c@LHCb>> / (<<B_s->phill::K_2c@LHCb>>) )"),


                make_expression_observable("B_s->phill::Q_3(q2)@LHCb", R"(Q_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::Z_3(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::Q_3@LHCb", R"(Q_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::Z_3@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::Q_4(q2)@LHCb", R"(Q_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::Z_4(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * (2.0 * <<B_s->phill::K_2s(q2)@LHCb>> - <<B_s->phill::K_3(q2)@LHCb>>) ) ^ 0.5) )"),
                make_expression_observable("B_s->phill::Q_4@LHCb", R"(Q_{4}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::Z_4@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * (2.0 * <<B_s->phill::K_2s@LHCb>> - <<B_s->phill::K_3@LHCb>>) ) ^ 0.5) )"),


                make_expression_observable("B_s->phill::Q_5(q2)@LHCb", R"(Q_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_5(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * (2.0 * <<B_s->phill::K_2s(q2)@LHCb>> - <<B_s->phill::K_3(q2)@LHCb>>) ) ^ 0.5) )"),
                make_expression_observable("B_s->phill::Q_5@LHCb", R"(Q_{5}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_5@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * (2.0 * <<B_s->phill::K_2s@LHCb>> - <<B_s->phill::K_3@LHCb>>) ) ^ 0.5) )"),


                make_expression_observable("B_s->phill::Q_6s(q2)@LHCb", R"(Q_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::Z_6s(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::Q_6s@LHCb", R"(Q_{6s}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::Z_6s@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::Q_7(q2)@LHCb", R"(Q_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::Z_7(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * (2.0 * <<B_s->phill::K_2s(q2)@LHCb>> - <<B_s->phill::K_3(q2)@LHCb>>) ) ^ 0.5) )"),
                make_expression_observable("B_s->phill::Q_7@LHCb", R"(Q_{7}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-1.0 * <<B_s->phill::Z_7@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * (2.0 * <<B_s->phill::K_2s@LHCb>> - <<B_s->phill::K_3@LHCb>>) ) ^ 0.5) )"),


                make_expression_observable("B_s->phill::Q_8(q2)@LHCb", R"(Q_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "((1.0/2.0^(0.5)) * <<B_s->phill::Z_8(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * (2.0 * <<B_s->phill::K_2s(q2)@LHCb>> - <<B_s->phill::K_3(q2)@LHCb>>) ) ^ 0.5) )"),
                make_expression_observable("B_s->phill::Q_8@LHCb", R"(Q_{8}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "((1.0/2.0^(0.5)) * <<B_s->phill::Z_8@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * (2.0 * <<B_s->phill::K_2s@LHCb>> - <<B_s->phill::K_3@LHCb>>) ) ^ 0.5) )"),


                make_expression_observable("B_s->phill::Q_9(q2)@LHCb", R"(Q_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::Z_9(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::Q_9@LHCb", R"(Q_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::Z_9@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::SP_1(q2)@LHCb", R"(SP_{1}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::S_3(q2)@LHCb>> / (<<B_s->phill::S_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::SP_1@LHCb", R"(SP_{1}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::S_3@LHCb>> / (<<B_s->phill::S_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::SP_2(q2)@LHCb", R"(SP_{2}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.125 * <<B_s->phill::S_6s(q2)@LHCb>> / (<<B_s->phill::S_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::SP_2@LHCb", R"(SP_{2}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.125 * <<B_s->phill::S_6s@LHCb>> / (<<B_s->phill::S_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::SP_3(q2)@LHCb", R"(SP_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.25 * <<B_s->phill::S_9(q2)@LHCb>> / (<<B_s->phill::S_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::SP_3@LHCb", R"(SP_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.25 * <<B_s->phill::S_9@LHCb>> / (<<B_s->phill::S_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::SP_4p(q2)@LHCb", R"(SP_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_4(q2)@LHCb>> / ((-1.0 * <<B_s->phill::S_2c(q2)@LHCb>> * <<B_s->phill::S_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::SP_4p@LHCb", R"(SP_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_4@LHCb>> / ((-1.0 * <<B_s->phill::S_2c@LHCb>> * <<B_s->phill::S_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::SP_5p(q2)@LHCb", R"(SP_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::S_5(q2)@LHCb>> / ((-1.0 * <<B_s->phill::S_2c(q2)@LHCb>> * <<B_s->phill::S_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::SP_5p@LHCb", R"(SP_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::S_5@LHCb>> / ((-1.0 * <<B_s->phill::S_2c@LHCb>> * <<B_s->phill::S_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::SP_6p(q2)@LHCb", R"(SP_{6p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::S_7(q2)@LHCb>> / ((-1.0 * <<B_s->phill::S_2c(q2)@LHCb>> * <<B_s->phill::S_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::SP_6p@LHCb", R"(SP_{6p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::S_7@LHCb>> / ((-1.0 * <<B_s->phill::S_2c@LHCb>> * <<B_s->phill::S_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::SP_8p(q2)@LHCb", R"(SP_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_8(q2)@LHCb>> / ((-1.0 * <<B_s->phill::S_2c(q2)@LHCb>> * <<B_s->phill::S_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::SP_8p@LHCb", R"(SP_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_8@LHCb>> / ((-1.0 * <<B_s->phill::S_2c@LHCb>> * <<B_s->phill::S_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::AP_1(q2)@LHCb", R"(AP_{1}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::A_3(q2)@LHCb>> / (<<B_s->phill::S_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::AP_1@LHCb", R"(AP_{1}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::A_3@LHCb>> / (<<B_s->phill::S_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::AP_2(q2)@LHCb", R"(AP_{2}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.125 * <<B_s->phill::A_6s(q2)@LHCb>> / (<<B_s->phill::S_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::AP_2@LHCb", R"(AP_{2}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.125 * <<B_s->phill::A_6s@LHCb>> / (<<B_s->phill::S_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::AP_3(q2)@LHCb", R"(AP_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.25 * <<B_s->phill::A_9(q2)@LHCb>> / (<<B_s->phill::S_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::AP_3@LHCb", R"(AP_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.25 * <<B_s->phill::A_9@LHCb>> / (<<B_s->phill::S_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::AP_4p(q2)@LHCb", R"(AP_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_4(q2)@LHCb>> / ((-1.0 * <<B_s->phill::S_2c(q2)@LHCb>> * <<B_s->phill::S_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::AP_4p@LHCb", R"(AP_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_4@LHCb>> / ((-1.0 * <<B_s->phill::S_2c@LHCb>> * <<B_s->phill::S_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::AP_5p(q2)@LHCb", R"(AP_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::A_5(q2)@LHCb>> / ((-1.0 * <<B_s->phill::S_2c(q2)@LHCb>> * <<B_s->phill::S_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::AP_5p@LHCb", R"(AP_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::A_5@LHCb>> / ((-1.0 * <<B_s->phill::S_2c@LHCb>> * <<B_s->phill::S_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::AP_6p(q2)@LHCb", R"(AP_{6p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::A_7(q2)@LHCb>> / ((-1.0 * <<B_s->phill::S_2c(q2)@LHCb>> * <<B_s->phill::S_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::AP_6p@LHCb", R"(AP_{6p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::A_7@LHCb>> / ((-1.0 * <<B_s->phill::S_2c@LHCb>> * <<B_s->phill::S_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::AP_8p(q2)@LHCb", R"(AP_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_8(q2)@LHCb>> / ((-1.0 * <<B_s->phill::S_2c(q2)@LHCb>> * <<B_s->phill::S_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::AP_8p@LHCb", R"(AP_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_8@LHCb>> / ((-1.0 * <<B_s->phill::S_2c@LHCb>> * <<B_s->phill::S_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::KP_1(q2)@LHCb", R"(KP_{1}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::K_3(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::KP_1@LHCb", R"(KP_{1}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::K_3@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::KP_2(q2)@LHCb", R"(KP_{2}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.125 * <<B_s->phill::K_6s(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::KP_2@LHCb", R"(KP_{2}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.125 * <<B_s->phill::K_6s@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::KP_3(q2)@LHCb", R"(KP_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.25 * <<B_s->phill::K_9(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::KP_3@LHCb", R"(KP_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.25 * <<B_s->phill::K_9@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::KP_4p(q2)@LHCb", R"(KP_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_4(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::KP_4p@LHCb", R"(KP_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_4@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::KP_5p(q2)@LHCb", R"(KP_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::K_5(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::KP_5p@LHCb", R"(KP_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::K_5@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::KP_6p(q2)@LHCb", R"(KP_{6p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::K_7(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::KP_6p@LHCb", R"(KP_{6p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::K_7@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::KP_8p(q2)@LHCb", R"(KP_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_8(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::KP_8p@LHCb", R"(KP_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_8@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::WP_1(q2)@LHCb", R"(WP_{1}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::W_3(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::WP_1@LHCb", R"(WP_{1}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::W_3@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::WP_2(q2)@LHCb", R"(WP_{2}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.125 * <<B_s->phill::W_6s(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::WP_2@LHCb", R"(WP_{2}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.125 * <<B_s->phill::W_6s@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::WP_3(q2)@LHCb", R"(WP_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.25 * <<B_s->phill::W_9(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::WP_3@LHCb", R"(WP_{3}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.25 * <<B_s->phill::W_9@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::WP_4p(q2)@LHCb", R"(WP_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_4(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::WP_4p@LHCb", R"(WP_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_4@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::WP_5p(q2)@LHCb", R"(WP_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::W_5(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::WP_5p@LHCb", R"(WP_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::W_5@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::WP_6p(q2)@LHCb", R"(WP_{6p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::W_7(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::WP_6p@LHCb", R"(WP_{6p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::W_7@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::WP_8p(q2)@LHCb", R"(WP_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_8(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::WP_8p@LHCb", R"(WP_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_8@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::M_4p(q2)@LHCb", R"(M_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_4(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::M_4p@LHCb", R"(M_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_4@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::M_5p(q2)@LHCb", R"(M_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::H_5(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::M_5p@LHCb", R"(M_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::H_5@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::M_7p(q2)@LHCb", R"(M_{7p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::H_7(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::M_7p@LHCb", R"(M_{7p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::H_7@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::M_8p(q2)@LHCb", R"(M_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_8(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::M_8p@LHCb", R"(M_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::H_8@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::Q_4p(q2)@LHCb", R"(Q_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_4(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::Q_4p@LHCb", R"(Q_{4p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_4@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::Q_5p(q2)@LHCb", R"(Q_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::Z_5(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::Q_5p@LHCb", R"(Q_{5p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(0.5 * <<B_s->phill::Z_5@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::Q_7p(q2)@LHCb", R"(Q_{7p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::Z_7(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::Q_7p@LHCb", R"(Q_{7p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(-0.5 * <<B_s->phill::Z_7@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::Q_8p(q2)@LHCb", R"(Q_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_8(q2)@LHCb>> / ((-1.0 * <<B_s->phill::K_2c(q2)@LHCb>> * <<B_s->phill::K_2s(q2)@LHCb>>) ^ 0.5) )"),
                make_expression_observable("B_s->phill::Q_8p@LHCb", R"(Q_{8p}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::Z_8@LHCb>> / ((-1.0 * <<B_s->phill::K_2c@LHCb>> * <<B_s->phill::K_2s@LHCb>>) ^ 0.5) )"),


                make_expression_observable("B_s->phill::SS(q2)@LHCb", R"(SS{}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_6c(q2)@LHCb>> / (<<B_s->phill::S_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::SS@LHCb", R"(SS{}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::S_6c@LHCb>> / (<<B_s->phill::S_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::KS(q2)@LHCb", R"(KS{}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_6c(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::KS@LHCb", R"(KS{}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::K_6c@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::WS(q2)@LHCb", R"(WS{}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_6c(q2)@LHCb>> / (<<B_s->phill::K_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::WS@LHCb", R"(WS{}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::W_6c@LHCb>> / (<<B_s->phill::K_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::AS(q2)@LHCb", R"(AS{}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_6c(q2)@LHCb>> / (<<B_s->phill::S_2s(q2)@LHCb>>) )"),
                make_expression_observable("B_s->phill::AS@LHCb", R"(AS{}(\bar{B}_s\to \phi\ell^+\ell^-))",
                                        Unit::None(),
                                        "(<<B_s->phill::A_6c@LHCb>> / (<<B_s->phill::S_2s@LHCb>>) )"),


                make_expression_observable("B_s->phill::A_FB@LHCb", R"(A_\mathrm{FB}^\mathrm{LHCb}(\bar{B}_s\to \phi\ell^+\ell^-))",
                        Unit::None(),
                        R"( -1.0 * <<B_s->phill::A_FB>> )"),

                make_expression_observable("B_s->phill::A_FB(q2)@LHCb", R"(A_\mathrm{FB}^\mathrm{LHCb}(\bar{B}_s\to \phi\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"( -1.0 * <<B_s->phill::A_FB(q2)>> )"),

                make_expression_observable("B_s->phill::R_phi", R"(R_\phi)",
                        Unit::None(),
                        R"(
                        <<B_s->phill::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        /
                        <<B_s->phill::BR;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                        )"),

                make_expression_observable("B_s->phill::NormalizedBR", R"(\mathcal{B}(\bar{B}_s\to \phi\ell^+\ell^-)/\mathcal{B}(\bar{B}_s\to\phi J/\psi))",
                        Unit::None(),
                        R"(
                        <<B_s->phill::BR>> / <<B_s->phipsi::BR;psi=J/psi>>
                        )"),

                make_observable("B_s->phill::A_para_left_real", R"(Re(A_\parallel^L)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_para_left_real,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_para_right_real", R"(Re(A_\parallel^R)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_para_right_real,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_perp_left_real", R"(Re(A_\perp^L)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_perp_left_real,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_perp_right_real", R"(Re(A_\perp^R)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_perp_right_real,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_0_left_real", R"(Re(A_0^L)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_long_left_real,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_0_right_real", R"(Re(A_0^R)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_long_right_real,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_time_real", R"(Re(A_t)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_time_real,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_scal_real", R"(Re(A_S)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_scal_real,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_para_left_imag", R"(Im(A_\parallel^L)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_para_left_imag,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_para_right_imag", R"(Im(A_\parallel^R)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_para_right_imag,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_perp_left_imag", R"(Im(A_\perp^L)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_perp_left_imag,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_perp_right_imag", R"(Im(A_\perp^R)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_perp_right_imag,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_0_left_imag", R"(Im(A_0^L)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_long_left_imag,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_0_right_imag", R"(Im(A_0^R)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_long_right_imag,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_time_imag", R"(Im(A_t)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_time_imag,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::A_scal_imag", R"(Im(A_S)(B_s^0 \rightarrow \phi \ell^+\ell^-))",
                        Unit::None(),
                        &BsToPhiDileptonAndConjugate::a_scal_imag,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),


                // make_expression_observable("B_s->phill::Q_8^-(q2)@LHCb", R"(Q_{8}^{-}(\bar{B}_s\to \phi\ell^+\ell^-))",
                //         Unit::None(),
                //         R"( 
                //         (<<B_s->phill::Z_8(q2)>>) / (
                //                 -2 * (<<B_s->phill::J_2c(q2);cp-conjugate=false>> + <<B_s->phill::J_2c(q2);cp-conjugate=true>>) * (
                //                 2 * (<<B_s->phill::J_2s(q2);cp-conjugate=false>> + <<B_s->phill::J_2s(q2);cp-conjugate=true>>) - 
                //                 (<<B_s->phill::J_3(q2);cp-conjugate=false>> + <<B_s->phill::J_3(q2);cp-conjugate=true>>)
                //                 )
                //             )^(0.5)
                //         )" ),

                // make_expression_observable("B_s->phill::Q_9(q2)@LHCb", R"(Q_{9}(\bar{B}_s\to \phi\ell^+\ell^-))",
                //         Unit::None(),
                //         R"( 
                //         -1. * (<<B_s->phill::Z_9(q2)>>) / (
                //                 -2. * (<<B_s->phill::J_2s(q2);cp-conjugate=false>> + <<B_s->phill::J_2s(q2);cp-conjugate=true>>) 
                //         )
                //         )" ),

                make_expression_observable("B_s->phill::expBR", R"(\langle\mathcal{B}\rangle(\bar{B}\to \bar{K}^*\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                         <<life_time::B_s>> / <<QM::hbar>>* 0.5 / (1.0 - <<B_s::ys>>^2.0) * (
                                        2.0 * (<<B_s->phill::J_1s;cp-conjugate=false>> + <<B_s->phill::J_1s;cp-conjugate=true>> - <<B_s::ys>> * <<B_s->phill::H_1s>>)
                                            + (<<B_s->phill::J_1c;cp-conjugate=false>> + <<B_s->phill::J_1c;cp-conjugate=true>> - <<B_s::ys>> * <<B_s->phill::H_1c>>)
                         - 1.0 / 3.0 * (2.0 * (<<B_s->phill::J_2s;cp-conjugate=false>> + <<B_s->phill::J_2s;cp-conjugate=true>> - <<B_s::ys>> * <<B_s->phill::H_2s>>)
                                            + (<<B_s->phill::J_2c;cp-conjugate=false>> + <<B_s->phill::J_2c;cp-conjugate=true>> - <<B_s::ys>> * <<B_s->phill::H_2c>>))
                        )
                        )"),

                make_expression_observable("B_s->phill::NormalizedexpBR", R"(\langle\mathcal{B}\rangle(B_s->\phi\ell\ell)/\mathcal{B}(B_s->\phi J/\psi))",
                        Unit::None(),
                        R"(
                        <<B_s->phill::expBR>> / <<B_s->phipsi::BR;psi=J/psi>>
                        )"),

                make_observable("B_s->phill::Re{C9_perp}(q2)",
                        Unit::None(),
                        &BsToPhiDilepton::real_C9_perp,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::Re{C9_para}(q2)",
                        Unit::None(),
                        &BsToPhiDilepton::real_C9_para,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::Im{C9_perp}(q2)",
                        Unit::None(),
                        &BsToPhiDilepton::imag_C9_perp,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),

                make_observable("B_s->phill::Im{C9_para}(q2)",
                        Unit::None(),
                        &BsToPhiDilepton::imag_C9_para,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}


    // Lambda_b -> Lambda l^+ l^-
    // {{{
    ObservableGroup
    make_lambdab_to_lambda_ll_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\Lambda_b \to \Lambda\ell^+\ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavor.)",
            {
                // Lambda_b -> Lambda l^+ l^-, Large Recoil
                make_observable("Lambda_b->Lambdall::dBR/dq2@LargeRecoil", R"(d\mathcal{B}(\Lambda_b\to\Lambda\ell^+\ell^-)/dq^2)",
                        Unit::InverseGeV2(),
                        &LambdaBToLambdaDilepton<LargeRecoil>::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^l(q2)@LargeRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_leptonic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^h(q2)@LargeRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_hadronic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^c(q2)@LargeRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_combined,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::F_0(q2)@LargeRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LargeRecoil>::differential_fzero,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::BR@LargeRecoil", R"(\mathcal{B}(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LargeRecoil>::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("Lambda_b->Lambdall::R_Lambda@LargeRecoil", R"(R_{\Lambda})",
                        Unit::None(),
                        R"(
                        <<Lambda_b->Lambdall::BR@LargeRecoil;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        /
                        <<Lambda_b->Lambdall::BR@LargeRecoil;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                        )"),

                make_observable("Lambda_b->Lambdall::A_FB^l@LargeRecoil", R"(A_\mathrm{FB}^\ell(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^h@LargeRecoil", R"(A_\mathrm{FB}^h(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_hadronic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^c@LargeRecoil", R"(A_\mathrm{FB}^{h,\ell}(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_combined,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::F_0@LargeRecoil", R"(F_0(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LargeRecoil>::integrated_fzero,
                        std::make_tuple("q2_min", "q2_max")),

                // Lambda_b -> Lambda l^+ l^-, Low Recoil
                make_observable("Lambda_b->Lambdall::dBR/dq2@LowRecoil", R"(d\mathcal{B}(\Lambda_b\to\Lambda\ell^+\ell^-)/dq^2)",
                        Unit::InverseGeV2(),
                        &LambdaBToLambdaDilepton<LowRecoil>::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^l(q2)@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_leptonic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^h(q2)@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_hadronic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^c(q2)@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_combined,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::F_0(q2)@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::differential_fzero,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::BR@LowRecoil", R"(\mathcal{B}(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^l@LowRecoil", R"(A_\mathrm{FB}^\ell(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^h@LowRecoil", R"(A_\mathrm{FB}^h(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_hadronic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^c@LowRecoil", R"(A_\mathrm{FB}^{h,\ell}(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_combined,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::F_0@LowRecoil", R"(F_0(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_fzero,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_1ss@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k1ss,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_1cc@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k1cc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_1c@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k1c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_2ss@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k2ss,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_2cc@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k2cc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_2c@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k2c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_3sc@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k3sc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_3s@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k3s,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_4sc@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k4sc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_4s@LowRecoil",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k4s,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_1@LowRecoil", R"(M_1)",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m1,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_2@LowRecoil", R"(M_2)",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m2,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_3@LowRecoil", R"(M_3)",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m3,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_4@LowRecoil", R"(M_4)",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_5@LowRecoil", R"(M_5)",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_6@LowRecoil", R"(M_6)",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m6,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_7@LowRecoil", R"(M_7)",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m7,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_8@LowRecoil", R"(M_8)",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m8,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_9@LowRecoil", R"(M_9)",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m9,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_10@LowRecoil", R"(M_{10})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m10,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_11@LowRecoil", R"(M_{11})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m11,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_12@LowRecoil", R"(M_{12})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m12,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_13@LowRecoil", R"(M_{13})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m13,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_14@LowRecoil", R"(M_{14})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m14,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_15@LowRecoil", R"(M_{15})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m15,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_16@LowRecoil", R"(M_{16})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m16,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_17@LowRecoil", R"(M_{17})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m17,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_18@LowRecoil", R"(M_{18})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m18,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_19@LowRecoil", R"(M_{19})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m19,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_20@LowRecoil", R"(M_{20})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m20,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_21@LowRecoil", R"(M_{21})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m21,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_22@LowRecoil", R"(M_{22})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m22,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_23@LowRecoil", R"(M_{23})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m23,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_24@LowRecoil", R"(M_{24})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m24,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_25@LowRecoil", R"(M_{25})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m25,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_26@LowRecoil", R"(M_{26})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m26,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_27@LowRecoil", R"(M_{27})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m27,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_28@LowRecoil", R"(M_{28})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m28,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_29@LowRecoil", R"(M_{29})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m29,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_30@LowRecoil", R"(M_{30})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m30,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_31@LowRecoil", R"(M_{31})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m31,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_32@LowRecoil", R"(M_{32})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m32,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_33@LowRecoil", R"(M_{33})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m33,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_34@LowRecoil", R"(M_{34})",
                        Unit::None(),
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m34,
                        std::make_tuple("q2_min", "q2_max")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}


    // Lambda_b -> Lambda(1520) l^+ l^-
    // {{{
    ObservableGroup
    make_lambdab_to_lambda1520_ll_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\Lambda_b \to \Lambda(1520))\ell^+\ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavor.)",
            {
                // Lambda_b -> Lambda(1520) l^+ l^-
                make_observable("Lambda_b->Lambda(1520)ll::dBR/dq2", R"(d\mathcal{B}(\bar{\Lambda}_b\to\bar{\Lambda}(1520)\ell^+\ell^-)/dq^2)",
                        Unit::InverseGeV2(),
                        &LambdaBToLambda1520Dilepton::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda(1520)ll::Gamma_CP_specific(q2)",
                        Unit::GeV(),
                        &LambdaBToLambda1520Dilepton::differential_decay_width,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda(1520)ll::A_FB^l(q2)",
                        Unit::None(),
                        &LambdaBToLambda1520Dilepton::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_expression_observable("Lambda_b->Lambda(1520)ll::Gamma(q2)", R"(\Gamma(\bar{\Lambda}_b\to\bar{\Lambda}(1520)\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"(
                        0.5 * (<<Lambda_b->Lambda(1520)ll::Gamma_CP_specific(q2);cp-conjugate=false>> + <<Lambda_b->Lambda(1520)ll::Gamma_CP_specific(q2);cp-conjugate=true>>)
                        )"),

                make_observable("Lambda_b->Lambda(1520)ll::L_1cc(q2)", R"(L_{1cc}(\bar{\Lambda}_b\to\bar{\Lambda}(1520)\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        &LambdaBToLambda1520Dilepton::differential_L_1cc,
                        std::make_tuple("q2")),

                make_expression_observable("Lambda_b->Lambda(1520)ll::S_1cc(q2)", R"(S_{1cc}(\Lambda_b\to\Lambda(1520)\ell^+\ell^-)(q^2))",
                        Unit::None(),
                        R"(
                        0.5 * (<<Lambda_b->Lambda(1520)ll::L_1cc(q2);cp-conjugate=false>> + <<Lambda_b->Lambda(1520)ll::L_1cc(q2);cp-conjugate=true>>)
                            / <<Lambda_b->Lambda(1520)ll::Gamma(q2)>>
                        )"),

                make_observable("Lambda_b->Lambda(1520)ll::BR", R"(\mathcal{B}(\bar{\Lambda}_b\to\bar{\Lambda}(1520)\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambda1520Dilepton::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda(1520)ll::Gamma_CP_specific",
                        Unit::GeV(),
                        &LambdaBToLambda1520Dilepton::integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("Lambda_b->Lambda(1520)ll::Gamma", R"(\Gamma(\bar{\Lambda}_b\to\bar{\Lambda}(1520)\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (<<Lambda_b->Lambda(1520)ll::Gamma_CP_specific;cp-conjugate=false>> + <<Lambda_b->Lambda(1520)ll::Gamma_CP_specific;cp-conjugate=true>>)
                        )"),

                make_observable("Lambda_b->Lambda(1520)ll::A_FB^l", R"(A_\mathrm{FB}^\ell(\bar{\Lambda}_b\to\bar{\Lambda}(1520)\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambda1520Dilepton::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda(1520)ll::L_1cc", R"(L_{1cc}(\bar{\Lambda}_b\to\bar{\Lambda}(1520)\ell^+\ell^-))",
                        Unit::None(),
                        &LambdaBToLambda1520Dilepton::integrated_L_1cc,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("Lambda_b->Lambda(1520)ll::S_1cc", R"(S_{1cc}(\Lambda_b\to\Lambda(1520)\ell^+\ell^-))",
                        Unit::None(),
                        R"(
                        0.5 * (<<Lambda_b->Lambda(1520)ll::L_1cc;cp-conjugate=false>> + <<Lambda_b->Lambda(1520)ll::L_1cc;cp-conjugate=true>>)
                            / <<Lambda_b->Lambda(1520)ll::Gamma>>
                        )"),

            }
        );

        return ObservableGroup(imp);
    }
    // }}}


    // Lambda_b -> Lambda(1520) gamma
    // {{{
    ObservableGroup
    make_lambdab_to_lambda1520_gamma_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\Lambda_b \to \Lambda(1520)) \gamma$ decays)",
            R"()",
            {
                // Lambda_b -> Lambda(1520) gamma
                make_observable("Lambda_b->Lambda(1520)gamma::BR", R"(\mathcal{B}(\Lambda_b\to\Lambda(1520)\gamma))",
                        Unit::None(),
                        &LambdaBToLambda1520Gamma::branching_ratio)
            }
        );

        return ObservableGroup(imp);
    }
    // }}}


    // B -> X_s {gamma, l^+ l^-}
    // {{{
    ObservableGroup
    make_b_to_xs_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B \to X_s \lbrace \gamma, \ell^+\ell^-\rbrace$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor.)",
            {
                // B->X_s gamma
                make_observable("B->X_sgamma::BR@Minimal",
                        Unit::None(),
                        &BToXsGamma<Minimal>::integrated_branching_ratio),

                // B->X_s gamma, NLO implementation
                make_observable("B->X_sgamma::BR(E_min)@NLO",
                        Unit::None(),
                        &BToXsGamma<NLO>::integrated_branching_ratio,
                        std::make_tuple("E_min")),

                make_observable("B->X_sgamma::E_1(E_min)@NLO",
                        Unit::GeV(),
                        &BToXsGamma<NLO>::photon_energy_moment_1,
                        std::make_tuple("E_min")),

                make_observable("B->X_sgamma::E_2(E_min)@NLO",
                        Unit::GeV2(),
                        &BToXsGamma<NLO>::photon_energy_moment_2,
                        std::make_tuple("E_min")),

                // B->X_s ll, HLMW2005
                make_observable("B->X_sll::dBR/dq2@HLMW2005",
                        Unit::InverseGeV2(),
                        &BToXsDilepton<HLMW2005>::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B->X_sll::BR@HLMW2005",
                        Unit::None(),
                        &BToXsDilepton<HLMW2005>::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}


    // B -> K nu nu
    // {{{
    ObservableGroup
    make_b_to_k_nu_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to K \nu\bar\nu$ decays)",
            R"()",
            {
                make_observable("B->Knunu::dBR/dq2", R"(d\mathcal{B}(\bar{B}\to \bar{K}\nu\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarDineutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "u" }, { "P"_ok, "K" } }),
                make_observable("B->Knunu::BR", R"(\mathcal{B}(\bar{B}\to \bar{K}\nu\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarDineutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "u" }, { "P"_ok, "K" } })
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> eta nu nu
    // {{{
    ObservableGroup
    make_bs_to_eta_nu_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_s\to \eta \nu\bar\nu$ decays)",
            R"()",
            {
                make_observable("B_s->etanunu::dBR/dq2", R"(d\mathcal{B}(\bar{B}_s\to\bar{\eta}\nu\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarDineutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" }, { "P"_ok, "eta" } }),
                make_observable("B_s->etanunu::BR", R"(\mathcal{B}(\bar{B}_s\to\bar{\eta}\nu\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarDineutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" }, { "P"_ok, "eta" } }),
                make_observable("B_s->eta_primenunu::dBR/dq2", R"(d\mathcal{B}(\bar{B}_s\to\bar{\eta_prime}\nu\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarDineutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "s" }, { "P"_ok, "eta_prime" } }),
                make_observable("B_s->eta_primenunu::BR", R"(\mathcal{B}(\bar{B}_s\to\bar{\eta_prime}\nu\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarDineutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "q"_ok, "s" }, { "P"_ok, "eta_prime" } })
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> K^* nu nu
    // {{{
    ObservableGroup
    make_b_to_kstar_nu_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to K^* \nu\bar\nu$ decays)",
            R"()",
            {
                make_observable("B->K^*nunu::dBR/dq2", R"(d\mathcal{B}(\bar{B}\to \bar{K}^*\nu\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToVectorDineutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "D"_ok, "s" }, { "I"_ok, "1/2" } }),
                make_observable("B->K^*nunu::BR", R"(\mathcal{B}(\bar{B}\to \bar{K}^*\nu\bar\nu))",
                        Unit::None(),
                        &BToVectorDineutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "D"_ok, "s" }, { "I"_ok, "1/2" } }),
                make_observable("B->K^*nunu::F_L(q2)", R"(F_L(\bar{B}\to \bar{K}^*\nu\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorDineutrino::differential_longitudinal_polarisation,
                        std::make_tuple("q2"),
                        Options{ { "D"_ok, "s" }, { "I"_ok, "1/2" } }),
                make_observable("B->K^*nunu::F_L", R"(F_L(\bar{B}\to \bar{K}^*\nu\bar\nu))",
                        Unit::None(),
                        &BToVectorDineutrino::integrated_longitudinal_polarisation,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "D"_ok, "s" }, { "I"_ok, "1/2" } })
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> phi nu nu
    // {{{
    ObservableGroup
    make_bs_to_phi_nu_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_s\to\phi\nu\bar\nu$ decays)",
            R"()",
            {
                make_observable("B_s->phinunu::dBR/dq2", R"(d\mathcal{B}(\bar{B}_s\to\phi\nu\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToVectorDineutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "D"_ok, "s" }, { "q"_ok, "s" }, { "I"_ok, "0" } }),
                make_observable("B_s->phinunu::BR", R"(\mathcal{B}(\bar{B}_s\to\phi\nu\bar\nu))",
                        Unit::None(),
                        &BToVectorDineutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "D"_ok, "s" }, { "q"_ok, "s" }, { "I"_ok, "0" } }),
                make_observable("B_s->phinunu::F_L(q2)", R"(F_L(\bar{B}_s\to\phi\nu\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorDineutrino::differential_longitudinal_polarisation,
                        std::make_tuple("q2"),
                        Options{ { "D"_ok, "s" }, { "q"_ok, "s" }, { "I"_ok, "0" } }),
                make_observable("B_s->phinunu::F_L", R"(F_L(\bar{B}_s\to\phi\nu\bar\nu))",
                        Unit::None(),
                        &BToVectorDineutrino::integrated_longitudinal_polarisation,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "D"_ok, "s" }, { "q"_ok, "s" }, { "I"_ok, "0" } })
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Lambda_b -> Lambda nu nu
    // {{{
    ObservableGroup
    make_lambda_b_to_lambda_nu_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\Lambda_b\to\Lambda\nu\bar\nu$ decays)",
            R"()",
            {
                make_observable("Lambda_b->Lambdanunu::dBR/dq2", R"(d\mathcal{B}(\bar{\Lambda}_b\to\Lambda\nu\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &LambdaBToLambdaDineutrino::differential_branching_ratio,
                        std::make_tuple("q2")),
                make_observable("Lambda_b->Lambdanunu::F_L(q^2)", R"(F_L(\bar{\Lambda}_b\to\Lambda\nu\bar\nu)(q^2))",
                        Unit::InverseGeV2(),
                        &LambdaBToLambdaDineutrino::differential_longitudinal_polarisation,
                        std::make_tuple("q2")),
                make_cacheable_observable("Lambda_b->Lambdanunu::BR", R"(\mathcal{B}(\bar{\Lambda}_b\to\Lambda\nu\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaDineutrino::prepare,
                        &LambdaBToLambdaDineutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),
                make_cacheable_observable("Lambda_b->Lambdanunu::F_L", R"(F_L(\bar{\Lambda}_b\to\Lambda\nu\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaDineutrino::prepare,
                        &LambdaBToLambdaDineutrino::integrated_longitudinal_polarisation,
                        std::make_tuple("q2_min", "q2_max"))
            }
        );

        return ObservableGroup(imp);
    }
    // }}}


    ObservableSection
    make_rare_b_decays_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in rare (semi)leptonic and radiative $b$-hadron decays",
            "",
            {
                // B_q -> l^+ l^-
                make_b_to_ll_group(),

                // B_q -> M charmonium
                make_b_to_p_charmonium_group(),
                make_b_to_v_charmonium_group(),

                // B_q -> V gamma
                make_b_to_v_gamma_group(),

                // B_q -> M l^+ l^-
                make_b_to_p_ll_group(),
                make_b_to_v_ll_group(),

                // Lambda_b -> Lambda l^+ l^-
                make_lambdab_to_lambda_ll_group(),

                // Lambda_b -> Lambda l^+ l^-
                make_lambdab_to_lambda1520_ll_group(),

                // Lambda_b -> Lambda gamma
                make_lambdab_to_lambda1520_gamma_group(),

                // B -> X_s {gamma, l^+ l^-}
                make_b_to_xs_group(),

                // B_{u,d,s} -> P nu nubar
                make_b_to_k_nu_nu_group(),
                make_bs_to_eta_nu_nu_group(),

                // B_{u,d} -> V nu nubar
                make_b_to_kstar_nu_nu_group(),
                make_bs_to_phi_nu_nu_group(),

                // Lambda_b -> Lambda nu nubar
                make_lambda_b_to_lambda_nu_nu_group()
            }
        );

        return ObservableSection(imp);
    }
}
