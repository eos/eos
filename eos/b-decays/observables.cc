/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

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

#include <eos/observable-impl.hh>
#include <eos/b-decays/b-to-d-pi-l-nu.hh>
#include <eos/b-decays/b-to-l-nu.hh>
#include <eos/b-decays/b-to-pi-pi-l-nu.hh>
#include <eos/b-decays/b-to-psd-l-nu.hh>
#include <eos/b-decays/b-to-vec-l-nu.hh>
#include <eos/b-decays/bs-to-kstar-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac2595-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac2625-l-nu.hh>
#include <eos/b-decays/inclusive-b-to-u.hh>
#include <eos/b-decays/properties.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Leptonic B decays
    // {{{
    ObservableGroup
    make_b_to_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B^-\to \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour.)",
            {
                make_observable("B_u->lnu::BR", R"(\mathcal{B}(B^- \to \ell^-\bar\nu))",
                        &BToLeptonNeutrino::branching_ratio),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Semileptonic B -> P(seudoscalar) decays
    // {{{

    // B -> pi l nu
    // {{{
    ObservableGroup
    make_b_to_pi_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \pi \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->pilnu::dBR/dq2", R"(d\mathcal{B}(B\to\pi\ell^-\bar\nu)/dq^2)",
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" } }),

                make_observable("B->pilnu::BR", R"(\mathcal{B}(B\to\pi\ell^-\bar\nu))",
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" } }),

                make_observable_ratio("B->pilnu::R_pi(q2)", R"(R_{\pi}(q^2))",
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, { "l", "tau" } },
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, { "l", "mu" } }),

                make_observable_ratio("B->pilnu::R_pi", R"(R_{\pi})",
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_tau_min", "q2_tau_max"),
                        Options{ { "U", "u" }, { "l", "tau"} },
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_tau_min", "q2_tau_max"),
                        Options{ { "U", "u" }, { "l", "mu"} }),

                make_observable("B->pilnu::P(q2)", R"(dP(B\to\pi\ell^-\bar\nu)/dq^2)",
                        &BToPseudoscalarLeptonNeutrino::differential_pdf_q2,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" } }),

                make_observable("B->pilnu::P(q2_min,q2_max)", R"(P(B\to\pi\ell^-\bar\nu))",
                        &BToPseudoscalarLeptonNeutrino::integrated_pdf_q2,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" } }),

                make_observable("B->pilnu::zeta",
                        &BToPseudoscalarLeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> D l nu
    // {{{
    ObservableGroup
    make_b_to_d_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \bar{D} \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->Dlnu::dBR/dq2", R"(d\mathcal{B}(B\to \bar{D}\ell^-\bar\nu)/dq^2)",
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" } }),

                make_observable("B->Dlnu::BR", R"(\mathcal{B}(B\to \bar{D}\ell^-\bar\nu))",
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" } }),

                make_observable("B->Dlnu::normdBR/ds",
                        &BToPseudoscalarLeptonNeutrino::normalized_differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" } }),

                make_observable("B->Dlnu::normBR",
                        &BToPseudoscalarLeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" } }),

                make_observable_ratio("B->Dlnu::R_D(q2)", R"(R_D(q^2))",
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, { "l", "tau" } },
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, { "l", "mu" } }),

                make_observable_ratio("B->Dlnu::R_D", R"(R_D)",
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_tau_min", "q2_tau_max"),
                        Options{ { "U", "c" }, { "l", "tau"} },
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_mu_min", "q2_mu_max"),
                        Options{ { "U", "c" }, { "l", "mu"} }),

                make_observable("B->Dlnu::A_FB(q2)", R"(A_{\text{FB}}(B\to \bar{D}\ell^-\bar\nu)(q^2))",
                        &BToPseudoscalarLeptonNeutrino::differential_a_fb_leptonic,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" } }),

                make_observable("B->Dlnu::A_FB", R"(A_{\text{FB}}(B\to \bar{D}\ell^-\bar\nu))",
                        &BToPseudoscalarLeptonNeutrino::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" } }),

                make_observable("B->Dlnu::P(w)",
                        &BToPseudoscalarLeptonNeutrino::differential_pdf_w,
                        std::make_tuple("w"),
                        Options{ { "U", "c" } }),

                make_observable("B->Dlnu::P(w_min,w_max)",
                        &BToPseudoscalarLeptonNeutrino::integrated_pdf_w,
                        std::make_tuple("w_min", "w_max"),
                        Options{ { "U", "c" } }),

                make_observable("B->Dlnu::A_l",
                        &BToPseudoscalarLeptonNeutrino::integrated_lepton_polarization,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> D_s l nu
    // {{{
    ObservableGroup
    make_bs_to_ds_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_s\to \bar{D_s} \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour.)"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B_s->D_slnu::dBR/dq2", R"(d\mathcal{B}(B_s\to \bar{D}_s\ell^-\bar\nu)/dq^2)",
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_slnu::BR", R"(\mathcal{B}(B_s\to \bar{D}_s\ell^-\bar\nu))",
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_slnu::normdBR/ds",
                        &BToPseudoscalarLeptonNeutrino::normalized_differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_slnu::normBR",
                        &BToPseudoscalarLeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, {"q", "s"} }),

                make_observable_ratio("B_s->D_slnu::R_D_s(q2)", R"(R_D_s(q^2))",
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, {"q", "s"}, { "l", "tau" } },
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, {"q", "s"}, { "l", "mu" } }),

                make_observable_ratio("B_s->D_slnu::R_D_s", R"(R_D_s)",
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_tau_min", "q2_tau_max"),
                        Options{ { "U", "c" }, {"q", "s"}, { "l", "tau"} },
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_mu_min", "q2_mu_max"),
                        Options{ { "U", "c" }, {"q", "s"}, { "l", "mu"} }),

                make_observable("B_s->D_slnu::A_FB(q2)", R"(A_{\text{FB}}(B_s\to \bar{D}_s\ell^-\bar\nu)(q^2))",
                        &BToPseudoscalarLeptonNeutrino::differential_a_fb_leptonic,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_slnu::A_FB", R"(A_{\text{FB}}(B_s\to \bar{D}_s\ell^-\bar\nu))",
                        &BToPseudoscalarLeptonNeutrino::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_slnu::P(w)",
                        &BToPseudoscalarLeptonNeutrino::differential_pdf_w,
                        std::make_tuple("w"),
                        Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_slnu::P(w_min,w_max)",
                        &BToPseudoscalarLeptonNeutrino::integrated_pdf_w,
                        std::make_tuple("w_min", "w_max"),
                        Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_slnu::A_l",
                        &BToPseudoscalarLeptonNeutrino::integrated_lepton_polarization,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, {"q", "s"} }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // Semileptonic B -> V(seudoscalar) decays
    // {{{

    // B -> D^* l nu
    // {{{
    ObservableGroup
    make_b_to_dstar_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \bar{D}^* \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // B -> D^* l nu
                make_observable("B->D^*lnu::dBR/dq2", R"(d\mathcal{B}(B\to \bar{D}^*\ell^-\bar\nu)/dq^2)",
                                &BToVectorLeptonNeutrino::differential_branching_ratio,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::normdBR/dq2",
                                &BToVectorLeptonNeutrino::normalized_differential_branching_ratio,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::A_FB(q2)", R"(A_{\text{FB}}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_a_fb_leptonic,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_1c(q2)", R"(J_{1c}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J1c_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_1s(q2)", R"(J_{1s}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J1s_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_2c(q2)", R"(J_{2c}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J2c_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_2s(q2)", R"(J_{2s}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J2s_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_3(q2)", R"(J_{3}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J3_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_4(q2)", R"(J_{4}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J4_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_5(q2)", R"(J_{5}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J5_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_6c(q2)", R"(J_{6c}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J6c_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_6s(q2)", R"(J_{6s}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J6s_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_7(q2)", R"(J_{7}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J7_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_8(q2)", R"(J_{8}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J8_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::J_9(q2)", R"(J_{9}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J9_normalized,
                                std::make_tuple("q2")),

                make_observable("B->D^*lnu::BR", R"(\mathcal{B}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_branching_ratio,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::normBR",
                                &BToVectorLeptonNeutrino::normalized_integrated_branching_ratio,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable_ratio("B->D^*lnu::R_D^*(q2)", R"(R_{D^*}(q^2))",
                                &BToVectorLeptonNeutrino::differential_branching_ratio,
                                std::make_tuple("q2"),
                                Options{ { "l", "tau" } },
                                &BToVectorLeptonNeutrino::differential_branching_ratio,
                                std::make_tuple("q2"),
                                Options{ { "l", "mu" } }),

                make_observable_ratio("B->D^*lnu::R_D^*", R"(R_{D^*})",
                                &BToVectorLeptonNeutrino::integrated_branching_ratio,
                                std::make_tuple("q2_tau_min", "q2_tau_max"),
                                Options{ { "l", "tau" } },
                                &BToVectorLeptonNeutrino::integrated_branching_ratio,
                                std::make_tuple("q2_mu_min", "q2_mu_max"),
                                Options{ { "l", "mu" } }),

                make_observable("B->D^*lnu::A_FB", R"(A_{\text{FB}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_fb_leptonic,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_L",
                                &BToVectorLeptonNeutrino::integrated_amplitude_polarization_L,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_T",
                                &BToVectorLeptonNeutrino::integrated_amplitude_polarization_T,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::F_L", R"(F_{\text{L}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_f_L,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::Ftilde_L", R"(\tilde{F}_{\text{L}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_ftilde_L,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_C^1", R"(A_{\text{C}}^1(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_c_1,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_C^2", R"(A_{\text{C}}^2(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_c_2,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_C^3", R"(A_{\text{C}}^3(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_c_3,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_T^1", R"(A_{\text{T}}^1(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_t_1,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_T^2", R"(A_{\text{T}}^2(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_t_2,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::A_T^3", R"(A_{\text{T}}^3(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_t_3,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_1c", R"(J_{1c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J1c_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_1s", R"(J_{1s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J1s_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_2c", R"(J_{2c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J2c_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_2s", R"(J_{2s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J2s_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_3", R"(J_{3}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J3_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_4", R"(J_{4}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J4_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_5", R"(J_{5}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J5_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_6c", R"(J_{6c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J6c_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_6s", R"(J_{6s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J6s_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_7", R"(J_{7}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J7_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_8", R"(J_{8}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J8_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::J_9", R"(J_{9}(B\to \bar{D}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J9_normalized,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->D^*lnu::P(w_min,w_max)",
                                &BToVectorLeptonNeutrino::integrated_pdf_w,
                                std::make_tuple("w_min", "w_max"),
                                Options{ { "U", "c" }, {"q", "d"} }),

                // B -> D pi l nu
                make_observable("B->Dpilnu::P(c_D)",
                                &BToDPiLeptonNeutrino::differential_pdf_d,
                                std::make_tuple("c_D")),

                make_observable("B->Dpilnu::P(c_l)",
                                &BToDPiLeptonNeutrino::differential_pdf_l,
                                std::make_tuple("c_l")),

                make_observable("B->Dpilnu::P(chi)",
                                &BToDPiLeptonNeutrino::differential_pdf_chi,
                                std::make_tuple("chi")),

                make_observable("B->Dpilnu::P(w)",
                                &BToDPiLeptonNeutrino::differential_pdf_w,
                                std::make_tuple("w")),

                make_observable("B->Dpilnu::P(q2)",
                                &BToDPiLeptonNeutrino::differential_pdf_q2,
                                std::make_tuple("q2")),

                make_observable("B->Dpilnu::A_l",
                                &BToDPiLeptonNeutrino::integrated_lepton_polarization,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Dpilnu::P(c_D_min,c_D_max)",
                                &BToDPiLeptonNeutrino::integrated_pdf_d,
                                std::make_tuple("c_D_min", "c_D_max")),

                make_observable("B->Dpilnu::P(c_l_min,c_l_max)",
                                &BToDPiLeptonNeutrino::integrated_pdf_l,
                                std::make_tuple("c_l_min", "c_l_max")),

                make_observable("B->Dpilnu::P(chi_min,chi_max)",
                                &BToDPiLeptonNeutrino::integrated_pdf_chi,
                                std::make_tuple("chi_min", "chi_max")),

                make_observable("B->Dpilnu::P(w_min,w_max)",
                                &BToDPiLeptonNeutrino::integrated_pdf_w,
                                std::make_tuple("w_min", "w_max")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> D_s^* l nu
    // {{{
    ObservableGroup
    make_bs_to_dsstar_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\bar{B}_s\to D_s^* \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour.)"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // B_s -> D_s^* l nu
                make_observable("B_s->D_s^*lnu::dBR/dq2", R"(d\mathcal{B}(B_s\to \bar{D}_s^*\ell^-\bar\nu)/dq^2)",
                                &BToVectorLeptonNeutrino::differential_branching_ratio,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::normdBR/dq2",
                                &BToVectorLeptonNeutrino::normalized_differential_branching_ratio,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::A_FB(q2)", R"(A_{\text{FB}}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_a_fb_leptonic,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_1c(q2)", R"(J_{1c}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J1c_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_1s(q2)", R"(J_{1s}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J1s_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_2c(q2)", R"(J_{2c}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J2c_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_2s(q2)", R"(J_{2s}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J2s_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_3(q2)", R"(J_{3}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J3_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_4(q2)", R"(J_{4}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J4_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_5(q2)", R"(J_{5}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J5_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_6c(q2)", R"(J_{6c}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J6c_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_6s(q2)", R"(J_{6s}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J6s_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_7(q2)", R"(J_{7}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J7_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_8(q2)", R"(J_{8}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J8_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_9(q2)", R"(J_{9}(B_s\to \bar{D_s}^*\ell^-\bar\nu)(q^2))",
                                &BToVectorLeptonNeutrino::differential_J9_normalized,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::BR", R"(\mathcal{B}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_branching_ratio,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::normBR",
                                &BToVectorLeptonNeutrino::normalized_integrated_branching_ratio,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable_ratio("B_s->D_s^*lnu::R_D_s^*(q2)", R"(R_{D_s^*}(q^2))",
                                &BToVectorLeptonNeutrino::differential_branching_ratio,
                                std::make_tuple("q2"),
                                Options{ { "l", "tau" }, { "U", "c" }, {"q", "s"} },
                                &BToVectorLeptonNeutrino::differential_branching_ratio,
                                std::make_tuple("q2"),
                                Options{ { "l", "mu" }, { "U", "c" }, {"q", "s"} }),

                make_observable_ratio("B_s->D_s^*lnu::R_D_s^*", R"(R_{D_s^*})",
                                &BToVectorLeptonNeutrino::integrated_branching_ratio,
                                std::make_tuple("q2_tau_min", "q2_tau_max"),
                                Options{ { "l", "tau" }, { "U", "c" }, {"q", "s"} },
                                &BToVectorLeptonNeutrino::integrated_branching_ratio,
                                std::make_tuple("q2_mu_min", "q2_mu_max"),
                                Options{ { "l", "mu" }, { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::A_FB", R"(A_{\text{FB}}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_fb_leptonic,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::A_L",
                                &BToVectorLeptonNeutrino::integrated_amplitude_polarization_L,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::A_T",
                                &BToVectorLeptonNeutrino::integrated_amplitude_polarization_T,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::F_L", R"(F_{\text{L}}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_f_L,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::A_C^1", R"(A_{\text{C}}^1(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_c_1,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::A_C^2", R"(A_{\text{C}}^2(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_c_2,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::A_C^3", R"(A_{\text{C}}^3(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_c_3,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::A_T^1", R"(A_{\text{T}}^1(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_t_1,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::A_T^2", R"(A_{\text{T}}^2(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_t_2,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::A_T^3", R"(A_{\text{T}}^3(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_a_t_3,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_1c", R"(J_{1c}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J1c_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_1s", R"(J_{1s}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J1s_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_2c", R"(J_{2c}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J2c_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_2s", R"(J_{2s}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J2s_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_3", R"(J_{3}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J3_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_4", R"(J_{4}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J4_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_5", R"(J_{5}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J5_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_6c", R"(J_{6c}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J6c_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_6s", R"(J_{6s}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J6s_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_7", R"(J_{7}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J7_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_8", R"(J_{8}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J8_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::J_9", R"(J_{9}(B_s\to \bar{D_s}^*\ell^-\bar\nu))",
                                &BToVectorLeptonNeutrino::integrated_J9_normalized,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),

                make_observable("B_s->D_s^*lnu::P(w_min,w_max)",
                                &BToVectorLeptonNeutrino::integrated_pdf_w,
                                std::make_tuple("w_min", "w_max"),
                                Options{ { "U", "c" }, {"q", "s"} }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> K^* l nu
    // {{{
    ObservableGroup
    make_bs_to_kstar_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_s\to \bar{K}^* \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // B_s -> K^* l nubar
                make_observable("B_s->K^*lnu::F_perp(q2)",
                        &BsToKstarLeptonNeutrino::Fperp,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::F_para(q2)",
                        &BsToKstarLeptonNeutrino::Fpara,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::F_long(q2)",
                        &BsToKstarLeptonNeutrino::Flong,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::d^4Gamma",
                        &BsToKstarLeptonNeutrino::four_differential_decay_width,
                        std::make_tuple("q2", "cos(theta_l)", "cos(theta_k)", "phi")),

                make_observable("B_s->K^*lnu::dBR/ds", R"(d\mathcal{B}(B_s\to \bar{K}^*\ell^-\bar\nu)/dq^2)",
                        &BsToKstarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_FB(q2)",
                        &BsToKstarLeptonNeutrino::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::BR",
                        &BsToKstarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::A_FB",
                        &BsToKstarLeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_1s",
                        &BsToKstarLeptonNeutrino::integrated_s_1s,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_1c",
                        &BsToKstarLeptonNeutrino::integrated_s_1c,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_2s",
                        &BsToKstarLeptonNeutrino::integrated_s_2s,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_2c",
                        &BsToKstarLeptonNeutrino::integrated_s_2c,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_3",
                        &BsToKstarLeptonNeutrino::integrated_s_3,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_4",
                        &BsToKstarLeptonNeutrino::integrated_s_4,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_5",
                        &BsToKstarLeptonNeutrino::integrated_s_5,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::Shat_6s",
                        &BsToKstarLeptonNeutrino::integrated_s_6s,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B_s->K^*lnu::A_FB(q2)",
                        &BsToKstarLeptonNeutrino::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^2(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_2,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^3(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_3,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^4(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_4,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^5(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_5,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^re(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_re,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_T^im(q2)",
                        &BsToKstarLeptonNeutrino::differential_transverse_asymmetry_im,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::F_L(q2)",
                        &BsToKstarLeptonNeutrino::differential_longitudinal_polarisation,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::F_T(q2)",
                        &BsToKstarLeptonNeutrino::differential_transversal_polarisation,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::H_T^1(q2)",
                        &BsToKstarLeptonNeutrino::differential_h_1,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::H_T^2(q2)",
                        &BsToKstarLeptonNeutrino::differential_h_2,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::H_T^3(q2)",
                        &BsToKstarLeptonNeutrino::differential_h_3,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::H_T^4(q2)",
                        &BsToKstarLeptonNeutrino::differential_h_4,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::H_T^5(q2)",
                        &BsToKstarLeptonNeutrino::differential_h_5,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*lnu::A_FB",
                        &BsToKstarLeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::BR",
                        &BsToKstarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::F_L",
                        &BsToKstarLeptonNeutrino::integrated_longitudinal_polarisation,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::F_T",
                        &BsToKstarLeptonNeutrino::integrated_transversal_polarisation,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^2",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_2,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^3",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^4",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^5",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^re",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_re,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::A_T^im",
                        &BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_im,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::H_T^1",
                        &BsToKstarLeptonNeutrino::integrated_h_1,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::H_T^2",
                        &BsToKstarLeptonNeutrino::integrated_h_2,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::H_T^3",
                        &BsToKstarLeptonNeutrino::integrated_h_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::H_T^4",
                        &BsToKstarLeptonNeutrino::integrated_h_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B_s->K^*lnu::H_T^5",
                        &BsToKstarLeptonNeutrino::integrated_h_5,
                        std::make_tuple("q2_min", "q2_max")),

                // B_s -> K^* l nubar Ratios
                make_observable("B_s->K^*lnu::R_long",
                        &BsToKstarLeptonNeutrinoRatios::ratio_long),

                make_observable("B_s->K^*lnu::R_para",
                        &BsToKstarLeptonNeutrinoRatios::ratio_para),

                make_observable("B_s->K^*lnu::R_perp",
                        &BsToKstarLeptonNeutrinoRatios::ratio_perp),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // Semileptonic B -> P(seudoscalar) P(seudoscalar) decays
    // {{{

    // B -> pi pi l nu
    // {{{
    ObservableGroup
    make_b_to_pi_pi_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \pi\pi \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->pipilnu::BR(q2,k2)", R"(d^2\mathcal{B}(B\to \pi\pi \ell^-\bar\nu)/(dq^2\,dk^2))",
                        &BToPiPiLeptonNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "k2")),

                make_observable("B->pipilnu::BR(q2,k2,cos(theta_pi))",
                        &BToPiPiLeptonNeutrino::triple_differential_branching_ratio,
                        std::make_tuple("q2", "k2", "cos(theta_pi)")),

                make_observable("B->pipilnu::A_FB(q2,k2)", R"(A_{\text{FB}}(B\to \pi\pi \ell^-\bar\nu)(q^2,k^2))",
                        &BToPiPiLeptonNeutrino::double_differential_forward_backward_asymmetry,
                        std::make_tuple("q2", "k2")),

                make_observable("B->pipilnu::P(cos(theta_pi))",
                        &BToPiPiLeptonNeutrino::partial_waves,
                        std::make_tuple("q2", "k2", "cos(theta_pi)")),

                make_observable("B->pipilnu::BR", R"(\mathcal{B}(B\to \pi\pi \ell^-\bar\nu))",
                        &BToPiPiLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max", "z_min", "z_max")),

                make_observable("B->pipilnu::A_FB", R"(A_{\text{FB}}(B\to \pi\pi \ell^-\bar\nu))",
                        &BToPiPiLeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // Semileptonic Lambda_b decays
    // {{{

    // Lambda_b -> Lambda_c l nu
    // {{{
    ObservableGroup
    make_lambdab_to_lambdac_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\Lambda_b\to \Lambda_c \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavour. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // Lambda_b -> Lambda_c l nu
                make_observable("Lambda_b->Lambda_clnu::dBR/dq2", R"(d\mathcal{B}(\Lambda_b\to\Lambda_c \ell^-\bar\nu)/dq^2)",
                        &LambdaBToLambdaCLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^l(q2)", R"(A_{\text{FB}}^\ell(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_leptonic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^h(q2)", R"(A_{\text{FB}}^h(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_hadronic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^c(q2)", R"(A_{\text{FB}}^{h\ell}(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_combined,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::F_0(q2)", R"(F_0(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        &LambdaBToLambdaCLeptonNeutrino::differential_fzero,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::BR", R"(\mathcal{B}(\Lambda_b\to\Lambda_c \ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable_ratio("Lambda_b->Lambda_clnu::R(Lambda_c)", R"(R(\Lambda_c))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_tau_min", "q2_tau_max"),
                        Options{ { "l", "tau" } },
                        &LambdaBToLambdaCLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_mu_min", "q2_mu_max"),
                        Options{ { "l", "mu" } }),

                make_observable("Lambda_b->Lambda_clnu::A_FB^l",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^h",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_hadronic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable_ratio("Lambda_b->Lambda_clnu::R(A_FB^h)(q2)", R"(R(A_{\text{FB}}^{\Lambda_c})(q^2))",
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_hadronic,
                        std::make_tuple("q2"),
                        Options{ { "l", "tau" } },
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_hadronic,
                        std::make_tuple("q2"),
                        Options{ { "l", "mu" } }),

                make_observable_ratio("Lambda_b->Lambda_clnu::R(A_FB^h)", R"(R(A_{\text{FB}}^{\Lambda_c}))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_hadronic,
                        std::make_tuple("q2_tau_min", "q2_tau_max"),
                        Options{ { "l", "tau" } },
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_hadronic,
                        std::make_tuple("q2_mu_min", "q2_mu_max"),
                        Options{ { "l", "mu" } }),

                make_observable("Lambda_b->Lambda_clnu::A_FB^c",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_combined,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::F_0",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_fzero,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_1ss", R"(K_{1ss}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k1ss,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_1cc", R"(K_{1cc}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k1cc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_1c", R"(K_{1c}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k1c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_2ss", R"(K_{2ss}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k2ss,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_2cc", R"(K_{2cc}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k2cc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_2c", R"(K_{2c}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k2c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_3sc", R"(K_{3sc}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k3sc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_3s", R"(K_{3s}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k3s,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_4sc", R"(K_{4sc}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k4sc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_4s", R"(K_{4s}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k4s,
                        std::make_tuple("q2_min", "q2_max")),

                // Lambda_b -> Lambda_c(2595) l nubar
                make_observable("Lambda_b->Lambda_c(2595)lnu::dBR/ds",
                        &LambdaBToLambdaC2595LeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::dBR/dsdtheta_l",
                        &LambdaBToLambdaC2595LeptonNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "theta_l")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::BR", R"(\mathcal{B}(\Lambda_b\to\Lambda_c(2595) \ell^-\bar\nu))",
                        &LambdaBToLambdaC2595LeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::A_FB",
                        &LambdaBToLambdaC2595LeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::Gamma_normalized(q2_min,q2_max)",
                        &LambdaBToLambdaC2595LeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable_ratio("Lambda_b->Lambda_c(2595)lnu::R_Lambda_c(2595)(q2)", R"(R_{\Lambda_c(2595}(q^2))",
                        &LambdaBToLambdaC2595LeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "l", "tau" } },
                        &LambdaBToLambdaC2595LeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "l", "mu" } }),

                make_observable_ratio("Lambda_b->Lambda_c(2595)lnu::R_Lambda_c(2595)", R"(R_{\Lambda_c(2595})",
                        &LambdaBToLambdaC2595LeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_tau_min", "q2_tau_max"),
                        Options{ { "l", "tau" } },
                        &LambdaBToLambdaC2595LeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_mu_min", "q2_mu_max"),
                        Options{ { "l", "mu" } }),

                // Lambda_b -> Lambda_c(2625) l nubar
                make_observable("Lambda_b->Lambda_c(2625)lnu::dBR/ds",
                        &LambdaBToLambdaC2625LeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::A_FB(q2)",
                        &LambdaBToLambdaC2625LeptonNeutrino::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::dBR/dsdtheta_l",
                        &LambdaBToLambdaC2625LeptonNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "theta_l")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::BR", R"(\mathcal{B}(\Lambda_b\to\Lambda_c(2625) \ell^-\bar\nu))",
                        &LambdaBToLambdaC2625LeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::A_FB",
                        &LambdaBToLambdaC2625LeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::Gamma_normalized(q2_min,q2_max)",
                        &LambdaBToLambdaC2625LeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable_ratio("Lambda_b->Lambda_c(2625)lnu::R_Lambda_c(2625)(q2)", R"(R_{\Lambda_c(2625}(q^2))",
                        &LambdaBToLambdaC2625LeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "l", "tau" } },
                        &LambdaBToLambdaC2625LeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "l", "mu" } }),

                make_observable_ratio("Lambda_b->Lambda_c(2625)lnu::R_Lambda_c(2625)", R"(R_{\Lambda_c(2625})",
                        &LambdaBToLambdaC2625LeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_tau_min", "q2_tau_max"),
                        Options{ { "l", "tau" } },
                        &LambdaBToLambdaC2625LeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_mu_min", "q2_mu_max"),
                        Options{ { "l", "mu" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // Misc.
    // {{{
    ObservableGroup
    make_b_to_xu_semileptonic_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Miscellaneous matrix elements)",
            R"()",
            {
                /* B Meson Properties */
                make_observable("B::M_B^*-M_B",
                        &BMesonProperties::mass_splitting_j1_j0),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_b_decays_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in (semi)leptonic $b$-hadron decays",
            "",
            {
                // B^- -> l^- nubar
                make_b_to_l_nu_group(),

                // B_{u,d} -> P l^- nubar
                make_b_to_pi_l_nu_group(),
                make_b_to_d_l_nu_group(),

                // B_s -> P l^- nubar
                make_bs_to_ds_l_nu_group(),

                // B_{u,d} -> V l^- nubar
                make_b_to_dstar_l_nu_group(),

                // B_s -> V l^- nubar
                make_bs_to_kstar_l_nu_group(),
                make_bs_to_dsstar_l_nu_group(),

                // Lambda_b
                make_lambdab_to_lambdac_l_nu_group(),

                // B -> X_u l^- nubar
                make_b_to_xu_semileptonic_group(),
            }
        );

        return ObservableSection(imp);
    }
}
