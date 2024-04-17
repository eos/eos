/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2019-2023 Danny van Dyk
 * Copyright (c) 2022 Philip Lüghausen
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
#include <eos/b-decays/b-to-gamma-l-nu.hh>
#include <eos/b-decays/b-to-l-nu.hh>
#include <eos/b-decays/b-to-3l-nu.hh>
#include <eos/b-decays/b-to-pi-pi-l-nu.hh>
#include <eos/b-decays/b-to-psd-l-nu.hh>
#include <eos/b-decays/b-to-psd-psd.hh>
#include <eos/b-decays/b-to-vec-l-nu.hh>
#include <eos/b-decays/b-to-vec-l-nu-impl.hh>
#include <eos/b-decays/bq-to-dq-psd.hh>
#include <eos/b-decays/bq-to-dstarq-psd.hh>
#include <eos/b-decays/lambdab-to-lambdac-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac2595-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac2625-l-nu.hh>
#include <eos/b-decays/lifetime.hh>
#include <eos/b-decays/inclusive-b-to-u.hh>
#include <eos/b-decays/properties.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
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
            R"(The option "l" selects the charged lepton flavor.)",
            {
                make_observable("B_u->lnu::BR", R"(\mathcal{B}(B^- \to \ell^-\bar\nu))",
                        Unit::None(),
                        &BToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q", "u" } }),

                make_observable("B_c->lnu::BR", R"(\mathcal{B}(B_c^- \to \ell^-\bar\nu))",
                        Unit::None(),
                        &BToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q", "c" } }),
            }
        );

        return ObservableGroup(imp);
    }

    ObservableGroup
    make_b_to_3l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B^-\to \ell^-\bar\nu\ell'^+\ell'^-$ decays)",
            R"(The option "l" selects the charged lepton flavour coming out of the )"
            R"(weak current, "lprime" selects the lepton flavour coming out of the photon.)",
            {
                make_observable("B_u->enumumu::d2BR/dq2/dk2", R"(\frac{d\mathcal{B}(B^- \to e^-\bar\nu\mu^+\mu^-)}{dq^2/dk^2})",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "k2"),
                        Options{ { "l", "e" }, { "lprime", "mu" } }),

                make_observable("B_u->munuee::d2BR/dq2/dk2", R"(\frac{d\mathcal{B}(B^- \to \mu^-\bar\nu e^+e^-)}{dq^2/dk^2})",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "k2"),
                        Options{ { "l", "mu" }, { "lprime", "e" } }),

                make_observable("B_u->taunuee::d2BR/dq2/dk2", R"(\frac{d\mathcal{B}(B^- \to \tau^-\bar\nu e^+e^-)}{dq^2/dk^2})",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "k2"),
                        Options{ { "l", "tau" }, { "lprime", "e" } }),

                make_observable("B_u->taunumumu::d2BR/dq2/dk2", R"(\frac{d\mathcal{B}(B^- \to \tau^-\bar\nu\mu^+\mu^-)}{dq^2/dk^2})",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "k2"),
                        Options{ { "l", "tau" }, { "lprime", "mu" } }),

                make_observable("B_u->enumumu::BR", R"(\mathcal{B}(B^- \to e^-\bar\nu\mu^+\mu^-))",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max"),
                        Options{ { "l", "e" }, { "lprime", "mu" } }),

                make_observable("B_u->munuee::BR", R"(\mathcal{B}(B^- \to \mu^-\bar\nu e^+e^-))",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max"),
                        Options{ { "l", "mu" }, { "lprime", "e" } }),

                make_observable("B_u->taunuee::BR", R"(\mathcal{B}(B^- \to \tau^-\bar\nu e^+e^-))",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max"),
                        Options{ { "l", "tau" }, { "lprime", "e" } }),

                make_observable("B_u->taunumumu::BR", R"(\mathcal{B}(B^- \to \tau^-\bar\nu\mu^+\mu^-))",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max"),
                        Options{ { "l", "tau" }, { "lprime", "mu" } }),

                make_observable("B_u->enumumu::A_FB", R"(A_{\mathrm{FB}}(B^- \to e^-\bar\nu\mu^+\mu^-))",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max"),
                        Options{ { "l", "e" }, { "lprime", "mu" } }),

                make_observable("B_u->munuee::A_FB", R"(A_{\mathrm{FB}}(B^- \to \mu^-\bar\nu e^+e^-))",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max"),
                        Options{ { "l", "mu" }, { "lprime", "e" } }),

                make_observable("B_u->taunumumu::A_FB", R"(A_{\mathrm{FB}}(B^- \to \tau^-\bar\nu\mu^+\mu^-))",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max"),
                        Options{ { "l", "tau" }, { "lprime", "mu" } }),

                make_observable("B_u->taunuee::A_FB", R"(A_{\mathrm{FB}}(B^- \to \tau^-\bar\nu e^+e^-))",
                        Unit::None(),
                        &BToThreeLeptonsNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max"),
                        Options{ { "l", "tau" }, { "lprime", "e" } }),
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
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->pilnu::dBR/dq2", R"(d\mathcal{B}(B\to\pi\ell^-\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::d^2BR/dq2/dcos(theta_l)", R"(d^2\mathcal{B}(B\to\pi\ell^-\bar\nu)/dq^2/d\cos(\theta_l))",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarLeptonNeutrino::two_differential_branching_ratio,
                        std::make_tuple("q2", "cos(theta_l)"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::BR", R"(\mathcal{B}(B\to\pi\ell^-\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::width", R"(\Gamma(B\to\pi\ell^-\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::width_p", R"(\Gamma(B\to\pi\ell^-\bar\nu)_p)",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width_p,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::width_0", R"(\Gamma(B\to\pi\ell^-\bar\nu)_0)",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width_0,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_expression_observable("B->pilnu::R_pi(q2)", R"(R_{\pi}(q^2))",
                        Unit::None(),
                        R"(
                        <<B->pilnu::dBR/dq2;l=tau>>
                        /
                        <<B->pilnu::dBR/dq2;l=mu>>
                        )"),

                make_expression_observable("B->pilnu::R_pi_p", R"(R_{\pi, P})",
                        Unit::None(),
                        R"(
                        <<B->pilnu::width_p;l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                        /
                        <<B->pilnu::width;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        )"),

                make_expression_observable("B->pilnu::R_pi_0", R"(R_{\pi, 0})",
                        Unit::None(),
                        R"(
                        <<B->pilnu::width_0;l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                        /
                        <<B->pilnu::width;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        )"),

                make_expression_observable("B->pilnu::R_pi", R"(R_{\pi})",
                        Unit::None(),
                        R"(
                        <<B->pilnu::BR;l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                        /
                        <<B->pilnu::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                       )"),

                make_observable("B->pilnu::A_FB(q2)", R"(A_{\mathrm{FB}}(B\to \pi\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::differential_a_fb_leptonic,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::A_FB", R"(A_{\mathrm{FB}}(B\to \pi\ell^-\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::P(q2)", R"(dP(B\to\pi\ell^-\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarLeptonNeutrino::differential_pdf_q2,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::P(q2_min,q2_max)", R"(P(B\to\pi\ell^-\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_pdf_q2,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::A_l",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_lepton_polarization,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::F_H",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_flat_term,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, { "I", "1" } }),

                make_observable("B->pilnu::zeta",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, { "I", "1" } }),
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
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->Dlnu::dBR/dq2", R"(d\mathcal{B}(B\to \bar{D}\ell^-\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->Dlnu::d^2BR/dq2/dcos(theta_l)", R"(d^2\mathcal{B}(B\to \bar{D}\ell^-\bar\nu)/dq^2/d\cos(\theta_l))",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarLeptonNeutrino::two_differential_branching_ratio,
                        std::make_tuple("q2", "cos(theta_l)"),
                        Options{ { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->Dlnu::BR", R"(\mathcal{B}(B\to \bar{D}\ell^-\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->Dlnu::normdBR/ds",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarLeptonNeutrino::normalized_differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->Dlnu::normBR",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, { "I", "1/2" } }),

                make_expression_observable("B->Dlnu::R_D(q2)", R"(R_D(q^2))",
                        Unit::None(),
                        R"(
                        <<B->Dlnu::dBR/dq2;l=tau>>
                        /
                        <<B->Dlnu::dBR/dq2;l=mu>>
                        )"),

                make_expression_observable("B->Dlnu::R_D", R"(R_D)",
                        Unit::None(),
                        R"(
                        <<B->Dlnu::BR;l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                        /
                        <<B->Dlnu::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        )"),

                make_observable("B->Dlnu::A_FB(q2)", R"(A_{\mathrm{FB}}(B\to \bar{D}\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::differential_a_fb_leptonic,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->Dlnu::A_FB", R"(A_{\mathrm{FB}}(B\to \bar{D}\ell^-\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->Dlnu::P(w)",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::differential_pdf_w,
                        std::make_tuple("w"),
                        Options{ { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->Dlnu::P(w_min,w_max)",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_pdf_w,
                        std::make_tuple("w_min", "w_max"),
                        Options{ { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->Dlnu::A_l",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_lepton_polarization,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, { "I", "1/2" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> K l nu
    // {{{
    ObservableGroup
    make_bs_to_k_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_s\to \bar{K} \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavor.)"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B_s->Klnu::dBR/dq2", R"(d\mathcal{B}(\bar{B}_s\to K\ell^-\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, { "I", "1/2" } }),

                make_observable("B_s->Klnu::BR", R"(\mathcal{B}(\bar{B}_s\to K\ell^-\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, {"q", "s"}, { "I", "1/2" } }),

                make_observable("B_s->Klnu::normdBR/ds",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarLeptonNeutrino::normalized_differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, { "I", "1/2" } }),

                make_observable("B_s->Klnu::normBR",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, {"q", "s"}, { "I", "1/2" } }),
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
            R"(The option "l" selects the charged lepton flavor.)"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B_s->D_slnu::dBR/dq2", R"(d\mathcal{B}(B_s\to \bar{D}_s\ell^-\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, {"q", "s"}, { "I", "0" } }),

                make_observable("B_s->D_slnu::BR", R"(\mathcal{B}(B_s\to \bar{D}_s\ell^-\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, {"q", "s"}, { "I", "0" } }),

                make_observable("B_s->D_slnu::normdBR/ds",
                        Unit::InverseGeV2(),
                        &BToPseudoscalarLeptonNeutrino::normalized_differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, {"q", "s"}, { "I", "0" } }),

                make_observable("B_s->D_slnu::normBR",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, {"q", "s"}, { "I", "0" } }),

                make_expression_observable("B_s->D_slnu::R_D_s(q2)", R"(R_{D_s}(q^2))",
                        Unit::None(),
                        R"(
                        <<B_s->D_slnu::dBR/dq2;U=c,q=s,l=tau>>
                        /
                        <<B_s->D_slnu::dBR/dq2;U=c,q=s,l=mu>>
                        )"),

                make_expression_observable("B_s->D_slnu::R_D_s", R"(R_{D_s})",
                        Unit::None(),
                        R"(
                        <<B_s->D_slnu::BR;U=c,q=s,l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                        /
                        <<B_s->D_slnu::BR;U=c,q=s,l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        )"),

                make_observable("B_s->D_slnu::A_FB(q2)", R"(A_{\mathrm{FB}}(B_s\to \bar{D}_s\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::differential_a_fb_leptonic,
                        std::make_tuple("q2"),
                        Options{ { "U", "c" }, {"q", "s"}, { "I", "0" } }),

                make_observable("B_s->D_slnu::A_FB", R"(A_{\mathrm{FB}}(B_s\to \bar{D}_s\ell^-\bar\nu))",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, {"q", "s"}, { "I", "0" } }),

                make_observable("B_s->D_slnu::P(w)",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::differential_pdf_w,
                        std::make_tuple("w"),
                        Options{ { "U", "c" }, {"q", "s"}, { "I", "0" } }),

                make_observable("B_s->D_slnu::P(w_min,w_max)",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_pdf_w,
                        std::make_tuple("w_min", "w_max"),
                        Options{ { "U", "c" }, {"q", "s"}, { "I", "0" } }),

                make_observable("B_s->D_slnu::A_l",
                        Unit::None(),
                        &BToPseudoscalarLeptonNeutrino::integrated_lepton_polarization,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "c" }, {"q", "s"}, { "I", "0" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // Semileptonic B -> V(pseudoscalar) decays
    // {{{

    // B -> gamma l nu
    // {{{
    ObservableGroup
    make_b_to_gamma_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_u \to \gamma \ell \nu_{\ell}$ decays)",
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B_u->gammalnu::BR(E_gamma_min)", R"(\mathcal{B}(B^- \to \gamma \ell^-\bar\nu))",
                        Unit::None(),
                        &BToGammaLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("E_gamma_min")),

                make_observable("B_u->gammalnu::A_FB(E_gamma_min)", R"(A_{\mathrm{FB}}(B^- \to \gamma \ell^-\bar\nu))",
                        Unit::None(),
                        &BToGammaLeptonNeutrino::forward_backward_asymmetry,
                        std::make_tuple("E_gamma_min")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> omega l nu
    // {{{
    ObservableGroup
    make_b_to_omega_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \omega \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->omegalnu::dBR/dq2", R"(d\mathcal{B}(B\to\omega\ell^-\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToVectorLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, { "q", "u" }, {"I", "0"} }),

                make_observable("B->omegalnu::BR", R"(\mathcal{B}(B\to\omega\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, { "q", "u" }, {"I", "0"}  }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> rho l nu
    // {{{
    ObservableGroup
    make_b_to_rho_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \rho \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->rholnu::dBR/dq2", R"(d\mathcal{B}(B\to\rho\ell^-\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToVectorLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"I", "1"} }),

                make_observable("B->rholnu::BR", R"(\mathcal{B}(B\to\rho\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, {"I", "1"}  }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> D^* l nu
    // {{{
    ObservableGroup
    make_b_to_dstar_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B\to \bar{D}^* \ell^-\bar\nu$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // B -> D^* l nu

                // q^2 - differential

                make_observable("B->D^*lnu::dBR/dq2", R"(d\mathcal{B}(B\to \bar{D}^*\ell^-\bar\nu)/dq^2)",
                                Unit::InverseGeV2(),
                                &BToVectorLeptonNeutrino::differential_branching_ratio,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::normdBR/dq2",
                                Unit::InverseGeV2(),
                                &BToVectorLeptonNeutrino::normalized_differential_branching_ratio,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_expression_observable("B->D^*lnu::R_{D^*}^{tau/mu}(q2)", R"(R_{D^*}^{\tau/\mu}(q^2))",
                                Unit::None(),
                                R"(
                                <<B->D^*lnu::dBR/dq2;l=tau>>
                                /
                                <<B->D^*lnu::dBR/dq2;l=mu>>
                                )"),

                make_expression_observable("B->D^*lnu::R_{D^*}^{e/mu}(q2)", R"(R_{D^*}^{e/\mu}(q^2))",
                                Unit::None(),
                                R"(
                                <<B->D^*lnu::dBR/dq2;l=e>>
                                /
                                <<B->D^*lnu::dBR/dq2;l=mu>>
                                )"),

                make_observable("B->D^*lnu::A_FB(q2)", R"(A_{\mathrm{FB}}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_a_fb_leptonic,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_1c(q2)", R"(J_{1c}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J1c,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_1s(q2)", R"(J_{1s}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J1s,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_2c(q2)", R"(J_{2c}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J2c,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_2s(q2)", R"(J_{2s}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J2s,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_3(q2)", R"(J_{3}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J3,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_4(q2)", R"(J_{4}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J4,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_5(q2)", R"(J_{5}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J5,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_6c(q2)", R"(J_{6c}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J6c,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_6s(q2)", R"(J_{6s}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J6s,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_7(q2)", R"(J_{7}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J7,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_8(q2)", R"(J_{8}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J8,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_observable("B->D^*lnu::J_9(q2)", R"(J_{9}(B\to \bar{D}^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J9,
                                std::make_tuple("q2"),
                                { { "U", "c" }, { "I", "1/2" } }),

                // q^2 - integrated

                make_observable("B->D^*lnu::normGamma_CP_specific", R"(\Gamma(B\to \bar{D}^*\ell^-\bar\nu)_{|V_{cb}|=1})",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::normalized_decay_width,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_expression_observable("B->D^*lnu::normGamma", R"(\bar{\Gamma}(B\to \bar{D}^*\ell^-\bar\nu)_{|V_{cb}|=1})",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::normGamma_CP_specific;cp-conjugate=false>>
                                +
                                0.5 * <<B->D^*lnu::normGamma_CP_specific;cp-conjugate=true>>
                                )"),

                make_observable("B->D^*lnu::BR_CP_specific", R"(\mathcal{B}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::integrated_branching_ratio,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_expression_observable("B->D^*lnu::BR", R"(\bar{\mathcal{B}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::BR_CP_specific;cp-conjugate=false>>
                                +
                                0.5 * <<B->D^*lnu::BR_CP_specific;cp-conjugate=true>>
                                )"),

                make_expression_observable("B->D^*lnu::BRbar", R"(\mathcal{B}(B\to \bar{D}^*\ell^-\bar\nu)_{\ell=e,\mu})",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                +
                                0.5 * <<B->D^*lnu::BR;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                                )"),

                make_expression_observable("B->D^*lnu::DeltaBR", R"(\Delta\mathcal{B}(B\to \bar{D}^*\ell^-\bar\nu)_{\ell=e,\mu})",
                                Unit::None(),
                                R"(
                                <<B->D^*lnu::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                -
                                <<B->D^*lnu::BR;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                                )"),

                make_observable("B->D^*lnu::normBR",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::normalized_integrated_branching_ratio,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_expression_observable("B->D^*lnu::R_D^*", R"(R_{D^*})",
                                Unit::None(),
                                R"(
                                <<B->D^*lnu::BR;l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                                /
                                <<B->D^*lnu::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                )"),

                make_cacheable_observable("B->D^*lnu::A_L", R"()",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_amplitude_polarization_L,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::A_T", R"()",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_amplitude_polarization_T,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::A_C^1", R"(A_{\mathrm{C}}^1(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_c_1,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::A_C^2", R"(A_{\mathrm{C}}^2(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_c_2,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::A_C^3", R"(A_{\mathrm{C}}^3(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_c_3,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::A_T^1", R"(A_{\mathrm{T}}^1(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_t_1,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::A_T^2", R"(A_{\mathrm{T}}^2(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_t_2,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::A_T^3", R"(A_{\mathrm{T}}^3(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_t_3,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_1c", R"(J_{1c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J1c,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_1s", R"(J_{1s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J1s,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_2c", R"(J_{2c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J2c,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_2s", R"(J_{2s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J2s,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_3", R"(J_{3}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J3,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_4", R"(J_{4}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J4,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_5", R"(J_{5}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J5,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_6c", R"(J_{6c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J6c,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_6s", R"(J_{6s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J6s,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_7", R"(J_{7}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J7,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_8", R"(J_{8}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J8,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::J_9", R"(J_{9}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J9,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_cacheable_observable("B->D^*lnu::A_FB_CP_specific", R"(A_{\mathrm{FB}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_fb_leptonic,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_expression_observable("B->D^*lnu::A_FB", R"(A_{\mathrm{FB}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::A_FB_CP_specific;cp-conjugate=false>>
                                +
                                0.5 * <<B->D^*lnu::A_FB_CP_specific;cp-conjugate=true>>
                                )"),

                make_expression_observable("B->D^*lnu::Abar_FB", R"(\bar{A}_{\mathrm{FB}}(B\to \bar{D}^*\ell^-\bar\nu)_{\ell=e,\mu})",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::A_FB;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                +
                                0.5 * <<B->D^*lnu::A_FB;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                                )"),

                make_expression_observable("B->D^*lnu::DeltaA_FB", R"()",
                                Unit::None(),
                                R"(
                                <<B->D^*lnu::A_FB;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                -
                                <<B->D^*lnu::A_FB;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                                )"),

                make_cacheable_observable("B->D^*lnu::F_L_CP_specific", R"(F_{\mathrm{L}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_f_L,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_expression_observable("B->D^*lnu::F_L", R"(F_{\mathrm{L}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::F_L_CP_specific;cp-conjugate=false>>
                                +
                                0.5 * <<B->D^*lnu::F_L_CP_specific;cp-conjugate=true>>
                                )"),

                make_expression_observable("B->D^*lnu::Fbar_L", R"(\bar{F}_{\mathrm{L}}(B\to \bar{D}^*\ell^-\bar\nu)_{\ell=e,\mu})",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::F_L;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                +
                                0.5 * <<B->D^*lnu::F_L;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                                )"),

                make_expression_observable("B->D^*lnu::DeltaF_L", R"()",
                                Unit::None(),
                                R"(
                                <<B->D^*lnu::F_L;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                -
                                <<B->D^*lnu::F_L;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                                )"),

                make_cacheable_observable("B->D^*lnu::Ftilde_L_CP_specific", R"(\tilde{F}_{\mathrm{L}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_ftilde_L,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, { "I", "1/2" } }),

                make_expression_observable("B->D^*lnu::Ftilde_L", R"(\tilde{F}_{\mathrm{L}}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::Ftilde_L_CP_specific;cp-conjugate=false>>
                                +
                                0.5 * <<B->D^*lnu::Ftilde_L_CP_specific;cp-conjugate=true>>
                                )"),

                make_expression_observable("B->D^*lnu::Ftildebar_L", R"(\bar{\tilde{F}}_{\mathrm{L}}(B\to \bar{D}^*\ell^-\bar\nu)_{\ell=e,\mu})",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::Ftilde_L;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                +
                                0.5 * <<B->D^*lnu::Ftilde_L;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                                )"),

                make_expression_observable("B->D^*lnu::DeltaFtilde_L", R"()",
                                Unit::None(),
                                R"(
                                <<B->D^*lnu::Ftilde_L;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                -
                                <<B->D^*lnu::Ftilde_L;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                                )"),

                /*  CP-symmetric normalized observables
                *
                *    S_i ~ (J_i + barJ_i) / (Gam_i + barGam_i)
                */
                make_expression_observable("B->D^*lnu::S_1c", R"(S_{1c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_1c;cp-conjugate=false>> + <<B->D^*lnu::J_1c;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::S_1s", R"(S_{1s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_1s;cp-conjugate=false>> + <<B->D^*lnu::J_1s;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::S_2c", R"(S_{2c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_2c;cp-conjugate=false>> + <<B->D^*lnu::J_2c;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::S_2s", R"(S_{2s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_2s;cp-conjugate=false>> + <<B->D^*lnu::J_2s;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::S_3", R"(S_3(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_3;cp-conjugate=false>> + <<B->D^*lnu::J_3;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::Sbar_3", R"()",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::S_3;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                +
                                0.5 * <<B->D^*lnu::S_3;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                                )"),

                make_expression_observable("B->D^*lnu::DeltaS_3", R"()",
                                Unit::None(),
                                R"(
                                <<B->D^*lnu::S_3;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                -
                                <<B->D^*lnu::S_3;l=e>>[q2_max=>q2_e_max,q2_min=>q2_e_min]
                                )"),

                make_expression_observable("B->D^*lnu::S_4", R"(S_4(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_4;cp-conjugate=false>> + <<B->D^*lnu::J_4;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::S_5", R"(S_5(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_5;cp-conjugate=false>> + <<B->D^*lnu::J_5;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::S_6c", R"(S_{6c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_6c;cp-conjugate=false>> + <<B->D^*lnu::J_6c;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::S_6s", R"(S_{6s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_6s;cp-conjugate=false>> + <<B->D^*lnu::J_6s;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::S_7", R"(S_7(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_7;cp-conjugate=false>> + <<B->D^*lnu::J_7;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::S_8", R"(S_8(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_8;cp-conjugate=false>> + <<B->D^*lnu::J_8;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::S_9", R"(S_9(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_9;cp-conjugate=false>> + <<B->D^*lnu::J_9;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                /*  CP-asymmetric normalized observables
                *
                *    A_i ~ (J_i - barJ_i) / (Gam_i + barGam_i)
                */
                make_expression_observable("B->D^*lnu::A_1c", R"(A_{1c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_1c;cp-conjugate=false>> - <<B->D^*lnu::J_1c;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_1s", R"(A_{1s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_1s;cp-conjugate=false>> - <<B->D^*lnu::J_1s;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_2c", R"(A_{2c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_2c;cp-conjugate=false>> - <<B->D^*lnu::J_2c;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_2s", R"(A_{2s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_2s;cp-conjugate=false>> - <<B->D^*lnu::J_2s;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_3", R"(A_3(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_3;cp-conjugate=false>> - <<B->D^*lnu::J_3;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_4", R"(A_4(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_4;cp-conjugate=false>> - <<B->D^*lnu::J_4;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_5", R"(A_5(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_5;cp-conjugate=false>> - <<B->D^*lnu::J_5;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_6c", R"(A_{6c}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_6c;cp-conjugate=false>> - <<B->D^*lnu::J_6c;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_6s", R"(A_{6s}(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_6s;cp-conjugate=false>> - <<B->D^*lnu::J_6s;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_7", R"(A_7(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_7;cp-conjugate=false>> - <<B->D^*lnu::J_7;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_8", R"(A_8(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_8;cp-conjugate=false>> - <<B->D^*lnu::J_8;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_expression_observable("B->D^*lnu::A_9", R"(A_9(B\to \bar{D}^*\ell^-\bar\nu))",
                                Unit::None(),
                                R"(
                                0.5 * (<<B->D^*lnu::J_9;cp-conjugate=false>> - <<B->D^*lnu::J_9;cp-conjugate=true>>)
                                /
                                <<B->D^*lnu::normGamma>>
                                )"),

                make_observable("B->D^*lnu::P(w_min,w_max)",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::integrated_pdf_w,
                                std::make_tuple("w_min", "w_max"),
                                Options{ { "U", "c" }, { "I", "1/2" } }),

                make_expression_observable("B->D^*lnu::Pbar(w_min,w_max)", R"()",
                                Unit::None(),
                                R"(
                                0.5 * <<B->D^*lnu::P(w_min,w_max);l=mu>>[w_mu_max=>w_max,w_mu_min=>w_min]
                                +
                                0.5 * <<B->D^*lnu::P(w_min,w_max);l=e>>[w_e_max=>w_max,w_e_min=>w_min]
                                )"),

                make_expression_observable("B->D^*lnu::DeltaP(w_min,w_max)", R"()",
                                Unit::None(),
                                R"(
                                <<B->D^*lnu::P(w_min,w_max);l=mu>>[w_mu_max=>w_max,w_mu_min=>w_min]
                                -
                                <<B->D^*lnu::P(w_min,w_max);l=e>>[w_e_max=>w_max,w_e_min=>w_min]
                                )"),

                // B -> D pi l nu
                make_observable("B->Dpilnu::P(cos(theta_D))",
                                Unit::None(),
                                &BToDPiLeptonNeutrino::differential_pdf_d,
                                std::make_tuple("cos(theta_D)")),

                make_observable("B->Dpilnu::P(cos(theta_l))",
                                Unit::None(),
                                &BToDPiLeptonNeutrino::differential_pdf_l,
                                std::make_tuple("cos(theta_l)")),

                make_observable("B->Dpilnu::P(phi)",
                                Unit::None(),
                                &BToDPiLeptonNeutrino::differential_pdf_chi,
                                std::make_tuple("phi")),

                make_observable("B->Dpilnu::P(w)",
                                Unit::None(),
                                &BToDPiLeptonNeutrino::differential_pdf_w,
                                std::make_tuple("w")),

                make_observable("B->Dpilnu::P(q2)",
                                Unit::InverseGeV2(),
                                &BToDPiLeptonNeutrino::differential_pdf_q2,
                                std::make_tuple("q2")),

                make_observable("B->Dpilnu::A_l",
                                Unit::None(),
                                &BToDPiLeptonNeutrino::integrated_lepton_polarization,
                                std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Dpilnu::P(cos(theta_D)_min,cos(theta_D)_max)",
                                Unit::None(),
                                &BToDPiLeptonNeutrino::integrated_pdf_d,
                                std::make_tuple("cos(theta_D)_min", "cos(theta_D)_max")),

                make_observable("B->Dpilnu::P(cos(theta_l)_min,cos(theta_l)_max)",
                                Unit::None(),
                                &BToDPiLeptonNeutrino::integrated_pdf_l,
                                std::make_tuple("cos(theta_l)_min", "cos(theta_l)_max")),

                make_observable("B->Dpilnu::P(phi_min,phi_max)",
                                Unit::None(),
                                &BToDPiLeptonNeutrino::integrated_pdf_chi,
                                std::make_tuple("phi_min", "phi_max")),

                make_observable("B->Dpilnu::P(w_min,w_max)",
                                Unit::None(),
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
            R"(The option "l" selects the charged lepton flavor.)"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // B_s -> D_s^* l nu
                make_observable("B_s->D_s^*lnu::dBR/dq2", R"(d\mathcal{B}(B_s\to \bar{D}_s^*\ell^-\bar\nu)/dq^2)",
                                Unit::InverseGeV2(),
                                &BToVectorLeptonNeutrino::differential_branching_ratio,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::normdBR/dq2",
                                Unit::InverseGeV2(),
                                &BToVectorLeptonNeutrino::normalized_differential_branching_ratio,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::A_FB(q2)", R"(A_{\mathrm{FB}}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_a_fb_leptonic,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_1c(q2)", R"(J_{1c}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J1c,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_1s(q2)", R"(J_{1s}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J1s,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_2c(q2)", R"(J_{2c}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J2c,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_2s(q2)", R"(J_{2s}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J2s,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_3(q2)", R"(J_{3}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J3,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_4(q2)", R"(J_{4}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J4,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_5(q2)", R"(J_{5}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J5,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_6c(q2)", R"(J_{6c}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J6c,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_6s(q2)", R"(J_{6s}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J6s,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_7(q2)", R"(J_{7}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J7,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_8(q2)", R"(J_{8}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J8,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::J_9(q2)", R"(J_{9}(B_s\to \bar{D}_s^*\ell^-\bar\nu)(q^2))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::differential_J9,
                                std::make_tuple("q2"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::BR", R"(\mathcal{B}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::integrated_branching_ratio,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::normBR",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::normalized_integrated_branching_ratio,
                                std::make_tuple("q2_min", "q2_max"),
                                Options{ { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_expression_observable("B_s->D_s^*lnu::R_D_s^*(q2)", R"(R_{D_s^*}(q^2))",
                                Unit::None(),
                                R"(
                                <<B_s->D_s^*lnu::dBR/dq2;U=c,q=s,l=tau>>
                                /
                                <<B_s->D_s^*lnu::dBR/dq2;U=c,q=s,l=mu>>
                                )"),

                make_expression_observable("B_s->D_s^*lnu::R_D_s^*", R"(R_{D_s^*})",
                                Unit::None(),
                                R"(
                                <<B_s->D_s^*lnu::BR;U=c,q=s,l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                                /
                                <<B_s->D_s^*lnu::BR;U=c,q=s,l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                                )"),

                make_cacheable_observable("B_s->D_s^*lnu::A_FB", R"(A_{\mathrm{FB}}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_fb_leptonic,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::A_L", R"(A_L(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_amplitude_polarization_L,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::A_T", R"(A_T(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_amplitude_polarization_T,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::F_L", R"(F_{\mathrm{L}}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_f_L,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::A_C^1", R"(A_{\mathrm{C}}^1(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_c_1,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::A_C^2", R"(A_{\mathrm{C}}^2(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_c_2,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::A_C^3", R"(A_{\mathrm{C}}^3(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_c_3,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::A_T^1", R"(A_{\mathrm{T}}^1(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_t_1,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::A_T^2", R"(A_{\mathrm{T}}^2(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_t_2,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::A_T^3", R"(A_{\mathrm{T}}^3(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_a_t_3,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_1c", R"(J_{1c}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J1c,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_1s", R"(J_{1s}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J1s,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_2c", R"(J_{2c}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J2c,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_2s", R"(J_{2s}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J2s,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_3", R"(J_{3}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J3,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_4", R"(J_{4}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J4,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_5", R"(J_{5}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J5,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_6c", R"(J_{6c}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J6c,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_6s", R"(J_{6s}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J6s,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_7", R"(J_{7}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J7,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_8", R"(J_{8}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J8,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_cacheable_observable("B_s->D_s^*lnu::J_9", R"(J_{9}(B_s\to \bar{D}_s^*\ell^-\bar\nu))",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::prepare,
                                &BToVectorLeptonNeutrino::integrated_J9,
                                std::make_tuple("q2_min", "q2_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),

                make_observable("B_s->D_s^*lnu::P(w_min,w_max)",
                                Unit::None(),
                                &BToVectorLeptonNeutrino::integrated_pdf_w,
                                std::make_tuple("w_min", "w_max"),
                                { { "U", "c" }, {"q", "s"}, {"I", "0"} }),
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
            R"(The option "l" selects the charged lepton flavor. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // B_s -> K^* l nubar
                make_observable("B_s->K^*lnu::dBR/dq2", R"(d\mathcal{B}(B_s\to \bar{K}^*\ell^-\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToVectorLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::normdBR/dq2", R"(d\mathcal{B}(B_s\to \bar{K}^*\ell^-\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &BToVectorLeptonNeutrino::normalized_differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::A_FB(q2)", R"(A_{\mathrm{FB}}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_a_fb_leptonic,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_1c(q2)", R"(J_{1c}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J1c,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_1s(q2)", R"(J_{1s}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J1s,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_2c(q2)", R"(J_{2c}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J2c,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_2s(q2)", R"(J_{2s}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J2s,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_3(q2)", R"(J_{3}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J3,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_4(q2)", R"(J_{4}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J4,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_5(q2)", R"(J_{5}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J5,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_6c(q2)", R"(J_{6c}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J6c,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_6s(q2)", R"(J_{6s}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J6s,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_7(q2)", R"(J_{7}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J7,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_8(q2)", R"(J_{8}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J8,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::J_9(q2)", R"(J_{9}(B_s\to \bar{K}^*\ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::differential_J9,
                        std::make_tuple("q2"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::BR", R"(\mathcal{B}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::normBR",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::A_FB", R"(A_{\mathrm{FB}}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::A_L", R"(A_L(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_amplitude_polarization_L,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::A_T", R"(A_T(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_amplitude_polarization_T,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::F_L", R"(F_{\mathrm{L}}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_f_L,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::A_C^1", R"(A_{\mathrm{C}}^1(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_a_c_1,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::A_C^2", R"(A_{\mathrm{C}}^2(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_a_c_2,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::A_C^3", R"(A_{\mathrm{C}}^3(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_a_c_3,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::A_T^1", R"(A_{\mathrm{T}}^1(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_a_t_1,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::A_T^2", R"(A_{\mathrm{T}}^2(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_a_t_2,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::A_T^3", R"(A_{\mathrm{T}}^3(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_a_t_3,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_1c", R"(J_{1c}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J1c,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_1s", R"(J_{1s}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J1s,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_2c", R"(J_{2c}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J2c,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_2s", R"(J_{2s}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J2s,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_3", R"(J_{3}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J3,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_4", R"(J_{4}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J4,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_5", R"(J_{5}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J5,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_6c", R"(J_{6c}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J6c,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_6s", R"(J_{6s}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J6s,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_7", R"(J_{7}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J7,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_8", R"(J_{8}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J8,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_cacheable_observable("B_s->K^*lnu::J_9", R"(J_{9}(B_s\to \bar{K}^*\ell^-\bar\nu))",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::prepare,
                        &BToVectorLeptonNeutrino::integrated_J9,
                        std::make_tuple("q2_min", "q2_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),

                make_observable("B_s->K^*lnu::P(w_min,w_max)",
                        Unit::None(),
                        &BToVectorLeptonNeutrino::integrated_pdf_w,
                        std::make_tuple("w_min", "w_max"),
                        { { "U", "u" }, {"q", "s"}, {"I", "1/2"} }),
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
            R"(The option "l" selects the charged lepton flavor. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("B->pipilnu::BR(q2,k2)", R"(d^2\mathcal{B}(B\to \pi\pi \ell^-\bar\nu)/(dq^2\,dk^2))",
                        Unit::InverseGeV4(),
                        &BToPiPiLeptonNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "k2")),

                make_observable("B->pipilnu::BR(q2,k2,cos(theta_pi))",
                        Unit::InverseGeV4(),
                        &BToPiPiLeptonNeutrino::triple_differential_branching_ratio,
                        std::make_tuple("q2", "k2", "cos(theta_pi)")),

                make_observable("B->pipilnu::A_FB(q2,k2)", R"(A_{\mathrm{FB}}(B\to \pi\pi \ell^-\bar\nu)(q^2,k^2))",
                        Unit::InverseGeV4(),
                        &BToPiPiLeptonNeutrino::double_differential_forward_backward_asymmetry,
                        std::make_tuple("q2", "k2")),

                make_observable("B->pipilnu::P(cos(theta_pi))",
                        Unit::None(),
                        &BToPiPiLeptonNeutrino::partial_waves,
                        std::make_tuple("q2", "k2", "cos(theta_pi)")),

                make_observable("B->pipilnu::BR", R"(\mathcal{B}(B\to \pi\pi \ell^-\bar\nu))",
                        Unit::None(),
                        &BToPiPiLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max", "k2_min", "k2_max", "z_min", "z_max")),

                make_observable("B->pipilnu::A_FB", R"(A_{\mathrm{FB}}(B\to \pi\pi \ell^-\bar\nu))",
                        Unit::None(),
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
            R"(The option "l" selects the charged lepton flavor. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                // Lambda_b -> Lambda_c l nu
                make_observable("Lambda_b->Lambda_clnu::dBR/dq2", R"(d\mathcal{B}(\Lambda_b\to\Lambda_c \ell^-\bar\nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &LambdaBToLambdaCLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^l(q2)", R"(A_{\mathrm{FB}}^\ell(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_leptonic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^h(q2)", R"(A_{\mathrm{FB}}^h(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_hadronic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^c(q2)", R"(A_{\mathrm{FB}}^{h\ell}(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::differential_a_fb_combined,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::F_0(q2)", R"(F_0(\Lambda_b\to\Lambda_c \ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::differential_fzero,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_clnu::BR", R"(\mathcal{B}(\Lambda_b\to\Lambda_c \ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("Lambda_b->Lambda_clnu::R(Lambda_c)", R"(R(\Lambda_c))",
                        Unit::None(),
                        R"(
                        <<Lambda_b->Lambda_clnu::BR;l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                        /
                        <<Lambda_b->Lambda_clnu::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        )"),

                make_observable("Lambda_b->Lambda_clnu::A_FB^l",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::A_FB^h",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_hadronic,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("Lambda_b->Lambda_clnu::R(A_FB^h)(q2)", R"(R(A_{\mathrm{FB}}^{\Lambda_c})(q^2))",
                        Unit::None(),
                        R"(
                        <<Lambda_b->Lambda_clnu::A_FB^h(q2);l=tau>>
                        /
                        <<Lambda_b->Lambda_clnu::A_FB^h(q2);l=mu>>
                        )"),

                make_expression_observable("Lambda_b->Lambda_clnu::R(A_FB^h)", R"(R(A_{\mathrm{FB}}^{\Lambda_c}))",
                        Unit::None(),
                        R"(
                        <<Lambda_b->Lambda_clnu::A_FB^h;l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                        /
                        <<Lambda_b->Lambda_clnu::A_FB^h;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        )"),

                make_observable("Lambda_b->Lambda_clnu::A_FB^c",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_combined,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::F_0",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_fzero,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_1ss", R"(K_{1ss}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k1ss,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_1cc", R"(K_{1cc}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k1cc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_1c", R"(K_{1c}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k1c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_2ss", R"(K_{2ss}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k2ss,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_2cc", R"(K_{2cc}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k2cc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_2c", R"(K_{2c}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k2c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_3sc", R"(K_{3sc}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k3sc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_3s", R"(K_{3s}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k3s,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_4sc", R"(K_{4sc}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k4sc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_clnu::K_4s", R"(K_{4s}(\Lambda_b\to\Lambda_c(\to \Lambda\pi)\ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaCLeptonNeutrino::integrated_k4s,
                        std::make_tuple("q2_min", "q2_max")),

                // Lambda_b -> Lambda_c(2595) l nubar
                make_observable("Lambda_b->Lambda_c(2595)lnu::dBR/ds", R"(d\mathcal{B}/dq^2(\Lambda_b\to\Lambda_c(2595) \ell^-\bar\nu))",
                        Unit::InverseGeV2(),
                        &LambdaBToLambdaC2595LeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::dBR/dsdcos(theta_l)",
                        Unit::InverseGeV2(),
                        &LambdaBToLambdaC2595LeptonNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "cos(theta_l)")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::BR", R"(\mathcal{B}(\Lambda_b\to\Lambda_c(2595) \ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaC2595LeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::A_FB", R"(A_\mathrm{FB}(\Lambda_b\to\Lambda_c(2595) \ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaC2595LeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2595)lnu::Gamma_normalized(q2_min,q2_max)",
                        Unit::None(),
                        &LambdaBToLambdaC2595LeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("Lambda_b->Lambda_c(2595)lnu::R_Lambda_c(2595)(q2)", R"(R_{\Lambda_c(2595)}(q^2))",
                        Unit::None(),
                        R"(
                        <<Lambda_b->Lambda_c(2595)lnu::dBR/ds;l=tau>>
                        /
                        <<Lambda_b->Lambda_c(2595)lnu::dBR/ds;l=mu>>
                        )"),

                make_expression_observable("Lambda_b->Lambda_c(2595)lnu::R_Lambda_c(2595)", R"(R_{\Lambda_c(2595)})",
                        Unit::None(),
                        R"(
                        <<Lambda_b->Lambda_c(2595)lnu::BR;l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                        /
                        <<Lambda_b->Lambda_c(2595)lnu::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        )"),

                // Lambda_b -> Lambda_c(2625) l nubar
                make_observable("Lambda_b->Lambda_c(2625)lnu::dBR/ds", R"(d\mathcal{B}/dq^2(\Lambda_b\to\Lambda_c(2625) \ell^-\bar\nu))",
                        Unit::InverseGeV2(),
                        &LambdaBToLambdaC2625LeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::A_FB(q2)", R"(A_\mathrm{FB}(\Lambda_b\to\Lambda_c(2625) \ell^-\bar\nu)(q^2))",
                        Unit::None(),
                        &LambdaBToLambdaC2625LeptonNeutrino::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::dBR/dsdcos(theta_l)",
                        Unit::InverseGeV2(),
                        &LambdaBToLambdaC2625LeptonNeutrino::double_differential_branching_ratio,
                        std::make_tuple("q2", "cos(theta_l)")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::BR", R"(\mathcal{B}(\Lambda_b\to\Lambda_c(2625) \ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaC2625LeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::A_FB", R"(A_\mathrm{FB}(\Lambda_b\to\Lambda_c(2625) \ell^-\bar\nu))",
                        Unit::None(),
                        &LambdaBToLambdaC2625LeptonNeutrino::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambda_c(2625)lnu::Gamma_normalized(q2_min,q2_max)",
                        Unit::None(),
                        &LambdaBToLambdaC2625LeptonNeutrino::normalized_integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_expression_observable("Lambda_b->Lambda_c(2625)lnu::R_Lambda_c(2625)(q2)", R"(R_{\Lambda_c(2625)}(q^2))",
                        Unit::None(),
                        R"(
                        <<Lambda_b->Lambda_c(2625)lnu::dBR/ds;l=tau>>
                        /
                        <<Lambda_b->Lambda_c(2625)lnu::dBR/ds;l=mu>>
                        )"),

                make_expression_observable("Lambda_b->Lambda_c(2625)lnu::R_Lambda_c(2625)", R"(R_{\Lambda_c(2625)})",
                        Unit::None(),
                        R"(
                        <<Lambda_b->Lambda_c(2625)lnu::BR;l=tau>>[q2_max=>q2_tau_max,q2_min=>q2_tau_min]
                        /
                        <<Lambda_b->Lambda_c(2625)lnu::BR;l=mu>>[q2_max=>q2_mu_max,q2_min=>q2_mu_min]
                        )")
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
            R"(Miscellaneous)",
            R"()",
            {
                /* B Meson Properties */
                make_observable("B::M_B^*-M_B", R"(M_{B^*} - M_B)",
                        Unit::GeV(),
                        &BMesonProperties::mass_splitting_j1_j0),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Class-I nonleptonic heavy-to-heavy
    // {{{
    ObservableGroup
    make_classI_nonleptonic_heavy_to_heavy_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Class-I Nonleptonic Heavy-to-Heavy Decays)",
            R"()",
            {
                /* B_s -> D_s pi */
                make_observable("B_s^0->D_s^+pi^-::BR", R"(\mathcal{B}(\bar{B}_s^0\to D_s^+\pi^-))",
                        Unit::None(),
                        &BqToDqPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s"} }),

                /* B -> D K */
                make_observable("B^0->D^+K^-::BR", R"(\mathcal{B}(\bar{B}^0\to D^+K^-))",
                        Unit::None(),
                        &BqToDqPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d"} }),

                /* B_s -> D_s^* pi */
                make_observable("B_s^0->D_s^*+pi^-::BR", R"(\mathcal{B}(\bar{B}_s^0\to D_s^{*+}\pi^-))",
                        Unit::None(),
                        &BqToDstarqPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s"} }),

                /* B -> D^* K */
                make_observable("B^0->D^*+K^-::BR", R"(\mathcal{B}(\bar{B}^0\to D^{*+}K^-))",
                        Unit::None(),
                        &BqToDstarqPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d"} }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Observables related to the B -> P(seudoscalar) P(seudoscalar) decays
    // {{{
    ObservableGroup
    make_nonleptonic_b_to_p_p_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables related to the $B\to PP$ decays)",
            R"()",
            {
                make_observable("B^+->pi^+pi^0::BR", R"(\mathcal{B}(B^+\to\pi^+\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "pi^+" }, { "P2", "pi^0" } }),

                make_observable("B^+->pi^+pi^0::BR_avg", R"(\bar{\mathcal{B}}(B^+\to\pi^+\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "pi^+" }, { "P2", "pi^0" } }),

                make_observable("B^+->pi^+pi^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^+\to\pi^+\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "pi^+" }, { "P2", "pi^0" } }),

                make_observable("B^+->pi^+pi^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^+\to\pi^+\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "pi^+" }, { "P2", "pi^0" } }),

                make_observable("B^+->K^+Kbar^0::BR", R"(\mathcal{B}(B^+\to K^+\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_u" }, { "P2", "Kbar_d" } }),

                make_observable("B^+->K^+Kbar^0::BR_avg", R"(\bar{\mathcal{B}}(B^+\to K^+\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_u" }, { "P2", "Kbar_d" } }),

                make_observable("B^+->K^+Kbar^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^+\to K^+\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_u" }, { "P2", "Kbar_d" } }),

                make_observable("B^+->K^+Kbar^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^+\to K^+\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_u" }, { "P2", "Kbar_d" } }),

                make_observable("B^0->pi^+pi^-::BR", R"(\mathcal{B}(B^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^+" }, { "P2", "pi^-" } }),

                make_observable("B^0->pi^+pi^-::BR_avg", R"(\bar{\mathcal{B}}(B^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^+" }, { "P2", "pi^-" } }),

                make_observable("B^0->pi^+pi^-::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^+" }, { "P2", "pi^-" } }),

                make_observable("B^0->pi^+pi^-::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^+" }, { "P2", "pi^-" } }),

                make_observable("B^0->pi^+pi^-::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^+" }, { "P2", "pi^-" }}),

                make_observable("B^0->pi^+pi^-::A_DG", R"(A^{\Delta\Gamma}_\mathrm{CP}(B^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::a_Delta_Gamma,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^+" }, { "P2", "pi^-" }}),

                make_observable("B^0->pi^0pi^0::BR", R"(\mathcal{B}(B^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B^0->pi^0pi^0::BR_avg", R"(\bar{\mathcal{B}}(B^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B^0->pi^0pi^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B^0->pi^0pi^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B^0->pi^0pi^0::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B^0->pi^0pi^0::A_DG", R"(A^{\Delta\Gamma}_\mathrm{CP}(B^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::a_Delta_Gamma,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "pi^0" }, { "P2", "pi^0" }}),

                make_observable("B^0->K^+K^-::BR", R"(\mathcal{B}(B^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B^0->K^+K^-::BR_avg", R"(\bar{\mathcal{B}}(B^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B^0->K^+K^-::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B^0->K^+K^-::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B^0->K^+K^-::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B^0->K^+K^-::A_DG", R"(A^{\Delta\Gamma}_\mathrm{CP}(B^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::a_Delta_Gamma,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B^0->K^0Kbar^0::BR", R"(\mathcal{B}(B^0\to K^0\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_d" }, { "P2", "Kbar_d" } }),

                make_observable("B^0->K^0Kbar^0::BR_avg", R"(\bar{\mathcal{B}}(B^0\to \bar{K}^0K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_d" }, { "P2", "Kbar_d" } }),

                make_observable("B^0->K^0Kbar^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to \bar{K}^0K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_d" }, { "P2", "Kbar_d" } }),

                make_observable("B^0->K^0Kbar^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^0\to \bar{K}^0K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_d" }, { "P2", "Kbar_d" } }),

                make_observable("B^0->K^0Kbar^0::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B^0\to \bar{K}^0K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "Kbar_d" }, { "P2", "K_d" } }),

                make_observable("B^0->K^0Kbar^0::A_DG", R"(A^{\Delta\Gamma}_\mathrm{CP}(B^0\to \bar{K}^0K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::a_Delta_Gamma,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "Kbar_d" }, { "P2", "K_d" } }),

                make_observable("B_s^0->K^-pi^+::BR", R"(\mathcal{B}(B_s^0\to K^-\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "Kbar_u" }, { "P2", "pi^+" } }),

                make_observable("B_s^0->K^-pi^+::BR_avg", R"(\bar{\mathcal{B}}(B_s^0\to K^-\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "Kbar_u" }, { "P2", "pi^+" } }),

                make_observable("B_s^0->K^-pi^+::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B_s^0\to K^-\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "Kbar_u" }, { "P2", "pi^+" } }),

                make_observable("B_s^0->K^-pi^+::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B_s^0\to K^-\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "Kbar_u" }, { "P2", "pi^+" } }),

                make_observable("B_s^0->K^-pi^+::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B_s^0\to K^-\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "Kbar_u" }, { "P2", "pi^+" } }),

                make_observable("B_s^0->Kbar^0pi^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B_s^0\to \bar{K}^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "Kbar_d" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->Kbar^0pi^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B_s^0\to \bar{K}^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "Kbar_d" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->K_Spi^0::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B_s^0\to K_S\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_S" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->K_Spi^0::A_DG", R"(A^{\Delta\Gamma}_\mathrm{CP}(B_s^0\to K_S\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::a_Delta_Gamma,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_S" }, { "P2", "pi^0" } }),

                make_observable("B^+->K^+pi^0::BR", R"(\mathcal{B}(B^+\to K^+\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_u" }, { "P2", "pi^0" } }),

                make_observable("B^+->K^+pi^0::BR_avg", R"(\bar{\mathcal{B}}(B^+\to K^+\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_u" }, { "P2", "pi^0" } }),

                make_observable("B^+->K^+pi^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^+\to K^+\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_u" }, { "P2", "pi^0" } }),

                make_observable("B^+->K^+pi^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^+\to K^+\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_u" }, { "P2", "pi^0" } }),

                make_observable("B^+->K^0pi^+::BR", R"(\mathcal{B}(B^+\to K^0\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_d" }, { "P2", "pi^+" } }),

                make_observable("B^+->K^0pi^+::BR_avg", R"(\bar{\mathcal{B}}(B^+\to K^0\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_d" }, { "P2", "pi^+" } }),

                make_observable("B^+->K^0pi^+::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^+\to K^0\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_d" }, { "P2", "pi^+" } }),

                make_observable("B^+->K^0pi^+::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^+\to K^0\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "K_d" }, { "P2", "pi^+" } }),

                make_observable("B^0->K^+pi^-::BR", R"(\mathcal{B}(B^0\to K^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "pi^-" } }),

                make_observable("B^0->K^+pi^-::BR_avg", R"(\bar{\mathcal{B}}(B^0\to K^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "pi^-" } }),

                make_observable("B^0->K^+pi^-::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to K^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "pi^-" } }),

                make_observable("B^0->K^+pi^-::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^0\to K^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "pi^-" } }),

                make_observable("B^0->K^+pi^-::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B^0\to K^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_u" }, { "P2", "pi^-" } }),

                make_observable("B^0->K^0pi^0::BR", R"(\mathcal{B}(B^0\to K^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_d" }, { "P2", "pi^0" } }),

                make_observable("B^0->K^0pi^0::BR_avg", R"(\bar{\mathcal{B}}(B^0\to K^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_d" }, { "P2", "pi^0" } }),

                make_observable("B^0->K^0pi^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to K^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_d" }, { "P2", "pi^0" } }),

                make_observable("B^0->K^0pi^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^0\to K^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_d" }, { "P2", "pi^0" } }),

                make_observable("B^0->K_Spi^0::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B^0\to K_S\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_S" }, { "P2", "pi^0" } }),

                make_observable("B^0->K_Spi^0::A_DG", R"(A^{\Delta\Gamma}_\mathrm{CP}(B^0\to K_S\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::a_Delta_Gamma,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "K_S" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->pi^+pi^-::BR", R"(\mathcal{B}(B_s^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^+" }, { "P2", "pi^-" } }),

                make_observable("B_s^0->pi^+pi^-::BR_avg", R"(\bar{\mathcal{B}}(B_s^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^+" }, { "P2", "pi^-" } }),

                make_observable("B_s^0->pi^+pi^-::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B_s^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^+" }, { "P2", "pi^-" } }),

                make_observable("B_s^0->pi^+pi^-::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B_s^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^+" }, { "P2", "pi^-" } }),

                make_observable("B_s^0->pi^+pi^-::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B_s^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^+" }, { "P2", "pi^-" } }),

                make_observable("B_s^0->pi^+pi^-::A_DG", R"(A^{\Delta\Gamma}_\mathrm{CP}(B_s^0\to\pi^+\pi^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::a_Delta_Gamma,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^+" }, { "P2", "pi^-" } }),

                make_observable("B_s^0->K^+K^-::BR", R"(\mathcal{B}(B_s^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B_s^0->K^+K^-::BR_avg", R"(\bar{\mathcal{B}}(B_s^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B_s^0->K^+K^-::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B_s^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B_s^0->K^+K^-::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B_s^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B_s^0->K^+K^-::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B_s^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B_s^0->K^+K^-::A_DG", R"(A^{\Delta\Gamma}_\mathrm{CP}(B_s^0\to K^+K^-))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::a_Delta_Gamma,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_u" }, { "P2", "Kbar_u" } }),

                make_observable("B_s^0->K^0Kbar^0::BR", R"(\mathcal{B}(B_s^0\to K^0\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_d" }, { "P2", "Kbar_d" } }),

                make_observable("B_s^0->K^0Kbar^0::BR_avg", R"(\bar{\mathcal{B}}(B_s^0\to K^0\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_d" }, { "P2", "Kbar_d" } }),

                make_observable("B_s^0->K^0Kbar^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B_s^0\to K^0\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_d" }, { "P2", "Kbar_d" } }),

                make_observable("B_s^0->K^0Kbar^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B_s^0\to K^0\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_d" }, { "P2", "Kbar_d" } }),

                make_observable("B_s^0->K^0Kbar^0::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B_s^0\to K^0\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_d" }, { "P2", "Kbar_d" } }),

                make_observable("B_s^0->K^0Kbar^0::A_DG", R"(A^{\Delta\Gamma}_\mathrm{CP}(B_s^0\to K^0\bar{K}^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::a_Delta_Gamma,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "K_d" }, { "P2", "Kbar_d" } }),

                make_observable("B_s^0->pi^0pi^0::BR", R"(\mathcal{B}(B_s^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->pi^0pi^0::BR_avg", R"(\bar{\mathcal{B}}(B_s^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->pi^0pi^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B_s^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->pi^0pi^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B_s^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->pi^0pi^0::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B_s^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->pi^0pi^0::A_DG", R"(A^{\Delta\Gamma}_\mathrm{CP}(B_s^0\to\pi^0\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::a_Delta_Gamma,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "pi^0" }, { "P2", "pi^0" } }),

                make_observable("B^+->etapi^+::BR", R"(\mathcal{B}(B^+\to\eta\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta" }, { "P2", "pi^+" } }),

                make_observable("B^+->etapi^+::BR_avg", R"(\bar{\mathcal{B}}(B^+\to\eta\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta" }, { "P2", "pi^+" } }),

                make_observable("B^+->etapi^+::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^+\to\eta\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta" }, { "P2", "pi^+" } }),

                make_observable("B^+->etapi^+::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^+\to\eta\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta" }, { "P2", "pi^+" } }),

                make_observable("B^0->etapi^0::BR", R"(\mathcal{B}(B^0\to\eta\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "pi^0" } }),

                make_observable("B^0->etapi^0::BR_avg", R"(\bar{\mathcal{B}}(B^0\to\eta\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "pi^0" } }),

                make_observable("B^0->etapi^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to\eta\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "pi^0" } }),

                make_observable("B^0->etapi^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^0\to\eta\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "pi^0" } }),

                make_observable("B^0->etapi^0::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B^0\to\eta\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "pi^0" } }),

                make_observable("B^+->etaK^+::BR", R"(\mathcal{B}(B^+\to\eta K^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta" }, { "P2", "K_u" } }),

                make_observable("B^+->etaK^+::BR_avg", R"(\bar{\mathcal{B}}(B^+\to\eta K^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta" }, { "P2", "K_u" } }),

                make_observable("B^+->etaK^+::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^+\to\eta K^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta" }, { "P2", "K_u" } }),

                make_observable("B^+->etaK^+::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^+\to\eta K^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta" }, { "P2", "K_u" } }),

                make_observable("B^0->etaK^0::BR", R"(\mathcal{B}(B^0\to\eta K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "K_d" } }),

                make_observable("B^0->etaK^0::BR_avg", R"(\bar{\mathcal{B}}(B^0\to\eta K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "K_d" } }),

                make_observable("B^0->etaK^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to\eta K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "K_d" } }),

                make_observable("B^0->etaK^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^0\to\eta K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "K_d" } }),

                make_observable("B^0->etaK_S::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B^0\to\eta K_S^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "K_S" } }),

                make_observable("B_s^0->etapi^0::BR", R"(\mathcal{B}(B_s^0\to\eta\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->etapi^0::BR_avg", R"(\bar{\mathcal{B}}(B_s^0\to\eta\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->etapi^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B_s^0\to\eta\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->etapi^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B_s^0\to\eta\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->etapi^0::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B_s^0\to\eta\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "pi^0" } }),

                make_observable("B_s^0->etaK^0::BR", R"(\mathcal{B}(B_s^0\to\eta\bar K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "Kbar_d" } }),

                make_observable("B_s^0->etaK^0::BR_avg", R"(\bar{\mathcal{B}}(B_s^0\to\eta\bar K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "Kbar_d" } }),

                make_observable("B_s^0->etaK^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B_s^0\to\eta\bar K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "Kbar_d" } }),

                make_observable("B_s^0->etaK^0::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B_s^0\to\eta\bar K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "Kbar_d" } }),

                make_observable("B_s^0->etaK_S::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B_s^0\to\eta\bar K_S^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "K_S" } }),

                make_observable("B^0->etaeta::BR", R"(\mathcal{B}(B^0\to\eta\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "eta" } }),

                make_observable("B^0->etaeta::BR_avg", R"(\bar{\mathcal{B}}(B^0\to\eta\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "eta" } }),

                make_observable("B^0->etaeta::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to\eta\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "eta" } }),

                make_observable("B^0->etaeta::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B^0\to\eta\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "eta" } }),

                make_observable("B^0->etaeta::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B^0\to\eta\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta" }, { "P2", "eta" } }),

                make_observable("B_s^0->etaeta::BR", R"(\mathcal{B}(B_s^0\to\eta\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "eta" } }),

                make_observable("B_s^0->etaeta::BR_avg", R"(\bar{\mathcal{B}}(B_s^0\to\eta\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "eta" } }),

                make_observable("B_s^0->etaeta::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B_s^0\to\eta\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "eta" } }),

                make_observable("B_s^0->etaeta::A_CP", R"(A^\mathrm{dir}_\mathrm{CP}(B_s^0\to\eta\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "eta" } }),

                make_observable("B_s^0->etaeta::S_CP", R"(A^\mathrm{mix}_\mathrm{CP}(B_s^0\to\eta\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta" }, { "P2", "eta" } }),

                make_observable("B^+->eta^primepi^+::BR", R"(\mathcal{B}(B^+\to\eta^\prime\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta_prime" }, { "P2", "pi^+" } }),

                make_observable("B^+->eta^primepi^+::BR_avg", R"(\bar{\mathcal{B}}(B^+\to\eta^\prime\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta_prime" }, { "P2", "pi^+" } }),

                make_observable("B^+->eta^primepi^+::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^+\to\eta^\prime\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta_prime" }, { "P2", "pi^+" } }),

                make_observable("B^+->eta^primepi^+::A_CP", R"(A\mathrm{CP}(B^+\to\eta^\prime\pi^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta_prime" }, { "P2", "pi^+" } }),

                make_observable("B^0->eta^primepi^0::BR", R"(\mathcal{B}(B^0\to\eta^\prime\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "pi^0" } }),

                make_observable("B^0->eta^primepi^0::BR_avg", R"(\bar{\mathcal{B}}(B^0\to\eta^\prime\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "pi^0" } }),

                make_observable("B^0->eta^primepi^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to\eta^\prime\pi^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "pi^0" } }),

                make_observable("B^+->eta^primeK^+::BR", R"(\mathcal{B}(B^+\to\eta^\prime K^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta_prime" }, { "P2", "K_u" } }),

                make_observable("B^+->eta^primeK^+::BR_avg", R"(\bar{\mathcal{B}}(B^+\to\eta^\prime K^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta_prime" }, { "P2", "K_u" } }),

                make_observable("B^+->eta^primeK^+::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^+\to\eta^\prime K^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta_prime" }, { "P2", "K_u" } }),

                make_observable("B^+->eta^primeK^+::A_CP", R"(A\mathrm{CP}(B^+\to\eta^\prime K^+))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "u" }, { "P1", "eta_prime" }, { "P2", "K_u" } }),

                make_observable("B^0->eta^primeK^0::BR", R"(\mathcal{B}(B^0\to\eta^\prime K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "K_d" } }),

                make_observable("B^0->eta^primeK^0::BR_avg", R"(\bar{\mathcal{B}}(B^0\to\eta^\prime K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "K_d" } }),

                make_observable("B^0->eta^primeK^0::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to\eta^\prime K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "K_d" } }),

                make_observable("B^0->eta^primeK^0::A_CP", R"(A\mathrm{CP}(B^0\to\eta^\prime K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "K_d" } }),

                make_observable("B^0->eta^primeK^0::S_CP", R"(S\mathrm{CP}(B^0\to\eta^\prime K^0))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "K_d" } }),

                make_observable("B^0->eta^primeeta^prime::BR", R"(\mathcal{B}(B^0\to\eta^\prime\eta^\prime))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "eta_prime" } }),

                make_observable("B^0->eta^primeeta^prime::BR_avg", R"(\bar{\mathcal{B}}(B^0\to\eta^\prime\eta^\prime))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "eta_prime" } }),

                make_observable("B^0->eta^primeeta^prime::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to\eta^\prime\eta^\prime))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "eta_prime" } }),

                make_observable("B^0->eta^primeeta::BR", R"(\mathcal{B}(B^0\to\eta^\prime\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "eta" } }),

                make_observable("B^0->eta^primeeta::BR_avg", R"(\bar{\mathcal{B}}(B^0\to\eta^\prime\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "eta" } }),

                make_observable("B^0->eta^primeeta::BR_exp", R"(\mathcal{B}_\mathrm{exp}(B^0\to\eta^\prime\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "d" }, { "P1", "eta_prime" }, { "P2", "eta" } }),

                make_observable("B_s^0->eta^primeeta::BR", R"(\mathcal{B}(B_s^0\to\eta^\prime\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta_prime" }, { "P2", "eta" } }),

                make_observable("B_s^0->eta^primeeta::BR_avg", R"(\bar{\mathcal{B}}(B_s^0\to\eta^\prime\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta_prime" }, { "P2", "eta" } }),

                make_observable("B_s^0->eta^primeeta::BR_exp", R"(\mathcal{B}\mathrm{exp}(B_s^0\to\eta^\prime\eta))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta_prime" }, { "P2", "eta" } }),

                make_observable("B_s^0->eta^primeeta^prime::BR", R"(\mathcal{B}(B_s^0\to\eta^\prime\eta^\prime))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta_prime" }, { "P2", "eta_prime" } }),

                make_observable("B_s^0->eta^primeeta^prime::BR_avg", R"(\bar{\mathcal{B}}(B_s^0\to\eta^\prime\eta^\prime))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::avg_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta_prime" }, { "P2", "eta_prime" } }),

                make_observable("B_s^0->eta^primeeta^prime::BR_exp", R"(\mathcal{B}\mathrm{exp}(B_s^0\to\eta^\prime\eta^\prime))",
                        Unit::None(),
                        &BToPseudoscalarPseudoscalar::exp_branching_ratio,
                        std::make_tuple(),
                        { { "q", "s" }, { "P1", "eta_prime" }, { "P2", "eta_prime" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // (Pseudo)Observables related to the B-meson lifetime
    // {{{
    ObservableGroup
    make_b_lifetime_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"((Pseudo)Observables related to the $B$-meson liftime)",
            R"()",
            {
                /* B^0 lifetime */
                make_observable("B^0::Gamma(dbcu)", R"(\Gamma(\bar{B}^0)^{dbcu})",
                        Unit::InversePicoSecond(),
                        &Lifetime::decay_width_dbcu_dim6_lo,
                        std::make_tuple(),
                        { { "q", "d"} }),

                make_observable("B^0::Gamma(sbcu)", R"(\Gamma(\bar{B}^0)^{sbcu})",
                        Unit::InversePicoSecond(),
                        &Lifetime::decay_width_sbcu_dim6_lo,
                        std::make_tuple(),
                        { { "q", "d"} }),

                /* B^- lifetime */
                make_observable("B^-::Gamma(dbcu)", R"(\Gamma(\bar{B}^-)^{dbcu})",
                        Unit::InversePicoSecond(),
                        &Lifetime::decay_width_dbcu_dim6_lo,
                        std::make_tuple(),
                        { { "q", "u"} }),

                make_observable("B^-::Gamma(sbcu)", R"(\Gamma(\bar{B}^-)^{sbcu})",
                        Unit::InversePicoSecond(),
                        &Lifetime::decay_width_sbcu_dim6_lo,
                        std::make_tuple(),
                        { { "q", "u"} }),

                /* B_s^0 lifetime */
                make_observable("B_s^0::Gamma(dbcu)", R"(\Gamma(\bar{B}_s^0)^{dbcu})",
                        Unit::InversePicoSecond(),
                        &Lifetime::decay_width_dbcu_dim6_lo,
                        std::make_tuple(),
                        { { "q", "s"} }),

                make_observable("B_s^0::Gamma(sbcu)", R"(\Gamma(\bar{B}_s^0)^{sbcu})",
                        Unit::InversePicoSecond(),
                        &Lifetime::decay_width_sbcu_dim6_lo,
                        std::make_tuple(),
                        { { "q", "s"} }),

                /* spectator independent LO B lifetime*/
                make_observable("B::Gamma(dbcu)", R"(\Gamma(\bar{B})_{\text{LO}}^{dbcu})",
                        Unit::InversePicoSecond(),
                        &Lifetime::decay_width_dbcu_dim3_lo,
                        std::make_tuple(),
                        { { "q", "s"} }),

                make_observable("B::Gamma(sbcu)", R"(\Gamma(\bar{B})_{\text{LO}}^{sbcu})",
                        Unit::InversePicoSecond(),
                        &Lifetime::decay_width_sbcu_dim3_lo,
                        std::make_tuple(),
                        { { "q", "s"} })
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

                // B^- -> l^- nubar lprime^+ lprime^-
                make_b_to_3l_nu_group(),

                // B_{u,d} -> P l^- nubar
                make_b_to_pi_l_nu_group(),
                make_b_to_d_l_nu_group(),

                // B_s -> P l^- nubar
                make_bs_to_k_l_nu_group(),
                make_bs_to_ds_l_nu_group(),

                // B_{u,d} -> V l^- nubar
                make_b_to_omega_l_nu_group(),
                make_b_to_rho_l_nu_group(),
                make_b_to_dstar_l_nu_group(),

                // B_u -> gamma l nu
                make_b_to_gamma_l_nu_group(),

                // B_s -> V l^- nubar
                make_bs_to_kstar_l_nu_group(),
                make_bs_to_dsstar_l_nu_group(),

                // B_{u,d} -> P P l^- nubar
                make_b_to_pi_pi_l_nu_group(),

                // Lambda_b
                make_lambdab_to_lambdac_l_nu_group(),

                // B -> X_u l^- nubar
                make_b_to_xu_semileptonic_group(),

                // class I nonleptonic heavy-to-heavy decays
                make_classI_nonleptonic_heavy_to_heavy_group(),

                // nonleptonic B -> PP decays
                make_nonleptonic_b_to_p_p_group(),

                // B-meson lifetime
                make_b_lifetime_group(),
            }
        );

        return ObservableSection(imp);
    }
}
