/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2023-2026 Danny van Dyk
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

#include <eos/signal-pdf-impl.hh>
#include <eos/b-decays/b-to-d-l-x-nu.hh>
#include <eos/b-decays/b-to-gamma-l-nu.hh>
#include <eos/b-decays/b-to-psd-l-nu.hh>
#include <eos/b-decays/b-to-pi-l-x-nu.hh>
#include <eos/b-decays/b-to-pi-pi-l-nu.hh>
#include <eos/b-decays/b-to-vec-l-nu.hh>
#include <eos/b-decays/b-to-3l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac2625-l-nu.hh>
#include <eos/utils/concrete-signal-pdf.hh>

using std::literals::string_literals::operator""s;

namespace eos
{
    // Leptonic and photoleptonic B decays
    // {{{
    SignalPDFGroup
    make_b_to_leptons_pdf_group()
    {
        auto imp = new Implementation<SignalPDFGroup>(
            R"(Signal PDFs in leptonic and photoleptonic $B$ decays)",
            R"()",
            {
                // B -> gamma l nu
                make_signal_pdf("B_u->gammalnu::P(E_gamma,cos(theta_l))",
                    R"(PDF for the decay $B^-\to \gamma \ell^- \bar\nu$ as a function of the photon energy $E_\gamma$
                    and the angle $\theta_\ell$ between the charged lepton and the photon in the $B$ rest frame.)",
                    Options{ },
                    "B_u->gammalnu::UnnormalizedPDF(E_gamma,cos(theta_l))",
                    std::make_tuple(
                        "E_gamma"s,
                        "cos(theta_l)"s
                    ),
                    "B_u->gammalnu::NormalizationPDF(E_gamma,cos(theta_l))",
                    std::make_tuple(
                        "E_gamma_min"s
                    )
                ),

                // B -> 3l nu
                make_signal_pdf("B_u->enumumu::P(q2,k2,z_gamma,z_w,phi)",
                    R"(PDF for the decay $B^-\to e^-\bar\nu \mu^+\mu^-$ as a function of the invariant $\mu^+\mu^-$ mass squared $q^2$,
                    the invariant $e^-\bar\nu$ mass squared $k^2$, the cosine of the photon helicity angle $z_\gamma$, the cosine
                    of the $W$ helicity angle $z_W$, and the azimuthal angle $\phi$.)",
                    Options{ { "l"_ok, "e"_ov }, { "lprime"_ok, "mu"_ov } },
                    "B_u->enumumu::UnnormalizedPDF(q2,k2,z_gamma,z_w,phi)",
                    std::make_tuple(
                        "q2"s,
                        "k2"s,
                        "z_gamma"s,
                        "z_w"s,
                        "phi"s
                    ),
                    "B_u->enumumu::NormalizationPDF(q2,k2,z_gamma,z_w,phi)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s,
                        "k2_min"s,
                        "k2_max"s
                    )
                ),

                make_signal_pdf("B_u->munuee::P(q2,k2,z_gamma,z_w,phi)",
                    R"(PDF for the decay $B^-\to \mu^-\bar\nu e^+e^-$ as a function of the invariant $e^+e^-$ mass squared $q^2$,
                    the invariant $\mu^-\bar\nu$ mass squared $k^2$, the cosine of the photon helicity angle $z_\gamma$, the cosine
                    of the $W$ helicity angle $z_W$, and the azimuthal angle $\phi$.)",
                    Options{ { "l"_ok, "mu"_ov }, { "lprime"_ok, "e"_ov } },
                    "B_u->munuee::UnnormalizedPDF(q2,k2,z_gamma,z_w,phi)",
                    std::make_tuple(
                        "q2"s,
                        "k2"s,
                        "z_gamma"s,
                        "z_w"s,
                        "phi"s
                    ),
                    "B_u->munuee::NormalizationPDF(q2,k2,z_gamma,z_w,phi)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s,
                        "k2_min"s,
                        "k2_max"s
                    )
                )
            }
        );

        return SignalPDFGroup(imp);
    }
    // }}}

    // Semileptonic B -> P(seudoscalar) decays
    // {{{

    SignalPDFGroup
    make_b_to_p_l_nu_pdf_group()
    {
        auto imp = new Implementation<SignalPDFGroup>(
            R"(Signal PDFs in semileptonic $B\to P \ell^-\bar\nu$ decays)",
            R"()",
            {
                // B -> pi l nu
                make_signal_pdf("B->pilnu::P(q2)",
                    R"(PDF for the decay $\bar{B}\to \pi \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "P"_ok, "pi"_ov } },
                    "B->pilnu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B->pilnu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                make_signal_pdf("B->pilnu::P(q2,cos(theta_l))",
                    R"(PDF for the decay $\bar{B}\to \pi \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$
                    and the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $\bar{B}$ flight direction in the $\ell^-\bar\nu$ rest frame.)",
                    Options{ { "P"_ok, "pi"_ov } },
                    "B->pilnu::UnnormalizedPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s
                    ),
                    "B->pilnu::NormalizationPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                // B -> D l nu
                make_signal_pdf("B->Dlnu::P(q2)",
                    R"(PDF for the decay $\bar{B}\to D \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "P"_ok, "D"_ov } },
                    "B->Dlnu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B->Dlnu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                make_signal_pdf("B->Dlnu::P(q2,cos(theta_l))",
                    R"(PDF for the decay $\bar{B}\to D \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$
                    and the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $\bar{B}$ flight direction in the $\ell^-\bar\nu$ rest frame.)",
                    Options{ { "P"_ok, "D"_ov } },
                    "B->Dlnu::UnnormalizedPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s
                    ),
                    "B->Dlnu::NormalizationPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                // B_s -> K l nu
                make_signal_pdf("B_s->Klnu::P(q2)",
                    R"(PDF for the decay $\bar{B}_s\to K \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "P"_ok, "K"_ov } },
                    "B_s->Klnu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B_s->Klnu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                make_signal_pdf("B_s->Klnu::P(q2,cos(theta_l))",
                    R"(PDF for the decay $\bar{B}_s\to K \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$
                    and the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $\bar{B}_s$ flight direction in the $\ell^-\bar\nu$ rest frame.)",
                    Options{ { "P"_ok, "K"_ov } },
                    "B_s->Klnu::UnnormalizedPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s
                    ),
                    "B_s->Klnu::NormalizationPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                // B_s -> D_s l nu
                make_signal_pdf("B_s->D_slnu::P(q2)",
                    R"(PDF for the decay $\bar{B}_s\to D_s \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "P"_ok, "D_s"_ov } },
                    "B_s->D_slnu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B_s->D_slnu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                make_signal_pdf("B_s->D_slnu::P(q2,cos(theta_l))",
                    R"(PDF for the decay $\bar{B}_s\to D_s \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$
                    and the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $\bar{B}_s$ flight direction in the $\ell^-\bar\nu$ rest frame.)",
                    Options{ { "P"_ok, "D_s"_ov } },
                    "B_s->D_slnu::UnnormalizedPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s
                    ),
                    "B_s->D_slnu::NormalizationPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                // B -> eta l nu
                make_signal_pdf("B->etalnu::P(q2)",
                    R"(PDF for the decay $\bar{B}\to \eta \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "P"_ok, "eta"_ov } },
                    "B->etalnu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B->etalnu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                make_signal_pdf("B->etalnu::P(q2,cos(theta_l))",
                    R"(PDF for the decay $\bar{B}\to \eta \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$
                    and the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $\bar{B}$ flight direction in the $\ell^-\bar\nu$ rest frame.)",
                    Options{ { "P"_ok, "eta"_ov } },
                    "B->etalnu::UnnormalizedPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s
                    ),
                    "B->etalnu::NormalizationPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                // B -> eta' l nu
                make_signal_pdf("B->eta_primelnu::P(q2)",
                    R"(PDF for the decay $\bar{B}\to \eta' \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "P"_ok, "eta_prime"_ov } },
                    "B->eta_primelnu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B->eta_primelnu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                make_signal_pdf("B->eta_primelnu::P(q2,cos(theta_l))",
                    R"(PDF for the decay $\bar{B}\to \eta' \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$
                    and the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $\bar{B}$ flight direction in the $\ell^-\bar\nu$ rest frame.)",
                    Options{ { "P"_ok, "eta_prime"_ov } },
                    "B->eta_primelnu::UnnormalizedPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s
                    ),
                    "B->eta_primelnu::NormalizationPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                )
            }
        );

        return SignalPDFGroup(imp);
    }
    // }}}

    // Semileptonic B -> V(pseudoscalar) decays
    // {{{
    SignalPDFGroup
    make_b_to_v_l_nu_pdf_group()
    {
        auto imp = new Implementation<SignalPDFGroup>(
            R"(Signal PDFs in semileptonic $B\to V \ell^-\bar\nu$ decays)",
            R"()",
            {
                // B -> D^* l nu
                make_signal_pdf("B->D^*lnu::P(q2)",
                    R"(PDF for the decay $\bar{B}\to D^*(\to D\pi) \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "V"_ok, "D^*"_ov } },
                    "B->D^*lnu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B->D^*lnu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                make_signal_pdf("B->D^*lnu::P(q2,cos(theta_l),cos(theta_D),phi)",
                    R"(PDF for the decay $\bar{B}\to D^*(\to D\pi) \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$,
                    the cosine of the angle $\theta_D$ between the $D$ and the negative $\ell^-\bar\nu$ flight direction in the $D^*$ rest frame,
                    the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $D^*$ flight direction in the $\ell^-\bar\nu$ rest frame,
                    and the azimuthal angle $\phi$ between the two decay planes.)",
                    Options{ { "V"_ok, "D^*"_ov } },
                    "B->D^*lnu::UnnormalizedPDF(q2,cos(theta_l),cos(theta_D),phi)",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s,
                        "cos(theta_D)"s,
                        "phi"s
                    ),
                    "B->D^*lnu::NormalizationPDF(q2,cos(theta_l),cos(theta_D),phi)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                // B_s -> K^* l nu
                make_signal_pdf("B_s->K^*lnu::P(q2)",
                    R"(PDF for the decay $\bar{B}_s\to K^*(\to K\pi) \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "U"_ok, "u"_ov }, {"q"_ok, "s"_ov}, { "I"_ok, "1/2"_ov } },
                    "B_s->K^*lnu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B_s->K^*lnu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                make_signal_pdf("B_s->K^*lnu::P(q2,cos(theta_l),cos(theta_K),phi)",
                    R"(PDF for the decay $\bar{B}_s\to K^*(\to K\pi) \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$,
                    the cosine of the angle $\theta_K$ between the $K$ and the negative $\ell^-\bar\nu$ flight direction in the $K^*$ rest frame,
                    the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $\bar{B}_s$ flight direction in the $\ell^-\bar\nu$ rest frame,
                    and the azimuthal angle $\phi$ between the two decay planes.)",
                    Options{ { "U"_ok, "u"_ov }, {"q"_ok, "s"_ov}, { "I"_ok, "1/2"_ov } },
                    "B_s->K^*lnu::UnnormalizedPDF(q2,cos(theta_l),cos(theta_K),phi)",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s,
                        "cos(theta_K)"s,
                        "phi"s
                    ),
                    "B_s->K^*lnu::NormalizationPDF(q2,cos(theta_l),cos(theta_K),phi)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                // B_s -> D_s^* l nu
                make_signal_pdf("B_s->D_s^*lnu::P(q2)",
                    R"(PDF for the decay $\bar{B}_s\to D_s^*(\to D_s\pi^0) \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "U"_ok, "c"_ov }, {"q"_ok, "s"_ov}, { "I"_ok, "0"_ov } },
                    "B_s->D_s^*lnu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B_s->D_s^*lnu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                make_signal_pdf("B_s->D_s^*lnu::P(q2,cos(theta_l),cos(theta_D_s),phi)",
                    R"(PDF for the decay $\bar{B}_s\to D_s^*(\to D_s\pi^0) \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$,
                    the cosine of the angle $\theta_{D_s}$ between the $D_s$ and the negative $\ell^-\bar\nu$ flight direction in the $D_s^*$ rest frame,
                    the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $\bar{B}_s$ flight direction in the $\ell^-\bar\nu$ rest frame,
                    and the azimuthal angle $\phi$ between the two decay planes.)",
                    Options{ { "U"_ok, "c"_ov }, {"q"_ok, "s"_ov}, { "I"_ok, "0"_ov } },
                    "B_s->D_s^*lnu::UnnormalizedPDF(q2,cos(theta_l),cos(theta_D_s),phi)",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s,
                        "cos(theta_D_s)"s,
                        "phi"s
                    ),
                    "B_s->D_s^*lnu::NormalizationPDF(q2,cos(theta_l),cos(theta_D_s),phi)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),
            }
        );

        return SignalPDFGroup(imp);
    }
    // }}}

    // Semileptonic B -> P(seudoscalar) P(seudoscalar) decays
    // {{{
    SignalPDFGroup
    make_b_to_p_p_l_nu_pdf_group()
    {
        auto imp = new Implementation<SignalPDFGroup>(
            R"(Signal PDFs in semileptonic $B\to PP \ell^-\bar\nu$ decays)",
            R"()",
            {
                make_signal_pdf("B^+->pi^+pi^-lnu::PDF(q2,k2,cos(theta_pi))",
                    R"(PDF for the decay $B^+\to \pi^+\pi^- \ell^+\nu$ as a function of the invariant $\pi^+\pi^-$ mass squared $k^2$,
                    the invariant $\ell^+\nu$ mass squared $q^2$, and the cosine of the angle $\theta_\pi$ between the $\pi^+$ and the
                    negative $B^+$ flight direction in the $\pi^+\pi^-$ rest frame.)",
                    Options{ { "U"_ok, "u"_ov }, {"q"_ok, "u"_ov}, {"I1"_ok, "1"_ov}, {"I2"_ok, "1"_ov}, {"C"_ok, "+-"_ov}, {"I"_ok, "0|1"_ov}, {"L"_ok, "S|P|D"_ov} },
                    "B^+->pi^+pi^-lnu::UnnormalizedPDF(q2,k2,cos(theta_pi))",
                    std::make_tuple(
                        "q2"s,
                        "k2"s,
                        "cos(theta_pi)"s
                    ),
                    "B^+->pi^+pi^-lnu::NormalizationPDF(q2,k2,cos(theta_pi))",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s,
                        "k2_min"s,
                        "k2_max"s,
                        "cos(theta)_min"s,
                        "cos(theta)_max"s
                    )
                ),
            }
        );

        return SignalPDFGroup(imp);
    }
    // }}}

    // Semileptonic Lambda_b -> 1/2^+ decays
    // {{{
    SignalPDFGroup
    make_lambdab_to_onehalfplus_l_nu_group()
    {
        auto imp = new Implementation<SignalPDFGroup>(
            R"(Signal PDFs in semileptonic $\Lambda_b\to 1/2^+ \ell^-\bar\nu$ decays)",
            R"()",
            {
                make_signal_pdf("Lambda_b->Lambda_clnu::P(q2,cos(theta_l),cos(theta_L),phi)",
                    R"(PDF for the decay $\Lambda_b\to \Lambda_c(\to \Lambda \pi) \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$,
                    the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $\Lambda_b$ flight direction in the $\ell^-\bar\nu$ rest frame,
                    the cosine of the angle $\theta_L$ between the $\Lambda$ and the negative $\Lambda_b$ flight direction in the $\Lambda_c$ rest frame,
                    and the azimuthal angle $\phi$ between the two decay planes.)",
                    Options{ },
                    "Lambda_b->Lambda_clnu::UnnormalizedPDF(q2,cos(theta_l),cos(theta_L),phi)",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s,
                        "cos(theta_L)"s,
                        "phi"s
                    ),
                    "Lambda_b->Lambda_clnu::NormalizationPDF(q2,cos(theta_l),cos(theta_L),phi)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),
            }
        );

        return SignalPDFGroup(imp);
    }
    // }}}

    // Semileptonic Lambda_b -> 3/2^- decays
    // {{{
    SignalPDFGroup
    make_lambdab_to_threehalfminus_l_nu_group()
    {
        auto imp = new Implementation<SignalPDFGroup>(
            R"(Signal PDFs in semileptonic $\Lambda_b\to 3/2^- \ell^-\bar\nu$ decays)",
            R"()",
            {
                // Lambda_b -> Lambda_c(2625) l nu
                make_signal_pdf("Lambda_b->Lambda_c(2625)lnu::P(q2,cos(theta_l))",
                    R"(PDF for the decay $\Lambda_b\to \Lambda_c(2625) \ell^-\bar\nu$ as a function of the invariant dilepton mass squared $q^2$
                    and the cosine of the angle $\theta_\ell$ between the charged lepton and the negative $\Lambda_b$ flight direction in the $\ell^-\bar\nu$ rest frame.)",
                    Options{ },
                    "Lambda_b->Lambda_c(2625)lnu::UnnormalizedPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s
                    ),
                    "Lambda_b->Lambda_c(2625)lnu::NormalizationPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),
            }
        );

        return SignalPDFGroup(imp);
    }
    // }}}

    SignalPDFSection
    make_b_decays_pdf_section()
    {
        auto imp = new Implementation<SignalPDFSection>(
            "Signal PDFs in (semi)leptonic $b$-hadron decays",
            "",
            {
                // Leptonic and photoleptonic B decays
                make_b_to_leptons_pdf_group(),

                // Semileptonic B_{u,d,s} -> P l^- nubar decays
                make_b_to_p_l_nu_pdf_group(),

                // Semileptonic B_{u,d,s} -> V l^- nubar decays
                make_b_to_v_l_nu_pdf_group(),

                // Semileptonic B_{u,d,s} -> P P l^- nubar decays
                make_b_to_p_p_l_nu_pdf_group(),

                // Semileptonic Lambda_b -> 1/2^+ l^- nubar decays
                make_lambdab_to_onehalfplus_l_nu_group(),

                // Semileptonic Lambda_b -> 3/2^- l^- nubar decays
                make_lambdab_to_threehalfminus_l_nu_group(),
            }
        );

        return SignalPDFSection(imp);
    }
}
