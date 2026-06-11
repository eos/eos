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
#include <eos/rare-b-decays/b-to-k-ll.hh>
#include <eos/rare-b-decays/b-to-kstar-ll.hh>
#include <eos/rare-b-decays/b-to-psd-nu-nu.hh>
#include <eos/rare-b-decays/b-to-vec-nu-nu.hh>
#include <eos/utils/concrete-signal-pdf.hh>

using std::literals::string_literals::operator""s;

namespace eos
{
    // Rare semileptonic B -> P(seudoscalar) decays
    // {{{
    SignalPDFGroup
    make_b_to_p_l_l_pdf_group()
    {
        auto imp = new Implementation<SignalPDFGroup>(
            R"(Signal PDFs in rare semileptonic $B\to P \lbrace \bar\nu\nu, \ell^+\ell^- \rbrace$ decays)",
            R"()",
            {
                // B -> K nubar nu
                make_signal_pdf("B^-->K^-nunu::P(q2)",
                    R"(PDF for the decay $\bar{B}\to \bar{K} \bar\nu\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "q"_ok, "u"_ov }, { "I"_ok, "1/2"_ov }, { "D"_ok, "s"_ov } },
                    "B->Knunu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B->Knunu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                // B -> K l^+ l^-
                make_signal_pdf("B->Kll::P(q2,cos(theta_l))",
                    R"(PDF for the decay $\bar{B}\to \bar{K} \ell^+\ell^-$ as a function of the invariant dilepton mass squared $q^2$
                    and the cosine of the angle $\theta_\ell$ between the positively charged lepton and the negative $B$
                    flight direction in the $\ell^+\ell^-$ rest frame.)",
                    Options{ {"tag"_ok, "BFS2004"_ov} },
                    "B->Kll::UnnormalizedPDF(q2,cos(theta_l))",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s
                    ),
                    "B->Kll::NormalizationPDF(q2,cos(theta_l))",
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

    // Rare semileptonic B -> V(pseudoscalar) decays
    // {{{
    SignalPDFGroup
    make_b_to_v_l_l_pdf_group()
    {
        auto imp = new Implementation<SignalPDFGroup>(
            R"(Signal PDFs in rare semileptonic $B\to V \lbrace \bar\nu\nu, \ell^+\ell^- \rbrace$ decays)",
            R"()",
            {
                // B -> K^* nu nu
                make_signal_pdf("B^-->K^*-nunu::P(q2)",
                    R"(PDF for the decay $\bar{B}\to \bar{K}^*(\to \bar{K}\pi) \bar\nu\nu$ as a function of the invariant dilepton mass squared $q^2$.)",
                    Options{ { "q"_ok, "u"_ov }, { "I"_ok, "1/2"_ov }, { "D"_ok, "s"_ov } },
                    "B->K^*nunu::UnnormalizedPDF(q2)",
                    std::make_tuple(
                        "q2"s
                    ),
                    "B->K^*nunu::NormalizationPDF(q2)",
                    std::make_tuple(
                        "q2_min"s,
                        "q2_max"s
                    )
                ),

                // B -> K^* l^+ l^-
                make_signal_pdf("B->K^*ll::P(q2,cos(theta_l),cos(theta_K),phi)",
                    R"(PDF for the decay $\bar{B}\to \bar{K}^*(\to \bar{K}\pi) \ell^+\ell^-$ as a function of the invariant dilepton mass squared $q^2$,
                    the cosine of the angle $\theta_K$ between the $\bar{K}$ and the negative $\ell^+\ell^-$ flight direction,
                    the cosine of the angle $\theta_\ell$ between the positively charged lepton and the negative $B$ flight direction in the $\ell^+\ell^-$ rest frame,
                    and the azimuthal angle $\phi$ between the two decay planes.)",
                    Options{ {"tag"_ok, "BFS2004"_ov} },
                    "B->K^*ll::UnnormalizedPDF(q2,cos(theta_l),cos(theta_K),phi)",
                    std::make_tuple(
                        "q2"s,
                        "cos(theta_l)"s,
                        "cos(theta_K)"s,
                        "phi"s
                    ),
                    "B->K^*ll::NormalizationPDF(q2,cos(theta_l),cos(theta_K),phi)",
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

    SignalPDFSection
    make_rare_b_decays_pdf_section()
    {
        auto imp = new Implementation<SignalPDFSection>(
            "Signal PDFs in rare (semi)leptonic $b$-hadron decays",
            "",
            {
                // Rare semileptonic B_{u,d,s} -> P {l^+ l^-, nu nubar} decays
                make_b_to_p_l_l_pdf_group(),

                // Rare semileptonic B_{u,d,s} -> V {l^+ l^-, nu nubar} decays
                make_b_to_v_l_l_pdf_group(),
            }
        );

        return SignalPDFSection(imp);
    }
}
