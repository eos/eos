/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2023-2025 Danny van Dyk
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
                make_signal_pdf("B^-->K^-nunu::dGamma/dq2",
                    Options{ { "q"_ok, "u" }, { "I"_ok, "1/2" }, { "D"_ok, "s" } },
                    &BToPseudoscalarDineutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 22.90, BToPseudoscalarDineutrino::kinematics_description_q2 }
                    ),
                    &BToPseudoscalarDineutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

                // B -> K l^+ l^-
                make_signal_pdf("B->Kll::d^2Gamma@LargeRecoil",
                    Options{ {"tag"_ok, "BFS2004"} },
                    &BToKDilepton::two_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "s",                  1.00,  6.00, BToKDilepton::kinematics_description_s },
                        KinematicRange{ "cos(theta_l)^LHCb", -1.0,  +1.0,  BToKDilepton::kinematics_description_c_theta_l }
                    ),
                    &BToKDilepton::integrated_decay_width,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

                make_signal_pdf("B->Kll::d^2Gamma@LowRecoil",
                    Options{ {"tag"_ok, "GP2004" }},
                    &BToKDilepton::two_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "s",                 15.00, 22.87, BToKDilepton::kinematics_description_s },
                        KinematicRange{ "cos(theta_l)^LHCb", -1.0,  +1.0,  BToKDilepton::kinematics_description_c_theta_l }
                    ),
                    &BToKDilepton::integrated_decay_width,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),
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
                make_signal_pdf("B^-->K^*-nunu::dGamma/dq2",
                    Options{ { "q"_ok, "u" }, { "I"_ok, "1/2" }, { "D"_ok, "s" } },
                    &BToVectorDineutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 19.25, BToVectorDineutrino::kinematics_description_q2 }
                    ),
                    &BToVectorDineutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

                // B -> K^* l^+ l^-
                make_signal_pdf("B->K^*ll::d^4Gamma@LargeRecoil",
                    Options{ {"tag"_ok, "BFS2004"} },
                    &BToKstarDilepton::decay_width_LHCb,
                    std::make_tuple(
                        KinematicRange{ "s",                  1.00,  6.00,       BToKstarDilepton::kinematics_description_s         },
                        KinematicRange{ "cos(theta_l)^LHCb", -1.0,  +1.0,        BToKstarDilepton::kinematics_description_c_theta_l },
                        KinematicRange{ "cos(theta_k)^LHCb", -1.0,  +1.0,        BToKstarDilepton::kinematics_description_c_theta_k },
                        KinematicRange{ "phi^LHCb",           0.0,   2.0 * M_PI, BToKstarDilepton::kinematics_description_phi       }
                    ),
                    std::function<double (const BToKstarDilepton *, const double &, const double &)>([] (const BToKstarDilepton * decay, const double & q2_min, const double & q2_max) -> double {
                        return decay->integrated_decay_width(decay->prepare(q2_min, q2_max));
                    }),
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

                make_signal_pdf("B->K^*ll::d^4Gamma@LowRecoil",
                    Options{ {"tag"_ok, "GP2004" }},
                    &BToKstarDilepton::decay_width_LHCb,
                    std::make_tuple(
                        KinematicRange{ "s",                  15.00,  19.21,      BToKstarDilepton::kinematics_description_s },
                        KinematicRange{ "cos(theta_l)^LHCb", -1.0,   +1.0,        BToKstarDilepton::kinematics_description_c_theta_l },
                        KinematicRange{ "cos(theta_k)^LHCb", -1.0,   +1.0,        BToKstarDilepton::kinematics_description_c_theta_k },
                        KinematicRange{ "phi^LHCb",           0.0,    2.0 * M_PI, BToKstarDilepton::kinematics_description_phi }
                    ),
                    std::function<double (const BToKstarDilepton *, const double &, const double &)>([] (const BToKstarDilepton * decay, const double & q2_min, const double & q2_max) -> double {
                        return decay->integrated_decay_width(decay->prepare(q2_min, q2_max));
                    }),
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),
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
