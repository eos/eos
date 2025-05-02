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
                make_signal_pdf("B->gammalnu::d^2Gamma/dEgamma/dcos(theta_l)",
                    Options{ },
                    &BToGammaLeptonNeutrino::fully_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "Egamma", 0.1, 2.64, BToGammaLeptonNeutrino::kinematics_description_Egamma },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0, BToGammaLeptonNeutrino::kinematics_description_c_theta_l}
                    ),
                    &BToGammaLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "E_gamma_min"
                    )
                ),

                // B -> 3l nu
                make_signal_pdf("B_u->enumumu::d^5Gamma",
                    Options{ { "l"_ok, "e" }, { "lprime"_ok, "mu" } },
                    &BToThreeLeptonsNeutrino::quintuple_differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0447, 27.8714, BToThreeLeptonsNeutrino::kinematics_description_q2 },
                        KinematicRange{ "k2", 0.00051, 25.6849, BToThreeLeptonsNeutrino::kinematics_description_k2 },
                        KinematicRange{ "z_gamma", -1.0, +1.0, BToThreeLeptonsNeutrino::kinematics_description_z_gamma },
                        KinematicRange{ "z_w", -1.0, +1.0, BToThreeLeptonsNeutrino::kinematics_description_z_w },
                        KinematicRange{ "phi", -M_PI, +M_PI, BToThreeLeptonsNeutrino::kinematics_description_phi }
                    ),
                    &BToThreeLeptonsNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max",
                        "k2_min",
                        "k2_max"
                    )
                ),

                make_signal_pdf("B_u->munuee::d^5Gamma",
                    Options{ { "l"_ok, "mu" }, { "lprime"_ok, "e" } },
                    &BToThreeLeptonsNeutrino::quintuple_differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 1.0e-6, 26.767, BToThreeLeptonsNeutrino::kinematics_description_q2 },
                        KinematicRange{ "k2", 0.011, 27.8606, BToThreeLeptonsNeutrino::kinematics_description_k2 },
                        KinematicRange{ "z_gamma", -1.0, +1.0, BToThreeLeptonsNeutrino::kinematics_description_z_gamma },
                        KinematicRange{ "z_w", -1.0, +1.0, BToThreeLeptonsNeutrino::kinematics_description_z_w },
                        KinematicRange{ "phi", -M_PI, +M_PI, BToThreeLeptonsNeutrino::kinematics_description_phi }
                    ),
                    &BToThreeLeptonsNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max",
                        "k2_min",
                        "k2_max"
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
                make_signal_pdf("B->pilnu::dGamma/dq2",
                    Options{ { "P"_ok, "pi" } },
                    &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 26.41, BToPseudoscalarLeptonNeutrino::kinematics_description_q2 }
                    ),
                    &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

                make_signal_pdf("B->pilnu::d^2Gamma/dq2/dcos(theta_l)",
                    Options{ { "P"_ok, "pi" } },
                    &BToPseudoscalarLeptonNeutrino::normalized_two_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 26.41, BToPseudoscalarLeptonNeutrino::kinematics_description_q2 },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0, BToPseudoscalarLeptonNeutrino::kinematics_description_c_theta_l}
                    ),
                    &BToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

                // B -> pi l X_nubar
                make_signal_pdf("B->pimu1nu::d^2Gamma",
                    Options{ },
                    &BToPiLeptonInclusiveNeutrinos::differential_decay_width_1nu,
                    std::make_tuple(
                        KinematicRange{ "s", 0.0, 26.41, BToPiLeptonInclusiveNeutrinos::kinematics_description_s},
                        KinematicRange{ "cos(theta)", -1.0, +1.0, BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta}
                    ),
                    &BToPiLeptonInclusiveNeutrinos::integrated_decay_width_1nu,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

                make_signal_pdf("B->pimu3nu::d^5Gamma",
                    Options{ },
                    &BToPiLeptonInclusiveNeutrinos::differential_decay_width_3nu,
                    std::make_tuple(
                        KinematicRange{ "s", 3.16, 26.41, BToPiLeptonInclusiveNeutrinos::kinematics_description_s },
                        KinematicRange{ "snunubar", 0.0, 3.16, BToPiLeptonInclusiveNeutrinos::kinematics_description_snunubar },
                        KinematicRange{ "cos(theta_tau)", -1.0, +1.0, BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta_tau },
                        KinematicRange{ "phi", 0.0, 2.0 * M_PI, BToPiLeptonInclusiveNeutrinos::kinematics_description_phi },
                        KinematicRange{ "cos(theta_mu^*)", -1.0, +1.0, BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta_mu_star }
                    ),
                    &BToPiLeptonInclusiveNeutrinos::integrated_decay_width_3nu,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

                // B -> D l nu
                make_signal_pdf("B->Dlnu::dGamma/dq2",
                    Options{ { "P"_ok, "D" } },
                    &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 11.62, BToPseudoscalarLeptonNeutrino::kinematics_description_q2 }
                    ),
                    &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

                make_signal_pdf("B->Dlnu::d^2Gamma/dq2/dcos(theta_l)",
                    Options{ { "P"_ok, "D" } },
                    &BToPseudoscalarLeptonNeutrino::normalized_two_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 11.62, BToPseudoscalarLeptonNeutrino::kinematics_description_q2 },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0, BToPseudoscalarLeptonNeutrino::kinematics_description_c_theta_l}
                    ),
                    &BToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

                // B -> D l X_nubar
                make_signal_pdf("B->Dmu1nu::d^2Gamma",
                    Options{ },
                    &BToDLeptonInclusiveNeutrinos::differential_decay_width_1nu,
                    std::make_tuple(
                        KinematicRange{ "s", 0.0, 19.71, BToDLeptonInclusiveNeutrinos::kinematics_description_s },
                        KinematicRange{ "cos(theta)", -1.0, +1.0, BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta}
                    ),
                    &BToDLeptonInclusiveNeutrinos::integrated_decay_width_1nu,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

                make_signal_pdf("B->Dmu3nu::d^5Gamma",
                    Options{ },
                    &BToDLeptonInclusiveNeutrinos::differential_decay_width_3nu,
                    std::make_tuple(
                        KinematicRange{ "s", 3.16, 19.71, BToDLeptonInclusiveNeutrinos::kinematics_description_s },
                        KinematicRange{ "snunubar", 0.0, 3.16, BToDLeptonInclusiveNeutrinos::kinematics_description_snunubar },
                        KinematicRange{ "cos(theta_tau)", -1.0, +1.0, BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta_tau },
                        KinematicRange{ "phi", 0.0, 2.0 * M_PI, BToDLeptonInclusiveNeutrinos::kinematics_description_phi },
                        KinematicRange{ "cos(theta_mu^*)", -1.0, +1.0, BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta_mu_star }
                    ),
                    &BToDLeptonInclusiveNeutrinos::integrated_decay_width_3nu,
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
                make_signal_pdf("B->D^*lnu::dBR",
                    Options{ { "V"_ok, "D^*" } },
                    &BToVectorLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 10.68, BToVectorLeptonNeutrino::kinematics_description_q2 }
                    ),
                    &BToVectorLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

                make_signal_pdf("B->D^*lnu::d^4Gamma",
                    Options{ { "V"_ok, "D^*" } },
                    &BToVectorLeptonNeutrino::normalized_four_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "q2",            0.0,  10.68,      BToVectorLeptonNeutrino::kinematics_description_q2        },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0,        BToVectorLeptonNeutrino::kinematics_description_c_theta_l },
                        KinematicRange{ "cos(theta_d)", -1.0, +1.0,        BToVectorLeptonNeutrino::kinematics_description_c_theta_d },
                        KinematicRange{ "phi",           0.0,  2.0 * M_PI, BToVectorLeptonNeutrino::kinematics_description_phi       }
                    ),
                    &BToVectorLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

                // B_s -> K^* l nu
                make_signal_pdf("B_s->K^*lnu::dBR",
                    Options{ { "U"_ok, "u" }, {"q"_ok, "s"}, { "I"_ok, "1/2" } },
                    &BToVectorLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 10.68, BToVectorLeptonNeutrino::kinematics_description_q2 }
                    ),
                    &BToVectorLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

                make_signal_pdf("B_s->K^*lnu::d^4Gamma",
                    Options{ { "U"_ok, "u" }, {"q"_ok, "s"}, { "I"_ok, "1/2" } },
                    &BToVectorLeptonNeutrino::normalized_four_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "q2",            0.0,  10.68,      BToVectorLeptonNeutrino::kinematics_description_q2        },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0,        BToVectorLeptonNeutrino::kinematics_description_c_theta_l },
                        KinematicRange{ "cos(theta_d)", -1.0, +1.0,        BToVectorLeptonNeutrino::kinematics_description_c_theta_d },
                        KinematicRange{ "phi",           0.0,  2.0 * M_PI, BToVectorLeptonNeutrino::kinematics_description_phi       }
                    ),
                    &BToVectorLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
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
                make_signal_pdf("B->pipimunu::d^3Gamma@QCDF",
                    Options{ },
                    &BToPiPiLeptonNeutrino::triple_differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.01, 0.93859, BToPiPiLeptonNeutrino::kinematics_description_q2 },
                        KinematicRange{ "k2", 18.582, 27.872, BToPiPiLeptonNeutrino::kinematics_description_k2 },
                        KinematicRange{ "cos(theta)", -1.0, +1.0, BToPiPiLeptonNeutrino::kinematics_description_z }
                    ),
                    &BToPiPiLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max",
                        "k2_min",
                        "k2_max",
                        "cos(theta)_min",
                        "cos(theta)_max"
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
                // Lambda_b -> Lambda_c l nu
                make_signal_pdf("Lambda_b->Lambda_clnu::dGamma",
                    Options{ },
                    &LambdaBToLambdaCLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.011, 11.1, LambdaBToLambdaCLeptonNeutrino::kinematics_description_q2 }
                    ),
                    &LambdaBToLambdaCLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

                make_signal_pdf("Lambda_b->Lambda_clnu::d^4Gamma",
                    Options{ },
                    &LambdaBToLambdaCLeptonNeutrino::four_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.011, 11.1, LambdaBToLambdaCLeptonNeutrino::kinematics_description_q2 },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0, LambdaBToLambdaCLeptonNeutrino::kinematics_description_c_theta_l },
                        KinematicRange{ "cos(theta_L)", -1.0, +1.0, LambdaBToLambdaCLeptonNeutrino::kinematics_description_c_theta_L },
                        KinematicRange{ "phi", 0.0, 2.0 * M_PI, LambdaBToLambdaCLeptonNeutrino::kinematics_description_phi }
                    ),
                    &LambdaBToLambdaCLeptonNeutrino::integrated_decay_width,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
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
                make_signal_pdf("Lambda_b->Lambda_c(2625)lnu::dGamma",
                    Options{ },
                    &LambdaBToLambdaC2625LeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "s", 0.011, 8.9478, LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_s }
                    ),
                    &LambdaBToLambdaC2625LeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

                make_signal_pdf("Lambda_b->Lambda_c(2625)lnu::d^2Gamma",
                    Options{ },
                    &LambdaBToLambdaC2625LeptonNeutrino::double_differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "s", 0.011, 8.9478, LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_s },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0, LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_c_theta_l }
                    ),
                    &LambdaBToLambdaC2625LeptonNeutrino::integrated_branching_ratio,
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
