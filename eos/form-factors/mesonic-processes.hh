/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020-2023 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_MESONIC_PROCESSES_HH
#define EOS_GUARD_EOS_FORM_FACTORS_MESONIC_PROCESSES_HH 1

namespace eos
{
    /* P -> P Processes */

    struct BToK {
        using Transition = PToP;
        static constexpr const char * label = "B->K";
        static constexpr const char * name_B = "mass::B_d";
        static constexpr const char * name_P = "mass::K_d";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange);
        static constexpr const double m_B = 5.279;
        static constexpr const double m_P = 0.492;
        // first resonances sorted by spin/parity
        static constexpr double mR2_0p = 5.630 * 5.630; // B_s scalar
        static constexpr double mR2_1m = 5.415 * 5.415; // B_s^*
        // Isospin-degeneracy factor
        static constexpr double eta  = 2.0;
        // OPE results for the unitarity bounds
        static constexpr double chi_0p_v  = 1.42e-2;
        static constexpr double chi_1m_v  = 1.20e-2 / (4.2 * 4.2);
        static constexpr double chi_1m_t  = 0.803e-2 / (4.2 * 4.2);

        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct BToPi {
        using Transition = PToP;
        static constexpr const char * label = "B->pi";
        static constexpr const char * name_B = "mass::B_d";
        static constexpr const char * name_P = "mass::pi^0";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down);
        static constexpr const double m_B = 5.279;
        static constexpr const double m_P = 0.135;
        static constexpr const double mR2_1m = 5.325 * 5.325; // B_{u,d}^*
        static constexpr const double mR2_0p = 5.540 * 5.540; // B_{u,d} scalar: M(B_s scalar) - M(B_s^*) + M(B_{u,d}^*)
        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct BsToK {
        using Transition = PToP;
        static constexpr const char * label = "B_s->K";
        static constexpr const char * name_B = "mass::B_s";
        static constexpr const char * name_P = "mass::K_d";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up);
        static constexpr const double m_B = 5.366;
        static constexpr const double m_P = 0.494;
        static constexpr const double mR2_1m = 5.325 * 5.325; // B_{u,d}^*
        static constexpr const double mR2_0p = 5.540 * 5.540; // B_{u,d} scalar: M(B_s scalar) - M(B_s^*) + M(B_{u,d}^*)
        // Isospin-degeneracy factor
        static constexpr double eta  = 1.0;
        // OPE results for the unitarity bounds --> isospin symmetry, using results for d quarks
        static constexpr double chi_0p_v  = 1.50e-2;
        static constexpr double chi_1m_v  = 1.16e-2 / (4.2 * 4.2);
        static constexpr double chi_1m_t  = 7.75e-3 / (4.2 * 4.2); // already divided by 4 due to different convention in [BFW:2010A]
        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct BToD {
        using Transition = PToP;
        static constexpr const char * label = "B->D";
        static constexpr const char * name_B = "mass::B_d";
        static constexpr const char * name_P = "mass::D_u";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm);
        static constexpr const double m_B = 5.279;
        static constexpr const double m_P = 1.870;
        // resonance masses from [HPQCD2015A]
        static constexpr const double mR2_1m = 6.330 * 6.330; // B_c^*
        static constexpr const double mR2_0p = 6.420 * 6.420; // B_c scalar
        static constexpr const bool uses_tensor_form_factors = true;
        static constexpr const char * hqe_prefix = "B(*)->D(*)";
    };

    struct BsToDs {
        using Transition = PToP;
        static constexpr const char * label = "B_s->D_s";
        static constexpr const char * name_B = "mass::B_s";
        static constexpr const char * name_P = "mass::D_s";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm);
        static constexpr const double m_B = 5.366;
        static constexpr const double m_P = 1.968;
        // resonance masses from [HPQCD2015A]
        static constexpr const double mR2_1m = 6.330 * 6.330; // B_c^*
        static constexpr const double mR2_0p = 6.420 * 6.420; // B_c scalar
        static constexpr const double tp = (m_B + m_P) * (m_B + m_P);
        static constexpr const bool uses_tensor_form_factors = true;
        static constexpr const char * hqe_prefix = "B_s(*)->D_s(*)";
    };

    struct DToPi {
        using Transition = PToP;
        static constexpr const char * label = "D->pi";
        static constexpr const char * name_B = "mass::D_u";
        static constexpr const char * name_P = "mass::pi^0";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::charm, QuarkFlavor::up);
        static constexpr const double m_B = 1.867;
        static constexpr const double m_P = 0.135;
        static constexpr const double mR2_1m = 2.007 * 2.007; // D^*0
        static constexpr const double mR2_0p = 2.300 * 2.300; // D^*0 scalar
        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct DToK {
        using Transition = PToP;
        static constexpr const char * label = "D->K";
        static constexpr const char * name_B = "mass::D_u";
        static constexpr const char * name_P = "mass::K_u";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::charm, QuarkFlavor::strange);
        static constexpr const double m_B = 1.867;
        static constexpr const double m_P = 0.492;
        static constexpr const double mR2_1m = 2.714 * 2.714; // Ds1
        static constexpr const double mR2_0p = 2.317 * 2.317; // Ds0
        static constexpr const bool uses_tensor_form_factors = true;
        // zero of the conformal mapping: z(t0, t0) = 0.0
        // This optimal value follows from z(0, t0) = - z(tm, t0)
        static constexpr double t0 = 1.04;
        // Isospin-degeneracy factor
        static constexpr double eta  = 2.0;
        // OPE results for the unitarity bounds (1103.1481)
        static constexpr double chi_0p_v  = 1.38e-2;
        static constexpr double chi_1m_v  = 9.35e-3;
        static constexpr double chi_1m_t  = 6.89e-03;
    };

    struct DsToK {
        using Transition = PToP;
        static constexpr const char * label = "D_s->K";
        static constexpr const char * name_B = "mass::D_s";
        static constexpr const char * name_P = "mass::K_u";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::charm, QuarkFlavor::up);
        static constexpr const double m_B = 1.968;
        static constexpr const double m_P = 0.492;
        static constexpr const double mR2_1m = 2.007 * 2.007; // D^*0
        static constexpr const double mR2_0p = 2.300 * 2.300; // D^*0 scalar
        static constexpr const bool uses_tensor_form_factors = true;
    };

    /* P -> V Processes */

    struct BToDstar {
        using Transition = PToV;
        static constexpr const char * label = "B->D^*";
        static constexpr const char * name_B = "mass::B_d";
        static constexpr const char * name_V = "mass::D_u^*";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm);
        static constexpr double m_B = 5.279;
        static constexpr double m_V = 2.0103;
        static constexpr double m_Bc = 6.2751;
        static constexpr double mR2_0m = (m_Bc + 0.000) * (m_Bc + 0.000);
        static constexpr double mR2_1m = (m_Bc + 0.056) * (m_Bc + 0.056);
        static constexpr double mR2_1p = (m_Bc + 0.492) * (m_Bc + 0.492);
        static constexpr const char * hqe_prefix = "B(*)->D(*)";
    };

    struct BsToDsstar {
        using Transition = PToV;
        static constexpr const char * label = "B_s->D_s^*";
        static constexpr const char * name_B = "mass::B_s";
        static constexpr const char * name_V = "mass::D_s^*";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm);
        static constexpr double m_B = 5.366;
        static constexpr double m_V = 2.1121;
        static constexpr double m_Bc = 6.2751;
        static constexpr double mR2_0m = (m_Bc + 0.000) * (m_Bc + 0.000);
        static constexpr double mR2_1m = (m_Bc + 0.056) * (m_Bc + 0.056);
        static constexpr double mR2_1p = (m_Bc + 0.492) * (m_Bc + 0.492);
        static constexpr const char * hqe_prefix = "B_s(*)->D_s(*)";
    };

    struct BToKstar {
        using Transition = PToV;
        static constexpr const char * label = "B->K^*";
        static constexpr const char * name_B = "mass::B_d";
        static constexpr const char * name_V = "mass::K_d^*";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange);
        static constexpr double m_B = 5.279;
        static constexpr double m_V = 0.896;
        // // first resonances sorted by spin/parity
        static constexpr double mR2_0m = 5.366 * 5.366; // B_s
        static constexpr double mR2_0p = 5.630 * 5.630; // B_s scalar
        static constexpr double mR2_1m = 5.415 * 5.415; // B_s^*
        static constexpr double mR2_1p = 5.829 * 5.829; // B_s,1
        // scalar pair production threshold: B + K
        static constexpr double tp_v = (5.279 + 0.492) * (5.279 + 0.492);
        // vector pair production threshold: B + K + pi
        static constexpr double tp_a = (5.279 + 0.492 + 0.135) * (5.279 + 0.492 + 0.135);
        // Isospin-degeneracy factor
        static constexpr double eta  = 2.0;
        // OPE results for the unitarity bounds
        static constexpr double chi_0m_a  = 1.57e-2;
        static constexpr double chi_0p_v  = 1.42e-2;
        static constexpr double chi_1m_v  = 1.20e-2 / (4.2 * 4.2);
        static constexpr double chi_1p_a  = 1.13e-2 / (4.2 * 4.2);
        static constexpr double chi_1m_t  = 0.803e-2 / (4.2 * 4.2);
        static constexpr double chi_1p_t5 = 0.748e-2 / (4.2 * 4.2);
    };

    struct BToOmega {
        using Transition = PToV;
        static constexpr const char * label = "B->omega";
        static constexpr const char * name_B = "mass::B_d";
        static constexpr const char * name_V = "mass::omega";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down);
        static constexpr double m_B = 5.279;
        static constexpr double m_V = 0.7827;
        static constexpr double mR2_0m = 5.279 * 5.279;
        static constexpr double mR2_1m = 5.325 * 5.325;
        static constexpr double mR2_1p = 5.724 * 5.724;
    };

    struct BToRho {
        using Transition = PToV;
        static constexpr const char * label = "B->rho";
        static constexpr const char * name_B = "mass::B_d";
        static constexpr const char * name_V = "mass::rho^0";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down);
        static constexpr double m_B = 5.279;
        static constexpr double m_V = 0.7751;
        static constexpr double mR2_0m = 5.279 * 5.279;
        static constexpr double mR2_1m = 5.325 * 5.325;
        static constexpr double mR2_1p = 5.724 * 5.724;
    };

    struct BsToPhi {
        using Transition = PToV;
        static constexpr const char * label = "B_s->phi";
        static constexpr const char * name_B = "mass::B_s";
        static constexpr const char * name_V = "mass::phi";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange);
        static constexpr double m_B = 5.366;
        static constexpr double m_V = 1.020;
        // // first resonances sorted by spin/parity
        static constexpr double mR2_0m = 5.366 * 5.366; // B_s
        static constexpr double mR2_0p = 5.630 * 5.630; // B_s scalar
        static constexpr double mR2_1m = 5.415 * 5.415; // B_s^*
        static constexpr double mR2_1p = 5.829 * 5.829; // B_s,1
        // scalar pair production threshold: B + K
        static constexpr double tp_v = (5.279 + 0.492) * (5.279 + 0.492);
        // vector pair production threshold: B + K + pi
        static constexpr double tp_a = (5.279 + 0.492 + 0.135) * (5.279 + 0.492 + 0.135);
        // Isospin-degeneracy factor
        static constexpr double eta  = 1.0;
        // OPE results for the unitarity bounds
        static constexpr double chi_0m_a  = 1.57e-2;
        static constexpr double chi_0p_v  = 1.42e-2;
        static constexpr double chi_1m_v  = 1.20e-2 / (4.2 * 4.2);
        static constexpr double chi_1p_a  = 1.13e-2 / (4.2 * 4.2);
        static constexpr double chi_1m_t  = 0.803e-2 / (4.2 * 4.2);
        static constexpr double chi_1p_t5 = 0.748e-2 / (4.2 * 4.2);
    };

    struct BsToKstar {
        using Transition = PToV;
        static constexpr const char * label = "B_s->K^*";
        static constexpr const char * name_B = "mass::B_s";
        static constexpr const char * name_V = "mass::K_d^*";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down);
        static constexpr double m_B = 5.366;
        static constexpr double m_V = 0.896;
        static constexpr double mR2_0m = 5.279 * 5.279;
        static constexpr double mR2_1m = 5.325 * 5.325;
        static constexpr double mR2_1p = 5.723 * 5.723;
    };

    /* P -> PP Processes */

    struct BToPiPi {
        using Transition = PToPP;
        static constexpr const char * label = "B->pipi";
        static constexpr double m_B  = 5.2795;
        static constexpr double m_P1 = 0.13957;
        static constexpr double m_P2 = 0.13957;

        // for pole and t_0 calculation in zhat
        static constexpr double m_Bst = 5.32465;

        // for pole calculation in z, depending on the current at hand
        static constexpr double mR2_1m = 5.32465;
        static constexpr double mR2_1p = 5.72590;
        static constexpr double mR2_0m = 5.27932;
    };

    struct BToKPi {
        using Transition = PToPP;
        static constexpr const char * label = "B->Kpi";
        static constexpr const char * name_B = "mass::B";
        static constexpr const char * name_P1 = "mass::K";
        static constexpr const char * name_P2 = "mass::Pi";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down);

        // for pole and t_0 calculation in zhat
        static constexpr double m_Bst = 5.32465;

        // for pole calculation in z, depending on the current at hand
        static constexpr double mR2_1m = 5.32465;
        static constexpr double mR2_1p = 5.72590;
        static constexpr double mR2_0m = 5.27932;
    };

    /* V -> P Processes */

    struct BstarToD {
        using Transition = VToP;
        static constexpr const char * label    = "B^*->D";
        static constexpr const char * name_Bst = "mass::B_d^*";
        static constexpr const char * name_P   = "mass::D_u";
        static constexpr const double m_Bc  = 6.2751;
        static constexpr const double mR2_0m = (m_Bc + 0.000) * (m_Bc + 0.000);
        static constexpr const double mR2_1m = (m_Bc + 0.056) * (m_Bc + 0.056);
        static constexpr const double mR2_1p = (m_Bc + 0.492) * (m_Bc + 0.492);
        static constexpr const char * hqe_prefix = "B(*)->D(*)";
    };

    struct BsstarToDs {
        using Transition = VToP;
        static constexpr const char * label = "B_s^*->D_s";
        static constexpr const char * name_Bst = "mass::B_s^*";
        static constexpr const char * name_P   = "mass::D_s";
        static constexpr const double m_B = 5.324;
        static constexpr const double m_P = 1.968;
        static constexpr const double m_Bc  = 6.2751;
        static constexpr const double mR2_0m = (m_Bc + 0.000) * (m_Bc + 0.000);
        static constexpr const double mR2_1m = (m_Bc + 0.056) * (m_Bc + 0.056);
        static constexpr const double mR2_1p = (m_Bc + 0.492) * (m_Bc + 0.492);
        static constexpr const char * hqe_prefix = "B_s(*)->D_s(*)";
    };

    /* V -> V Processes */

    struct BstarToDstar {
        using Transition = VToV;
        static constexpr const char * label = "B^*->D^*";
        static constexpr const char * name_Bst = "mass::B_d^*";
        static constexpr const char * name_V   = "mass::D_u^*";
        static constexpr const double m_V1 = 5.324;
        static constexpr const double m_V2 = 2.010;
        static constexpr const double m_Bc = 6.2751;
        static constexpr const double mR2_0m = (m_Bc + 0.000) * (m_Bc + 0.000);
        static constexpr const double mR2_1m = (m_Bc + 0.056) * (m_Bc + 0.056);
        static constexpr const double mR2_1p = (m_Bc + 0.492) * (m_Bc + 0.492);
        static constexpr const char * hqe_prefix = "B(*)->D(*)";
    };

    struct BsstarToDsstar {
        using Transition = VToV;
        static constexpr const char * label = "B_s^*->D_s^*";
        static constexpr const char * name_Bst = "mass::B_s^*";
        static constexpr const char * name_P   = "mass::D_s";
        static constexpr const double m_B = 5.324;
        static constexpr const double m_V = 2.010;
        static constexpr const double m_Bc  = 6.2751;
        static constexpr const double mR2_0m = (m_Bc + 0.000) * (m_Bc + 0.000);
        static constexpr const double mR2_1m = (m_Bc + 0.056) * (m_Bc + 0.056);
        static constexpr const double mR2_1p = (m_Bc + 0.492) * (m_Bc + 0.492);
        static constexpr const char * hqe_prefix = "B_s(*)->D_s(*)";
    };
}

#endif
