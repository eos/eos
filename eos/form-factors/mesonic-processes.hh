/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Danny van Dyk
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

#ifndef ISSUE239_GUARD_EOS_FORM_FACTORS_MESONIC_PROCESSES_HH
#define ISSUE239_GUARD_EOS_FORM_FACTORS_MESONIC_PROCESSES_HH 1

namespace eos
{
    /* P -> P Processes */

    struct BToK {
        using Transition = PToP;
        static constexpr const char * label = "B->K";
        static constexpr const double m_B = 5.279;
        static constexpr const double m_P = 0.492;
        static constexpr const double m2_Br1m = 5.415 * 5.415; // B_s^*
        static constexpr const double m2_Br0p = 5.630 * 5.630; // B_s scalar
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct BToPi {
        using Transition = PToP;
        static constexpr const char * label = "B->pi";
        static constexpr const double m_B = 5.279;
        static constexpr const double m_P = 0.135;
        static constexpr const double m2_Br1m = 5.325 * 5.325; // B_{u,d}^*
        static constexpr const double m2_Br0p = 5.540 * 5.540; // B_{u,d} scalar: M(B_s scalar) - M(B_s^*) + M(B_{u,d}^*)
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct BsToK {
        using Transition = PToP;
        static constexpr const char * label = "B_s->K";
        static constexpr const double m_B = 5.366;
        static constexpr const double m_P = 0.494;
        static constexpr const double m2_Br1m = 5.325 * 5.325; // B_{u,d}^*
        static constexpr const double m2_Br0p = 5.540 * 5.540; // B_{u,d} scalar: M(B_s scalar) - M(B_s^*) + M(B_{u,d}^*)
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct BToD {
        using Transition = PToP;
        static constexpr const char * label = "B->D";
        static constexpr const char * name_B = "mass::B_d";
        static constexpr const char * name_P = "mass::D_u";
        static constexpr const double m_B = 5.279;
        static constexpr const double m_P = 1.870;
        // resonance masses from [HPQCD2015A]
        static constexpr const double m2_Br1m = 6.330 * 6.330; // B_c^*
        static constexpr const double m2_Br0p = 6.420 * 6.420; // B_c scalar
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = true;
        static constexpr const char * hqe_prefix = "B(*)->D(*)";
    };

    struct BsToDs {
        using Transition = PToP;
        static constexpr const char * label = "B_s->D_s";
        static constexpr const char * name_B = "mass::B_s";
        static constexpr const char * name_P = "mass::D_s";
        static constexpr const double m_B = 5.366;
        static constexpr const double m_P = 1.968;
        // resonance masses from [HPQCD2015A]
        static constexpr const double m2_Br1m = 6.330 * 6.330; // B_c^*
        static constexpr const double m2_Br0p = 6.420 * 6.420; // B_c scalar
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = true;
        static constexpr const char * hqe_prefix = "B_s(*)->D_s(*)";
    };

    struct DToPi {
        using Transition = PToP;
        static constexpr const char * label = "D->pi";
        static constexpr const double m_B = 1.867;
        static constexpr const double m_P = 0.135;
        static constexpr const double m2_Br1m = 2.007 * 2.007; // D^*0
        static constexpr const double m2_Br0p = 2.300 * 2.300; // D^*0 scalar
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct DToK {
        using Transition = PToP;
        static constexpr const char * label = "D->K";
        static constexpr const double m_B = 1.867;
        static constexpr const double m_P = 0.492;
        static constexpr const double m2_Br1m = 2.714 * 2.714; // Ds1
        static constexpr const double m2_Br0p = 2.317 * 2.317; // Ds0
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct DsToK {
        using Transition = PToP;
        static constexpr const char * label = "D_s->K";
        static constexpr const double m_B = 1.968;
        static constexpr const double m_P = 0.492;
        static constexpr const double m2_Br1m = 2.007 * 2.007; // D^*0
        static constexpr const double m2_Br0p = 2.300 * 2.300; // D^*0 scalar
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = true;
    };


    /* P -> V Processes */

    struct BToDstar {
        using Transition = PToV;
        static constexpr const char * label = "B->D^*";
        static constexpr const char * name_B = "mass::B_d";
        static constexpr const char * name_V = "mass::D_u^*";
        static constexpr double mB = 5.279;
        static constexpr double mV = 2.0103;
        static constexpr double mBc = 6.2751;
        static constexpr double mR2_0m = (mBc + 0.000) * (mBc + 0.000);
        static constexpr double mR2_1m = (mBc + 0.056) * (mBc + 0.056);
        static constexpr double mR2_1p = (mBc + 0.492) * (mBc + 0.492);
        static constexpr const char * hqe_prefix = "B(*)->D(*)";
    };

    struct BsToDsstar {
        using Transition = PToV;
        static constexpr const char * label = "B_s->D_s^*";
        static constexpr const char * name_B = "mass::B_s";
        static constexpr const char * name_V = "mass::D_s^*";
        static constexpr double mB = 5.366;
        static constexpr double mV = 2.1121;
        static constexpr double mBc = 6.2751;
        static constexpr double mR2_0m = (mBc + 0.000) * (mBc + 0.000);
        static constexpr double mR2_1m = (mBc + 0.056) * (mBc + 0.056);
        static constexpr double mR2_1p = (mBc + 0.492) * (mBc + 0.492);
        static constexpr const char * hqe_prefix = "B_s(*)->D_s(*)";
    };

    struct BToKstar {
        using Transition = PToV;
        static constexpr const char * label = "B->K^*";
        static constexpr double mB = 5.279;
        static constexpr double mV = 0.896;
        static constexpr double mR2_0m = 5.366 * 5.366;
        static constexpr double mR2_1m = 5.415 * 5.415;
        static constexpr double mR2_1p = 5.829 * 5.829;
    };

    struct BToOmega {
        using Transition = PToV;
        static constexpr const char * label = "B->omega";
        static constexpr double mB = 5.279;
        static constexpr double mV = 0.7827;
        static constexpr double mR2_0m = 5.279 * 5.279;
        static constexpr double mR2_1m = 5.325 * 5.325;
        static constexpr double mR2_1p = 5.724 * 5.724;
    };

    struct BToRho {
        using Transition = PToV;
        static constexpr const char * label = "B->rho";
        static constexpr double mB = 5.279;
        static constexpr double mV = 0.7751;
        static constexpr double mR2_0m = 5.279 * 5.279;
        static constexpr double mR2_1m = 5.325 * 5.325;
        static constexpr double mR2_1p = 5.724 * 5.724;
    };

    struct BsToPhi {
        using Transition = PToV;
        static constexpr const char * label = "B_s->phi";
        static constexpr double mB = 5.366;
        static constexpr double mV = 1.020;
        static constexpr double mR2_0m = 5.366 * 5.366;
        static constexpr double mR2_1m = 5.415 * 5.415;
        static constexpr double mR2_1p = 5.829 * 5.829;
    };

    struct BsToKstar {
        using Transition = PToV;
        static constexpr const char * label = "B_s->K^*";
        static constexpr double mB = 5.366;
        static constexpr double mV = 0.896;
        static constexpr double mR2_0m = 5.279 * 5.279;
        static constexpr double mR2_1m = 5.325 * 5.325;
        static constexpr double mR2_1p = 5.723 * 5.723;
    };

    /* P -> PP Processes */

    struct BToPiPi {
        using Transition = PToPP;
        static constexpr const char * label = "B->pipi";
        static constexpr double mB  = 5.2795;
        static constexpr double mP1 = 0.13957;
        static constexpr double mP2 = 0.13957;

        // for pole and t_0 calculation in zhat
        static constexpr double mBst = 5.32465;

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
        static constexpr const double mBc  = 6.2751;
        static constexpr const double mR2_0m = (mBc + 0.000) * (mBc + 0.000);
        static constexpr const double mR2_1m = (mBc + 0.056) * (mBc + 0.056);
        static constexpr const double mR2_1p = (mBc + 0.492) * (mBc + 0.492);
        static constexpr const char * hqe_prefix = "B(*)->D(*)";
    };

    struct BsstarToDs {
        using Transition = VToP;
        static constexpr const char * label = "B_s^*->D_s";
        static constexpr const char * name_Bst = "mass::B_s^*";
        static constexpr const char * name_P   = "mass::D_s";
        static constexpr const double m_B = 5.324;
        static constexpr const double m_P = 1.968;
        static constexpr const double mBc  = 6.2751;
        static constexpr const double mR2_0m = (mBc + 0.000) * (mBc + 0.000);
        static constexpr const double mR2_1m = (mBc + 0.056) * (mBc + 0.056);
        static constexpr const double mR2_1p = (mBc + 0.492) * (mBc + 0.492);
        static constexpr const char * hqe_prefix = "B_s(*)->D_s(*)";
    };

    /* V -> V Processes */

    struct BstarToDstar {
        using Transition = VToV;
        static constexpr const char * label = "B^*->D^*";
        static constexpr const char * name_Bst = "mass::B_d^*";
        static constexpr const char * name_V   = "mass::D_u^*";
        static constexpr const double mV1 = 5.324;
        static constexpr const double mV2 = 2.010;
        static constexpr const double mBc = 6.2751;
        static constexpr const double mR2_0m = (mBc + 0.000) * (mBc + 0.000);
        static constexpr const double mR2_1m = (mBc + 0.056) * (mBc + 0.056);
        static constexpr const double mR2_1p = (mBc + 0.492) * (mBc + 0.492);
        static constexpr const char * hqe_prefix = "B(*)->D(*)";
    };

    struct BsstarToDsstar {
        using Transition = VToV;
        static constexpr const char * label = "B_s^*->D_s^*";
        static constexpr const char * name_Bst = "mass::B_s^*";
        static constexpr const char * name_P   = "mass::D_s";
        static constexpr const double m_B = 5.324;
        static constexpr const double m_V = 2.010;
        static constexpr const double mBc  = 6.2751;
        static constexpr const double mR2_0m = (mBc + 0.000) * (mBc + 0.000);
        static constexpr const double mR2_1m = (mBc + 0.056) * (mBc + 0.056);
        static constexpr const double mR2_1p = (mBc + 0.492) * (mBc + 0.492);
        static constexpr const char * hqe_prefix = "B_s(*)->D_s(*)";
    };
}

#endif
