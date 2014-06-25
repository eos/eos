/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013, 2014 Danny van Dyk
 * Copyright (c) 2010, 2011 Christian Wacker
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

#include <eos/form-factors/analytic-b-to-pi.hh>
#include <eos/form-factors/mesonic-impl.hh>
#include <eos/utils/destringify.hh>

#include <map>
#include <cmath>

namespace eos
{
    /* P -> V Processes */

    /* B_{u,d} -> K^* */

    /* For the values below, cf. [BZ2004], Table 8, p. 28 */
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_v_r1 = +0.923;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_v_r2 = -0.511;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_v_m2r = 5.32 * 5.32;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_v_m2fit = 49.40;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a0_r1 = +1.364;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a0_r2 = -0.990;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a0_m2r = 5.28 * 5.28;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a0_m2fit = 36.78;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a1_r2 = +0.290;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a1_m2fit = 40.38;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a2_r1 = -0.084;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a2_r2 = +0.342;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a2_m2fit = 52.00;
    template class BZ2004FormFactors<BToKstar, PToV>;


    // mass B_d, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToV>::_m_B       = 5.2795;
    // mass K^*0, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToV>::_m_Kstar   = 0.89594;
    // mass B_s (0-), cf. [KMPW2010]
    const double KMPW2010FormFactors<PToV>::_m_Bs2_0m  = 5.366 * 5.366;
    // mass B_s (1-), cf. [KMPW2010]
    const double KMPW2010FormFactors<PToV>::_m_Bs2_1m  = 5.412 * 5.412;
    // mass B_s (1+), cf. [KMPW2010]
    const double KMPW2010FormFactors<PToV>::_m_Bs2_1p  = 5.829 * 5.829;
    const double KMPW2010FormFactors<PToV>::_tau_p     = (_m_B + _m_Kstar) * (_m_B + _m_Kstar);
    const double KMPW2010FormFactors<PToV>::_tau_m     = (_m_B - _m_Kstar) * (_m_B - _m_Kstar);
    const double KMPW2010FormFactors<PToV>::_tau_0     = _tau_p - std::sqrt(_tau_p * _tau_p - _tau_m * _tau_p);
    template class KMPW2010FormFactors<PToV>;


    /* B_s -> phi */
    /* For the values below, cf. [BZ2004], Table 8, p. 28 */
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_v_r1 = +1.484;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_v_r2 = -1.049;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_v_m2r = 5.42 * 5.42;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_v_m2fit = 39.52;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a0_r1 = +3.310;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a0_r2 = -2.835;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a0_m2r = 5.37 * 5.37;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a0_m2fit = 31.57;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a1_r2 = +0.308;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a1_m2fit = 36.54;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a2_r1 = -0.054;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a2_r2 = +0.288;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a2_m2fit = 48.94;
    template class BZ2004FormFactors<BsToPhi, PToV>;


    FormFactors<PToV>::~FormFactors()
    {
    }

    std::shared_ptr<FormFactors<PToV>>
    FormFactorFactory<PToV>::create(const std::string & label, const Parameters & parameters)
    {
        std::shared_ptr<FormFactors<PToV>> result;

        typedef std::tuple<std::string, std::string> KeyType;
        typedef std::function<FormFactors<PToV> * (const Parameters &, unsigned)> ValueType;
        static const std::map<KeyType, ValueType> form_factors
        {
            { KeyType("B->K^*",     "BZ2004"),   &BZ2004FormFactors<BToKstar, PToV>::make   },
            { KeyType("B->K^*",     "KMPW2010"), &KMPW2010FormFactors<PToV>::make           },
            { KeyType("Bs->phi",    "BZ2004"),   &BZ2004FormFactors<BsToPhi, PToV>::make    },
        };

        /*
         * Labels have the form
         *
         *   PROCESS@NAME[:SET]
         *
         * The brackets indicate the latter part to be optional.
         */

        std::string process, name, input(label);
        unsigned set(0);

        std::string::size_type sep_at(input.find('@')), sep_colon(input.find(':'));
        if (std::string::npos == sep_at)
            return result;

        if (std::string::npos != sep_colon)
        {
            set = destringify<unsigned>(input.substr(sep_colon + 1));
            input.erase(sep_colon + 1);
        }

        name = input.substr(sep_at + 1);
        process = input.substr(0, sep_at);

        auto i = form_factors.find(KeyType(process, name));
        if (form_factors.cend() == i)
            return result;

        result = std::shared_ptr<FormFactors<PToV>>(i->second(parameters, set));

        return result;
    }

    /* P -> P Processes */

    /* B_{u,d} -> K */

    // [BCL2008]
    template class BCL2008FormFactors<BToK>;

    /* For the values below, cf. [BZ2004v2], Table 1, p. 8 */
    template <> const double BZ2004FormFactors<BToK, PToP>::_r1_p     = 0.162;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r2_p     = 0.173;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r1_t     = 0.161;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r2_t     = 0.198;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r2_0     = 0.330;
    template <> const double BZ2004FormFactors<BToK, PToP>::_mfit2    = 37.46;
    template <> const double BZ2004FormFactors<BToK, PToP>::_m12      = 5.41 * 5.41;
    template class BZ2004FormFactors<BToK, PToP>;


    // [BZ2004v3], Table B, p. 26
    template <> const double BZ2004FormFactorsSplit<BToK>::_r1_p_asymptotic    = 0.0541;
    template <> const double BZ2004FormFactorsSplit<BToK>::_r2_p_asymptotic    = 0.2166;
    template <> const double BZ2004FormFactorsSplit<BToK>::_r2_0_asymptotic    = 0.2719;
    template <> const double BZ2004FormFactorsSplit<BToK>::_r1_t_asymptotic    = 0.0244;
    template <> const double BZ2004FormFactorsSplit<BToK>::_r2_t_asymptotic    = 0.2590;
    template <> const double BZ2004FormFactorsSplit<BToK>::_mfit2_0_asymptotic = 30.33;

    // [BZ2004v3], Table D, p. 28
    template <> const double BZ2004FormFactorsSplit<BToK>::_m12_asymptotic     = 5.41 * 5.41;

    // [BZ2004v3], Table C, set 2
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_a_1 =  0.310;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_b_1 =  0.930e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_c_1 =  0.139e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_d_1 = -0.083e-3;

    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_a_2 =  0.228;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_b_2 = -0.632e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_c_2 =  0.017e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_d_2 = -0.143e-3;

    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_a_4 = -0.173;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_b_4 = -0.947e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_c_4 =  0.005e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_p_d_4 =  0.196e-3;

    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_a_1 =  0.308;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_b_1 =  0.106e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_c_1 =  0.026e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_d_1 = -0.048e-3;

    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_a_2 =  0.226;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_b_2 = -1.031e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_c_2 = -0.092e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_d_2 = -0.005e-3;

    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_a_4 = -0.170;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_b_4 = -0.838e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_c_4 =  0.209e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_0_d_4 =  0.001e-3;

    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_a_1 =  0.381;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_b_1 =  1.056e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_c_1 =  0.167e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_d_1 = -0.108e-3;

    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_a_2 =  0.264;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_b_2 = -0.858e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_c_2 = -0.011e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_d_2 = -0.153e-3;

    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_a_4 = -0.217;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_b_4 = -1.165e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_c_4 =  0.101e-2;
    template <> const double BZ2004FormFactorsSplit<BToK>::_f_t_d_4 =  0.187e-3;

    template class BZ2004FormFactorsSplit<BToK>;


    // mass B_u, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToP>::_m_B      = 5.27917;
    // mass K^+, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToP>::_m_K      = 0.493677;
    // mass B_s^* (1-), cf. [KMPW2010]
    const double KMPW2010FormFactors<PToP>::_m_Bs2    = 5.412 * 5.412;
    const double KMPW2010FormFactors<PToP>::_tau_p    = (_m_B + _m_K) * (_m_B + _m_K);
    const double KMPW2010FormFactors<PToP>::_tau_m    = (_m_B - _m_K) * (_m_B - _m_K);
    const double KMPW2010FormFactors<PToP>::_tau_0    = _tau_p - std::sqrt(_tau_p * _tau_p - _tau_m * _tau_p);
    template class KMPW2010FormFactors<PToP>;


    // mass B_u, cf. [PDG 2010]
    const double BFW2010FormFactors<PToP>::_m_B       = 5.27917;
    // mass K^+, cf. [PDG 2010]
    const double BFW2010FormFactors<PToP>::_m_K       = 0.493677;
    // mass B_s^* (1-)
    const double BFW2010FormFactors<PToP>::_m_Bs2     = 5.412 * 5.412;
    const double BFW2010FormFactors<PToP>::_tau_p     = (_m_B + _m_K) * (_m_B + _m_K);
    const double BFW2010FormFactors<PToP>::_tau_m     = (_m_B - _m_K) * (_m_B - _m_K);
    const double BFW2010FormFactors<PToP>::_tau_0     = _tau_p - std::sqrt(_tau_p * _tau_p - _tau_m * _tau_p);
    template class BFW2010FormFactors<PToP>;

    /* B_{u,d} -> pi */

    // [BCL2008]
    template class BCL2008FormFactors<BToPi>;


    FormFactors<PToP>::~FormFactors()
    {
    }

    std::shared_ptr<FormFactors<PToP>>
    FormFactorFactory<PToP>::create(const std::string & label, const Parameters & parameters)
    {
        std::shared_ptr<FormFactors<PToP>> result;

        typedef std::tuple<std::string, std::string> KeyType;
        typedef std::function<FormFactors<PToP> * (const Parameters &, unsigned)> ValueType;
        static const std::map<KeyType, ValueType> form_factors
        {
            // parametrizations
            { KeyType("B->K",     "BZ2004v2"),      &BZ2004FormFactors<BToK, PToP>::make     },
            { KeyType("B->K",     "BZ2004v2Split"), &BZ2004FormFactorsSplit<BToK>::make      },
            { KeyType("B->K",     "KMPW2010"),      &KMPW2010FormFactors<PToP>::make         },
            { KeyType("B->K",     "BFW2010"),       &BFW2010FormFactors<PToP>::make          },
            { KeyType("B->pi",    "BCL2008"),       &BCL2008FormFactors<BToPi>::make         },
            // analytic computations
            { KeyType("B->pi",    "DKMMO2008"), &AnalyticFormFactorBToPiDKMMO2008::make      },
        };

        /*
         * Labels have the form
         *
         *   PROCESS@NAME[:SET]
         *
         * The brackets indicate the latter part to be optional.
         */

        std::string process, name, input(label);
        unsigned set(0);

        std::string::size_type sep_at(input.find('@')), sep_colon(input.find(':'));
        if (std::string::npos == sep_at)
            return result;

        if (std::string::npos != sep_colon)
        {
            set = destringify<unsigned>(input.substr(sep_colon + 1));
            input.erase(sep_colon + 1);
        }

        name = input.substr(sep_at + 1);
        process = input.substr(0, sep_at);

        auto i = form_factors.find(KeyType(process, name));
        if (form_factors.cend() == i)
            return result;

        result = std::shared_ptr<FormFactors<PToP>>(i->second(parameters, set));

        return result;
    }
}
