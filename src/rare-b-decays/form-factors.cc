/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
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

#include <src/rare-b-decays/form-factors-impl.hh>
#include <src/utils/destringify.hh>

#include <map>
#include <cmath>

namespace eos
{
    /* P -> V Processes */

    /* B_{u,d} -> K^* */
    template class BZ2004FormFactors<BToKstar, PToV>;

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

    /* B_s -> phi */
    template class BZ2004FormFactors<BsToPhi, PToV>;

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
            { KeyType("B->K^*",     "BZ2004"), &BZ2004FormFactors<BToKstar, PToV>::make   },
            { KeyType("Bs->phi",    "BZ2004"), &BZ2004FormFactors<BsToPhi, PToV>::make    },
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
    template class BZ2004FormFactors<BToK, PToP>;

    /* For the values below, cf. [BZ2004v2], Table 1, p. 8 */
    template <> const double BZ2004FormFactors<BToK, PToP>::_r1_p     = 0.162;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r2_p     = 0.173;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r1_t     = 0.161;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r2_t     = 0.198;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r2_0     = 0.330;
    template <> const double BZ2004FormFactors<BToK, PToP>::_mfit2    = 37.46;
    template <> const double BZ2004FormFactors<BToK, PToP>::_m12      = 5.41 * 5.41;


    /* B_{u,d} -> K */
    template class KMPW2010FormFactors<BToK>;

    /* For the values below, cf. [KMPW2010], Table 4, p. 31 */
    template <> const double KMPW2010FormFactors<BToK>::_f0_p     =  0.34;
    template <> const double KMPW2010FormFactors<BToK>::_f0_0     =  0.34;
    template <> const double KMPW2010FormFactors<BToK>::_f0_t     =  0.39;
    template <> const double KMPW2010FormFactors<BToK>::_b1_p     = -2.1;
    template <> const double KMPW2010FormFactors<BToK>::_b1_0     = -4.3;
    template <> const double KMPW2010FormFactors<BToK>::_b1_t     = -2.2;

    template <> const double KMPW2010FormFactors<BToK>::_tau_p    = (_m_B + _m_K) * (_m_B + _m_K);
    template <> const double KMPW2010FormFactors<BToK>::_tau_m    = (_m_B - _m_K) * (_m_B - _m_K);
    template <> const double KMPW2010FormFactors<BToK>::_tau_0    = _tau_p - std::sqrt(_tau_p * _tau_p - _tau_m * _tau_p);

    // Masses, cf. PDG 2008
    template <> const double KMPW2010FormFactors<BToK>::_m_B      = 5.280;
    template <> const double KMPW2010FormFactors<BToK>::_m_K      = 0.498;
    template <> const double KMPW2010FormFactors<BToK>::_m_Bs2    = 5.325 * 5.325;

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
            { KeyType("B->K",     "BZ2004v2"), &BZ2004FormFactors<BToK, PToP>::make   },
            { KeyType("B->K",     "KMPW2010"), &KMPW2010FormFactors<BToK>::make   }
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
