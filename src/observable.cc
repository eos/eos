/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <src/observable.hh>
#include <src/rare-b-decays/exclusive-b-to-dilepton.hh>
#include <src/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <src/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil.hh>
#include <src/rare-b-decays/exclusive-b-to-s-gamma.hh>
#include <src/rare-b-decays/form-factors-impl.hh>
#include <src/rare-b-decays/inclusive-b-to-s-dilepton.hh>
#include <src/rare-b-decays/inclusive-b-to-s-gamma.hh>
#include <src/utils/concrete_observable.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <map>

namespace eos
{
    ObservableNameError::ObservableNameError(const std::string & name) :
        Exception("Observable name '" + name + "' is malformed")
    {
    }

    ObservableFactory::ObservableFactory()
    {
    }

    ObservableFactory::~ObservableFactory()
    {
    }

    template <typename Decay_, typename ... Args_>
    std::pair<std::string, ObservableFactory *> make_observable(const char * name,
            double (Decay_::* function)(const Args_ & ...) const)
    {
        std::string sname(name);

        return std::make_pair(sname, make_concrete_observable_factory(sname, function, std::make_tuple()));
    }

    template <typename Decay_, typename Tuple_, typename ... Args_>
    std::pair<std::string, ObservableFactory *> make_observable(const char * name,
            double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        std::string sname(name);

        return std::make_pair(sname, make_concrete_observable_factory(sname, function, kinematics_names));
    }

    ObservablePtr
    Observable::make(const std::string & _name, const Parameters & parameters, const Kinematics & kinematics, const Options & _options)
    {
        static const std::map<std::string, ObservableFactory *> simple_observables
        {
            /* Exclusive Decays */

            // B_q -> ll
            make_observable("B_q->ll::BR",
                    &BToDilepton::branching_ratio),

            // B -> K, cf. [BZ2004v2]
            make_observable("B->K::f_+(s)@BZ2004v2",
                    &BZ2004FormFactors<BToK, PToP>::f_p,
                    std::make_tuple("s")),

            make_observable("B->K::f_T(s)@BZ2004v2",
                    &BZ2004FormFactors<BToK, PToP>::f_t,
                    std::make_tuple("s")),

            make_observable("B->K::f_0(s)@BZ2004v2",
                    &BZ2004FormFactors<BToK, PToP>::f_0,
                    std::make_tuple("s")),

            // B -> K, cf. [BZ2004v3]
            make_observable("B->K::f_+(s)@BZ2004v2Split",
                    &BZ2004FormFactorsSplit<BToK>::f_p,
                    std::make_tuple("s")),

            make_observable("B->K::f_T(s)@BZ2004v2Split",
                    &BZ2004FormFactorsSplit<BToK>::f_t,
                    std::make_tuple("s")),

            make_observable("B->K::f_0(s)@BZ2004v2Split",
                    &BZ2004FormFactorsSplit<BToK>::f_0,
                    std::make_tuple("s")),

            // B -> K, cf. [KMPW2010]
            make_observable("B->K::f_+(s)@KMPW2010",
                    &KMPW2010FormFactors<PToP>::f_p,
                    std::make_tuple("s")),

            make_observable("B->K::f_T(s)@KMPW2010",
                    &KMPW2010FormFactors<PToP>::f_t,
                    std::make_tuple("s")),

            make_observable("B->K::f_0(s)@KMPW2010",
                    &KMPW2010FormFactors<PToP>::f_0,
                    std::make_tuple("s")),

            // B -> K^*, cf. [BZ2004]
            make_observable("B->K^*::V(s)@BZ2004",
                    &BZ2004FormFactors<BToKstar, PToV>::v,
                    std::make_tuple("s")),

            make_observable("B->K^*::A_0(s)@BZ2004",
                    &BZ2004FormFactors<BToKstar, PToV>::a_0,
                    std::make_tuple("s")),

            make_observable("B->K^*::A_1(s)@BZ2004",
                    &BZ2004FormFactors<BToKstar, PToV>::a_1,
                    std::make_tuple("s")),

            make_observable("B->K^*::A_2(s)@BZ2004",
                    &BZ2004FormFactors<BToKstar, PToV>::a_2,
                    std::make_tuple("s")),

            // B -> K^*, cf. [KMPW2010]
            make_observable("B->K^*::V(s)@KMPW2010",
                    &KMPW2010FormFactors<PToV>::v,
                    std::make_tuple("s")),

            make_observable("B->K^*::A_0(s)@KMPW2010",
                    &KMPW2010FormFactors<PToV>::a_0,
                    std::make_tuple("s")),

            make_observable("B->K^*::A_1(s)@KMPW2010",
                    &KMPW2010FormFactors<PToV>::a_1,
                    std::make_tuple("s")),

            make_observable("B->K^*::A_2(s)@KMPW2010",
                    &KMPW2010FormFactors<PToV>::a_2,
                    std::make_tuple("s")),

            // B_s -> phi, cf. [BZ2004]
            make_observable("B_s->phi::V(s)@BZ2004",
                    &BZ2004FormFactors<BToKstar, PToV>::v,
                    std::make_tuple("s")),

            make_observable("B_s->phi::A_0(s)@BZ2004",
                    &BZ2004FormFactors<BToKstar, PToV>::a_0,
                    std::make_tuple("s")),

            make_observable("B_s->phi::A_1(s)@BZ2004",
                    &BZ2004FormFactors<BToKstar, PToV>::a_1,
                    std::make_tuple("s")),

            make_observable("B_s->phi::A_2(s)@BZ2004",
                    &BZ2004FormFactors<BToKstar, PToV>::a_2,
                    std::make_tuple("s")),

            // B -> K^* gamma
            make_observable("B->K^*gamma::BR",
                    &BToKstarGamma::branching_ratio),

            make_observable("B->K^*gamma::S_K^*gamma",
                    &BToKstarGamma::s_kstar_gamma),

            // B -> K ll, Large Recoil
            make_observable("B->Kll::dBR/ds@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->Kll::F_H(s)@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::differential_flat_term,
                    std::make_tuple("s")),

            make_observable("B->Kll::BR@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::F_H@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_flat_term,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::R_K@LargeRecoil",
                    &BToKDilepton<LowRecoil>::integrated_ratio_muons_electrons,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::a_l@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::a_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::c_l@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::c_l,
                    std::make_tuple("s")),

            // B -> K ll, Low Recoil
            make_observable("B->Kll::dBR/ds@LowRecoil",
                    &BToKDilepton<LowRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->Kll::F_H(s)@LowRecoil",
                    &BToKDilepton<LowRecoil>::differential_flat_term,
                    std::make_tuple("s")),

            make_observable("B->Kll::BR@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::F_H@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_flat_term,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::R_K@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_ratio_muons_electrons,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::a_l@LowRecoil",
                    &BToKDilepton<LowRecoil>::a_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::c_l@LowRecoil",
                    &BToKDilepton<LowRecoil>::c_l,
                    std::make_tuple("s")),

            // B -> K^* ll, Large Recoil
            make_observable("B->K^*ll::dBR/ds@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_FB(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^2(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^3(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^4(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^5(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::F_L(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_longitudinal_polarisation,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_FB@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::Abar_FB@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_unnormalized_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::BR@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_L@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_longitudinal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::Fbar_L@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_unnormalized_longitudinal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::s_0^A_FB@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::a_fb_zero_crossing),


            // B -> K^* ll, Low Recoil
            make_observable("B->K^*ll::dBR/ds@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_FB(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^2(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^3(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^4(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::F_L(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_longitudinal_polarisation,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::rho_1(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::rho_1,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::rho_2(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::rho_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_FB@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::Abar_FB@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_unnormalized_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nA_FB@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::BR@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_L@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_longitudinal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nF_L@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_longitudinal_polarisation_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^2@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nA_T^2@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_2_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nA_T^3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_3_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^4@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nA_T^4@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_4_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^1@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_1,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nH_T^1@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_1_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^2@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nH_T^2@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_2_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nH_T^3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_3_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::Re{Y}(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::real_y,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::Im{Y}(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::imag_y,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::Re{C_9^eff}(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::real_c9eff,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::Im{C_9^eff}(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::imag_c9eff,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::a_CP^1(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_1,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::a_CP^2(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::a_CP^3(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::a_CP^mix(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_mix,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::a_CP^1@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_cp_asymmetry_1,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::a_CP^2@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_cp_asymmetry_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::a_CP^3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_cp_asymmetry_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::Gamma+Gammabar@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_cp_summed_decay_width,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::Gamma-Gammabar@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_unnormalized_cp_asymmetry_1,
                    std::make_tuple("s_min", "s_max")),

            /* Inclusive Decays */

            // B->X_s ll, HLMW2005
            make_observable("B->X_sll::dBR/ds@HLMW2005",
                    &BToXsDilepton<HLMW2005>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->X_sll::BR@HLMW2005",
                    &BToXsDilepton<HLMW2005>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            // B->X_s gamma
            make_observable("B->X_sgamma::BR@Minimal",
                    &BToXsGamma<Minimal>::integrated_branching_ratio),
        };

        Options options(_options);
        std::string name(_name);

        std::string::size_type pos;
        while (std::string::npos != (pos = name.rfind(',')))
        {
            std::string::size_type sep(name.find('=', pos + 1));
            if (std::string::npos == sep)
                throw ObservableNameError(_name);

            std::string key(name.substr(pos + 1, sep - pos - 1));
            std::string value(name.substr(sep + 1));

            options.set(key, value);
            name.erase(pos);
        }

        auto i(simple_observables.find(name));
        if (simple_observables.end() == i)
            return ObservablePtr();

        return i->second->make(parameters, kinematics, options);
    }
}
