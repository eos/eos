/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014, 2015 Danny van Dyk
 * Copyright (c) 2017 Thomas Blake
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_BARYONIC_B_TO_S_DILEPTON_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_BARYONIC_B_TO_S_DILEPTON_HH 1

#include <eos/rare-b-decays/decays.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /*
     * Decay: Lambda_b -> Lambda l^+ l^- at large recoil, cf. [BFvD:2014]
     */
    template <> class LambdaBToLambdaDilepton<LargeRecoil> :
        public ParameterUser,
        public PrivateImplementationPattern<LambdaBToLambdaDilepton<LargeRecoil>>
    {
        public:
            LambdaBToLambdaDilepton(const Parameters &, const Options &);
            ~LambdaBToLambdaDilepton();

            double differential_branching_ratio(const double & s) const;
            double differential_a_fb_leptonic(const double & s) const;
            double differential_a_fb_hadronic(const double & s) const;
            double differential_a_fb_combined(const double & s) const;
            double differential_fzero(const double & s) const;

            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_a_fb_leptonic(const double & s_min, const double & s_max) const;
            double integrated_a_fb_hadronic(const double & s_min, const double & s_max) const;
            double integrated_a_fb_combined(const double & s_min, const double & s_max) const;
            double integrated_fzero(const double & s_min, const double & s_max) const;

            double integrated_m1(const double & s_min, const double & s_max) const;
            double integrated_m2(const double & s_min, const double & s_max) const;
            double integrated_m3(const double & s_min, const double & s_max) const;
            double integrated_m4(const double & s_min, const double & s_max) const;
            double integrated_m5(const double & s_min, const double & s_max) const;
            double integrated_m6(const double & s_min, const double & s_max) const;
            double integrated_m7(const double & s_min, const double & s_max) const;
            double integrated_m8(const double & s_min, const double & s_max) const;
            double integrated_m9(const double & s_min, const double & s_max) const;
            double integrated_m10(const double & s_min, const double & s_max) const;
            double integrated_m11(const double & s_min, const double & s_max) const;
            double integrated_m12(const double & s_min, const double & s_max) const;
            double integrated_m13(const double & s_min, const double & s_max) const;
            double integrated_m14(const double & s_min, const double & s_max) const;
            double integrated_m15(const double & s_min, const double & s_max) const;
            double integrated_m16(const double & s_min, const double & s_max) const;
            double integrated_m17(const double & s_min, const double & s_max) const;
            double integrated_m18(const double & s_min, const double & s_max) const;
            double integrated_m19(const double & s_min, const double & s_max) const;
            double integrated_m20(const double & s_min, const double & s_max) const;
            double integrated_m21(const double & s_min, const double & s_max) const;
            double integrated_m22(const double & s_min, const double & s_max) const;
            double integrated_m23(const double & s_min, const double & s_max) const;
            double integrated_m24(const double & s_min, const double & s_max) const;
            double integrated_m25(const double & s_min, const double & s_max) const;
            double integrated_m26(const double & s_min, const double & s_max) const;
            double integrated_m27(const double & s_min, const double & s_max) const;
            double integrated_m28(const double & s_min, const double & s_max) const;
            double integrated_m29(const double & s_min, const double & s_max) const;
            double integrated_m30(const double & s_min, const double & s_max) const;
            double integrated_m31(const double & s_min, const double & s_max) const;
            double integrated_m32(const double & s_min, const double & s_max) const;
            double integrated_m33(const double & s_min, const double & s_max) const;
            double integrated_m34(const double & s_min, const double & s_max) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };

    /*
     * Decay: Lambda_b -> Lambda l^+ l^- at low recoil, cf. [BFvD:2014]
     */
    template <> class LambdaBToLambdaDilepton<LowRecoil> :
        public ParameterUser,
        public PrivateImplementationPattern<LambdaBToLambdaDilepton<LowRecoil>>
    {
        public:
            LambdaBToLambdaDilepton(const Parameters &, const Options &);
            ~LambdaBToLambdaDilepton();

            double differential_branching_ratio(const double & s) const;
            double differential_a_fb_leptonic(const double & s) const;
            double differential_a_fb_hadronic(const double & s) const;
            double differential_a_fb_combined(const double & s) const;
            double differential_fzero(const double & s) const;

            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_a_fb_leptonic(const double & s_min, const double & s_max) const;
            double integrated_a_fb_hadronic(const double & s_min, const double & s_max) const;
            double integrated_a_fb_combined(const double & s_min, const double & s_max) const;
            double integrated_fzero(const double & s_min, const double & s_max) const;

            double integrated_k1ss(const double & s_min, const double & s_max) const;
            double integrated_k1cc(const double & s_min, const double & s_max) const;
            double integrated_k1c(const double & s_min, const double & s_max) const;
            double integrated_k2ss(const double & s_min, const double & s_max) const;
            double integrated_k2cc(const double & s_min, const double & s_max) const;
            double integrated_k2c(const double & s_min, const double & s_max) const;
            double integrated_k3sc(const double & s_min, const double & s_max) const;
            double integrated_k3s(const double & s_min, const double & s_max) const;
            double integrated_k4sc(const double & s_min, const double & s_max) const;
            double integrated_k4s(const double & s_min, const double & s_max) const;

            double integrated_m1(const double & s_min, const double & s_max) const;
            double integrated_m2(const double & s_min, const double & s_max) const;
            double integrated_m3(const double & s_min, const double & s_max) const;
            double integrated_m4(const double & s_min, const double & s_max) const;
            double integrated_m5(const double & s_min, const double & s_max) const;
            double integrated_m6(const double & s_min, const double & s_max) const;
            double integrated_m7(const double & s_min, const double & s_max) const;
            double integrated_m8(const double & s_min, const double & s_max) const;
            double integrated_m9(const double & s_min, const double & s_max) const;
            double integrated_m10(const double & s_min, const double & s_max) const;
            double integrated_m11(const double & s_min, const double & s_max) const;
            double integrated_m12(const double & s_min, const double & s_max) const;
            double integrated_m13(const double & s_min, const double & s_max) const;
            double integrated_m14(const double & s_min, const double & s_max) const;
            double integrated_m15(const double & s_min, const double & s_max) const;
            double integrated_m16(const double & s_min, const double & s_max) const;
            double integrated_m17(const double & s_min, const double & s_max) const;
            double integrated_m18(const double & s_min, const double & s_max) const;
            double integrated_m19(const double & s_min, const double & s_max) const;
            double integrated_m20(const double & s_min, const double & s_max) const;
            double integrated_m21(const double & s_min, const double & s_max) const;
            double integrated_m22(const double & s_min, const double & s_max) const;
            double integrated_m23(const double & s_min, const double & s_max) const;
            double integrated_m24(const double & s_min, const double & s_max) const;
            double integrated_m25(const double & s_min, const double & s_max) const;
            double integrated_m26(const double & s_min, const double & s_max) const;
            double integrated_m27(const double & s_min, const double & s_max) const;
            double integrated_m28(const double & s_min, const double & s_max) const;
            double integrated_m29(const double & s_min, const double & s_max) const;
            double integrated_m30(const double & s_min, const double & s_max) const;
            double integrated_m31(const double & s_min, const double & s_max) const;
            double integrated_m32(const double & s_min, const double & s_max) const;
            double integrated_m33(const double & s_min, const double & s_max) const;
            double integrated_m34(const double & s_min, const double & s_max) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };
}

#endif
