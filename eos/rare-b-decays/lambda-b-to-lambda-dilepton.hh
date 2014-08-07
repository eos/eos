/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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

namespace eos
{
    /*
     * Decay: Lambda_b -> Lambda l^+ l^- at large recoil, cf. [BFvD2014]
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
    };

    /*
     * Decay: Lambda_b -> Lambda l^+ l^- at low recoil, cf. [BFvD2014]
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
    };
}

#endif
