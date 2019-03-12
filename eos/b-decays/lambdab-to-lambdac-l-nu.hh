/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Ahmet Kokulu
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

#ifndef EOS_GUARD_EOS_B_DECAYS_BARYONIC_B_TO_C_LEPTON_NEUTRINO_HH
#define EOS_GUARD_EOS_B_DECAYS_BARYONIC_B_TO_C_LEPTON_NEUTRINO_HH 1

#include <eos/decays.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*
     * Decay: Lambda_b -> Lambda_c lnu
     */
    class LambdaBToLambdaCLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<LambdaBToLambdaCLeptonNeutrino>
    {
        public:
            LambdaBToLambdaCLeptonNeutrino(const Parameters &, const Options &);
            ~LambdaBToLambdaCLeptonNeutrino();

            double differential_branching_ratio(const double & s) const;
            double differential_a_fb_leptonic(const double & s) const;
            double differential_a_fb_hadronic(const double & s) const;
            double differential_a_fb_combined(const double & s) const;
            double differential_fzero(const double & s) const;
            double differential_ratio_tau_mu(const double & s) const;
            double differential_ratio_a_fb_hadronic_tau_mu(const double & s) const;

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
            double integrated_ratio_tau_mu(const double & s_min_mu, const double & s_min_tau, const double & s_max_mu, const double & s_max_tau) const;
            double integrated_ratio_a_fb_hadronic_tau_mu(const double & s_min_mu, const double & s_min_tau, const double & s_max_mu, const double & s_max_tau) const;

    /*!
    * Descriptions of the process and its kinematics.
    */
        static const std::string description;
        static const std::string kinematics_description_s;
    };
}

#endif
