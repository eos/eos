/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Elena Graverini
 * Copyright (c) 2017 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_B_DECAYS_LAMBDAB_TO_LAMBDAC2625_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_LAMBDAB_TO_LAMBDAC2625_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*
     * Decay: Lambda_b -> Lambda_c(2625) l nu
     */
    class LambdaBToLambdaC2625LeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<LambdaBToLambdaC2625LeptonNeutrino>
    {
        public:
            LambdaBToLambdaC2625LeptonNeutrino(const Parameters & parameters, const Options & options);
            ~LambdaBToLambdaC2625LeptonNeutrino();

            // [BBGIOvD] parametrization for the differential decay width
            double a_l(const double & s) const;
            double b_l(const double & s) const;
            double c_l(const double & s) const;
            double gamma_0(const double & s) const;

            // Differential Observables
            double differential_branching_ratio(const double & s) const;
            double differential_forward_backward_asymmetry(const double & s) const;
            double double_differential_branching_ratio(const double & s, const double & theta_l) const;

            // Integrated Observables
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const;

            // Integrated Observables (normalized to 1)
            double normalized_integrated_branching_ratio(const double & s_min, const double & s_max) const;

            // R_Lambda_c_2625
            double differential_r_lambdac2625(const double & s) const;
            double integrated_r_lambdac2625() const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_s;
            static const std::string kinematics_description_c_theta_l;
    };
}

#endif
