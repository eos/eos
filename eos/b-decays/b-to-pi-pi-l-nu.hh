/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_B_DECAYS_EXCLUSIVE_B_TO_PI_PI_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_EXCLUSIVE_B_TO_PI_PI_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*
     * Decay: B -> pi pi l nu
     */
    class BToPiPiLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToPiPiLeptonNeutrino>
    {
        public:
            BToPiPiLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToPiPiLeptonNeutrino();

            // Differential observables
            double double_differential_branching_ratio(const double & q2, const double & k2) const;
            double double_differential_forward_backward_asymmetry(const double & q2, const double & k2) const;
            double triple_differential_branching_ratio(const double & q2, const double & k2, const double & z) const;
            double partial_waves(const double & q2, const double & k2, const double & z) const;

            // Integrated observables
            double integrated_branching_ratio(const double & q2min, const double & q2max,
                    const double & k2min, const double & k2max,
                    const double & zmin, const double & zmax) const;
            double integrated_forward_backward_asymmetry(const double & q2min, const double & q2max,
                    const double & k2min, const double & k2max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_q2;
            static const std::string kinematics_description_k2;
            static const std::string kinematics_description_z;
    };
}

#endif
