/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013, 2014, 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_D_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_D_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*
     * Decay: B -> D l nu
     */
    class BToDLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToDLeptonNeutrino>
    {
        public:
            BToDLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToDLeptonNeutrino();

            // Differential Observables
            double differential_branching_ratio(const double & s) const;

            // Differential Observables - normalized(|Vcb|=1)
            double normalized_differential_branching_ratio(const double & s) const;

            // Integrated Observables
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;

            // Integrated Observables - normalized(|Vcb|=1)
            double normalized_integrated_branching_ratio(const double & s_min, const double & s_max) const;

            // R_D
            double differential_r_d(const double & s) const;
            double integrated_r_d() const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_s;
    };
}

#endif
