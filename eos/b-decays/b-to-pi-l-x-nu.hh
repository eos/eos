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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_PI_L_X_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_PI_L_X_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /*
     * Decay: B -> pi l X_nu
     */
    class BToPiLeptonInclusiveNeutrinos :
        public ParameterUser,
        public PrivateImplementationPattern<BToPiLeptonInclusiveNeutrinos>
    {
        public:
            BToPiLeptonInclusiveNeutrinos(const Parameters & parameters, const Options & options);
            ~BToPiLeptonInclusiveNeutrinos();

            double differential_decay_width_1nu(const double & s, const double & c_theta_tau) const;

            double differential_decay_width_3nu(const double & s, const double & snunubar,
                    const double & c_theta_tau, const double & phi, const double & c_theta_mu_star) const;

            double integrated_decay_width_1nu(const double & s_min, const double & s_max) const;

            double integrated_decay_width_3nu(const double & s_min, const double & s_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_s;
            static const std::string kinematics_description_snunubar;
            static const std::string kinematics_description_c_theta;
            static const std::string kinematics_description_c_theta_tau;
            static const std::string kinematics_description_c_theta_mu_star;
            static const std::string kinematics_description_phi;

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
