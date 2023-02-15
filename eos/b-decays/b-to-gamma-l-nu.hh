/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Philip Lüghausen
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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_GAMMA_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_GAMMA_L_NU_HH 1

#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    class BToGammaLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToGammaLeptonNeutrino>
    {
        public:
            BToGammaLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToGammaLeptonNeutrino();

            // Observables
            double forward_backward_asymmetry(const double & E_gamma_min) const;
            double fully_differential_decay_width(const double & E_gamma, const double & costheta) const;
            double integrated_branching_ratio(const double & E_gamma_min) const;
            double differential_decay_width_dEgamma(const double & E_gamma) const;
            double integrated_decay_width(const double & E_gamma_min) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_Egamma;
            static const std::string kinematics_description_c_theta_l;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();

            /* Diagnostics for unit tests */
            Diagnostics diagnostics() const;
    };
}

#endif

