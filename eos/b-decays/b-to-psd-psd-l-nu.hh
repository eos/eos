/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_PSD_PSD_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_PSD_PSD_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /*
     * Decay: B -> P P' l nu
     */
    class BToPPLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToPPLeptonNeutrino>
    {
        public:
            BToPPLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToPPLeptonNeutrino();

            /*!
             * Observables
             */
            // Differential in k2 and q2
            double double_differential_decay_width(const double & q2, const double & k2) const;
            double double_differential_mesonic_afb(const double & q2, const double & k2) const;
            double double_differential_branching_ratio(const double & q2, const double & k2) const;

            // Double differential observables integrated in a q2-sqrts bin
            double integrated_branching_ratio(const double & q2_min, const double & q2_max, const double & sqrt_k2_min, const double & sqrt_k2_max) const;
            double integrated_mesonic_afb(const double & q2_min, const double & q2_max, const double & sqrt_k2_min, const double & sqrt_k2_max) const;

            // Fully integrated branching ratio
            double fully_integrated_branching_ratio() const;
            double fully_integrated_mesonic_afb() const;

            // Branching ratio in 1D windows
            double q2_integrated_branching_ratio(const double & sqrt_k2_min, const double & sqrt_k2_max) const;
            double q2_integrated_mesonic_afb(const double & sqrt_k2_min, const double & sqrt_k2_max) const;
            double s_integrated_branching_ratio(const double & q2_min, const double & q2_max) const;
            double s_integrated_mesonic_afb(const double & q2_min, const double & q2_max) const;

            // Single differential branching ratios
            double integrated_branching_ratio_q2(const double & q2) const;
            double integrated_mesonic_afb_q2(const double & q2) const;
            double integrated_branching_ratio_sqrt_k2(const double & sqrt_k2) const;
            double integrated_mesonic_afb_sqrt_k2(const double & sqrt_k2) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_q2;
            static const std::string kinematics_description_k2;
            static const std::string kinematics_description_c_theta_l;
            static const std::string kinematics_description_c_theta_nu;
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
