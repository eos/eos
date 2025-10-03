/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Matthew Kirk
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

#ifndef EOS_GUARD_EOS_S_DECAYS_K_TO_PI_L_NU_HH
#define EOS_GUARD_EOS_S_DECAYS_K_TO_PI_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /*
     * Decay: K -> pi l nu
     */
    class KToPiLeptonNeutrino : public ParameterUser, public PrivateImplementationPattern<KToPiLeptonNeutrino>
    {
        public:
            KToPiLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~KToPiLeptonNeutrino();

            // Two-fold differential observables - normalized to |V_us|=1
            double normalized_two_differential_decay_width(const double & q2, const double & c_theta_l) const;

            // Single-differential Observables
            double differential_branching_ratio(const double & q2) const;
            double differential_a_fb_leptonic(const double & q2) const;
            double differential_flat_term(const double & q2) const;
            double differential_lepton_polarization(const double & q2) const;

            // Single-differential Observables - normalized to |V_us|=1
            double normalized_differential_branching_ratio(const double & q2) const;

            // Integrated Observables
            double total_branching_ratio() const;
            double integrated_branching_ratio(const double & q2_min, const double & q2_max) const;
            double integrated_a_fb_leptonic(const double & q2_min, const double & q2_max) const;
            double integrated_flat_term(const double & q2_min, const double & q2_max) const;
            double integrated_lepton_polarization(const double & q2_min, const double & q2_max) const;

            // Integrated Observables - normalized to |V_us|=1
            double normalized_integrated_branching_ratio(const double & q2_min, const double & q2_max) const;
            double normalized_integrated_decay_width(const double & q2_min, const double & q2_max) const;
            double normalized_integrated_decay_width_0(const double & q2_min, const double & q2_max) const;
            double normalized_integrated_decay_width_p(const double & q2_min, const double & q2_max) const;

            // PDF
            double differential_pdf_q2(const double & q2) const;
            double integrated_pdf_q2(const double & q2_min, const double & q2_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_q2;
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
    };
} // namespace eos

#endif
