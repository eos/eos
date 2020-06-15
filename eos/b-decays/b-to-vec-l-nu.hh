/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018, 2019 Ahmet Kokulu
 * Copyright (c) 2019 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_Dstar_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_Dstar_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /*
     * Decay: B -> D^* l nu
     */
    class BToVectorLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToVectorLeptonNeutrino>
    {
        public:
            BToVectorLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToVectorLeptonNeutrino();

            // Differential Observables
            double differential_branching_ratio(const double & q2) const;
            double differential_a_fb_leptonic(const double & q2) const;
            double differential_J1c_normalized(const double & q2) const;
            double differential_J1s_normalized(const double & q2) const;
            double differential_J2c_normalized(const double & q2) const;
            double differential_J2s_normalized(const double & q2) const;
            double differential_J3_normalized(const double & q2) const;
            double differential_J4_normalized(const double & q2) const;
            double differential_J5_normalized(const double & q2) const;
            double differential_J6c_normalized(const double & q2) const;
            double differential_J6s_normalized(const double & q2) const;
            double differential_J7_normalized(const double & q2) const;
            double differential_J8_normalized(const double & q2) const;
            double differential_J9_normalized(const double & q2) const;
            double differential_ratio_tau_mu(const double & q2) const;

            // Differential Observables - normalized(|Vcb|=1)
            double normalized_differential_branching_ratio(const double & q2) const;

            // Four Differential Observables - normalized(|Vcb|=1)
            double normalized_four_differential_decay_width(const double & q2, const double & c_theta_l, const double & c_theta_d, const double & phi) const;

            // Integrated Observables
            double integrated_branching_ratio(const double & q2_min, const double & q2_max) const;
            double integrated_a_fb_leptonic(const double & q2_min, const double & q2_max) const;

            double integrated_amplitude_polarization_L(const double & q2_min, const double & q2_max) const;
            double integrated_amplitude_polarization_T(const double & q2_min, const double & q2_max) const;
            double integrated_f_L(const double & q2_min, const double & q2_max) const;
            double integrated_ftilde_L(const double & q2_min, const double & q2_max) const;
            double integrated_a_c_1(const double & q2_min, const double & q2_max) const;
            double integrated_a_c_2(const double & q2_min, const double & q2_max) const;
            double integrated_a_c_3(const double & q2_min, const double & q2_max) const;
            double integrated_a_t_1(const double & q2_min, const double & q2_max) const;
            double integrated_a_t_2(const double & q2_min, const double & q2_max) const;
            double integrated_a_t_3(const double & q2_min, const double & q2_max) const;
            double integrated_J1c_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J1s_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J2c_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J2s_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J3_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J4_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J5_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J6c_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J6s_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J7_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J8_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_J9_normalized(const double & q2_min, const double & q2_max) const;

            // Integrated Observables - normalized(|Vcb|=1)
            double normalized_integrated_branching_ratio(const double & q2_min, const double & q2_max) const;

            // PDF
            double differential_pdf_q2(const double & q2) const;
            double differential_pdf_w(const double & w) const;
            double integrated_pdf_q2(const double & q2_min, const double & q2_max) const;
            double integrated_pdf_w(const double & w_min, const double & w_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_q2;
            static const std::string kinematics_description_c_theta_l;
            static const std::string kinematics_description_c_theta_d;
            static const std::string kinematics_description_phi;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;
    };
}

#endif
