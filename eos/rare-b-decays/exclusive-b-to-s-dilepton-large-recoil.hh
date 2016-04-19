/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_LARGE_RECOIL_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_LARGE_RECOIL_HH 1

#include <eos/decays.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*
     * Decay: B -> K^* l lbar at Large Recoil, cf. [BHP2008]
     */
    template <>
    class BToKstarDilepton<LargeRecoil> :
        public ParameterUser,
        public PrivateImplementationPattern<BToKstarDilepton<LargeRecoil>>
    {
        public:
            BToKstarDilepton(const Parameters & parameters, const Options & options);
            ~BToKstarDilepton();

            // Soft form factors
            double xi_perp(const double & s) const;
            double xi_para(const double & s) const;

            // [BHP2008], Appendix C
            complex<double> a_long(const Helicity & h, const double & s) const;
            complex<double> a_perp(const Helicity & h, const double & s) const;
            complex<double> a_par(const Helicity & h, const double & s) const;
            complex<double> a_timelike(const double & s) const;
            complex<double> a_scalar(const double & s) const;
            complex<double> a_par_perp(const double & s) const;
            complex<double> a_t_long(const double & s) const;
            complex<double> a_t_perp(const double & s) const;
            complex<double> a_t_par(const double & s) const;
            complex<double> a_long_par(const double & s) const;
            complex<double> a_long_perp(const double & s) const;

            // Inverse Observables
            double a_fb_zero_crossing() const;

            // Four Differential Observables
            double four_differential_decay_width(const double & s, const double & c_theta_l, const double & c_theta_k, const double & phi) const;

            // Single Differential Observables
            double differential_branching_ratio(const double & s) const;
            double differential_decay_width(const double & s) const;
            double differential_isospin_asymmetry(const double & s) const;
            double differential_forward_backward_asymmetry(const double & s) const;
            double differential_longitudinal_polarisation(const double & s) const;
            double differential_transversal_polarisation(const double & s) const;
            double differential_transverse_asymmetry_2(const double & s) const;
            double differential_transverse_asymmetry_3(const double & s) const;
            double differential_transverse_asymmetry_4(const double & s) const;
            double differential_transverse_asymmetry_5(const double & s) const;
            double differential_transverse_asymmetry_re(const double & s) const;
            double differential_transverse_asymmetry_im(const double & s) const;
            double differential_p_prime_4(const double & s) const;
            double differential_p_prime_5(const double & s) const;
            double differential_p_prime_6(const double & s) const;
            double differential_h_1(const double & s) const;
            double differential_h_2(const double & s) const;
            double differential_h_3(const double & s) const;
            double differential_h_4(const double & s) const;
            double differential_h_5(const double & s) const;
            double differential_j_1s(const double & s) const;
            double differential_j_1c(const double & s) const;
            double differential_j_2s(const double & s) const;
            double differential_j_2c(const double & s) const;
            double differential_j_3(const double & s) const;
            double differential_j_3_normalized(const double & s) const;
            double differential_j_3_normalized_cp_averaged(const double & s) const;
            double differential_j_4(const double & s) const;
            double differential_j_5(const double & s) const;
            double differential_j_6s(const double & s) const;
            double differential_j_6c(const double & s) const;
            double differential_j_6c_cp_averaged(const double & s) const;
            double differential_j_7(const double & s) const;
            double differential_j_8(const double & s) const;
            double differential_j_9(const double & s) const;
            double differential_j_9_normalized(const double & s) const;
            double differential_j_9_normalized_cp_averaged(const double & s) const;
            double differential_j_1c_plus_j_2c_cp_averaged(const double & s) const;
            double differential_j_1s_minus_3j_2s_cp_averaged(const double & s) const;
            double differential_d_4(const double & s) const;
            double differential_d_5(const double & s) const;
            double differential_d_6s(const double & s) const;
            double differential_ratio_muons_electrons(const double & s) const;

            // Integrated Observables
            double integrated_decay_width(const double & s_min, const double & s_max) const;
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_branching_ratio_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_cp_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_isospin_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const;
            double integrated_longitudinal_polarisation_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_transversal_polarisation(const double & s_min, const double & s_max) const;
            double integrated_transversal_polarisation_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_2_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_5(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_re(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_im(const double & s_min, const double & s_max) const;
            double integrated_p_prime_4(const double & s_min, const double & s_max) const;
            double integrated_p_prime_5(const double & s_min, const double & s_max) const;
            double integrated_p_prime_6(const double & s_min, const double & s_max) const;
            double integrated_h_1(const double & s_min, const double & s_max) const;
            double integrated_h_2(const double & s_min, const double & s_max) const;
            double integrated_h_3(const double & s_min, const double & s_max) const;
            double integrated_h_4(const double & s_min, const double & s_max) const;
            double integrated_h_5(const double & s_min, const double & s_max) const;
            double integrated_j_1s(const double & s_min, const double & s_max) const;
            double integrated_j_1c(const double & s_min, const double & s_max) const;
            double integrated_j_2s(const double & s_min, const double & s_max) const;
            double integrated_j_2c(const double & s_min, const double & s_max) const;
            double integrated_j_3(const double & s_min, const double & s_max) const;
            double integrated_j_3_normalized(const double & s_min, const double & s_max) const;
            double integrated_j_3_normalized_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_j_4(const double & s_min, const double & s_max) const;
            double integrated_j_4_normalized_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_j_5(const double & s_min, const double & s_max) const;
            double integrated_j_5_normalized_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_j_6s(const double & s_min, const double & s_max) const;
            double integrated_j_6c(const double & s_min, const double & s_max) const;
            double integrated_j_7(const double & s_min, const double & s_max) const;
            double integrated_j_7_normalized_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_j_8(const double & s_min, const double & s_max) const;
            double integrated_j_8_normalized_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_j_9(const double & s_min, const double & s_max) const;
            double integrated_j_9_normalized(const double & s_min, const double & s_max) const;
            double integrated_j_9_normalized_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_a_9(const double & s_min, const double & s_max) const;
            double integrated_d_4(const double & s_min, const double & s_max) const;
            double integrated_d_5(const double & s_min, const double & s_max) const;
            double integrated_d_6s(const double & s_min, const double & s_max) const;
            double integrated_ratio_muons_electrons(const double & s_min, const double & s_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_s;
            static const std::string kinematics_description_c_theta_l;
            static const std::string kinematics_description_c_theta_k;
            static const std::string kinematics_description_phi;
    };

    /*
     * Decay: B -> K l lbar at Large Recoil, cf. [BFS2001], [BHP2007]
     */
    template <>
    class BToKDilepton<LargeRecoil> :
        public ParameterUser,
        public PrivateImplementationPattern<BToKDilepton<LargeRecoil>>
    {
        public:
            BToKDilepton(const Parameters & parameters, const Options & options);
            ~BToKDilepton();

            // Amplitudes
            std::complex<double> F_A(const double & s) const;
            std::complex<double> F_V(const double & s) const;
            std::complex<double> F_S(const double & s) const;
            std::complex<double> F_P(const double & s) const;
            std::complex<double> F_T(const double & s) const;
            std::complex<double> F_T5(const double & s) const;

            // Angular Observables
            double a_l(const double & s) const;
            double b_l(const double & s) const;
            double c_l(const double & s) const;

            // Two Differential Observables
            double two_differential_decay_width(const double & s, const double & c_theta_l) const;

            // Differential Observables
            double differential_branching_ratio(const double & s) const;
            double differential_flat_term(const double & s) const;
            double differential_forward_backward_asymmetry(const double & s) const;
            double differential_ratio_muons_electrons(const double & s) const;

            // Integrated Observables
            double integrated_decay_width(const double & s_min, const double & s_max) const;
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_branching_ratio_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_cp_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_flat_term(const double & s_min, const double & s_max) const;
            double integrated_flat_term_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_ratio_muons_electrons(const double & s_min, const double & s_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_s;
            static const std::string kinematics_description_c_theta_l;
    };
}

#endif
