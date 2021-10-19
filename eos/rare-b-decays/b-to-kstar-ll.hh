/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2014 Christoph Bobeth
 * Copyright (c) 2016, 2017 Danny van Dyk
 * Copyright (c) 2021 Méril Reboud
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_B_TO_KSTAR_LL_HH 1

#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

#include <array>

namespace eos
{
    /*
     * Decay: B -> K^* l lbar.
     */
    class BToKstarDilepton :
        public ParameterUser,
        public PrivateImplementationPattern<BToKstarDilepton>
    {
        public:
            BToKstarDilepton(const Parameters & parameters, const Options & options);
            ~BToKstarDilepton();

            struct AngularCoefficients;
            struct Amplitudes;
            class AmplitudeGenerator;
            struct DipoleFormFactors;

            /*!
             * @name Signal PDFs
             */
            // @{
            double decay_width(const double & s, const double & c_theta_l, const double & c_theta_k, const double & phi) const;
            double decay_width_LHCb(const double & s, const double & c_theta_l_LHCb, const double & c_theta_k_LHCb, const double & phi_LHCb) const;
            // @}

            /*!
             * @name Inverse observables
             *
             * These observables have to be computed as inverse problems,
             * e.g. zero crossings.
             */
            // @{
            double a_fb_zero_crossing() const;
            // @}

            /*!
             * @name Simple observables (@f$q^2@f$-differential)
             *
             * These observables are differential in q^2, the dilepton
             * invariant mass. They arise form one decay mode,
             * e.g. from @f$\bar{B}^0 \to @f$\bar{K}^{*0} \ell^+ \ell^-@f$, only.
             */
            // @{
            double differential_decay_width(const double & q2) const;
            double differential_branching_ratio(const double & q2) const;
            double differential_forward_backward_asymmetry(const double & q2) const;
            double differential_longitudinal_polarisation(const double & q2) const;
            double differential_transversal_polarisation(const double & q2) const;
            // @}

            /*!
             * @name Transverse asymmetries (@f$q^2@f$-differential)
             *
             * These transverse asymmetries partially cancel hadronic
             * matrix elements at small @f$q^2@f$.
             */
            // @{
            double differential_transverse_asymmetry_2(const double & q2) const;
            double differential_transverse_asymmetry_3(const double & q2) const;
            double differential_transverse_asymmetry_4(const double & q2) const;
            double differential_transverse_asymmetry_5(const double & q2) const;
            double differential_transverse_asymmetry_re(const double & q2) const;
            double differential_transverse_asymmetry_im(const double & q2) const;
            // @}

            /*!
             * @name Optimized observables for large recoil (@f$q^2@f$-differential)
             *
             * These transverse asymmetries partially cancel hadronic
             * matrix elements at small @f$q^2@f$.
             */
            // @{
            double differential_p_prime_4(const double & q2) const;
            double differential_p_prime_5(const double & q2) const;
            double differential_p_prime_6(const double & q2) const;
            // @}

            /*!
             * @name Optimized observables for low recoil (@f$q^2@f$-differential)
             *
             * These transverse asymmetries partially cancel hadronic
             * matrix elements at large @f$q^2@f$.
             */
            // @{
            double differential_h_1(const double & q2) const;
            double differential_h_2(const double & q2) const;
            double differential_h_3(const double & q2) const;
            double differential_h_4(const double & q2) const;
            double differential_h_5(const double & q2) const;
            // @}

            /*!
             * @name Angular observables (@f$q^2@f$-differential)
             */
            // @{
            double differential_j_1s(const double & s) const;
            double differential_j_1c(const double & s) const;
            double differential_j_2s(const double & s) const;
            double differential_j_2c(const double & s) const;
            double differential_j_3(const double & s) const;
            double differential_j_4(const double & s) const;
            double differential_j_5(const double & s) const;
            double differential_j_6s(const double & s) const;
            double differential_j_6c(const double & s) const;
            double differential_j_6c_cp_averaged(const double & s) const;
            double differential_j_7(const double & s) const;
            double differential_j_8(const double & s) const;
            double differential_j_9(const double & s) const;
            double differential_j_1c_plus_j_2c_cp_averaged(const double & s) const;
            double differential_j_1s_minus_3j_2s_cp_averaged(const double & s) const;
            // @}

            /*!
             * @name Simple observables (@f$q^2@f$-integrated)
             *
             * These observables are integrated over q^2, the dilepton
             * invariant mass. They arise from one decay mode,
             * e.g. from @f$\bar{B}^0 \to \bar{K}^{*0} \ell^+ \ell^-@f$, only.
             */
            // @{
            double integrated_decay_width(const double & q2_min, const double & q2_max) const;
            double integrated_branching_ratio(const double & q2_min, const double & q2_max) const;
            double integrated_branching_ratio_cp_averaged(const double & q2_min, const double & q2_max) const;
            double integrated_unnormalized_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & q2_min, const double & q2_max) const;
            double integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_longitudinal_polarisation(const double & q2_min, const double & q2_max) const;
            double integrated_longitudinal_polarisation_cp_averaged(const double & q2_min, const double & q2_max) const;
            double integrated_transversal_polarisation(const double & q2_min, const double & q2_max) const;
            double integrated_transversal_polarisation_cp_averaged(const double & q2_min, const double & q2_max) const;
            // @}

            /*!
             * @name Direct CP-asymmetry (@f$q^2@f$-integrated)
             *
             * These observables are integrated over q^2, the dilepton
             * invariant mass. They arise from the CP-antisymmetrization
             * of the decay modes @f$\bar{B}^0 \to \bar{K}^{*0} \ell^+ \ell^-@f$
             * and @f$B^0 \to K^{*0} \ell^+ \ell^-@f$.
             */
            // @{
            double integrated_cp_asymmetry(const double & q2_min, const double & q2_max) const;
            // @}

            /*!
             * @name Transverse asymmetries (@f$q^2@f$-integrated)
             *
             * These transverse asymmetries partially cancel hadronic
             * matrix elements at small @f$q^2@f$.
             */
            // @{
            double integrated_transverse_asymmetry_2(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_2_cp_averaged(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_3(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_4(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_5(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_re(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_im(const double & q2_min, const double & q2_max) const;
            // @}

            /*!
             * @name Optimized observables for large recoil (@f$q^2@f$-integrated)
             *
             * These transverse asymmetries partially cancel hadronic
             * matrix elements at small @f$q^2@f$.
             */
            // @{
            double integrated_p_prime_4(const double & q2_min, const double & q2_max) const;
            double integrated_p_prime_5(const double & q2_min, const double & q2_max) const;
            double integrated_p_prime_6(const double & q2_min, const double & q2_max) const;
            // @}

            /*!
             * @name Optimized observables for low recoil (@f$q^2@f$-integrated)
             *
             * These transverse asymmetries partially cancel hadronic
             * matrix elements at large @f$q^2@f$.
             */
            // @{
            double integrated_h_1(const double & q2_min, const double & q2_max) const;
            double integrated_h_2(const double & q2_min, const double & q2_max) const;
            double integrated_h_3(const double & q2_min, const double & q2_max) const;
            double integrated_h_4(const double & q2_min, const double & q2_max) const;
            double integrated_h_5(const double & q2_min, const double & q2_max) const;
            // @}

            /*!
             * @name Angular observables (@f$q^2@f$-integrated)
             */
            // @{
            double integrated_j_1s(const double & q2_min, const double & q2_max) const;
            double integrated_j_1c(const double & q2_min, const double & q2_max) const;
            double integrated_j_2s(const double & q2_min, const double & q2_max) const;
            double integrated_j_2c(const double & q2_min, const double & q2_max) const;
            double integrated_j_3(const double & q2_min, const double & q2_max) const;
            double integrated_j_3_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_j_3_normalized_cp_averaged(const double & q2_min, const double & q2_max) const;
            double integrated_j_4(const double & q2_min, const double & q2_max) const;
            double integrated_j_5(const double & q2_min, const double & q2_max) const;
            double integrated_j_6s(const double & q2_min, const double & q2_max) const;
            double integrated_j_6c(const double & q2_min, const double & q2_max) const;
            double integrated_j_7(const double & q2_min, const double & q2_max) const;
            double integrated_j_8(const double & q2_min, const double & q2_max) const;
            double integrated_j_9(const double & q2_min, const double & q2_max) const;
            double integrated_j_9_normalized(const double & q2_min, const double & q2_max) const;
            double integrated_j_9_normalized_cp_averaged(const double & q2_min, const double & q2_max) const;
            // @}

            /*!
             * @name CP-symmetrized angular observables (@f$q^2@f$-integrated), in theory convention
             */
            // @{
            double integrated_s_1s(const double & q2_min, const double & q2_max) const;
            double integrated_s_1c(const double & q2_min, const double & q2_max) const;
            double integrated_s_2s(const double & q2_min, const double & q2_max) const;
            double integrated_s_2c(const double & q2_min, const double & q2_max) const;
            double integrated_s_3(const double & q2_min, const double & q2_max) const;
            double integrated_s_4(const double & q2_min, const double & q2_max) const;
            double integrated_s_5(const double & q2_min, const double & q2_max) const;
            double integrated_s_6s(const double & q2_min, const double & q2_max) const;
            double integrated_s_6c(const double & q2_min, const double & q2_max) const;
            double integrated_s_7(const double & q2_min, const double & q2_max) const;
            double integrated_s_8(const double & q2_min, const double & q2_max) const;
            double integrated_s_9(const double & q2_min, const double & q2_max) const;
            // @}

            /*!
             * @name CP-symmetrized angular observables (@f$q^2@f$-integrated), in LHCb convention
             */
            // @{
            double integrated_s_1s_LHCb(const double & q2_min, const double & q2_max) const
            {
                return integrated_s_1s(q2_min,q2_max);
            }
            double integrated_s_1c_LHCb(const double & q2_min, const double & q2_max) const
            {
                return integrated_s_1c(q2_min,q2_max);
            }
            double integrated_s_2s_LHCb(const double & q2_min, const double & q2_max) const
            {
                return integrated_s_2s(q2_min,q2_max);
            }
            double integrated_s_2c_LHCb(const double & q2_min, const double & q2_max) const
            {
                return integrated_s_2c(q2_min,q2_max);
            }
            double integrated_s_3_LHCb(const double & q2_min, const double & q2_max) const
            {
                return integrated_s_3(q2_min,q2_max);
            }
            double integrated_s_4_LHCb(const double & q2_min, const double & q2_max) const
            {
                return -integrated_s_4(q2_min,q2_max);
            }
            double integrated_s_5_LHCb(const double & q2_min, const double & q2_max) const
            {
                return integrated_s_5(q2_min,q2_max);
            }
            double integrated_s_6s_LHCb(const double & q2_min, const double & q2_max) const
            {
                return -integrated_s_6s(q2_min,q2_max);
            }
            double integrated_s_6c_LHCb(const double & q2_min, const double & q2_max) const
            {
                return  -integrated_s_6c(q2_min,q2_max);
            }
            double integrated_s_7_LHCb(const double & q2_min, const double & q2_max) const
            {
                return -integrated_s_7(q2_min,q2_max);
            }
            double integrated_s_8_LHCb(const double & q2_min, const double & q2_max) const
            {
                return integrated_s_8(q2_min,q2_max);
            }
            double integrated_s_9_LHCb(const double & q2_min, const double & q2_max) const
            {
                return -integrated_s_9(q2_min,q2_max);
            }
            double integrated_forward_backward_asymmetry_LHCb(const double & q2_min, const double & q2_max) const
            {
                return -integrated_forward_backward_asymmetry(q2_min,q2_max); 
            }

            /*!
             * @name CP-antisymmetrized angular observables (@f$q^2@f$-integrated)
             */
            // @{
            double integrated_a_1s(const double & q2_min, const double & q2_max) const;
            double integrated_a_1c(const double & q2_min, const double & q2_max) const;
            double integrated_a_2s(const double & q2_min, const double & q2_max) const;
            double integrated_a_2c(const double & q2_min, const double & q2_max) const;
            double integrated_a_3(const double & q2_min, const double & q2_max) const;
            double integrated_a_4(const double & q2_min, const double & q2_max) const;
            double integrated_a_5(const double & q2_min, const double & q2_max) const;
            double integrated_a_6s(const double & q2_min, const double & q2_max) const;
            double integrated_a_6c(const double & q2_min, const double & q2_max) const;
            double integrated_a_7(const double & q2_min, const double & q2_max) const;
            double integrated_a_8(const double & q2_min, const double & q2_max) const;
            double integrated_a_9(const double & q2_min, const double & q2_max) const;
            // @}

            /*!
             * Probes of symmetry relations in the large-energy limit (q^2 << m_b^2)
             */
            double differential_symrel_le_a1v(const double & q2) const;
            double differential_symrel_le_t1v(const double & q2) const;
            double differential_symrel_le_t2v(const double & q2) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            // @{
            static const std::string description;
            static const std::string kinematics_description_s;
            static const std::string kinematics_description_c_theta_l;
            static const std::string kinematics_description_c_theta_k;
            static const std::string kinematics_description_phi;
            // @}

            /*!
             * Auxilliary methods for unit tests and diagnostic purposes.
             */
            Amplitudes amplitudes(const double & q2) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;
    };

    /*!
     * Amplitudes for the decay B -> K^* l lbar.
     */
    struct BToKstarDilepton::Amplitudes
    {
        complex<double> a_long_right, a_long_left;
        complex<double> a_perp_right, a_perp_left;
        complex<double> a_para_right, a_para_left;
        complex<double> a_time, a_scal;
        complex<double> a_para_perp, a_time_long;
        complex<double> a_time_perp, a_long_perp;
        complex<double> a_time_para, a_long_para;
    };

    namespace btovll
    {
        struct Amplitudes
        {
            complex<double> a_long_right, a_long_left;
            complex<double> a_perp_right, a_perp_left;
            complex<double> a_par_right, a_par_left;
            complex<double> a_timelike;
            complex<double> a_scalar;
            complex<double> a_par_perp, a_t_long;
            complex<double> a_t_perp, a_long_perp;
            complex<double> a_t_par, a_long_par;
        };

        struct AngularCoefficients
        {
            double j1s, j1c;
            double j2s, j2c;
            double j3;
            double j4;
            double j5;
            double j6s, j6c;
            double j7;
            double j8;
            double j9;
        };

        inline AngularCoefficients array_to_angular_coefficients(const std::array<double, 12> & arr)
        {
            AngularCoefficients a_c = { arr[0], arr[1], arr[2], arr[3], arr[4],  arr[5],
                arr[6], arr[7], arr[8], arr[9], arr[10], arr[11] };

            return a_c;
        }

        inline double decay_width(const AngularCoefficients & a_c)
        {
            // cf. [BHvD2010], p. 6, eq. (2.7)
            return 2.0 * a_c.j1s + a_c.j1c - 1.0 / 3.0 * (2.0 * a_c.j2s + a_c.j2c);
        }

        inline std::array<double, 12> angular_coefficients_array(const Amplitudes & A, const double & s, const double & m_l)
        {
            // cf. [BHvD2010], p. 26, eqs. (A1)-(A11)
            // cf. [BHvD2012], app B, eqs. (B1)-(B12)
            std::array<double, 12> result;

            double z = 4.0 * power_of<2>(m_l) / s;
            double y = m_l / std::sqrt(s);
            double beta2 = 1.0 - z;
            double beta = std::sqrt(beta2);

            // j1s
            result[0] = 3.0 / 4.0 * (
                  (2.0 + beta2) / 4.0 * (norm(A.a_perp_left) + norm(A.a_perp_right) + norm(A.a_par_left) + norm(A.a_par_right))
                  + z * real(A.a_perp_left * conj(A.a_perp_right) + A.a_par_left * conj(A.a_par_right))
                  + 4.0 * beta2 * (norm(A.a_long_perp) + norm(A.a_long_par))
                  + 4.0 * (4.0 - 3.0 * beta2) * (norm(A.a_t_perp) + norm(A.a_t_par))
                  + 8.0 * std::sqrt(2.0) * y * real(
                       (A.a_par_left + A.a_par_right)   * conj(A.a_t_par)
                     + (A.a_perp_left + A.a_perp_right) * conj(A.a_t_perp)
                  )
               );
            // j1c
            result[1] = 3.0 / 4.0 * (
                  norm(A.a_long_left) + norm(A.a_long_right)
                  + z * (norm(A.a_timelike) + 2.0 * real(A.a_long_left * conj(A.a_long_right)))
                  + beta2 * norm(A.a_scalar)
                  + 8.0 * (2.0 - beta2) * norm(A.a_t_long)
                  + 8.0 * beta2 * norm(A.a_par_perp)
                  + 16.0 * y * real((A.a_long_left + A.a_long_right) * conj(A.a_t_long))
               );
            // j2s
            result[2] = 3.0 * beta2 / 16.0 * (
                  norm(A.a_perp_left) + norm(A.a_perp_right) + norm(A.a_par_left) + norm(A.a_par_right)
                  - 16.0 * (norm(A.a_t_perp) + norm(A.a_t_par) + norm(A.a_long_perp) + norm(A.a_long_par))
               );
            // j2c
            result[3] = -3.0 * beta2 / 4.0 * (
                  norm(A.a_long_left) + norm(A.a_long_right)
                  - 8.0 * (norm(A.a_t_long) + norm(A.a_par_perp))
               );
            // j3
            result[4] = 3.0 / 8.0 * beta2 * (
                  norm(A.a_perp_left) + norm(A.a_perp_right) - norm(A.a_par_left) - norm(A.a_par_right)
                  + 16.0 * (norm(A.a_t_par) - norm(A.a_t_perp) + norm(A.a_long_par) - norm(A.a_long_perp))
               );
            // j4
            result[5] = 3.0 / (4.0 * std::sqrt(2.0)) * beta2 * real(
                  A.a_long_left * conj(A.a_par_left) + A.a_long_right * conj(A.a_par_right)
                  - 8.0 * std::sqrt(2.0) * (A.a_t_long * conj(A.a_t_par) + A.a_par_perp * conj(A.a_long_par))
               );
            // j5
            result[6] = 3.0 * std::sqrt(2.0) / 4.0 * beta * real(
                  A.a_long_left * conj(A.a_perp_left) - A.a_long_right * conj(A.a_perp_right)
                  - 2.0 * std::sqrt(2.0) * A.a_t_par * conj(A.a_scalar)
                  - y * (
                     (A.a_par_left + A.a_par_right) * conj(A.a_scalar)
                     + 4.0 * std::sqrt(2.0) * A.a_long_par * conj(A.a_timelike)
                     - 4.0 * std::sqrt(2.0) * (A.a_long_left - A.a_long_right) * conj(A.a_t_perp)
                     - 4.0 * (A.a_perp_left - A.a_perp_right) * conj(A.a_t_long)
                  )
               );
            // j6s
            result[7] = 3.0 / 2.0 * beta * real(
                  A.a_par_left * conj(A.a_perp_left) - A.a_par_right * conj(A.a_perp_right)
                  + 4.0 * std::sqrt(2.0) * y * (
                     (A.a_perp_left - A.a_perp_right) * conj(A.a_t_par)
                     + (A.a_par_left - A.a_par_right) * conj(A.a_t_perp)
                  )
               );
            // j6c
            result[8] = 3.0 * beta * real(
                  2.0 * A.a_t_long * conj(A.a_scalar)
                  + y * (
                     (A.a_long_left + A.a_long_right) * conj(A.a_scalar)
                     + 4.0 * A.a_par_perp * conj(A.a_timelike)
                  )
               );
            // j7
            result[9] = 3.0 * std::sqrt(2.0) / 4.0 * beta * imag(
                  A.a_long_left * conj(A.a_par_left) - A.a_long_right * conj(A.a_par_right)
                  + 2.0 * std::sqrt(2.0) * A.a_t_perp * conj(A.a_scalar)
                  + y * (
                     (A.a_perp_left + A.a_perp_right) * conj(A.a_scalar)
                     + 4.0 * std::sqrt(2.0) * A.a_long_perp * conj(A.a_timelike)
                     + 4.0 * std::sqrt(2.0) * (A.a_long_left - A.a_long_right) * conj(A.a_t_par)
                     - 4.0 * (A.a_par_left - A.a_par_right) * conj(A.a_t_long)
                  )
               );
            // j8
            result[10] = 3.0 / 4.0 / std::sqrt(2.0) * beta2 * imag(
                  A.a_long_left * conj(A.a_perp_left) + A.a_long_right * conj(A.a_perp_right)
               );
            // j9
            result[11] = 3.0 / 4.0 * beta2 * imag(
                  conj(A.a_par_left) * A.a_perp_left + conj(A.a_par_right) * A.a_perp_right
               );

            return result;
        }
    }
}

#endif
