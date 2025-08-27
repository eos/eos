/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_BS_TO_PHI_LL_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_BS_TO_PHI_LL_HH 1

#define DECLARE_DIFFERENTIAL_H_FUNCTION(suffix) \
    double differential_H_##suffix(const double & s) const { \
        return differential_H(s, #suffix); \
    }
#define DECLARE_INTEGRATED_H_FUNCTION(suffix) \
    double integrated_H_##suffix(const double & s_min, const double & s_max) const { \
        return integrated_H(s_min, s_max, #suffix); \
    }

#define DECLARE_DIFFERENTIAL_Z_FUNCTION(suffix) \
    double differential_Z_##suffix(const double & s) const { \
        return differential_Z(s, #suffix); \
    }
#define DECLARE_INTEGRATED_Z_FUNCTION(suffix) \
    double integrated_Z_##suffix(const double & s_min, const double & s_max) const { \
        return integrated_Z(s_min, s_max, #suffix); \
    }

#define DECLARE_DIFFERENTIAL_S_FUNCTION(suffix) \
    double differential_S_##suffix(const double & s) const { \
        return differential_S(s, #suffix); \
    }
#define DECLARE_INTEGRATED_S_FUNCTION(suffix) \
    double integrated_S_##suffix(const double & s_min, const double & s_max) const { \
        return integrated_S(s_min, s_max, #suffix); \
    }

#define DECLARE_DIFFERENTIAL_K_FUNCTION(suffix) \
    double differential_K_##suffix(const double & s) const { \
        return differential_K(s, #suffix); \
    }
#define DECLARE_INTEGRATED_K_FUNCTION(suffix) \
    double integrated_K_##suffix(const double & s_min, const double & s_max) const { \
        return integrated_K(s_min, s_max, #suffix); \
    }

#define DECLARE_DIFFERENTIAL_A_FUNCTION(suffix) \
    double differential_A_##suffix(const double & s) const { \
        return differential_A(s, #suffix); \
    }
#define DECLARE_INTEGRATED_A_FUNCTION(suffix) \
    double integrated_A_##suffix(const double & s_min, const double & s_max) const { \
        return integrated_A(s_min, s_max, #suffix); \
    }

#define DECLARE_DIFFERENTIAL_W_FUNCTION(suffix) \
    double differential_W_##suffix(const double & s) const { \
        return differential_W(s, #suffix); \
    }
#define DECLARE_INTEGRATED_W_FUNCTION(suffix) \
    double integrated_W_##suffix(const double & s_min, const double & s_max) const { \
        return integrated_W(s_min, s_max, #suffix); \
    }

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

#include <array>

namespace eos
{
    /*
     * Decay: B_s -> phi l lbar.
     */
    class BsToPhiDilepton :
        public ParameterUser,
        public PrivateImplementationPattern<BsToPhiDilepton>
    {
        public:
            BsToPhiDilepton(const Parameters & parameters, const Options & options);
            ~BsToPhiDilepton();

            struct AngularCoefficients;
            struct Amplitudes;
            class AmplitudeGenerator;
            struct DipoleFormFactors;
            struct FormFactorCorrections;

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
             * e.g. from @f$\bar{B}_s^0 \to @f$\phi \ell^+ \ell^-@f$, only.
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
            double differential_j_7(const double & s) const;
            double differential_j_8(const double & s) const;
            double differential_j_9(const double & s) const;
            // @}

            /*!
             * @name Simple observables (@f$q^2@f$-integrated)
             *
             * These observables are integrated over q^2, the dilepton
             * invariant mass. They arise from one decay mode,
             * e.g. from @f$\bar{B}_s^0 \to \phi \ell^+ \ell^-@f$, only.
             */
            // @{
            double integrated_decay_width(const double & q2_min, const double & q2_max) const;
            double integrated_branching_ratio(const double & q2_min, const double & q2_max) const;
            double integrated_unnormalized_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & q2_min, const double & q2_max) const;
            double integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_longitudinal_polarisation(const double & q2_min, const double & q2_max) const;
            double integrated_transversal_polarisation(const double & q2_min, const double & q2_max) const;
            // @}

            /*!
             * @name Transverse asymmetries (@f$q^2@f$-integrated)
             *
             * These transverse asymmetries partially cancel hadronic
             * matrix elements at small @f$q^2@f$.
             */
            // @{
            double integrated_transverse_asymmetry_2(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_3(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_4(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_5(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_re(const double & q2_min, const double & q2_max) const;
            double integrated_transverse_asymmetry_im(const double & q2_min, const double & q2_max) const;
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
            double integrated_j_4(const double & q2_min, const double & q2_max) const;
            double integrated_j_5(const double & q2_min, const double & q2_max) const;
            double integrated_j_6s(const double & q2_min, const double & q2_max) const;
            double integrated_j_6c(const double & q2_min, const double & q2_max) const;
            double integrated_j_7(const double & q2_min, const double & q2_max) const;
            double integrated_j_8(const double & q2_min, const double & q2_max) const;
            double integrated_j_9(const double & q2_min, const double & q2_max) const;
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
             * Auxiliary methods for unit tests and diagnostic purposes.
             */
            Amplitudes amplitudes(const double & q2) const;
            double m_l() const;
            double phiBs() const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();

            AngularCoefficients pub_differential_angular_coefficients(const double & s) const;
            AngularCoefficients pub_integrated_angular_coefficients(const double & s_min, const double & s_max) const;

    	    /*!
             * Test functions, [BFS2001] eqs. (40-41)
             */
            double real_C9_perp(const double & s) const;
            double real_C9_para(const double & s) const;
            double imag_C9_perp(const double & s) const;
            double imag_C9_para(const double & s) const;
    };

    /*!
     * Amplitudes for the decay B_s -> phi l lbar.
     */
    struct BsToPhiDilepton::Amplitudes
    {
        complex<double> a_long_right, a_long_left;
        complex<double> a_perp_right, a_perp_left;
        complex<double> a_para_right, a_para_left;
        complex<double> a_time, a_scal;
        complex<double> a_para_perp, a_time_long;
        complex<double> a_time_perp, a_long_perp;
        complex<double> a_time_para, a_long_para;
    };


    class BsToPhiDileptonAndConjugate:
        public ParameterUser
    {
        public:
            BsToPhiDileptonAndConjugate(const Parameters & parameters, const Options & options);
            ~BsToPhiDileptonAndConjugate();

            BsToPhiDilepton bstophidilepton;
            BsToPhiDilepton bstophidilepton_conjugate;

            struct AngularhCoefficients;
            struct AngularsCoefficients;

            // Parameters related to mixing
            inline double decay_width(const BsToPhiDilepton::AngularCoefficients & a_c,
                                      const BsToPhiDilepton::AngularCoefficients & a_cc,
                                      const BsToPhiDileptonAndConjugate::AngularhCoefficients & a_h) const;
            double integrated_decay_width(const double & q2_min, const double & q2_max) const;
            double differential_decay_width(const double & q2) const;

            inline std::array<double, 12> angular_h_coefficients_array(const BsToPhiDilepton::Amplitudes & A, const BsToPhiDilepton::Amplitudes & Atilda, const double & s) const;
            inline std::array<double, 12> differential_angular_h_coefficients_array(const double & s) const;
            BsToPhiDileptonAndConjugate::AngularhCoefficients differential_angular_h_coefficients(const double & s) const;
            BsToPhiDileptonAndConjugate::AngularhCoefficients integrated_angular_h_coefficients(const double & s_min, const double & s_max) const;
            
            inline std::array<double, 12> angular_s_coefficients_array(const BsToPhiDilepton::Amplitudes & A, const BsToPhiDilepton::Amplitudes & Atilda, const double & s) const;
            inline std::array<double, 12> differential_angular_s_coefficients_array(const double & s) const;
            BsToPhiDileptonAndConjugate::AngularsCoefficients differential_angular_s_coefficients(const double & s) const;
            BsToPhiDileptonAndConjugate::AngularsCoefficients integrated_angular_s_coefficients(const double & s_min, const double & s_max) const;

            complex<double> a_long_right(const double & s) const;
            complex<double> a_long_left(const double & s) const;
            complex<double> a_perp_right(const double & s) const;
            complex<double> a_perp_left(const double & s) const;
            complex<double> a_para_right(const double & s) const;
            complex<double> a_para_left(const double & s) const;
            complex<double> a_time(const double & s) const;
            complex<double> a_scal(const double & s) const;

            double a_long_right_real(const double & s) const;
            double a_long_left_real(const double & s) const;
            double a_perp_right_real(const double & s) const;
            double a_perp_left_real(const double & s) const;
            double a_para_right_real(const double & s) const;
            double a_para_left_real(const double & s) const;
            double a_time_real(const double & s) const;
            double a_scal_real(const double & s) const;

            double a_long_right_imag(const double & s) const;
            double a_long_left_imag(const double & s) const;
            double a_perp_right_imag(const double & s) const;
            double a_perp_left_imag(const double & s) const;
            double a_para_right_imag(const double & s) const;
            double a_para_left_imag(const double & s) const;
            double a_time_imag(const double & s) const;
            double a_scal_imag(const double & s) const;

            double integrated_S(const double & s_min, const double & s_max, const std::string & name) const;
            DECLARE_INTEGRATED_S_FUNCTION(1s)
            DECLARE_INTEGRATED_S_FUNCTION(1c)
            DECLARE_INTEGRATED_S_FUNCTION(2s)
            DECLARE_INTEGRATED_S_FUNCTION(2c)
            DECLARE_INTEGRATED_S_FUNCTION(3)
            DECLARE_INTEGRATED_S_FUNCTION(4)
            DECLARE_INTEGRATED_S_FUNCTION(5)
            DECLARE_INTEGRATED_S_FUNCTION(6s)
            DECLARE_INTEGRATED_S_FUNCTION(6c)
            DECLARE_INTEGRATED_S_FUNCTION(7)
            DECLARE_INTEGRATED_S_FUNCTION(8)
            DECLARE_INTEGRATED_S_FUNCTION(9)

            double differential_S(const double & s, const std::string & name) const;
            DECLARE_DIFFERENTIAL_S_FUNCTION(1s)
            DECLARE_DIFFERENTIAL_S_FUNCTION(1c)
            DECLARE_DIFFERENTIAL_S_FUNCTION(2s)
            DECLARE_DIFFERENTIAL_S_FUNCTION(2c)
            DECLARE_DIFFERENTIAL_S_FUNCTION(3)
            DECLARE_DIFFERENTIAL_S_FUNCTION(4)
            DECLARE_DIFFERENTIAL_S_FUNCTION(5)
            DECLARE_DIFFERENTIAL_S_FUNCTION(6s)
            DECLARE_DIFFERENTIAL_S_FUNCTION(6c)
            DECLARE_DIFFERENTIAL_S_FUNCTION(7)
            DECLARE_DIFFERENTIAL_S_FUNCTION(8)
            DECLARE_DIFFERENTIAL_S_FUNCTION(9)

            double integrated_K(const double & s_min, const double & s_max, const std::string & name) const;
            DECLARE_INTEGRATED_K_FUNCTION(1s)
            DECLARE_INTEGRATED_K_FUNCTION(1c)
            DECLARE_INTEGRATED_K_FUNCTION(2s)
            DECLARE_INTEGRATED_K_FUNCTION(2c)
            DECLARE_INTEGRATED_K_FUNCTION(3)
            DECLARE_INTEGRATED_K_FUNCTION(4)
            DECLARE_INTEGRATED_K_FUNCTION(5)
            DECLARE_INTEGRATED_K_FUNCTION(6s)
            DECLARE_INTEGRATED_K_FUNCTION(6c)
            DECLARE_INTEGRATED_K_FUNCTION(7)
            DECLARE_INTEGRATED_K_FUNCTION(8)
            DECLARE_INTEGRATED_K_FUNCTION(9)

            double differential_K(const double & s, const std::string & name) const;
            DECLARE_DIFFERENTIAL_K_FUNCTION(1s)
            DECLARE_DIFFERENTIAL_K_FUNCTION(1c)
            DECLARE_DIFFERENTIAL_K_FUNCTION(2s)
            DECLARE_DIFFERENTIAL_K_FUNCTION(2c)
            DECLARE_DIFFERENTIAL_K_FUNCTION(3)
            DECLARE_DIFFERENTIAL_K_FUNCTION(4)
            DECLARE_DIFFERENTIAL_K_FUNCTION(5)
            DECLARE_DIFFERENTIAL_K_FUNCTION(6s)
            DECLARE_DIFFERENTIAL_K_FUNCTION(6c)
            DECLARE_DIFFERENTIAL_K_FUNCTION(7)
            DECLARE_DIFFERENTIAL_K_FUNCTION(8)
            DECLARE_DIFFERENTIAL_K_FUNCTION(9)

            double integrated_A(const double & s_min, const double & s_max, const std::string & name) const;
            DECLARE_INTEGRATED_A_FUNCTION(1s)
            DECLARE_INTEGRATED_A_FUNCTION(1c)
            DECLARE_INTEGRATED_A_FUNCTION(2s)
            DECLARE_INTEGRATED_A_FUNCTION(2c)
            DECLARE_INTEGRATED_A_FUNCTION(3)
            DECLARE_INTEGRATED_A_FUNCTION(4)
            DECLARE_INTEGRATED_A_FUNCTION(5)
            DECLARE_INTEGRATED_A_FUNCTION(6s)
            DECLARE_INTEGRATED_A_FUNCTION(6c)
            DECLARE_INTEGRATED_A_FUNCTION(7)
            DECLARE_INTEGRATED_A_FUNCTION(8)
            DECLARE_INTEGRATED_A_FUNCTION(9)

            double differential_A(const double & s, const std::string & name) const;
            DECLARE_DIFFERENTIAL_A_FUNCTION(1s)
            DECLARE_DIFFERENTIAL_A_FUNCTION(1c)
            DECLARE_DIFFERENTIAL_A_FUNCTION(2s)
            DECLARE_DIFFERENTIAL_A_FUNCTION(2c)
            DECLARE_DIFFERENTIAL_A_FUNCTION(3)
            DECLARE_DIFFERENTIAL_A_FUNCTION(4)
            DECLARE_DIFFERENTIAL_A_FUNCTION(5)
            DECLARE_DIFFERENTIAL_A_FUNCTION(6s)
            DECLARE_DIFFERENTIAL_A_FUNCTION(6c)
            DECLARE_DIFFERENTIAL_A_FUNCTION(7)
            DECLARE_DIFFERENTIAL_A_FUNCTION(8)
            DECLARE_DIFFERENTIAL_A_FUNCTION(9)

            double integrated_W(const double & s_min, const double & s_max, const std::string & name) const;
            DECLARE_INTEGRATED_W_FUNCTION(1s)
            DECLARE_INTEGRATED_W_FUNCTION(1c)
            DECLARE_INTEGRATED_W_FUNCTION(2s)
            DECLARE_INTEGRATED_W_FUNCTION(2c)
            DECLARE_INTEGRATED_W_FUNCTION(3)
            DECLARE_INTEGRATED_W_FUNCTION(4)
            DECLARE_INTEGRATED_W_FUNCTION(5)
            DECLARE_INTEGRATED_W_FUNCTION(6s)
            DECLARE_INTEGRATED_W_FUNCTION(6c)
            DECLARE_INTEGRATED_W_FUNCTION(7)
            DECLARE_INTEGRATED_W_FUNCTION(8)
            DECLARE_INTEGRATED_W_FUNCTION(9)

            double differential_W(const double & s, const std::string & name) const;
            DECLARE_DIFFERENTIAL_W_FUNCTION(1s)
            DECLARE_DIFFERENTIAL_W_FUNCTION(1c)
            DECLARE_DIFFERENTIAL_W_FUNCTION(2s)
            DECLARE_DIFFERENTIAL_W_FUNCTION(2c)
            DECLARE_DIFFERENTIAL_W_FUNCTION(3)
            DECLARE_DIFFERENTIAL_W_FUNCTION(4)
            DECLARE_DIFFERENTIAL_W_FUNCTION(5)
            DECLARE_DIFFERENTIAL_W_FUNCTION(6s)
            DECLARE_DIFFERENTIAL_W_FUNCTION(6c)
            DECLARE_DIFFERENTIAL_W_FUNCTION(7)
            DECLARE_DIFFERENTIAL_W_FUNCTION(8)
            DECLARE_DIFFERENTIAL_W_FUNCTION(9)

            double integrated_H(const double & q2_min, const double & q2_max, const std::string & name) const;
            DECLARE_INTEGRATED_H_FUNCTION(1s)
            DECLARE_INTEGRATED_H_FUNCTION(1c)
            DECLARE_INTEGRATED_H_FUNCTION(2s)
            DECLARE_INTEGRATED_H_FUNCTION(2c)
            DECLARE_INTEGRATED_H_FUNCTION(3)
            DECLARE_INTEGRATED_H_FUNCTION(4)
            DECLARE_INTEGRATED_H_FUNCTION(5)
            DECLARE_INTEGRATED_H_FUNCTION(6s)
            DECLARE_INTEGRATED_H_FUNCTION(6c)
            DECLARE_INTEGRATED_H_FUNCTION(7)
            DECLARE_INTEGRATED_H_FUNCTION(8)
            DECLARE_INTEGRATED_H_FUNCTION(9)

            double differential_H(const double & s, const std::string & name) const;
            DECLARE_DIFFERENTIAL_H_FUNCTION(1s)
            DECLARE_DIFFERENTIAL_H_FUNCTION(1c)
            DECLARE_DIFFERENTIAL_H_FUNCTION(2s)
            DECLARE_DIFFERENTIAL_H_FUNCTION(2c)
            DECLARE_DIFFERENTIAL_H_FUNCTION(3)
            DECLARE_DIFFERENTIAL_H_FUNCTION(4)
            DECLARE_DIFFERENTIAL_H_FUNCTION(5)
            DECLARE_DIFFERENTIAL_H_FUNCTION(6s)
            DECLARE_DIFFERENTIAL_H_FUNCTION(6c)
            DECLARE_DIFFERENTIAL_H_FUNCTION(7)
            DECLARE_DIFFERENTIAL_H_FUNCTION(8)
            DECLARE_DIFFERENTIAL_H_FUNCTION(9)

            double integrated_Z(const double & q2_min, const double & q2_max, const std::string & name) const;
            DECLARE_INTEGRATED_Z_FUNCTION(1s)
            DECLARE_INTEGRATED_Z_FUNCTION(1c)
            DECLARE_INTEGRATED_Z_FUNCTION(2s)
            DECLARE_INTEGRATED_Z_FUNCTION(2c)
            DECLARE_INTEGRATED_Z_FUNCTION(3)
            DECLARE_INTEGRATED_Z_FUNCTION(4)
            DECLARE_INTEGRATED_Z_FUNCTION(5)
            DECLARE_INTEGRATED_Z_FUNCTION(6c)
            DECLARE_INTEGRATED_Z_FUNCTION(6s)
            DECLARE_INTEGRATED_Z_FUNCTION(7)
            DECLARE_INTEGRATED_Z_FUNCTION(8)
            DECLARE_INTEGRATED_Z_FUNCTION(9)

            double differential_Z(const double & q2, const std::string & name) const;
            DECLARE_DIFFERENTIAL_Z_FUNCTION(1s)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(1c)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(2s)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(2c)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(3)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(4)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(5)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(6c)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(6s)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(7)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(8)
            DECLARE_DIFFERENTIAL_Z_FUNCTION(9)

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static const std::vector<OptionSpecification> options;
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();

        private:
            const double m_y;       // ΔΓ_s / Γ_s
            const double m_x;       // ΔM_s / Γ_s
            const double m_gamma_s; // Γ_s
    };

}

#endif
