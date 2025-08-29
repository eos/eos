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
#define EOS_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_HH 1

#include <eos/observable.hh>

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
             * e.g. from @f$\bar{B}^0 \to \bar{K}^{*0} \ell^+ \ell^-@f$, only.
             */
            // @{
            class IntermediateResult;
            const IntermediateResult * prepare(const double & q2_min, const double & q2_max) const;
            double integrated_decay_width(const IntermediateResult * ir) const;
            double integrated_branching_ratio(const IntermediateResult * ir) const;
            double integrated_unnormalized_forward_backward_asymmetry(const IntermediateResult * ir) const;
            double integrated_forward_backward_asymmetry(const IntermediateResult * ir) const;
            double integrated_longitudinal_polarisation(const IntermediateResult * ir) const;
            double integrated_transversal_polarisation(const IntermediateResult * ir) const;
            // @}

            /*!
             * @name Transverse asymmetries (@f$q^2@f$-integrated)
             *
             * These transverse asymmetries partially cancel hadronic
             * matrix elements at small @f$q^2@f$.
             */
            // @{
            double integrated_transverse_asymmetry_2(const IntermediateResult * ir) const;
            double integrated_transverse_asymmetry_3(const IntermediateResult * ir) const;
            double integrated_transverse_asymmetry_4(const IntermediateResult * ir) const;
            double integrated_transverse_asymmetry_5(const IntermediateResult * ir) const;
            double integrated_transverse_asymmetry_re(const IntermediateResult * ir) const;
            double integrated_transverse_asymmetry_im(const IntermediateResult * ir) const;
            // @}

            /*!
             * @name Optimized observables for low recoil (@f$q^2@f$-integrated)
             *
             * These transverse asymmetries partially cancel hadronic
             * matrix elements at large @f$q^2@f$.
             */
            // @{
            double integrated_h_1(const IntermediateResult * ir) const;
            double integrated_h_2(const IntermediateResult * ir) const;
            double integrated_h_3(const IntermediateResult * ir) const;
            double integrated_h_4(const IntermediateResult * ir) const;
            double integrated_h_5(const IntermediateResult * ir) const;
            // @}

            /*!
             * @name Angular observables (@f$q^2@f$-integrated)
             */
            // @{
            double integrated_j_1s(const IntermediateResult * ir) const;
            double integrated_j_1c(const IntermediateResult * ir) const;
            double integrated_j_2s(const IntermediateResult * ir) const;
            double integrated_j_2c(const IntermediateResult * ir) const;
            double integrated_j_3(const IntermediateResult * ir) const;
            double integrated_j_4(const IntermediateResult * ir) const;
            double integrated_j_5(const IntermediateResult * ir) const;
            double integrated_j_6s(const IntermediateResult * ir) const;
            double integrated_j_6c(const IntermediateResult * ir) const;
            double integrated_j_7(const IntermediateResult * ir) const;
            double integrated_j_8(const IntermediateResult * ir) const;
            double integrated_j_9(const IntermediateResult * ir) const;
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

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();

	         /*!
	          * Test functions, [BFS:2001] eqs. (40-41)
             */
            double real_C9_perp(const double & s) const;
            double real_C9_para(const double & s) const;
            double imag_C9_perp(const double & s) const;
            double imag_C9_para(const double & s) const;

            /*!
             * Test functions, return H_sb / H_c for each helicity
             */
            double H_perp_corrections(const double & s) const;
            double H_para_corrections(const double & s) const;
            double H_long_corrections(const double & s) const;
    };
}

#endif
