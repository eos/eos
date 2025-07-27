/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2014 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_QCD_HH
#define EOS_GUARD_EOS_UTILS_QCD_HH 1

#include <array>

namespace eos
{
    class QCD
    {
        public:
            /// The four leading coefficients of the \f$\alpha_s\f$ expansion of the quark mass' anomalous dimension \f$\gamma_m\f$.
            using AnomalousMassDimension = std::array<double, 4>;

            /// The four leading coefficient of the \f$\alpha_s\f$ expansion of the beta function of QCD.
            using BetaFunction = std::array<double, 4>;

            /*!
             * Calculate RGE running of strong coupling alpha_s from scale mu_0 down to scale mu in the MSbar scheme.
             *
             * Calculation according to [CKS2000].
             *
             * @param mu          scale at which alpha_s shall be evaluated
             * @param alpha_s_0   alpha_s at the initial scale mu_0
             * @param mu_0        initial scale mu_0
             * @param beta        parameters of QCD beta function that control the running
             */
            static double alpha_s(const double & mu, const double & alpha_s_0, const double & mu_0, const BetaFunction & beta, unsigned int loop_order = 4);

            /*!
             * Calculate RGE running of quark mass m_q in the MSbar scheme.
             *
             * Calculation according to [CKS2000].
             *
             * @param m_q         quark mass in the MSbar scheme at the scale m_q
             * @param alpha_s_0   alpha_s at the scale mu_0 = m_q
             * @param alpha_s_mu  alpha_s at the scale mu
             * @param beta        parameters of QCD beta function that control the running
             * @param gamma_m     parameters of QCD anomalous mass dimension that control the running
             */
            static double m_q_msbar(const double & m_q, const double & alpha_s_0, const double & alpha_s_mu, const BetaFunction & beta, const AnomalousMassDimension & gamma_m,
                                    unsigned int loop_order = 4);

            /*!
             * Calculate the shift from pole mass scheme to MSbar.
             *
             * Calculation according to [MvR1999].
             *
             * @param m_q_pole    quark mass in the pole mass scheme
             * @param alpha_s     alpha_s at scale mu = m_q_pole
             * @param nf          number of active QCD flavors that control the calculation
             */
            static double m_q_msbar(const double & m_q_pole, const double & alpha_s, const double & nf);

            /*!
             * Calculate the shift from MSbar scheme to pole mass.
             *
             * Calculation according to [CERN2003-002].
             *
             * @param m_q         quark mass in the MSbar scheme at the scale m_q
             * @param alpha_s     alpha_s at the scale m_q
             * @param nf          number of active QCD flavors that control the calculation
             */
            static double m_q_pole(const double & m_q, const double & alpha_s, const double & nf, unsigned int loop_order = 3);

            /*!
             * Calculate the shift from MSbar scheme to potential-Subtracted mass (PS mass).
             *
             * Calculation according to [B1998].
             *
             * @param m_q         quark mass in the MSbar scheme at the scale m_q
             * @param alpha_s     alpha_s at the scale m_q
             * @param mu_f        Factorization scale
             * @param nf          number of active QCD flavors that control the calculation
             * @param beta        parameters of QCD beta function that control the calculation
             */
            static double m_q_ps(const double & m_q, const double & alpha_s, const double & mu_f, const double & nf, const BetaFunction & beta);

            /*!
             * Calculate the shift from MSbar scheme to kinetic scheme (kin mass)
             *
             * Calculation according to [BBMU2003].
             *
             * @param m_q         quark mass in the MSbar scheme at the scale m_q
             * @param alpha_s     alpha_s at the scale m_q
             * @param mu          scale at which the mass shall be evaluated
             * @param beta        parameters of QCD beta function that control the calculation
             */
            static double m_q_kin(const double & m_q, const double & alpha_s, const double & mu, const BetaFunction & beta);

            /*!
             * The quadratic casimir operator for the fundamental representation of SU(3).
             */
            static const double casimir_f;

            /*!
             * The quadratic casimir operator for the adjoint representation SU(3).
             */
            static const double casimir_a;

            /*!
             * The 6-flavor-QCD constants.
             */
            ///@{
            static const BetaFunction beta_function_nf_6;
            ///@}

            /*!
             * The 5-flavor-QCD constants.
             */
            ///@{
            static const BetaFunction           beta_function_nf_5;
            static const AnomalousMassDimension gamma_m_nf_5;
            ///@}

            /*!
             * The 4-flavor-QCD constants.
             */
            ///@{
            static const BetaFunction           beta_function_nf_4;
            static const AnomalousMassDimension gamma_m_nf_4;
            ///@}

            /*!
             * The 3-flavor QCD constants.
             */
            ///@{
            static const BetaFunction           beta_function_nf_3;
            static const AnomalousMassDimension gamma_m_nf_3;
            ///@}
    };
} // namespace eos

#endif
