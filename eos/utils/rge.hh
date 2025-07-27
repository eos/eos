/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_RGE_HH
#define EOS_GUARD_EOS_UTILS_RGE_HH 1

#include <eos/maths/gsl-interface.hh>

namespace eos
{

    template <typename Accuracy_, unsigned nf_, unsigned dim_> class MultiplicativeRenormalizationGroupEvolution;

    namespace accuracy
    {
        struct LL;
        struct NLL;
    } // namespace accuracy

    template <unsigned nf_, unsigned dim_> class MultiplicativeRenormalizationGroupEvolution<accuracy::LL, nf_, dim_>
    {
        private:
            // gamma_0 = V^-1,T . diag(gamma_0_ev) . V^T, see [BBL:1995A], p. 34, eq. (III.95)
            std::array<double, dim_> _gamma_0_ev;
            GSLMatrixPtr             _V, _Vinv;

            // evolution matrix
            GSLMatrixPtr _U_0;

            // initial conditions
            GSLVectorPtr _c_0_0;

            // temporary storage objects
            GSLMatrixPtr _tmp_matrix;
            GSLVectorPtr _tmp_vector;

        public:
            /*!
             * Constructor.
             *
             * This class expects provision with the anomalous mass dimension matrix (ADM) at LO only,
             * to provide RGE evolution to leading logarithmic accuracy. The LO gamma_0 matrix is
             * diagonalized by the matrix V, see [BBL:1995A], p. 34, eq. (III.95):
             *
             *   gamma_0 = V^-1,T . diag(gamma_0_ev) . V^T
             *
             * Note that, as in [BBL:1995A], the ADM for the operators is expected, not the ADM for
             * the Wilson coefficients, which is related to the operator ADM by transposition.
             *
             * @param gamma_0_ev The eigenvalues of the LO term for the anomalous dimension matrix.
             * @param V The matrix that diagonalizes the LO term for the anomalous dimension matrix.
             */
            MultiplicativeRenormalizationGroupEvolution(const std::array<double, dim_> & gamma_0_ev, const std::array<std::array<double, dim_>, dim_> & V);

            /*!
             * Evolve the Wilson coefficients from the scale mu_0 to the scale mu at leading-logarithmic accuracy.
             *
             * Expects the Wilson coefficients as a series in powers of alpha_(mu_0) / (4 pi):
             *
             *   c_0 = c_0_0 + O(alpha_(mu_0))
             *
             * @param alpha_s_mu The value of the strong coupling constant at the scale mu.
             * @param alpha_s_0 The value of the strong coupling constant at the scale mu_0.
             * @param c_0_0 The initial conditions for the Wilson coefficients at the scale mu_0 at order alpha_s^0
             */
            std::array<double, dim_> evolve(const double & alpha_s_mu, const double & alpha_s_0, const std::array<double, dim_> & c_0_0) const;
    };

    // Next-to-leading logarithmic accuracy, see [BBL:1995A], p. 34, eq. (III.93)
    template <unsigned nf_, unsigned dim_> class MultiplicativeRenormalizationGroupEvolution<accuracy::NLL, nf_, dim_>
    {
        private:
            // gamma_0 = V^-1,T . diag(gamma_0_ev) . V^T, see [BBL:1995A], p. 34, eq. (III.95)
            std::array<double, dim_> _gamma_0_ev;
            GSLMatrixPtr             _V, _Vinv;

            // gamma_1 = V^-1,T . G . V^T, see [BBL:1995A], p. 34, eq. (III.96)
            GSLMatrixPtr _G;

            // evolution matrices
            GSLMatrixPtr _U_0;
            GSLMatrixPtr _J;

            // initial conditions
            GSLVectorPtr _c_0_0;
            GSLVectorPtr _c_0_1;

            // temporary storage objects
            GSLMatrixPtr _tmp_matrix, _tmp_matrix_2;
            GSLVectorPtr _tmp_vector, _tmp_vector_2;

        public:
            /*!
             * Constructor.
             *
             * This class expects provision with the anomalous mass dimension matrix (ADM) at LO and NLO,
             * to provide RGE evolution to next-to-leading logarithmic accuracy. The LO gamma_0 matrix is
             * diagonalized by the matrix V, see [BBL:1995A], p. 34, eq. (III.95):
             *
             *   gamma_0 = V^-1,T . diag(gamma_0_ev) . V^T
             *
             * Note that, as in [BBL:1995A], the ADM for the operators is expected, not the ADM for the
             * Wilson coefficients, which is related to the operator ADM by transposition.
             *
             * @param gamma_0_ev The eigenvalues of the LO term for the anomalous dimension matrix.
             * @param V The matrix that diagonalizes the LO term for the anomalous dimension matrix.
             * @param gamma_1 The NLO term for the anomalous dimension matrix.
             */
            MultiplicativeRenormalizationGroupEvolution(const std::array<double, dim_> & gamma_0_ev, const std::array<std::array<double, dim_>, dim_> & V,
                                                        const std::array<std::array<double, dim_>, dim_> & gamma_1);

            /*!
             * Evolve the Wilson coefficients from the scale mu_0 to the scale mu at next-to-leading logarithmic accuracy.
             *
             * Expects the Wilson coefficients as a series in powers of alpha_(mu_0) / (4 pi):
             *
             *   c_0 = c_0_0 + alpha_(mu_0) / (4 pi) c_0_1 + O(alpha_(mu_0)^2)
             *
             * @param alpha_s_mu The value of the strong coupling constant at the scale mu.
             * @param alpha_s_0 The value of the strong coupling constant at the scale mu_0.
             * @param c_0_0 The initial conditions for the Wilson coefficients at the scale mu_0 at order alpha_s^0
             * @param c_0_1 The initial conditions for the Wilson coefficients at the scale mu_0 at order alpha_s^1,
             *              reduced by r^T . c_0_0, cf. [BBL:1995A], p. 34, eqs. (III.84) & (III.99).
             */
            std::array<double, dim_> evolve(const double & alpha_s_mu, const double & alpha_s_0, const std::array<double, dim_> & c_0_0,
                                            const std::array<double, dim_> & c_0_1) const;
    };
} // namespace eos

#endif /* EOS_GUARD_EOS_UTILS_RGE_HH */
