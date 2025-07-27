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

#ifndef EOS_GUARD_EOS_UTILS_RGE_IMPL_HH
#define EOS_GUARD_EOS_UTILS_RGE_IMPL_HH 1

#include <eos/utils/rge.hh>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

namespace eos
{
    template <unsigned nf_> struct QCDBetaFunction;

    template <> struct QCDBetaFunction<5u>
    {
            static constexpr double beta_0 = 23.0 / 3.0;
            static constexpr double beta_1 = 116.0 / 3.0;
    };

    template <unsigned nf_, unsigned dim_>
    MultiplicativeRenormalizationGroupEvolution<accuracy::LL, nf_, dim_>::MultiplicativeRenormalizationGroupEvolution(const std::array<double, dim_> &                   gamma_0_ev,
                                                                                                                      const std::array<std::array<double, dim_>, dim_> & V) :
        _gamma_0_ev(gamma_0_ev),
        _V(make_gsl_matrix(dim_, dim_)),
        _Vinv(make_gsl_matrix(dim_, dim_)),
        _U_0(make_gsl_matrix(dim_, dim_)),
        _c_0_0(make_gsl_vector(dim_)),
        _tmp_matrix(make_gsl_matrix(dim_, dim_)),
        _tmp_vector(make_gsl_vector(dim_))
    {
        for (unsigned i = 0; i < dim_; ++i)
        {
            for (unsigned j = 0; j < dim_; ++j)
            {
                gsl_matrix_set(_V.get(), i, j, V[i][j]);
            }
        }

        // tmp_matrix <- V
        gsl_matrix_memcpy(_tmp_matrix.get(), _V.get());

        gsl_permutation * p = gsl_permutation_alloc(dim_);
        int               signum;
        // tmp_matrix stores L, U from LU decomposition: P * V = L * U
        gsl_linalg_LU_decomp(_tmp_matrix.get(), p, &signum);
        // Vinv <- V^-1
        gsl_linalg_LU_invert(_tmp_matrix.get(), p, _Vinv.get());

        gsl_permutation_free(p);
    }

    template <unsigned nf_, unsigned dim_>
    std::array<double, dim_>
    MultiplicativeRenormalizationGroupEvolution<accuracy::LL, nf_, dim_>::evolve(const double & alpha_s_mu, const double & alpha_s_0, const std::array<double, dim_> & c_0_0) const
    {
        // LL evolution:
        //   c(mu) = U_0 . c(mu_0),
        // where
        //   U_0 = V . diag[ eta^(gamma_0_ev / (2 * beta_0)) ] . V^-1
        // since
        //   gamma_0 = V^-1,T . diag[ gamma_0_ev ] . V^T,

        const double eta    = alpha_s_0 / alpha_s_mu;
        const double beta_0 = QCDBetaFunction<nf_>::beta_0;

        gsl_matrix_set_identity(_U_0.get());

        // U_0 <- diag[ eta^(gamma_0_ev / (2 * beta_0)) ]
        // cf. [BBL:1995A], p. 34, eq. (III.94)
        for (unsigned i = 0; i < dim_; ++i)
        {
            gsl_matrix_set(_U_0.get(), i, i, std::pow(eta, _gamma_0_ev[i] / (2.0 * beta_0)));
            gsl_vector_set(_c_0_0.get(), i, c_0_0[i]);
        }
        // tmp_matrix <- V . U_0
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, _V.get(), _U_0.get(), 0.0, _tmp_matrix.get());
        // U_0 <- tmp_matrix . V
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, _tmp_matrix.get(), _Vinv.get(), 0.0, _U_0.get());

        // tmp_vector <- U_0 * c_0_0
        gsl_blas_dgemv(CblasNoTrans, 1.0, _U_0.get(), _c_0_0.get(), 0.0, _tmp_vector.get());

        std::array<double, dim_> result;
        for (unsigned i = 0; i < dim_; ++i)
        {
            result[i] = gsl_vector_get(_tmp_vector.get(), i);
        }

        return result;
    }

    template <unsigned nf_, unsigned dim_>
    MultiplicativeRenormalizationGroupEvolution<accuracy::NLL, nf_, dim_>::MultiplicativeRenormalizationGroupEvolution(const std::array<double, dim_> & gamma_0_ev,
                                                                                                                       const std::array<std::array<double, dim_>, dim_> & V,
                                                                                                                       const std::array<std::array<double, dim_>, dim_> & gamma_1) :
        _gamma_0_ev(gamma_0_ev),
        _V(make_gsl_matrix(dim_, dim_)),
        _Vinv(make_gsl_matrix(dim_, dim_)),
        _G(make_gsl_matrix(dim_, dim_)),
        _U_0(make_gsl_matrix(dim_, dim_)),
        _J(make_gsl_matrix(dim_, dim_)),
        _c_0_0(make_gsl_vector(dim_)),
        _c_0_1(make_gsl_vector(dim_)),
        _tmp_matrix(make_gsl_matrix(dim_, dim_)),
        _tmp_matrix_2(make_gsl_matrix(dim_, dim_)),
        _tmp_vector(make_gsl_vector(dim_)),
        _tmp_vector_2(make_gsl_vector(dim_))
    {
        // V <- V
        // G <- gamma_1
        for (unsigned i = 0; i < dim_; ++i)
        {
            for (unsigned j = 0; j < dim_; ++j)
            {
                gsl_matrix_set(_V.get(), i, j, V[i][j]);
                gsl_matrix_set(_G.get(), i, j, gamma_1[i][j]);
            }
        }

        // tmp_matrix <- V
        gsl_matrix_memcpy(_tmp_matrix.get(), _V.get());

        gsl_permutation * p = gsl_permutation_alloc(dim_);
        int               signum;
        // tmp_matrix stores L, U from LU decomposition: P . V = L . U
        gsl_linalg_LU_decomp(_tmp_matrix.get(), p, &signum);
        // Vinv <- V^-1
        gsl_linalg_LU_invert(_tmp_matrix.get(), p, _Vinv.get());

        gsl_permutation_free(p);

        // G = V^-1 . gamma_1^T . V, cf. [BBL:1995A], p. 34, eq. (III.96)
        // tmp_matrix <- V^-1 . gamma_1^T
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, _Vinv.get(), _G.get(), 0.0, _tmp_matrix.get());
        // G <- tmp_matrix . V
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, _tmp_matrix.get(), _V.get(), 0.0, _G.get());

        // J <- H = delta_ij gamma_0_ev_i beta_1 / (2 beta_0^2)
        //        - G_ij / (2 beta_0 + gamma_0_ev_i - gamma_0_ev_j)
        // cf. [BBL:1995A], p. 34, eq. (III.97)
        const double beta_0 = QCDBetaFunction<nf_>::beta_0;
        const double beta_1 = QCDBetaFunction<nf_>::beta_1;

        for (unsigned i = 0; i < dim_; ++i)
        {
            for (unsigned j = 0; j < dim_; ++j)
            {
                double value = -1.0 * gsl_matrix_get(_G.get(), i, j) / (2.0 * beta_0 + gamma_0_ev[i] - gamma_0_ev[j]);
                if (i == j)
                {
                    value += gamma_0_ev[i] * beta_1 / (2.0 * beta_0 * beta_0);
                }
                gsl_matrix_set(_J.get(), i, j, value);
            }
        }
        // tmp <- H . V^-1
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, _J.get(), _Vinv.get(), 0.0, _tmp_matrix.get());
        // J <- V . tmp
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, _V.get(), _tmp_matrix.get(), 0.0, _J.get());
    }

    template <unsigned nf_, unsigned dim_>
    std::array<double, dim_>
    MultiplicativeRenormalizationGroupEvolution<accuracy::NLL, nf_, dim_>::evolve(const double & alpha_s_mu, const double & alpha_s_0, const std::array<double, dim_> & c_0_0,
                                                                                  const std::array<double, dim_> & c_0_1) const
    {
        // NLL evolution:
        //   U = (1 + a_s_mu J) . U_0 . (1 - a_s_0 J)
        //   c(mu) = U . c_0(mu_0) + a_s_0 * U_0 . c_0_1(mu_0)
        // where
        //   a_s(x) = alpha_s(x) / (4 pi),
        //   U_0 = V . diag[ eta^(gamma_0_ev / (2 * beta_0)) ] . V^-1
        //   J = V . H . V^-1,
        // since
        //   gamma_0 = V^-1,T . diag[ gamma_0_ev ] . V^T,
        //   H = delta_ij gamma_0_ev_i beta_1 / (2 beta_0^2) - G_ij / (2 beta_0 + gamma_0_ev_i - gamma_0_ev_j),
        //   G = V^-1 . gamma_1^T . V

        const double eta    = alpha_s_0 / alpha_s_mu;
        const double beta_0 = QCDBetaFunction<nf_>::beta_0;

        gsl_matrix_set_identity(_U_0.get());

        // U_0 <- diag[ eta^(gamma_0_ev / (2 * beta_0)) ]
        // cf. [BBL:1995A], p. 34, eq. (III.94)
        for (unsigned i = 0; i < dim_; ++i)
        {
            gsl_matrix_set(_U_0.get(), i, i, std::pow(eta, _gamma_0_ev[i] / (2.0 * beta_0)));
            gsl_vector_set(_c_0_0.get(), i, c_0_0[i]);
            gsl_vector_set(_c_0_1.get(), i, c_0_1[i]);
        }
        // tmp_matrix <- V . U_0
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, _V.get(), _U_0.get(), 0.0, _tmp_matrix.get());
        // U_0 <- tmp_matrix * V^-1
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, _tmp_matrix.get(), _Vinv.get(), 0.0, _U_0.get());

        // tmp_vector := c_0_0 + alpha_s_0 / (4 pi) * (c_0_1 - J * c_0_0), cf. [BBL:1995A], p. 34, eq. (III.99)
        // tmp_vector <- c_0_1
        gsl_vector_memcpy(_tmp_vector.get(), _c_0_1.get());
        // tmp_vector <- -1.0 * J * c_0_0 + 1.0 * tmp_vector
        gsl_blas_dgemv(CblasNoTrans, -1.0, _J.get(), _c_0_0.get(), 1.0, _tmp_vector.get());
        // tmp_vector <- c_0_0 + alpha_s(mu_0) / (4 pi) * tmp_vector
        gsl_vector_axpby(1.0, _c_0_0.get(), alpha_s_0 / (4.0 * M_PI), _tmp_vector.get());

        // tmp_matrix <- unit_matrix
        gsl_matrix_set_zero(_tmp_matrix.get());
        gsl_matrix_set_identity(_tmp_matrix.get());
        // tmp_matrix2 <- alpha_s(mu) / (4 pi) * J
        gsl_matrix_memcpy(_tmp_matrix_2.get(), _J.get());
        gsl_matrix_scale(_tmp_matrix_2.get(), alpha_s_mu / (4.0 * M_PI));
        // tmp_matrix <- tmp_matrix + tmp_matrix_2
        gsl_matrix_add(_tmp_matrix.get(), _tmp_matrix_2.get());

        // tmp_vector_2 <- U_0 . tmp_vector
        gsl_blas_dgemv(CblasNoTrans, 1.0, _U_0.get(), _tmp_vector.get(), 0.0, _tmp_vector_2.get());
        // tmp_vector <- tmp_matrix . tmp_vector_2
        gsl_blas_dgemv(CblasNoTrans, 1.0, _tmp_matrix.get(), _tmp_vector_2.get(), 0.0, _tmp_vector.get());

        std::array<double, dim_> result;
        for (unsigned i = 0; i < dim_; ++i)
        {
            result[i] = gsl_vector_get(_tmp_vector.get(), i);
        }

        return result;
    }
} // namespace eos

#endif /* EOS_GUARD_EOS_UTILS_RGE_IMPL_HH */
