/*
 * Copyright (c) 2024 Florian Herren
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

#ifndef EOS_GUARD_EOS_UTILS_OMNES_FACTOR_IMPL_HH
#define EOS_GUARD_EOS_UTILS_OMNES_FACTOR_IMPL_HH 1

#include <eos/maths/complex.hh>
#include <eos/maths/legendre-polynomial-vector.hh>
#include <eos/maths/omnes-factor.hh>
#include <eos/utils/exception.hh>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <array>
#include <functional>

namespace eos
{

    // Allocate memory and initialise zeros & weights
    template <unsigned order_, unsigned nints_>
    OmnesFactor<order_, nints_>::OmnesFactor(const std::array<double, nints_> & intervals, std::function<double(const double &)> scattering_phase) :
        _intervals(intervals),
        _sys(gsl_matrix_calloc(nints_ * order_ + 1, nints_ * order_)),
        _sys2(gsl_matrix_calloc(nints_ * order_ + 1, nints_ * order_)),
        _V(gsl_matrix_calloc(nints_ * order_, nints_ * order_)),
        _S(gsl_vector_calloc(nints_ * order_)),
        _work(gsl_vector_calloc(nints_ * order_)),
        _x(gsl_vector_calloc(nints_ * order_)),
        _b(gsl_vector_calloc(nints_ * order_ + 1)),
        _err(-1.0),
        _eps(1.0e-5),
        _scattering_phase(scattering_phase)
    {
        // Perform pointer check
        if ((_sys == nullptr) || (_sys2 == nullptr) || (_V == nullptr) || (_S == nullptr) || (_work == nullptr) || (_x == nullptr) || (_b == nullptr)) { throw std::bad_alloc(); }

        // Obtain zeros of Legendre Polynomials and compute Gauss-Legendre weights
        LegendrePVector<order_> lp;
        lp.gauss_legendre(_zeros, _weights);

        // Compute values for evaluation
        for (unsigned i = 0 ; i < nints_ - 1 ; i++)
        {
            for (unsigned j = 0 ; j < order_ ; j++)
            {
                _slist[i * order_ + j] = (_intervals[i] + _intervals[i + 1] + (_intervals[i + 1] - _intervals[i]) * _zeros[j]) / 2.0;
                _tanvals[i][j] = std::tan(_scattering_phase(_slist[i * order_ + j]));
            }
        }

        for (unsigned j = 0 ; j < order_ ; j++)
        {
            _slist[(nints_ - 1) * order_ + j] = (2.0 * _intervals[nints_ - 1] / (1 - _zeros[j]));
            _tanvals[nints_ - 1][j] = std::tan(_scattering_phase(_slist[(nints_ - 1) * order_ + j]));
        }

        LegendrePVector<order_ - 1> lp_v;
        for(unsigned i = 0 ; i < order_ ; i++)
        {
            _p_j_u[i] = lp_v(_zeros[i]);
            for (unsigned j = 0 ; j < order_ ; j++)
            {
                _p_j_u[i][j] *= (2 * j + 1);
            }
            _weights[i] /= M_PI;
        }
    }

    // Deallocate memory
    template <unsigned order_, unsigned nints_>
    OmnesFactor<order_, nints_>::~OmnesFactor()
    {
        if (_sys) { gsl_matrix_free(_sys); }
        _sys = nullptr;

        if (_sys2) { gsl_matrix_free(_sys2); }
        _sys2 = nullptr;

        if (_V) { gsl_matrix_free(_V); }
        _V = nullptr;

        if (_S) { gsl_vector_free(_S); }
        _S = nullptr;

        if (_work) { gsl_vector_free(_work); }
        _work = nullptr;

        if (_x) { gsl_vector_free(_x); }
        _x = nullptr;

        if (_b) { gsl_vector_free(_b); }
        _b = nullptr;
    }

    // Prepare system of equations
    template <unsigned order_, unsigned nints_>
    std::array<double, nints_ * order_> OmnesFactor<order_, nints_>::solve_sys(const double & bc_pos)
    {
        std::array<double, (nints_ + 1) * order_> bc  {};
        std::array<double, nints_ * order_> sol {};

        for (unsigned i = 0 ; i < nints_ - 1 ; i++)
        {
            const std::array<double, order_> bcslice = rr_ab(bc_pos, _intervals[i], _intervals[i + 1]);
            for (unsigned j = 0 ; j < order_ ; j++)
            {
                bc[i * order_ + j] = bcslice[j];
            }
        }

        std::array<double, order_> bc_end = rr_inf(bc_pos, _intervals[nints_ - 1]);
        for (unsigned j = 0 ; j < order_ ; j++)
        {
            bc[(nints_ - 1) * order_ + j] = bc_end[j];
        }

        for (unsigned i = 0 ; i < nints_ * order_; i++)
        {
            for (unsigned j = 0 ; j < nints_ ; j++)
            {
                std::array<double, order_> row  {};
                if (j != nints_ - 1)
                {
                    row = rr_ab(_slist[i], _intervals[j], _intervals[j + 1]);
                }
                else
                {
                    row = rr_inf(_slist[i], _intervals[j]);
                }
                for (unsigned k = 0 ; k < order_ ; k++)
                {
                    if (i == j * order_ + k)
                    {
                        gsl_matrix_set(_sys, i, j * order_ + k, 1.0 - row[k]);
                    }
                    else
                    {
                        gsl_matrix_set(_sys, i, j * order_ + k, -row[k]);
                    }
                }
            }
        }

        for (unsigned j = 0 ; j < nints_ * order_ ; j++)
        {
            gsl_matrix_set(_sys, nints_ * order_, j, bc[j]);
        }

        // Initialize b & sys2
        gsl_vector_set(_b, nints_ * order_, 1.0);
        gsl_matrix_memcpy(_sys2, _sys);

        // Solve things, TODO: Check for error output....
        gsl_linalg_SV_decomp(_sys2, _V, _S, _work);
        gsl_linalg_SV_solve(_sys2, _V, _S, _b, _x);

        // Convert gsl -> vector
        for (unsigned i = 0 ; i < nints_ * order_ ; i++)
        {
            sol[i] = gsl_vector_get(_x, i);
        }

        gsl_blas_dgemv(CblasNoTrans, 1.0, _sys, _x, -1.0, _b);
        _err = gsl_blas_dnrm2(_b);

        return sol;
    }

    // Compute sums of Eq. 58 in [M:1999A]
    template <unsigned order_, unsigned nints_>
    std::array<double, order_> OmnesFactor<order_, nints_>::lq_sum(const double & z) const
    {
        std::array<double, order_> ret_vec  {};
        LegendreReQVector<order_ - 1> lq_v;
        const std::array<double, order_> q_j_z = lq_v(z);

        for (unsigned i = 0 ; i < order_ ; i++)
        {
            ret_vec[i] = 0.0;
            for (unsigned j = 0 ; j < order_ ; j++)
            {
                ret_vec[i] += _p_j_u[i][j] * q_j_z[j];
            }
        }
        return ret_vec;
    }

    // Compute generic part of integrand of Eq. 57 in [M:1999A]
    template <unsigned order_, unsigned nints_>
    std::array<double, order_> OmnesFactor<order_, nints_>::p_ab(const double & z, const double & a, const double & b) const
    {
        std::array<double, order_> ret_vec = lq_sum((2 * z - a - b) / (b - a));
        for (unsigned i = 0 ; i < order_ ; i++)
        {
            ret_vec[i] *= -_weights[i];
        }
        return ret_vec;
    }

    // Compute generic part of integrand of Eq. 60 in [M:1999A]
    template <unsigned order_, unsigned nints_>
    std::array<double, order_> OmnesFactor<order_, nints_>::p_inf(const double & z, const double & a) const
    {
        std::array<double, order_> ret_vec  {};
        if (std::abs(z) > 1e-10)
        {
            ret_vec = lq_sum(1.0 - 2.0 * a / z);
            for (unsigned i = 0 ; i < order_ ; i++)
            {
                ret_vec[i] *= -2.0 * a / z * _weights[i] / (1 - _zeros[i]);
            }
        }
        else
        {
            // Special case for z = 0
            for (unsigned i = 0 ; i < order_ ; i++)
            {
                ret_vec[i] = _p_j_u[i][0] * _weights[i] / (1 - _zeros[i]);
            }
        }
        return ret_vec;
    }

    // Compute full integrand of Eq. 57 in [M:1999A]
    template <unsigned order_, unsigned nints_>
    std::array<double, order_> OmnesFactor<order_, nints_>::rr_ab(const double & z, const double & a, const double & b) const
    {
        std::array<double, order_> ret_vec  {};
        const std::array<double, order_> pab = p_ab(z, a, b);

        for (unsigned i = 0 ; i < order_ ; i++)
        {
            ret_vec[i] = pab[i] * std::tan(_scattering_phase((a + b + (b - a) * _zeros[i]) / 2));
        }
        return ret_vec;
    }

    // Compute full integrand of Eq. 60 in [M:1999A]
    template <unsigned order_, unsigned nints_>
    std::array<double, order_> OmnesFactor<order_, nints_>::rr_inf(const double & z, const double & a) const
    {
        std::array<double, order_> ret_vec  {};
        const std::array<double, order_> pinf = p_inf(z, a);

        for (unsigned i = 0 ; i < order_ ; i++)
        {
            ret_vec[i] = pinf[i] * std::tan(_scattering_phase(2 * a / (1 - _zeros[i])));
        }
        return ret_vec;
    }

    template <unsigned order_, unsigned nints_>
    double OmnesFactor<order_, nints_>::omnes_abs(const double & s) const
    {
        double res = 0.0;

        for (unsigned i = 0 ; i < nints_ - 1 ; i++)
        {
            const std::array<double, order_> pab = p_ab(s, _intervals[i], _intervals[i + 1]);
            for (unsigned j = 0 ; j < order_ ; j++)
            {
                res += pab[j] * _tanvals[i][j] * _sol[i * order_ + j];
            }
        }

        const std::array<double, order_> pinf = p_inf(s, _intervals[nints_ - 1]);
        for (unsigned j = 0 ; j < order_ ; j++)
        {
            res += pinf[j] * _tanvals[nints_ - 1][j] * _sol[(nints_ - 1) * order_ + j];
        }

        return res;
    }

    template <unsigned order_, unsigned nints_>
    complex<double> OmnesFactor<order_, nints_>::evaluate_omnes(const double & s) const
    {
        if (s > _intervals[0])
        {
            // TODO: Can we do better?
            if (std::abs(_scattering_phase(s) - M_PI / 2.0) < 1e-8)
            {
                // Here we use that omnes_abs crosses 0 when delta becomes Pi/2 and cancel the simple pole from the tangent with the s-dependence of the derivative
                // We drop the denominators (12 * _eps) of abs_deriv and phase_deriv since they anyways cancel
                double abs_deriv = (-omnes_abs(s + 2.0 * _eps) + 8.0 * omnes_abs(s + _eps) - 8.0 * omnes_abs(s - _eps) + omnes_abs(s - 2.0 * _eps));
                double phase_deriv = (- _scattering_phase(s + 2.0 * _eps) + 8.0 * _scattering_phase(s + _eps) - 8.0 * _scattering_phase(s - _eps) + _scattering_phase(s - 2.0 * _eps));
                return complex<double>(0.0, -abs_deriv / phase_deriv);
            }
            else
            {
                return omnes_abs(s) * complex<double>(1.0, std::tan(_scattering_phase(s)));
            }
        }
        else
        {
            return complex<double>(omnes_abs(s), 0);
        }
    }
}

#endif
