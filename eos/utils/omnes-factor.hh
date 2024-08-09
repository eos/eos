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

#ifndef EOS_GUARD_EOS_UTILS_OMNES_FACTOR_HH
#define EOS_GUARD_EOS_UTILS_OMNES_FACTOR_HH 1

#include <eos/maths/complex.hh>
#include <eos/utils/exception.hh>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <array>
#include <functional>

namespace eos
{
    // Abstract class implementing the algorithm of [M:1999A] to solve the Omnes integral equation
    template <unsigned order_, unsigned nints_>
    class OmnesFactor
    {
        protected:
            // Vector containing integral borders
            std::array<double, nints_> _intervals;

            // Matrix holding the system of equations
            gsl_matrix * _sys;

            // Ingredients for SVD
            gsl_matrix * _sys2;
            gsl_matrix * _V;
            gsl_vector * _S;
            gsl_vector * _work;
            gsl_vector * _x;
            gsl_vector * _b;

            // Array containing Gauss-Legendre weights
            std::array<double, order_> _weights {};

            // Array containing zeros of Legendre Polynomials
            std::array<double, order_> _zeros {};

            // Array containing solution weights
            std::array<double, nints_ * order_> _sol {};

            // Error of approximation
            double _err;

            // Phase input
            std::function<double(const double &)> _scattering_phase;

            // Base constructor
            OmnesFactor(const std::array<double, nints_> & intervals, std::function<double(const double &)> scattering_phase);

            // Helper routines to compute ingredients for solution
            std::array<double, order_> lq_sum(const double & z) const;
            std::array<double, order_> p_ab(const double & z, const double & a, const double & b) const;
            std::array<double, order_> p_inf(const double & z, const double & a) const;
            std::array<double, order_> rr_ab(const double & z, const double & a, const double & b) const;
            std::array<double, order_> rr_inf(const double & z, const double & a) const;

            // Solve the system of equations
            std::array<double, nints_ * order_> solve_sys(const double & bc_pos);

            // Evaluate results
            complex<double> evaluate_omnes(const double & s) const;

        public:
            // Constructors
            OmnesFactor(const std::array<double, nints_> & intervals, std::function<double(const double &)> scattering_phase, const std::array<double, nints_ * order_> & sol) :
                OmnesFactor(intervals, scattering_phase) { _sol = sol; }

            OmnesFactor(const std::array<double, nints_> & intervals, std::function<double(const double &)> scattering_phase, const double & bcpos) :
                OmnesFactor(intervals, scattering_phase) { _sol = solve_sys(bcpos); }

            // Destructor
            ~OmnesFactor();

            // Return weights
            std::array<double, nints_ * order_> get_weights() { return _sol; }

            // Return Omnes factor evaluated at s
            complex<double> constexpr operator() (const double & s) const
            {
                double eps = 1e-7;
                for (unsigned i = 0 ; i < nints_ ; i++)
                {
                    if (std::abs(_intervals[i] - s) < eps)
                    {
                        // Can we be better here?
                        return (evaluate_omnes(s - eps) + evaluate_omnes(s + eps)) / 2.0;
                    }
                }
                return evaluate_omnes(s);
            }
    };

}

#endif
