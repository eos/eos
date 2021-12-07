/*
 * Copyright (c) 2021 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_SZEGO_POLYNOMIALS_HH
#define EOS_GUARD_SRC_UTILS_SZEGO_POLYNOMIALS_HH 1

#include <eos/utils/complex.hh>
#include <eos/utils/power_of.hh>

#include <array>

namespace eos
{
    /*
     * Reprensentation of a Szego polynomial in the form of its
     * Verblunsky coefficients. We assume real-valued Vernblunsky
     * coefficients.
     */
    template <unsigned order_>
    class SzegoPolynomial
    {
        private:
            const double _norm_measure;

            const std::array<double, order_> _verblunsky_coefficients;

            const std::array<double, order_ + 1> _norms;

            std::array<double, order_ + 1> _calculate_norms() const
            {
                std::array<double, order_ + 1> result;

                result[0] = _norm_measure;

                for (unsigned i = 1 ; i <= order_ ; ++i)
                {
                    result[i] = result[i - 1] * (1.0 - power_of<2>(_verblunsky_coefficients[i - 1]));
                }

                for (unsigned i = 0 ; i <= order_ ; ++i)
                {
                    result[i] = std::sqrt(result[i]);
                }

                return result;
            }

        public:
            SzegoPolynomial(const double & norm_measure, std::array<double, order_> && verblunsky_coefficients) :
                _norm_measure(norm_measure),
                _verblunsky_coefficients(verblunsky_coefficients),
                _norms(_calculate_norms())
            {
            }

            SzegoPolynomial(const SzegoPolynomial &) = default;
            SzegoPolynomial(SzegoPolynomial &&) = default;
            ~SzegoPolynomial() = default;
            SzegoPolynomial & operator= (SzegoPolynomial &&) = default;
            SzegoPolynomial & operator= (const SzegoPolynomial &) = default;

            // Evaluate the normalized polynomials at on the real z axis, in the interval [-1, +1].
            // Note that, contrary to the literature [S:2004B], we use an integral measure dmu, which
            // yields \int dmu = `_norm_measure`, rather than the usual \int dmu = 1.
            std::array<double, order_ + 1> operator() (const double & z) const
            {
                std::array<double, order_ + 1> phi;
                std::array<double, order_ + 1> phi_star;

                phi[0]      = 1.0;
                phi_star[0] = 1.0;

                // we use real-valued Verblunsky coefficients only.
                for (unsigned n = 1 ; n <= order_ ; ++n)
                {
                    // cf. [S:2004B], eq. (1.4), p.2
                    phi[n]      = z * phi[n - 1]  - _verblunsky_coefficients[n - 1] * phi_star[n - 1];

                    // cf. [S:2004B], eqs. (1.4) and (1.5) in combination;
                    phi_star[n] = phi_star[n - 1] - _verblunsky_coefficients[n - 1] * z * phi[n - 1];
                }

                std::array<double, order_ + 1> result;
                for (unsigned n = 0 ; n <= order_ ; ++n)
                {
                    result[n] = phi[n] / _norms[n];
                };

                return result;
            }

            // Trivial generalization to real Verblunsky coefficients and complex z
            std::array<complex<double>, order_ + 1> operator() (const complex<double> & z) const
            {
                std::array<complex<double>, order_ + 1> phi;
                std::array<complex<double>, order_ + 1> phi_star;

                phi[0]      = 1.0;
                phi_star[0] = 1.0;

                // we use real-valued Verblunsky coefficients only.
                for (unsigned n = 1 ; n <= order_ ; ++n)
                {
                    // cf. [S:2004B], eq. (1.4), p.2
                    phi[n]      = z * phi[n - 1]  - _verblunsky_coefficients[n - 1] * phi_star[n - 1];

                    // cf. [S:2004B], eqs. (1.4) and (1.5) in combination;
                    phi_star[n] = phi_star[n - 1] - _verblunsky_coefficients[n - 1] * z * phi[n - 1];
                }

                std::array<complex<double>, order_ + 1> result;
                for (unsigned n = 0 ; n <= order_ ; ++n)
                {
                    result[n] = phi[n] / _norms[n];
                };

                return result;
            }
    };
}

#endif
