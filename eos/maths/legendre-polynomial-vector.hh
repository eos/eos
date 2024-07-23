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

#ifndef EOS_GUARD_EOS_MATHS_LEGENDRE_POLYNOMIAL_VECTOR_HH
#define EOS_GUARD_EOS_MATHS_LEGENDRE_POLYNOMIAL_VECTOR_HH 1

#include <array>
#include <vector>
#include <complex>

#include <boost/math/special_functions/legendre.hpp>

namespace eos
{

    /*
     * Representation of a vector of Legendre polynomials for m = 0.
     *
     */
    template <unsigned order_>
    class LegendrePVector
    {
        public:
            LegendrePVector()
            {
            }

            LegendrePVector(const LegendrePVector &) = default;
            LegendrePVector(LegendrePVector &&) = default;
            ~LegendrePVector() = default;
            LegendrePVector & operator= (LegendrePVector &&) = default;
            LegendrePVector & operator= (const LegendrePVector &) = default;

            // Evaluate the vector of Legendre polynomials
            std::array<double, order_ + 1> constexpr operator() (const double & z) const
            {
                std::array<double, order_ + 1> ret_vec;

                if (order_ == 0)
                {
                    ret_vec[0] = 1;
                }
                else if (order_ == 1)
                {
                    ret_vec[0] = 1;
                    ret_vec[1] = z;
                }
                else
                {
                    ret_vec[0] = 1;
                    ret_vec[1] = z;
                    for (unsigned i = 2 ; i <= order_ ; i++)
                    {
                        ret_vec[i] = (z * ret_vec[i - 1] * (2 * i - 1) - ret_vec[i - 2] * (i - 1)) / i;
                    }
                }
                return ret_vec;
            }

            // Return zeros and compute Gauss-Legendre weights
            void gauss_legendre(std::array<double, order_> & zeros, std::array<double, order_> & weights)
            {
                std::vector<double> zerosp = boost::math::legendre_p_zeros<double>(order_);
                std::array<double, order_ / 2 + order_ % 2> weightsp;
                const unsigned len = order_ / 2 + order_ % 2;

                for (unsigned i = 0 ; i < len ; i++)
                {
                    double lp = boost::math::legendre_p<double>(order_ + 1, zerosp[i]);
                    weightsp[i] = 2 * (1 - zerosp[i] * zerosp[i]) / ((order_ + 1) * (order_ + 1) * lp * lp);
                }

                // Flip order for first half
                for (unsigned i = 0 ; i < len ; i++)
                {
                    zeros[i] = -1 * zerosp[len - 1 - i];
                    weights[i] = weightsp[len - 1 - i];
                }

                // Double them, take care in case of odd order_
                if (order_ % 2 == 1)
                {
                    for (unsigned i = 1 ; i < len ; i++)
                    {
                        zeros[len + i - 1] = zerosp[i];
                        weights[len + i - 1] = weightsp[i];
                    }
                }
                else
                {
                    for (unsigned i = 0 ; i < len ; i++)
                    {
                        zeros[len + i] = zerosp[i];
                        weights[len + i] = weightsp[i];
                    }
                }
            }
    };

    /*
     * Representation of a vector of the real part of the associated Legendre functions for m = 0.
     * Implementation is based on Zhang, Shanjie and Jin, Jianming. “Computation of Special Functions”, John Wiley and Sons, 1996.
     *
     */
    template <unsigned order_>
    class LegendreReQVector
    {
        private:
            const double _eps;
            const unsigned _cut;

        public:
            // 10^-14 chosen to stay away from double precision noise
            // 1000 iterations are enough to cover all practical cases,
            // but keep us from going into an infinite loop if we really
            // hit a point where things do not converge.
            LegendreReQVector() :
                _eps(1e-14),
                _cut(1000)
            {
            }

            LegendreReQVector(const double & eps, const unsigned & cut) :
                _eps(eps),
                _cut(cut)
            {
            }

            LegendreReQVector(const LegendreReQVector &) = default;
            LegendreReQVector(LegendreReQVector &&) = default;
            ~LegendreReQVector() = default;
            LegendreReQVector & operator= (LegendreReQVector &&) = default;
            LegendreReQVector & operator= (const LegendreReQVector &) = default;

            // Evaluate the vector of associated Legendre functions
            // Throws exception if algorithm doesn't terminate
            // Returns real part if abs(z) > 1
            std::array<double, order_ + 1> constexpr operator() (const double & z) const
            {
                std::array<double, order_ + 1> ret_vec;

                // Standard forward recursion
                if (std::abs(z) < 1)
                {
                    if (order_ == 0)
                    {
                        ret_vec[0] = std::log((1 + z) / (1 - z)) / 2;
                    }
                    else if (order_ == 1)
                    {
                        ret_vec[0] = std::log((1 + z) / (1 - z)) / 2;
                        ret_vec[1] = z * ret_vec[0] - 1;
                    }
                    else
                    {
                        ret_vec[0] = std::log((1 + z) / (1 - z)) / 2;
                        ret_vec[1] = z * ret_vec[0] - 1;
                        for (unsigned i = 2 ; i <= order_ ; i++)
                        {
                            ret_vec[i] = (z * ret_vec[i - 1] * (2 * i - 1) - ret_vec[i - 2] * (i - 1)) / i;
                        }
                    }
                }
                // Heuristic parameter for when we need the backward recursion
                else if (std::abs(z) >= 1.021)
                {
                    if (order_ == 0)
                    {
                        ret_vec[0] = std::log(std::abs((1 + z) / (1 - z))) / 2;
                    }
                    else if (order_ == 1)
                    {
                        ret_vec[0] = std::log(std::abs((1 + z) / (1 - z))) / 2;
                        ret_vec[1] = z * ret_vec[0] - 1;
                    }
                    else
                    {
                        double q2 = 1 / z;
                        double q1 = 1;
                        for (unsigned i = 1 ; i <= order_ ; i++)
                        {
                            q2 *= i / z / (2 * i + 1);
                            if (i == order_ - 1)
                            {
                                q1 = q2;
                            }
                        }

                        double t1 = 1;
                        double qr = 1;
                        for (unsigned k = 1 ; k <= _cut ; k++)
                        {
                            qr *= (order_ / 2.0 + k - 1) * (k + (order_ - 1) / 2.0) / (z * z * k * ((order_ + k) - 1.0 / 2.0));
                            t1 += qr;
                            if (std::abs(qr / t1) < _eps)
                            {
                                break;
                            }
                            if (k == _cut)
                            {
                                throw InternalError("Maximum number of iterations reached in LegendreReQVector!");
                            }
                        }
                        ret_vec[order_ - 1] = t1 * q1;

                        double t2 = 1;
                        qr = 1;
                        for (unsigned k = 1 ; k <= _cut ; k++)
                        {
                            qr *= ((order_ + 1) / 2.0 + k - 1) * (k + order_ / 2.0) / (z * z * k * ((order_ + 1 + k) - 1.0 / 2.0));
                            t2 += qr;
                            if (std::abs(qr / t2) < _eps)
                            {
                                break;
                            }
                            if (k == _cut)
                            {
                                throw InternalError("Maximum number of iterations reached in LegendreReQVector!");
                            }
                        }
                        ret_vec[order_] = t2 * q2;

                        for (unsigned i = order_ ; i >= 4 ; i--)
                        {
                            ret_vec[i - 2] = ((2 * i - 1) * z * ret_vec[i - 1] - i * ret_vec[i]) / (i - 1);
                        }
                        ret_vec[0] = std::log(std::abs((1 + z) / (1 - z))) / 2;
                        ret_vec[1] = z * ret_vec[0] - 1;
                    }
                }
                else
                {
                    if (order_ == 0)
                    {
                        ret_vec[0] = std::log(std::abs((1 + z) / (1 - z))) / 2;
                    }
                    else if (order_ == 1)
                    {
                        ret_vec[0] = std::log(std::abs((1 + z) / (1 - z))) / 2;
                        ret_vec[1] = z * ret_vec[0] - 1;
                    }
                    else
                    {
                        ret_vec[0] = std::log(std::abs((1 + z) / (1 - z))) / 2;
                        ret_vec[1] = z * ret_vec[0] - 1;
                        for (unsigned i = 2 ; i <= order_ ; i++)
                        {
                            ret_vec[i] = (z * ret_vec[i - 1] * (2 * i - 1) - ret_vec[i - 2] * (i - 1)) / i;
                        }
                    }
                }
                return ret_vec;
            }
    };
}

#endif
