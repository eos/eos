/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Fatemeh Nouri
 * Copyright (c) 2026 Méril Reboud
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BHKMNR2026_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BHKMNR2026_HH 1

#include <eos/form-factors/mesonic-processes.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/reference-name.hh>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_poly.h>

#include <array>

namespace eos
{
    template <typename Process_> class BHKMNR2026FormFactors;

    template <> class BHKMNR2026FormFactors<VacuumToPiPi> : public FormFactors<VacuumToPP>
    {
        private:
            // parameters for form factor f_+ (I=1 projection)
            std::array<UsedParameter, 9u> _a_fp_I1; // unconstrained expansion coefficients
            std::array<UsedParameter, 3u> _M_fp_I1; // masses of the rho, rho', etc.
            std::array<UsedParameter, 3u> _G_fp_I1; // widths of the rho, rho', etc.

            RestrictedOption _n_resonances_I1; // number of used I = 1 resonances

            // parameters for form factor f_+ (I=0 projection)
            std::array<UsedParameter, 5u> _a_fp_I0; // unconstrained expansion coefficients
            std::array<UsedParameter, 2u> _M_fp_I0; // masses of the omega, phi, etc.
            std::array<UsedParameter, 2u> _G_fp_I0; // widths of the omega, phi, etc.

            RestrictedOption _n_resonances_I0; // number of used I = 0 resonances

            UsedParameter _m_pi; // pion mass

            UsedParameter _s_0;  // free parameter
            UsedParameter _s_in; // inelastic threshold

            UsedParameter _hbar;

            // Isospin option
            std::array<bool, 2> _switch_I;
            IsospinOption       _opt_I;

            // matrices and vectors needed for the linear system of equations to determine the constrained coefficients
            gsl_matrix *                 _M;
            gsl_matrix *                 _inv_M;
            gsl_vector *                 _L;
            gsl_permutation *            _perm;
            gsl_vector *                 _constrained_coefficents;
            gsl_poly_complex_workspace * _poly_workspace;

            inline std::string
            _par_name(const std::string & ff, const std::string & isospin, const std::string & index) const
            {
                return "0->pipi::a_(" + ff + "," + isospin + ")^" + index + "@BHKMNR2026";
            }

            inline complex<double>
            _s_p() const
            {
                return complex<double>(4.0 * power_of<2>(_m_pi()), 0.0); // pair-production threshold s_plus
            }

            inline complex<double>
            _s_m() const
            {
                return complex<double>(0.0, 0.0); // start of the left-hand cut s_minus
            }

            inline complex<double>
            _phi_to_s(const complex<double> & phi, const complex<double> & s_in) const
            {
                const complex<double> s_p = _s_p();

                return (s_p * power_of<2>(1.0 + power_of<2>(phi)) - 4.0 * s_in * power_of<2>(phi)) / power_of<2>(1.0 - power_of<2>(phi));
            }

            inline complex<double>
            _s_to_phi_11(const complex<double> & s, const complex<double> & s_in) const
            {
                const complex<double> s_p = _s_p();
                const double          eps = 1e-14;

                // phi_11(s = s_+) = 0 by construction
                if (std::abs(s - s_p) < eps)
                {
                    return complex<double>(0.0, 0.0);
                }

                return (std::sqrt(s_in - s) - std::sqrt(s_in - s_p)) / std::sqrt(s_p - s);
            }

            inline complex<double>
            _s_to_phi_21(const complex<double> & s, const complex<double> & s_in) const
            {
                const complex<double> s_p = _s_p();
                static const double   eps = 1.0e-14;

                // phi_21(s = s_+) = 0 by construction
                if (std::abs(s - s_p) < eps)
                {
                    return complex<double>(0.0, 0.0);
                }

                return -std::sqrt(s_p - s) / (std::sqrt(s_in - s_p) + std::sqrt(s_in - s));
            }

            inline complex<double>
            _s_to_phi_22(const complex<double> & s, const complex<double> & s_in) const
            {
                const complex<double> s_p = _s_p();
                const double          eps = 1e-14;

                // phi_22(s = s_+) -> infinity
                if (std::abs(s - s_p) < eps)
                {
                    throw InternalError("phi_22 is not finite at s = s_plus");
                }

                return (std::sqrt(s_in - s) + std::sqrt(s_in - s_p)) / std::sqrt(s_p - s);
            }

            inline complex<double>
            _s_to_phi_12(const complex<double> & s, const complex<double> & s_in) const
            {
                const complex<double> s_p = _s_p();
                const double          eps = 1e-14;

                // phi_12(s = s_+) -> infinity
                if (std::abs(s - s_p) < eps)
                {
                    throw InternalError("phi_12 is not finite at s = s_plus");
                }

                return (-std::sqrt(s_in - s) - std::sqrt(s_in - s_p)) / std::sqrt(s_p - s);
            }

            // map that deals with the LHC of the RS 21, c.f. eq. (4.3) of arxiv:2510.25584
            inline complex<double>
            _chi(const complex<double> & x, const complex<double> & x_L, const complex<double> & x_0) const
            {
                const complex<double> A = (x * power_of<2>(x_L - 1.0) - x_L * power_of<2>(x - 1.0)) * power_of<2>(x_0 - 1.0);
                const complex<double> B = (x_0 * power_of<2>(x_L - 1.0) - x_L * power_of<2>(x_0 - 1.0)) * power_of<2>(x - 1.0);

                return (std::sqrt(A) - std::sqrt(B)) / (std::sqrt(A) + std::sqrt(B));
            }

            inline complex<double>
            _chi_inverse(const complex<double> & y, const complex<double> & x_L, const complex<double> & x_0) const
            {
                const complex<double> A = power_of<2>(y - 1.0) * power_of<2>(x_L - 1.0) * power_of<2>(x_0 - 1.0);
                const complex<double> B = 4.0 * x_0 * power_of<2>(y + 1.0) * power_of<2>(x_L - 1.0) - 16.0 * y * x_L * power_of<2>(x_0 - 1.0);

                return (std::sqrt(A + B) - std::sqrt(A)) / (std::sqrt(A + B) + std::sqrt(A));
            }

            inline complex<double>
            _s_to_psi_11(const complex<double> & s) const
            {
                const complex<double> s_0   = _s_0();
                const complex<double> s_in  = _s_in();
                const complex<double> phi_L = _s_to_phi_21(_s_m(), s_in);

                return _chi(_s_to_phi_11(s, s_in), phi_L, _s_to_phi_11(s_0, s_in));
            }

            inline complex<double>
            _s_to_psi_21(const complex<double> & s) const
            {
                const complex<double> s_0   = _s_0();
                const complex<double> s_in  = _s_in();
                const complex<double> phi_L = _s_to_phi_21(_s_m(), s_in);

                return _chi(_s_to_phi_21(s, s_in), phi_L, _s_to_phi_11(s_0, s_in));
            }

            inline complex<double>
            _s_to_psi_22(const complex<double> & s) const
            {
                const complex<double> s_0   = _s_0();
                const complex<double> s_in  = _s_in();
                const complex<double> phi_L = _s_to_phi_21(_s_m(), s_in);

                return 1.0 / _chi(1.0 / _s_to_phi_22(s, s_in), phi_L, _s_to_phi_11(s_0, s_in));
            }

            inline complex<double>
            _s_to_psi_12(const complex<double> & s) const
            {
                const complex<double> s_0   = _s_0();
                const complex<double> s_in  = _s_in();
                const complex<double> phi_L = _s_to_phi_21(_s_m(), s_in);

                return 1.0 / _chi(1.0 / _s_to_phi_12(s, s_in), phi_L, _s_to_phi_11(s_0, s_in));
            }

            inline complex<double>
            _psi_r(const double & M, const double & Gamma) const
            {
                if (M * M < _s_in()) // the resonance is below the inelastic threshold, the pole is on the 21 Riemann sheet
                {
                    return _s_to_psi_21(power_of<2>(complex<double>(M, -Gamma / 2.0)));
                }
                else // the resonance is above the inelastic threshold, we only consider the pole on the 22 Riemann sheet
                {
                    return _s_to_psi_22(power_of<2>(complex<double>(M, -Gamma / 2.0)));
                }
            }

            // c.f. eq. (4.10) of arxiv:2510.25584, but we also factorize (psi - 1.0)^2 * (psi + 1.0) to
            // fix the threshold and high energy behaviors of the form factor
            inline complex<double>
            _P(const complex<double> & psi) const
            {
                complex<double>   psi_r;
                const std::size_t num_resonances = stoi(_n_resonances_I1.value());
                complex<double>   result         = power_of<2>(psi - 1.0) * (psi + 1.0);

                for (auto i = 0u; i < num_resonances; i++)
                {
                    psi_r   = _psi_r(_M_fp_I1[i](), _G_fp_I1[i]());
                    result /= (psi - psi_r) * (psi - std::conj(psi_r));
                }

                return result;
            }

            inline complex<double>
            _Q(const complex<double> & psi) const
            {
                complex<double>   psi_r;
                const std::size_t num_resonances = stoi(_n_resonances_I0.value());
                complex<double>   result         = 1.0;

                for (auto i = 0u; i < num_resonances; i++)
                {
                    psi_r   = _psi_r(_M_fp_I0[i](), _G_fp_I0[i]());
                    result /= (psi - psi_r) * (psi - std::conj(psi_r));
                }

                return result;
            }

            // Residue of P (in psi) at the k^th resonance
            inline complex<double>
            _P_residue(const unsigned & k) const
            {
                const std::size_t num_resonances_I1 = stoi(_n_resonances_I1.value());

                if (k >= num_resonances_I1)
                {
                    throw InternalError("The residue index must be smaller than the number of used resonances.");
                }

                const complex<double> psi_rk = _psi_r(_M_fp_I1[k](), _G_fp_I1[k]());

                complex<double> result = power_of<2>(psi_rk - 1.0) * (psi_rk + 1.0) / (psi_rk - std::conj(psi_rk));

                for (auto i = 0u; i < num_resonances_I1; ++i)
                {
                    if (i != k)
                    {
                        const complex<double> psi_r = _psi_r(_M_fp_I1[i](), _G_fp_I1[i]());

                        result /= (psi_rk - psi_r) * (psi_rk - std::conj(psi_r));
                    }
                }

                return result;
            }

            inline complex<double>
            _Q_residue(const unsigned & k) const
            {
                const std::size_t num_resonances_I0 = stoi(_n_resonances_I0.value());

                if (k >= num_resonances_I0)
                {
                    throw InternalError("The residue index must be smaller than the number of used resonances.");
                }

                const complex<double> psi_rk = _psi_r(_M_fp_I0[k](), _G_fp_I0[k]());

                complex<double> result = 1.0 / (psi_rk - std::conj(psi_rk));

                for (auto i = 0u; i < num_resonances_I0; ++i)
                {
                    if (i != k)
                    {
                        const complex<double> psi_r = _psi_r(_M_fp_I0[i](), _G_fp_I0[i]());

                        result /= (psi_rk - psi_r) * (psi_rk - std::conj(psi_r));
                    }
                }

                return result;
            }

            // The four following functions are used to compute the expansion parameters a constrained by the threshold and q2 = 0 behaviours
            inline complex<double>
            _dPdpsi(const complex<double> & psi) const
            {
                const std::size_t num_resonances_I1 = stoi(_n_resonances_I1.value());

                complex<double> pole_factor = 1.0;
                complex<double> sum         = 0.0;

                for (auto i = 0u; i < num_resonances_I1; ++i)
                {
                    const complex<double> psi_r = _psi_r(_M_fp_I1[i](), _G_fp_I1[i]());

                    const complex<double> denom = (psi - psi_r) * (psi - std::conj(psi_r));

                    pole_factor /= denom;

                    sum += (2.0 * psi - psi_r - std::conj(psi_r)) / denom;
                }

                const complex<double> numerator = power_of<2>(psi - 1.0) * (psi + 1.0);

                const complex<double> numerator_derivative = (psi - 1.0) * (3.0 * psi + 1.0);

                return pole_factor * (numerator_derivative - numerator * sum);
            }

            inline complex<double>
            _dQdpsi(const complex<double> & psi) const
            {
                const std::size_t     num_resonances_I0 = stoi(_n_resonances_I0.value());
                const complex<double> Q_val             = _Q(psi);
                complex<double>       sum               = complex<double>(0.0, 0.0);

                for (auto i = 0u; i < num_resonances_I0; ++i)
                {
                    const complex<double> psi_r = _psi_r(_M_fp_I0[i](), _G_fp_I0[i]());
                    const complex<double> denom = (psi - psi_r) * (psi - std::conj(psi_r));

                    sum += (2.0 * psi - psi_r - std::conj(psi_r)) / denom;
                }

                return -Q_val * sum;
            }

            inline complex<double>
            _dfdpsi_terms_I1(const unsigned k, const complex<double> & psi) const
            {
                const complex<double> dP_val = _dPdpsi(psi);

                switch (k)
                {
                    case 0: return dP_val;
                    case 1: return dP_val * psi + _P(psi);
                    default:
                        complex<double> psi_km1 = std::pow(psi, k - 1);
                        complex<double> psi_k   = psi_km1 * psi;
                        return dP_val * psi_k + static_cast<double>(k) * _P(psi) * psi_km1;
                }
            }

            inline complex<double>
            _dfdpsi_terms_I0(const unsigned k, const complex<double> & psi) const
            {
                const complex<double> dQ_val = _dQdpsi(psi);

                switch (k)
                {
                    case 0: return dQ_val;
                    case 1: return dQ_val * psi + _Q(psi);
                    default:
                        complex<double> psi_km1 = std::pow(psi, k - 1);
                        complex<double> psi_k   = psi_km1 * psi;
                        return dQ_val * psi_k + static_cast<double>(k) * _Q(psi) * psi_km1;
                }
            }

            // Value and derivatives of P with respect to psi. This function will be used to find scattering lenght parameters.
            struct PDerivatives
            {
                    complex<double> P0;
                    complex<double> P1;
                    complex<double> P2;
                    complex<double> P3;
                    complex<double> P4;
                    complex<double> P5;
            };

            inline PDerivatives
            _P_derivatives(const complex<double> & psi) const
            {
                // This implementation must not be evaluated at psi = +/-1.
                if ((std::abs(psi - 1.0) < 1e-14) || (std::abs(psi + 1.0) < 1e-14))
                {
                    throw InternalError("Some derivatives where evaluated outside their region of validity.");
                }

                const std::size_t num_resonances = stoi(_n_resonances_I1.value());

                const complex<double> P_val = _P(psi);

                complex<double> L = 2.0 / (psi - 1.0) + 1.0 / (psi + 1.0), L1 = -2.0 / power_of<2>(psi - 1.0) - 1.0 / power_of<2>(psi + 1.0),
                                L2 = 4.0 / power_of<3>(psi - 1.0) + 2.0 / power_of<3>(psi + 1.0), L3 = -12.0 / power_of<4>(psi - 1.0) - 6.0 / power_of<4>(psi + 1.0),
                                L4 = 48.0 / power_of<5>(psi - 1.0) + 24.0 / power_of<5>(psi + 1.0);

                for (auto i = 0u; i < num_resonances; ++i)
                {
                    const complex<double> psi_r = _psi_r(_M_fp_I1[i](), _G_fp_I1[i]());

                    const complex<double> inverse_a = 1.0 / (psi - psi_r);

                    const complex<double> inverse_b = 1.0 / (psi - std::conj(psi_r));

                    const complex<double> inverse_a2 = inverse_a * inverse_a;
                    const complex<double> inverse_b2 = inverse_b * inverse_b;

                    const complex<double> inverse_a3 = inverse_a2 * inverse_a;
                    const complex<double> inverse_b3 = inverse_b2 * inverse_b;

                    const complex<double> inverse_a4 = inverse_a3 * inverse_a;
                    const complex<double> inverse_b4 = inverse_b3 * inverse_b;

                    const complex<double> inverse_a5 = inverse_a4 * inverse_a;
                    const complex<double> inverse_b5 = inverse_b4 * inverse_b;

                    L  -= inverse_a + inverse_b;
                    L1 += inverse_a2 + inverse_b2;
                    L2 -= 2.0 * (inverse_a3 + inverse_b3);
                    L3 += 6.0 * (inverse_a4 + inverse_b4);
                    L4 -= 24.0 * (inverse_a5 + inverse_b5);
                }

                PDerivatives result;

                result.P0 = P_val;
                result.P1 = P_val * L;
                result.P2 = P_val * (L * L + L1);
                result.P3 = P_val * (L * L * L + 3.0 * L * L1 + L2);
                result.P4 = P_val * (power_of<4>(L) + 6.0 * L * L * L1 + 3.0 * L1 * L1 + 4.0 * L * L2 + L3);
                result.P5 = P_val * (power_of<5>(L) + 10.0 * power_of<3>(L) * L1 + 15.0 * L * L1 * L1 + 10.0 * L * L * L2 + 10.0 * L1 * L2 + 5.0 * L * L3 + L4);

                return result;
            }

            struct QDerivatives
            {
                    complex<double> Q0;
                    complex<double> Q1;
                    complex<double> Q2;
                    complex<double> Q3;
                    complex<double> Q4;
                    complex<double> Q5;
            };

            inline QDerivatives
            _Q_derivatives(const complex<double> & psi) const
            {
                // This implementation must not be evaluated at psi = +/-1.
                if ((std::abs(psi - 1.0) < 1e-14) || (std::abs(psi + 1.0) < 1e-14))
                {
                    throw InternalError("Some derivatives where evaluated outside their region of validity.");
                }

                const std::size_t num_resonances = stoi(_n_resonances_I0.value());

                const complex<double> Q_val = _Q(psi);

                complex<double> L = 0.0, L1 = 0.0, L2 = 0.0, L3 = 0.0, L4 = 0.0;

                for (auto i = 0u; i < num_resonances; ++i)
                {
                    const complex<double> psi_r = _psi_r(_M_fp_I0[i](), _G_fp_I0[i]());

                    const complex<double> inverse_a = 1.0 / (psi - psi_r);

                    const complex<double> inverse_b = 1.0 / (psi - std::conj(psi_r));

                    const complex<double> inverse_a2 = inverse_a * inverse_a;
                    const complex<double> inverse_b2 = inverse_b * inverse_b;

                    const complex<double> inverse_a3 = inverse_a2 * inverse_a;
                    const complex<double> inverse_b3 = inverse_b2 * inverse_b;

                    const complex<double> inverse_a4 = inverse_a3 * inverse_a;
                    const complex<double> inverse_b4 = inverse_b3 * inverse_b;

                    const complex<double> inverse_a5 = inverse_a4 * inverse_a;
                    const complex<double> inverse_b5 = inverse_b4 * inverse_b;

                    L  -= inverse_a + inverse_b;
                    L1 += inverse_a2 + inverse_b2;
                    L2 -= 2.0 * (inverse_a3 + inverse_b3);
                    L3 += 6.0 * (inverse_a4 + inverse_b4);
                    L4 -= 24.0 * (inverse_a5 + inverse_b5);
                }

                QDerivatives result;

                result.Q0 = Q_val;
                result.Q1 = Q_val * L;
                result.Q2 = Q_val * (L * L + L1);
                result.Q3 = Q_val * (L * L * L + 3.0 * L * L1 + L2);
                result.Q4 = Q_val * (power_of<4>(L) + 6.0 * L * L * L1 + 3.0 * L1 * L1 + 4.0 * L * L2 + L3);
                result.Q5 = Q_val * (power_of<5>(L) + 10.0 * power_of<3>(L) * L1 + 15.0 * L * L1 * L1 + 10.0 * L * L * L2 + 10.0 * L1 * L2 + 5.0 * L * L3 + L4);

                return result;
            }

            /*
             * Stores a function value and its first five ordinary derivatives with respect to psi:
             *
             *     Fk = d^k F / d psi^k.
             *
             * These are derivatives, not Taylor coefficients; the corresponding Taylor term is Fk / k!.
             */
            struct FormFactorDerivatives
            {
                    complex<double> F0;
                    complex<double> F1;
                    complex<double> F2;
                    complex<double> F3;
                    complex<double> F4;
                    complex<double> F5;
            };

            /*
             * Computes the value and first five derivatives of the series
             *
             *     S(psi) = sum_k a[k] psi^k.
             *
             * This is used for both the I = 1 and I = 0 expansions.
             */
            inline FormFactorDerivatives
            _series_derivatives(const complex<double> & psi, const std::array<double, 13u> & a) const
            {
                std::array<complex<double>, 13u> psi_power;

                psi_power[0] = 1.0;

                for (auto k = 1u; k < psi_power.size(); ++k)
                {
                    psi_power[k] = psi_power[k - 1u] * psi;
                }

                FormFactorDerivatives result{};

                for (auto k = 0u; k < a.size(); ++k)
                {
                    result.F0 += a[k] * psi_power[k];

                    if (k >= 1u)
                    {
                        result.F1 += a[k] * static_cast<double>(k) * psi_power[k - 1u];
                    }

                    if (k >= 2u)
                    {
                        result.F2 += a[k] * static_cast<double>(k * (k - 1u)) * psi_power[k - 2u];
                    }

                    if (k >= 3u)
                    {
                        result.F3 += a[k] * static_cast<double>(k * (k - 1u) * (k - 2u)) * psi_power[k - 3u];
                    }

                    if (k >= 4u)
                    {
                        result.F4 += a[k] * static_cast<double>(k * (k - 1u) * (k - 2u) * (k - 3u)) * psi_power[k - 4u];
                    }

                    if (k >= 5u)
                    {
                        result.F5 += a[k] * static_cast<double>(k * (k - 1u) * (k - 2u) * (k - 3u) * (k - 4u)) * psi_power[k - 5u];
                    }
                }

                return result;
            }

            /*
             * Computes the value and first five derivatives of the product
             * A(psi) B(psi) using the Leibniz rule.
             */
            inline FormFactorDerivatives
            _product_derivatives(const FormFactorDerivatives & A, const FormFactorDerivatives & B) const
            {
                FormFactorDerivatives result;

                result.F0 = A.F0 * B.F0;
                result.F1 = A.F1 * B.F0 + A.F0 * B.F1;
                result.F2 = A.F2 * B.F0 + 2.0 * A.F1 * B.F1 + A.F0 * B.F2;
                result.F3 = A.F3 * B.F0 + 3.0 * A.F2 * B.F1 + 3.0 * A.F1 * B.F2 + A.F0 * B.F3;
                result.F4 = A.F4 * B.F0 + 4.0 * A.F3 * B.F1 + 6.0 * A.F2 * B.F2 + 4.0 * A.F1 * B.F3 + A.F0 * B.F4;
                result.F5 = A.F5 * B.F0 + 5.0 * A.F4 * B.F1 + 10.0 * A.F3 * B.F2 + 10.0 * A.F2 * B.F3 + 5.0 * A.F1 * B.F4 + A.F0 * B.F5;

                return result;
            }

            /*
             * Computes the derivatives of
             *
             *     f_I1(psi) = P(psi) * S_I1(psi).
             */

            inline FormFactorDerivatives
            _f_I1_derivatives(const complex<double> & psi, const std::array<double, 13u> & coefficients) const
            {
                const auto                  P = this->_P_derivatives(psi);
                const FormFactorDerivatives P_derivatives{ P.P0, P.P1, P.P2, P.P3, P.P4, P.P5 };
                const FormFactorDerivatives S_derivatives = this->_series_derivatives(psi, coefficients);

                return this->_product_derivatives(P_derivatives, S_derivatives);
            }

            /*
             * Computes the derivatives of
             *
             *     delta_I0(psi) = Q(psi) * S_I0(psi).
             */

            inline FormFactorDerivatives
            _f_I0_derivatives(const complex<double> & psi, const std::array<double, 13u> & coefficients) const
            {
                const auto                  Q             = _Q_derivatives(psi);
                const FormFactorDerivatives Q_derivatives = { Q.Q0, Q.Q1, Q.Q2, Q.Q3, Q.Q4, Q.Q5 };
                const FormFactorDerivatives S_derivatives = this->_series_derivatives(psi, coefficients);

                return this->_product_derivatives(Q_derivatives, S_derivatives);
            }

            /*
             * Computes the derivatives of the complete form factor
             *
             *     f_p(psi) = f_I1(psi) * (switch_I1 + switch_I0 * delta_I0(psi)).
             *
             * The I = 0 contribution is included only when enabled. These
             * derivatives are used for the charge radius and for the derivative
             * ratios entering the threshold parameters a11 and b11.
             */

            inline FormFactorDerivatives
            _f_p_derivatives(const complex<double> & psi) const
            {
                // Prepare I = 1 expansion coefficients
                std::array<double, 13u> a_I1{};

                const auto constrained_a_I1 = constrained_a_fp_I1();
                std::copy(constrained_a_I1.cbegin(), constrained_a_I1.cend(), a_I1.begin());

                for (auto i = 0u; i < _a_fp_I1.size(); ++i)
                {
                    a_I1[i + 4u] = _a_fp_I1[i]();
                }

                const FormFactorDerivatives f_I1 = _f_I1_derivatives(psi, a_I1);

                /*
                 * Derivatives of the second factor:
                 *
                 *     switch_I1 + switch_I0 * delta_I0(psi)
                 *
                 * If I = 0 is disabled, this is simply the constant
                 * switch_I1, whose higher derivatives vanish.
                 */
                FormFactorDerivatives multiplier{ static_cast<double>(_switch_I[1]), 0.0, 0.0, 0.0, 0.0, 0.0 };

                if (_switch_I[0])
                {
                    // Prepare I = 0 expansion coefficients
                    std::array<double, 13u> a_I0{};

                    const auto constrained_a_I0 = constrained_a_fp_I0();

                    std::copy(constrained_a_I0.cbegin(), constrained_a_I0.cend(), a_I0.begin());

                    for (auto i = 0u; i < _a_fp_I0.size(); ++i)
                    {
                        a_I0[i + 4u] = _a_fp_I0[i]();
                    }

                    const FormFactorDerivatives f_I0 = _f_I0_derivatives(psi, a_I0);

                    multiplier.F0 += f_I0.F0;
                    multiplier.F1  = f_I0.F1;
                    multiplier.F2  = f_I0.F2;
                    multiplier.F3  = f_I0.F3;
                    multiplier.F4  = f_I0.F4;
                    multiplier.F5  = f_I0.F5;
                }

                return _product_derivatives(f_I1, multiplier);
            }


        public:
            BHKMNR2026FormFactors(const Parameters & p, const Options & o);
            ~BHKMNR2026FormFactors();

            static FormFactors<VacuumToPP> * make(const Parameters & p, const Options & o);

            /* auxiliary functions */
            std::array<double, 4u> constrained_a_fp_I1() const;

            double
            a_0() const
            {
                return this->constrained_a_fp_I1()[0];
            }

            double
            a_1() const
            {
                return this->constrained_a_fp_I1()[1];
            }

            double
            a_2() const
            {
                return this->constrained_a_fp_I1()[2];
            }

            double
            a_3() const
            {
                return this->constrained_a_fp_I1()[3];
            }

            std::array<double, 4u> constrained_a_fp_I0() const;

            complex<double> psi(const complex<double> & s) const;
            complex<double> psi21(const complex<double> & s) const;
            complex<double> P(const complex<double> & psi) const;
            complex<double> Q(const complex<double> & psi) const;
            complex<double> dPdpsi(const complex<double> & psi) const;
            complex<double> dQdpsi(const complex<double> & psi) const;
            complex<double> dfdpsi_terms_I1(const unsigned k, const complex<double> & psi) const;
            complex<double> dfdpsi_terms_I0(const unsigned k, const complex<double> & psi) const;
            complex<double> series(const complex<double> & psi, const std::array<double, 13> & a) const;
            complex<double> f_p_of_psi(const complex<double> & psi) const;

            double
            abs2_f_p_of_psi(const double & re_psi, const double & im_psi) const
            {
                return std::norm(f_p_of_psi(complex<double>(re_psi, im_psi)));
            }

            double
            arg_f_p_of_psi(const double & re_psi, const double & im_psi) const
            {
                return std::arg(f_p_of_psi(complex<double>(re_psi, im_psi)));
            }

            /* form factors on the real axis */
            virtual complex<double> f_p(const double & s) const override;
            virtual complex<double> f_0(const double & s) const override;
            virtual complex<double> f_t(const double & s) const override;

            /* form factor in the complex s plane */
            virtual complex<double> f_p(const complex<double> & s) const override;
            virtual complex<double> f_0(const complex<double> & s) const override;
            virtual complex<double> f_t(const complex<double> & s) const override;

            /* form factors on the 21 Rieman sheet */
            complex<double> f_p_21(const double & s) const;
            complex<double> f_p_21(const complex<double> & s) const;

            /* Isospin 1, P-wave partial wave */
            complex<double> partial_wave(const double & s) const;
            complex<double> partial_wave(const complex<double> & s) const;

            double
            re_partial_wave(const double & s) const
            {
                return std::real(partial_wave(s));
            }

            double
            im_partial_wave(const double & s) const
            {
                return std::imag(partial_wave(s));
            }

            std::array<complex<double>, 4u> scattering_length_parameters() const;

            double
            d2fdpsi2_over_f() const
            {
                return std::real(scattering_length_parameters()[0]);
            }

            double
            d3fdpsi3_over_f() const
            {
                return std::real(scattering_length_parameters()[1]);
            }

            double
            d4fdpsi4_over_f() const
            {
                return std::real(scattering_length_parameters()[2]);
            }

            double
            d5fdpsi5_over_f() const
            {
                return std::real(scattering_length_parameters()[3]);
            }

            complex<double> dfdpsi_11(const complex<double> & s) const;

            double
            dfdpsi_11_at_0() const
            {
                return std::real(dfdpsi_11(0.0));
            }

            double dispersive_integrand(const double & s) const;
            double saturation() const;

            // residue functions
            complex<double> residue_I1(const unsigned & k) const;
            double          re_residue_rho() const;
            double          im_residue_rho() const;


            // Test of roots on the first RS
            // Returns the sum of the inverse of the modulus of the roots in the first Riemann sheet
            double root_penalty() const;

            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification>           option_specifications;

            static const std::set<ReferenceName> references;
    };

    extern template class BHKMNR2026FormFactors<VacuumToPiPi>;
} // namespace eos

#endif
