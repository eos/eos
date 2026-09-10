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

#include <eos/form-factors/parametric-bhkmnr2026.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_poly.h>

#include <functional>
#include <numeric>

namespace eos
{
    using namespace std::literals::string_literals;

    BHKMNR2026FormFactors<VacuumToPiPi>::BHKMNR2026FormFactors(const Parameters & p, const Options & o) :
        _a_fp_I1{
            { UsedParameter(p[_par_name("+", "1", "4")], *this),
             UsedParameter(p[_par_name("+", "1", "5")], *this),
             UsedParameter(p[_par_name("+", "1", "6")], *this),
             UsedParameter(p[_par_name("+", "1", "7")], *this),
             UsedParameter(p[_par_name("+", "1", "8")], *this),
             UsedParameter(p[_par_name("+", "1", "9")], *this),
             UsedParameter(p[_par_name("+", "1", "10")], *this),
             UsedParameter(p[_par_name("+", "1", "11")], *this),
             UsedParameter(p[_par_name("+", "1", "12")], *this) }
    },
        _M_fp_I1{ { UsedParameter(p["0->pipi::M_(+,1,0)@BHKMNR2026"], *this),
                    UsedParameter(p["0->pipi::M_(+,1,1)@BHKMNR2026"], *this),
                    UsedParameter(p["0->pipi::M_(+,1,2)@BHKMNR2026"], *this) } },
        _G_fp_I1{ { UsedParameter(p["0->pipi::Gamma_(+,1,0)@BHKMNR2026"], *this),
                    UsedParameter(p["0->pipi::Gamma_(+,1,1)@BHKMNR2026"], *this),
                    UsedParameter(p["0->pipi::Gamma_(+,1,2)@BHKMNR2026"], *this) } },
        _n_resonances_I1(o, option_specifications, "n-resonances-I1"_ok),
        _a_fp_I0{ { UsedParameter(p[_par_name("+", "0", "4")], *this),
                    UsedParameter(p[_par_name("+", "0", "5")], *this),
                    UsedParameter(p[_par_name("+", "0", "6")], *this),
                    UsedParameter(p[_par_name("+", "0", "7")], *this),
                    UsedParameter(p[_par_name("+", "0", "8")], *this) } },
        _M_fp_I0{ { UsedParameter(p["0->pipi::M_(+,0,0)@BHKMNR2026"], *this), UsedParameter(p["0->pipi::M_(+,0,1)@BHKMNR2026"], *this) } },
        _G_fp_I0{ { UsedParameter(p["0->pipi::Gamma_(+,0,0)@BHKMNR2026"], *this), UsedParameter(p["0->pipi::Gamma_(+,0,1)@BHKMNR2026"], *this) } },
        _n_resonances_I0(o, option_specifications, "n-resonances-I0"_ok), _m_pi(p["mass::pi^+"], *this), _s_0(p["0->pipi::s_0@BHKMNR2026"], *this),
        _s_in(p["0->pipi::s_in@BHKMNR2026"], *this), _hbar(p["QM::hbar"], *this), _opt_I(o, option_specifications, "I"_ok), _M(gsl_matrix_alloc(4, 4)),
        _inv_M(gsl_matrix_alloc(4, 4)), _L(gsl_vector_alloc(4)), _perm(gsl_permutation_calloc(4)), _constrained_coefficents(gsl_vector_alloc(4)),
        _poly_workspace(gsl_poly_complex_workspace_alloc(13))
    {
        // Perform pointer checks
        if (_M == nullptr)
        {
            throw std::bad_alloc();
        }
        if (_inv_M == nullptr)
        {
            throw std::bad_alloc();
        }
        if (_L == nullptr)
        {
            throw std::bad_alloc();
        }
        if (_perm == nullptr)
        {
            throw std::bad_alloc();
        }
        if (_constrained_coefficents == nullptr)
        {
            throw std::bad_alloc();
        }
        if (_poly_workspace == nullptr)
        {
            throw std::bad_alloc();
        }

        _switch_I[0] = (_opt_I.value() && Isospin::zero);
        _switch_I[1] = (_opt_I.value() && Isospin::one);
    }

    BHKMNR2026FormFactors<VacuumToPiPi>::~BHKMNR2026FormFactors()
    {
        if (_perm)
        {
            gsl_permutation_free(_perm);
        }
        _perm = nullptr;
        if (_constrained_coefficents)
        {
            gsl_vector_free(_constrained_coefficents);
        }
        _constrained_coefficents = nullptr;
        if (_L)
        {
            gsl_vector_free(_L);
        }
        _L = nullptr;
        if (_inv_M)
        {
            gsl_matrix_free(_inv_M);
        }
        _inv_M = nullptr;
        if (_M)
        {
            gsl_matrix_free(_M);
        }
        _M = nullptr;
        if (_poly_workspace)
        {
            gsl_poly_complex_workspace_free(_poly_workspace);
        }
        _poly_workspace = nullptr;
    }

    FormFactors<VacuumToPP> *
    BHKMNR2026FormFactors<VacuumToPiPi>::make(const Parameters & p, const Options & o)
    {
        return new BHKMNR2026FormFactors<VacuumToPiPi>(p, o);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::psi(const complex<double> & s) const
    {
        return this->_s_to_psi_11(s);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::psi21(const complex<double> & s) const
    {
        return this->_s_to_psi_21(s);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::P(const complex<double> & psi) const
    {
        return this->_P(psi);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::dPdpsi(const complex<double> & psi) const
    {
        return this->_dPdpsi(psi);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::Q(const complex<double> & psi) const
    {
        return this->_Q(psi);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::dQdpsi(const complex<double> & psi) const
    {
        return this->_dQdpsi(psi);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::dfdpsi_terms_I1(const unsigned k, const complex<double> & psi) const
    {
        return this->_dfdpsi_terms_I1(k, psi);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::dfdpsi_terms_I0(const unsigned k, const complex<double> & psi) const
    {
        return this->_dfdpsi_terms_I0(k, psi);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::series(const complex<double> & psi, const std::array<double, 13> & a) const
    {
        complex<double> series    = 0.0;
        complex<double> psi_power = 1.0;

        for (auto k = 0u; k < a.size(); k++)
        {
            series    += a[k] * psi_power;
            psi_power *= psi;
        }

        return series;
    }

    std::array<double, 4u>
    BHKMNR2026FormFactors<VacuumToPiPi>::constrained_a_fp_I1() const
    {
        const complex<double> psi_p  = _s_to_psi_11(_s_p());
        const complex<double> psi_in = _s_to_psi_11(_s_in());
        const complex<double> psi_0  = _s_to_psi_11(complex<double>(0.0, 0.0));
        const complex<double> P0     = _P(psi_0);

        // Solve the system M.{a0, a1, a2, a3} = L
        // Fill M
        complex<double> psi0_pow = complex<double>(1.0, 0.0);

        for (auto k = 0u; k < 4; k++)
        {
            const complex<double> val_p  = _dfdpsi_terms_I1(k, psi_p);
            const complex<double> val_in = _dfdpsi_terms_I1(k, psi_in);

            gsl_matrix_set(_M, 0, k, val_p.real());
            gsl_matrix_set(_M, 1, k, val_in.real());
            gsl_matrix_set(_M, 2, k, val_in.imag());
            gsl_matrix_set(_M, 3, k, (P0 * psi0_pow).real());

            psi0_pow *= psi_0;
        }

        // Fill L
        const unsigned n = 3 + _a_fp_I1.size();

        complex<double> sum_p  = complex<double>(0.0, 0.0);
        complex<double> sum_in = complex<double>(0.0, 0.0);
        complex<double> sum_0  = complex<double>(0.0, 0.0);

        for (auto k = 4u; k <= n; k++)
        {
            const double ak = _a_fp_I1[k - 4]();

            sum_p  += ak * _dfdpsi_terms_I1(k, psi_p);
            sum_in += ak * _dfdpsi_terms_I1(k, psi_in);
            sum_0  += ak * psi0_pow;

            psi0_pow *= psi_0;
        }

        gsl_vector_set(_L, 0, -sum_p.real());
        gsl_vector_set(_L, 1, -sum_in.real());
        gsl_vector_set(_L, 2, -sum_in.imag());

        const complex<double> entry = 1.0 - P0 * sum_0; // f_p_I1(s = 0) = 1.0
        gsl_vector_set(_L, 3, entry.real());

        // Invert M and solve for the constrained coefficients
        int signum = 0;
        gsl_permutation_init(_perm);
        gsl_linalg_LU_decomp(_M, _perm, &signum);
        gsl_linalg_LU_invert(_M, _perm, _inv_M);

        gsl_blas_dgemv(CblasNoTrans, 1.0, _inv_M, _L, 0.0, _constrained_coefficents);

        std::array<double, 4u> result;
        for (auto i = 0u; i < 4; ++i)
        {
            result[i] = gsl_vector_get(_constrained_coefficents, i);
        }
        return result;
    }

    std::array<double, 4u>
    BHKMNR2026FormFactors<VacuumToPiPi>::constrained_a_fp_I0() const
    {
        const complex<double> psi_p  = _s_to_psi_11(_s_p());
        const complex<double> psi_in = _s_to_psi_11(_s_in());
        const complex<double> psi_0  = _s_to_psi_11(complex<double>(0.0, 0.0));
        const complex<double> Q0     = _Q(psi_0);

        // Solve the system M.{a0, a1, a2, a3} = L
        // Fill M
        complex<double> psi0_pow = complex<double>(1.0, 0.0);

        for (auto k = 0u; k < 4; k++)
        {
            const complex<double> val_p  = _dfdpsi_terms_I0(k, psi_p);
            const complex<double> val_in = _dfdpsi_terms_I0(k, psi_in);

            gsl_matrix_set(_M, 0, k, val_p.real());
            gsl_matrix_set(_M, 1, k, val_in.real());
            gsl_matrix_set(_M, 2, k, val_in.imag());
            gsl_matrix_set(_M, 3, k, (Q0 * psi0_pow).real());

            psi0_pow *= psi_0;
        }

        // Fill L
        const unsigned n = 3 + _a_fp_I0.size();

        complex<double> sum_p  = complex<double>(0.0, 0.0);
        complex<double> sum_in = complex<double>(0.0, 0.0);
        complex<double> sum_0  = complex<double>(0.0, 0.0);

        for (auto k = 4u; k <= n; k++)
        {
            const double ak = _a_fp_I0[k - 4]();

            sum_p  += ak * _dfdpsi_terms_I0(k, psi_p);
            sum_in += ak * _dfdpsi_terms_I0(k, psi_in);
            sum_0  += ak * psi0_pow;

            psi0_pow *= psi_0;
        }

        gsl_vector_set(_L, 0, -sum_p.real());
        gsl_vector_set(_L, 1, -sum_in.real());
        gsl_vector_set(_L, 2, -sum_in.imag());

        const complex<double> entry = -Q0 * sum_0; // delta_I0(s = 0) = 0.0
        gsl_vector_set(_L, 3, entry.real());

        // Invert M and solve for the constrained coefficients
        int signum = 0;
        gsl_permutation_init(_perm);
        gsl_linalg_LU_decomp(_M, _perm, &signum);
        gsl_linalg_LU_invert(_M, _perm, _inv_M);

        gsl_blas_dgemv(CblasNoTrans, 1.0, _inv_M, _L, 0.0, _constrained_coefficents);

        std::array<double, 4u> result;
        for (auto i = 0u; i < 4; ++i)
        {
            result[i] = gsl_vector_get(_constrained_coefficents, i);
        }
        return result;
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_p_of_psi(const complex<double> & psi) const
    {
        // prepare I=1 expansion coefficients
        std::array<double, 13> a_I1;
        const auto             constrained_a_I1 = this->constrained_a_fp_I1();
        std::copy(constrained_a_I1.cbegin(), constrained_a_I1.cend(), a_I1.begin()); // copy constrained coefficients
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a_I1.begin() + 4);             // copy unconstrained coefficients

        const complex<double> f_I1 = this->_P(psi) * this->series(psi, a_I1);

        complex<double> delta_I0 = 0.0;

        if (_switch_I[0])
        {
            // prepare I=1 expansion coefficients
            std::array<double, 13> a_I0;
            a_I0.fill(0.0);
            const auto constrained_a_I0 = this->constrained_a_fp_I0();
            std::copy(constrained_a_I0.cbegin(), constrained_a_I0.cend(), a_I0.begin()); // copy constrained coefficients
            std::copy(_a_fp_I0.cbegin(), _a_fp_I0.cend(), a_I0.begin() + 4);             // copy unconstrained coefficients

            delta_I0 = this->_Q(psi) * this->series(psi, a_I0);
        }

        return f_I1 * (static_cast<double>(_switch_I[1]) + delta_I0);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_p(const complex<double> & s) const
    {
        return f_p_of_psi(this->_s_to_psi_11(s));
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_p(const double & s) const
    {
        static const double eps = 1.0e-14;
        return f_p(complex<double>(s, eps));
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_p_21(const complex<double> & s) const
    {
        return f_p_of_psi(this->_s_to_psi_21(s));
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_p_21(const double & s) const
    {
        static const double eps = 1.0e-14;
        return f_p_21(complex<double>(s, eps));
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::partial_wave(const complex<double> & s) const
    {
        const complex<double> s_p = this->_s_p();
        static const double   eps = 1.0e-14;

        if (std::abs(s - s_p) < eps)
        {
            return complex<double>(0.0, 0.0);
        }

        // epsilon prescription
        complex<double> s_evaluation = s;
        if (s.real() == 0.0 && s.imag() == 0.0)
        {
            s_evaluation = complex<double>(0.0, 1.0e-14);
        }

        const complex<double> f_p_11 = this->f_p(s_evaluation);

        // Evaluate directly to avoid the near-zero cutoff in f_p_21().
        const complex<double> f_p_21 = this->f_p_of_psi(this->_s_to_psi_21(s_evaluation));

        return std::sqrt(s_evaluation) / (2.0 * std::sqrt(s_p - s_evaluation)) * (1.0 - f_p_11 / f_p_21);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::partial_wave(const double & s) const
    {
        static const double eps = 1.0e-14;
        return partial_wave(complex<double>(s, eps));
    }

    std::array<complex<double>, 4u>
    BHKMNR2026FormFactors<VacuumToPiPi>::scattering_length_parameters() const
    {
        const complex<double> psi_p = this->_s_to_psi_11(this->_s_p());

        // prepare I=1 expansion coefficients
        std::array<double, 13> a_I1;
        const auto             constrained_a_I1 = this->constrained_a_fp_I1();
        std::copy(constrained_a_I1.cbegin(), constrained_a_I1.cend(), a_I1.begin()); // copy constrained coefficients
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a_I1.begin() + 4);             // copy unconstrained coefficients

        const FormFactorDerivatives f_I1       = this->_f_I1_derivatives(psi_p, a_I1);
        FormFactorDerivatives       multiplier = { static_cast<double>(_switch_I[1]), 0.0, 0.0, 0.0, 0.0, 0.0 };

        if (_switch_I[0])
        {
            // prepare I=1 expansion coefficients
            std::array<double, 13> a_I0;
            a_I0.fill(0.0);
            const auto constrained_a_I0 = this->constrained_a_fp_I0();
            std::copy(constrained_a_I0.cbegin(), constrained_a_I0.cend(), a_I0.begin()); // copy constrained coefficients
            std::copy(_a_fp_I0.cbegin(), _a_fp_I0.cend(), a_I0.begin() + 4);             // copy unconstrained coefficients

            const FormFactorDerivatives f_I0 = this->_f_I0_derivatives(psi_p, a_I0);

            multiplier.F0 += f_I0.F0;
            multiplier.F1  = f_I0.F1;
            multiplier.F2  = f_I0.F2;
            multiplier.F3  = f_I0.F3;
            multiplier.F4  = f_I0.F4;
            multiplier.F5  = f_I0.F5;
        }

        const FormFactorDerivatives f_full = this->_product_derivatives(f_I1, multiplier);

        return { f_full.F2 / f_full.F0, f_full.F3 / f_full.F0, f_full.F4 / f_full.F0, f_full.F5 / f_full.F0 };
    }

    // This function is used to compute charge radius of pion
    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::dfdpsi_11(const complex<double> & s) const
    {
        const complex<double> psi = _s_to_psi_11(s);

        return _f_p_derivatives(psi).F1;
    }

    double
    BHKMNR2026FormFactors<VacuumToPiPi>::dispersive_integrand(const double & x) const
    {
        // change of variable s = s_p + x / (1 - x), then integral over s from s_p to infinity will become integral over x from 0 to 1
        const double Q2    = 1.0;
        const double chi   = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A
        const double s_p   = real(this->_s_p());
        const double denom = 48.0 * power_of<2>(M_PI) * chi;

        complex<double> f_p = this->f_p(s_p + x / (1.0 - x));

        return std::pow(x, 1.5) * std::norm(f_p) / std::sqrt(s_p * (1 - x) + x) / power_of<3>((s_p + Q2) * (1 - x) + x) / denom;
    }

    double
    BHKMNR2026FormFactors<VacuumToPiPi>::saturation() const
    {
        std::function<double(const double &)> f = [this](const double & x) -> double { return this->dispersive_integrand(x); };
        return integrate<1, 1>(f, 0, 1, cubature::Config().epsrel(1.0e-5));
    }

    // This function is used to compute the rho-pipi and rho-gamma couplings
    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::residue_I1(const unsigned & k) const
    {
        // prepare expansion coefficients
        std::array<double, 13> a;
        const auto             constrained_a = this->constrained_a_fp_I1();
        std::copy(constrained_a.cbegin(), constrained_a.cend(), a.begin()); // copy constrained coefficients
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a.begin() + 4);       // copy unconstrained coefficients

        const complex<double> psi_r    = this->_psi_r(this->_M_fp_I1[k](), this->_G_fp_I1[k]());
        const complex<double> series_r = series(psi_r, a);
        const complex<double> P_r      = this->_P_residue(k);

        return P_r * series_r;
    }

    double
    BHKMNR2026FormFactors<VacuumToPiPi>::re_residue_rho() const
    {
        return std::real(this->residue_I1(0u));
    }

    double
    BHKMNR2026FormFactors<VacuumToPiPi>::im_residue_rho() const
    {
        return std::imag(this->residue_I1(0u));
    }

    double
    BHKMNR2026FormFactors<VacuumToPiPi>::root_penalty() const
    {
        // prepare expansion coefficients
        std::array<double, 13> reversed_a;
        const auto             constrained_a = this->constrained_a_fp_I1();

        if (abs(constrained_a[0]) < 1e-14)
        {
            throw InternalError("Not implemented!");
            return 0.0;
        }

        std::copy(constrained_a.cbegin(), constrained_a.cend(), reversed_a.rbegin()); // reverse copy constrained coefficients
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), reversed_a.rbegin() + 4);       // reverse copy unconstrained coefficients

        // Fill the root array with values on the second Riemann sheet
        std::array<complex<double>, 12> roots;
        roots.fill(-1.0);
        // Find the roots of the reciprocal adjoint polynomial since a_0 != 0 is less likely
        double adjoint_roots[24];
        gsl_poly_complex_solve(reversed_a.data(), 13, _poly_workspace, adjoint_roots);

        // Convert the roots back to our polynomial roots (no need for complex conjugation since the roots are conjugated)
        for (size_t i = 0; i < 12; ++i)
        {
            complex<double> adjoint_root = complex<double>(adjoint_roots[2 * i], adjoint_roots[2 * i + 1]);
            if (abs(adjoint_root) > 1e-12)
            {
                roots[i] = 1.0 / adjoint_root;
            }
        }

        double penalty = 0.0;
        for (auto r : roots)
        {
            if (abs(r) < 1.0)
            {
                const complex<double> phi_root = _chi_inverse(r, _s_to_phi_21(_s_m(), _s_in()), _s_to_phi_11(_s_0(), _s_in()));
                if (std::real(phi_root) > 0.0 && std::abs(phi_root) < 1.0) // Root on RS11
                {
                    penalty += 1.0 / std::abs(r) - 1.0;
                }
                else if (std::real(phi_root) < 0.0 && std::abs(phi_root) < 1.0 && std::abs(std::imag(r)) < 1e-12) // Root on the real axis of RS21
                {
                    penalty += 1.0 / std::abs(r) - 1.0;
                }
            }
        }

        return 1.0 + penalty;
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_0(const double & /*s*/) const
    {
        return 0.0; // vanishes in our approximation
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_0(const complex<double> & /*s*/) const
    {
        return 0.0; // vanishes in our approximation
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_t(const double & /*s*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_t(const complex<double> & /*s*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    const std::vector<OptionSpecification> BHKMNR2026FormFactors<VacuumToPiPi>::option_specifications{
        { "n-resonances-I1"_ok, { "1"_ov, "2"_ov, "3"_ov },   "1"_ov },
        { "n-resonances-I0"_ok,         { "1"_ov, "2"_ov },   "1"_ov },
        {               "I"_ok,               { "0|1"_ov }, "0|1"_ov }
    };

    std::vector<OptionSpecification>::const_iterator
    BHKMNR2026FormFactors<VacuumToPiPi>::begin_options()
    {
        return option_specifications.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BHKMNR2026FormFactors<VacuumToPiPi>::end_options()
    {
        return option_specifications.cend();
    }

    const std::set<ReferenceName> BHKMNR2026FormFactors<VacuumToPiPi>::references{ "BHKMR:2025A"_rn };

} // namespace eos
