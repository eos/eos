/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2023 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2018, 2021 Christoph Bobeth
 * Copyright (c) 2022 Philip Lüghausen
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

#include <eos/models/top-loops.hh>
#include <eos/models/standard-model.hh>
#include <eos/maths/matrix.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/log.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/rge-impl.hh>
#include <eos/utils/stringify.hh>

#include <array>
#include <cmath>

#include <gsl/gsl_sf_clausen.h>
#include <gsl/gsl_sf_dilog.h>

namespace eos
{
    using std::sqrt;

    SMComponent<components::CKM>::SMComponent(const Parameters & p, ParameterUser & u) :
        _A__ckm(p["CKM::A"], u),
        _lambda__ckm(p["CKM::lambda"], u),
        _rhobar__ckm(p["CKM::rhobar"], u),
        _etabar__ckm(p["CKM::etabar"], u)
    {
    }

namespace implementation
{
    // return rho + i eta, cf. [CKMfitter04], Eq. (17), p. 12
    complex<double> rho_eta(const double & A, const double & lambda, const double & rhobar, const double & etabar)
    {
        double A2 = power_of<2>(A), lambda2 = power_of<2>(lambda), lambda4 = power_of<2>(lambda2);

        complex<double> result = complex<double>(rhobar, etabar) * std::sqrt(1.0 - A2 * lambda4)
            / std::sqrt(1.0 - lambda2) / (1.0 - A2 * lambda4 * complex<double>(rhobar, etabar));

        return result;
    }
}

    /*
     * For the parametrisation of all CKM matrix elements, cf. [CKMfitter04], Footnote 4, p. 10
     */
    complex<double>
    SMComponent<components::CKM>::ckm_cd() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A__ckm, _lambda__ckm, _rhobar__ckm, _etabar__ckm);
        double A2 = power_of<2>(_A__ckm()), lambda4 = power_of<4>(_lambda__ckm()), lambda6 = power_of<6>(_lambda__ckm());

        complex<double> result = _lambda__ckm() * (1.0 - A2 * lambda4 * (1.0 - 2.0 * rho_eta) / 2.0 + A2 * lambda6 * rho_eta / 2.0);

        return result;
    }

    complex<double>
    SMComponent<components::CKM>::ckm_cs() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A__ckm, _lambda__ckm, _rhobar__ckm, _etabar__ckm);
        double A2 = power_of<2>(_A__ckm()), A4 = power_of<2>(A2);
        double lambda2 = power_of<2>(_lambda__ckm()), lambda4 = power_of<2>(lambda2), lambda6 = lambda4 * lambda2, lambda8 = lambda4 * lambda4;;

        complex<double> result = 1.0 - lambda2 / 2.0 - lambda4 * (1.0 + 4.0 * A2) / 8.0
            - lambda6 * (1.0 - 4.0 * A2 + 16.0 * A2 * rho_eta) / 16.0 - lambda8 * (5.0 - 8.0 * A2 + 16.0 * A4) / 128.0;

        return result;
    }

    complex<double>
    SMComponent<components::CKM>::ckm_cb() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A__ckm, _lambda__ckm, _rhobar__ckm, _etabar__ckm);
        double A2 = power_of<2>(_A__ckm()), lambda2 = power_of<2>(_lambda__ckm()), lambda6 = power_of<3>(lambda2);

        double result = _A__ckm * lambda2 * (1.0 - 0.5 * A2 * lambda6 * std::norm(rho_eta));

        return complex<double>(result, 0.0);
    }

    complex<double>
    SMComponent<components::CKM>::ckm_ud() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A__ckm, _lambda__ckm, _rhobar__ckm, _etabar__ckm);
        double A2 = power_of<2>(_A__ckm()), lambda2 = power_of<2>(_lambda__ckm()), lambda4 = lambda2 * lambda2,
               lambda6 = lambda2 * lambda4, lambda8 = lambda4 * lambda4;

        double result = 1.0 - lambda2 / 2.0 - lambda4 / 8.0 - lambda6 * (1.0 + 8.0 * A2 * std::norm(rho_eta)) / 16.0
            - lambda8 * (5.0 - 32.0 * A2 * std::norm(rho_eta)) / 128.0;

        return complex<double>(result, 0.0);
    }

    complex<double>
    SMComponent<components::CKM>::ckm_us() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A__ckm, _lambda__ckm, _rhobar__ckm, _etabar__ckm);
        double A2 = power_of<2>(_A__ckm()), lambda6 = power_of<6>(_lambda__ckm());

        double result = _lambda__ckm * (1.0 - 0.5 * A2 * lambda6 * std::norm(rho_eta));

        return complex<double>(result, 0.0);
    }

    complex<double>
    SMComponent<components::CKM>::ckm_ub() const
    {
        complex<double> rho_eta_conj = std::conj(implementation::rho_eta(_A__ckm, _lambda__ckm, _rhobar__ckm, _etabar__ckm));

        complex<double> result = _A__ckm * power_of<3>(_lambda__ckm()) * rho_eta_conj;

        return result;
    }

    complex<double>
    SMComponent<components::CKM>::ckm_td() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A__ckm, _lambda__ckm, _rhobar__ckm, _etabar__ckm);
        double A2 = power_of<2>(_A__ckm()), lambda2 = power_of<2>(_lambda__ckm()), lambda3 = _lambda__ckm() * lambda2, lambda4 = lambda2 * lambda2;

        complex<double> result = _A__ckm * lambda3 *
            ((1.0 - rho_eta) + lambda2 * rho_eta / 2.0 + lambda4 * (1.0 + 4.0 * A2) * rho_eta / 8.0);

        return result;
    }

    complex<double>
    SMComponent<components::CKM>::ckm_ts() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A__ckm, _lambda__ckm, _rhobar__ckm, _etabar__ckm);
        double A2 = power_of<2>(_A__ckm()), lambda2 = power_of<2>(_lambda__ckm()), lambda4 = lambda2 * lambda2, lambda6 = lambda2 * lambda4;

        complex<double> result = -1.0 * _A__ckm * lambda2 *
            (1.0 - lambda2 * (1.0 - 2.0 * rho_eta) / 2.0 - lambda4 / 8.0 - lambda6 * (1.0 + 8.0 * A2 * rho_eta) / 16.0);

        return result;
    }

    complex<double>
    SMComponent<components::CKM>::ckm_tb() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A__ckm, _lambda__ckm, _rhobar__ckm, _etabar__ckm);
        double A2 = power_of<2>(_A__ckm()), A4 = A2 * A2;
        double lambda4 = power_of<4>(_lambda__ckm()), lambda6 = power_of<6>(_lambda__ckm()), lambda8 = lambda4 * lambda4;

        double result = 1.0 - A2 * lambda4 / 2.0 - A2 * lambda6 * std::norm(rho_eta) / 2.0 - A4 * lambda8 / 8.0;

        return complex<double>(result, 0.0);
    }

    SMComponent<components::QCD>::SMComponent(const Parameters & p, ParameterUser & u) :
        _alpha_s_Z__qcd(p["QCD::alpha_s(MZ)"], u),
        _mu_t__qcd(p["QCD::mu_t"], u),
        _mu_b__qcd(p["QCD::mu_b"], u),
        _mu_c__qcd(p["QCD::mu_c"], u),
        _lambda_qcd__qcd(p["QCD::Lambda"], u),
        _m_t_pole__qcd(p["mass::t(pole)"], u),
        _m_b_MSbar__qcd(p["mass::b(MSbar)"], u),
        _m_c_MSbar__qcd(p["mass::c"], u),
        _m_s_MSbar__qcd(p["mass::s(2GeV)"], u),
        _m_d_MSbar__qcd(p["mass::d(2GeV)"], u),
        _m_u_MSbar__qcd(p["mass::u(2GeV)"], u),
        _m_Z__qcd(p["mass::Z"], u)
    {
    }

    double
    SMComponent<components::QCD>::alpha_s(const double & mu) const
    {
        double alpha_s_0 = _alpha_s_Z__qcd, mu_0 = _m_Z__qcd;

        if (mu >= _m_Z__qcd)
        {
            if (mu < _mu_t__qcd)
                return QCD::alpha_s(mu, alpha_s_0, mu_0, QCD::beta_function_nf_5);

            alpha_s_0 = QCD::alpha_s(_mu_t__qcd, alpha_s_0, mu_0, QCD::beta_function_nf_5);
            mu_0 = _mu_t__qcd;

            return QCD::alpha_s(mu, alpha_s_0, mu_0, QCD::beta_function_nf_6);
        }

        if (mu >= _mu_b__qcd)
            return QCD::alpha_s(mu, alpha_s_0, mu_0, QCD::beta_function_nf_5);

        alpha_s_0 = QCD::alpha_s(_mu_b__qcd, alpha_s_0, mu_0, QCD::beta_function_nf_5);
        mu_0 = _mu_b__qcd;

        if (mu >= _mu_c__qcd)
            return QCD::alpha_s(mu, alpha_s_0, mu_0, QCD::beta_function_nf_4);

        alpha_s_0 = QCD::alpha_s(_mu_c__qcd, alpha_s_0, mu_0, QCD::beta_function_nf_4);
        mu_0 = _mu_c__qcd;

        if (mu >= _lambda_qcd__qcd)
            return QCD::alpha_s(mu, alpha_s_0, mu_0, QCD::beta_function_nf_3);

        throw InternalError("SMComponent<components::QCD>::alpha_s: Cannot run alpha_s to mu < lambda_qcd");
    }

    double
    SMComponent<components::QCD>::m_t_msbar(const double & mu) const
    {
        double alpha_s_m_t_pole = this->alpha_s(_m_t_pole__qcd);
        double m_t_msbar_m_t_pole = QCD::m_q_msbar(_m_t_pole__qcd, alpha_s_m_t_pole, 5.0);

        if ((_mu_b__qcd <= mu) && (mu < _mu_t__qcd))
            return QCD::m_q_msbar(m_t_msbar_m_t_pole, alpha_s_m_t_pole, this->alpha_s(mu), QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

        throw InternalError("SMComponent<components::QCD>::m_t_msbar: Running of m_t_MSbar to mu >= mu_t or to mu < m_b not yet implemented");
    }

    double
    SMComponent<components::QCD>::m_t_pole() const
    {
        return _m_t_pole__qcd();
    }

    double
    SMComponent<components::QCD>::m_b_kin(const double & mu_kin) const
    {
        double m_b_MSbar = _m_b_MSbar__qcd();
        double alpha_mu_0 = alpha_s(m_b_MSbar);

        return QCD::m_q_kin(m_b_MSbar, alpha_mu_0, mu_kin, QCD::beta_function_nf_5);
    }

    double
    SMComponent<components::QCD>::m_b_msbar(const double & mu) const
    {
        double m_b_MSbar = _m_b_MSbar__qcd();
        double alpha_mu_0 = alpha_s(m_b_MSbar);

        if (mu > m_b_MSbar)
        {
            if (mu < _mu_t__qcd)
                return QCD::m_q_msbar(m_b_MSbar, alpha_mu_0, alpha_s(mu), QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

            throw InternalError("SMComponent<components::QCD>::m_b_msbar: Running of m_b_MSbar to mu > mu_t not yet implemented");
        }
        else
        {
            if (mu >= _mu_c__qcd)
                return QCD::m_q_msbar(m_b_MSbar, alpha_mu_0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            throw InternalError("SMComponent<components::QCD>::m_b_msbar: Running of m_b_MSbar to mu < mu_c not yet implemented");
        }
    }

    double
    SMComponent<components::QCD>::m_b_pole(unsigned int loop_order) const
    {
        // The true (central) pole mass of the bottom is very close to the values
        // that can be calculated by the following quadratic polynomial.
        // This holds vor 4.13 <= m_b_MSbar <= 4.37, which corresponds to the values from [PDG2010].
        using Coefficients = std::array<double, 4>;
        static const std::array<Coefficients, 4> c = {{
            // m0,                a,                 b,                   c
            { 0.0,                0.0,               1.0,                 0.0                 }, // trivial order
            { 3.8870091768922093, 4.156247812901621, 1.2735213574282815, -0.25935468202619605 }, // loop order 1
            { 3.962932009714688,  4.38323050264802,  1.2544893957664187, -0.26527600396378315 }, // loop order 2
            { 4.19,               4.7266,            1.14485,            -0.168099            }  // loop order 3
        }};
        if (loop_order > c.size() - 1) {
            throw InternalError("SMComponent<components::QCD>::m_b_pole: maximum loop order (" + stringify(c.size() - 1) + ") exceeded (" + stringify(loop_order) + ")");
        }
        double m_b_MSbar = _m_b_MSbar__qcd();

        // Initial guess
        //                a                               m0                  b                                          m0                  c
        double m_b_pole = c[loop_order][1] + (m_b_MSbar - c[loop_order][0]) * c[loop_order][2] + power_of<2>(m_b_MSbar - c[loop_order][0]) * c[loop_order][3];

        // Iterative fixed-point procedure
        for (int i = 0 ; i < 10 ; ++i)
        {
            m_b_MSbar = m_b_msbar(m_b_pole);
            // Neglect the dependence of alpha_s on the loop order
            double next = QCD::m_q_pole(m_b_MSbar, alpha_s(m_b_pole), 5.0, loop_order);

            double delta = (m_b_pole - next) / m_b_pole;
            m_b_pole = next;

            if (std::abs(delta) < 1e-3)
                return m_b_pole;
        }

        throw InternalError("SMComponent<components::QCD>::m_b_pole: fixed-point procedure did not converge");
    }

    double
    SMComponent<components::QCD>::m_b_ps(const double & mu_f) const
    {
        double m_b_MSbar = _m_b_MSbar__qcd();

        return QCD::m_q_ps(m_b_MSbar, alpha_s(m_b_MSbar), mu_f, 5.0, QCD::beta_function_nf_5);
    }

    /* Charm */
    double
    SMComponent<components::QCD>::m_c_kin(const double & mu_kin) const
    {
        double m_c_MSbar = _m_c_MSbar__qcd();
        double alpha_mu_0 = alpha_s(m_c_MSbar);

        return QCD::m_q_kin(m_c_MSbar, alpha_mu_0, mu_kin, QCD::beta_function_nf_4);
    }

    double
    SMComponent<components::QCD>::m_c_msbar(const double & mu) const
    {
        double m_c_0 = _m_c_MSbar__qcd();
        double alpha_s_mu0 = alpha_s(m_c_0);

        if (mu >= _mu_c__qcd)
        {
            if (mu <= _mu_b__qcd)
                return QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            double alpha_s_b = alpha_s(_mu_b__qcd);
            m_c_0 = QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s_b, QCD::beta_function_nf_4, QCD::gamma_m_nf_4);
            alpha_s_mu0 = alpha_s_b;

            if (mu <= _mu_t__qcd)
                return QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

            throw InternalError("SMComponent<components::QCD>::m_c_msbar: Running of m_c_MSbar to mu > mu_t not yet implemented");
        }
        else
        {
            throw InternalError("SMComponent<components::QCD>::m_c_msbar: Running of m_c_MSbar to mu < mu_c not yet implemented");
        }
    }

    double
    SMComponent<components::QCD>::m_c_pole() const
    {
        // The true (central) pole mass of the charm is very close to the values
        // that can be calculated by the following quadratic polynomial.
        // This holds vor 1.16 <= m_c_MSbar <= 1.34, which corresponds to the values from [PDG2010].
        static const double m0 = 1.27, a = 1.59564, b = 1.13191, c = -0.737165;

        double m_c_MSbar = _m_c_MSbar__qcd();
        double m_c_pole = a + (m_c_MSbar - m0) * b + power_of<2>(m_c_MSbar - m0) * c;

        for (int i = 0 ; i < 10 ; ++i)
        {
            m_c_MSbar = m_c_msbar(m_c_pole);
            double next = QCD::m_q_pole(m_c_MSbar, alpha_s(m_c_pole), 4.0);

            double delta = (m_c_pole - next) / m_c_pole;
            m_c_pole = next;

            if (std::abs(delta) < 1e-3)
                break;
        }

        return m_c_pole;
    }

    double
    SMComponent<components::QCD>::m_s_msbar(const double & mu) const
    {
        double m_s_0 = _m_s_MSbar__qcd();
        double alpha_s_mu0 = alpha_s(2.0);

        if (mu >= 2.0)
        {
            if (mu <= _mu_b__qcd)
                return QCD::m_q_msbar(m_s_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            double alpha_s_b = alpha_s(_mu_b__qcd);
            m_s_0 = QCD::m_q_msbar(m_s_0, alpha_s_mu0, alpha_s_b, QCD::beta_function_nf_4, QCD::gamma_m_nf_4);
            alpha_s_mu0 = alpha_s_b;

            if (mu <= _mu_t__qcd)
                return QCD::m_q_msbar(m_s_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

            throw InternalError("SMComponent<components::QCD>::m_s_msbar: Running of m_s_MSbar to mu > mu_t not yet implemented");
        }
        else
        {
            if (mu >= _mu_c__qcd)
                return QCD::m_q_msbar(m_s_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            double alpha_s_c = alpha_s(_mu_c__qcd);
            double m_s_c = QCD::m_q_msbar(m_s_0, alpha_s_mu0, alpha_s_c, QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            if (mu >= 0.5)
                return QCD::m_q_msbar(m_s_c, alpha_s_c, alpha_s(mu), QCD::beta_function_nf_3, QCD::gamma_m_nf_3);

            throw InternalError("SMComponent<components::QCD>::m_s_msbar: Running of m_s_MSbar to mu < 0.5 GeV not yet implemented");
        }
    }

    double
    SMComponent<components::QCD>::m_ud_msbar(const double & mu) const
    {
        double m_ud_0 = _m_u_MSbar__qcd() + _m_d_MSbar__qcd();
        double alpha_s_mu0 = alpha_s(2.0);

        if (mu >= 2.0)
        {
            if (mu <= _mu_b__qcd)
                return QCD::m_q_msbar(m_ud_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            double alpha_s_b = alpha_s(_mu_b__qcd);
            m_ud_0 = QCD::m_q_msbar(m_ud_0, alpha_s_mu0, alpha_s_b, QCD::beta_function_nf_4, QCD::gamma_m_nf_4);
            alpha_s_mu0 = alpha_s_b;

            if (mu <= _mu_t__qcd)
                return QCD::m_q_msbar(m_ud_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

            throw InternalError("SMComponent<components::QCD>::m_ud_msbar: Running of m_ud_MSbar to mu > mu_t not yet implemented");
        }
        else
        {
            if (mu >= _mu_c__qcd)
                return QCD::m_q_msbar(m_ud_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            double alpha_s_c = alpha_s(_mu_c__qcd);
            m_ud_0 = QCD::m_q_msbar(m_ud_0, alpha_s_mu0, alpha_s_c, QCD::beta_function_nf_4, QCD::gamma_m_nf_4);
            alpha_s_mu0 = alpha_s_c;

            if (mu >= 1.0)
                return QCD::m_q_msbar(m_ud_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_3, QCD::gamma_m_nf_3);

            throw InternalError("SMComponent<components::QCD>::m_ud_msbar: Running of m_ud_MSbar to mu < 1.0 GeV not yet implemented");
        }
    }

    double
    SMComponent<components::QCD>::m_u_msbar(const double & mu) const
    {
        double m_u_0 = _m_u_MSbar__qcd();
        double alpha_s_mu0 = alpha_s(2.0);

        if (mu >= 2.0)
        {
            if (mu <= _mu_b__qcd)
                return QCD::m_q_msbar(m_u_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            double alpha_s_b = alpha_s(_mu_b__qcd);
            m_u_0 = QCD::m_q_msbar(m_u_0, alpha_s_mu0, alpha_s_b, QCD::beta_function_nf_4, QCD::gamma_m_nf_4);
            alpha_s_mu0 = alpha_s_b;

            if (mu <= _mu_t__qcd)
                return QCD::m_q_msbar(m_u_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

            throw InternalError("SMComponent<components::QCD>::m_u_msbar: Running of m_u_MSbar to mu > mu_t not yet implemented");
        }
        else
        {
            if (mu >= _mu_c__qcd)
                return QCD::m_q_msbar(m_u_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            double alpha_s_c = alpha_s(_mu_c__qcd);
            m_u_0 = QCD::m_q_msbar(m_u_0, alpha_s_mu0, alpha_s_c, QCD::beta_function_nf_4, QCD::gamma_m_nf_4);
            alpha_s_mu0 = alpha_s_c;

            if (mu >= 1.0)
                return QCD::m_q_msbar(m_u_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_3, QCD::gamma_m_nf_3);

            throw InternalError("SMComponent<components::QCD>::m_u_msbar: Running of m_u_MSbar to mu < 1.0 GeV not yet implemented");
        }
    }

    double
    SMComponent<components::QCD>::m_d_msbar(const double & mu) const
    {
        double m_d_0 = _m_d_MSbar__qcd();
        double alpha_s_mu0 = alpha_s(2.0);

        if (mu >= 2.0)
        {
            if (mu <= _mu_b__qcd)
                return QCD::m_q_msbar(m_d_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            double alpha_s_b = alpha_s(_mu_b__qcd);
            m_d_0 = QCD::m_q_msbar(m_d_0, alpha_s_mu0, alpha_s_b, QCD::beta_function_nf_4, QCD::gamma_m_nf_4);
            alpha_s_mu0 = alpha_s_b;

            if (mu <= _mu_t__qcd)
                return QCD::m_q_msbar(m_d_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

            throw InternalError("SMComponent<components::QCD>::m_d_msbar: Running of m_d_MSbar to mu > mu_t not yet implemented");
        }
        else
        {
            if (mu >= _mu_c__qcd)
                return QCD::m_q_msbar(m_d_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            double alpha_s_c = alpha_s(_mu_c__qcd);
            m_d_0 = QCD::m_q_msbar(m_d_0, alpha_s_mu0, alpha_s_c, QCD::beta_function_nf_4, QCD::gamma_m_nf_4);
            alpha_s_mu0 = alpha_s_c;

            if (mu >= 1.0)
                return QCD::m_q_msbar(m_d_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_3, QCD::gamma_m_nf_3);

            throw InternalError("SMComponent<components::QCD>::m_d_msbar: Running of m_d_MSbar to mu < 1.0 GeV not yet implemented");
        }
    }

    SMComponent<components::DeltaBS1>::SMComponent(const Parameters & p, ParameterUser & u) :
        _alpha_s_Z__deltabs1(p["QCD::alpha_s(MZ)"], u),
        _mu_t__deltabs1(p["QCD::mu_t"], u),
        _mu_b__deltabs1(p["QCD::mu_b"], u),
        _mu_c__deltabs1(p["QCD::mu_c"], u),
        _sw2__deltabs1(p["GSW::sin^2(theta)"], u),
        _m_t_pole__deltabs1(p["mass::t(pole)"], u),
        _m_W__deltabs1(p["mass::W"], u),
        _m_Z__deltabs1(p["mass::Z"], u),
        _mu_0c__deltabs1(p["b->s::mu_0c"], u),
        _mu_0t__deltabs1(p["b->s::mu_0t"], u)
    {
    }

    /* b->s Wilson coefficients */
namespace implementation
{
    /*
     * Initial scale Wilson coefficients from the charm sector, cf. [BMU1999], between Eqs. (4) and (5), pp. 4-5
     *
     * x_c = m_t(mu_0c)^2 / m_W^2
     * log_c = ln(mu_0c^2 / m_W^2)
     * sw2 = sin^2(theta_Weinberg)
     */
    std::array<complex<double>, 15>
    initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd0()
    {
        std::array<complex<double>, 15> result
        {{
            0.0, -1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0
        }};

        return result;
    }

    std::array<complex<double>, 15>
    initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd1(const double & log_c, const double & sw2)
    {
        std::array<complex<double>, 15> result
        {{
            -15.0 - 6.0 * log_c, 0.0, 0.0, 7.0/9.0 - 2.0 / 3.0 * log_c, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            23.0/36.0, 1.0/3.0, -0.25 / sw2 - 38.0/27.0, 0.25 / sw2
        }};

        return result;
    }

    std::array<complex<double>, 15>
    initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd2(const double & x_c, const double & log_c, const double & sw2)
    {
        std::array<complex<double>, 15> result;
        result.fill(0.0);
        result[0] = -(16.0 * x_c + 8.0) * std::sqrt(4.0 * x_c - 1.0) * gsl_sf_clausen(2.0 * std::asin(1.0 / 2.0 / std::sqrt(x_c)))
            + (16.0 * x_c + 20.0 / 3.0) * std::log(x_c) + 32.0 * x_c + 112.0 / 9.0
            - 7987.0 / 72.0 - 17.0 / 3.0 * M_PI * M_PI - 475.0 / 6.0 * log_c - 17.0 * log_c * log_c;
        result[1] = -127.0 / 18.0 - 4.0 / 3.0 * M_PI * M_PI - 46.0 / 3.0 * log_c - 4.0 * log_c * log_c;
        result[2] = 680.0 / 243.0 + 20.0 / 81.0 * M_PI * M_PI + 68.0 / 81.0 * log_c + 20.0 / 27.0 * log_c * log_c;
        result[3] = -950.0 / 243.0 - 10.0 / 81.0 * M_PI * M_PI - 124.0 / 27.0 * log_c - 10.0 / 27.0  * log_c * log_c;
        result[4] = -68.0 / 243.0 - 2.0 / 81.0 * M_PI * M_PI - 14.0 / 81.0 * log_c - 2.0 / 27.0 * log_c * log_c;
        result[5] = -85.0 / 162.0 - 5.0 / 108.0 * M_PI * M_PI - 35.0 / 108.0 * log_c - 5.0 / 36.0 * log_c * log_c;
        result[11] = -713.0 / 243.0 - 4.0 / 81.0 * log_c;
        result[12] = -91.0 / 324.0 + 4.0 / 27.0 * log_c;
        result[13] = -1.0 / sw2 - 524.0 / 729.0 + 128.0 / 243.0 * M_PI * M_PI + 16.0 / 3.0 * log_c + 128.0 / 81.0 * log_c * log_c;
        result[14] = 1.0 / sw2;

        return result;
    }

    /*
     * Initial scale Wilson coefficients from the top sector, cf. [BMU1999], between Eqs. (4) and (5), pp. 4-5
     *
     * x_t = m_t(mu_0t)^2 / m_W^2
     * log_t = ln(mu_0t / m_t(mu_0t))
     * sw2 = sin^2(theta_Weinberg)
     */
    std::array<complex<double>, 15>
    initial_scale_wilson_coefficients_b_to_s_top_sector_qcd0()
    {
        std::array<complex<double>, 15> result;
        result.fill(0.0);

        return result;
    }

    std::array<complex<double>, 15>
    initial_scale_wilson_coefficients_b_to_s_top_sector_qcd1(const double & x_t, const double & sw2)
    {
        std::array<complex<double>, 15> result;
        result.fill(0.0);
        result[3] = TopLoops::E0(x_t);
        result[11] = -0.5 * TopLoops::A0(x_t);
        result[12] = -0.5 * TopLoops::F0(x_t);
        result[13] = (1.0 - 4.0 * sw2) / sw2 * TopLoops::C0(x_t) - TopLoops::B0(x_t) / sw2 - TopLoops::D0(x_t);
        result[14] = (TopLoops::B0(x_t) - TopLoops::C0(x_t)) / sw2;

        return result;
    }

    std::array<complex<double>, 15>
    initial_scale_wilson_coefficients_b_to_s_top_sector_qcd2(const double & x_t, const double & log_t, const double & sw2)
    {
        std::array<complex<double>, 15> result;
        result.fill(0.0);
        result[2] = TopLoops::G1(x_t, log_t);
        result[3] = TopLoops::E1(x_t, log_t);
        result[4] = -0.1 * TopLoops::G1(x_t, log_t) + 2.0 / 15.0 * TopLoops::E0(x_t);
        result[5] = -3.0 / 16.0 * TopLoops::E1(x_t, log_t) + 0.25 * TopLoops::E0(x_t);
        result[11] = -0.5 * TopLoops::A1(x_t, log_t);
        result[12] = -0.5 * TopLoops::F1(x_t, log_t);
        result[13] = (1.0 - 4.0 * sw2) / sw2 * TopLoops::C1(x_t, log_t) - TopLoops::B1(x_t, log_t) / sw2 - TopLoops::D1(x_t, log_t);
        result[14] = (TopLoops::B1(x_t, log_t) - TopLoops::C1(x_t, log_t)) / sw2;

        return result;
    }
}

    WilsonCoefficients<BToS>
    SMComponent<components::DeltaBS1>::wilson_coefficients_b_to_s(const double & mu, const LeptonFlavor & /*lepton_flavor*/, const bool & /*cp_conjugate*/) const
    {
        /*
         * In the SM all Wilson coefficients are real-valued -> all weak phases are zero.
         * Therefore, CP conjugation leaves the Wilson coefficients invariant.
         *
         * In the SM there is lepton flavor universality.
         */

        // Calculation according to [BMU1999], Eq. (25), p. 7

        if (mu >= _mu_t__deltabs1)
            throw InternalError("SMComponent<components::DeltaB1>::wilson_coefficients_b_to_s: Evolution to mu >= mu_t is not yet implemented!");

        if (mu <= _mu_c__deltabs1)
            throw InternalError("SMComponent<components::DeltaB1>::wilson_coefficients_b_to_s: Evolution to mu <= mu_c is not yet implemented!");

        // only evolve the wilson coefficients for 5 active flavors
        static const double nf = 5.0;

        // calculate all alpha_s values
        const double alpha_s_mu_0c = QCD::alpha_s(_mu_0c__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
        const double alpha_s_mu_0t = QCD::alpha_s(_mu_0t__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);

        double alpha_s = 0.0;
        if (mu < _mu_b__deltabs1)
        {
            alpha_s = QCD::alpha_s(_mu_b__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
            alpha_s = QCD::alpha_s(mu, alpha_s, _mu_b__deltabs1, QCD::beta_function_nf_4);
        }
        else
        {
            alpha_s = QCD::alpha_s(mu, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
        }

        double alpha_s_m_t_pole = 0.0;
        if (_mu_t__deltabs1 <= _m_t_pole__deltabs1)
        {
            alpha_s_m_t_pole = QCD::alpha_s(_mu_t__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
            alpha_s_m_t_pole = QCD::alpha_s(_m_t_pole__deltabs1, alpha_s_m_t_pole, _mu_t__deltabs1, QCD::beta_function_nf_6);
        }
        else
        {
            Log::instance()->message("sm_component<deltab1>.wc", ll_error)
                << "mu_t > m_t_pole!";

            alpha_s_m_t_pole = QCD::alpha_s(_m_t_pole__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
        }

        // calculate m_t at the matching scales in the MSbar scheme
        const double m_t_msbar_m_t_pole = QCD::m_q_msbar(_m_t_pole__deltabs1, alpha_s_m_t_pole, 5.0);
        const double m_t_mu_0c = QCD::m_q_msbar(m_t_msbar_m_t_pole, alpha_s_m_t_pole, alpha_s_mu_0c, QCD::beta_function_nf_5, QCD::gamma_m_nf_5);
        const double m_t_mu_0t = QCD::m_q_msbar(m_t_msbar_m_t_pole, alpha_s_m_t_pole, alpha_s_mu_0t, QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

        // calculate dependent inputs
        const double log_c = 2.0 * std::log(_mu_0c__deltabs1 / _m_W__deltabs1), log_t = std::log(_mu_0t__deltabs1 / m_t_mu_0t);
        const double x_c = power_of<2>(m_t_mu_0c / _m_W__deltabs1), x_t = power_of<2>(m_t_mu_0t / _m_W__deltabs1);

        WilsonCoefficients<BToS> downscaled_charm = evolve(implementation::initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd0(),
                implementation::initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd1(log_c, _sw2__deltabs1),
                implementation::initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd2(x_c, log_c, _sw2__deltabs1),
                alpha_s_mu_0c, alpha_s, nf, QCD::beta_function_nf_5);
        WilsonCoefficients<BToS> downscaled_top = evolve(implementation::initial_scale_wilson_coefficients_b_to_s_top_sector_qcd0(),
                implementation::initial_scale_wilson_coefficients_b_to_s_top_sector_qcd1(x_t, _sw2__deltabs1),
                implementation::initial_scale_wilson_coefficients_b_to_s_top_sector_qcd2(x_t, log_t, _sw2__deltabs1),
                alpha_s_mu_0t, alpha_s, nf, QCD::beta_function_nf_5);

        WilsonCoefficients<BToS> wc = downscaled_top;
        wc._sm_like_coefficients = wc._sm_like_coefficients + complex<double>(-1.0, 0.0) * downscaled_charm._sm_like_coefficients;

        return wc;
    }

    SMComponent<components::WET::SBSB>::SMComponent(const Parameters & p, ParameterUser & u) :
        _G_Fermi__deltabs2(p["WET::G_Fermi"], u),
        _alpha_s_Z__deltabs2(p["QCD::alpha_s(MZ)"], u),
        _mu_t__deltabs2(p["QCD::mu_t"], u),
        _mu_b__deltabs2(p["QCD::mu_b"], u),
        _mu_c__deltabs2(p["QCD::mu_c"], u),
        _sw2__deltabs2(p["GSW::sin^2(theta)"], u),
        _m_t_pole__deltabs2(p["mass::t(pole)"], u),
        _m_W__deltabs2(p["mass::W"], u),
        _m_Z__deltabs2(p["mass::Z"], u),
        _mu_0__deltabs2(p["sbsb::mu_0"], u),
        _mu__deltabs2(p["sbsb::mu"], u)
    {
    }

    WilsonCoefficients<wc::SBSB>
    SMComponent<components::WET::SBSB>::wet_sbsb() const
    {
        if (_mu__deltabs2 >= _mu_t__deltabs2)
            throw InternalError("SMComponent<components::DeltaB1>::wilson_coefficients_sbsb: Evolution to mu >= mu_t is illdefined!");

        if (_mu__deltabs2 <= _mu_c__deltabs2)
            throw InternalError("SMComponent<components::DeltaB1>::wilson_coefficients_sbsb: Evolution to mu <= mu_c is not implemented!");

        // only evolve the wilson coefficients for 5 active flavors
        static const double nf    = 5.0;
        static const auto & beta4 = QCD::beta_function_nf_4;
        static const auto & beta5 = QCD::beta_function_nf_5;
        static const auto & beta6 = QCD::beta_function_nf_6;

        // calculate all alpha_s values
        const double alpha_s_mu_0 = QCD::alpha_s(_mu_0__deltabs2, _alpha_s_Z__deltabs2, _m_Z__deltabs2, beta5);

        double alpha_s = 0.0;
        if (_mu__deltabs2 < _mu_b__deltabs2)
        {
            alpha_s = QCD::alpha_s(_mu_b__deltabs2, _alpha_s_Z__deltabs2, _m_Z__deltabs2, beta5);
            alpha_s = QCD::alpha_s(_mu__deltabs2, alpha_s, _mu_b__deltabs2, beta4);
        }
        else
        {
            alpha_s = QCD::alpha_s(_mu__deltabs2, _alpha_s_Z__deltabs2, _m_Z__deltabs2, beta5);
        }

        double alpha_s_m_t_pole = 0.0;
        if (_mu_t__deltabs2 <= _m_t_pole__deltabs2)
        {
            alpha_s_m_t_pole = QCD::alpha_s(_mu_t__deltabs2, _alpha_s_Z__deltabs2, _m_Z__deltabs2, beta5);
            alpha_s_m_t_pole = QCD::alpha_s(_m_t_pole__deltabs2, alpha_s_m_t_pole, _mu_t__deltabs2, beta6);
        }
        else
        {
            Log::instance()->message("sm_component<deltabs2>.wc", ll_error)
                << "mu_t > m_t_pole!";

            alpha_s_m_t_pole = QCD::alpha_s(_m_t_pole__deltabs2, _alpha_s_Z__deltabs2, _m_Z__deltabs2, beta5);
        }

        // calculate m_t at the matching scale in the MSbar scheme
        const double m_t_msbar_m_t_pole = QCD::m_q_msbar(_m_t_pole__deltabs2, alpha_s_m_t_pole, nf);
        const double m_t_mu_0 = QCD::m_q_msbar(m_t_msbar_m_t_pole, alpha_s_m_t_pole, alpha_s_mu_0, beta5, QCD::gamma_m_nf_5);

        // calculate dependent inputs
        const double log_t = std::log(_mu_0__deltabs2 / _m_W__deltabs2);
        const double xt = power_of<2>(m_t_mu_0 / _m_W__deltabs2);
        const double xt2 = xt  * xt,
                     xt3 = xt2 * xt,
                     xt4 = xt2 * xt2;
        const double lnxt = std::log(xt), ln2xt = lnxt * lnxt;
        // GSL and [BBL:1995A] convention for the dilogarithm agree:
        //   gsl_sf_dilog(1.0 - x) = L_2(1.0 - x)
        const double L2 = gsl_sf_dilog(1.0 - xt);

        /*
        * Initial scale Wilson coefficients from the top sector, cf. [BBL:1995A], Eqs. (XIII.1) to (XIII.5), pp. 118.
        */

        // anomalous mass dimension, eq. (XII.7), in the five-flavor scheme
        const double Nc      = 3.0;
        const double gamma_0 = 6.0 * (Nc - 1.0) / Nc;
        const double gamma_1 = (-21.0 + 57.0 / Nc - 19.0 / 3.0 * Nc + 4.0 / 3.0 * nf) * (Nc - 1.0) /(2.0 * Nc);
        const double d5      = gamma_0  / (2.0 * beta5[0]);
        const double J5      = d5 * beta5[1] / beta5[0] - gamma_1 / (2.0 * beta5[0]);

        // one-loop (Inami-Lim) function S_0 = S_0(x_t, x_t), cf. [BBL:1995A], Eq. (XII.4), p. 101
        const double S_0 = (4.0 * xt - 11.0 * xt2 + xt3) / (4.0 * power_of<2>(1.0 - xt))
                         - 3.0 * xt3 * std::log(xt) / (2.0 * power_of<3>(1.0 - xt));
        // derivative of S_0 w.r.t. to xt
        const double S_0_d1 = (4.0 - 18.0 * xt - 3.0 * xt2 - xt3) / (4.0 * power_of<3>(1.0 - xt))
                            - 9.0 * xt2 * std::log(xt) / (2.0 * power_of<4>(1.0 - xt));

        // two-loop function, eqs. (XII.11)-(XII.14)
        const double B_t   = 5.0 * (Nc - 1.0) / (2.0 * Nc) + 3.0 * (Nc * Nc - 1.0) / (2.0 * Nc);
        // two-loop function S_1 (color singlet part)
        const double S_1_1 = -xt * (4.0 - 39.0 * xt + 168.0 * xt2 + 11.0 * xt3) / (4.0 * power_of<3>(1.0 - xt))
                           - 3.0 * xt * (4.0 - 24.0 * xt + 36.0 * xt2 + 7.0 * xt3 + xt4) / (2.0 * power_of<4>(1.0 - xt)) * lnxt
                           + 3.0 * xt3 * (13.0 + 4.0 * xt + xt2) / (2.0 * power_of<4>(1.0 - xt)) * ln2xt
                           - 3.0 * xt3 * (5.0 + xt) / power_of<3>(1.0 - xt) * L2;
        // two-loop function S_1 (color octet part)
        const double S_1_8 = -(64.0 - 68.0 * xt - 17.0 * xt2 + 11.0 * xt3) / (4.0 * power_of<2>(1.0 - xt))
                           + (32.0 - 68.0 * xt + 32.0 * xt2 - 28.0 * xt3 + 3.0 * xt4) / (2.0 * power_of<3>(1.0 - xt)) * lnxt
                           + xt2 * (4.0 - 7.0 * xt + 7.0 * xt2 - 2.0 * xt3) / (2.0 * power_of<4>(1.0 - xt)) * ln2xt
                           + 2.0 * xt * (4.0 - 7.0 * xt - 7.0 * xt2 + xt3) / power_of<3>(1.0 - xt) * L2
                           + 16.0 / xt * (M_PI * M_PI / 6.0 - L2);
        // two-loop function S_1 (full result)
        const double S_1   = (Nc - 1.0) / (2.0 * Nc) * S_1_8 + (Nc * Nc - 1.0) / (2.0 * Nc) * S_1_1;

        // auxilliary quantities
        const double eta   = std::pow(alpha_s_mu_0 / alpha_s, 6.0 / 23.0);
        // eta2B from (XIII.3), except for
        //  - a factor alpha_s(mu)^(-6/23), which has been absorbed into eta
        const double eta2B = 1.0 + alpha_s_mu_0 / (4.0 * M_PI) * (
            S_1 / S_0 + B_t - J5 + gamma_0 * log_t + 8.0 * xt * S_0_d1 / (S_0) * 2.0 * log_t
        );
        // U5(mu, mu_0) corresponds to the square brackets in (XIII.1) and (XIII.5)
        const double U5    = 1.0 + alpha_s / (4.0 * M_PI) * J5;

        // We use an effective Hamiltonian
        //   H^eff = 4 GF / sqrt(2) lambda_q^2 * C_1 * O_1,
        // where 4 * O_1 coincides with the operator Q in eq. (XIII.2).
        // We can obtain C_i from eq. (XIII.1)
        WilsonCoefficients<wc::SBSB> wc;
        wc._coefficients[0] = _G_Fermi__deltabs2 * power_of<2>(_m_W__deltabs2()) * std::sqrt(2.0) / (16.0 * M_PI * M_PI)
                  * S_0 * eta2B * eta * U5;

        return wc;
    }

    SMComponent<components::WET::UBLNu>::SMComponent(const Parameters & /* p */, ParameterUser & /* u */)
    {
    }

    WilsonCoefficients<ChargedCurrent>
    SMComponent<components::WET::UBLNu>::wet_ublnu(LeptonFlavor lepton_flavor, const bool & /* cp_conjugate */) const
    {
        // universal electroweak correction, cf. [S1982]
        // etaEW = 1 + alpha_e/pi log(m_Z/mu_b)
        // TODO: provide this to b->ulv and b->clv
        const double etaEW = 1.0066;

        WilsonCoefficients<ChargedCurrent> wc;
        wc._coefficients.fill(complex<double>(0.0));
        wc._coefficients[0] = complex<double>(etaEW);

        return wc;
    }

    SMComponent<components::WET::CBLNu>::SMComponent(const Parameters & /* p */, ParameterUser & /* c */)
    {
    }

    WilsonCoefficients<ChargedCurrent>
    SMComponent<components::WET::CBLNu>::wet_cblnu(LeptonFlavor lepton_flavor, const bool & /* cp_conjugate */) const
    {
        // universal electroweak correction, cf. [S:1982A]
        // etaEW = 1 + alpha_e/pi log(m_Z/mu_b)
        // TODO: provide this to b->ulv and b->clv
        const double etaEW = 1.0066;

        WilsonCoefficients<ChargedCurrent> wc;
        wc._coefficients.fill(complex<double>(0.0));
        wc._coefficients[0] = complex<double>(etaEW);

        return wc;
    }

    SMComponent<components::WET::SBNuNu>::SMComponent(const Parameters &  p , ParameterUser &  u) :
        _alpha_s_Z__sbnunu(p["QCD::alpha_s(MZ)"], u),
        _mu_t__sbnunu(p["QCD::mu_t"], u),
        _sw2__sbnunu(p["GSW::sin^2(theta)"], u),
        _m_t_pole__sbnunu(p["mass::t(pole)"], u),
        _m_W__sbnunu(p["mass::W"], u),
        _m_Z__sbnunu(p["mass::Z"], u),
        _mu_0__sbnunu(p["sbnunu::mu_0"], u)
    {
    }

    WilsonCoefficients<wc::SBNuNu>
    SMComponent<components::WET::SBNuNu>::wet_sbnunu(const bool & /* cp_conjugate */) const
    {
        // SM Wilson coefficients are real so cp conjugation has no effect

        // calculate alpha_s
        static const double nf    = 5.0;
        static const auto & beta5 = QCD::beta_function_nf_5;
        static const auto & beta6 = QCD::beta_function_nf_6;

        const double alpha_s_mu_0 = QCD::alpha_s(_mu_0__sbnunu, _alpha_s_Z__sbnunu, _m_Z__sbnunu, beta5);

        double alpha_s_m_t_pole = 0.0;
        if (_mu_t__sbnunu <= _m_t_pole__sbnunu)
        {
            alpha_s_m_t_pole = QCD::alpha_s(_mu_t__sbnunu, _alpha_s_Z__sbnunu, _m_Z__sbnunu, beta5);
            alpha_s_m_t_pole = QCD::alpha_s(_m_t_pole__sbnunu, alpha_s_m_t_pole, _mu_t__sbnunu, beta6);
        }
        else
        {
            Log::instance()->message("sm_component<sbnunu>.wc", ll_error)
                << "mu_t > m_t_pole!";

            alpha_s_m_t_pole = QCD::alpha_s(_m_t_pole__sbnunu, _alpha_s_Z__sbnunu, _m_Z__sbnunu, beta5);
        }

        // calculate m_t at the matching scale in the MSbar scheme
        const double m_t_msbar_m_t_pole = QCD::m_q_msbar(_m_t_pole__sbnunu, alpha_s_m_t_pole, nf);
        const double m_t_mu_0 = QCD::m_q_msbar(m_t_msbar_m_t_pole, alpha_s_m_t_pole, alpha_s_mu_0, beta5, QCD::gamma_m_nf_5);
        const double x_t = power_of<2>(m_t_mu_0 / _m_W__sbnunu);
        const double x_t2 = power_of<2>(x_t), x_t3 = power_of<3>(x_t), x_t4 = power_of<4>(x_t);

        // [BGS:2010A] TODO: implement EW corrections
        const double X_t0 = x_t / 8.0 * (
            (x_t + 2.0) / (x_t - 1.0) + (3.0 * x_t - 6.0) / power_of<2>(x_t - 1.0) * std::log(x_t)
        );
        const double X_t1 =
            - (29.0 * x_t - x_t2 - 4.0 * x_t3) / 3.0 / power_of<2>(1.0 - x_t)
            - (x_t + 9.0 * x_t2 - x_t3 - x_t4) / power_of<3>(1.0 - x_t) * std::log(x_t)
            + (8.0 * x_t + 4.0 * x_t2 + x_t3 - x_t4) / 2.0 / power_of<3>(1.0 - x_t) * power_of<2>(std::log(x_t))
            - (4.0 * x_t - x_t3) / power_of<2>(1.0 - x_t) * gsl_sf_dilog(1.0 - x_t)
            + 16.0 * x_t * std::log(_mu_t__sbnunu / _m_W__sbnunu) * (
                (8.0 - 9.0 * x_t + x_t3 + 6.0 * std::log(x_t)) / (8.0 * power_of<3>(x_t - 1.0)) // dX_t0 / dx_t
            );

        const double X_t = X_t0 + alpha_s_mu_0 / 4.0 / M_PI * X_t1;

        // [BGNS:2014A] eq. 3
        WilsonCoefficients<wc::SBNuNu> wc;
        wc._coefficients.fill(complex<double>(0.0));
        wc._coefficients[0] = X_t / _sw2__sbnunu;

        return wc;
    }

    SMComponent<components::WET::SBCU>::SMComponent(const Parameters & p , ParameterUser & u) :
        _alpha_s_Z__sbcu(p["QCD::alpha_s(MZ)"], u),
        _m_Z__sbcu(p["mass::Z"], u),
        _m_W__sbcu(p["mass::W"], u),
        _mu_0__sbcu(p["sbcu::mu_0"], u),
        _mu__sbcu(p["sbcu::mu"], u)
    {
    }

    WilsonCoefficients<wc::SBCU>
    SMComponent<components::WET::SBCU>::wet_sbcu(const bool & /* cp_conjugate */) const
    {
        // SM Wilson coefficients are real so cp conjugation has no effect

        // RGE
        static const MultiplicativeRenormalizationGroupEvolution<accuracy::NLL, 5u, 10u> rge
        {
            // gamma_0: eigenvalues
            (2.0 / 3.0) * std::array<double, 10u>{{
                -24.0, -12.0, 6.0, 3.0, (-17.0 - sqrt(241.0)), -24.0, (+1.0 + sqrt(241.0)), (+1.0 - sqrt(241.0)), 3.0, (-17 + sqrt(241.0))
            }},
            // gamma_0: V
            {{
                {{  -8.0 / 3.0, +4.0 / 3.0,  -8.0 / 3.0, +64.0 / 3.0,  0.0,                                      0.0,   0.0,                                    0.0,                                       0.0,   0.0                                   }},
                {{ -16.0,       -4.0,        -4.0,       -16.0,        0.0,                                      0.0,   0.0,                                    0.0,                                       0.0,   0.0                                   }},
                {{  +1.0 / 6.0, -1.0 / 3.0,  +2.0 / 3.0,  -4.0 / 3.0,  0.0,                                      0.0,   0.0,                                    0.0,                                       0.0,   0.0                                   }},
                {{  +1.0,       +1.0,        +1.0,        +1.0,        0.0,                                      0.0,   0.0,                                    0.0,                                       0.0,   0.0                                   }},
                {{   0.0,        0.0,         0.0,         0.0,       -53.0 / 3.0 - sqrt(241.0),               -64.0,  86.0 / 15.0 - 2.0 / 5.0 * sqrt(241.0),   2.0 / 15.0 * (43.0 + 3.0 * sqrt(241.0)),   0.0, -53.0 / 3.0 + sqrt(241.0)               }},
                {{   0.0,        0.0,         0.0,         0.0,       -16.0,                                     0.0, -16.0,                                  -16.0,                                     -64.0, -16.0                                   }},
                {{   0.0,        0.0,         0.0,         0.0,       +79.0 / 4.0 + 11.0 / 12.0 * sqrt(241.0), +16.0, (-207.0 + 7.0 * sqrt(241.0)) / 30.0,    (-207.0 - 7.0 * sqrt(241.0)) / 30.0,         0.0, +79.0 / 4.0 - 11.0 / 12.0 * sqrt(241.0) }},
                {{   0.0,        0.0,         0.0,         0.0,       +27.0 + sqrt(241.0),                       0.0,  2.0 * (51.0 - sqrt(241.0)) / 5.0,      2.0 * (51.0 + sqrt(241.0)) / 5.0,          +16.0, +27.0 - sqrt(241.0)                     }},
                {{   0.0,        0.0,         0.0,         0.0,       (53.0 + 3.0 * sqrt(241.0)) / 48.0,        +1.0, (-43.0 + 3.0 * sqrt(241.0)) / 120.0,    (-43.0 - 3.0 * sqrt(241.0)) / 120.0,         0.0,  (53.0 - 3.0 * sqrt(241.0)) / 48.0      }},
                {{   0.0,        0.0,         0.0,         0.0,       +1.0,                                      0.0,  +1.0,                                   +1.0,                                      +1.0,  +1.0                                   }}
            }},
            // gamma_1
            {{
                {{    44.0 /  9.0,  -899.0 / 3.0,  -32.0 /  9.0,  245.0 / 12.0,      0.0,              0.0,            0.0,            0.0,              0.0,            0.0         }},
                {{  -646.0 / 27.0, -2072.0 / 9.0, -115.0 / 54.0,  739.0 / 72.0,      0.0,              0.0,            0.0,            0.0,              0.0,            0.0         }},
                {{ -6848.0 /  9.0, -1344.0,        524.0 /  9.0,  178.0,             0.0,              0.0,            0.0,            0.0,              0.0,            0.0         }},
                {{     0.0,        -8468.0 / 9.0, -172.0 /  9.0,  367.0 / 18.0,      0.0,              0.0,            0.0,            0.0,              0.0,            0.0         }},
                {{     0.0,            0.0,          0.0,           0.0,         -1832.0 /  9.0,     -64.0 / 3.0,   -104.0 /  9.0,  -296.0 /   9.0,     -7.0 /  18.0,   11.0 /  48.0 }},
                {{     0.0,            0.0,          0.0,           0.0,          -128.0 / 27.0,     608.0 / 9.0,    -52.0 / 81.0, -1783.0 / 108.0,     11.0 / 216.0,   59.0 / 144.0 }},
                {{     0.0,            0.0,          0.0,           0.0,         -9488.0 / 27.0,    7108.0 / 9.0,   3052.0 /  9.0,   -31.0 /   9.0,    521.0 /  27.0, -217.0 /  36.0 }},
                {{     0.0,            0.0,          0.0,           0.0,        -25528.0 / 81.0,     896.0 / 3.0,  -6974.0 / 81.0, -4727.0 /  27.0,    863.0 / 162.0,   38.0 /   9.0 }},
                {{     0.0,            0.0,          0.0,           0.0,        -26368.0 / 27.0, -249088.0 / 9.0, -91456.0 /  9.0, 68192.0 /   9.0,  -8912.0 /  27.0, 8143.0 /   9.0 }},
                {{     0.0,            0.0,          0.0,           0.0,        510976.0 / 81.0,  -14080.0 / 9.0,   1600.0 / 81.0, 46960.0 /  27.0, -11794.0 /  81.0,  -41.0 /   6.0 }},
            }}
        };

        WilsonCoefficients<wc::SBCU> wc;
        wc._unprimed.fill(complex<double>(0.0));
        wc._primed.fill(complex<double>(0.0));

        // only unprimed WCs are non-zero in the SM
        // leading order in alpha_s
        const std::array<double, 10u> lo_unprimed = {
            -1.0 / 9.0, -2.0 / 3.0, +1.0 / 36.0, +1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        };
        // next-to-leading order in alpha_s
        const double L = 2.0 * std::log(_mu_0__sbcu / _m_W__sbcu);
        const std::array<double, 10u> nlo_unprimed = {
            +52.0 / 27.0 - 8.0 / 9.0 * L,
            -85.0 /  9.0 + 2.0 / 3.0 * L,
             -1.0 / 27.0 + 2.0 / 9.0 * L,
            +19.0 / 36.0 - 1.0 / 6.0 * L,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        };

        static const auto & beta5 = QCD::beta_function_nf_5;
        const double alpha_s_mu_0 = QCD::alpha_s(_mu_0__sbcu, _alpha_s_Z__sbcu, _m_Z__sbcu, beta5);
        const double alpha_s_mu   = QCD::alpha_s(_mu__sbcu,   _alpha_s_Z__sbcu, _m_Z__sbcu, beta5);

        const auto _unprimed = rge.evolve(alpha_s_mu, alpha_s_mu_0, lo_unprimed, nlo_unprimed);

        std::copy(_unprimed.begin(), _unprimed.end(), wc._unprimed.begin());

        return wc;
    }

    SMComponent<components::WET::DBCU>::SMComponent(const Parameters & p , ParameterUser & u) :
        _alpha_s_Z__dbcu(p["QCD::alpha_s(MZ)"], u),
        _m_Z__dbcu(p["mass::Z"], u),
        _m_W__dbcu(p["mass::W"], u),
        _mu_0__dbcu(p["dbcu::mu_0"], u),
        _mu__dbcu(p["dbcu::mu"], u)
    {
    }

    WilsonCoefficients<wc::DBCU>
    SMComponent<components::WET::DBCU>::wet_dbcu(const bool & /* cp_conjugate */) const
    {
        // SM Wilson coefficients are real so cp conjugation has no effect

        // RGE
        static const MultiplicativeRenormalizationGroupEvolution<accuracy::NLL, 5u, 10u> rge
        {
            // gamma_0: eigenvalues
            (2.0 / 3.0) * std::array<double, 10u>{{
                -24.0, -12.0, 6.0, 3.0, (-17.0 - sqrt(241.0)), -24.0, (+1.0 + sqrt(241.0)), (+1.0 - sqrt(241.0)), 3.0, (-17 + sqrt(241.0))
            }},
            // gamma_0: V
            {{
                {{  -8.0 / 3.0, +4.0 / 3.0,  -8.0 / 3.0, +64.0 / 3.0,  0.0,                                      0.0,   0.0,                                    0.0,                                       0.0,   0.0                                   }},
                {{ -16.0,       -4.0,        -4.0,       -16.0,        0.0,                                      0.0,   0.0,                                    0.0,                                       0.0,   0.0                                   }},
                {{  +1.0 / 6.0, -1.0 / 3.0,  +2.0 / 3.0,  -4.0 / 3.0,  0.0,                                      0.0,   0.0,                                    0.0,                                       0.0,   0.0                                   }},
                {{  +1.0,       +1.0,        +1.0,        +1.0,        0.0,                                      0.0,   0.0,                                    0.0,                                       0.0,   0.0                                   }},
                {{   0.0,        0.0,         0.0,         0.0,       -53.0 / 3.0 - sqrt(241.0),               -64.0,  86.0 / 15.0 - 2.0 / 5.0 * sqrt(241.0),   2.0 / 15.0 * (43.0 + 3.0 * sqrt(241.0)),   0.0, -53.0 / 3.0 + sqrt(241.0)               }},
                {{   0.0,        0.0,         0.0,         0.0,       -16.0,                                     0.0, -16.0,                                  -16.0,                                     -64.0, -16.0                                   }},
                {{   0.0,        0.0,         0.0,         0.0,       +79.0 / 4.0 + 11.0 / 12.0 * sqrt(241.0), +16.0, (-207.0 + 7.0 * sqrt(241.0)) / 30.0,    (-207.0 - 7.0 * sqrt(241.0)) / 30.0,         0.0, +79.0 / 4.0 - 11.0 / 12.0 * sqrt(241.0) }},
                {{   0.0,        0.0,         0.0,         0.0,       +27.0 + sqrt(241.0),                       0.0,  2.0 * (51.0 - sqrt(241.0)) / 5.0,      2.0 * (51.0 + sqrt(241.0)) / 5.0,          +16.0, +27.0 - sqrt(241.0)                     }},
                {{   0.0,        0.0,         0.0,         0.0,       (53.0 + 3.0 * sqrt(241.0)) / 48.0,        +1.0, (-43.0 + 3.0 * sqrt(241.0)) / 120.0,    (-43.0 - 3.0 * sqrt(241.0)) / 120.0,         0.0,  (53.0 - 3.0 * sqrt(241.0)) / 48.0      }},
                {{   0.0,        0.0,         0.0,         0.0,       +1.0,                                      0.0,  +1.0,                                   +1.0,                                      +1.0,  +1.0                                   }}
            }},
            // gamma_1
            {{
                {{    44.0 /  9.0,  -899.0 / 3.0,  -32.0 /  9.0,  245.0 / 12.0,      0.0,              0.0,            0.0,            0.0,              0.0,            0.0         }},
                {{  -646.0 / 27.0, -2072.0 / 9.0, -115.0 / 54.0,  739.0 / 72.0,      0.0,              0.0,            0.0,            0.0,              0.0,            0.0         }},
                {{ -6848.0 /  9.0, -1344.0,        524.0 /  9.0,  178.0,             0.0,              0.0,            0.0,            0.0,              0.0,            0.0         }},
                {{     0.0,        -8468.0 / 9.0, -172.0 /  9.0,  367.0 / 18.0,      0.0,              0.0,            0.0,            0.0,              0.0,            0.0         }},
                {{     0.0,            0.0,          0.0,           0.0,         -1832.0 /  9.0,     -64.0 / 3.0,   -104.0 /  9.0,  -296.0 /   9.0,     -7.0 /  18.0,   11.0 /  48.0 }},
                {{     0.0,            0.0,          0.0,           0.0,          -128.0 / 27.0,     608.0 / 9.0,    -52.0 / 81.0, -1783.0 / 108.0,     11.0 / 216.0,   59.0 / 144.0 }},
                {{     0.0,            0.0,          0.0,           0.0,         -9488.0 / 27.0,    7108.0 / 9.0,   3052.0 /  9.0,   -31.0 /   9.0,    521.0 /  27.0, -217.0 /  36.0 }},
                {{     0.0,            0.0,          0.0,           0.0,        -25528.0 / 81.0,     896.0 / 3.0,  -6974.0 / 81.0, -4727.0 /  27.0,    863.0 / 162.0,   38.0 /   9.0 }},
                {{     0.0,            0.0,          0.0,           0.0,        -26368.0 / 27.0, -249088.0 / 9.0, -91456.0 /  9.0, 68192.0 /   9.0,  -8912.0 /  27.0, 8143.0 /   9.0 }},
                {{     0.0,            0.0,          0.0,           0.0,        510976.0 / 81.0,  -14080.0 / 9.0,   1600.0 / 81.0, 46960.0 /  27.0, -11794.0 /  81.0,  -41.0 /   6.0 }},
            }}
        };

        WilsonCoefficients<wc::DBCU> wc;
        wc._unprimed.fill(complex<double>(0.0));
        wc._primed.fill(complex<double>(0.0));

        // only unprimed WCs are non-zero in the SM
        // leading order in alpha_s
        const std::array<double, 10u> lo_unprimed = {
            -1.0 / 9.0, -2.0 / 3.0, +1.0 / 36.0, +1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        };
        // next-to-leading order in alpha_s
        const double L = 2.0 * std::log(_mu_0__dbcu / _m_W__dbcu);
        const std::array<double, 10u> nlo_unprimed = {
            +52.0 / 27.0 - 8.0 / 9.0 * L,
            -85.0 /  9.0 + 2.0 / 3.0 * L,
             -1.0 / 27.0 + 2.0 / 9.0 * L,
            +19.0 / 36.0 - 1.0 / 6.0 * L,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        };

        static const auto & beta5 = QCD::beta_function_nf_5;
        const double alpha_s_mu_0 = QCD::alpha_s(_mu_0__dbcu, _alpha_s_Z__dbcu, _m_Z__dbcu, beta5);
        const double alpha_s_mu   = QCD::alpha_s(_mu__dbcu,   _alpha_s_Z__dbcu, _m_Z__dbcu, beta5);

        const auto _unprimed = rge.evolve(alpha_s_mu, alpha_s_mu_0, lo_unprimed, nlo_unprimed);

        std::copy(_unprimed.begin(), _unprimed.end(), wc._unprimed.begin());

        return wc;
    }

    SMComponent<components::WET::SCNuL>::SMComponent(const Parameters & p, ParameterUser & u) :
        _alpha_e__scnul(p["QED::alpha_e(m_c)"], u),
        _m_Z__scnul(p["mass::Z"], u),
        _mu__scnul{ UsedParameter(p["scnuee::mu"], u), UsedParameter(p["scnumumu::mu"], u), UsedParameter(p["scnutautau::mu"], u) }
    {
    }

    WilsonCoefficients<bern::ClassII>
    SMComponent<components::WET::SCNuL>::wet_scnul(LeptonFlavor lepton_flavor, const bool & /* cp_conjugate */) const
    {
        // determine renormalization scale
        const double mu = _mu__scnul[static_cast<size_t>(lepton_flavor)];

        // compute universal electroweak correction, cf. [S:1982A], eq. (1) with Qbar = 1 / 6.
        const double etaEW = 1.0 + _alpha_e__scnul / M_PI * std::log(_m_Z__scnul / mu);

        WilsonCoefficients<bern::ClassII> wc;
        wc._coefficients.fill(complex<double>(0.0));
        wc._coefficients[0] = complex<double>(etaEW);

        return wc;
    }

    SMComponent<components::WET::SB>::SMComponent(const Parameters & p , ParameterUser & u) :
        _alpha_s_Z__sbqq(p["QCD::alpha_s(MZ)"], u),
        _m_Z__sbqq(p["mass::Z"], u),
        _m_W__sbqq(p["mass::W"], u),
        _mu_0__sbqq(p["sbqq::mu_0"], u),
        _mu__sb(p["sb::mu"], u),
        _mu_t__sbqq(p["QCD::mu_t"], u),
        _m_t_pole__sbqq(p["mass::t(pole)"], u)
    {
    }

    WilsonCoefficients<wc::SBQQ>
    SMComponent<components::WET::SB>::wet_sbqq(const bool & /* cp_conjugate */) const
    {
        // SM Wilson coefficients are real so cp conjugation has no effect

        // RGE
        static const MultiplicativeRenormalizationGroupEvolution<accuracy::NLL, 5u, 11u> rge
        {
            // gamma_0: eigenvalues
            std::array<double, 11u>{{
                -19.494503, -13.790721, -12.819708, 12.093029, -8.0000000, -6.4858258, 6.2654910, 4.0000000, 4.0000000, 2.2332779, 2.2211815
            }},
            // gamma_0: V
            {{
                {{0,0,0,0,-1.00000000,0,0,0,1.00000000,0,0}},
                {{0,0,0,0,0.333333333,0,0,0,0.666666667,0,0}},
                {{-0.00888718779,-0.0107035682,1.51315026,0.135698334,-0.0370370370,-0.0200744548,-0.113763345,0.0634920635,0.0317460317,0.0810493047,1.06027001}},
                {{-0.0565215183,-0.0655927704,9.74852994,-0.241743497,0.111111111,0.00583079358,0.0408844728,0.0952380952,0.0476190476,-0.0763605912,-0.986766989}},
                {{0.000663783899,0.000483358181,-0.0463214633,-0.0162052853,0.00925925926,0.00336138495,0.0161441871,-0.0158730159,-0.00793650794,-0.00411591439,-0.0546038906}},
                {{0.00274062101,0.00505945904,-0.841241774,-0.0199191229,-0.0277777778,-0.00548838474,0.0186705861,-0.0238095238,-0.0119047619,0.00556786345,0.0713875393}},
                {{-0.0322597457,0,1.82885075,-0.165870115,0,0,0,0,0,0,-1.63072088}},
                {{-0.417566380,0,2.32733294,-2.77515573,0,0,0,0,0,0,0.865389169}},
                {{-0.0199827510,0,-0.559744385,0.453850259,0,0,0,0,0,0,0.125876878}},
                {{0.0431354138,0,0.460604211,0.504137026,0,0,0,0,0,0,-0.00787665082}},
                {{0,0,0,0,0,0,0,1.00000000,0,0,0}}
            }},
            // gamma_1
            {{
                {{-355. / 9., -502. / 27., -1412. / 243., -1369. / 243., 134. / 243., -35. / 162., 0., 0., 0., 0., 0}},
                {{-35. / 3., -28. / 3., -416. / 81., 1280. / 81., 56. / 81., 35. / 27., 0., 0., 0., 0., 0}},
                {{0., 0., -4468. / 81., -31469. / 81., 400. / 81., 3373. / 108., 0., 0., 0., 0., 0}},
                {{0., 0., -8158. / 243., -59399. / 243., 269. / 483., 12899. / 648., 0., 0., 0., 0., 0}},
                {{0., 0., -251680. / 81., -128648. / 81., 23836. / 81., 6106. / 27., 0., 0., 0., 0., 0}},
                {{0., 0., 58640. / 243., -26348. / 243., -14324. / 243., -2551. / 162., 0., 0., 0., 0., 0}},
                {{0., 0., 832. / 243., -4000. / 243., -112. / 243., -70. / 81., -404. / 9., -3077. / 9., 32. / 9., 1031. / 36., 0}},
                {{0., 0., 3376. / 729., 6344. / 729., -280. / 729., 55. / 486., -2698. / 81., -8035. / 27., -49. / 162., 4493. / 216., 0}},
                {{0., 0., 2272. / 243., -72088. / 243., -688. / 243., -1240. / 81., -19072. / 9., -14096. / 9., 1708. / 9., 1622. / 9., 0}},
                {{0., 0., 45424. / 729., 84236. / 729., -3880. / 729., 1220. / 243., 32288. / 81., -15976. / 27., -6692. / 81., -2437. / 54., 0}},
                {{0., 0., -1576. / 81., 446. / 27., 172. / 81., 40. / 27., 0., 0., 0., 0., 325. / 9.}}
            }}
        };

        WilsonCoefficients<wc::SBQQ> wc;
        wc._unprimed.fill(complex<double>(0.0));
        wc._primed.fill(complex<double>(0.0));

        // only unprimed WCs are non-zero in the SM
        // leading order in alpha_s
        const std::array<double, 11u> lo_unprimed = {
            0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.
        };
        // next-to-leading order in alpha_s

        // calculate alpha_s
        static const double nf    = 5.0;
        static const auto & beta5 = QCD::beta_function_nf_5;
        static const auto & beta6 = QCD::beta_function_nf_6;
        const double alpha_s_mu_0 = QCD::alpha_s(_mu_0__sbqq, _alpha_s_Z__sbqq, _m_Z__sbqq, beta5);

        double alpha_s_m_t_pole = 0.0;
        if (_mu_t__sbqq <= _m_t_pole__sbqq)
        {
            alpha_s_m_t_pole = QCD::alpha_s(_mu_t__sbqq, _alpha_s_Z__sbqq, _m_Z__sbqq, beta5);
            alpha_s_m_t_pole = QCD::alpha_s(_m_t_pole__sbqq, alpha_s_m_t_pole, _mu_t__sbqq, beta6);
        }
        else
        {
            Log::instance()->message("sm_component<sbqq>.wc", ll_error)
                << "mu_t > m_t_pole!";

            alpha_s_m_t_pole = QCD::alpha_s(_m_t_pole__sbqq, _alpha_s_Z__sbqq, _m_Z__sbqq, beta5);
        }

        // calculate m_t at the matching scale in the MSbar scheme
        const double m_t_msbar_m_t_pole = QCD::m_q_msbar(_m_t_pole__sbqq, alpha_s_m_t_pole, nf);
        const double m_t_mu_0 = QCD::m_q_msbar(m_t_msbar_m_t_pole, alpha_s_m_t_pole, alpha_s_mu_0, beta5, QCD::gamma_m_nf_5);
        const double xt = power_of<2>(m_t_mu_0/_m_W__sbqq);

        const double E0 = ( (8 - 42*xt + 35*power_of<2>(xt) - 7*power_of<3>(xt)) / (12*power_of<3>(xt-1)) ) - (4 - 16*xt+ 9*power_of<2>(xt))*std::log(xt) / (6*power_of<4>(xt-1));
        const double L = 2.0 * std::log(_mu_0__sbqq / _m_W__sbqq);
        const std::array<double, 11u> nlo_unprimed = {
            +15.0 + 6.0 * L,
            0.,
            0.,
            E0    + 2.0 / 3.0 * L,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
        };

        const double alpha_s_mu   = QCD::alpha_s(_mu__sb,     _alpha_s_Z__sbqq, _m_Z__sbqq, beta5);

        const auto _unprimed = rge.evolve(alpha_s_mu, alpha_s_mu_0, lo_unprimed, nlo_unprimed);

        std::copy(_unprimed.begin(), _unprimed.end(), wc._unprimed.begin());

        return wc;
    }

    StandardModel::StandardModel(const Parameters & p) :
        SMComponent<components::CKM>(p, *this),
        SMComponent<components::QCD>(p, *this),
        SMComponent<components::WET::SBSB>(p, *this),
        SMComponent<components::DeltaBS1>(p, *this),
        SMComponent<components::WET::CBLNu>(p, *this),
        SMComponent<components::WET::UBLNu>(p, *this),
        SMComponent<components::WET::SBNuNu>(p, *this),
        SMComponent<components::WET::SBCU>(p, *this),
        SMComponent<components::WET::DBCU>(p, *this),
        SMComponent<components::WET::SCNuL>(p, *this),
        SMComponent<components::WET::SB>(p, *this)
    {
    }

    StandardModel::~StandardModel()
    {
    }

    std::shared_ptr<Model>
    StandardModel::make(const Parameters & parameters, const Options &)
    {
        return std::shared_ptr<Model>(new StandardModel(parameters));
    }
}
