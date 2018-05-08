/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2017 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2018 Christoph Bobeth
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

#include <eos/utils/log.hh>
#include <eos/utils/matrix.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/top-loops.hh>
#include <eos/utils/standard-model.hh>

#include <cmath>

#include <gsl/gsl_sf_clausen.h>

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
        _m_ud_MSbar__qcd(p["mass::ud(2GeV)"], u),
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
    SMComponent<components::QCD>::m_b_pole() const
    {
        // The true (central) pole mass of the bottom is very close to the values
        // that can be calculated by the following quadratic polynomial.
        // This holds vor 4.13 <= m_b_MSbar <= 4.37, which corresponds to the values from [PDG2010].
        static const double m0 = 4.19, a = 4.7266, b = 1.14485, c = -0.168099;
        double m_b_MSbar = _m_b_MSbar__qcd();
        double m_b_pole = a + (m_b_MSbar - m0) * b + power_of<2>(m_b_MSbar - m0) * c;

        for (int i = 0 ; i < 10 ; ++i)
        {
            m_b_MSbar = m_b_msbar(m_b_pole);
            double next = QCD::m_q_pole(m_b_MSbar, alpha_s(m_b_pole), 5.0);

            double delta = (m_b_pole - next) / m_b_pole;
            m_b_pole = next;

            if (std::abs(delta) < 1e-3)
                break;
        }

        return m_b_pole;
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
            throw InternalError("SMComponent<components::QCD>::m_s_msbar: Running of m_s_MSbar to mu < 2.0 GeV not yet implemented");
        }
    }

    double
    SMComponent<components::QCD>::m_ud_msbar(const double & mu) const
    {
        double m_ud_0 = _m_ud_MSbar__qcd();
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
        _mu_0t__deltabs1(p["b->s::mu_0t"], u),
        _mu__deltabs1(p["mu"], u)
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
    SMComponent<components::DeltaBS1>::wilson_coefficients_b_to_s(const std::string & /*lepton_flavour*/, const bool & /*cp_conjugate*/) const
    {
        /*
         * In the SM all Wilson coefficients are real-valued -> all weak phases are zero.
         * Therefore, CP conjugation leaves the Wilson coefficients invariant.
         *
         * In the SM there is lepton flavour universality.
         */

        // Calculation according to [BMU1999], Eq. (25), p. 7

        if (_mu__deltabs1 >= _mu_t__deltabs1)
            throw InternalError("SMComponent<components::DeltaB1>::wilson_coefficients_b_to_s: Evolution to mu >= mu_t is not yet implemented!");

        if (_mu__deltabs1 <= _mu_c__deltabs1)
            throw InternalError("SMComponent<components::DeltaB1>::wilson_coefficients_b_to_s: Evolution to mu <= mu_c is not yet implemented!");

        // only evolve the wilson coefficients for 5 active flavors
        static const double nf = 5.0;

        // calculate all alpha_s values
        const double alpha_s_mu_0c = QCD::alpha_s(_mu_0c__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
        const double alpha_s_mu_0t = QCD::alpha_s(_mu_0t__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);

        double alpha_s = 0.0;
        if (_mu__deltabs1 < _mu_b__deltabs1)
        {
            alpha_s = QCD::alpha_s(_mu_b__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
            alpha_s = QCD::alpha_s(_mu__deltabs1, alpha_s, _mu_b__deltabs1, QCD::beta_function_nf_4);
        }
        else
        {
            alpha_s = QCD::alpha_s(_mu__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
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

    SMComponent<components::DeltaBU1>::SMComponent(const Parameters & /* p */, ParameterUser & /* u */)
    {
    }

    WilsonCoefficients<BToU>
    SMComponent<components::DeltaBU1>::wilson_coefficients_b_to_u(const std::string & lepton_flavour, const bool & /* cp_conjugate */) const
    {
        WilsonCoefficients<BToU> wc;
        wc._coefficients.fill(complex<double>(0.0));
        wc._coefficients[0] = complex<double>(1.0);

        return wc;
    }
    
    SMComponent<components::DeltaBC1>::SMComponent(const Parameters & /* p */, ParameterUser & /* c */)
    {
    }
    
    WilsonCoefficients<BToC>
    SMComponent<components::DeltaBC1>::wilson_coefficients_b_to_c(const std::string & lepton_flavour, const bool & /* cp_conjugate */) const
    {
        WilsonCoefficients<BToC> wc;
        wc._coefficients.fill(complex<double>(0.0));
        wc._coefficients[0] = complex<double>(1.0);
        
        return wc;
    }

    StandardModel::StandardModel(const Parameters & p) :
        SMComponent<components::CKM>(p, *this),
        SMComponent<components::QCD>(p, *this),
        SMComponent<components::DeltaBS1>(p, *this),
        SMComponent<components::DeltaBU1>(p, *this),
        SMComponent<components::DeltaBC1>(p, *this)
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
