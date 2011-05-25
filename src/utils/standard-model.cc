/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <src/utils/matrix.hh>
#include <src/utils/power_of.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>
#include <src/utils/top-loops.hh>
#include <src/utils/standard-model.hh>

#include <cmath>

#include <gsl/gsl_sf_clausen.h>

namespace eos
{
    using std::sqrt;

    StandardModel::StandardModel(const Parameters & p) :
        _A(p["CKM::A"], *this),
        _lambda(p["CKM::lambda"], *this),
        _rhobar(p["CKM::rhobar"], *this),
        _etabar(p["CKM::etabar"], *this),
        _alpha_s_Z(p["QCD::alpha_s(MZ)"], *this),
        _mu_t(p["QCD::mu_t"], *this),
        _mu_b(p["QCD::mu_b"], *this),
        _mu_c(p["QCD::mu_c"], *this),
        _lambda_qcd(p["QCD::Lambda"], *this),
        _m_t_pole(p["mass::t(pole)"], *this),
        _m_b_MSbar(p["mass::b(MSbar)"], *this),
        _m_c_MSbar(p["mass::c"], *this),
        _m_W(p["mass::W"], *this),
        _m_Z(p["mass::Z"], *this),
        _sw2(p["GSW::sin^2(theta)"], *this),
        _mu(p["mu"], *this),
        _mu_0c(p["b->s::mu_0c"], *this),
        _mu_0t(p["b->s::mu_0t"], *this)
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

    double
    StandardModel::alpha_s(const double & mu) const
    {
        double alpha_s_0 = _alpha_s_Z, mu_0 = _m_Z;

        if (mu >= _m_Z)
        {
            if (mu < _mu_t)
                return QCD::alpha_s(mu, alpha_s_0, mu_0, QCD::beta_function_nf_5);

            alpha_s_0 = QCD::alpha_s(_mu_t, alpha_s_0, mu_0, QCD::beta_function_nf_5);
            mu_0 = _mu_t;

            return QCD::alpha_s(mu, alpha_s_0, mu_0, QCD::beta_function_nf_6);
        }

        if (mu >= _mu_b)
            return QCD::alpha_s(mu, alpha_s_0, mu_0, QCD::beta_function_nf_5);

        alpha_s_0 = QCD::alpha_s(_mu_b, alpha_s_0, mu_0, QCD::beta_function_nf_5);
        mu_0 = _mu_b;

        if (mu >= _mu_c)
            return QCD::alpha_s(mu, alpha_s_0, mu_0, QCD::beta_function_nf_4);

        alpha_s_0 = QCD::alpha_s(_mu_c, alpha_s_0, mu_0, QCD::beta_function_nf_4);
        mu_0 = _mu_c;

        if (mu >= _lambda_qcd)
            return QCD::alpha_s(mu, alpha_s_0, mu_0, QCD::beta_function_nf_3);

        throw InternalError("StandardModel::alpha_s: Cannot run alpha_s to mu < lambda_qcd");
    }

    double
    StandardModel::m_t_msbar(const double & mu) const
    {
        double alpha_s_m_t_pole = this->alpha_s(_m_t_pole);
        double m_t_msbar_m_t_pole = QCD::m_q_msbar(_m_t_pole, alpha_s_m_t_pole, 5.0);

        if ((_mu_b <= mu) && (mu < _mu_t))
            return QCD::m_q_msbar(m_t_msbar_m_t_pole, alpha_s_m_t_pole, this->alpha_s(mu), QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

        throw InternalError("StandardModel::m_t_msbar: Running of m_t_MSbar to mu >= mu_t or to mu < m_b not yet implemented");
    }

    double
    StandardModel::m_t_pole() const
    {
        return _m_t_pole();
    }

    double
    StandardModel::m_b_msbar(const double & mu) const
    {
        double m_b_MSbar = _m_b_MSbar();
        double alpha_mu_0 = alpha_s(m_b_MSbar);

        if (mu > _m_b_MSbar)
        {
            if (mu < _mu_t)
                return QCD::m_q_msbar(m_b_MSbar, alpha_mu_0, alpha_s(mu), QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

            throw InternalError("StandardModel::m_b_msbar: Running of m_b_MSbar to mu > mu_t not yet implemented");
        }
        else
        {
            if (mu >= _mu_c)
                return QCD::m_q_msbar(m_b_MSbar, alpha_mu_0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            throw InternalError("StandardModel::m_b_msbar: Running of m_b_MSbar to mu < mu_c not yet implemented");
        }
    }

    double
    StandardModel::m_b_pole() const
    {
        double m_b_MSbar = _m_b_MSbar();

        return QCD::m_q_pole(m_b_MSbar, alpha_s(m_b_MSbar), 5.0);
    }

    double
    StandardModel::m_b_ps(const double & mu_f) const
    {
        double m_b_MSbar = _m_b_MSbar();

        return QCD::m_q_ps(m_b_MSbar, alpha_s(m_b_MSbar), mu_f, 5.0, QCD::beta_function_nf_5);
    }

    /* Charm */
    double
    StandardModel::m_c_msbar(const double & mu) const
    {
        double m_c_0 = _m_c_MSbar();
        double alpha_s_mu0 = alpha_s(m_c_0);

        if (mu >= _m_c_MSbar)
        {
            if (mu <= _mu_b)
                return QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_4, QCD::gamma_m_nf_4);

            double alpha_s_b = alpha_s(_mu_b);
            m_c_0 = QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s_b, QCD::beta_function_nf_4, QCD::gamma_m_nf_4);
            alpha_s_mu0 = alpha_s_b;

            if (mu <= _mu_t)
                return QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s(mu), QCD::beta_function_nf_5, QCD::gamma_m_nf_5);

            throw InternalError("StandardModel::m_c_msbar: Running of m_c_MSbar to mu > mu_t not yet implemented");
        }
        else
        {
            throw InternalError("StandardModel::m_c_msbar: Running of m_c_MSbar to mu < mu_c not yet implemented");
        }
    }

    double
    StandardModel::m_c_pole() const
    {
        double m_c_MSbar = _m_c_MSbar();

        return QCD::m_q_pole(m_c_MSbar, alpha_s(m_c_MSbar), 4.0);
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
    StandardModel::ckm_cd() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A, _lambda, _rhobar, _etabar);
        double A2 = power_of<2>(_A()), lambda4 = power_of<4>(_lambda()), lambda6 = power_of<6>(_lambda());

        complex<double> result = _lambda() * (1.0 - A2 * lambda4 * (1.0 - 2.0 * rho_eta) / 2.0 + A2 * lambda6 * rho_eta / 2.0);

        return result;
    }

    complex<double>
    StandardModel::ckm_cs() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A, _lambda, _rhobar, _etabar);
        double A2 = power_of<2>(_A()), A4 = power_of<2>(A2);
        double lambda2 = power_of<2>(_lambda()), lambda4 = power_of<2>(lambda2), lambda6 = lambda4 * lambda2, lambda8 = lambda4 * lambda4;;

        complex<double> result = 1.0 - lambda2 / 2.0 - lambda4 * (1.0 + 4.0 * A2) / 8.0
            - lambda6 * (1.0 - 4.0 * A2 + 16.0 * A2 * rho_eta) / 16.0 - lambda8 * (5.0 - 8.0 * A2 + 16.0 * A4) / 128.0;

        return result;
    }

    complex<double>
    StandardModel::ckm_cb() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A, _lambda, _rhobar, _etabar);
        double A2 = power_of<2>(_A()), lambda2 = power_of<2>(_lambda()), lambda6 = power_of<3>(lambda2);

        double result = _A * lambda2 * (1.0 - 0.5 * A2 * lambda6 * std::norm(rho_eta));

        return complex<double>(result, 0.0);
    }

    complex<double>
    StandardModel::ckm_ud() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A, _lambda, _rhobar, _etabar);
        double A2 = power_of<2>(_A()), lambda2 = power_of<2>(_lambda()), lambda4 = lambda2 * lambda2,
               lambda6 = lambda2 * lambda4, lambda8 = lambda4 * lambda4;

        double result = 1.0 - lambda2 / 2.0 - lambda4 / 8.0 - lambda6 * (1.0 + 8.0 * A2 * std::norm(rho_eta)) / 16.0
            - lambda8 * (5.0 - 32.0 * A2 * std::norm(rho_eta)) / 128.0;

        return complex<double>(result, 0.0);
    }

    complex<double>
    StandardModel::ckm_us() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A, _lambda, _rhobar, _etabar);
        double A2 = power_of<2>(_A()), lambda6 = power_of<6>(_lambda());

        double result = _lambda * (1.0 - 0.5 * A2 * lambda6 * std::norm(rho_eta));

        return complex<double>(result, 0.0);
    }

    complex<double>
    StandardModel::ckm_ub() const
    {
        complex<double> rho_eta_conj = std::conj(implementation::rho_eta(_A, _lambda, _rhobar, _etabar));

        complex<double> result = _A * power_of<3>(_lambda()) * rho_eta_conj;

        return result;
    }

    complex<double>
    StandardModel::ckm_td() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A, _lambda, _rhobar, _etabar);
        double A2 = power_of<2>(_A()), lambda2 = power_of<2>(_lambda()), lambda3 = _lambda() * lambda2, lambda4 = lambda2 * lambda2;

        complex<double> result = _A * lambda3 *
            ((1.0 - rho_eta) + lambda2 * rho_eta / 2.0 + lambda4 * (1.0 + 4.0 * A2) * rho_eta / 8.0);

        return result;
    }

    complex<double>
    StandardModel::ckm_ts() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A, _lambda, _rhobar, _etabar);
        double A2 = power_of<2>(_A()), lambda2 = power_of<2>(_lambda()), lambda4 = lambda2 * lambda2, lambda6 = lambda2 * lambda4;

        complex<double> result = -1.0 * _A * lambda2 *
            (1.0 - lambda2 * (1.0 - 2.0 * rho_eta) / 2.0 - lambda4 / 8.0 - lambda6 * (1.0 + 8.0 * A2 * rho_eta) / 16.0);

        return result;
    }

    complex<double>
    StandardModel::ckm_tb() const
    {
        complex<double> rho_eta = implementation::rho_eta(_A, _lambda, _rhobar, _etabar);
        double A2 = power_of<2>(_A()), A4 = A2 * A2;
        double lambda4 = power_of<4>(_lambda()), lambda6 = power_of<6>(_lambda()), lambda8 = lambda4 * lambda4;

        double result = 1.0 - A2 * lambda4 / 2.0 - A2 * lambda6 * std::norm(rho_eta) / 2.0 - A4 * lambda8 / 8.0;

        return complex<double>(result, 0.0);
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
    std::array<double, 15>
    initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd0()
    {
        std::array<double, 15> result
        {{
            0.0, -1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0
        }};

        return result;
    }

    std::array<double, 15>
    initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd1(const double & log_c, const double & sw2)
    {
        std::array<double, 15> result
        {{
            -15.0 - 6.0 * log_c, 0.0, 0.0, 7.0/9.0 - 2.0 / 3.0 * log_c, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            23.0/36.0, 1.0/3.0, -0.25 / sw2 - 38.0/27.0, 0.25 / sw2
        }};

        return result;
    }

    std::array<double, 15>
    initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd2(const double & x_c, const double & log_c, const double & sw2)
    {
        std::array<double, 15> result;
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
    std::array<double, 15>
    initial_scale_wilson_coefficients_b_to_s_top_sector_qcd0()
    {
        std::array<double, 15> result;
        result.fill(0.0);

        return result;
    }

    std::array<double, 15>
    initial_scale_wilson_coefficients_b_to_s_top_sector_qcd1(const double & x_t, const double & sw2)
    {
        std::array<double, 15> result;
        result.fill(0.0);
        result[3] = TopLoops::E0(x_t);
        result[11] = -0.5 * TopLoops::A0(x_t);
        result[12] = -0.5 * TopLoops::F0(x_t);
        result[13] = (1.0 - 4.0 * sw2) / sw2 * TopLoops::C0(x_t) - TopLoops::B0(x_t) / sw2 - TopLoops::D0(x_t);
        result[14] = (TopLoops::B0(x_t) - TopLoops::C0(x_t)) / sw2;

        return result;
    }

    std::array<double, 15>
    initial_scale_wilson_coefficients_b_to_s_top_sector_qcd2(const double & x_t, const double & log_t, const double & sw2)
    {
        std::array<double, 15> result;
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
    StandardModel::wilson_coefficients_b_to_s() const
    {
        // Calculation according to [BMU1999], Eq. (25), p. 7

        if (_mu >= _mu_t)
            throw InternalError("StandardModel::wilson_coefficients_b_to_s: Evolution to mu >= mu_t is not yet implemented!");

        if (_mu <= _mu_c)
            throw InternalError("StandardModel::wilson_coefficients_b_to_s: Evolution to mu <= mu_c is not yet implemented!");

        static const double nf = 5.0;
        const double m_t_mu_0c = this->m_t_msbar(_mu_0c), m_t_mu_0t = this->m_t_msbar(_mu_0t);
        const double log_c = 2.0 * std::log(_mu_0c / _m_W), log_t = std::log(_mu_0t / m_t_mu_0t);
        const double x_c = power_of<2>(m_t_mu_0c / _m_W), x_t = power_of<2>(m_t_mu_0t / _m_W);

        double alpha_s_0 = this->alpha_s(_m_W);
        double alpha_s = this->alpha_s(_mu);

        WilsonCoefficients<BToS> downscaled_charm = evolve(implementation::initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd0(),
                implementation::initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd1(log_c, _sw2),
                implementation::initial_scale_wilson_coefficients_b_to_s_charm_sector_qcd2(x_c, log_c, _sw2),
                alpha_s_0, alpha_s, nf, QCD::beta_function_nf_5);
        WilsonCoefficients<BToS> downscaled_top = evolve(implementation::initial_scale_wilson_coefficients_b_to_s_top_sector_qcd0(),
                implementation::initial_scale_wilson_coefficients_b_to_s_top_sector_qcd1(x_t, _sw2),
                implementation::initial_scale_wilson_coefficients_b_to_s_top_sector_qcd2(x_t, log_t, _sw2),
                alpha_s_0, alpha_s, nf, QCD::beta_function_nf_5);

        WilsonCoefficients<BToS> wc = downscaled_top;
        wc._coefficients = wc._coefficients + (-1.0) * downscaled_charm._coefficients;

        return wc;
    }
}
