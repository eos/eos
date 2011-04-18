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

#include <src/utils/standard-model.hh>
#include <src/utils/power_of.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>

namespace eos
{
    using std::sqrt;

    StandardModel::StandardModel(const Parameters & p) :
        _A(p["CKM::A"]),
        _lambda(p["CKM::lambda"]),
        _rhobar(p["CKM::rhobar"]),
        _etabar(p["CKM::etabar"]),
        _alpha_s_Z(p["QCD::alpha_s(MZ)"]),
        _mu_t(p["QCD::mu_t"]),
        _mu_b(p["QCD::mu_b"]),
        _mu_c(p["QCD::mu_c"]),
        _lambda_qcd(p["QCD::Lambda"]),
        _m_b_MSbar(p["mass::b(MSbar)"]),
        _m_c_MSbar(p["mass::c"]),
        _m_Z(p["mass::Z"]),
        _mu(p["mu"])
    {
    }

    StandardModel::~StandardModel()
    {
    }

    std::shared_ptr<Model>
    StandardModel::make(const Parameters & parameters)
    {
        return std::shared_ptr<Model>(new StandardModel(parameters));
    }

    double
    StandardModel::alpha_s(const double & mu) const
    {
        if (mu > _m_Z)
            throw InternalError("StandardModel::alpha_s: Running of alpha_s to mu > m_Z not yet implemented");

        double alpha_s_0 = _alpha_s_Z, mu_0 = _m_Z;

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
}
