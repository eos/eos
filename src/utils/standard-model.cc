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
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>

namespace eos
{
    using std::sqrt;

    struct StandardModelImplementation
    {
        // 5 flavor QCD
        static const QCD::BetaFunction qcd_nf5_beta_function;
        static const QCD::AnomalousMassDimension qcd_nf5_gamma_m;

        // 4 flavor QCD
        static const QCD::BetaFunction qcd_nf4_beta_function;
        static const QCD::AnomalousMassDimension qcd_nf4_gamma_m;

        // 3 flavor QCD
        static const QCD::BetaFunction qcd_nf3_beta_function;
    };

    /* 5 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 5
    const QCD::BetaFunction StandardModelImplementation::qcd_nf5_beta_function
    {{
        23.0 / 3.0,
        116.0 / 3.0,
        9769.0 / 54.0,
        4826.1563287908967,
    }};

    // cf. [CKS2000], Eq. (7), p. 5 with n_f = 5
    const QCD::AnomalousMassDimension StandardModelImplementation::qcd_nf5_gamma_m
    {{
        1.0,
        506.0 / 9.0,
        474.87124557719461,
        2824.7862379694232,
    }};

    /* 4 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 4
    const QCD::BetaFunction StandardModelImplementation::qcd_nf4_beta_function
    {{
        25.0 / 3.0,
        154.0 / 3.0,
        21943.0 / 54.0,
        8035.1864197901160,
    }};

    // cf. [CKS2000], Eq. (7), p. 5 with n_f = 4
    const QCD::AnomalousMassDimension StandardModelImplementation::qcd_nf4_gamma_m
    {{
        1.0,
        526.0 / 9.0,
        636.61057670866927,
        6989.5510103599477,
    }};

    /* 3 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 3
    const QCD::BetaFunction StandardModelImplementation::qcd_nf3_beta_function
    {{
        9.0,
        64.0,
        3863.0 / 6.0,
        12090.378130803711,
    }};

    StandardModel::StandardModel(const Parameters & p) :
        _A(p["CKM::A"]),
        _lambda(p["CKM::lambda"]),
        _rhobar(p["CKM::rhobar"]),
        _etabar(p["CKM::etabar"]),
        _ckm_cb(p["CKM::|V_cb|"]),
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
            return QCD::alpha_s(mu, alpha_s_0, mu_0, StandardModelImplementation::qcd_nf5_beta_function);

        alpha_s_0 = QCD::alpha_s(_mu_b, alpha_s_0, mu_0, StandardModelImplementation::qcd_nf5_beta_function);
        mu_0 = _mu_b;

        if (mu >= _mu_c)
            return QCD::alpha_s(mu, alpha_s_0, mu_0, StandardModelImplementation::qcd_nf4_beta_function);

        alpha_s_0 = QCD::alpha_s(_mu_c, alpha_s_0, mu_0, StandardModelImplementation::qcd_nf4_beta_function);
        mu_0 = _mu_c;

        if (mu >= _lambda_qcd)
            return QCD::alpha_s(mu, alpha_s_0, mu_0, StandardModelImplementation::qcd_nf3_beta_function);

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
                return QCD::m_q_msbar(m_b_MSbar, alpha_mu_0, alpha_s(mu), StandardModelImplementation::qcd_nf5_beta_function, StandardModelImplementation::qcd_nf5_gamma_m);

            throw InternalError("StandardModel::m_b_msbar: Running of m_b_MSbar to mu > mu_t not yet implemented");
        }
        else
        {
            if (mu >= _mu_c)
                return QCD::m_q_msbar(m_b_MSbar, alpha_mu_0, alpha_s(mu), StandardModelImplementation::qcd_nf4_beta_function, StandardModelImplementation::qcd_nf4_gamma_m);

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

        return QCD::m_q_ps(m_b_MSbar, alpha_s(m_b_MSbar), mu_f, 5.0, StandardModelImplementation::qcd_nf5_beta_function);
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
                return QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s(mu), StandardModelImplementation::qcd_nf4_beta_function, StandardModelImplementation::qcd_nf4_gamma_m);

            double alpha_s_b = alpha_s(_mu_b);
            m_c_0 = QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s_b, StandardModelImplementation::qcd_nf4_beta_function, StandardModelImplementation::qcd_nf4_gamma_m);
            alpha_s_mu0 = alpha_s_b;

            if (mu <= _mu_t)
                return QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s(mu), StandardModelImplementation::qcd_nf5_beta_function, StandardModelImplementation::qcd_nf5_gamma_m);

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

    complex<double>
    StandardModel::ckm_cb() const
    {
        return complex<double>(_ckm_cb, 0.0);
    }

    complex<double>
    StandardModel::ckm_us() const
    {
        return complex<double>(_lambda, 0.0);
    }

    complex<double>
    StandardModel::ckm_ub() const
    {
        double lambda2 = _lambda * _lambda, lambda3 = _lambda * lambda2, lambda4 = lambda2 * lambda3;
        double A2 = _A * _A;
        complex<double> num = complex<double>(_rhobar, -1.0 * _etabar);
        complex<double> denom = 1.0 - num * A2 * lambda4;
        double rest =  _A * lambda3 * sqrt((1.0 - A2 * lambda4) / (1.0 - lambda2));

        return (num / denom) * rest;
    }

    complex<double>
    StandardModel::ckm_ts() const
    {
        return complex<double>(-1.0 * _A * _lambda * _lambda, 0.0);
    }

    complex<double>
    StandardModel::ckm_tb() const
    {
        return complex<double>(1.0, 0.0);
    }
}
