/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/standard-model.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>

namespace eos
{
    using std::sqrt;

    template <> struct Implementation<StandardModel>
    {
        /* CKM Wolfenstein parameters */
        Parameter A;
        Parameter lambda;
        Parameter rhobar;
        Parameter etabar;
        // Measurement of V_cb is quite good, so use that
        Parameter ckm_cb;

        /* QCD parameters */
        Parameter alpha_s_Z;
        Parameter mu_t;
        Parameter mu_b;
        Parameter mu_c;

        /* Masses */
        Parameter m_b_MSbar;
        Parameter m_c_MSbar;
        Parameter m_Z;

        /* Energy Scale */
        Parameter mu;

        Implementation(const Parameters & p) :
            A(p["CKM::A"]),
            lambda(p["CKM::lambda"]),
            rhobar(p["CKM::rhobar"]),
            etabar(p["CKM::etabar"]),
            ckm_cb(p["CKM::|V_cb|"]),
            alpha_s_Z(p["QCD::alpha_s(MZ)"]),
            mu_t(p["QCD::mu_t"]),
            mu_b(p["QCD::mu_b"]),
            mu_c(p["QCD::mu_c"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c_MSbar(p["mass::c"]),
            m_Z(p["mass::Z"]),
            mu(p["mu"])
        {
        }

        // 5 flavor QCD
        static const QCD::BetaFunction qcd_nf5_beta_function;
        static const QCD::AnomalousMassDimension qcd_nf5_gamma_m;

        // 4 flavor QCD
        static const QCD::BetaFunction qcd_nf4_beta_function;
        static const QCD::AnomalousMassDimension qcd_nf4_gamma_m;
    };

    /* 5 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 5
    const QCD::BetaFunction Implementation<StandardModel>::qcd_nf5_beta_function
    {{
        23.0 / 3.0,
        116.0 / 3.0,
        9769.0 / 54.0,
        4826.1563287908967,
    }};

    // cf. [CKS2000], Eq. (7), p. 5 with n_f = 5
    const QCD::AnomalousMassDimension Implementation<StandardModel>::qcd_nf5_gamma_m
    {{
        1.0,
        506.0 / 9.0,
        474.87124557719461,
        2824.7862379694232,
    }};

    /* 4 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 4
    const QCD::BetaFunction Implementation<StandardModel>::qcd_nf4_beta_function
    {{
        25.0 / 3.0,
        154.0 / 3.0,
        21943.0 / 54.0,
        8035.1864197901160,
    }};

    // cf. [CKS2000], Eq. (7), p. 5 with n_f = 4
    const QCD::AnomalousMassDimension Implementation<StandardModel>::qcd_nf4_gamma_m
    {{
        1.0,
        526.0 / 9.0,
        636.61057670866927,
        6989.5510103599477,
    }};

    StandardModel::StandardModel(const Parameters & p) :
        PrivateImplementationPattern<StandardModel>(new Implementation<StandardModel>(p))
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
        if (mu > _imp->m_Z)
            throw InternalError("StandardModel::alpha_s: Running of alpha_s to mu > m_Z not yet implemented");

        double alpha_s_0 = _imp->alpha_s_Z, mu_0 = _imp->m_Z;

        if (mu >= _imp->mu_b)
            return QCD::alpha_s(mu, alpha_s_0, mu_0, Implementation<StandardModel>::qcd_nf5_beta_function);

        alpha_s_0 = QCD::alpha_s(_imp->mu_b, alpha_s_0, mu_0, Implementation<StandardModel>::qcd_nf5_beta_function);
        mu_0 = _imp->mu_b;

        if (mu >= _imp->mu_c)
            return QCD::alpha_s(mu, alpha_s_0, mu_0, Implementation<StandardModel>::qcd_nf4_beta_function);

        throw InternalError("StandardModel::alpha_s: Running of alpha_s to mu < mu_c not yet implemented");
    }

    double
    StandardModel::m_b_msbar(const double & mu) const
    {
        double m_b_MSbar = _imp->m_b_MSbar();
        double alpha_mu_0 = alpha_s(m_b_MSbar);

        if (mu > _imp->m_b_MSbar)
        {
            if (mu < _imp->mu_t)
                return QCD::m_q_msbar(m_b_MSbar, alpha_mu_0, alpha_s(mu), Implementation<StandardModel>::qcd_nf5_beta_function, Implementation<StandardModel>::qcd_nf5_gamma_m);

            throw InternalError("StandardModel::m_b_msbar: Running of m_b_MSbar to mu > mu_t not yet implemented");
        }
        else
        {
            if (mu >= _imp->mu_c)
                return QCD::m_q_msbar(m_b_MSbar, alpha_mu_0, alpha_s(mu), Implementation<StandardModel>::qcd_nf4_beta_function, Implementation<StandardModel>::qcd_nf4_gamma_m);

            throw InternalError("StandardModel::m_b_msbar: Running of m_b_MSbar to mu < mu_c not yet implemented");
        }
    }

    double
    StandardModel::m_b_pole() const
    {
        double m_b_MSbar = _imp->m_b_MSbar();

        return QCD::m_q_pole(m_b_MSbar, alpha_s(m_b_MSbar), 5.0);
    }

    double
    StandardModel::m_b_ps(const double & mu_f) const
    {
        double m_b_MSbar = _imp->m_b_MSbar();

        return QCD::m_q_ps(m_b_MSbar, alpha_s(m_b_MSbar), mu_f, 5.0, Implementation<StandardModel>::qcd_nf5_beta_function);
    }

    /* Charm */
    double
    StandardModel::m_c_msbar(const double & mu) const
    {
        double m_c_0 = _imp->m_c_MSbar();
        double alpha_s_mu0 = alpha_s(m_c_0);

        if (mu >= _imp->m_c_MSbar)
        {
            if (mu <= _imp->mu_b)
                return QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s(mu), Implementation<StandardModel>::qcd_nf4_beta_function, Implementation<StandardModel>::qcd_nf4_gamma_m);

            double alpha_s_b = alpha_s(_imp->mu_b);
            m_c_0 = QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s_b, Implementation<StandardModel>::qcd_nf4_beta_function, Implementation<StandardModel>::qcd_nf4_gamma_m);
            alpha_s_mu0 = alpha_s_b;

            if (mu <= _imp->mu_t)
                return QCD::m_q_msbar(m_c_0, alpha_s_mu0, alpha_s(mu), Implementation<StandardModel>::qcd_nf5_beta_function, Implementation<StandardModel>::qcd_nf5_gamma_m);

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
        double m_c_MSbar = _imp->m_c_MSbar();

        return QCD::m_q_pole(m_c_MSbar, alpha_s(m_c_MSbar), 4.0);
    }

    complex<double>
    StandardModel::ckm_cb() const
    {
        return complex<double>(_imp->ckm_cb, 0.0);
    }

    complex<double>
    StandardModel::ckm_us() const
    {
        return complex<double>(_imp->lambda, 0.0);
    }

    complex<double>
    StandardModel::ckm_ub() const
    {
        double lambda2 = _imp->lambda * _imp->lambda, lambda3 = _imp->lambda * lambda2, lambda4 = lambda2 * lambda3;
        double A2 = _imp->A * _imp->A;
        complex<double> num = complex<double>(_imp->rhobar, -1.0 * _imp->etabar);
        complex<double> denom = 1.0 - num * A2 * lambda4;
        double rest =  _imp->A * lambda3 * sqrt((1.0 - A2 * lambda4) / (1.0 - lambda2));

        return (num / denom) * rest;
    }

    complex<double>
    StandardModel::ckm_ts() const
    {
        return complex<double>(-1.0 * _imp->A * _imp->lambda * _imp->lambda, 0.0);
    }

    complex<double>
    StandardModel::ckm_tb() const
    {
        return complex<double>(1.0, 0.0);
    }
}
