/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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

#include <eos/form-factors/baryonic.hh>
#include <eos/rare-b-decays/charm-loops.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-dilepton.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/integrate-impl.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <functional>
#include <map>
#include <string>

namespace eos
{
    using std::norm;

    namespace lambdab_to_lambda_dilepton
    {
        struct Amplitudes
        {
            complex<double> a_perp_0_L, a_perp_0_R;
            complex<double> a_para_0_L, a_para_0_R;
            complex<double> a_perp_1_L, a_perp_1_R;
            complex<double> a_para_1_L, a_para_1_R;
            double alpha;
        };

        struct AngularObservables
        {
            std::array<double, 10> _k;

            AngularObservables(const Amplitudes & a)
            {
                // cf. [BFvD2014], eqs. (3.29) - (3.32), pp. 8-9
                _k[0] = (std::norm(a.a_perp_1_R) + std::norm(a.a_para_1_R) + 2.0 * std::norm(a.a_perp_0_R) + 2.0 * std::norm(a.a_para_0_R)
                      +  std::norm(a.a_perp_1_L) + std::norm(a.a_para_1_L) + 2.0 * std::norm(a.a_perp_0_L) + 2.0 * std::norm(a.a_para_0_L)) / 4.0;
                _k[1] = (std::norm(a.a_perp_1_R) + std::norm(a.a_para_1_R)
                      +  std::norm(a.a_perp_1_L) + std::norm(a.a_para_1_L)) / 2.0;
                _k[2] = -std::real(a.a_perp_1_R * std::conj(a.a_para_1_R) - a.a_perp_1_L * std::conj(a.a_para_1_L));
                _k[3] = (std::real(a.a_perp_1_R * std::conj(a.a_para_1_R) + 2.0 * a.a_perp_0_R * std::conj(a.a_para_0_R)
                      +            a.a_perp_1_L * std::conj(a.a_para_1_L) + 2.0 * a.a_perp_0_L * std::conj(a.a_para_0_L))) * a.alpha / 2.0;
                _k[4] = std::real(a.a_perp_1_R * std::conj(a.a_para_1_R)
                      +           a.a_perp_1_L * std::conj(a.a_para_1_L)) * a.alpha;
                _k[5] = ((std::norm(a.a_perp_1_R) + std::norm(a.a_para_1_R))
                      -  (std::norm(a.a_perp_1_L) + std::norm(a.a_para_1_L))) * a.alpha * -0.5;
                _k[6] = std::imag(a.a_perp_1_R * std::conj(a.a_perp_0_R) - a.a_para_1_R * std::conj(a.a_para_0_R)
                      +           a.a_perp_1_L * std::conj(a.a_perp_0_L) - a.a_para_1_L * std::conj(a.a_para_0_L)) * a.alpha / std::sqrt(2.0);
                _k[7] = std::imag(a.a_perp_1_R * std::conj(a.a_para_0_R) - a.a_para_1_R * std::conj(a.a_perp_0_R)
                      -          (a.a_perp_1_L * std::conj(a.a_para_0_L) - a.a_para_1_L * std::conj(a.a_perp_0_L))) * a.alpha / std::sqrt(2.0);
                _k[8] = std::real(a.a_perp_1_R * std::conj(a.a_para_0_R) - a.a_para_1_R * std::conj(a.a_perp_0_R)
                      +           a.a_perp_1_L * std::conj(a.a_para_0_L) - a.a_para_1_L * std::conj(a.a_perp_0_L)) * a.alpha / std::sqrt(2.0);
                _k[9] = std::real(a.a_perp_1_R * std::conj(a.a_perp_0_R) - a.a_para_1_R * std::conj(a.a_para_0_R)
                      -          (a.a_perp_1_L * std::conj(a.a_perp_0_L) - a.a_para_1_L * std::conj(a.a_para_0_L))) * a.alpha / std::sqrt(2.0);
            }

            AngularObservables(const std::array<double, 10> & k) :
                _k(k)
            {
            }

            inline double k1ss() const { return _k[0]; }
            inline double k1cc() const { return _k[1]; }
            inline double k1c()  const { return _k[2]; }
            inline double k2ss() const { return _k[3]; }
            inline double k2cc() const { return _k[4]; }
            inline double k2c()  const { return _k[5]; }
            inline double k3sc() const { return _k[6]; }
            inline double k3s()  const { return _k[7]; }
            inline double k4sc() const { return _k[8]; }
            inline double k4s()  const { return _k[9]; }

            inline double decay_width() const
            {
                return 2.0 * k1ss() + k1cc();
            }

            inline double a_fb_leptonic() const
            {
                return 3.0 / 2.0 * k1c() / decay_width();
            }

            inline double a_fb_hadronic() const
            {
                return 1.0 / 2.0 * (2.0 * k2ss() + k2cc()) / decay_width();
            }

            inline double a_fb_combined() const
            {
                return 3.0 / 4.0 * k2c() / decay_width();
            }

            inline double f_zero() const
            {
                return (2.0 * k1ss() - k1cc()) / decay_width();
            }
        };
    }

    /* Large Recoil */

    template <> struct Implementation<LambdaBToLambdaDilepton<LargeRecoil>>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;
        UsedParameter tau_Lambda_b;

        UsedParameter g_fermi;

        UsedParameter m_Lambda_b;
        UsedParameter m_Lambda;
        UsedParameter alpha;

        UsedParameter alpha_e;
        UsedParameter mu;

        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["hbar"], u),
            tau_Lambda_b(p["life_time::Lambda_b"], u),
            g_fermi(p["G_Fermi"], u),
            m_Lambda_b(p["mass::Lambda_b"], u),
            m_Lambda(p["mass::Lambda"], u),
            alpha(p["Lambda::alpha"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            mu(p["mu"], u)
        {
            form_factors = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda@" + o.get("form-factors", "BFvD2014"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);

        }

        double norm(const double & s) const
        {
            return g_fermi() * alpha_e() * abs(model->ckm_tb() * conj(model->ckm_ts()))
                * std::sqrt(s / 3.0 / 2048 / std::pow(M_PI, 5.0) / power_of<3>(m_Lambda_b())
                * std::sqrt(lambda(m_Lambda_b * m_Lambda_b, m_Lambda * m_Lambda, s))); // cf. [BFvD2014], Eq. (?), p. ??
        }

        double kappa() const
        {
            return (1.0 - 2.0 * model->alpha_s(mu) / (3.0 * M_PI) * std::log(mu / model->m_b_msbar(mu)));
        }

        lambdab_to_lambda_dilepton::Amplitudes amplitudes(const double & s)
        {
            lambdab_to_lambda_dilepton::Amplitudes result;

            double alpha_s = model->alpha_s(mu());
            double m_b_MSbar = model->m_b_msbar(mu), m_b_PS = model->m_b_ps(2.0), m_b_PS2 = m_b_PS * m_b_PS;
            double m_c_pole = model->m_c_pole();

            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu());

            complex<double> lambda_hat_u = model->ckm_ub() * conj(model->ckm_us()) / std::abs(model->ckm_tb() * conj(model->ckm_ts()));
            double sqrtsminus = std::sqrt(power_of<2>(m_Lambda_b - m_Lambda) - s), sqrtsplus = std::sqrt(power_of<2>(m_Lambda_b + m_Lambda) - s), sqrts = std::sqrt(s);
            double N = norm(s);

            /* Y(s) for the up and the top sector */
            // cf. [BFS2001], Eq. (10), p. 4
            complex<double> Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2() + 6.0 * wc.c3() + 60.0 * wc.c5();
            complex<double> Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
            complex<double> Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64 / 3.0 * wc.c6());
            complex<double> Y_top_ = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());
            // Use b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses,
            // then replace b pole mass by the PS mass.
            complex<double> Y_top = Y_top_c * CharmLoops::h(mu, s, m_c_pole)
                 + Y_top_b * CharmLoops::h(mu, s, m_b_PS)
                 + Y_top_0 * CharmLoops::h(mu, s)
                 + Y_top_;
            // cf. [BFS2004], Eq. (43), p. 24
            complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (CharmLoops::h(mu, s, m_c_pole) - CharmLoops::h(mu, s));

            // calculate effective wilson coefficients
            // cf. [BFS2001], below Eq. (9), p. 4
            complex<double> c7eff = wc.c7() - 1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6();
            // cf. [BFS2001], below Eq. (26), p. 8
            complex<double> c8eff = wc.c8() + wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

            // two loop virtual corrections, cf. [AAGW2001]
            // charm quarks
            complex<double> F27c = CharmLoops::F27_massive(mu(), s, m_b_PS, m_c_pole);
            complex<double> F17c = -F27c / 6.0;
            complex<double> F19c = CharmLoops::F19_massive(mu(), s, m_b_PS, m_c_pole);
            complex<double> F29c = CharmLoops::F29_massive(mu(), s, m_b_PS, m_c_pole);
            // up quarks
            complex<double> F27u = CharmLoops::F27_massless(mu(), s, m_b_PS);
            complex<double> F17u = -F27u / 6.0;
            complex<double> F19u = CharmLoops::F19_massless(mu(), s, m_b_PS);
            complex<double> F29u = CharmLoops::F29_massless(mu(), s, m_b_PS);
            // gluon
            complex<double> F87  = CharmLoops::F87_massless(mu(), s, m_b_PS);
            complex<double> F89  = CharmLoops::F89_massless(s, m_b_PS);

            // integredients for form factor relations
            // cf. [FY2011]
            double L = -1.0 * (m_b_PS2 - s) / s * std::log(1.0 - s / m_b_PS2);

            // ratio of tensor to vector form factors
            // cf. [BFvD2014], eqs. (??)-(??)
            double R1p = 1.0 + alpha_s / (3.0 * M_PI) * (2.0 * std::log(m_b_PS / mu()) - 2.0 - L);
            double R1m = 1.0 + alpha_s / (3.0 * M_PI) * (2.0 * std::log(m_b_PS / mu()) - 2.0 - L);
            double R0p = 1.0 + alpha_s / (3.0 * M_PI) * (2.0 * std::log(m_b_PS / mu()) - 2.0 + 2.0 * L);
            double R0m = 1.0 + alpha_s / (3.0 * M_PI) * (2.0 * std::log(m_b_PS / mu()) - 2.0 + 2.0 * L);

            // functions tau
            // cf. [BFvD2014], eqs. (??)-(??)
            complex<double> tau_1p = (m_Lambda_b + m_Lambda) / m_Lambda_b * (
                    c7eff + wc.c7prime()
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F17c + wc.c2() * F27c + c8eff * F87)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F17c - F17u) + wc.c2() * (F27c - F27u) + c8eff * F87) * lambda_hat_u
                ) * R1p
                + s / (2.0 * m_b_MSbar * m_Lambda_b) * (
                    Y_top + lambda_hat_u * Y_up
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F19c + wc.c2() * F29c + wc.c8() * F89)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F19c - F19u) + wc.c2() * (F29c - F29u) + wc.c8() * F89) * lambda_hat_u
                );
            complex<double> tau_1m = (m_Lambda_b - m_Lambda) / m_Lambda_b * (
                    c7eff - wc.c7prime()
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F17c + wc.c2() * F27c)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F17c - F17u) + wc.c2() * (F27c - F27u) + c8eff * F87) * lambda_hat_u
                ) * R1m
                + s / (2.0 * m_b_MSbar * m_Lambda_b) * (
                    Y_top + lambda_hat_u * Y_up
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F19c + wc.c2() * F29c + c8eff * F87)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F19c - F19u) + wc.c2() * (F29c - F29u) + wc.c8() * F89) * lambda_hat_u
                );
            complex<double> tau_0p = m_Lambda_b / (m_Lambda_b + m_Lambda) * (
                    c7eff + wc.c7prime()
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F17c + wc.c2() * F27c)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F17c - F17u) + wc.c2() * (F27c - F27u) + c8eff * F87) * lambda_hat_u
                ) * R0p
                + m_Lambda_b / (2.0 * m_b_MSbar) * (
                    Y_top + lambda_hat_u * Y_up
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F19c + wc.c2() * F29c)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F19c - F19u) + wc.c2() * (F29c - F29u) + wc.c8() * F89) * lambda_hat_u
                );
            complex<double> tau_0m = m_Lambda_b / (m_Lambda_b - m_Lambda) * (
                    c7eff - wc.c7prime()
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F17c + wc.c2() * F27c)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F17c - F17u) + wc.c2() * (F27c - F27u) + c8eff * F87) * lambda_hat_u
                ) * R0m
                + m_Lambda_b / (2.0 * m_b_MSbar) * (
                    Y_top + lambda_hat_u * Y_up
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F19c + wc.c2() * F29c)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F19c - F19u) + wc.c2() * (F29c - F29u) + wc.c8() * F89) * lambda_hat_u
                );

            // cf. [BFvD2014], eqs. (??)-(??)
            result.a_perp_1_R = -2.0 *       N * (wc.c9() + wc.c9prime() + (wc.c10() + wc.c10prime()) + 2.0 * m_b_MSbar * m_Lambda_b / s * tau_1p) * form_factors->f_perp_v(s) * sqrtsminus;
            result.a_perp_1_L = -2.0 *       N * (wc.c9() + wc.c9prime() - (wc.c10() + wc.c10prime()) + 2.0 * m_b_MSbar * m_Lambda_b / s * tau_1p) * form_factors->f_perp_v(s) * sqrtsminus;

            result.a_para_1_R = +2.0 *       N * (wc.c9() - wc.c9prime() + (wc.c10() - wc.c10prime()) + 2.0 * m_b_MSbar * m_Lambda_b / s * tau_1m) * form_factors->f_perp_a(s) * sqrtsplus;
            result.a_para_1_L = +2.0 *       N * (wc.c9() - wc.c9prime() - (wc.c10() - wc.c10prime()) + 2.0 * m_b_MSbar * m_Lambda_b / s * tau_1m) * form_factors->f_perp_a(s) * sqrtsplus;

            result.a_perp_0_R = +sqrt(2.0) * N * (wc.c9() + wc.c9prime() + (wc.c10() + wc.c10prime()) + 2.0 * m_b_MSbar / m_Lambda_b * tau_0p) * form_factors->f_long_v(s) * (m_Lambda_b + m_Lambda)  / sqrts * sqrtsminus;
            result.a_perp_0_L = +sqrt(2.0) * N * (wc.c9() + wc.c9prime() - (wc.c10() + wc.c10prime()) + 2.0 * m_b_MSbar / m_Lambda_b * tau_0p) * form_factors->f_long_v(s) * (m_Lambda_b + m_Lambda)  / sqrts * sqrtsminus;

            result.a_para_0_R = -sqrt(2.0) * N * (wc.c9() - wc.c9prime() + (wc.c10() - wc.c10prime()) + 2.0 * m_b_MSbar / m_Lambda_b * tau_0m) * form_factors->f_long_a(s) * (m_Lambda_b - m_Lambda)  / sqrts * sqrtsplus;
            result.a_para_0_L = -sqrt(2.0) * N * (wc.c9() - wc.c9prime() - (wc.c10() - wc.c10prime()) + 2.0 * m_b_MSbar / m_Lambda_b * tau_0m) * form_factors->f_long_a(s) * (m_Lambda_b - m_Lambda)  / sqrts * sqrtsplus;

            result.alpha = this->alpha();

            return result;
        }

        std::array<double, 10> _differential_angular_observables(const double & s)
        {
            return lambdab_to_lambda_dilepton::AngularObservables(this->amplitudes(s))._k;
        }

        std::array<double, 10> _integrated_angular_observables(const double & s_min, const double & s_max)
        {
            std::function<std::array<double, 10> (const double &)> integrand(std::bind(&Implementation::_differential_angular_observables, this, std::placeholders::_1));

            return integrate(integrand, 32, s_min, s_max);
        }

        inline lambdab_to_lambda_dilepton::AngularObservables differential_angular_observables(const double & s)
        {
            return lambdab_to_lambda_dilepton::AngularObservables{ _differential_angular_observables(s) };
        }

        inline lambdab_to_lambda_dilepton::AngularObservables integrated_angular_observables(const double & s_min, const double & s_max)
        {
            return lambdab_to_lambda_dilepton::AngularObservables{ _integrated_angular_observables(s_min, s_max) };
        }
    };

    LambdaBToLambdaDilepton<LargeRecoil>::LambdaBToLambdaDilepton(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<LambdaBToLambdaDilepton<LargeRecoil>>(new Implementation<LambdaBToLambdaDilepton<LargeRecoil>>(p, o, *this))
    {
    }

    LambdaBToLambdaDilepton<LargeRecoil>::~LambdaBToLambdaDilepton()
    {
    }

    /* q^2-differential observables */
    double
    LambdaBToLambdaDilepton<LargeRecoil>::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_angular_observables(s).decay_width() * _imp->tau_Lambda_b / _imp->hbar;
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_leptonic();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_hadronic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_hadronic();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_combined(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_combined();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::differential_fzero(const double & s) const
    {
        return _imp->differential_angular_observables(s).f_zero();
    }

    /* q^2-integrated observables */
    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).decay_width() * _imp->tau_Lambda_b / _imp->hbar;
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_hadronic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_hadronic();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_combined(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_combined();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_fzero(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_zero();
    }

    /* Low Recoil */

    namespace lambdab_to_lambda_dilepton
    {
        // Handle ratio of tensor to vector form factors
        // model-independently.
        class TensorToVectorRatios :
            public ParameterUser
        {
            public:
                virtual ~TensorToVectorRatios()
                {
                }

                virtual double R_perp_V(const double & s) const = 0;
                virtual double R_perp_A(const double & s) const = 0;
                virtual double R_long_V(const double & s) const = 0;
                virtual double R_long_A(const double & s) const = 0;

                static TensorToVectorRatios * make(const std::string & name, const Parameters & p, const Options & o);
        };

        // HQET estimates based on [BFvD2014]
        class EstimatedTensorToVectorRatios :
            public TensorToVectorRatios
        {
            private:
                UsedParameter _sigma_perp_V;
                UsedParameter _sigma_perp_A;
                UsedParameter _sigma_long_V;
                UsedParameter _sigma_long_A;

            public:
                EstimatedTensorToVectorRatios(const Parameters & p, const Options &) :
                    _sigma_perp_V(p["Lambda_b->Lambda::sigma_perp^V@BFvD2014"], *this),
                    _sigma_perp_A(p["Lambda_b->Lambda::sigma_perp^A@BFvD2014"], *this),
                    _sigma_long_V(p["Lambda_b->Lambda::sigma_long^V@BFvD2014"], *this),
                    _sigma_long_A(p["Lambda_b->Lambda::sigma_long^A@BFvD2014"], *this)
                {
                }

                ~EstimatedTensorToVectorRatios()
                {
                }

                virtual double R_perp_V(const double & s) const
                {
                    // cf. [BFvD2014], eq. (C.2) and table 4.
                    return 1.205 - 1.00e-2 * s + _sigma_perp_V() * (5.01 - 0.195 * s);
                }

                virtual double R_perp_A(const double & s) const
                {
                    // cf. [BFvD2014], eq. (C.2) and table 4.
                    return 1.342 - 1.57e-2 * s + _sigma_perp_A() * (4.51 - 0.171 * s);
                }

                virtual double R_long_V(const double & s) const
                {
                    // cf. [BFvD2014], eq. (C.2) and table 4.
                    return 1.235 - 1.00e-2 * s + _sigma_long_V() * (5.74 - 0.231 * s);
                }

                virtual double R_long_A(const double & s) const
                {
                    // cf. [BFvD2014], eq. (C.2) and table 4.
                    return 1.294 - 1.29e-2 * s + _sigma_long_A() * (4.61 - 0.176 * s);
                }

                static EstimatedTensorToVectorRatios * make(const Parameters & p, const Options & o)
                {
                    return new EstimatedTensorToVectorRatios(p, o);
                }
        };

        TensorToVectorRatios *
        TensorToVectorRatios::make(const std::string & name, const Parameters & p, const Options & o)
        {
            typedef std::function<TensorToVectorRatios * (const Parameters &, const Options &)> RatiosMaker;
            static const std::map<std::string, RatiosMaker> ratios_makers
            {
                std::make_pair("Estimate",   &EstimatedTensorToVectorRatios::make)
            };

            auto i = ratios_makers.find(name);

            if (ratios_makers.cend() == i)
                throw InternalError("No such ratio '" + name + "'");

            return i->second(p, o);
        }
    }

    template <> struct Implementation<LambdaBToLambdaDilepton<LowRecoil>>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;
        UsedParameter tau_Lambda_b;

        UsedParameter g_fermi;

        UsedParameter m_Lambda_b;
        UsedParameter m_Lambda;
        UsedParameter alpha;

        UsedParameter alpha_e;
        UsedParameter mu;

        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

        std::shared_ptr<lambdab_to_lambda_dilepton::TensorToVectorRatios> ratios;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["hbar"], u),
            tau_Lambda_b(p["life_time::Lambda_b"], u),
            g_fermi(p["G_Fermi"], u),
            m_Lambda_b(p["mass::Lambda_b"], u),
            m_Lambda(p["mass::Lambda"], u),
            alpha(p["Lambda::alpha"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            mu(p["mu"], u)
        {
            form_factors = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda@" + o.get("form-factors", "BFvD2014"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            using namespace lambdab_to_lambda_dilepton;
            ratios = std::shared_ptr<TensorToVectorRatios>(TensorToVectorRatios::make(o.get("ratios", "Estimate"), p, o));

            u.uses(*ratios);
            u.uses(*form_factors);
            u.uses(*model);
        }

        double norm(const double & s) const
        {
            return g_fermi() * alpha_e() * abs(model->ckm_tb() * conj(model->ckm_ts()))
                * std::sqrt(s / 3.0 / 2048 / std::pow(M_PI, 5.0) / power_of<3>(m_Lambda_b())
                * std::sqrt(lambda(m_Lambda_b * m_Lambda_b, m_Lambda * m_Lambda, s))); // cf. [BFvD2014], Eq. (3.18), p. 6
        }

        double kappa() const
        {
            return (1.0 - 2.0 * model->alpha_s(mu) / (3.0 * M_PI) * std::log(mu / model->m_b_msbar(mu)));
        }

        lambdab_to_lambda_dilepton::Amplitudes amplitudes(const double & s)
        {
            lambdab_to_lambda_dilepton::Amplitudes result;

            double alpha_s = model->alpha_s(mu()), m_b = model->m_b_ps(2.0), m_c = model->m_c_msbar(mu());
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu());
            complex<double> lambda_hat_u = model->ckm_ub() * conj(model->ckm_us()) / std::abs(model->ckm_tb() * conj(model->ckm_ts()));
            double sqrtsminus = std::sqrt(power_of<2>(m_Lambda_b - m_Lambda) - s), sqrtsplus = std::sqrt(power_of<2>(m_Lambda_b + m_Lambda) - s), sqrts = std::sqrt(s);
            double N = norm(s), kappa = this->kappa();

            // calculate effective wilson coefficients
            complex<double> c7eff = ShortDistanceLowRecoil::c7eff(s, mu(), alpha_s, m_b, true, wc);
            complex<double> c9eff = ShortDistanceLowRecoil::c9eff(s, mu(), alpha_s, m_b, m_c, true, false, lambda_hat_u, wc);

            // cf. [BFvD2014], eq.s (??), p. ??
            double zeta_perp_V = (m_Lambda_b + m_Lambda) / m_Lambda_b * (ratios->R_perp_V(s) - 1.0) + m_Lambda / m_Lambda_b;
            double zeta_perp_A = (m_Lambda_b - m_Lambda) / m_Lambda_b * (ratios->R_perp_A(s) - 1.0) - m_Lambda / m_Lambda_b;
            double zeta_long_V = s / ((m_Lambda_b + m_Lambda) * m_Lambda_b) * ratios->R_long_V(s) - 1.0;
            double zeta_long_A = s / ((m_Lambda_b - m_Lambda) * m_Lambda_b) * ratios->R_long_A(s) - 1.0;

            // cf. [BFvD2014], eqs. (4.9)-(4.10), p. 11
            result.a_perp_1_R = -2.0 *       N * (c9eff + wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff + wc.c7prime()) * zeta_perp_V + (wc.c10() + wc.c10prime())) * form_factors->f_perp_v(s) * sqrtsminus;
            result.a_perp_1_L = -2.0 *       N * (c9eff + wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff + wc.c7prime()) * zeta_perp_V - (wc.c10() + wc.c10prime())) * form_factors->f_perp_v(s) * sqrtsminus;

            result.a_para_1_R = +2.0 *       N * (c9eff - wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff - wc.c7prime()) * zeta_perp_A + (wc.c10() - wc.c10prime())) * form_factors->f_perp_a(s) * sqrtsplus;
            result.a_para_1_L = +2.0 *       N * (c9eff - wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff - wc.c7prime()) * zeta_perp_A - (wc.c10() - wc.c10prime())) * form_factors->f_perp_a(s) * sqrtsplus;

            result.a_perp_0_R = +sqrt(2.0) * N * (c9eff + wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff + wc.c7prime()) * zeta_long_V + (wc.c10() + wc.c10prime())) * form_factors->f_long_v(s)
                * (m_Lambda_b + m_Lambda)  / sqrts * sqrtsminus;
            result.a_perp_0_L = +sqrt(2.0) * N * (c9eff + wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff + wc.c7prime()) * zeta_long_V - (wc.c10() + wc.c10prime())) * form_factors->f_long_v(s)
                * (m_Lambda_b + m_Lambda)  / sqrts * sqrtsminus;

            result.a_para_0_R = -sqrt(2.0) * N * (c9eff - wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff - wc.c7prime()) * zeta_long_A + (wc.c10() - wc.c10prime())) * form_factors->f_long_a(s)
                * (m_Lambda_b - m_Lambda)  / sqrts * sqrtsplus;
            result.a_para_0_L = -sqrt(2.0) * N * (c9eff - wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff - wc.c7prime()) * zeta_long_A - (wc.c10() - wc.c10prime())) * form_factors->f_long_a(s)
                * (m_Lambda_b - m_Lambda)  / sqrts * sqrtsplus;

            result.alpha = this->alpha();

            return result;
        }

        std::array<double, 10> _differential_angular_observables(const double & s)
        {
            return lambdab_to_lambda_dilepton::AngularObservables(this->amplitudes(s))._k;
        }

        std::array<double, 10> _integrated_angular_observables(const double & s_min, const double & s_max)
        {
            std::function<std::array<double, 10> (const double &)> integrand(std::bind(&Implementation::_differential_angular_observables, this, std::placeholders::_1));

            return integrate(integrand, 32, s_min, s_max);
        }

        inline lambdab_to_lambda_dilepton::AngularObservables differential_angular_observables(const double & s)
        {
            return lambdab_to_lambda_dilepton::AngularObservables{ _differential_angular_observables(s) };
        }

        inline lambdab_to_lambda_dilepton::AngularObservables integrated_angular_observables(const double & s_min, const double & s_max)
        {
            return lambdab_to_lambda_dilepton::AngularObservables{ _integrated_angular_observables(s_min, s_max) };
        }
    };

    LambdaBToLambdaDilepton<LowRecoil>::LambdaBToLambdaDilepton(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<LambdaBToLambdaDilepton<LowRecoil>>(new Implementation<LambdaBToLambdaDilepton<LowRecoil>>(p, o, *this))
    {
    }

    LambdaBToLambdaDilepton<LowRecoil>::~LambdaBToLambdaDilepton()
    {
    }

    /* q^2-differential observables */
    double
    LambdaBToLambdaDilepton<LowRecoil>::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_angular_observables(s).decay_width() * _imp->tau_Lambda_b / _imp->hbar;
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_leptonic();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_hadronic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_hadronic();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_combined(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_combined();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::differential_fzero(const double & s) const
    {
        return _imp->differential_angular_observables(s).f_zero();
    }

    /* q^2-integrated observables */
    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).decay_width() * _imp->tau_Lambda_b / _imp->hbar;
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_hadronic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_hadronic();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_combined(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_combined();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_fzero(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_zero();
    }
}
