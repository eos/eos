/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019      Ahmet Kokulu
 * Copyright (c) 2019-2025 Danny van Dyk
 * Copyright (c) 2023      Méril Reboud
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

#include <eos/c-decays/lambdac-to-lambda-l-nu.hh>
#include <eos/form-factors/baryonic.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <functional>
#include <map>
#include <string>

namespace eos
{
    using std::norm;

    namespace lambdac_to_lambda_l_nu
    {
        struct Amplitudes
        {
            complex<double> perp_0_L;
            complex<double> para_0_L;
            complex<double> perp_1_L;
            complex<double> para_1_L;
            complex<double> perp_t_L;
            complex<double> para_t_L;

            complex<double> perp_0_T;
            complex<double> para_0_T;
            complex<double> perp_1_T;
            complex<double> para_1_T;

            double alpha;
            double beta;
        };

        struct AngularObservables
        {
            std::array<double, 10> _k;

            AngularObservables(const Amplitudes & a)
            {
                using std::norm;
                using std::real;
                using std::imag;
                using std::conj;
                using std::sqrt;

                double beta       = a.beta;
                double sqrt1mbeta = sqrt(1.0 - beta);

                // cf. [BKTvD2019], eqs. (2.5)-(2.14), pp. 4-5.

                // K_{1ss}
                _k[0] = (
                            + 2.0                * (norm(a.para_0_L) + norm(a.perp_0_L))
                            + (2.0 - beta)       * (norm(a.para_1_L) + norm(a.perp_1_L) + norm(a.para_1_T) + norm(a.perp_1_T))
                            + 2.0 * (1.0 - beta) * (norm(a.para_t_L) + norm(a.perp_t_L) + norm(a.para_0_T) + norm(a.perp_0_T))
                            - 4.0 * sqrt1mbeta   * real(
                                  a.para_0_T * conj(a.para_0_L) + a.perp_0_T * conj(a.perp_0_L)
                                + a.para_1_T * conj(a.para_1_L) + a.perp_1_T * conj(a.perp_1_L)
                            )
                        ) / 4.0;

                // K_{1cc}
                _k[1] = (
                            +                    (norm(a.para_1_L) + norm(a.perp_1_L) + norm(a.para_0_T) + norm(a.perp_0_T))
                            + (1.0 - beta)     * (norm(a.para_0_L) + norm(a.perp_0_L) + norm(a.para_t_L) + norm(a.perp_t_L) + norm(a.para_1_T) + norm(a.perp_1_T))
                            - 2.0 * sqrt1mbeta * real(
                                  a.para_0_T * conj(a.para_0_L) + a.perp_0_T * conj(a.perp_0_L)
                                + a.para_1_T * conj(a.para_1_L) + a.perp_1_T * conj(a.perp_1_L)
                            )
                        ) / 2.0;

                // K_{1c}
                _k[2] = real(
                            +                 a.perp_1_L * conj(a.para_1_L)
                            + (1.0 - beta) * (a.para_0_L * conj(a.para_t_L) + a.perp_0_L * conj(a.perp_t_L) + a.perp_1_T * conj(a.para_1_T))
                            - sqrt1mbeta   * (
                                  a.perp_1_T * conj(a.para_1_L) + a.para_1_T * conj(a.perp_1_L)
                                + a.para_0_T * conj(a.para_t_L) + a.perp_0_T * conj(a.perp_t_L)
                            )
                        );

                // K_{2ss}
                _k[3] = a.alpha * real(
                            + 2.0                *  a.perp_0_L * conj(a.para_0_L)
                            + (2.0 - beta)       * (a.perp_1_L * conj(a.para_1_L) - a.perp_1_T * conj(a.para_1_T))
                            + 2.0 * (1.0 - beta) * (a.perp_t_L * conj(a.para_t_L) - a.perp_0_T * conj(a.para_0_T))
                            - 2.0 * sqrt1mbeta   * (
                                  a.perp_0_T * conj(a.para_0_L) + a.para_0_T * conj(a.perp_0_L)
                                + a.perp_1_T * conj(a.para_1_L) + a.para_1_T * conj(a.perp_1_L)
                            )
                        ) / 2.0;

                // K_{2cc}
                _k[4] = a.alpha * real(
                            +                (a.perp_1_L * conj(a.para_1_L) + a.perp_0_T * conj(a.para_0_T))
                            + (1.0 - beta) * (a.perp_0_L * conj(a.para_0_L) + a.perp_t_L * conj(a.para_t_L) + a.perp_1_T * conj(a.para_1_T))
                            - sqrt1mbeta   * (
                                  a.perp_0_T * conj(a.para_0_L) + a.para_0_T * conj(a.perp_0_L)
                                + a.perp_1_T * conj(a.para_1_L) + a.para_1_T * conj(a.perp_1_L)
                            )
                        );

                // K_{2c}
                _k[5] = a.alpha * (
                            +                      (norm(a.para_1_L) + norm(a.perp_1_L))
                            + (1.0 - beta)       * (norm(a.para_1_T) + norm(a.perp_1_T))
                            + 2.0 * (1.0 - beta) * real(
                                a.perp_0_L * conj(a.para_t_L) + a.para_0_L * conj(a.perp_t_L)
                            )
                            - 2.0 * sqrt1mbeta   * real(
                                  a.para_1_T * conj(a.para_1_L) + a.perp_1_T * conj(a.perp_1_L)
                                + a.perp_0_T * conj(a.para_t_L) + a.para_0_T * conj(a.perp_t_L)
                            )
                        ) / 2.0;

                // K_{3sc}
                _k[6] = a.alpha * beta * imag(
                            + a.perp_1_L * conj(a.perp_0_L) - a.para_1_L * conj(a.para_0_L)
                            + a.para_1_T * conj(a.para_0_T) - a.perp_1_T * conj(a.perp_0_T)
                        ) / sqrt(2.0);

                // K_{3s}
                _k[7] = a.alpha * imag(
                            +                (a.para_1_L * conj(a.perp_0_L) - a.perp_1_L * conj(a.para_0_L))
                            + (1.0 - beta) * (
                                + a.para_1_L * conj(a.para_t_L) - a.perp_1_L * conj(a.perp_t_L)
                                + a.para_1_T * conj(a.perp_0_T) - a.perp_1_T * conj(a.para_0_T)
                            )
                            + sqrt1mbeta   * (
                                + a.perp_0_T * conj(a.para_1_L) + a.perp_1_T * conj(a.para_0_L) + a.perp_1_T * conj(a.perp_t_L)
                                - a.para_0_T * conj(a.perp_1_L) - a.para_1_T * conj(a.perp_0_L) - a.para_1_T * conj(a.para_t_L)
                            )
                        ) / sqrt(2.0);

                // K_{4sc}
                _k[8] = a.alpha * beta * real(
                              a.perp_1_L * conj(a.para_0_L) - a.para_1_L * conj(a.perp_0_L)
                            + a.perp_0_T * conj(a.para_1_T) - a.perp_1_T * conj(a.para_0_T)
                        ) / sqrt(2.0);

                // K_{4s}
                _k[9] = a.alpha * real(
                            +                (a.para_1_L * conj(a.para_0_L) - a.perp_1_L * conj(a.perp_0_L))
                            + (1.0 - beta) * (
                                + a.para_1_L * conj(a.perp_t_L) - a.perp_1_L * conj(a.para_t_L)
                                + a.para_1_T * conj(a.para_0_T) - a.perp_1_T * conj(a.perp_0_T)
                            )
                            + sqrt1mbeta   * (
                                + a.perp_0_T * conj(a.perp_1_L) + a.perp_1_T * conj(a.perp_0_L) + a.perp_1_T * conj(a.para_t_L)
                                - a.para_0_T * conj(a.para_1_L) - a.para_1_T * conj(a.para_0_L) - a.para_1_T * conj(a.perp_t_L)
                            )
                        ) / sqrt(2.0);
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

            inline double k1()  const { return _k[0]; }
            inline double k2()  const { return _k[1]; }
            inline double k3()  const { return _k[2]; }
            inline double k4()  const { return _k[3]; }
            inline double k5()  const { return _k[4]; }
            inline double k6()  const { return _k[5]; }
            inline double k7()  const { return _k[6]; }
            inline double k8()  const { return _k[7]; }
            inline double k9()  const { return _k[8]; }
            inline double k10() const { return _k[9]; }

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

            inline double d4gamma(const double & c_lep, const double & c_lam, const double & phi) const
            {
                const double c2_lep = c_lep * c_lep;
                const double s2_lep = 1.0 - c2_lep;
                const double s_lep  = std::sqrt(s2_lep);
                const double s_lam  = std::sqrt(1.0 - c_lam * c_lam);
                const double c_phi  = std::cos(phi), s_phi = std::sin(phi);

                // cf. [BKTvD2019], p. 2, eqs. (2.3) and (2.4)
                return 3.0 / (8.0 * M_PI) * (
                       k1ss() * s2_lep        + k1cc() * c2_lep + k1c() * c_lep
                    + (k2ss() * s2_lep        + k2cc() * c2_lep + k2c() * c_lep) * c_lam
                    + (k3sc() * s_lep * c_lep + k3s()  * s_lep) * s_lam * s_phi
                    + (k4sc() * s_lep * c_lep + k4s()  * s_lep) * s_lam * c_phi
                );
            }
        };
    }

    template <> struct Implementation<LambdaCToLambdaLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        Parameters parameters;

        UsedParameter hbar;
        UsedParameter tau_Lambda_c;

        UsedParameter g_fermi;

        LeptonFlavorOption opt_l;
        UsedParameter m_l;

        UsedParameter m_Lambda_c;
        UsedParameter m_Lambda;
        UsedParameter alpha;

        UsedParameter mu;

        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            parameters(p),
            hbar(p["QM::hbar"], u),
            tau_Lambda_c(p["life_time::Lambda_c"], u),
            g_fermi(p["WET::G_Fermi"], u),
            opt_l(o, options, "l"_ok),
            m_l(p["mass::" + opt_l.str()], u),
            m_Lambda_c(p["mass::Lambda_c"], u),
            m_Lambda(p["mass::Lambda"], u),
            alpha(p["Lambda::alpha"], u),
            mu(p["scnu" + opt_l.str() + opt_l.str() + "::mu"], u),
            form_factors(FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_c->Lambda::" + o.get("form-factors"_ok, "BMRvD2022"), p, o))
        {
            u.uses(*form_factors);
            u.uses(*model);
        }

        const complex<double> norm(const double & q2) const
        {
            double lam = lambda(m_Lambda_c * m_Lambda_c, m_Lambda * m_Lambda, q2);

            if ((lam <= 0) || (q2 <= m_l * m_l))
                return 0.0;

            return g_fermi() * model->ckm_cs() * (1.0 - m_l * m_l / q2) * std::sqrt(q2 / 3.0 / 128 / power_of<3>(M_PI * m_Lambda_c()) * std::sqrt(lam));
        }

        lambdac_to_lambda_l_nu::Amplitudes amplitudes(const double & s)
        {
            using std::sqrt;

            lambdac_to_lambda_l_nu::Amplitudes result;

            const auto wc = model->wet_scnul(opt_l.value(), false);
            const complex<double> cvl = wc.cvl();
            const complex<double> cvr = wc.cvr();
            const complex<double> csl = wc.csl();
            const complex<double> csr = wc.csr();
            const complex<double> ct  = wc.ct();

            // baryonic form factors (10)
            const double fftV  = form_factors->f_time_v(s);
            const double ff0V  = form_factors->f_long_v(s);
            const double ffpV  = form_factors->f_perp_v(s);
            const double fftA  = form_factors->f_time_a(s);
            const double ff0A  = form_factors->f_long_a(s);
            const double ffpA  = form_factors->f_perp_a(s);
            const double ff0T  = form_factors->f_long_t(s);
            const double ff0T5 = form_factors->f_long_t5(s);
            const double ffpT  = form_factors->f_perp_t(s);
            const double ffpT5 = form_factors->f_perp_t5(s);
            // running quark masses
            const double mcatmu = model->m_c_msbar(mu);
            const double msatmu = model->m_s_msbar(mu);

            // kinematics
            const double beta = (1.0 - m_l * m_l / s);
            const double m_l_hat = std::sqrt(1.0 - beta);
            const double sqrtsminus = std::sqrt(power_of<2>(m_Lambda_c - m_Lambda) - s);
            const double sqrtsplus  = std::sqrt(power_of<2>(m_Lambda_c + m_Lambda) - s);
            const double sqrts      = std::sqrt(s);

            // normalization
            const complex<double> N = norm(s);

            // b->c case transversity amplitudes
            // cf. [BKTvD2019], eqs. (2.18)-(2.23), p. 6, including contributions from the vector and scalar operators.
            result.perp_1_L = -2.0 * N * ffpV * (cvl + cvr) * sqrtsminus;
            result.para_1_L = +2.0 * N * ffpA * (cvl - cvr) * sqrtsplus;
            result.perp_0_L = +std::sqrt(2.0) * N * ff0V * ((m_Lambda_c + m_Lambda) / sqrts) * (cvl + cvr) * sqrtsminus;
            result.para_0_L = -std::sqrt(2.0) * N * ff0A * ((m_Lambda_c - m_Lambda) / sqrts) * (cvl - cvr) * sqrtsplus;
            result.perp_t_L = +std::sqrt(2.0) * N * sqrtsplus  * fftV * ( ((m_Lambda_c - m_Lambda) / sqrts) * (cvl + cvr) + ((m_Lambda_c - m_Lambda) / (mcatmu - msatmu)) * (csl + csr) / m_l_hat );
            result.para_t_L = -std::sqrt(2.0) * N * sqrtsminus * fftA * ( ((m_Lambda_c + m_Lambda) / sqrts) * (cvl - cvr) - ((m_Lambda_c + m_Lambda) / (mcatmu + msatmu)) * (csl - csr) / m_l_hat );

            // cf. [BKTvD2019], eqs. (2.26)-(2.29), p. 6, including contributions from the tensor operator.
            result.para_0_T = -sqrt(8.0) * N * ff0T5 * sqrtsplus  * ct;
            result.perp_0_T = -sqrt(8.0) * N * ff0T  * sqrtsminus * ct;
            result.para_1_T = +sqrt(4.0) * N * ffpT5 * sqrtsplus  * ct * (m_Lambda_c - m_Lambda) / sqrts;
            result.perp_1_T = +sqrt(4.0) * N * ffpT  * sqrtsminus * ct * (m_Lambda_c + m_Lambda) / sqrts;

            result.alpha = this->alpha();
            result.beta  = beta;

            return result;
        }

        std::array<double, 10> _differential_angular_observables(const double & q2)
        {
            return lambdac_to_lambda_l_nu::AngularObservables(this->amplitudes(q2))._k;
        }

        std::array<double, 10> _integrated_angular_observables(const double & q2_min, const double & q2_max)
        {
            std::function<std::array<double, 10> (const double &)> integrand(std::bind(&Implementation::_differential_angular_observables, this, std::placeholders::_1));

            return integrate1D(integrand, 64, q2_min, q2_max);
        }

        inline lambdac_to_lambda_l_nu::AngularObservables differential_angular_observables(const double & q2)
        {
            return lambdac_to_lambda_l_nu::AngularObservables{ _differential_angular_observables(q2) };
        }

        inline lambdac_to_lambda_l_nu::AngularObservables integrated_angular_observables(const double & q2_min, const double & q2_max)
        {
            return lambdac_to_lambda_l_nu::AngularObservables{ _integrated_angular_observables(q2_min, q2_max) };
        }
    };

    const std::vector<OptionSpecification>
    Implementation<LambdaCToLambdaLeptonNeutrino>::options
    {
        { "l"_ok, { "e", "mu", "tau" }, "mu" }
    };

    LambdaCToLambdaLeptonNeutrino::LambdaCToLambdaLeptonNeutrino(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<LambdaCToLambdaLeptonNeutrino>(new Implementation<LambdaCToLambdaLeptonNeutrino>(p, o, *this))
    {
    }

    LambdaCToLambdaLeptonNeutrino::~LambdaCToLambdaLeptonNeutrino()
    {
    }

    /* for four-differential signal PDF */
    double
    LambdaCToLambdaLeptonNeutrino::four_differential_decay_width(const double & q2, const double & c_lep,
            const double & c_lam, const double & phi) const
    {
        return _imp->differential_angular_observables(q2).d4gamma(c_lep, c_lam, phi);
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_decay_width(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).decay_width();
    }

    /* q^2-differential observables */

    double
    LambdaCToLambdaLeptonNeutrino::differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).decay_width() * _imp->tau_Lambda_c / _imp->hbar;
    }

    double
    LambdaCToLambdaLeptonNeutrino::differential_a_fb_leptonic(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).a_fb_leptonic();
    }

    double
    LambdaCToLambdaLeptonNeutrino::differential_a_fb_hadronic(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).a_fb_hadronic();
    }

    double
    LambdaCToLambdaLeptonNeutrino::differential_a_fb_combined(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).a_fb_combined();
    }

    double
    LambdaCToLambdaLeptonNeutrino::differential_fzero(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).f_zero();
    }

    /* q^2-integrated observables */

    double
    LambdaCToLambdaLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).decay_width() * _imp->tau_Lambda_c / _imp->hbar;
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_a_fb_hadronic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_hadronic();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_a_fb_combined(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_combined();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_fzero(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_zero();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_k1ss(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1ss() / o.decay_width();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_k1cc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1cc() / o.decay_width();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_k1c(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1c() / o.decay_width();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_k2ss(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2ss() / o.decay_width();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_k2cc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2cc() / o.decay_width();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_k2c(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2c() / o.decay_width();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_k3sc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k3sc() / o.decay_width();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_k3s(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k3s() / o.decay_width();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_k4sc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k4sc() / o.decay_width();
    }

    double
    LambdaCToLambdaLeptonNeutrino::integrated_k4s(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k4s() / o.decay_width();
    }

    const std::string
    LambdaCToLambdaLeptonNeutrino::description = "\
    The decay Lambda_c -> Lambda lbar nu, where lbar=e^+,mu^+,tau^+ is a charged antilepton.";

    const std::string
    LambdaCToLambdaLeptonNeutrino::kinematics_description_q2 = "\
    The invariant mass of the lbar-nu pair in GeV^2.";

    const std::string
    LambdaCToLambdaLeptonNeutrino::kinematics_description_c_theta_l = "\
    The cosine of the helicity angle between the direction of flight of the charged antilepton and of the Lambda_c in the lbar-nu rest frame.";

    const std::string
    LambdaCToLambdaLeptonNeutrino::kinematics_description_c_theta_L = "\
    The cosine of the helicity angle between the direction of flight of the Lambda and of the pion in the Lambda_c rest frame.";

    const std::string
    LambdaCToLambdaLeptonNeutrino::kinematics_description_phi = "\
    The azimuthal angle between the two decay planes.";

    const std::set<ReferenceName>
    LambdaCToLambdaLeptonNeutrino::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    LambdaCToLambdaLeptonNeutrino::begin_options()
    {
        return Implementation<LambdaCToLambdaLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    LambdaCToLambdaLeptonNeutrino::end_options()
    {
        return Implementation<LambdaCToLambdaLeptonNeutrino>::options.cend();
    }
}
