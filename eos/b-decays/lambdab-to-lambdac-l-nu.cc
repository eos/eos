/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Ahmet Kokulu
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
#include <eos/b-decays/lambdab-to-lambdac-l-nu.hh>
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
#include <eos/utils/integrate.hh>
#include <eos/utils/save.hh>
#include <eos/utils/options-impl.hh>

#include <cmath>
#include <functional>
#include <map>
#include <string>

namespace eos
{
    using std::norm;

    namespace lambdab_to_lambdac_l_nu
    {
        struct Amplitudes
        {
            complex<double> a_perp_0_L;
            complex<double> a_para_0_L;
            complex<double> a_perp_1_L;
            complex<double> a_para_1_L;
            complex<double> a_perp_t_L;
            complex<double> a_para_t_L;
            double alpha;
            double v;
        };

        struct AngularObservables
        {
            std::array<double, 10> _k;

            AngularObservables(const Amplitudes & a)
            {
                // charged lepton velocity in the dilepton rest frame
                double v = a.v;

                _k[0] = ( (2.0 - v) * (std::norm(a.a_perp_1_L) + std::norm(a.a_para_1_L)) +
                          2.0 * std::norm(a.a_perp_0_L) +
                         2.0 * std::norm(a.a_para_0_L) + 2.0 * (1.0 - v) * (std::norm(a.a_perp_t_L) + std::norm(a.a_para_t_L)) ) / 4.0;

                _k[1] = ( std::norm(a.a_perp_1_L) +
                         std::norm(a.a_para_1_L) + (1.0 - v) * ( std::norm(a.a_perp_0_L) + std::norm(a.a_para_0_L) + std::norm(a.a_perp_t_L) + std::norm(a.a_para_t_L) ) ) / 2.0;

                _k[2] = std::real( a.a_perp_1_L*std::conj(a.a_para_1_L) + (1.0 - v) * ( a.a_perp_0_L*std::conj(a.a_perp_t_L) + a.a_para_0_L*std::conj(a.a_para_t_L) ) );

                _k[3] = a.alpha * std::real( 2.0 * (1.0 - v) * a.a_perp_t_L*std::conj(a.a_para_t_L) + (2.0 - v) * a.a_perp_1_L*std::conj(a.a_para_1_L) + 2.0 * a.a_perp_0_L*std::conj(a.a_para_0_L) ) / 2.0;

                _k[4] = a.alpha * std::real( a.a_perp_1_L*std::conj(a.a_para_1_L) + (1.0 - v) * (a.a_perp_0_L*std::conj(a.a_para_0_L) + a.a_perp_t_L*std::conj(a.a_para_t_L)) );

                _k[5] = ( std::norm(a.a_perp_1_L) +
                         std::norm(a.a_para_1_L) + 2.0 * (1.0 - v) * std::real( a.a_perp_0_L*std::conj(a.a_para_t_L) + a.a_para_0_L*std::conj(a.a_perp_t_L) ) ) * a.alpha / 2.0;

                _k[6] = v * std::imag( a.a_perp_1_L*std::conj(a.a_perp_0_L) -
                                      a.a_para_1_L*std::conj(a.a_para_0_L) ) * a.alpha / std::sqrt(2.0);

                _k[7] = std::imag( - a.a_perp_1_L*std::conj(a.a_para_0_L) +
                                  a.a_para_1_L*std::conj(a.a_perp_0_L) + (1.0 - v) * ( a.a_para_1_L*std::conj(a.a_para_t_L) - a.a_perp_1_L*std::conj(a.a_perp_t_L) ) ) * a.alpha / std::sqrt(2.0);

                _k[8] = v * std::real( a.a_perp_1_L*std::conj(a.a_para_0_L) -
                                      a.a_para_1_L*std::conj(a.a_perp_0_L) ) * a.alpha / std::sqrt(2.0);

                _k[9] = std::real( -a.a_perp_1_L*std::conj(a.a_perp_0_L) +
                                  a.a_para_1_L*std::conj(a.a_para_0_L) + (1.0 - v) * ( a.a_para_1_L*std::conj(a.a_perp_t_L) - a.a_perp_1_L*std::conj(a.a_para_t_L) ) ) * a.alpha / std::sqrt(2.0);
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
        };
    }

    /**/

    template <> struct Implementation<LambdaBToLambdaCLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        Parameters parameters;

        UsedParameter hbar;
        UsedParameter tau_Lambda_b;

        UsedParameter g_fermi;

        SwitchOption opt_l;
        UsedParameter m_l;

        UsedParameter m_Lambda_b;
        UsedParameter m_Lambda_c;
        UsedParameter alpha;

        UsedParameter mu;

        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            hbar(p["hbar"], u),
            tau_Lambda_b(p["life_time::Lambda_b"], u),
            g_fermi(p["G_Fermi"], u),
            opt_l(o, "l", {"e", "mu", "tau"}, "mu"),
            m_l(p["mass::" + opt_l.value()], u),
            m_Lambda_b(p["mass::Lambda_b"], u),
            m_Lambda_c(p["mass::Lambda_c"], u),
            alpha(p["Lambda_c::alpha"], u),
            mu(p["mu"], u)
        {
            form_factors = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda_c::" + o.get("form-factors", "DKMR2017"), p, o);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);

        }

        const complex<double> norm(const double & s) const
        {
            // charged lepton velocity in the dilepton rest frame
            double v = (1.0 - m_l * m_l / s);
            double lam = lambda(m_Lambda_b * m_Lambda_b, m_Lambda_c * m_Lambda_c, s);

            return g_fermi() * 4.0 * M_PI * model->ckm_cb() * v * std::sqrt(s / 3.0 / 2048 / std::pow(M_PI, 5.0) / power_of<3>(m_Lambda_b()) * std::sqrt(lam));
        }

        lambdab_to_lambdac_l_nu::Amplitudes amplitudes(const double & s)
        {
            lambdab_to_lambdac_l_nu::Amplitudes result;

            // define below the b->c WCs in EOS basis
            const WilsonCoefficients<BToC> wc = model->wilson_coefficients_b_to_c(opt_l.value(), false);

            const complex<double> cvl        = wc.cvl();
            const complex<double> cvr        = wc.cvr();
            const complex<double> csl        = wc.csl();
            const complex<double> csr        = wc.csr();
            const complex<double> ct         = wc.ct();

            // baryonic form factors (10)
            double fftV  = form_factors->f_time_v(s);
            double ff0V  = form_factors->f_long_v(s);
            double ffpV  = form_factors->f_perp_v(s);
            double fftA  = form_factors->f_time_a(s);
            double ff0A  = form_factors->f_long_a(s);
            double ffpA  = form_factors->f_perp_a(s);
            double ff0T  = form_factors->f_long_t(s);
            double ff0T5 = form_factors->f_long_t5(s);
            double ffpT  = form_factors->f_perp_t(s);
            double ffpT5 = form_factors->f_perp_t5(s);
            // running quark masses
            double mbatmu = model->m_b_msbar(mu);
            double mcatmu = model->m_c_msbar(mu);
            // charged lepton velocity in the dilepton rest frame
            double v = (1.0 - m_l * m_l / s);
            double m_l_hat = std::sqrt(1.0 - v);

            double sqrtsminus = std::sqrt(power_of<2>(m_Lambda_b - m_Lambda_c) - s), sqrtsplus = std::sqrt(power_of<2>(m_Lambda_b + m_Lambda_c) - s), sqrts = std::sqrt(s);
            const complex<double> N = norm(s);

            // b->c case transversity amplitudes A's. cf. from [BKvD2019]
            // VA & SP operators' contributions
            result.a_perp_1_L = -2.0 * N * ffpV * (cvl + cvr) * sqrtsminus;
            result.a_para_1_L = +2.0 * N * ffpA * (cvl - cvr) * sqrtsplus;
            result.a_perp_0_L = +std::sqrt(2.0) * N * ff0V * ((m_Lambda_b + m_Lambda_c) / sqrts) * (cvl + cvr) * sqrtsminus;
            result.a_para_0_L = -std::sqrt(2.0) * N * ff0A * ((m_Lambda_b - m_Lambda_c) / sqrts) * (cvl - cvr) * sqrtsplus;
            result.a_perp_t_L = +std::sqrt(2.0) * N * sqrtsplus * fftV * ( ((m_Lambda_b - m_Lambda_c) / sqrts) * (cvl + cvr) + ((m_Lambda_b - m_Lambda_c) / (mbatmu - mcatmu)) * (csl + csr) / m_l_hat );
            result.a_para_t_L = -std::sqrt(2.0) * N * sqrtsminus * fftA * ( ((m_Lambda_b + m_Lambda_c) / sqrts) * (cvl - cvr) - ((m_Lambda_b + m_Lambda_c) / (mbatmu + mcatmu)) * (csl - csr) / m_l_hat );

            result.alpha = this->alpha();
            result.v     = v;

            return result;
        }

        std::array<double, 10> _differential_angular_observables(const double & s)
        {
            return lambdab_to_lambdac_l_nu::AngularObservables(this->amplitudes(s))._k;
        }

        // define below integrated observables in generic form
        std::array<double, 10> _integrated_angular_observables(const double & s_min, const double & s_max)
        {
            std::function<std::array<double, 10> (const double &)> integrand(std::bind(&Implementation::_differential_angular_observables, this, std::placeholders::_1));
            // second argument of integrate1D is some power of 2
            return integrate1D(integrand, 32, s_min, s_max);
        }

        inline lambdab_to_lambdac_l_nu::AngularObservables differential_angular_observables(const double & s)
        {
            return lambdab_to_lambdac_l_nu::AngularObservables{ _differential_angular_observables(s) };
        }

        inline lambdab_to_lambdac_l_nu::AngularObservables integrated_angular_observables(const double & s_min, const double & s_max)
        {
            return lambdab_to_lambdac_l_nu::AngularObservables{ _integrated_angular_observables(s_min, s_max) };
        }
    };

    LambdaBToLambdaCLeptonNeutrino::LambdaBToLambdaCLeptonNeutrino(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<LambdaBToLambdaCLeptonNeutrino>(new Implementation<LambdaBToLambdaCLeptonNeutrino>(p, o, *this))
    {
    }

    LambdaBToLambdaCLeptonNeutrino::~LambdaBToLambdaCLeptonNeutrino()
    {
    }

    /* q^2-differential observables */
    double
    LambdaBToLambdaCLeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_angular_observables(s).decay_width() * _imp->tau_Lambda_b / _imp->hbar;
    }

    double
    LambdaBToLambdaCLeptonNeutrino::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_leptonic();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::differential_a_fb_hadronic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_hadronic();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::differential_a_fb_combined(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_combined();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::differential_fzero(const double & s) const
    {
        return _imp->differential_angular_observables(s).f_zero();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::differential_ratio_tau_mu(const double &s) const
    {
        double br_tau;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            br_tau = this->differential_branching_ratio(s);
        }

        double br_mu;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            br_mu = this->differential_branching_ratio(s);
        }

        return br_tau / br_mu;
    }

    double
    LambdaBToLambdaCLeptonNeutrino::differential_ratio_a_fb_hadronic_tau_mu(const double &s) const
    {
        double a_fb_hadronic_tau;
        {
            Save<Parameter,     double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            a_fb_hadronic_tau = this->differential_a_fb_hadronic(s);
        }

        double a_fb_hadronic_mu;
        {
            Save<Parameter,    double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            a_fb_hadronic_mu = this->differential_a_fb_hadronic(s);
        }

        return a_fb_hadronic_tau / a_fb_hadronic_mu;
    }

    /* q^2-integrated observables */
    double
    LambdaBToLambdaCLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).decay_width() * _imp->tau_Lambda_b / _imp->hbar;
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_hadronic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_hadronic();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_a_fb_combined(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_combined();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_fzero(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_zero();
    }

    //*
    double
    LambdaBToLambdaCLeptonNeutrino::integrated_ratio_tau_mu(const double & s_min_mu, const double & s_min_tau, const double & s_max_mu, const double & s_max_tau) const
    {

        double br_mu;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "mu");
            
            br_mu = this->integrated_branching_ratio(s_min_mu, s_max_mu);
        }

        double br_tau;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "tau");
            
            br_tau = this->integrated_branching_ratio(s_min_tau, s_max_tau);
        }
        
        return br_tau / br_mu;
    }

    //*
    double
    LambdaBToLambdaCLeptonNeutrino::integrated_ratio_a_fb_hadronic_tau_mu(const double & s_min_mu, const double & s_min_tau, const double & s_max_mu, const double & s_max_tau) const
    {

        double integrated_a_fb_hadronic_mu;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "mu");
            
            integrated_a_fb_hadronic_mu = this->integrated_a_fb_hadronic(s_min_mu, s_max_mu);
        }

        double integrated_a_fb_hadronic_tau;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "tau");
            
            integrated_a_fb_hadronic_tau = this->integrated_a_fb_hadronic(s_min_tau, s_max_tau);
        }

        return integrated_a_fb_hadronic_tau / integrated_a_fb_hadronic_mu;
    }

    //*
    double
    LambdaBToLambdaCLeptonNeutrino::integrated_k1ss(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1ss() / o.decay_width();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_k1cc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1cc() / o.decay_width();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_k1c(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1c() / o.decay_width();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_k2ss(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2ss() / o.decay_width();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_k2cc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2cc() / o.decay_width();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_k2c(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2c() / o.decay_width();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_k3sc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k3sc() / o.decay_width();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_k3s(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k3s() / o.decay_width();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_k4sc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k4sc() / o.decay_width();
    }

    double
    LambdaBToLambdaCLeptonNeutrino::integrated_k4s(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k4s() / o.decay_width();
    }
    //*


    const std::string
    LambdaBToLambdaCLeptonNeutrino::description = "\
    The decay Lambda_b -> Lambda_c l nu, where l=e,mu,tau is a lepton.";

    const std::string
    LambdaBToLambdaCLeptonNeutrino::kinematics_description_s = "\
    The invariant mass of the l-nubar pair in GeV^2.";

}
