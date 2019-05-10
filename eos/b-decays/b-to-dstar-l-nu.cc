/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018, 2019 Ahmet Kokulu
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

#include <eos/form-factors/form-factors.hh>
#include <eos/b-decays/b-to-dstar-l-nu.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/integrate-impl.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/save.hh>

#include <cmath>
#include <functional>
#include <map>
#include <string>

namespace eos
{
    using std::norm;

    namespace b_to_dstar_l_nu
    {
        struct Amplitudes
        {
            complex<double> a_0;
            complex<double> a_0_T;
            complex<double> a_plus;
            complex<double> a_plus_T;
            complex<double> a_minus;
            complex<double> a_minus_T;
            complex<double> a_P;
            complex<double> a_t;
            complex<double> a_t_P;
            complex<double> a_para;
            complex<double> a_para_T;
            complex<double> a_perp;
            complex<double> a_perp_T;
            double v;
            double NF;
        };

        // angular observables V's. cf. from [DSD2014], p. 16, redifined V's in order to include NF
        struct AngularObservables
        {
            std::array<double, 12> _vv;

            AngularObservables(const Amplitudes & a)
            {
                // charged lepton velocity in the dilepton rest frame
                double v  = a.v;
                double NF = a.NF;

                _vv[0] = NF * 2.0 * ( (2.0 - v) * (std::norm(a.a_0) + 16.0 * std::norm(a.a_0_T)) + 2.0 * (1.0 - v) * std::norm(a.a_t_P) - 16.0 * std::sqrt(1.0 - v) * std::real(a.a_0_T * std::conj(a.a_0)) );

                _vv[1] = NF * 2.0 * v * ( - std::norm(a.a_0) + 16.0 * std::norm(a.a_0_T) );

                _vv[2] = - NF * 8.0 * std::real( (1.0 - v) * a.a_t_P * std::conj(a.a_0) - 4.0 * std::sqrt(1.0 - v) * a.a_t_P * std::conj(a.a_0_T) );

                _vv[3] = NF * ( (4.0 - v) * (std::norm(a.a_para) + std::norm(a.a_perp)) / 2.0 + 8.0 * (4.0 - 3.0 * v) * (std::norm(a.a_para_T) + std::norm(a.a_perp_T)) - 16.0 * std::sqrt(1.0 - v) * std::real( a.a_para_T * std::conj(a.a_para) + a.a_perp_T * std::conj(a.a_perp) ) );

                _vv[4] = NF * v * ( (std::norm(a.a_para) + std::norm(a.a_perp)) / 2.0 - 8.0 * (std::norm(a.a_para_T) + std::norm(a.a_perp_T)) );

                _vv[5] = NF * 4.0 * std::real( - a.a_para * std::conj(a.a_perp) - 16.0 * (1.0 - v) * a.a_para_T * std::conj(a.a_perp_T) + 4.0 * std::sqrt(1.0 - v) * ( a.a_perp_T * std::conj(a.a_para) + a.a_para_T * std::conj(a.a_perp) ) );

                _vv[6] = NF * v * ( - (std::norm(a.a_para) - std::norm(a.a_perp)) + 16.0 * (std::norm(a.a_para_T) - std::norm(a.a_perp_T)) );

                _vv[7] = NF * 2.0 * v * std::imag( a.a_para * std::conj(a.a_perp) );

                _vv[8] = NF * std::sqrt(2.0) * v * std::real( a.a_para * std::conj(a.a_0) - 16.0 * a.a_para_T * std::conj(a.a_0_T) );

                _vv[9] = NF * 2.0 * std::sqrt(2.0) * std::real( - a.a_perp * std::conj(a.a_0) + (1.0 - v) * ( a.a_para * std::conj(a.a_t_P) - 16.0 * a.a_perp_T * std::conj(a.a_0_T) ) + 4.0 * std::sqrt(1.0 - v) * ( a.a_0_T * std::conj(a.a_perp) + a.a_perp_T * std::conj(a.a_0) - a.a_para_T * std::conj(a.a_t_P) ) );

                _vv[10] = NF * 2.0 * std::sqrt(2.0) * std::imag( - a.a_para * std::conj(a.a_0) + (1.0 - v) * a.a_perp * std::conj(a.a_t_P) + 4.0 * std::sqrt(1.0 - v) * ( a.a_0_T * std::conj(a.a_para) - a.a_para_T * std::conj(a.a_0) + a.a_perp_T * std::conj(a.a_t_P) ) );

                _vv[11] = NF * std::sqrt(2.0) * v * std::imag( a.a_perp * std::conj(a.a_0) );
            }

            AngularObservables(const std::array<double, 12> & vv) :
                _vv(vv)
            {
            }

            inline double vv10()    const  { return _vv[0]; }
            inline double vv20()    const  { return _vv[1]; }
            inline double vv30()    const  { return _vv[2]; }
            inline double vv1T()    const  { return _vv[3]; }
            inline double vv2T()    const  { return _vv[4]; }
            inline double vv3T()    const  { return _vv[5]; }
            inline double vv4T()    const  { return _vv[6]; }
            inline double vv5T()    const  { return _vv[7]; }
            inline double vv10T()   const  { return _vv[8]; }
            inline double vv20T()   const  { return _vv[9]; }
            inline double vv30T()   const  { return _vv[10];}
            inline double vv40T()   const  { return _vv[11];}

            // longitudinal polarization amplitude
            inline double normalized_amplitude_polarization_L() const
            {
                return  vv10() - vv20() / 3.0;
            }

            // transverse polarization amplitude
            inline double normalized_amplitude_polarization_T() const
            {
                return  2.0 * (vv1T() - vv2T() / 3.0);
            }

            // redefined decay width
            inline double normalized_decay_width() const
            {
                return 3.0 / 4.0 * ( normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T() );
            }

            // polarization fraction
            inline double f_L() const
            {
                return  normalized_amplitude_polarization_L() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            // a_fb leptonic
            inline double a_fb_leptonic() const
            {
                return  (vv3T() + vv30() / 2.0) / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            // transverse azimuthal asymmetries
            inline double a_c_1() const
            {
                return  4.0 * vv4T() / (3.0 * (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T()));
            }

            inline double a_c_2() const
            {
                return  vv20T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            inline double a_c_3() const
            {
                return  vv10T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            // T-odd CP asymmetries
            inline double a_t_1() const
            {
                return  4.0 * vv5T() / (3.0 * (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T()));
            }

            inline double a_t_2() const
            {
                return  vv30T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            inline double a_t_3() const
            {
                return  vv40T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }
        };
    }

    /**/
    template <> struct Implementation<BToDstarLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        Parameters parameters;

        SwitchOption opt_q;

        UsedParameter hbar;

        UsedParameter tau_B;

        UsedParameter g_fermi;

        SwitchOption opt_l;

        UsedParameter m_l;

        UsedParameter m_B;

        UsedParameter m_Dstar;

        UsedParameter mu;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            opt_q(o, "q", { "u", "d" }, "d"),
            hbar(p["hbar"], u),
            tau_B(p["life_time::B_" + opt_q.value()], u),
            g_fermi(p["G_Fermi"], u),
            opt_l(o, "l", {"e", "mu", "tau"}, "mu"),
            m_l(p["mass::" + opt_l.value()], u),
            m_B(p["mass::B_" + opt_q.value()], u),
            m_Dstar(p["mass::D^*_" + opt_q.value()], u),
            mu(p["mu"], u)
        {
            form_factors = FormFactorFactory<PToV>::create("B->D^*::" + o.get("form-factors", "BSZ2015"), p, o);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);

        }

        // normalization cf. [DSD2014] eq. (7), p. 5
        double norm(const double & s) const
        {
            // charged lepton velocity in the dilepton rest frame
            double v    = (1.0 - m_l * m_l / s);
            double lam  = lambda(m_B * m_B, m_Dstar * m_Dstar, s);
            double p    = std::sqrt(lam) / (2.0 * m_B);
            // normalized prefactor (|Vcb|^2=1)
            return power_of<2>(g_fermi()) * p * s * power_of<2>(v) / (3.0 * 64.0 * power_of<3>(M_PI) * m_B * m_B);
        }

        b_to_dstar_l_nu::Amplitudes amplitudes(const double & s)
        {
            b_to_dstar_l_nu::Amplitudes result;

            // NP contributions in EFT including tensor operator cf. [DSD2014], p. 3
            const WilsonCoefficients<BToC> wc = model->wilson_coefficients_b_to_c(opt_l.value(), false);
            const complex<double> VL = wc.cvl() - 1.0;
            const complex<double> VR = wc.cvr();
            const complex<double> SL = wc.csl();
            const complex<double> SR = wc.csr();
            const complex<double> gV = VR + VL;
            const complex<double> gA = VR - VL;
            const complex<double> gS = SR + SL;
            const complex<double> gP = SR - SL;
            const complex<double> TL = wc.ct();

            // form factors
            double aff0  = form_factors->a_0(s);
            double aff1  = form_factors->a_1(s);
            double aff2  = form_factors->a_2(s);
            double vff   = form_factors->v(s);
            double tff1  = form_factors->t_1(s);
            double tff2  = form_factors->t_2(s);
            double tff3  = form_factors->t_3(s);
            // running quark masses
            double mbatmu = model->m_b_msbar(mu);
            double mcatmu = model->m_c_msbar(mu);
            // charged lepton velocity in the dilepton rest frame
            double v        = (1.0 - m_l * m_l / s);
            double m_l_hat  = std::sqrt(1.0 - v);
            double lam      = lambda(m_B * m_B, m_Dstar * m_Dstar, s);

            double sqrts = std::sqrt(s);
            double NF = norm(s);

            // transversity amplitudes A's. cf. [DSD2014], p.17
            result.a_0          = (1.0 - gA) * (m_B + m_Dstar) / (2.0 * m_Dstar * sqrts) * ( (m_B * m_B - m_Dstar * m_Dstar - s) * aff1 - lam * aff2 / power_of<2>(m_B + m_Dstar) );
            result.a_0_T        = TL / (2.0 * m_Dstar) * ( (m_B * m_B + 3.0 * m_Dstar * m_Dstar - s) * tff2 - lam * tff3 / (m_B * m_B - m_Dstar * m_Dstar) );
            result.a_plus       = ( (m_B + m_Dstar) * aff1 * (1.0 - gA) - std::sqrt(lam) * vff * (1.0 + gV) / (m_B + m_Dstar) );
            result.a_minus      = ( (m_B + m_Dstar) * aff1 * (1.0 - gA) + std::sqrt(lam) * vff * (1.0 + gV) / (m_B + m_Dstar) );
            result.a_plus_T     = TL * ( (m_B * m_B - m_Dstar * m_Dstar) * tff2 / sqrts + std::sqrt(lam) * tff1 / sqrts );
            result.a_minus_T    = TL * ( (m_B * m_B - m_Dstar * m_Dstar) * tff2 / sqrts - std::sqrt(lam) * tff1 / sqrts );
            result.a_t          =  std::sqrt(lam) * aff0 * (1.0 - gA) / sqrts;
            result.a_P          =  std::sqrt(lam) * aff0 * gP / (mbatmu + mcatmu);
            result.a_t_P        =  result.a_t + result.a_P / m_l_hat;
            result.a_para       =  (result.a_plus + result.a_minus) / std::sqrt(2.0);
            result.a_para_T     =  (result.a_plus_T + result.a_minus_T) / std::sqrt(2.0);
            result.a_perp       =  (result.a_plus - result.a_minus) / std::sqrt(2.0);
            result.a_perp_T     =  (result.a_plus_T - result.a_minus_T) / std::sqrt(2.0);

            result.v     = v;
            result.NF    = NF;

            return result;
        }

        std::array<double, 12> _differential_angular_observables(const double & s)
        {
            return b_to_dstar_l_nu::AngularObservables(this->amplitudes(s))._vv;
        }

        // define below integrated observables in generic form
        std::array<double, 12> _integrated_angular_observables(const double & s_min, const double & s_max)
        {
            std::function<std::array<double, 12> (const double &)> integrand(std::bind(&Implementation::_differential_angular_observables, this, std::placeholders::_1));
            // second argument of integrate1D is some power of 2
            return integrate1D(integrand, 64, s_min, s_max);
        }

        inline b_to_dstar_l_nu::AngularObservables differential_angular_observables(const double & s)
        {
            return b_to_dstar_l_nu::AngularObservables{ _differential_angular_observables(s) };
        }

        inline b_to_dstar_l_nu::AngularObservables integrated_angular_observables(const double & s_min, const double & s_max)
        {
            return b_to_dstar_l_nu::AngularObservables{ _integrated_angular_observables(s_min, s_max) };
        }
    };

    BToDstarLeptonNeutrino::BToDstarLeptonNeutrino(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BToDstarLeptonNeutrino>(new Implementation<BToDstarLeptonNeutrino>(p, o, *this))
    {
    }

    BToDstarLeptonNeutrino::~BToDstarLeptonNeutrino()
    {
    }

    /* q^2-differential observables */

    // |Vcb|=1
    double
    BToDstarLeptonNeutrino::normalized_differential_branching_ratio(const double & s) const
    {
        return _imp->differential_angular_observables(s).normalized_decay_width() * _imp->tau_B / _imp->hbar;
    }

    double
    BToDstarLeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_angular_observables(s).normalized_decay_width() * std::norm(_imp->model->ckm_cb()) * _imp->tau_B / _imp->hbar;
    }

    double
    BToDstarLeptonNeutrino::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_leptonic();
    }

    double
    BToDstarLeptonNeutrino::differential_ratio_tau_mu(const double &s) const
    {
        double br_tau;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "tau");
            br_tau = this->normalized_differential_branching_ratio(s);
        }

        double br_mu;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "mu");
            br_mu = this->normalized_differential_branching_ratio(s);
        }

        return br_tau / br_mu;
    }

    /* q^2-integrated observables */

    // |Vcb|=1
        double
    BToDstarLeptonNeutrino::normalized_integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).normalized_decay_width() * _imp->tau_B / _imp->hbar;
    }
    
    double
    BToDstarLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).normalized_decay_width() * std::norm(_imp->model->ckm_cb()) * _imp->tau_B / _imp->hbar;
    }

    double
    BToDstarLeptonNeutrino::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic();
    }

    double
    BToDstarLeptonNeutrino::integrated_ratio_tau_mu(const double & s_min_mu, const double & s_min_tau, const double & s_max_mu, const double & s_max_tau) const
    {

        double br_mu;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "mu");
            
            br_mu = this->normalized_integrated_branching_ratio(s_min_mu, s_max_mu);
        }

        double br_tau;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "tau");
            
            br_tau = this->normalized_integrated_branching_ratio(s_min_tau, s_max_tau);
        }

        return br_tau / br_mu;
    }

    double
    BToDstarLeptonNeutrino::integrated_amplitude_polarization_L(const double & s_min, const double & s_max) const
    {
        auto   o = _imp->integrated_angular_observables(s_min, s_max);
        return o.normalized_amplitude_polarization_L() * std::norm(_imp->model->ckm_cb());
    }

    double
    BToDstarLeptonNeutrino::integrated_amplitude_polarization_T(const double & s_min, const double & s_max) const
    {
        auto   o = _imp->integrated_angular_observables(s_min, s_max);
        return o.normalized_amplitude_polarization_T() * std::norm(_imp->model->ckm_cb());
    }

    double
    BToDstarLeptonNeutrino::integrated_f_L(const double & s_min, const double & s_max) const
    {
        auto   o = _imp->integrated_angular_observables(s_min, s_max);
        return o.f_L();
    }

    double
    BToDstarLeptonNeutrino::integrated_a_c_1(const double & s_min, const double & s_max) const
    {
        auto   o = _imp->integrated_angular_observables(s_min, s_max);
        return o.a_c_1();
    }

    double
    BToDstarLeptonNeutrino::integrated_a_c_2(const double & s_min, const double & s_max) const
    {
        auto   o = _imp->integrated_angular_observables(s_min, s_max);
        return o.a_c_2();
    }

    double
    BToDstarLeptonNeutrino::integrated_a_c_3(const double & s_min, const double & s_max) const
    {
        auto   o = _imp->integrated_angular_observables(s_min, s_max);
        return o.a_c_3();
    }

    double
    BToDstarLeptonNeutrino::integrated_a_t_1(const double & s_min, const double & s_max) const
    {
        auto   o = _imp->integrated_angular_observables(s_min, s_max);
        return o.a_t_1();
    }

    double
    BToDstarLeptonNeutrino::integrated_a_t_2(const double & s_min, const double & s_max) const
    {
        auto   o = _imp->integrated_angular_observables(s_min, s_max);
        return o.a_t_2();
    }

    double
    BToDstarLeptonNeutrino::integrated_a_t_3(const double & s_min, const double & s_max) const
    {
        auto   o = _imp->integrated_angular_observables(s_min, s_max);
        return o.a_t_3();
    }

    //* cf. [DSD2014], eq. (6), p. 5 - normalized(|Vcb|=1)
    double
    BToDstarLeptonNeutrino::normalized_four_differential_decay_width(const double & s, const double & c_theta_l, const double & c_theta_d, const double & phi) const
    {
        // compute d^4 Gamma, cf. [DSD2014], p. 5, eq. (6)
        // Trigonometric identities: Cosine squared of the angles
        double c_theta_d_2 = c_theta_d * c_theta_d;
        double c_theta_l_2 = c_theta_l * c_theta_l;
        double c_phi = std::cos(phi);
        // Sine squared of the angles
        double s_theta_d_2 = 1.0 - c_theta_d_2;
        double s_theta_l_2 = 1.0 - c_theta_l_2;
        // Sine of the angles
        double s_theta_d = std::sqrt(s_theta_d_2);
        double s_theta_l = std::sqrt(s_theta_l_2);
        double s_phi = std::sin(phi);
        // Cosine of twice the angle
        //double c_2_theta_d = 2.0 * c_theta_d_2 - 1.0;
        double c_2_theta_l = 2.0 * c_theta_l_2 - 1.0;
        double c_2_phi = std::cos(2.0 * phi);
        // Sine of twice the angle
        double s_2_theta_d = 2.0 * s_theta_d * c_theta_d;
        double s_2_theta_l = 2.0 * s_theta_l * c_theta_l;
        double s_2_phi = std::sin(2.0 * phi);

        b_to_dstar_l_nu::AngularObservables a_o = _imp->differential_angular_observables(s);

        double result = 9.0 / 32.0 / M_PI * (
                                            (a_o.vv10() + a_o.vv20() * c_2_theta_l + a_o.vv30() * c_theta_l) * c_theta_d_2
                                            +  (a_o.vv1T() + a_o.vv2T() * c_2_theta_l + a_o.vv3T() * c_theta_l) * s_theta_d_2
                                            +  a_o.vv4T() * s_theta_d_2 * s_theta_l_2 * c_2_phi + a_o.vv10T() * s_2_theta_d * s_2_theta_l * c_phi
                                            +  a_o.vv20T() * s_2_theta_d * s_theta_l * c_phi + a_o.vv5T() * s_theta_d_2 * s_theta_l_2 * s_2_phi
                                            +  a_o.vv30T() * s_2_theta_d * s_theta_l * s_phi +  a_o.vv40T() * s_2_theta_d * s_2_theta_l * s_phi
                                            );

        return result;
    }

    const std::string
    BToDstarLeptonNeutrino::description = "\
    The decay B->D^* l nu, where l=e,mu,tau is a lepton.";

    const std::string
    BToDstarLeptonNeutrino::kinematics_description_s = "\
    The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    BToDstarLeptonNeutrino::kinematics_description_c_theta_l = "\
    The cosine of the charged lepton's helicity angle theta_l in the l-nubar rest frame.";

    const std::string
    BToDstarLeptonNeutrino::kinematics_description_c_theta_d = "\
    The cosine of the D's helicity angle theta_d in the D-pi rest frame.";

    const std::string
    BToDstarLeptonNeutrino::kinematics_description_phi = "\
    The azimuthal angle between the D-pi plane and the l-nubar  plane.";

}
