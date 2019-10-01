/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2015 Marzia Bordone
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
#include <eos/b-decays/b-to-psd-l-nu.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/save.hh>

namespace eos
{
    namespace b_to_psd_l_nu
    {
        struct Amplitudes
        {
            // helicity amplitudes, cf. [DSD2014] eqs. 13-14
            complex<double> h_0;
            complex<double> h_t;
            complex<double> h_S;
            complex<double> h_T;
            complex<double> h_tS;
            double v;
            double p;
            double NF;
        };
    }

    /* */
    template <>
    struct Implementation<BToPseudoscalarLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        Parameters parameters;

        SwitchOption opt_U;
        SwitchOption opt_q;

        UsedParameter m_B;

        UsedParameter mu;

        UsedParameter tau_B;

        UsedParameter m_P;

        SwitchOption opt_l;

        UsedParameter m_l;

        UsedParameter g_fermi;

        UsedParameter hbar;

        std::function<double (const double &)> m_U_msbar;
        std::function<complex<double> ()> v_Ub;
        std::function<WilsonCoefficients<ChargedCurrent> (const std::string &, bool)> wc;

        inline std::string _mass_P() const
        {
            switch (opt_U.value()[0])
            {
                case 'c':
                    return std::string("mass::D_") + opt_q.value();
                    break;

                case 'u':
                    switch (opt_q.value()[0])
                    {
                        case 'd':
                            return std::string("mass::pi^+");
                            break;

                        case 's':
                            return std::string("mass::K_u");
                            break;

                        case 'u':
                            return std::string("mass::pi^0");
                            break;

                        default:
                            throw InternalError("Should never reach this part, either!");
                    }
                    break;

                default:
                    throw InternalError("Should never reach this part!");
            };

            return "";
        }

        inline std::string _process() const
        {
            switch (opt_U.value()[0])
            {
                case 'c':
                    switch (opt_q.value()[0])
                    {
                        case 'd':
                        case 'u':
                            return std::string("B->D");
                            break;

                        case 's':
                            return std::string("B_s->D_s");
                            break;

                        default:
                            throw InternalError("Should never reach this part, either!");
                    }
                    break;

                case 'u':
                    switch (opt_q.value()[0])
                    {
                        case 'd':
                        case 'u':
                            return std::string("B->pi");
                            break;

                        case 's':
                            return std::string("B_s->K");
                            break;

                        default:
                            throw InternalError("Should never reach this part, either!");
                    }
                    break;

                default:
                    throw InternalError("Should never reach this part!");
            };

            return "";
        }

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            opt_U(o, "U", { "c", "u" }),
            opt_q(o, "q", { "u", "d", "s" }, "d"),
            m_B(p["mass::B_" + opt_q.value()], u),
            mu(p["mu"], u),
            tau_B(p["life_time::B_" + opt_q.value()], u),
            m_P(p[_mass_P()], u),
            opt_l(o, "l", {"e", "mu", "tau"}, "mu"),
            m_l(p["mass::" + opt_l.value()], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u)
        {
            form_factors = FormFactorFactory<PToP>::create(_process() + "::" + o.get("form-factors", "BSZ2015"), p, o);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            using std::placeholders::_1;
            using std::placeholders::_2;
            if ('u' == opt_U.value()[0])
            {
                m_U_msbar = std::bind(&ModelComponent<components::QCD>::m_u_msbar, model.get(), _1);
                v_Ub      = std::bind(&ModelComponent<components::CKM>::ckm_ub, model.get());
                wc        = std::bind(&ModelComponent<components::DeltaBU1>::wilson_coefficients_b_to_u, model.get(), _1, _2);
            }
            else
            {
                m_U_msbar = std::bind(&ModelComponent<components::QCD>::m_c_msbar, model.get(), _1);
                v_Ub      = std::bind(&ModelComponent<components::CKM>::ckm_cb, model.get());
                wc        = std::bind(&ModelComponent<components::DeltaBC1>::wilson_coefficients_b_to_c, model.get(), _1, _2);
            }

            u.uses(*form_factors);
            u.uses(*model);
        }

        b_to_psd_l_nu::Amplitudes amplitudes(const double & s) const
        {           
            // NP contributions in EFT including tensor operator (cf. [DSD2014]). (in SM cvl=1, and all other couplings are zero)
            auto wc = this->wc(opt_l.value(), false);
            const complex<double> gV = wc.cvr() + (wc.cvl() - 1.0);
            const complex<double> gS = wc.csr() + wc.csl();
            const complex<double> gT = wc.ct();

            // form factors
            const double fp = form_factors->f_p(s);
            const double f0 = form_factors->f_0(s);
            const double fT = form_factors->f_t(s);
            
            // running quark masses
            const double mbatmu = model->m_b_msbar(mu);
            const double mUatmu = m_U_msbar(mu);

            const double lam = lambda(m_B * m_B, m_P * m_P, s);
            const double p = sqrt(lam) / (2.0 * m_B);
            
            // make sure we return NaN if s < m_l^2, v = lepton velocity in the dilepton rest frame
            const double v = (1.0 - m_l * m_l / s);
            const double v2 = v * v;
            const double ml_hat = sqrt(1.0 - v);
            // universal electroweak correction, cf. [S1982]
            const double etaEW = 1.0066;
            const double nP = v2 * m_B * s * power_of<2>(g_fermi() * etaEW) / (256.0 * power_of<3>(M_PI * m_B));
          
            // helicity amplitudes, cf. [DSD2014] eqs. 13-14
            b_to_psd_l_nu::Amplitudes result;

            if (s >= pow(m_l, 2)) {
                result.h_0 = 2.0 * m_B * p * fp * (1.0 + gV) / sqrt(s) ;
                result.h_t = (1.0 + gV) * (m_B * m_B - m_P * m_P) * f0 / sqrt(s);
                result.h_S = - gS * (m_B * m_B - m_P * m_P) * f0 / (mbatmu - mUatmu);
                result.h_T = - 2.0 * m_B * p * fT * gT / (m_B + m_P);
                result.h_tS = result.h_t - result.h_S / ml_hat;
            
                result.v  = v;
                result.p  = p;
                result.NF = nP;
            }
            else { // TODO: return NaN
                                
            }
  
            return result;
        }

        // normalized (|V_Ub| = 1) two-fold-distribution, cf. [DSD2014], eq. (12), p. 6
        double normalized_two_differential_decay_width(const double & s, const double & c_theta_l) const
        {
            
            // TODO: remove this once implemented that amplitudes set to NaN 
            if (s <= pow(m_l, 2))
                return 0.0;

            //  d^2 Gamma, cf. [DSD2014], p. 6, eq. (13)
            // Trigonometric identities
            double c_thl_2 = c_theta_l * c_theta_l;
            double s_thl_2 = 1.0 - c_thl_2;
            double c_2_thl = 2.0 * c_thl_2 - 1.0;

            b_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return 2.0 * amp.NF * amp.p * (
                       std::norm(amp.h_0) * s_thl_2 
                       + (1.0 - amp.v) * power_of<2>(std::abs(amp.h_0) * c_theta_l - std::abs(amp.h_tS))
                       + 8.0 * ( ((2.0 - amp.v) + amp.v * c_2_thl) * std::norm(amp.h_T) 
                           - sqrt(1.0 - amp.v) * std::real(amp.h_T * (std::conj(amp.h_0) - std::conj(amp.h_tS) * c_theta_l)))
                   );
        }

        // normalized to |V_Ub = 1|, obtained using cf. [DSD2014], eq. (12), agrees with Sakaki'13 et al cf. [STTW2013]
        double normalized_differential_decay_width(const double & s) const
        {
            // TODO: remove this once implemented that amplitudes set to NaN 
            if (s <= pow(m_l, 2))
                return 0.0;

           b_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            // normalized(|V_cb|=1) differential decay width
            return 4.0 / 3.0 * amp.NF * amp.p * ( 
                       std::norm(amp.h_0) * (3.0 - amp.v) 
                       + 3.0 * std::norm(amp.h_tS) * (1.0 - amp.v) 
                       + 16.0 * std::norm(amp.h_T) * (3.0 - 2.0 * amp.v)
                       - 24.0 * sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_0))
                   );
        }

        // obtained using cf. [DSD2014], eq. (12), defined as int_1^0 d^2Gamma - int_0^-1 d^2Gamma
        double numerator_differential_a_fb_leptonic(const double & s) const
        {
            // TODO: remove this once implemented that amplitudes set to NaN 
            if (s <= pow(m_l, 2))
                return 0.0;

            b_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return - 4.0 * amp.NF * amp.p * (
                       std::abs(amp.h_0) * std::abs(amp.h_tS) * (1.0 - amp.v) 
                       - 4.0 * sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_tS))
                   );
        }

        // differential decay width
        double differential_decay_width(const double & s) const
        {
            return normalized_differential_decay_width(s) * std::norm(v_Ub());
        }

        // differential branching_ratio
        double differential_branching_ratio(const double & s) const
        {
            return differential_decay_width(s) * tau_B / hbar;
        }

        // "normalized" (|V_Ub|=1) differential branching_ratio
        double normalized_differential_branching_ratio(const double & s) const
        {
            return normalized_differential_decay_width(s) * tau_B / hbar;
        }

        double pdf_q2(const double & q2) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = power_of<2>(m_B() - m_P());

            std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::normalized_differential_branching_ratio, this, std::placeholders::_1);
            const double num   = normalized_differential_branching_ratio(q2);
            const double denom = integrate1D(f, 32, q2_min, q2_max);

            return num / denom;
        }

        double pdf_w(const double & w) const
        {
            const double m_B = this->m_B(), m_B2 = m_B * m_B;
            const double m_P = this->m_P(), m_P2 = m_P * m_P;
            const double q2  = m_B2 + m_P2 - 2.0 * m_B * m_P * w;

            return 2.0 * m_B * m_P * pdf_q2(q2);
        }

        double integrated_pdf_q2(const double & q2_min, const double & q2_max) const
        {
            const double q2_abs_min = power_of<2>(m_l());
            const double q2_abs_max = power_of<2>(m_B() - m_P());

            std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::normalized_differential_branching_ratio, this, std::placeholders::_1);
            const double num   = integrate<GSL::QAGS>(f, q2_min,     q2_max);
            const double denom = integrate<GSL::QAGS>(f, q2_abs_min, q2_abs_max);

            return num / denom;
        }

        double integrated_pdf_w(const double & w_min, const double & w_max) const
        {
            const double m_B    = this->m_B(), m_B2 = m_B * m_B;
            const double m_P    = this->m_P(), m_P2 = m_P * m_P;
            const double q2_max = m_B2 + m_P2 - 2.0 * m_B * m_P * w_min;
            const double q2_min = m_B2 + m_P2 - 2.0 * m_B * m_P * w_max;

            return integrated_pdf_q2(q2_min, q2_max) / (w_max - w_min);
        }

        double lepton_polarization_numerator(const double & q2) const
        {
            const double m_l2 = m_l() * m_l();
            const double m_B  = this->m_B(), m_B2 = m_B * m_B;
            const double m_P  = this->m_P(), m_P2 = m_P * m_P;
            const double p_P  = sqrt(eos::lambda(m_B2, m_P2, q2)) / (2.0 * m_B);
            const double sqrt_q2  = sqrt(q2);

            const double f_p = form_factors->f_p(q2);
            const double f_0 = form_factors->f_0(q2);

            // cf. [CJLP2012]
            const double H_0 = 2.0 * m_B * p_P / sqrt_q2 * f_p;
            const double H_t = (m_B2 - m_P2) / sqrt_q2 * f_0;

            const double H_02 = H_0 * H_0;
            const double H_t2 = H_t * H_t;

            const double nf = p_P * q2 * pow(1.0 - m_l2 / q2, 2);

            // cf. [CJLP2012]], eq. (20), p. 13
            const double num  = H_02 * (1.0 - 0.5 * m_l2 / q2) - 3.0 / 2.0 * m_l2 / q2 * H_t2;

            return nf * num;
        }

        double lepton_polarization_denominator(const double & q2) const
        {
            const double m_l2 = m_l() * m_l();
            const double m_B  = this->m_B(), m_B2 = m_B * m_B;
            const double m_P  = this->m_P(), m_P2 = m_P * m_P;
            const double p_P  = sqrt(eos::lambda(m_B2, m_P2, q2)) / (2.0 * m_B);
            const double sqrt_q2  = sqrt(q2);

            const double f_p = form_factors->f_p(q2);
            const double f_0 = form_factors->f_0(q2);

            // cf. [CJLP2012]
            const double H_0 = 2.0 * m_B * p_P / sqrt_q2 * f_p;
            const double H_t = (m_B2 - m_P2) / sqrt_q2 * f_0;

            const double H_02 = H_0 * H_0;
            const double H_t2 = H_t * H_t;

            const double nf = p_P * q2 * pow(1.0 - m_l2 / q2, 2);

            // cf. [CJLP2012]], eq. (20), p. 13
            const double denom  = H_02 * (1.0 + 0.5 * m_l2 / q2) + 3.0 / 2.0 * m_l2 / q2 * H_t2;

            return nf * denom;
        }

        double lepton_polarization(const double & q2_min, const double & q2_max) const
        {
            std::function<double (const double &)> integrand_num   = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::lepton_polarization_numerator,   this, std::placeholders::_1);
            std::function<double (const double &)> integrand_denom = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::lepton_polarization_denominator, this, std::placeholders::_1);
            const auto num   = integrate<GSL::QAGS>(integrand_num,   q2_min, q2_max);
            const auto denom = integrate<GSL::QAGS>(integrand_denom, q2_min, q2_max);

            return num / denom;
        }
    };

    BToPseudoscalarLeptonNeutrino::BToPseudoscalarLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToPseudoscalarLeptonNeutrino>(new Implementation<BToPseudoscalarLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToPseudoscalarLeptonNeutrino::~BToPseudoscalarLeptonNeutrino()
    {
    }

    //* normalized (|V_Ub|=1) two-fold-distribution, cf. [DSD2014], eq. (13), p. 6
    double
    BToPseudoscalarLeptonNeutrino::normalized_two_differential_decay_width(const double & s, const double & c_theta_l) const
    {
        return _imp->normalized_two_differential_decay_width(s, c_theta_l);
    }

    double
    BToPseudoscalarLeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    BToPseudoscalarLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max);
    }

    // normalized_differential_branching_ratio (|V_Ub|=1)
    double
    BToPseudoscalarLeptonNeutrino::normalized_differential_branching_ratio(const double & s) const
    {
        return _imp->normalized_differential_branching_ratio(s);
    }

    // normalized (|V_Ub|=1) integrated branching_ratio
    double
    BToPseudoscalarLeptonNeutrino::normalized_integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::normalized_differential_branching_ratio,
                                                             _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max);
    }

    double
    BToPseudoscalarLeptonNeutrino::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->numerator_differential_a_fb_leptonic(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    BToPseudoscalarLeptonNeutrino::integrated_lepton_polarization(const double & q2_min, const double & q2_max) const
    {
        return _imp->lepton_polarization(q2_min, q2_max);
    }

    double
    BToPseudoscalarLeptonNeutrino::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::numerator_differential_a_fb_leptonic,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_numerator = integrate<GSL::QAGS>(f, s_min, s_max);
        }

        double integrated_denominator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_denominator = integrate<GSL::QAGS>(f, s_min, s_max);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    BToPseudoscalarLeptonNeutrino::differential_pdf_w(const double & w) const
    {
        return _imp->pdf_w(w);
    }

    double
    BToPseudoscalarLeptonNeutrino::integrated_pdf_w(const double & w_min, const double & w_max) const
    {
        return _imp->integrated_pdf_w(w_min, w_max);
    }

    const std::string
    BToPseudoscalarLeptonNeutrino::description = "\
    The decay B->P l nu, where P=(U qbar), and l=e,mu,tau is a lepton.";

    const std::string
    BToPseudoscalarLeptonNeutrino::kinematics_description_s = "\
    The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    BToPseudoscalarLeptonNeutrino::kinematics_description_c_theta_l = "\
    The cosine of the polar angle theta_l between the charged lepton and the direction opposite to P(seudoscalar) meson in the l-nubar rest frame.";

}
