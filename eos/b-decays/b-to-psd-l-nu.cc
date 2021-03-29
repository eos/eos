/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2015 Marzia Bordone
 * Copyright (c) 2018, 2019 Ahmet Kokulu
 * Copyright (c) 2021 Christoph Bobeth
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
#include <eos/utils/destringify.hh>
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
            // helicity amplitudes, cf. [DDS:2014A] eqs. 13-14
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

        GSL::QAGS::Config int_config;

        bool cp_conjugate;

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
            hbar(p["QM::hbar"], u),
            int_config(GSL::QAGS::Config().epsrel(0.5e-3)),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false")))
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
            // NP contributions in EFT including tensor operator (cf. [DDS:2014A]).
            auto wc = this->wc(opt_l.value(), cp_conjugate);
            const complex<double> gV = wc.cvr() + (wc.cvl() - 1.0); // in SM cvl=1 => gV contains NP contribution of cvl
            const complex<double> gS = wc.csr() + wc.csl();
            const complex<double> gT = wc.ct();

            // form factors
            const double fp = form_factors->f_p(s);
            const double f0 = form_factors->f_0(s);
            const double fT = form_factors->f_t(s);

            // running quark masses
            const double mbatmu = model->m_b_msbar(mu);
            const double mUatmu = m_U_msbar(mu);

            const double m_B = this->m_B(), m_B2 = m_B * m_B;
            const double m_P = this->m_P(), m_P2 = m_P * m_P;
            const double lam = eos::lambda(m_B2, m_P2, s);
            const double p = std::sqrt(lam) / (2.0 * m_B);

            // v = lepton velocity in the dilepton rest frame
            const double m_l = this->m_l();
            const double v = (1.0 - m_l * m_l / s);
            const double ml_hat = std::sqrt(1.0 - v);
            const double NF = v * v * s * power_of<2>(g_fermi()) / (256.0 * power_of<3>(M_PI) * m_B2);

            // helicity amplitudes, cf. [DDS:2014A] eqs. 13-14
            b_to_psd_l_nu::Amplitudes result;

            if (s >= power_of<2>(m_l) && s <= power_of<2>(m_B - m_P))
            {
                result.h_0  = 2.0 * m_B * p * fp * (1.0 + gV) / std::sqrt(s);
                result.h_t  = (1.0 + gV) * (m_B2 - m_P2) * f0 / std::sqrt(s);
                result.h_S  = - gS * (m_B2 - m_P2) * f0 / (mbatmu - mUatmu);
                result.h_T  = - 2.0 * m_B * p * fT * gT / (m_B + m_P);
                result.h_tS = result.h_t - result.h_S / ml_hat;

                result.v  = v;
                result.p  = p;
                result.NF = NF;
            }
            else // set amplitudes to zero outside of physical phase space
            {
                result.h_0  = 0.0;
                result.h_t  = 0.0;
                result.h_S  = 0.0;
                result.h_T  = 0.0;
                result.h_tS = 0.0;

                result.v  = 0.99; // avoid NaN in std::sqrt(1.0 - v);
                result.p  = 0.0;
                result.NF = 0.0;
            }

            return result;
        }

        // normalized (|V_Ub| = 1) two-fold-distribution, cf. [DDS:2014A], eq. (12), p. 6
        double normalized_two_differential_decay_width(const double & s, const double & c_theta_l) const
        {
            //  d^2 Gamma, cf. [DDS:2014A], p. 6, eq. (13)
            double c_thl_2 = c_theta_l * c_theta_l;
            double s_thl_2 = 1.0 - c_thl_2;
            double c_2_thl = 2.0 * c_thl_2 - 1.0;

            b_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return 2.0 * amp.NF * amp.p * (
                       std::norm(amp.h_0) * s_thl_2
                       + (1.0 - amp.v) * std::norm(amp.h_0 * c_theta_l - amp.h_tS)
                       + 8.0 * ( ((2.0 - amp.v) + amp.v * c_2_thl) * std::norm(amp.h_T)
                           - std::sqrt(1.0 - amp.v) * std::real(amp.h_T * (std::conj(amp.h_0) - std::conj(amp.h_tS) * c_theta_l)))
                   );
        }

        // normalized to |V_Ub = 1|, obtained using cf. [DSD:2014A], eq. (12), agrees with Sakaki'13 et al cf. [STTW:2013A]
        double normalized_differential_decay_width(const double & s) const
        {
            b_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return 4.0 / 3.0 * amp.NF * amp.p * (
                       std::norm(amp.h_0) * (3.0 - amp.v)
                       + 3.0 * std::norm(amp.h_tS) * (1.0 - amp.v)
                       + 16.0 * std::norm(amp.h_T) * (3.0 - 2.0 * amp.v)
                       - 24.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_0))
                   );
        }

        double normalized_differential_decay_width_p(const double & s) const
        {
            b_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return 4.0 / 3.0 * amp.NF * amp.p * (
                       std::norm(amp.h_0) * (3.0 - amp.v)
                       );
        }

        double normalized_differential_decay_width_0(const double & s) const
        {
            b_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return 4.0 / 3.0 * amp.NF * amp.p * (
                       3.0 * std::norm(amp.h_t) * (1.0 - amp.v)
                   );
        }

        // obtained using cf. [DDS:2014A], eq. (12), defined as int_1^0 d^2Gamma - int_0^-1 d^2Gamma
        // in eq. (12) from cf. [DDS:2014A], (H0 * cos(theta) - HtS)^2 we interpret as |H0 * cos(theta) - HtS|^2
        // crosschecked against [BFNT:2019A] and [STTW:2013A]
        double numerator_differential_a_fb_leptonic(const double & s) const
        {
            b_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return - 4.0 * amp.NF * amp.p * (
                       std::real(amp.h_0 * std::conj(amp.h_tS)) * (1.0 - amp.v)
                       - 4.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_tS))
                   );
        }

        // obtained using cf. [DDS:2014A], eq. (12) and [BHP2007] eq.(1.2)
        double numerator_differential_flat_term(const double & s) const
        {
            b_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return amp.NF * amp.p * (
                       (std::norm(amp.h_0) + std::norm(amp.h_tS)) * (1.0 - amp.v)
                       + 16.0 * std::norm(amp.h_T)
                       - 8.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_0))
                   );
        }

        // obtained using cf. [STTW2013], eq. (49a - 49b)
        double numerator_differential_lepton_polarization(const double & s) const
        {
            b_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            const double dGplus = (std::norm(amp.h_0) + 3.0 * std::norm(amp.h_t)) * (1.0 - amp.v) / 2.0
                                + 3.0 / 2.0 * std::norm(amp.h_S)
                                + 8.0 * std::norm(amp.h_T)
                                - std::sqrt(1.0 - amp.v) * std::real(3.0 * amp.h_t * std::conj(amp.h_S)
                                                                   + 4.0 * amp.h_0 * std::conj(amp.h_T));
            const double dGminus = std::norm(amp.h_0)
                                 + 16.0 * std::norm(amp.h_T) * (1.0 - amp.v)
                                 - 8.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_0 * std::conj(amp.h_T));

            return 8.0 / 3.0 * amp.NF * amp.p * (dGplus - dGminus);
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
            const double denom = integrate<GSL::QAGS>(f, q2_min, q2_max, int_config);

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
            const double num   = integrate<GSL::QAGS>(f, q2_min,     q2_max,     int_config);
            const double denom = integrate<GSL::QAGS>(f, q2_abs_min, q2_abs_max, int_config);

            return num / denom / (q2_max - q2_min);
        }

        double integrated_pdf_w(const double & w_min, const double & w_max) const
        {
            const double m_B    = this->m_B(), m_B2 = m_B * m_B;
            const double m_P    = this->m_P(), m_P2 = m_P * m_P;
            const double q2_max = m_B2 + m_P2 - 2.0 * m_B * m_P * w_min;
            const double q2_min = m_B2 + m_P2 - 2.0 * m_B * m_P * w_max;

            return integrated_pdf_q2(q2_min, q2_max) * (q2_max - q2_min) / (w_max - w_min);
        }
    };

    BToPseudoscalarLeptonNeutrino::BToPseudoscalarLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToPseudoscalarLeptonNeutrino>(new Implementation<BToPseudoscalarLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToPseudoscalarLeptonNeutrino::~BToPseudoscalarLeptonNeutrino()
    {
    }

    // normalized (|V_Ub|=1) two-fold-distribution, cf. [DDS:2014A], eq. (13), p. 6
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

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
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

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    // normalized (|V_Ub|=1) integrated decay_width
    double
    BToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width_p(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width_p,
                                                             _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    double
    BToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width_0(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width_0,
                                                             _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    double
    BToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width,
                                                             _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    double
    BToPseudoscalarLeptonNeutrino::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->numerator_differential_a_fb_leptonic(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    BToPseudoscalarLeptonNeutrino::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::numerator_differential_a_fb_leptonic,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_numerator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        double integrated_denominator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_denominator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    BToPseudoscalarLeptonNeutrino::differential_flat_term(const double & s) const
    {
        return _imp->numerator_differential_flat_term(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    BToPseudoscalarLeptonNeutrino::integrated_flat_term(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::numerator_differential_flat_term,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_numerator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        double integrated_denominator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_denominator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    BToPseudoscalarLeptonNeutrino::differential_lepton_polarization(const double & s) const
    {
        return _imp->numerator_differential_lepton_polarization(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    BToPseudoscalarLeptonNeutrino::integrated_lepton_polarization(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::numerator_differential_lepton_polarization,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_numerator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        double integrated_denominator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_denominator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    BToPseudoscalarLeptonNeutrino::differential_pdf_q2(const double & q2) const
    {
        return _imp->pdf_q2(q2);
    }

    double
    BToPseudoscalarLeptonNeutrino::differential_pdf_w(const double & w) const
    {
        return _imp->pdf_w(w);
    }

    double
    BToPseudoscalarLeptonNeutrino::integrated_pdf_q2(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_pdf_q2(q2_min, q2_max);
    }

    double
    BToPseudoscalarLeptonNeutrino::integrated_pdf_w(const double & w_min, const double & w_max) const
    {
        return _imp->integrated_pdf_w(w_min, w_max);
    }

    const std::string
    BToPseudoscalarLeptonNeutrino::description = "\
    The decay B->P l nu, where both B=(B qbar) and P=(U qbar) are pseudoscalars, and l=e,mu,tau is a lepton.";

    const std::string
    BToPseudoscalarLeptonNeutrino::kinematics_description_q2 = "\
    The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    BToPseudoscalarLeptonNeutrino::kinematics_description_w = "\
    The recoil parameter of the B and P states, with w=1 corresponding to zero recoil.";

    const std::string
    BToPseudoscalarLeptonNeutrino::kinematics_description_c_theta_l = "\
    The cosine of the polar angle theta_l between the charged lepton and the direction opposite to P(seudoscalar) meson in the l-nubar rest frame.";

    const std::set<ReferenceName>
    BToPseudoscalarLeptonNeutrino::references
    {
        "S:1982A"_rn,
        "DDS:2014A"_rn,
        "STTW:2013A"_rn
    };
}
