/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Matthew Kirk
 * Copyright (c) 2025 Danny van Dyk
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
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/s-decays/k-to-pi-l-nu.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <map>
#include <string>

namespace eos
{
    using namespace std::literals::string_literals;

    namespace k_to_pi_l_nu
    {
        struct Amplitudes
        {
                // helicity amplitudes, cf. [DDS:2014A] eqs. 13-14
                complex<double> h_0;
                complex<double> h_t;
                complex<double> h_S;
                complex<double> h_T;
                complex<double> h_tS;
                double          v;
                double          p;
                double          NF;
        };
    } // namespace k_to_pi_l_nu

    template <> struct Implementation<KToPiLeptonNeutrino>
    {
            std::shared_ptr<Model> model;

            Parameters parameters;

            RestrictedOption opt_K;

            UsedParameter m_K;

            UsedParameter tau_K;

            UsedParameter m_pi;

            LeptonFlavorOption opt_l;

            UsedParameter m_l;

            UsedParameter g_fermi;

            UsedParameter hbar;

            const double FF_normalisation_factor;

            UsedParameter mu;

            static const std::vector<OptionSpecification> options;

            GSL::QAGS::Config int_config;

            BooleanOption opt_cp_conjugate;

            std::shared_ptr<FormFactors<VacuumToPP>> form_factors;

            Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
                model(Model::make(o.get("model"_ok, "SM"), p, o)),
                parameters(p),
                opt_K(o, options, "K"_ok),
                m_K(p["mass::" + opt_K.value()], u),
                tau_K(p["life_time::" + opt_K.value()], u),
                m_pi(p["mass::pi^" + std::string(opt_K.value() == "K_u" ? "0" : "-")], u),
                opt_l(o, options, "l"_ok),
                m_l(p["mass::" + opt_l.str()], u),
                g_fermi(p["WET::G_Fermi"], u),
                hbar(p["QM::hbar"], u),
                FF_normalisation_factor(opt_K.value() == "K_L" ? std::sqrt(2.0) : -std::sqrt(2.0)), // K_L -> pi+ vs K_S -> pi+ and K- -> pi^0
                mu(p["us" + opt_l.str() + "nu" + opt_l.str() + "::mu"], u),
                int_config(GSL::QAGS::Config().epsrel(0.5e-3)),
                opt_cp_conjugate(o, options, "cp-conjugate"_ok),
                form_factors(FormFactorFactory<VacuumToPP>::create("0->Kpi::KSvD2025", p, o))
            {
                Context ctx("When constructing K->pilnu observable");

                u.uses(*form_factors);
                u.uses(*model);
            }

            k_to_pi_l_nu::Amplitudes
            amplitudes(const double & s) const
            {
                // NP contributions in EFT including tensor operator (cf. [DDS:2014A]).
                const WilsonCoefficients<ChargedCurrent> wc = model->wet_uslnu(opt_l.value(), opt_cp_conjugate.value());
                const complex<double>                    gV = wc.cvr() + (wc.cvl() - 1.0); // in SM cvl=1 => gV contains NP contribution of cvl
                const complex<double>                    gS = wc.csr() + wc.csl();
                const complex<double>                    gT = wc.ct();

                // form factors
                const complex<double> fp = form_factors->f_p(s) / FF_normalisation_factor;
                const complex<double> f0 = form_factors->f_0(s) / FF_normalisation_factor;
                // const complex<double> fT = form_factors->f_t(s) / FF_normalisation_factor;
                const complex<double> fT = 0; // tensor form factor not implemented yet!

                // running quark masses
                const double msatmu = model->m_s_msbar(mu);
                const double muatmu = model->m_u_msbar(mu);

                const double m_K = this->m_K(), m_K2 = m_K * m_K;
                const double m_pi = this->m_pi(), m_pi2 = m_pi * m_pi;
                const double lam = eos::lambda(m_K2, m_pi2, s);
                const double p   = std::sqrt(lam) / (2.0 * m_K);

                // v = lepton velocity in the dilepton rest frame
                const double m_l    = this->m_l();
                const double v      = (1.0 - m_l * m_l / s);
                const double ml_hat = std::sqrt(1.0 - v);
                const double NF     = v * v * s * power_of<2>(g_fermi()) / (256.0 * power_of<3>(M_PI) * m_K2);

                // helicity amplitudes, cf. [DDS:2014A] eqs. 13-14
                k_to_pi_l_nu::Amplitudes result;

                if (s >= power_of<2>(m_l) && s <= power_of<2>(m_K - m_pi))
                {
                    result.h_0 = 2.0 * m_K * p * fp * (1.0 + gV) / std::sqrt(s);
                    result.h_t = (1.0 + gV) * (m_K2 - m_pi2) * f0 / std::sqrt(s);
                    result.h_S = gS * (m_K2 - m_pi2) * f0 / (msatmu - muatmu);
                    result.h_T = 2.0 * m_K * p * fT * gT / (m_K + m_pi);

                    result.h_tS = result.h_t - result.h_S / ml_hat;

                    result.v  = v;
                    result.p  = p;
                    result.NF = NF;
                }
                else // set amplitudes to zero outside of physical phase space
                {
                    result.h_0 = 0.0;
                    result.h_t = 0.0;
                    result.h_S = 0.0;
                    result.h_T = 0.0;

                    result.h_tS = 0.0;

                    result.v  = 0.99; // avoid NaN in std::sqrt(1.0 - v);
                    result.p  = 0.0;
                    result.NF = 0.0;
                }

                return result;
            }

            // normalized (|V_us| = 1) two-fold-distribution, cf. [DDS:2014A], eq. (12), p. 6
            double
            normalized_two_differential_decay_width(const double & s, const double & c_theta_l) const
            {
                //  d^2 Gamma, cf. [DDS:2014A], p. 6, eq. (13)
                double c_thl_2 = c_theta_l * c_theta_l;
                double s_thl_2 = 1.0 - c_thl_2;
                double c_2_thl = 2.0 * c_thl_2 - 1.0;

                k_to_pi_l_nu::Amplitudes amp(this->amplitudes(s));

                return 2.0 * amp.NF * amp.p
                       * (std::norm(amp.h_0) * s_thl_2 + (1.0 - amp.v) * std::norm(amp.h_0 * c_theta_l - amp.h_tS)
                          + 8.0
                                    * (((2.0 - amp.v) + amp.v * c_2_thl) * std::norm(amp.h_T)
                                       - std::sqrt(1.0 - amp.v) * std::real(amp.h_T * (std::conj(amp.h_0) - std::conj(amp.h_tS) * c_theta_l))));
            }

            // normalized to |V_us = 1|, obtained using cf. [DDS:2014A], eq. (12), agrees with Sakaki'13 et al cf. [STTW:2013A]
            double
            normalized_differential_decay_width(const double & s) const
            {
                k_to_pi_l_nu::Amplitudes amp(this->amplitudes(s));

                return 4.0 / 3.0 * amp.NF * amp.p
                       * (std::norm(amp.h_0) * (3.0 - amp.v) + 3.0 * std::norm(amp.h_tS) * (1.0 - amp.v) + 16.0 * std::norm(amp.h_T) * (3.0 - 2.0 * amp.v)
                          - 24.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_0)));
            }

            double
            normalized_differential_decay_width_p(const double & s) const
            {
                k_to_pi_l_nu::Amplitudes amp(this->amplitudes(s));

                return 4.0 / 3.0 * amp.NF * amp.p * (std::norm(amp.h_0) * (3.0 - amp.v));
            }

            double
            normalized_differential_decay_width_0(const double & s) const
            {
                k_to_pi_l_nu::Amplitudes amp(this->amplitudes(s));

                return 4.0 / 3.0 * amp.NF * amp.p * (3.0 * std::norm(amp.h_t) * (1.0 - amp.v));
            }

            // obtained using cf. [DDS:2014A], eq. (12), defined as int_1^0 d^2Gamma - int_0^-1 d^2Gamma
            // in eq. (12) from cf. [DDS:2014A], (H0 * cos(theta) - HtS)^2 we interpret as |H0 * cos(theta) - HtS|^2
            // crosschecked against [BFNT:2019A] and [STTW:2013A]
            double
            numerator_differential_a_fb_leptonic(const double & s) const
            {
                k_to_pi_l_nu::Amplitudes amp(this->amplitudes(s));

                return -4.0 * amp.NF * amp.p * (std::real(amp.h_0 * std::conj(amp.h_tS)) * (1.0 - amp.v) - 4.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_tS)));
            }

            // obtained using cf. [DDS:2014A], eq. (12) and [BHP:2007A] eq.(1.2)
            double
            numerator_differential_flat_term(const double & s) const
            {
                k_to_pi_l_nu::Amplitudes amp(this->amplitudes(s));

                return amp.NF * amp.p
                       * ((std::norm(amp.h_0) + std::norm(amp.h_tS)) * (1.0 - amp.v) + 16.0 * std::norm(amp.h_T)
                          - 8.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_0)));
            }

            // obtained using cf. [STTW:2013A], eq. (49a - 49b)
            double
            numerator_differential_lepton_polarization(const double & s) const
            {
                k_to_pi_l_nu::Amplitudes amp(this->amplitudes(s));

                const double dGplus = (std::norm(amp.h_0) + 3.0 * std::norm(amp.h_t)) * (1.0 - amp.v) / 2.0 + 3.0 / 2.0 * std::norm(amp.h_S) + 8.0 * std::norm(amp.h_T)
                                      - std::sqrt(1.0 - amp.v) * std::real(3.0 * amp.h_t * std::conj(amp.h_S) + 4.0 * amp.h_0 * std::conj(amp.h_T));
                const double dGminus = std::norm(amp.h_0) + 16.0 * std::norm(amp.h_T) * (1.0 - amp.v) - 8.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_0 * std::conj(amp.h_T));

                return 8.0 / 3.0 * amp.NF * amp.p * (dGplus - dGminus);
            }

            // differential decay width
            double
            differential_decay_width(const double & s) const
            {
                return normalized_differential_decay_width(s) * std::norm(model->ckm_us());
            }

            // differential branching_ratio
            double
            differential_branching_ratio(const double & s) const
            {
                return differential_decay_width(s) * tau_K / hbar;
            }

            // "normalized" (|V_us|=1) differential branching_ratio
            double
            normalized_differential_branching_ratio(const double & s) const
            {
                return normalized_differential_decay_width(s) * tau_K / hbar;
            }

            double
            total_branching_ratio() const
            {
                const double q2_min = power_of<2>(m_l());
                const double q2_max = power_of<2>(m_K() - m_pi());

                std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::differential_branching_ratio, this, std::placeholders::_1);
                return integrate<GSL::QAGS>(f, q2_min, q2_max, int_config);
            }

            double
            pdf_q2(const double & q2) const
            {
                const double q2_min = power_of<2>(m_l());
                const double q2_max = power_of<2>(m_K() - m_pi());

                std::function<double(const double &)> f     = std::bind(&Implementation<KToPiLeptonNeutrino>::normalized_differential_branching_ratio, this, std::placeholders::_1);
                const double                          num   = normalized_differential_branching_ratio(q2);
                const double                          denom = integrate<GSL::QAGS>(f, q2_min, q2_max, int_config);

                return num / denom;
            }

            double
            integrated_pdf_q2(const double & q2_min, const double & q2_max) const
            {
                const double q2_abs_min = power_of<2>(m_l());
                const double q2_abs_max = power_of<2>(m_K() - m_pi());

                std::function<double(const double &)> f     = std::bind(&Implementation<KToPiLeptonNeutrino>::normalized_differential_branching_ratio, this, std::placeholders::_1);
                const double                          num   = integrate<GSL::QAGS>(f, q2_min, q2_max, int_config);
                const double                          denom = integrate<GSL::QAGS>(f, q2_abs_min, q2_abs_max, int_config);

                return num / denom / (q2_max - q2_min);
            }
    };

    const std::vector<OptionSpecification> Implementation<KToPiLeptonNeutrino>::options{
        Model::option_specification(),      FormFactorFactory<PToP>::option_specification(), { "cp-conjugate"_ok,      { "true"s, "false"s }, "false"s },
        {            "l"_ok,            { "e"s, "mu"s },    "mu"s },
               {            "K"_ok, { "K_u"s, "K_S"s, "K_L"s },   "K_u"s },
    };

    KToPiLeptonNeutrino::KToPiLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<KToPiLeptonNeutrino>(new Implementation<KToPiLeptonNeutrino>(parameters, options, *this))
    {
    }

    KToPiLeptonNeutrino::~KToPiLeptonNeutrino() {}

    // normalized (|V_us|=1) two-fold-distribution, cf. [DDS:2014A], eq. (13), p. 6
    double
    KToPiLeptonNeutrino::normalized_two_differential_decay_width(const double & s, const double & c_theta_l) const
    {
        return _imp->normalized_two_differential_decay_width(s, c_theta_l);
    }

    double
    KToPiLeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    KToPiLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::differential_branching_ratio, _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    double
    KToPiLeptonNeutrino::total_branching_ratio() const
    {
        return _imp->total_branching_ratio();
    }

    // normalized_differential_branching_ratio (|V_us|=1)
    double
    KToPiLeptonNeutrino::normalized_differential_branching_ratio(const double & s) const
    {
        return _imp->normalized_differential_branching_ratio(s);
    }

    // normalized (|V_us|=1) integrated branching_ratio
    double
    KToPiLeptonNeutrino::normalized_integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::normalized_differential_branching_ratio, _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    // normalized (|V_us|=1) integrated decay_width
    double
    KToPiLeptonNeutrino::normalized_integrated_decay_width_p(const double & s_min, const double & s_max) const
    {
        std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::normalized_differential_decay_width_p, _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    double
    KToPiLeptonNeutrino::normalized_integrated_decay_width_0(const double & s_min, const double & s_max) const
    {
        std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::normalized_differential_decay_width_0, _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    double
    KToPiLeptonNeutrino::normalized_integrated_decay_width(const double & s_min, const double & s_max) const
    {
        std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::normalized_differential_decay_width, _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    double
    KToPiLeptonNeutrino::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->numerator_differential_a_fb_leptonic(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    KToPiLeptonNeutrino::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::numerator_differential_a_fb_leptonic, _imp.get(), std::placeholders::_1);
            integrated_numerator                    = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        double integrated_denominator;
        {
            std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::normalized_differential_decay_width, _imp.get(), std::placeholders::_1);
            integrated_denominator                  = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    KToPiLeptonNeutrino::differential_flat_term(const double & s) const
    {
        return _imp->numerator_differential_flat_term(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    KToPiLeptonNeutrino::integrated_flat_term(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::numerator_differential_flat_term, _imp.get(), std::placeholders::_1);
            integrated_numerator                    = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        double integrated_denominator;
        {
            std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::normalized_differential_decay_width, _imp.get(), std::placeholders::_1);
            integrated_denominator                  = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    KToPiLeptonNeutrino::differential_lepton_polarization(const double & s) const
    {
        return _imp->numerator_differential_lepton_polarization(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    KToPiLeptonNeutrino::integrated_lepton_polarization(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double(const double &)> f =
                    std::bind(&Implementation<KToPiLeptonNeutrino>::numerator_differential_lepton_polarization, _imp.get(), std::placeholders::_1);
            integrated_numerator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        double integrated_denominator;
        {
            std::function<double(const double &)> f = std::bind(&Implementation<KToPiLeptonNeutrino>::normalized_differential_decay_width, _imp.get(), std::placeholders::_1);
            integrated_denominator                  = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    KToPiLeptonNeutrino::differential_pdf_q2(const double & q2) const
    {
        return _imp->pdf_q2(q2);
    }

    double
    KToPiLeptonNeutrino::integrated_pdf_q2(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_pdf_q2(q2_min, q2_max);
    }

    const std::string KToPiLeptonNeutrino::description = "\
    The decay K->pi l nu, where both K and pi are pseudoscalars, and l=e,mu is a lepton.";

    const std::string KToPiLeptonNeutrino::kinematics_description_q2 = "\
    The invariant mass of the l-nubar pair in GeV^2.";

    const std::string KToPiLeptonNeutrino::kinematics_description_c_theta_l = "\
    The cosine of the polar angle theta_l between the charged lepton and the direction opposite to pi meson in the l-nubar rest frame.";

    const std::set<ReferenceName> KToPiLeptonNeutrino::references{ "S:1982A"_rn, "DDS:2014A"_rn, "STTW:2013A"_rn };

    std::vector<OptionSpecification>::const_iterator
    KToPiLeptonNeutrino::begin_options()
    {
        return Implementation<KToPiLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    KToPiLeptonNeutrino::end_options()
    {
        return Implementation<KToPiLeptonNeutrino>::options.cend();
    }
} // namespace eos
