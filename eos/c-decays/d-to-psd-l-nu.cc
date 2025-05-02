/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Carolina Bolognani
 * Copyright (c) 2023 Méril Reboud
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

#include <eos/c-decays/d-to-psd-l-nu.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <map>
#include <string>

namespace eos
{
    namespace d_to_psd_l_nu
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
    struct Implementation<DToPseudoscalarLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        Parameters parameters;

        QuarkFlavorOption opt_Q;
        QuarkFlavorOption opt_q;
        SwitchOption opt_I;

        UsedParameter m_D;

        UsedParameter tau_D;

        UsedParameter m_P;

        LeptonFlavorOption opt_l;

        UsedParameter m_l;

        UsedParameter g_fermi;

        UsedParameter hbar;

        const double isospin_factor;

        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        std::function<double (const double &)> m_Q_msbar;
        std::function<complex<double> ()> v_cQ;
        std::function<WilsonCoefficients<ChargedCurrent> (LeptonFlavor, bool)> wc;

        GSL::QAGS::Config int_config;

        SpecifiedOption opt_cp_conjugate;

        bool cp_conjugate;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        // { Q, q, I } -> { process, m_D, m_P, c_I }
        // Q: d, s; the quark flavor in the weak transition
        // q: u, d, s: the spectar quark flavor
        // I: 1, 0, 1/2: the total isospin of the daughter meson
        // process: string that can be used to obtain the form factor
        // D: name of the D meson
        // P: name of the daughter meson
        // c_I: isospin factor by which the amplitudes are multiplied
        static const std::map<std::tuple<QuarkFlavor, QuarkFlavor, std::string>, std::tuple<std::string, std::string, std::string, double>> process_map;

        inline std::string _process() const
        {
            const QuarkFlavor Q = opt_Q.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I = opt_I.value();
            const auto p = process_map.find(std::make_tuple(Q, q, I));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of Q=" + stringify(Q) + ", q=" + stringify(q) + ", I=" + I);

            return std::get<0>(p->second);
        }

        inline std::string _D() const
        {
            const QuarkFlavor Q = opt_Q.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I = opt_I.value();
            const auto p = process_map.find(std::make_tuple(Q, q, I));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of Q=" + stringify(Q) + ", q=" + stringify(q) + ", I=" + I);

            return std::get<1>(p->second);
        }

        inline std::string _P() const
        {
            const QuarkFlavor Q = opt_Q.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I = opt_I.value();
            const auto p = process_map.find(std::make_tuple(Q, q, I));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of Q=" + stringify(Q) + ", q=" + stringify(q) + ", I=" + I);

            return std::get<2>(p->second);
        }

        inline double _isospin_factor() const
        {
            const QuarkFlavor Q = opt_Q.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I = opt_I.value();
            const auto p = process_map.find(std::make_tuple(Q, q, I));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of Q=" + stringify(Q) + ", q=" + stringify(q) + ", I=" + I);

            return std::get<3>(p->second);
        }


        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            parameters(p),
            opt_Q(o, options, "Q"_ok),
            opt_q(o, options, "q"_ok),
            opt_I(o, "I"_ok, { "1", "0", "1/2" }),
            m_D(p["mass::" + _D()], u),
            tau_D(p["life_time::" + _D()], u),
            m_P(p["mass::" + _P()], u),
            opt_l(o, options, "l"_ok),
            m_l(p["mass::" + opt_l.str()], u),
            g_fermi(p["WET::G_Fermi"], u),
            hbar(p["QM::hbar"], u),
            isospin_factor(_isospin_factor()),
            mu(p[opt_Q.str() + "cnu" + opt_l.str() + opt_l.str() + "::mu"], u),
            int_config(GSL::QAGS::Config().epsrel(0.5e-3)),
            opt_cp_conjugate(o, options, "cp-conjugate"_ok),
            cp_conjugate(destringify<bool>(opt_cp_conjugate.value())),
            form_factors(FormFactorFactory<PToP>::create(_process() + "::" + o.get("form-factors"_ok, "BSZ2015"), p, o))
        {
            Context ctx("When constructing D->Plnu observable");

            switch (opt_Q.value())
            {
                case QuarkFlavor::strange:
                    m_Q_msbar = [this](const double & mu) -> double { return model->m_s_msbar(mu); };
                    v_cQ      = [this]() -> complex<double> { return model->ckm_cs(); };
                    wc        = [this](LeptonFlavor l, bool cp) -> WilsonCoefficients<ChargedCurrent> { return model->wet_scnul(l, cp); };
                    break;
                case QuarkFlavor::down:
                    throw InternalError("c -> d nu_l l decays are not yet implemented.");
                default:
                    throw InternalError("Invalid quark flavor: " + stringify(opt_Q.value()));
            }

            u.uses(*form_factors);
            u.uses(*model);
        }

        d_to_psd_l_nu::Amplitudes amplitudes(const double & s) const
        {
            // NP contributions in EFT including tensor operator (cf. [DDS:2014A]).
            const WilsonCoefficients<ChargedCurrent> wc = this->wc(opt_l.value(), cp_conjugate);
            const complex<double> gV = wc.cvr() + (wc.cvl() - 1.0); // in SM cvl=1 => gV contains NP contribution of cvl
            const complex<double> gS = wc.csr() + wc.csl();
            const complex<double> gT = wc.ct();

            // form factors
            const double fp = form_factors->f_p(s);
            const double f0 = form_factors->f_0(s);
            const double fT = form_factors->f_t(s);

            // running quark masses
            const double mcatmu = model->m_c_msbar(mu);
            const double mDatmu = m_Q_msbar(mu);

            const double m_D = this->m_D(), m_D2 = m_D * m_D;
            const double m_P = this->m_P(), m_P2 = m_P * m_P;
            const double lam = eos::lambda(m_D2, m_P2, s);
            const double p = std::sqrt(lam) / (2.0 * m_D);

            // v = lepton velocity in the dilepton rest frame
            const double m_l = this->m_l();
            const double v = (1.0 - m_l * m_l / s);
            const double ml_hat = std::sqrt(1.0 - v);
            const double NF = v * v * s * power_of<2>(g_fermi()) / (256.0 * power_of<3>(M_PI) * m_D2);

            // isospin factor
            const double isospin = this->isospin_factor;

            // helicity amplitudes, cf. [DDS:2014A] eqs. 13-14
            d_to_psd_l_nu::Amplitudes result;

            if (s >= power_of<2>(m_l) && s <= power_of<2>(m_D - m_P))
            {
                result.h_0  =   isospin * 2.0 * m_D * p * fp * (1.0 + gV) / std::sqrt(s);
                result.h_t  =   isospin * (1.0 + gV) * (m_D2 - m_P2) * f0 / std::sqrt(s);
                result.h_S  = - isospin * gS * (m_D2 - m_P2) * f0 / (mcatmu - mDatmu);
                result.h_T  = - isospin * 2.0 * m_D * p * fT * gT / (m_D + m_P);

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

        // normalized (|V_cQ| = 1) two-fold-distribution, cf. [DDS:2014A], eq. (12), p. 6
        double normalized_two_differential_decay_width(const double & s, const double & c_theta_l) const
        {
            //  d^2 Gamma, cf. [DDS:2014A], p. 6, eq. (13)
            double c_thl_2 = c_theta_l * c_theta_l;
            double s_thl_2 = 1.0 - c_thl_2;
            double c_2_thl = 2.0 * c_thl_2 - 1.0;

            d_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return 2.0 * amp.NF * amp.p * (
                       std::norm(amp.h_0) * s_thl_2
                       + (1.0 - amp.v) * std::norm(amp.h_0 * c_theta_l - amp.h_tS)
                       + 8.0 * ( ((2.0 - amp.v) + amp.v * c_2_thl) * std::norm(amp.h_T)
                           - std::sqrt(1.0 - amp.v) * std::real(amp.h_T * (std::conj(amp.h_0) - std::conj(amp.h_tS) * c_theta_l)))
                   );
        }

        // normalized to |V_cQ = 1|, obtained using cf. [DSD:2014A], eq. (12), agrees with Sakaki'13 et al cf. [STTW:2013A]
        double normalized_differential_decay_width(const double & s) const
        {
            d_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return 4.0 / 3.0 * amp.NF * amp.p * (
                       std::norm(amp.h_0) * (3.0 - amp.v)
                       + 3.0 * std::norm(amp.h_tS) * (1.0 - amp.v)
                       + 16.0 * std::norm(amp.h_T) * (3.0 - 2.0 * amp.v)
                       - 24.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_0))
                   );
        }

        double normalized_differential_decay_width_p(const double & s) const
        {
            d_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return 4.0 / 3.0 * amp.NF * amp.p * (
                       std::norm(amp.h_0) * (3.0 - amp.v)
                       );
        }

        double normalized_differential_decay_width_0(const double & s) const
        {
            d_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return 4.0 / 3.0 * amp.NF * amp.p * (
                       3.0 * std::norm(amp.h_t) * (1.0 - amp.v)
                   );
        }

        // obtained using cf. [DDS:2014A], eq. (12), defined as int_1^0 d^2Gamma - int_0^-1 d^2Gamma
        // in eq. (12) from cf. [DDS:2014A], (H0 * cos(theta) - HtS)^2 we interpret as |H0 * cos(theta) - HtS|^2
        // crosschecked against [BFNT:2019A] and [STTW:2013A]
        double numerator_differential_a_fb_leptonic(const double & s) const
        {
            d_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return - 4.0 * amp.NF * amp.p * (
                       std::real(amp.h_0 * std::conj(amp.h_tS)) * (1.0 - amp.v)
                       - 4.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_tS))
                   );
        }

        // obtained using cf. [DDS:2014A], eq. (12) and [BHP2007] eq.(1.2)
        double numerator_differential_flat_term(const double & s) const
        {
            d_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

            return amp.NF * amp.p * (
                       (std::norm(amp.h_0) + std::norm(amp.h_tS)) * (1.0 - amp.v)
                       + 16.0 * std::norm(amp.h_T)
                       - 8.0 * std::sqrt(1.0 - amp.v) * std::real(amp.h_T * std::conj(amp.h_0))
                   );
        }

        // obtained using cf. [STTW2013], eq. (49a - 49b)
        double numerator_differential_lepton_polarization(const double & s) const
        {
            d_to_psd_l_nu::Amplitudes amp(this->amplitudes(s));

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
            return normalized_differential_decay_width(s) * std::norm(v_cQ());
        }

        // differential branching_ratio
        double differential_branching_ratio(const double & s) const
        {
            return differential_decay_width(s) * tau_D / hbar;
        }

        // "normalized" (|V_cQ|=1) differential branching_ratio
        double normalized_differential_branching_ratio(const double & s) const
        {
            return normalized_differential_decay_width(s) * tau_D / hbar;
        }

        double pdf_q2(const double & q2) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = power_of<2>(m_D() - m_P());

            std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::normalized_differential_branching_ratio, this, std::placeholders::_1);
            const double num   = normalized_differential_branching_ratio(q2);
            const double denom = integrate<GSL::QAGS>(f, q2_min, q2_max, int_config);

            return num / denom;
        }

        double pdf_w(const double & w) const
        {
            const double m_D = this->m_D(), m_D2 = m_D * m_D;
            const double m_P = this->m_P(), m_P2 = m_P * m_P;
            const double q2  = m_D2 + m_P2 - 2.0 * m_D * m_P * w;

            return 2.0 * m_D * m_P * pdf_q2(q2);
        }

        double integrated_pdf_q2(const double & q2_min, const double & q2_max) const
        {
            const double q2_abs_min = power_of<2>(m_l());
            const double q2_abs_max = power_of<2>(m_D() - m_P());

            std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::normalized_differential_branching_ratio, this, std::placeholders::_1);
            const double num   = integrate<GSL::QAGS>(f, q2_min,     q2_max,     int_config);
            const double denom = integrate<GSL::QAGS>(f, q2_abs_min, q2_abs_max, int_config);

            return num / denom / (q2_max - q2_min);
        }

        double integrated_pdf_w(const double & w_min, const double & w_max) const
        {
            const double m_D    = this->m_D(), m_D2 = m_D * m_D;
            const double m_P    = this->m_P(), m_P2 = m_P * m_P;
            const double q2_max = m_D2 + m_P2 - 2.0 * m_D * m_P * w_min;
            const double q2_min = m_D2 + m_P2 - 2.0 * m_D * m_P * w_max;

            return integrated_pdf_q2(q2_min, q2_max) * (q2_max - q2_min) / (w_max - w_min);
        }
    };

    const std::map<std::tuple<QuarkFlavor, QuarkFlavor, std::string>, std::tuple<std::string, std::string, std::string, double>>
    Implementation<DToPseudoscalarLeptonNeutrino>::Implementation::process_map
    {
        { { QuarkFlavor::down,    QuarkFlavor::up,      "1"   }, { "D->pi",  "D_u", "pi^+", 1.0                   } },
        { { QuarkFlavor::down,    QuarkFlavor::down,    "1"   }, { "D->pi",  "D_d", "pi^0", 1.0 / std::sqrt(2.0)  } },
        { { QuarkFlavor::down,    QuarkFlavor::strange, "1/2" }, { "D_s->K", "D_s", "K_d",  1.0                   } },
        { { QuarkFlavor::strange, QuarkFlavor::up,      "1/2" }, { "D->K",   "D_u", "K_u",  1.0                   } },
        { { QuarkFlavor::strange, QuarkFlavor::down,    "1/2" }, { "D->K",   "D_d", "K_d",  1.0                   } },
    };

    const std::vector<OptionSpecification>
    Implementation<DToPseudoscalarLeptonNeutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToP>::option_specification(),
        { "cp-conjugate"_ok, { "true", "false" },  "false" },
        { "l"_ok, { "e", "mu", "tau" }, "mu" },
        { "Q"_ok, { "d", "s" }, "s" },
        { "q"_ok, { "u", "d", "s" }, "u" },
        { "I"_ok, { "1", "0", "1/2" }, "1" },
    };

    DToPseudoscalarLeptonNeutrino::DToPseudoscalarLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<DToPseudoscalarLeptonNeutrino>(new Implementation<DToPseudoscalarLeptonNeutrino>(parameters, options, *this))
    {
    }

    DToPseudoscalarLeptonNeutrino::~DToPseudoscalarLeptonNeutrino()
    {
    }

    // normalized (|V_cQ|=1) two-fold-distribution, cf. [DDS:2014A], eq. (13), p. 6
    double
    DToPseudoscalarLeptonNeutrino::normalized_two_differential_decay_width(const double & s, const double & c_theta_l) const
    {
        return _imp->normalized_two_differential_decay_width(s, c_theta_l);
    }

    double
    DToPseudoscalarLeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    DToPseudoscalarLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    // normalized_differential_branching_ratio (|V_cQ|=1)
    double
    DToPseudoscalarLeptonNeutrino::normalized_differential_branching_ratio(const double & s) const
    {
        return _imp->normalized_differential_branching_ratio(s);
    }

    // normalized (|V_cQ|=1) integrated branching_ratio
    double
    DToPseudoscalarLeptonNeutrino::normalized_integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::normalized_differential_branching_ratio,
                                                             _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    // normalized (|V_cQ|=1) integrated decay_width
    double
    DToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width_p(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width_p,
                                                             _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    double
    DToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width_0(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width_0,
                                                             _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    double
    DToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width,
                                                             _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
    }

    double
    DToPseudoscalarLeptonNeutrino::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->numerator_differential_a_fb_leptonic(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    DToPseudoscalarLeptonNeutrino::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::numerator_differential_a_fb_leptonic,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_numerator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        double integrated_denominator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_denominator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    DToPseudoscalarLeptonNeutrino::differential_flat_term(const double & s) const
    {
        return _imp->numerator_differential_flat_term(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    DToPseudoscalarLeptonNeutrino::integrated_flat_term(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::numerator_differential_flat_term,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_numerator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        double integrated_denominator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_denominator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    DToPseudoscalarLeptonNeutrino::differential_lepton_polarization(const double & s) const
    {
        return _imp->numerator_differential_lepton_polarization(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    DToPseudoscalarLeptonNeutrino::integrated_lepton_polarization(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::numerator_differential_lepton_polarization,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_numerator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        double integrated_denominator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<DToPseudoscalarLeptonNeutrino>::normalized_differential_decay_width,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_denominator = integrate<GSL::QAGS>(f, s_min, s_max, _imp->int_config);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    DToPseudoscalarLeptonNeutrino::differential_pdf_q2(const double & q2) const
    {
        return _imp->pdf_q2(q2);
    }

    double
    DToPseudoscalarLeptonNeutrino::differential_pdf_w(const double & w) const
    {
        return _imp->pdf_w(w);
    }

    double
    DToPseudoscalarLeptonNeutrino::integrated_pdf_q2(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_pdf_q2(q2_min, q2_max);
    }

    double
    DToPseudoscalarLeptonNeutrino::integrated_pdf_w(const double & w_min, const double & w_max) const
    {
        return _imp->integrated_pdf_w(w_min, w_max);
    }

    const std::string
    DToPseudoscalarLeptonNeutrino::description = "\
    The decay D->P l nu, where both D=(c qbar) and P=(Q qbar) are pseudoscalars, and l=e,mu,tau is a lepton.";

    const std::string
    DToPseudoscalarLeptonNeutrino::kinematics_description_q2 = "\
    The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    DToPseudoscalarLeptonNeutrino::kinematics_description_w = "\
    The recoil parameter of the D and P states, with w=1 corresponding to zero recoil.";

    const std::string
    DToPseudoscalarLeptonNeutrino::kinematics_description_c_theta_l = "\
    The cosine of the polar angle theta_l between the charged lepton and the direction opposite to P(seudoscalar) meson in the l-nubar rest frame.";

    const std::set<ReferenceName>
    DToPseudoscalarLeptonNeutrino::references
    {
        "S:1982A"_rn,
        "DDS:2014A"_rn,
        "STTW:2013A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    DToPseudoscalarLeptonNeutrino::begin_options()
    {
        return Implementation<DToPseudoscalarLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    DToPseudoscalarLeptonNeutrino::end_options()
    {
        return Implementation<DToPseudoscalarLeptonNeutrino>::options.cend();
    }
}
