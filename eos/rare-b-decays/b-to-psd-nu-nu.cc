/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021-2025 Danny van Dyk
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

#include <eos/rare-b-decays/b-to-psd-nu-nu.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/quantum-numbers.hh>

#include <map>
#include <string>

namespace eos
{
    template <>
    struct Implementation<BToPseudoscalarDineutrino>
    {
        std::shared_ptr<Model> model;

        Parameters parameters;

        LightMesonOption opt_P;
        QuarkFlavorOption opt_q;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_P;

        UsedParameter g_fermi;

        UsedParameter alpha_e;

        UsedParameter hbar;

        const double isospin_factor;

        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        std::function<complex<double> ()> lambda_t;
        std::function<WilsonCoefficients<wc::SBNuNu> ()> wc;

        GSL::QAGS::Config int_config;

        BooleanOption opt_cp_conjugate;

        bool cp_conjugate;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        // { q, P } -> { process, U, B_name, P_name, c_I }
        // q: u, d, s: the spectar quark flavor
        // P: K,eta, eta_prime: the type of daughter meson
        // process: string that can be used to obtain the form factor
        // U: the quark flavor in the weak transition
        // B_name: name of the B meson
        // P_name: name of the daughter meson
        // c_I: isospin factor by which the amplitudes are multiplied
        static const std::map<std::tuple<QuarkFlavor, std::string>, std::tuple<std::string, QuarkFlavor, std::string, std::string, double>> process_map;

        inline std::string _process() const
        {
            const QuarkFlavor q = opt_q.value();
            const auto p = process_map.find(std::make_tuple(q, opt_P.str()));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + opt_P.str());

            return std::get<0>(p->second);
        }

        inline QuarkFlavor _D() const
        {
            const QuarkFlavor q = opt_q.value();
            const auto p = process_map.find(std::make_tuple(q, opt_P.str()));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + opt_P.str());

            return std::get<1>(p->second);
        }

        inline std::string _B() const
        {
            const QuarkFlavor q = opt_q.value();
            const auto p = process_map.find(std::make_tuple(q, opt_P.str()));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + opt_P.str());

            return std::get<2>(p->second);
        }

        inline std::string _P() const
        {
            const QuarkFlavor q = opt_q.value();
            const auto p = process_map.find(std::make_tuple(q, opt_P.str()));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + opt_P.str());

            return std::get<3>(p->second);
        }

        inline double _isospin_factor() const
        {
            const QuarkFlavor q = opt_q.value();
            const auto p = process_map.find(std::make_tuple(q, opt_P.str()));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + opt_P.str());

            return std::get<4>(p->second);
        }


        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            parameters(p),
            opt_P(o, options, "P"_ok),
            opt_q(o, options, "q"_ok),
            m_B(p["mass::" + _B()], u),
            tau_B(p["life_time::" + _B()], u),
            m_P(p["mass::" + _P()], u),
            g_fermi(p["WET::G_Fermi"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            hbar(p["QM::hbar"], u),
            isospin_factor(_isospin_factor()),
            mu(p[stringify(_D()) + "b" + "nunu::mu"], u),
            int_config(GSL::QAGS::Config().epsrel(0.5e-3)),
            opt_cp_conjugate(o, options, "cp-conjugate"_ok),
            cp_conjugate(opt_cp_conjugate.value()),
            form_factors(FormFactorFactory<PToP>::create(_process() + "::" + o.get("form-factors"_ok, "BSZ2015"), p, o))
        {
            Context ctx("When constructing B->Pnunu observables");

            switch (_D())
            {
                case QuarkFlavor::strange:
                    lambda_t  = [*this] () { return model->ckm_tb() * std::conj(model->ckm_ts()); };
                    wc        = [*this] () { return model->wet_sbnunu(cp_conjugate); };
                    break;

                default:
                    throw InternalError("Unexpected quark flavor: '" + stringify(_D()) + "'");
            }

            u.uses(*form_factors);
            u.uses(*model);
        }

        // differential decay width
        double differential_decay_width(const double & q2) const
        {
            const double m_B = this->m_B(), m_B2 = m_B * m_B;
            const double m_P = this->m_P(), m_P2 = m_P * m_P;
            const double m_b = model->m_b_msbar(mu);
            const double m_s = model->m_s_msbar(mu);
            const double lambda = eos::lambda(m_B2, m_P2, q2);

            if (q2 < 0 || q2 >= power_of<2>(m_B - m_P))
                return 0.0;

            double sqrt_lambda = std::sqrt(lambda);

            const auto wc = this->wc();

            const double f_p = form_factors->f_p(q2);
            const double f_0 = form_factors->f_0(q2);
            const double f_t = form_factors->f_t(q2);


            // using different normalization than [FLS:2021A], eq. (1)
            // note that eq. (1) is a Lagrangian, while we use the Hamiltonian definition
            const double norm = power_of<2>(4.0 * g_fermi() * alpha_e / (2.0 * M_PI)) / 2.0 * std::norm(lambda_t())
                // remainder as in [FLS:2021A], eq. (8), except for moving the q2 factor into the square brackets
                * sqrt_lambda / (power_of<3>(4.0 * M_PI * m_B));

            // first term in square brackets in [FLS:2021A], eq. (8)
            const double contr_vector = lambda / 24.0 * f_p * f_p * std::norm(wc.cVL() +  wc.cVR());
            // second line in [FLS:2021A], eq. (8) (ignoring the {bs} Wilson coefficients)
            const double contr_scalar = q2 * power_of<2>((m_B2 - m_P2) / (m_b - m_s)) / 8 * f_0 * f_0 * std::norm(wc.cSL() +  wc.cSR());
            // third line in [FLS:2021A], eq. (8) (ignoring the {bs} Wilson coefficients)
            const double contr_tensor = q2 * 2/3 * lambda / power_of<2>(m_B + m_P) * f_t * f_t * std::norm(wc.cTL());

            // assume the production of 3 diagonal neutrino flavors (nu_i nubar_i)
            return 3.0 * norm * (contr_vector + contr_scalar + contr_tensor);
        }

        // differential branching_ratio
        double differential_branching_ratio(const double & s) const
        {
            return differential_decay_width(s) * tau_B / hbar;
        }
    };

    const std::map<std::tuple<QuarkFlavor, std::string>, std::tuple<std::string, QuarkFlavor, std::string, std::string, double>>
    Implementation<BToPseudoscalarDineutrino>::Implementation::process_map
    {
        { { QuarkFlavor::up,      "K" },         { "B->K",           QuarkFlavor::strange,   "B_u", "K_u",       1.0  } },
        { { QuarkFlavor::down,    "K" },         { "B->K",           QuarkFlavor::strange,   "B_d", "K_d",       1.0  } },
        { { QuarkFlavor::strange, "eta" },       { "B_s->eta",       QuarkFlavor::strange,   "B_s", "eta",       1.0  } },
        { { QuarkFlavor::strange, "eta_prime" }, { "B_s->eta_prime", QuarkFlavor::strange,   "B_s", "eta_prime", 1.0  } },
    };

    const std::vector<OptionSpecification>
    Implementation<BToPseudoscalarDineutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToP>::option_specification(),
        { "cp-conjugate"_ok, { "true", "false" },          "false" },
        { "P"_ok,            { "K", "eta", "eta_prime" },  ""      },
        { "q"_ok,            { "u", "d", "s" },            "u"     },
    };

    BToPseudoscalarDineutrino::BToPseudoscalarDineutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToPseudoscalarDineutrino>(new Implementation<BToPseudoscalarDineutrino>(parameters, options, *this))
    {
    }

    BToPseudoscalarDineutrino::~BToPseudoscalarDineutrino()
    {
    }

    double
    BToPseudoscalarDineutrino::differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_branching_ratio(q2);
    }

    double
    BToPseudoscalarDineutrino::integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToPseudoscalarDineutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, q2_min, q2_max, _imp->int_config);
    }

    const std::string
    BToPseudoscalarDineutrino::description = "\
    The decay B->P nu nu, where both B=(B qbar) and P=(U qbar) are pseudoscalars, and l=e,mu,tau is a lepton.";

    const std::string
    BToPseudoscalarDineutrino::kinematics_description_q2 = "\
    The invariant mass of the nu-nubar pair in GeV^2.";

    const std::set<ReferenceName>
    BToPseudoscalarDineutrino::references
    {
        "BGNS:2014A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BToPseudoscalarDineutrino::begin_options()
    {
        return Implementation<BToPseudoscalarDineutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToPseudoscalarDineutrino::end_options()
    {
        return Implementation<BToPseudoscalarDineutrino>::options.cend();
    }
}
