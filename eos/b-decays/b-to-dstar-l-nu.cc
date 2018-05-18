/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Ahmet Kokulu
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
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/save.hh>

namespace eos
{
    template <>
    struct Implementation<BToDstarLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        Parameters parameters;

        SwitchOption opt_q;

        UsedParameter m_B;
        
        UsedParameter mu;

        UsedParameter tau_B;

        UsedParameter m_Dstar;

        SwitchOption opt_l;

        UsedParameter m_l;

        UsedParameter g_fermi;

        UsedParameter hbar;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
        model(Model::make(o.get("model", "SM"), p, o)),
        parameters(p),
        opt_q(o, "q", { "u", "d" }, "d"),
        m_B(p["mass::B_" + opt_q.value()], u),
        // mu is the renormalization scale
        mu(p["mu"], u),
        tau_B(p["life_time::B_" + opt_q.value()], u),
        m_Dstar(p["mass::D^*_" + opt_q.value()], u),
        opt_l(o, "l", {"e", "mu", "tau"}, "mu"),
        m_l(p["mass::" + opt_l.value()], u),
        g_fermi(p["G_Fermi"], u),
        hbar(p["hbar"], u)
        {
            form_factors = FormFactorFactory<PToV>::create("B->D^*@" + o.get("form-factors", "BSZ2015"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);
        }

        // normalized to V_cb = 1
        double normalized_differential_decay_width(const double & s) const
        {
            // form factors
            double aff0 = form_factors->a_0(s);
            double aff1 = form_factors->a_1(s);
            double aff2 = form_factors->a_2(s);
            double vff  = form_factors->v(s);
            // running quark masses
            double mbatmu = model->m_b_msbar(mu);
            double mcatmu = model->m_c_msbar(mu);
            double lam = lambda(m_B * m_B, m_Dstar * m_Dstar, s);
            double p = sqrt(lam) / (2.0 * m_B);
            // make sure we return NaN if s < m_l^2, vv = lepton velocity in the dilepton rest frame
            double v = (1.0 - m_l * m_l / s);
            double mvm1 = m_l * m_l / s;
            double v2 = v * v;
            double norm = 8.0 * v2 * m_B * s * power_of<2>(g_fermi()) / (3.0 * 256.0 * power_of<3>(M_PI * m_B));
            // define helicity amplitudes
            const double aa0 = ((-4.0 * aff2 * m_B * m_B * p * p) / (m_B + m_Dstar) + aff1 * (m_B + m_Dstar) * (m_B * m_B - m_Dstar * m_Dstar - s)) / (2.0 * m_Dstar * sqrt(s));
            const double aapar = 2.0 * (m_B + m_Dstar) * aff1 / sqrt(2.0);
            const double aaperp = - 4.0 * m_B * p * vff / (sqrt(2.0) * (m_B + m_Dstar));
            const double aat = 2.0 * m_B * p * aff0 / sqrt(s);
            const double aap = - 2.0 * m_B * p * aff0 / (mcatmu + mbatmu);

            // NP contributions in EFT, cf. e.g. [DBG2013]
            const WilsonCoefficients<BToC> wc = model->wilson_coefficients_b_to_c(opt_l.value(), false);
            const complex<double> gv = wc.cvl() + wc.cvr();
            const complex<double> ga = wc.cvl() - wc.cvr();
            const complex<double> gp = wc.csl() - wc.csr();
            // further definitions
            const complex<double> aatP = aat * ga + (1.0 / sqrt(mvm1)) * aap * gp;
            const double aaAV2 = power_of<2>(aa0) * std::norm(ga) + power_of<2>(aapar) * std::norm(ga) + power_of<2>(aaperp) * std::norm(gv);

            // normalized(|V_cb|=1) differential decay width for B->Dstarlnu including NP (in SM gv=1, and all other couplings are zero)
            return norm * p * (aaAV2 + (mvm1 / 2.0) * (aaAV2 + 3.0 * std::norm(aatP)));
        }

        // differential decay width including NP
        double differential_decay_width(const double & s) const
        {
            return normalized_differential_decay_width(s) * std::norm(model->ckm_cb());
        }

        // differential branching_ratio including NP
        double differential_branching_ratio(const double & s) const
        {
            return differential_decay_width(s) * tau_B / hbar;
        }

        // "normalized"(|Vcb|=1) differential branching_ratio including NP
        double normalized_differential_branching_ratio(const double & s) const
        {
            return normalized_differential_decay_width(s) * tau_B / hbar;
        }
    };

    BToDstarLeptonNeutrino::BToDstarLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDstarLeptonNeutrino>(new Implementation<BToDstarLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToDstarLeptonNeutrino::~BToDstarLeptonNeutrino()
    {
    }

    double
    BToDstarLeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    BToDstarLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDstarLeptonNeutrino>::differential_branching_ratio, _imp.get(), std::placeholders::_1);
        return integrate<GSL::QAGS>(f, s_min, s_max);
    }

    // normalized_differential_branching_ratio (|V_cb|=1)
    double
    BToDstarLeptonNeutrino::normalized_differential_branching_ratio(const double & s) const
    {
        return _imp->normalized_differential_branching_ratio(s);
    }

    // normalized(|Vcb|=1) integrated branching_ratio including NP
    double
    BToDstarLeptonNeutrino::normalized_integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDstarLeptonNeutrino>::normalized_differential_branching_ratio, _imp.get(), std::placeholders::_1);
        return integrate<GSL::QAGS>(f, s_min, s_max);
    }

    double
    BToDstarLeptonNeutrino::differential_r_d(const double & s) const
    {
        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "mu");
            br_muons = _imp->differential_branching_ratio(s);
        }

        double br_taus;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "tau");
            br_taus = _imp->differential_branching_ratio(s);
        }
        return br_taus / br_muons;
    }

    double
    BToDstarLeptonNeutrino::integrated_r_d() const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDstarLeptonNeutrino>::differential_branching_ratio, _imp.get(), std::placeholders::_1);

        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "mu");
            // note that upper s limit is now less than the B->Dlnu case
            br_muons = integrate<GSL::QAGS>(f, 0.02, 10.68);
        }

        double br_taus;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "tau");
            // note that upper s limit is now less than the B->Dlnu case
            br_taus = integrate<GSL::QAGS>(f, 3.16, 10.68);
        }
        return br_taus / br_muons;
    }

    const std::string
    BToDstarLeptonNeutrino::description = "\
The decay B->D^* l nu, where l=e,mu,tau is a lepton.";

    const std::string
    BToDstarLeptonNeutrino::kinematics_description_s = "\
The invariant mass of the l-nubar pair in GeV^2.";

}
