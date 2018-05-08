/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2015 Marzia Bordone
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
#include <eos/b-decays/b-to-d-l-nu.hh>
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
    struct Implementation<BToDLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        Parameters parameters;

        SwitchOption opt_q;

        UsedParameter m_B;

        UsedParameter mu;

        UsedParameter tau_B;

        UsedParameter m_D;

        SwitchOption opt_l;

        UsedParameter m_l;

        UsedParameter g_fermi;

        UsedParameter hbar;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            opt_q(o, "q", { "u", "d" }, "d"),
            m_B(p["mass::B_" + opt_q.value()], u),
            mu(p["mu"], u),
            tau_B(p["life_time::B_" + opt_q.value()], u),
            m_D(p["mass::D_" + opt_q.value()], u),
            opt_l(o, "l", {"e", "mu", "tau"}, "mu"),
            m_l(p["mass::" + opt_l.value()], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u)
        {
            form_factors = FormFactorFactory<PToP>::create("B->D@" + o.get("form-factors", "BCL2008"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);
        }

        // normalized to V_cb = 1
        double normalized_differential_decay_width(const double & s) const
        {
            double fp = form_factors->f_p(s);
            double f0 = form_factors->f_0(s);
            // running quark masses
            double mbatmu = model->m_b_msbar(mu);
            double mcatmu = model->m_c_msbar(mu);
            double lam = lambda(m_B * m_B, m_D * m_D, s);
            double p = sqrt(lam) / (2.0 * m_B);
            // make sure we return NaN if s < m_l^2, v = lepton velocity in the dilepton rest frame
            double v = (1.0 - m_l * m_l / s);
            double mvm1 = m_l * m_l / s;
            double v2 = v * v;
            double norm = 8.0 * v2 * m_B * s * power_of<2>(g_fermi())
                / (3.0 * 256.0 * power_of<3>(M_PI * m_B));
            // helicity amplitudes
            const double hh0 = 2.0 * m_B * p * fp / (sqrt(s));
            const double hht = (m_B * m_B - m_D * m_D) * f0 / (sqrt(s));
            const double hhs = (m_B * m_B - m_D * m_D) * f0 / (mbatmu - mcatmu);
            // NP contributions in EFT, cf. e.g. [DBG2013]
            const WilsonCoefficients<BToC> wc = model->wilson_coefficients_b_to_c(opt_l.value(), false);
            const complex<double> gv = wc.cvl() + wc.cvr();
            const complex<double> gs = wc.csl() + wc.csr();
            // normalized(|V_cb|=1) differential decay width including NP (in SM gv=1, and all other couplings are zero)
            return norm * p * (power_of<2>(hh0) * std::norm(gv) * (1.0 + mvm1 / 2.0) + (3.0 * mvm1 / 2.0) * std::norm(hht * gv + hhs * gs / sqrt(mvm1)));
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

    BToDLeptonNeutrino::BToDLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDLeptonNeutrino>(new Implementation<BToDLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToDLeptonNeutrino::~BToDLeptonNeutrino()
    {
    }

    double
    BToDLeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    BToDLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonNeutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max);
    }

    // normalized_differential_branching_ratio (|V_cb|=1)
    double
    BToDLeptonNeutrino::normalized_differential_branching_ratio(const double & s) const
    {
        return _imp->normalized_differential_branching_ratio(s);
    }
    
    // normalized(|Vcb|=1) integrated branching_ratio including NP
    double
    BToDLeptonNeutrino::normalized_integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonNeutrino>::normalized_differential_branching_ratio,
                                                             _imp.get(), std::placeholders::_1);
        
        return integrate<GSL::QAGS>(f, s_min, s_max);
    }

    double
    BToDLeptonNeutrino::differential_r_d(const double & s) const
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
    BToDLeptonNeutrino::integrated_r_d() const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonNeutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        double br_muons;
        {

            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "mu");

            br_muons = integrate<GSL::QAGS>(f, 0.02, 11.62);
        }

        double br_taus;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "tau");

            br_taus = integrate<GSL::QAGS>(f, 3.16, 11.62);
        }

        return br_taus / br_muons;
    }

    const std::string
    BToDLeptonNeutrino::description = "\
The decay B->D l nu, where l=e,mu,tau is a lepton.";

    const std::string
    BToDLeptonNeutrino::kinematics_description_s = "\
The invariant mass of the l-nubar pair in GeV^2.";

}
