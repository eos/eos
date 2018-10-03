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
            double lam = lambda(m_B * m_B, m_Dstar * m_Dstar, s);
            double p = sqrt(lam) / (2.0 * m_B);
            // make sure we return NaN if s < m_l^2, vv = lepton velocity in the dilepton rest frame
            double v = (1.0 - m_l * m_l / s);
            double mvm1 = m_l * m_l / s;
            double v2 = v * v;
            double norm = 8.0 * v2 * m_B * s * power_of<2>(g_fermi()) / (3.0 * 256.0 * power_of<3>(M_PI * m_B));
            // helicity amplitudes, cf. e.g. Sakaki:2013bfa
            const double hhVp = (m_B + m_Dstar) * aff1 - 2.0 * m_B * p * vff / (m_B + m_Dstar);
            const double hhVm = (m_B + m_Dstar) * aff1 + 2.0 * m_B * p * vff / (m_B + m_Dstar);
            const double hhV0 = -((m_B + m_Dstar) / (2.0 * m_Dstar * sqrt(s))) * ((m_B * m_B - m_Dstar * m_Dstar - s) * aff1 - (4.0 * m_B * m_B * p * p * aff2) / ((m_B + m_Dstar) * (m_B + m_Dstar)));
            const double hhVt = -2.0 * m_B * p * aff0 / sqrt(s);
            const double hhS  = -2.0 * m_B * p * aff0 / (mcatmu + mbatmu);
            const double hhTp =  (1.0 / sqrt(s)) * ((m_B * m_B - m_Dstar * m_Dstar) * tff2 + 2.0 * m_B * p * tff1);
            const double hhTm =  (1.0 / sqrt(s)) * (-(m_B * m_B - m_Dstar * m_Dstar) * tff2 + 2.0 * m_B * p * tff1);
            const double hhT0 =  (1.0 / (2.0 * m_Dstar)) * (-(m_B * m_B + 3.0 * m_Dstar * m_Dstar - s) * tff2 + 4.0 * m_B * m_B * p * p * tff3 / (m_B * m_B - m_Dstar * m_Dstar));

            // NP contributions in EFT including tensor operator
            const WilsonCoefficients<BToC> wc = model->wilson_coefficients_b_to_c(opt_l.value(), false);
            const complex<double> cv1 = wc.cvl() - 1.0;
            const complex<double> cv2 = wc.cvr();
            const complex<double> cs1 = wc.csr();
            const complex<double> cs2 = wc.csl();
            const complex<double> cT  = wc.ct();

            // normalized(|V_cb|=1) differential decay width for B->Dstarlnu including tensor operators cont. (in SM cvl=1, and all other couplings are zero), cf. e.g. Sakaki:2013bfa
            return norm * p * ((std::norm(1.0 + cv1) + std::norm(cv2)) * ((1.0 + mvm1 / 2.0) * (power_of<2>(hhVp) + power_of<2>(hhVm) + power_of<2>(hhV0)) + (3.0 * mvm1 / 2.0) * power_of<2>(hhVt)) - 2.0 * std::real((1.0 + cv1) * std::conj(cv2)) * ((1.0 + mvm1 / 2.0) * (power_of<2>(hhV0) + 2.0 * hhVp * hhVm) + (3.0 * mvm1 / 2.0) * power_of<2>(hhVt)) + (3.0 / 2.0) * std::norm(cs1 - cs2) * power_of<2>(hhS) + 8.0 * std::norm(cT) * (1.0 + 2.0 * mvm1) * (power_of<2>(hhTp) + power_of<2>(hhTm) + power_of<2>(hhT0)) + 3.0 * std::real((1.0 + cv1 - cv2) * (std::conj(cs1) - std::conj(cs2))) * sqrt(mvm1) * hhS * hhVt - 12.0 * std::real((1.0 + cv1) * std::conj(cT)) * sqrt(mvm1) * (hhT0 * hhV0 + hhTp * hhVp - hhTm * hhVm) + 12.0 * std::real(cv2 * std::conj(cT)) * sqrt(mvm1) * (hhT0 * hhV0 + hhTp * hhVm - hhTm * hhVp));
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
    BToDstarLeptonNeutrino::integrated_r_d(const double & s_min_mu, const double & s_min_tau, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDstarLeptonNeutrino>::differential_branching_ratio, _imp.get(), std::placeholders::_1);

        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "mu");
            // note that upper s limit is now less than the B->Dlnu case
            br_muons = integrate<GSL::QAGS>(f, s_min_mu, s_max);
        }

        double br_taus;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "tau");
            // note that upper s limit is now less than the B->Dlnu case
            br_taus = integrate<GSL::QAGS>(f, s_min_tau, s_max);
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
