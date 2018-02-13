/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2015 Marzia Bordone
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

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_D;

        UsedParameter m_l;

        UsedParameter g_fermi;

        UsedParameter hbar;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            tau_B(p["life_time::B_" + o.get("q", "d")], u),
            m_D(p["mass::D^" + std::string(o.get("q", "d") == "d" ? "+" : "0")], u),
            m_l(p["mass::" + o.get("l", "mu")], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u)
        {
            if ((o.get("q", "d") != "d") && (o.get("q", "d") != "u")) // q = d is the default
            {
                // only B_{d,u} mesons can decay in this channel
                throw InternalError("BToDLeptonNeutrino: q = '" + o["q"] + "' is not a valid option for this decay channel");
            }

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
            double lam = lambda(m_B * m_B, m_D * m_D, s);
            double norm = power_of<2>(g_fermi())
                / (384.0 * power_of<3>(M_PI * m_B));
            // make sure we return NaN if s < m_l^2
            double sqrtv = sqrt(1.0 - m_l * m_l / s);
            double v = sqrtv * sqrtv, v2 = v * v;

            // correct result
            return norm * sqrt(lam) * v2 * ((3.0 - v) * fp * fp * lam + 3.0 * f0 * f0 * (1.0 - v) * power_of<2>(m_B * m_B - m_D * m_D));
        }

        double differential_decay_width(const double & s) const
        {
            return normalized_differential_decay_width(s) * std::norm(model->ckm_cb());
        }

        double differential_branching_ratio(const double & s) const
        {
            return differential_decay_width(s) * tau_B / hbar;
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

        return integrate1D(f, 128, s_min, s_max);
    }

    double
    BToDLeptonNeutrino::differential_r_d(const double & s) const
    {
        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, 0 /*_imp->parameters["mass::mu"]()*/);
            br_muons = _imp->differential_branching_ratio(s);
        }

        double br_taus;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
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
            Save<Parameter, double> save_m_l(_imp->m_l, 0 /*_imp->parameters["mass::mu"]()*/);
            br_muons = integrate1D(f, 128, 0.02, 11.62);
        }

        double br_taus;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            br_taus = integrate1D(f, 128, 3.16, 11.62);
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
