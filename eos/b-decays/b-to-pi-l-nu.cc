/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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
#include <eos/b-decays/b-to-pi-l-nu.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    template <>
    struct Implementation<BToPiLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_pi;

        UsedParameter g_fermi;

        UsedParameter hbar;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            tau_B(p["life_time::B_" + o.get("q", "d")], u),
            m_pi(p["mass::pi^" + std::string(o.get("q", "d") == "d" ? "+" : "0")], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u)
        {
            if (o.get("l", "mu") == "tau")
            {
                throw InternalError("BToPiLeptonNeutrino: l == 'tau' is not a valid option for this decay channel");
            }

            if ((o.get("q", "d") != "d") && (o.get("q", "u") != "u"))
            {
                // only B_{d,u} mesons can decay in this channel
                throw InternalError("BToPiLeptonNeutrino: q = '" + o["q"] + "' is not a valid option for this decay channel");
            }

            form_factors = FormFactorFactory<PToP>::create("B->pi@" + o.get("form-factors", "BCL2008"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);
        }

        double differential_branching_ratio(const double & s) const
        {
            // cf. e.g. [BCL2008], eq. (2), p. 1
            double fp = form_factors->f_p(s);
            double lam = lambda(m_B * m_B, m_pi * m_pi, s);
            double norm = power_of<2>(g_fermi * std::abs(model->ckm_ub()))
                / (192.0 * power_of<3>(M_PI * m_B));

            return norm * lam * std::sqrt(lam) * fp * fp * tau_B / hbar;
        }

        double differential_zeta(const double & s) const
        {
            // cf. e.g. [BCL2008], eq. (2), p. 1
            double fp = form_factors->f_p(s);
            double lam = lambda(m_B * m_B, m_pi * m_pi, s);
            double norm = power_of<2>(g_fermi())
                / (192.0 * power_of<3>(M_PI * m_B()));

            return norm * lam * std::sqrt(lam) * fp * fp / hbar;
        }
    };

    BToPiLeptonNeutrino::BToPiLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToPiLeptonNeutrino>(new Implementation<BToPiLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToPiLeptonNeutrino::~BToPiLeptonNeutrino()
    {
    }

    double
    BToPiLeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    BToPiLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToPiLeptonNeutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        return integrate(f, 64, s_min, s_max);
    }

    double
    BToPiLeptonNeutrino::differential_zeta(const double & s) const
    {
        return _imp->differential_zeta(s);
    }

    double
    BToPiLeptonNeutrino::integrated_zeta(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToPiLeptonNeutrino>::differential_zeta,
                _imp.get(), std::placeholders::_1);

        return integrate(f, 64, s_min, s_max);
    }
}
