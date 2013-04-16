/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013 Danny van Dyk
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

#include <eos/rare-b-decays/exclusive-b-to-dilepton.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    template <>
    struct Implementation<BToDilepton>
    {
        std::shared_ptr<Model> model;

        UsedParameter f_B;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter delta_gamma_B;

        UsedParameter mu;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        UsedParameter hbar;

        UsedParameter m_l;

        std::function<complex<double> (const Model *)> lambda;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            f_B(p["decay-constant::B_" + o.get("q", "d")], u),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            tau_B(p["life_time::B_" + o.get("q", "d")], u),
            delta_gamma_B(p["life_time::Delta_B_" + o.get("q", "d")], u),
            mu(p["mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u),
            m_l(p["mass::" + o.get("l", "mu")], u)
        {
            if (o.get("q", "d") == "d")
            {
                lambda = &lambda_t_d;
            }
            else if (o.get("q", "d") == "s")
            {
                lambda = &lambda_t_s;
            }
            else
            {
                // only neutral B mesons can decay in this channel
                throw InternalError("ExclusiveBToDilepton: q = '" + o["q"] + "' is not a valid option for a neutral decay channel");
            }

            u.uses(*model);
        }

        // CKM factors
        static complex<double> lambda_t_d(const Model * model) { return model->ckm_tb() * conj(model->ckm_td()); }
        static complex<double> lambda_t_s(const Model * model) { return model->ckm_tb() * conj(model->ckm_ts()); }

        double branching_ratio_time_zero() const
        {
            double lambda_t = abs(lambda(model.get()));
            double beta_l = std::sqrt(1.0 - 4.0 * power_of<2>(m_l / m_B()));

            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s();

            // cf. [BEKU2002], Eq. (3.6) with c_S,P(') = 0
            return power_of<2>(g_fermi() * alpha_e() * lambda_t) / 64.0 / power_of<3>(M_PI)
                * beta_l * m_B() * power_of<2>(f_B() * 2.0 * m_l) * std::norm(wc.c10() - wc.c10prime()) * tau_B / hbar;
        }

        double branching_ratio_untagged_integrated() const
        {
            double y_s = tau_B() * delta_gamma_B / 2.0;

            // In the case of vanishing C_S(') and C_P(') we find simply
            return branching_ratio_time_zero() / (1.0 - y_s);
        }
    };

    BToDilepton::BToDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDilepton>(new Implementation<BToDilepton>(parameters, options, *this))
    {
    }

    BToDilepton::~BToDilepton()
    {
    }

    double
    BToDilepton::branching_ratio_time_zero() const
    {
        return _imp->branching_ratio_time_zero();
    }

    double
    BToDilepton::branching_ratio_untagged_integrated() const
    {
        return _imp->branching_ratio_untagged_integrated();
    }
}
