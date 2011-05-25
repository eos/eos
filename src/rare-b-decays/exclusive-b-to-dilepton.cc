/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <src/rare-b-decays/exclusive-b-to-dilepton.hh>
#include <src/utils/destringify.hh>
#include <src/utils/model.hh>
#include <src/utils/power_of.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    template <>
    struct Implementation<BToDilepton>
    {
        std::shared_ptr<Model> model;

        UsedParameter abs_c10;

        UsedParameter arg_c10;

        UsedParameter abs_c10prime;

        UsedParameter arg_c10prime;

        UsedParameter f_B;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter mu;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        UsedParameter hbar;

        double m_l;

        std::function<complex<double> (const Model *)> lambda;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            abs_c10(p["Abs{c10}"], u),
            arg_c10(p["Arg{c10}"], u),
            abs_c10prime(p["Abs{c10'}"], u),
            arg_c10prime(p["Arg{c10'}"], u),
            f_B(p["decay-constant::B_" + o.get("q", "d")], u),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            tau_B(p["life_time::B_" + o.get("q", "d")], u),
            mu(p["mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u)
        {
            static const double m_mu = 0.10565836; // (GeV), cf. [PDG2008], p. 13
            static const double m_e = 0.00051099892; // (GeV), cf. [PDG2008], p. 13

            m_l = m_mu;

            if (o.has("l") && ("e" == o["l"]))
            {
                m_l = m_e;
            }

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

        inline complex<double> c10() const { return std::polar(abs_c10(), arg_c10()); }
        inline complex<double> c10prime() const { return std::polar(abs_c10prime(), arg_c10prime()); }

        double branching_ratio() const
        {
            double lambda_t = abs(lambda(model.get()));
            double beta_l = std::sqrt(1.0 - 4.0 * power_of<2>(m_l / m_B()));

            // cf. [BEKU2002], Eq. (3.6) with c_S,P(') = 0
            return power_of<2>(g_fermi() * alpha_e() * lambda_t) / 64.0 / power_of<3>(M_PI)
                * beta_l * m_B() * power_of<2>(f_B() * 2.0 * m_l) * std::norm(c10() - c10prime()) * tau_B / hbar;
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
    BToDilepton::branching_ratio() const
    {
        return _imp->branching_ratio();
    }
}
