/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
 * Copyright (c) 2014      Frederik Beaujean
 * Copyright (c) 2014      Christoph Bobeth
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

#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/rare-b-decays/b-to-ll.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    using namespace std::literals::string_literals;

    template <>
    struct Implementation<BToDilepton>
    {
        std::shared_ptr<Model> model;

        LeptonFlavorOption opt_l;
        QuarkFlavorOption opt_q;

        UsedParameter f_B;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter delta_gamma_B;

        UsedParameter mu;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        UsedParameter hbar;

        UsedParameter m_l;

        UsedParameter m_b;

        UsedParameter m_q;

        static const std::vector<OptionSpecification> options;

        std::function<complex<double> (const Model *)> lambda;

        using xi_t = std::array<complex<double>, 4>;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            opt_l(o, options, "l"_ok),
            opt_q(o, options, "q"_ok),
            f_B(p["decay-constant::B_" + opt_q.str()], u),
            m_B(p["mass::B_" + opt_q.str()], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            delta_gamma_B(p["life_time::Delta_B_" + opt_q.str()], u),
            mu(p["sb" + opt_l.str() + opt_l.str() + "::mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["WET::G_Fermi"], u),
            hbar(p["QM::hbar"], u),
            m_l(p["mass::" + opt_l.str()], u),
            m_b(p["mass::b(MSbar)"], u),
            m_q(p["mass::" + opt_q.str() + "(2GeV)"], u)
        {
            Context ctx("When constructing B->ll observables");

            switch (opt_q.value())
            {
                case QuarkFlavor::strange:
                    lambda = &lambda_t_s;
                    break;
                default:
                    // only neutral B mesons can decay in this channel
                    throw InternalError("ExclusiveBToDilepton: q = '" + opt_q.str() + "' is not a valid option for a neutral decay channel");
            }

            u.uses(*model);
        }

        // CKM factors
        static complex<double> lambda_t_d(const Model * model) { return model->ckm_tb() * conj(model->ckm_td()); }
        static complex<double> lambda_t_s(const Model * model) { return model->ckm_tb() * conj(model->ckm_ts()); }

        xi_t calc_amplitudes() const
        {
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), opt_l.value());

            double factor = power_of<2>(m_B()) / 2.0 / m_l / (m_b + m_q);
            complex<double> S = std::sqrt(1.0 - 4.0 * power_of<2>(m_l / m_B)) * factor * (wc.cS() - wc.cSprime());
            complex<double> P = (wc.c10() - wc.c10prime()) + factor * (wc.cP() - wc.cPprime());

            xi_t xi;
            xi[0] = -1.0 * (P + S) / std::conj(S - P);
            xi[1] = -1.0 * (S - P) / std::conj(P + S);
            xi[2] = std::norm(P) + std::norm(S);
            xi[3] = P * P - S * S;

            return xi;
        }

        double y_q() const
        {
            return tau_B() * delta_gamma_B / 2.0;
        }

        // cf. [BEKU2002], Eq. (3.6)
        double branching_ratio_time_zero() const
        {
            double lambda_t = abs(lambda(model.get()));
            double beta_l = std::sqrt(1.0 - 4.0 * power_of<2>(m_l / m_B()));

            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), opt_l.value());

            return power_of<2>(g_fermi() * alpha_e() * lambda_t * f_B()) / 64.0 / power_of<3>(M_PI) * tau_B / hbar
                * beta_l * power_of<3>(m_B()) * (
                        power_of<2>(beta_l) * std::norm(m_B / (m_b + m_q) * (wc.cS() - wc.cSprime()))
                        + std::norm(m_B / (m_b + m_q) * (wc.cP() - wc.cPprime()) + 2.0 * m_l / m_B * (wc.c10() - wc.c10prime())));
        }

        // [F2012], Eq. (29), (30)
        double branching_ratio_untagged_integrated() const
        {
            xi_t xi = calc_amplitudes();
            double factor = power_of<2>(g_fermi() * alpha_e() * f_B * 2.0 * m_l) * tau_B / hbar * m_B * std::norm(lambda(model.get())) / (64.0 * power_of<3>(M_PI))
                            * std::sqrt(1 - 4* power_of<2>(m_l/m_B));
            return factor / (1.0 - power_of<2>(y_q())) * (std::real(xi[2]) + std::real(xi[3]) * y_q());
        }

        // [F2012], Eq. (25)
        double cp_asymmetry_del_gamma() const
        {
            xi_t xi = calc_amplitudes();
            return 2.0 * std::real(xi[0]) / (1.0 + std::norm(xi[0]));
        }

        // [F2012], Eq. (24)
        double cp_asymmetry_mixing_S() const
        {
            xi_t xi = calc_amplitudes();
            return 2.0 * std::imag(xi[0]) / (1.0 + std::norm(xi[0]));
        }

        // [F2012], Eq. (8)
        double effective_lifetime() const
        {
            const double cp_asym = cp_asymmetry_del_gamma();
            const double y = y_q();

            return tau_B() / hbar / (1.0 - power_of<2>(y)) * (1.0 + 2.0 * cp_asym * y + power_of<2>(y)) / (1.0 + cp_asym * y);
        }
    };

    BToDilepton::BToDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDilepton>(new Implementation<BToDilepton>(parameters, options, *this))
    {
    }

    BToDilepton::~BToDilepton()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<BToDilepton>::options
    {
        Model::option_specification(),
        {"l"_ok, { "e"s, "mu"s, "tau"s }, "mu"s},
        {"q"_ok, { "s"s }, "s"s}
    };

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

    double
    BToDilepton::cp_asymmetry_del_gamma() const
    {
        return _imp->cp_asymmetry_del_gamma();
    }

    double
    BToDilepton::cp_asymmetry_mixing_S() const
    {
        return _imp->cp_asymmetry_mixing_S();
    }

    double
    BToDilepton::effective_lifetime() const
    {
        return _imp->effective_lifetime();
    }

    const std::set<ReferenceName>
    BToDilepton::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BToDilepton::begin_options()
    {
        return Implementation<BToDilepton>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToDilepton::end_options()
    {
        return Implementation<BToDilepton>::options.cend();
    }
}
