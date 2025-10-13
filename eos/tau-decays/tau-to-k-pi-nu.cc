/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Danny van Dyk
 * Copyright (c) 2025 Matthew Kirk
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
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/tau-decays/tau-to-k-pi-nu.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    using std::norm;

    /*
     * Decay: tau^+ -> [K pi]^+ nubar, cf. [CCH:2017A]
     */
    template <> struct Implementation<TauToKPiNeutrino>
    {
            std::shared_ptr<FormFactors<VacuumToPP>> form_factors;

            SpecifiedOption opt_model;

            std::shared_ptr<Model> model;

            UsedParameter hbar;

            UsedParameter g_fermi;

            UsedParameter m_tau;

            RestrictedOption opt_K;

            UsedParameter m_K;

            UsedParameter m_pi;

            UsedParameter m_s;

            UsedParameter m_u;

            UsedParameter tau_tau;

            const double isospin_factor;

            UsedParameter mu;

            GSL::QAGS::Config int_config;

            static const std::vector<OptionSpecification> options;

            Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
                form_factors(FormFactorFactory<VacuumToPP>::create("0->Kpi::KSvD2025", p, o)),
                opt_model(o, options, "model"_ok),
                opt_K(o, options, "K"_ok),
                model(Model::make(opt_model.value(), p, o)),
                hbar(p["QM::hbar"], u),
                g_fermi(p["WET::G_Fermi"], u),
                m_tau(p["mass::tau"], u),
                m_K(p["mass::" + opt_K.value()], u),
                m_pi(p["mass::pi^" + std::string(opt_K.value() == "K_u" ? "0" : "-")], u),
                m_s(p["mass::s(2GeV)"], u),
                m_u(p["mass::u(2GeV)"], u),
                tau_tau(p["life_time::tau"], u),
                isospin_factor(opt_K.value() == "K_L" ? +1.0 : -1.0), // K_L -> pi+ vs K_S -> pi+ and K- -> pi^0
                mu(p["ustaunutau::mu"], u),
                int_config(GSL::QAGS::Config().epsrel(0.5e-3))
            {
                Context ctx("When constructing tau^+ -> [K pi]^+ nubar observable");

                u.uses(*model);
                u.uses(*form_factors);
            }

            double
            differential_decay_width(const double & k2) const
            {
                // Return zero outside of physical phase space
                if (k2 < power_of<2>(m_K + m_pi) || k2 > power_of<2>(m_tau))
                {
                    return 0.0;
                }

                // Expression taken from [CCH:2017A], page 2, eqs. (12-14)
                // Except the missing factor of S_EW, as this is a RG evolution for the vector coefficient
                // alone, and is taken care of by our RG

                // Compare EOS basis with [CCH:2017A], page 2, eq. (9)
                const WilsonCoefficients<ChargedCurrent> wc  = model->wet_uslnu(LeptonFlavor::tauon, false);
                const complex<double>                    cV  = std::conj(wc.cvl() + wc.cvr());
                const complex<double>                    cA  = -std::conj(wc.cvl() + wc.cvr());
                const complex<double>                    cS  = std::conj(wc.csl() + wc.csr());
                const complex<double>                    icP = -std::conj(wc.csl() + wc.csr());
                const complex<double>                    cT  = 2.0 * std::conj(wc.ct());

                const auto m_tau2 = power_of<2>(m_tau);
                const auto m_K2   = power_of<2>(m_K);
                const auto m_pi2  = power_of<2>(m_pi);

                const auto fp = form_factors->f_p(k2);
                const auto f0 = form_factors->f_0(k2);
                // const auto BT = ( (-2.0 * m_K) / (m_K + m_pi) ) *  form_factors->f_t(k2);
                const auto BT = 0.0; // neglect tensor form factor as not implemented yet!

                const double lambda_piK = lambda(k2, m_pi2, m_K2);
                const auto   xi         = (m_tau2 + 2.0 * k2) * lambda_piK / (3.0 * m_tau2 * power_of<2>(m_K2 - m_pi2));

                const auto T = 3.0 * k2 * m_tau * cT * BT / ((m_tau2 + 2.0 * k2) * m_K);
                const auto V = fp * cV - T;
                const auto A = fp * cA + T;
                const auto S = f0 * (cV + k2 * cS / (m_tau * (m_s - m_u)));
                const auto P = f0 * (cA - k2 * icP / (m_tau * (m_s - m_u)));

                return power_of<2>(g_fermi * std::abs(model->ckm_us())) * sqrt(lambda_piK) * power_of<2>((m_tau2 - k2) * (m_K2 - m_pi2))
                       / (1024.0 * pow(M_PI, 3.0) * m_tau * pow(k2, 3.0))
                       * (xi * (std::norm(V) + std::norm(A) + 4.0 * power_of<2>(m_tau2 - k2) * std::norm(T) / (9.0 * k2 * m_tau2)) + std::norm(S) + std::norm(P));
            }

            double
            differential_branching_ratio(const double & k2) const
            {
                return differential_decay_width(k2) * tau_tau / hbar;
            }

            double
            total_branching_ratio() const
            {
                const double                          q2_min = power_of<2>(m_K() + m_pi());
                const double                          q2_max = power_of<2>(m_tau());
                std::function<double(const double &)> f      = std::bind(&Implementation<TauToKPiNeutrino>::differential_branching_ratio, this, std::placeholders::_1);
                return integrate<GSL::QAGS>(f, q2_min, q2_max, int_config);
            }

            double
            differential_pdf_q2(const double & k2) const
            {
                const double num   = differential_branching_ratio(k2);
                const double denom = total_branching_ratio();
                return num / denom;
            }

            double
            integrated_pdf_q2(const double & q2_min, const double & q2_max) const
            {
                std::function<double(const double &)> f     = std::bind(&Implementation<TauToKPiNeutrino>::differential_branching_ratio, this, std::placeholders::_1);
                const double                          num   = integrate<GSL::QAGS>(f, q2_min, q2_max, int_config);
                const double                          denom = total_branching_ratio();
                return num / denom / (q2_max - q2_min);
            }
    };

    const std::vector<OptionSpecification> Implementation<TauToKPiNeutrino>::options{
        Model::option_specification(),
        { "K"_ok, { "K_u", "K_S", "K_L" }, "K_u" },
    };

    TauToKPiNeutrino::TauToKPiNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<TauToKPiNeutrino>(new Implementation<TauToKPiNeutrino>(parameters, options, *this))
    {
    }

    TauToKPiNeutrino::~TauToKPiNeutrino() {}

    double
    TauToKPiNeutrino::differential_branching_ratio(const double & k2) const
    {
        return _imp->differential_branching_ratio(k2);
    }

    double
    TauToKPiNeutrino::differential_decay_width(const double & k2) const
    {
        return _imp->differential_decay_width(k2);
    }

    double
    TauToKPiNeutrino::decay_width(const double & q2_min, const double & q2_max) const
    {
        // Integrate the differential decay width over the specified kinematic range
        std::function<double(const double &)> f = std::bind(&Implementation<TauToKPiNeutrino>::differential_decay_width, _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, q2_min, q2_max, _imp->int_config);
    }

    double
    TauToKPiNeutrino::branching_ratio(const double & q2_min, const double & q2_max) const
    {
        return decay_width(q2_min, q2_max) * _imp->tau_tau / _imp->hbar;
    }

    double
    TauToKPiNeutrino::total_branching_ratio() const
    {
        return _imp->total_branching_ratio();
    }

    double
    TauToKPiNeutrino::differential_pdf_q2(const double & q2) const
    {
        return _imp->differential_pdf_q2(q2);
    }

    double
    TauToKPiNeutrino::integrated_pdf_q2(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_pdf_q2(q2_min, q2_max);
    }

    const std::set<ReferenceName> TauToKPiNeutrino::references{ "CCH:2017A"_rn };

    std::vector<OptionSpecification>::const_iterator
    TauToKPiNeutrino::begin_options()
    {
        return Implementation<TauToKPiNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    TauToKPiNeutrino::end_options()
    {
        return Implementation<TauToKPiNeutrino>::options.cend();
    }
} // namespace eos
