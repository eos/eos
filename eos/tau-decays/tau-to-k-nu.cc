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

#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/tau-decays/tau-to-k-nu.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    using std::norm;

    /*
     * Decay: tau^+ -> K^+ nubar, based on [S:2025A] and [DBG:2013A]
     */
    template <> struct Implementation<TauToKNeutrino>
    {
            SpecifiedOption opt_model;

            std::shared_ptr<Model> model;

            UsedParameter hbar;

            UsedParameter g_fermi;

            UsedParameter m_K;

            UsedParameter f_K;

            UsedParameter tau_tau;

            UsedParameter m_tau;

            BooleanOption opt_cp_conjugate;

            UsedParameter mu;

            static const std::vector<OptionSpecification> options;

            Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
                opt_model(o, options, "model"_ok),
                model(Model::make(opt_model.value(), p, o)),
                hbar(p["QM::hbar"], u),
                g_fermi(p["WET::G_Fermi"], u),
                m_K(p["mass::K_u"], u),
                f_K(p["decay-constant::K_u"], u),
                tau_tau(p["life_time::tau"], u),
                m_tau(p["mass::tau"], u),
                opt_cp_conjugate(o, options, "cp-conjugate"_ok),
                mu(p["ustaunutau::mu"], u)
            {
                Context ctx("When constructing tau->K-nu observable");

                u.uses(*model);
            }

            inline double
            beta_K() const
            {
                return (1.0 - power_of<2>(m_K() / m_tau()));
            }

            double
            decay_width() const
            {
                const WilsonCoefficients<ChargedCurrent> wc = model->wet_uslnu(LeptonFlavor::tauon, opt_cp_conjugate.value());

                // cf. [DBG2013], eq. (5), p. 5
                const complex<double> ga = wc.cvl() - wc.cvr();
                const complex<double> gp = wc.csl() - wc.csr();

                // masses
                const double m_tau = this->m_tau(), m_tau3 = power_of<3>(m_tau);
                const double m_K2   = power_of<2>(this->m_K());
                const double msatmu = model->m_s_msbar(mu);
                const double muatmu = model->m_u_msbar(mu);

                return power_of<2>(g_fermi * std::abs(model->ckm_us()) * f_K * beta_K()) * m_tau3 / (16.0 * M_PI) * norm(ga - gp * m_K2 / (msatmu + muatmu) / m_tau);
            }

            double
            branching_ratio() const
            {
                return decay_width() * tau_tau / hbar;
            }
    };

    const std::vector<OptionSpecification> Implementation<TauToKNeutrino>::options{
        Model::option_specification(),
        { "cp-conjugate"_ok, { "true", "false" }, "false" },
    };

    TauToKNeutrino::TauToKNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<TauToKNeutrino>(new Implementation<TauToKNeutrino>(parameters, options, *this))
    {
    }

    TauToKNeutrino::~TauToKNeutrino() {}

    double
    TauToKNeutrino::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    TauToKNeutrino::decay_width() const
    {
        return _imp->decay_width();
    }

    const std::set<ReferenceName> TauToKNeutrino::references{ "DBG:2013A"_rn, "S:2025A"_rn };

    std::vector<OptionSpecification>::const_iterator
    TauToKNeutrino::begin_options()
    {
        return Implementation<TauToKNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    TauToKNeutrino::end_options()
    {
        return Implementation<TauToKNeutrino>::options.cend();
    }
} // namespace eos
