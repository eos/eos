/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Danny van Dyk
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
#include <eos/s-decays/k-to-l-nu.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    using std::norm;

    /*
     * Decay: K^- -> l^- nubar, adapter from cf. [DBG:2013A]
     */
    template <> struct Implementation<KToLeptonNeutrino>
    {
            SpecifiedOption opt_model;

            std::shared_ptr<Model> model;

            UsedParameter hbar;

            UsedParameter g_fermi;

            UsedParameter m_K;

            UsedParameter f_K;

            UsedParameter tau_K;

            LeptonFlavorOption opt_l;

            UsedParameter m_l;

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
                tau_K(p["life_time::K_u"], u),
                opt_l(o, options, "l"_ok),
                m_l(p["mass::" + opt_l.str()], u),
                opt_cp_conjugate(o, options, "cp-conjugate"_ok),
                mu(p["us" + opt_l.str() + "nu" + opt_l.str() + "::mu"], u)
            {
                Context ctx("When constructing K_u->lnu observable");

                u.uses(*model);
            }

            inline double
            beta_l() const
            {
                return (1.0 - power_of<2>(m_l() / m_K()));
            }

            double
            decay_width() const
            {
                const WilsonCoefficients<ChargedCurrent> wc = model->wet_uslnu(opt_l.value(), opt_cp_conjugate.value());

                // cf. [DBG2013], eq. (5), p. 5
                const complex<double> ga = wc.cvl() - wc.cvr();
                const complex<double> gp = wc.csl() - wc.csr();

                // masses
                const double m_K = this->m_K(), m_K2 = m_K * m_K;
                const double m_l    = this->m_l();
                const double msatmu = model->m_s_msbar(mu);
                const double muatmu = model->m_u_msbar(mu);

                return power_of<2>(g_fermi * std::abs(model->ckm_us()) * f_K * beta_l()) * m_K / (8.0 * M_PI) * norm(ga * m_l - gp * m_K2 / (msatmu + muatmu));
            }

            double
            branching_ratio() const
            {
                return decay_width() * tau_K / hbar;
            }
    };

    const std::vector<OptionSpecification> Implementation<KToLeptonNeutrino>::options{
        Model::option_specification(),
        { "cp-conjugate"_ok, { "true", "false" }, "false" },
        {            "l"_ok,       { "e", "mu" },    "mu" },
    };

    KToLeptonNeutrino::KToLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<KToLeptonNeutrino>(new Implementation<KToLeptonNeutrino>(parameters, options, *this))
    {
    }

    KToLeptonNeutrino::~KToLeptonNeutrino() {}

    double
    KToLeptonNeutrino::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    KToLeptonNeutrino::decay_width() const
    {
        return _imp->decay_width();
    }

    const std::set<ReferenceName> KToLeptonNeutrino::references{ "DBG:2013A"_rn };

    std::vector<OptionSpecification>::const_iterator
    KToLeptonNeutrino::begin_options()
    {
        return Implementation<KToLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    KToLeptonNeutrino::end_options()
    {
        return Implementation<KToLeptonNeutrino>::options.cend();
    }
} // namespace eos
