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
#include <eos/tau-decays/tau-to-k-pi-nu.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    using std::norm;

    /*
     * Decay: tau^- -> [K pi]^- nubar, cf. [CCH:2017A]
     */
    template <> struct Implementation<TauToKPiNeutrino>
    {
            SpecifiedOption opt_model;

            std::shared_ptr<Model> model;

            UsedParameter hbar;

            UsedParameter g_fermi;

            UsedParameter m_tau;

            QuarkFlavorOption opt_q;

            UsedParameter m_K;

            UsedParameter m_pi;

            UsedParameter tau_tau;

            UsedParameter mu;

            static const std::vector<OptionSpecification> options;

            Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
                opt_model(o, options, "model"),
                model(Model::make(opt_model.value(), p, o)),
                hbar(p["QM::hbar"], u),
                g_fermi(p["WET::G_Fermi"], u),
                m_tau(p["mass::tau"], u),
                opt_q(o, options, "q"),
                m_K(p["mass::K_" + opt_q.str()], u),
                m_pi(p["mass::pi^" + std::string(opt_q.value() == QuarkFlavor::up ? "0" : "-")], u),
                tau_tau(p["life_time::tau"], u),
                mu(p["ustaunutau::mu"], u)
            {
                Context ctx("When constructing tau^-->[K pi]^- nubar observable");

                u.uses(*model);
            }

            double
            differential_decay_width(const double & k2) const
            {
                const WilsonCoefficients<ChargedCurrent> wc = model->wet_uslnu(LeptonFlavor::tauon, false);

                return 0.0;
            }

            double
            differential_branching_ratio(const double & k2) const
            {
                return differential_decay_width(k2) * tau_tau / hbar;
            }
    };

    const std::vector<OptionSpecification> Implementation<TauToKPiNeutrino>::options{
        Model::option_specification(),
        { "q", { "u", "d" }, "u" }, // CHECK w.r.t. notes
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
