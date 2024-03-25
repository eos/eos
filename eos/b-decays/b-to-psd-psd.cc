/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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

#include <eos/b-decays/b-to-psd-psd.hh>
#include <eos/models/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{

    template <>
    struct Implementation<BToPseudoscalarPseudoscalar>
    {
        QuarkFlavorOption opt_q;
        LightMesonOption opt_p1;
        LightMesonOption opt_p2;
        UsedParameter hbar;
        UsedParameter tau;
        std::shared_ptr<NonleptonicAmplitudes<PToPP>> nl_amplitudes;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_q(o, options, "q"),
            opt_p1(o, options, "P1"),
            opt_p2(o, options, "P2"),
            hbar(p["QM::hbar"], u),
            tau(p["life_time::B_" + opt_q.str()], u),
            nl_amplitudes(NonleptonicAmplitudeFactory<PToPP>::create("B->PP::" + o.get("representation", "topological"), p, o))
        {
            Context ctx("When constructing B->PP observable");

            u.uses(*nl_amplitudes);
        }

        double decay_width() const
        {
            return std::norm(nl_amplitudes->amplitude());
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToPseudoscalarPseudoscalar>::options
    {
        Model::option_specification(),
        NonleptonicAmplitudeFactory<PToPP>::option_specification(),
        { "q", { "u", "d", "s" }, "d" },
        { "P1", { "pi^0", "pi^+", "K_d", "K_u", "eta_q", "eta'_q", "eta_s", "eta'_s" }, "pi^0" },
        { "P2", { "pi^0", "pi^+", "K_d", "K_u", "eta_q", "eta'_q", "eta_s", "eta'_s" }, "pi^0" },
    };

    BToPseudoscalarPseudoscalar::BToPseudoscalarPseudoscalar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToPseudoscalarPseudoscalar>(new Implementation<BToPseudoscalarPseudoscalar>(parameters, options, *this))
    {

    }

    BToPseudoscalarPseudoscalar::~BToPseudoscalarPseudoscalar()
    {
    }

    double
    BToPseudoscalarPseudoscalar::decay_width() const
    {
        return _imp->decay_width();
    }

    double
    BToPseudoscalarPseudoscalar::branching_ratio() const
    {
        return _imp->decay_width() * _imp->tau() / _imp->hbar();
    }

    const std::string
    BToPseudoscalarPseudoscalar::description = "\
    The decay B->PP, where all states are pseudoscalars.";

    const std::set<ReferenceName>
    BToPseudoscalarPseudoscalar::references
    {
        "HTX:2021A"_rn,
    };

    std::vector<OptionSpecification>::const_iterator
    BToPseudoscalarPseudoscalar::begin_options()
    {
        return Implementation<BToPseudoscalarPseudoscalar>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToPseudoscalarPseudoscalar::end_options()
    {
        return Implementation<BToPseudoscalarPseudoscalar>::options.cend();
    }
}
