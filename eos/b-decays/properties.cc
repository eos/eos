/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012-2024 Danny van Dyk
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

#include <eos/b-decays/properties.hh>
#include <eos/models/model.hh>
#include <eos/utils/options.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <algorithm>
#include <array>

namespace eos
{
    template <>
    struct Implementation<BMesonProperties>
    {
        std::shared_ptr<Model> model;

        QuarkFlavorOption opt_q;

        // HQE non-perturbative matrix elements
        UsedParameter mu2_g;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            opt_q(o, options, "q"),
            mu2_g(p["B->B::mu_G^2@1GeV"], u)
        {
            Context ctx("When constructing the B meson properties");

            u.uses(*model);
        }

        double mass_splitting_j1_j0() const
        {
            // We use a kinetic scale of 1 GeV
            static const double mu_kin = 1.0;

            // The NLO contribution is given by [U2001], eq. (9), p. 5
            double c = (1.0 + 3.0 * model->alpha_s(4.6) / (2.0 * M_PI) * (2.0 + std::log(mu_kin / 4.6))) / model->m_b_kin(mu_kin);

            // Corrections of order 1/mb^2 can be estimates, cf. [U2001], eqs. (12) and (18)
            double sum_rho3 = -0.5; // Assumption for the sum of 1/mb^2 operators from eq. (18)
            double delta = sum_rho3 / (3.0 * power_of<2>(model->m_b_kin(mu_kin)));

            // cf. [N1997], p. 15, eq. (32), we also consider alpha_s corrections
            // to the rate, cf. [U2001], eq. (11), p. 5.
            return 2.0 / 3.0 * c * mu2_g + delta;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BMesonProperties>::options
    {
        Model::option_specification(),
        { "q", { "u", "d", "s", "c"}, "d" }
    };

    BMesonProperties::BMesonProperties(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BMesonProperties>(new Implementation<BMesonProperties>(parameters, options, *this))
    {
    }

    BMesonProperties::~BMesonProperties()
    {
    }

    double
    BMesonProperties::mass_splitting_j1_j0() const
    {
        return _imp->mass_splitting_j1_j0();
    }

    const std::set<ReferenceName>
    BMesonProperties::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BMesonProperties::begin_options()
    {
        return Implementation<BMesonProperties>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BMesonProperties::end_options()
    {
        return Implementation<BMesonProperties>::options.cend();
    }
}
