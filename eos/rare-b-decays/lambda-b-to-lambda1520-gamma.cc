/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Méril Reboud
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
#include <eos/rare-b-decays/lambda-b-to-lambda1520-gamma-base.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda1520-gamma-naive.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>

#include <gsl/gsl_sf.h>

namespace eos
{
    using namespace std::placeholders;
    using std::norm;

    /*!
     * Implementation for the decay @f$\bar{\Lambda_b} \to \bar{\Lambda}(1520) \gamma@f$.
     */
    template <>
    struct Implementation<LambdaBToLambda1520Gamma>
    {
        std::shared_ptr<LambdaBToLambda1520Gamma::AmplitudeGenerator> amplitude_generator;

        std::shared_ptr<Model> model;

        UsedParameter hbar;
        UsedParameter tau;
        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "WET"), p, o)),
            hbar(p["QM::hbar"], u),
            tau(p["life_time::Lambda_b"], u),
            mu(p["sb::mu"], u)
        {
            Context ctx("When constructing Lb->L(1520)gamma observables");

            std::string tag = o.get("tag", "");

            if ("Naive" == tag)
            {
                amplitude_generator.reset(new LambdaBToLambda1520GammaAmplitudes<tag::Naive>(p, o));
            }
            else
            {
                throw InternalError("LambdaBToLambda1520Gamma: Unknown tag or no valid tag specified (tag = '" + tag + "')!");
            }

            u.uses(*model);
            u.uses(*amplitude_generator);
        }

        double decay_rate()
        {
            auto amps = amplitude_generator->amplitudes();

            return norm(amps.a_perp12) + norm(amps.a_para12) + 3.0 * (norm(amps.a_perp32) + norm(amps.a_para32));
        }
    };

    LambdaBToLambda1520Gamma::LambdaBToLambda1520Gamma(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<LambdaBToLambda1520Gamma>(new Implementation<LambdaBToLambda1520Gamma>(parameters, options, *this))
    {
    }

    LambdaBToLambda1520Gamma::~LambdaBToLambda1520Gamma() = default;

    const std::vector<OptionSpecification>
    Implementation<LambdaBToLambda1520Gamma>::options
    {
    };

    double
    LambdaBToLambda1520Gamma::decay_rate() const
    {
        return _imp->decay_rate();
    }

    double
    LambdaBToLambda1520Gamma::branching_ratio() const
    {
        return _imp->decay_rate() * _imp->tau() / _imp->hbar();
    }

    const std::set<ReferenceName>
    LambdaBToLambda1520Gamma::references
    {
        "ABR:2022A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambda1520Gamma::begin_options()
    {
        return Implementation<LambdaBToLambda1520Gamma>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambda1520Gamma::end_options()
    {
        return Implementation<LambdaBToLambda1520Gamma>::options.cend();
    }
}
