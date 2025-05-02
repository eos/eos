/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2025 Danny van Dyk
 * Copyright (c) 2018      Ahmet Kokulu
 * Copyright (c) 2018      Christoph Bobeth
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

#include <eos/b-decays/inclusive-b-to-u.hh>
#include <eos/maths/complex.hh>
#include <eos/models/model.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    template <>
    struct Implementation<BToXuLeptonNeutrino<Naive>>
    {
        std::shared_ptr<Model> model;

        LeptonFlavorOption opt_l;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            opt_l(o, options, "l"_ok)
        {
            Context ctx("When constructing B->X_ulnu observable");

            u.uses(*model);
        }

        ~Implementation()
        {
        }

        double v_ub() const
        {
            // inclusive |V_ub|^2 = |V_ub^eff|^2 (|C_V,LL|^2 + |C_V,RL|^2)

            double v_ub_eff_squared = std::norm(model->ckm_ub());
            auto wc = model->wet_ublnu(opt_l.value(), false);

            return std::sqrt(v_ub_eff_squared * (std::norm(wc.cvl()) + std::norm(wc.cvr())));
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToXuLeptonNeutrino<Naive>>::options
    {
        Model::option_specification(),
        { "l"_ok, { "e", "mu", "tau" }, "mu" }
    };

    BToXuLeptonNeutrino<Naive>::BToXuLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToXuLeptonNeutrino<Naive>>(new Implementation<BToXuLeptonNeutrino<Naive>>(parameters, options, *this))
    {
    }

    BToXuLeptonNeutrino<Naive>::~BToXuLeptonNeutrino()
    {
    }

    double
    BToXuLeptonNeutrino<Naive>::v_ub() const
    {
        return _imp->v_ub();
    }

    const std::set<ReferenceName>
    BToXuLeptonNeutrino<Naive>::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BToXuLeptonNeutrino<Naive>::begin_options()
    {
        return Implementation<BToXuLeptonNeutrino<Naive>>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToXuLeptonNeutrino<Naive>::end_options()
    {
        return Implementation<BToXuLeptonNeutrino<Naive>>::options.cend();
    }
}
