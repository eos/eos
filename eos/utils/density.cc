/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013, 2015, 2016 Danny van Dyk
 * Copyright (c) 2013 Frederik Beaujean
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

#include <eos/utils/density-impl.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <vector>

namespace eos
{
    Density::~Density()
    {
    }

    void
    Density::dump_descriptions(hdf5::File & file, const std::string & data_set_base) const
    {
        auto data_set = file.create_data_set(data_set_base + "/parameters", Density::Output::description_type());

        auto record = Density::Output::description_record();

        // loop over descriptions
        for (auto & d : *this)
        {
            std::get<0>(record) = d.parameter->name().c_str();
            std::get<1>(record) = d.min;
            std::get<2>(record) = d.max;
            std::get<3>(record) = int(d.nuisance);

            data_set << record;
        }
    }

    Density::Output::DescriptionType
    Density::Output::description_type()
    {
        return
        DescriptionType
        {
         "parameter description",
         hdf5::Scalar<const char *>("name"),
         hdf5::Scalar<double>("min"),
         hdf5::Scalar<double>("max"),
         hdf5::Scalar<int>("nuisance"),
        };
    }

    std::tuple<const char *, double, double, int>
    Density::Output::description_record()
    {
        return std::make_tuple("name", 1.0, 2.0, 3);
    }

    template class WrappedForwardIterator<Density::IteratorTag, const ParameterDescription>;

    template <>
    struct Implementation<ProductDensity>
    {
        DensityPtr x, y;

        std::vector<ParameterDescription> descriptions;

        Implementation(const DensityPtr & x, const DensityPtr & y) :
            x(x),
            y(y)
        {
            descriptions.insert(descriptions.end(), x->begin(), x->end());
            descriptions.insert(descriptions.end(), y->begin(), y->end());
        }
    };

    ProductDensity::ProductDensity(const DensityPtr & x, const DensityPtr & y) :
        PrivateImplementationPattern<ProductDensity>(new Implementation<ProductDensity>(x, y))
    {
    }

    ProductDensity::~ProductDensity()
    {
    }

    DensityPtr
    ProductDensity::clone() const
    {
        return DensityPtr(new ProductDensity(_imp->x->clone(), _imp->y->clone()));
    }

    double
    ProductDensity::evaluate() const
    {
        // since densities are evaluates in the log scale, the product turns
        // into a sum.
        return _imp->x->evaluate() + _imp->y->evaluate();
    }

    Density::Iterator
    ProductDensity::begin() const
    {
        return Density::Iterator(_imp->descriptions.cbegin());
    }

    Density::Iterator
    ProductDensity::end() const
    {
        return Density::Iterator(_imp->descriptions.cend());
    }
}
