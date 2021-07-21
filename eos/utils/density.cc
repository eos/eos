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
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <vector>

namespace eos
{
    Density::~Density()
    {
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
