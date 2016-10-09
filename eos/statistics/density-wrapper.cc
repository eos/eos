/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#include <eos/statistics/density-wrapper.hh>
#include <eos/utils/density-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

namespace eos
{
    template <>
    struct WrappedForwardIteratorTraits<SimpleParameters::IteratorTag>
    {
        typedef std::vector<ParameterDescription>::const_iterator UnderlyingIterator;
    };

    DensityWrapper::DensityWrapper(const WrappedDensity & density) :
        _density(density)
    {
    }

    DensityWrapper::DensityWrapper(RawDensity func) :
        _density(WrappedDensity(func))
    {
    }

    DensityWrapper::~DensityWrapper()
    {
    }

    void
    DensityWrapper::add_parameter(const std::string & name, const double & min,
            const double & max, bool nuisance)
    {
        _parameters.declare(name, min, max, nuisance);
    }

    SimpleParameters &
    DensityWrapper::parameters()
    {
        return _parameters;
    }

    double
    DensityWrapper::evaluate() const
    {
        return _density(_parameters.values());
    }

    DensityPtr
    DensityWrapper::clone() const
    {
        DensityWrapper * density = new DensityWrapper(_density);
        density->_parameters = this->_parameters.clone();
        return DensityPtr(density);
    }

    Density::Iterator
    DensityWrapper::begin() const
    {
        return Density::Iterator(_parameters.begin().underlying_iterator<std::vector<ParameterDescription>::const_iterator>());
    }

    Density::Iterator
    DensityWrapper::end() const
    {
        return Density::Iterator(_parameters.end().underlying_iterator<std::vector<ParameterDescription>::const_iterator>());
    }
}
