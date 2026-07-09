/*
 * Copyright (c) 2026 Danny van Dyk
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

#include <eos/maths/resolution-convolution-impl.hh>

#include <format>

namespace eos
{
    template class ConcreteResolutionConvolution<1>;
    template class ConcreteResolutionConvolution<2>;
    template class ConcreteResolutionConvolution<3>;
    template class ConcreteResolutionConvolution<4>;

    ResolutionConvolution::ResolutionConvolution(const std::vector<AxisGeometry> & axes) :
        _axes(axes),
        _size(1)
    {
        for (const auto & axis : axes)
            _size *= axis.points;
    }

    ResolutionConvolution::~ResolutionConvolution() = default;

    std::unique_ptr<ResolutionConvolution>
    ResolutionConvolution::make(const std::vector<AxisGeometry> & axes)
    {
        switch (axes.size())
        {
            case 1:
                return std::make_unique<ConcreteResolutionConvolution<1>>(axes);
            case 2:
                return std::make_unique<ConcreteResolutionConvolution<2>>(axes);
            case 3:
                return std::make_unique<ConcreteResolutionConvolution<3>>(axes);
            case 4:
                return std::make_unique<ConcreteResolutionConvolution<4>>(axes);
            default:
                throw InternalError(std::format("ResolutionConvolution: dimensionality D = {} is not supported (must be in [1, 4])", axes.size()));
        }
    }
} // namespace eos
