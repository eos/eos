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

#ifndef EOS_GUARD_EOS_MATHS_RESOLUTION_CONVOLUTION_IMPL_HH
#define EOS_GUARD_EOS_MATHS_RESOLUTION_CONVOLUTION_IMPL_HH 1

#include <eos/maths/dft-container-impl.hh>
#include <eos/maths/dft-plan-impl.hh>
#include <eos/maths/resolution-convolution.hh>
#include <eos/utils/exception.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <numeric>

namespace eos
{
    template <std::size_t rank_>
    class ConcreteResolutionConvolution :
        public ResolutionConvolution
    {
        private:
            // Row-major grid geometry.
            std::array<std::size_t, rank_> _dimensions;
            std::array<std::size_t, rank_> _strides;

            // The DFT plans and the pre-computed spectrum of the resolution kernel. The forward plan
            // transforms the signal; the backward plan turns the product spectrum back into the
            // smeared grid.
            dft::Plan<rank_, dft::Direction::Forward>   _forward_plan;
            dft::Plan<rank_, dft::Direction::Backward>  _backward_plan;
            dft::Container<std::complex<double>, rank_> _resolution_dft;

            bool _have_resolution;

            // The convolved-grid buffer returned by convolve().
            std::vector<double> _result;

            // Validate the axes and extract the per-axis point counts.
            static std::array<std::size_t, rank_> to_dimensions(const std::vector<AxisGeometry> & axes)
            {
                if (axes.size() != rank_)
                    throw InternalError(std::format("ConcreteResolutionConvolution<{}>: expected {} axes but got {}", rank_, rank_, axes.size()));

                std::array<std::size_t, rank_> result;
                for (std::size_t d = 0 ; d < rank_ ; ++d)
                {
                    // The real DFT requires even dimensions; interpolation requires at least two points per axis.
                    if ((axes[d].points < 2) || (axes[d].points % 2 != 0))
                        throw InternalError(std::format("ConcreteResolutionConvolution<{}>: axis {} has {} points (must be even and >= 2)", rank_, d, axes[d].points));

                    if (! (axes[d].spacing > 0.0))
                        throw InternalError(std::format("ConcreteResolutionConvolution<{}>: axis {} has non-positive spacing {}", rank_, d, axes[d].spacing));

                    result[d] = axes[d].points;
                }

                return result;
            }

            // Row-major strides derived from the dimensions.
            static std::array<std::size_t, rank_> to_strides(const std::array<std::size_t, rank_> & dimensions)
            {
                std::array<std::size_t, rank_> result;
                result[rank_ - 1] = 1;
                for (std::size_t d = rank_ - 1 ; d-- > 0 ; )
                    result[d] = result[d + 1] * dimensions[d + 1];

                return result;
            }

            // Decompose a flat row-major index into its per-axis multi-index.
            std::array<std::size_t, rank_> multi_index(std::size_t flat) const
            {
                std::array<std::size_t, rank_> result;
                for (std::size_t d = 0 ; d < rank_ ; ++d)
                    result[d] = (flat / _strides[d]) % _dimensions[d];

                return result;
            }

        public:
            explicit ConcreteResolutionConvolution(const std::vector<AxisGeometry> & axes) :
                ResolutionConvolution(axes),
                _dimensions(to_dimensions(axes)),
                _strides(to_strides(_dimensions)),
                _forward_plan(_dimensions),
                _backward_plan(_dimensions),
                _resolution_dft(dft::impl::frequency_dimensions<rank_>(_dimensions)),
                _have_resolution(false),
                _result(_size, 0.0)
            {
            }

            virtual ~ConcreteResolutionConvolution() = default;

            virtual void set_resolution(const std::vector<double> & kernel_centred) override
            {
                if (kernel_centred.size() != _size)
                    throw InternalError(std::format("ResolutionConvolution: resolution kernel has {} entries but the grid has {} points", kernel_centred.size(), _size));

                const double total = std::accumulate(kernel_centred.begin(), kernel_centred.end(), 0.0);
                if (! (total > 0.0))
                    throw InternalError("ResolutionConvolution: resolution kernel must have positive total weight");

                const double inv_total = 1.0 / total;

                // Write the kernel into the forward plan's time-domain buffer, simultaneously
                // renormalising it to unit sum and converting from centred to wrap-around order.
                // For the (necessarily even) dimensions used here, the shift along each axis is by
                // points / 2, so the centred sample at index k maps to the wrap-around index
                // (k + points / 2) mod points -- this is a per-axis np.fft.ifftshift.
                double * time_data = _forward_plan.time_domain_container().data();
                for (std::size_t flat = 0 ; flat < _size ; ++flat)
                {
                    const std::array<std::size_t, rank_> centred = multi_index(flat);

                    std::size_t shifted = 0;
                    for (std::size_t d = 0 ; d < rank_ ; ++d)
                        shifted += ((centred[d] + _dimensions[d] / 2) % _dimensions[d]) * _strides[d];

                    time_data[shifted] = kernel_centred[flat] * inv_total;
                }

                _forward_plan.transform();
                _resolution_dft = _forward_plan.frequency_domain_container();

                _have_resolution = true;
            }

            virtual const std::vector<double> & convolve(const std::vector<double> & signal_grid) override
            {
                if (! _have_resolution)
                    throw InternalError("ResolutionConvolution::convolve: called before set_resolution");

                if (signal_grid.size() != _size)
                    throw InternalError(std::format("ResolutionConvolution::convolve: signal has {} entries but the grid has {} points", signal_grid.size(), _size));

                // Forward DFT of the signal values.
                double * time_data = _forward_plan.time_domain_container().data();
                std::copy(signal_grid.begin(), signal_grid.end(), time_data);
                _forward_plan.transform();

                // Multiply the signal spectrum by the pre-computed resolution spectrum, leaving the
                // result in the backward plan's frequency-domain buffer. The data must be copied in
                // place: assigning the container (operator=) would reallocate its buffer and
                // invalidate the pointer the backward FFTW plan was created with.
                const auto  freq_dims = dft::impl::frequency_dimensions<rank_>(_dimensions);
                std::size_t freq_size = 1;
                for (const auto & d : freq_dims)
                    freq_size *= d;

                const std::complex<double> * forward_freq_data  = _forward_plan.frequency_domain_container().data();
                std::complex<double> *       backward_freq_data = _backward_plan.frequency_domain_container().data();
                std::copy(forward_freq_data, forward_freq_data + freq_size, backward_freq_data);
                _backward_plan.frequency_domain_container() *= _resolution_dft;

                // Backward DFT. The DFT is unnormalised, so the result is N times the circular
                // convolution; we divide by N here.
                _backward_plan.transform();

                const double   norm   = 1.0 / static_cast<double>(_size);
                const double * result = _backward_plan.time_domain_container().data();
                for (std::size_t i = 0 ; i < _size ; ++i)
                    _result[i] = result[i] * norm;

                return _result;
            }

            virtual double interpolate(const std::vector<double> & grid, const std::vector<double> & point) const override
            {
                if (grid.size() != _size)
                    throw InternalError(std::format("ResolutionConvolution::interpolate: grid has {} entries but the geometry has {} points", grid.size(), _size));

                if (point.size() != rank_)
                    throw InternalError(std::format("ResolutionConvolution::interpolate: point has {} coordinates but the grid has rank {}", point.size(), rank_));

                std::size_t               base = 0;
                std::array<double, rank_> fraction;
                for (std::size_t d = 0 ; d < rank_ ; ++d)
                {
                    const double p = (point[d] - _axes[d].origin) / _axes[d].spacing;

                    if ((p < 0.0) || (p > static_cast<double>(_dimensions[d] - 1)))
                        throw InternalError(std::format("ResolutionConvolution::interpolate: point lies outside the grid along axis {}", d));

                    // Lower corner index, clamped so that the upper corner (index + 1) is always valid.
                    std::size_t j = static_cast<std::size_t>(std::floor(p));
                    double      t = p - static_cast<double>(j);
                    if (j >= _dimensions[d] - 1)
                    {
                        j = _dimensions[d] - 2;
                        t = 1.0;
                    }

                    base       += j * _strides[d];
                    fraction[d] = t;
                }

                double value = 0.0;
                for (std::size_t corner = 0 ; corner < (std::size_t(1) << rank_) ; ++corner)
                {
                    double      weight = 1.0;
                    std::size_t offset = 0;
                    for (std::size_t d = 0 ; d < rank_ ; ++d)
                    {
                        const bool upper = (corner >> d) & 1u;
                        weight *= upper ? fraction[d] : (1.0 - fraction[d]);
                        offset += upper ? _strides[d] : 0;
                    }
                    value += weight * grid[base + offset];
                }

                return value;
            }
    };
} // namespace eos

#endif
