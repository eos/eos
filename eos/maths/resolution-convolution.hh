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

#ifndef EOS_GUARD_EOS_MATHS_RESOLUTION_CONVOLUTION_HH
#define EOS_GUARD_EOS_MATHS_RESOLUTION_CONVOLUTION_HH 1

#include <cstddef>
#include <memory>
#include <span>
#include <vector>

namespace eos
{
    /*!
     * Circular convolution of a signal sampled on a uniform tensor-product grid with a fixed,
     * non-negative resolution kernel, evaluated via the discrete Fourier transform.
     *
     * This is the numerical core shared by DetectorLevelPDF and (in the future) the unbinned
     * likelihood. It is pure numerics and knows nothing about SignalPDF: a caller samples its
     * signal and resolution onto the grid and hands the values in as flat, row-major arrays.
     *
     * The class is rank-erased: a runtime dimensionality D in [1, 4] is dispatched by make() to a
     * compile-time rank, which is the range the underlying dft::Plan is instantiated for.
     *
     * The convolution is circular. The caller is responsible for padding the grid with a
     * sufficiently large region in which both the signal and the resolution are negligible, so
     * that no appreciable density wraps across the boundary.
     */
    class ResolutionConvolution
    {
        public:
            // Geometry of one axis: grid index i has coordinate origin + i * spacing.
            struct AxisGeometry
            {
                    double      origin;  // coordinate of grid index 0
                    double      spacing; // uniform step (> 0)
                    std::size_t points;  // number of points (even and >= 2)
            };

            /*!
             * A set of query points bound to one grid at construction: the flat index of the lower
             * corner and the per-axis fractional offsets of every point are precomputed once, so
             * operator() takes only a point index and cannot be aimed at a grid of a different
             * geometry.
             */
            class Interpolator
            {
                public:
                    virtual ~Interpolator();

                    /// Number of query points.
                    virtual std::size_t size() const = 0;

                    /// Multilinear interpolation of the bound grid at query point i.
                    virtual double operator() (std::size_t i) const = 0;
            };

        protected:
            std::vector<AxisGeometry> _axes;
            std::size_t               _size; // total number of grid points, i.e. the product of the axis point counts

            explicit ResolutionConvolution(const std::vector<AxisGeometry> & axes);

        public:
            virtual ~ResolutionConvolution();

            /*!
             * Install (or refresh) the resolution kernel, sampled on the *centred* offset grid in
             * natural row-major order -- i.e. the zero offset (the peak of a symmetric kernel) sits
             * at the central index (points / 2) of each axis. The engine renormalises the kernel to
             * unit sum (so that the convolution preserves the signal's normalisation), shifts it to
             * the wrap-around order the circular convolution expects, and caches its forward DFT.
             * Call only when the resolution actually changes.
             *
             * @param kernel_centred The resolution kernel, in centred order, one entry per grid point.
             *
             * Throws if the kernel size does not match the grid, or if its total weight is not
             * positive.
             */
            virtual void set_resolution(const std::vector<double> & kernel_centred) = 0;

            /*!
             * Convolve a signal sampled on the grid (flat, row-major) with the cached resolution.
             *
             * Returns a reference to the internal result buffer, allocated once and thereafter only
             * overwritten elementwise, never resized or reassigned, so the reference (and any span
             * taken over it, see result()) stays valid for the lifetime of the engine.
             *
             * @param signal_grid The signal, sampled on the grid, flat and row-major.
             *
             * Throws if set_resolution() has not been called, or if the signal size does not match
             * the grid.
             */
            virtual const std::vector<double> & convolve(const std::vector<double> & signal_grid) = 0;

            /*!
             * The buffer written by convolve(), without triggering a new convolution. Stable for the
             * lifetime of the engine, per the contract on convolve() above.
             */
            virtual std::span<const double> result() const = 0;

            /*!
             * Multilinear interpolation of a grid (as returned by convolve() or result()) at a
             * physical point.
             *
             * @param grid  The grid to interpolate, flat and row-major.
             * @param point One coordinate per axis, in axis order.
             *
             * Throws if the point lies outside the grid, or if the sizes do not match.
             */
            virtual double interpolate(std::span<const double> grid, const std::vector<double> & point) const = 0;

            /*!
             * Build an Interpolator bound to grid and a fixed set of query points.
             *
             * @param grid   The grid to interpolate, as returned by convolve() or result().
             * @param points The query points, one coordinate vector per point.
             *
             * Throws under the same conditions and message style as interpolate(): a grid size
             * mismatch, a point with the wrong number of coordinates, or a point outside the grid.
             */
            virtual std::unique_ptr<Interpolator> make_interpolator(std::span<const double> grid, const std::vector<std::vector<double>> & points) const = 0;

            // Dimensionality D of the grid.
            std::size_t
            rank() const
            {
                return _axes.size();
            }

            // Total number of grid points N = product of the axis point counts.
            std::size_t
            size() const
            {
                return _size;
            }

            const std::vector<AxisGeometry> &
            geometry() const
            {
                return _axes;
            }

            /*!
             * Construct an engine for the given axes.
             *
             * @param axes The per-axis grid geometry; the runtime dimensionality D = axes.size() is
             *             dispatched to a compile-time rank.
             *
             * Throws for D outside [1, 4].
             */
            static std::unique_ptr<ResolutionConvolution> make(const std::vector<AxisGeometry> & axes);
    };
} // namespace eos

#endif
