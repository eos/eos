/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#ifndef EOS_GUARD_EOS_UTILS_DETECTOR_LEVEL_PDF_HH
#define EOS_GUARD_EOS_UTILS_DETECTOR_LEVEL_PDF_HH 1

#include <eos/maths/resolution-convolution.hh>
#include <eos/signal-pdf.hh>
#include <eos/utils/observable_cache.hh>

#include <memory>
#include <span>
#include <string>
#include <vector>

namespace eos
{
    /*
     * A DetectorLevelPDF represents a truth-level SignalPDF after convolution with a detector
     * resolution function. It is itself a SignalPDF, so it can be evaluated, plotted, and inspected
     * point-by-point like any other SignalPDF -- the same convolution that the unbinned likelihood
     * performs internally is exposed here for debugging.
     *
     * The truth PDF is sampled on a uniform tensor-product grid, convolved with the resolution via
     * eos::ResolutionConvolution, and the smeared grid is interpolated at the point given by
     * kinematics(). The convolution is circular, so the grid must be padded with a region in which
     * both the PDF and the resolution are negligible (see eos::ResolutionConvolution).
     *
     * The resolution can be supplied in two ways:
     *   - as a SignalPDF over the per-axis *offset* variables (the general, D-dimensional case), or
     *   - as a pre-computed grid of kernel values in centred (pre-ifftshift) order, ranks 1-4.
     *
     * As for every SignalPDF, evaluate() and evaluate_linear() return the *unnormalized* smeared
     * density and normalization() returns the logarithm of its normalization; a normalized log
     * density is evaluate() - normalization().
     *
     * The per-grid-point truth PDFs are registered with the ObservableCache handed in at
     * construction as a single batch, and the truth normalization as one further cached observable.
     * The PDF is tied to that cache for its lifetime and never creates one of its own, and it never
     * calls update() on it -- whoever owns the cache does. grid() and update_grid() separate reading
     * the smeared grid from recomputing it, so calling grid() inside a per-event loop is cheap.
     */
    class DetectorLevelPDF : public SignalPDF
    {
        public:
            // One sampling axis of the convolution grid. Along this axis the grid has `points`
            // (even, >= 2) points spanning the inclusive range [min, max]; the resolution PDF is
            // sampled over the corresponding displacement variable `offset_variable` (defaults to
            // `variable` when empty).
            struct Axis
            {
                    std::string variable;
                    double      min;
                    double      max;
                    std::size_t points;
                    std::string offset_variable;

                    Axis(const std::string & variable, double min, double max, std::size_t points, const std::string & offset_variable = "") :
                        variable(variable),
                        min(min),
                        max(max),
                        points(points),
                        offset_variable(offset_variable)
                    {
                    }
            };

            // General case: the resolution is a SignalPDF over the offset variables of `axes`.
            DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const QualifiedName & resolution_name, const Options & options,
                             const std::vector<Axis> & axes);

            // Sampled resolution, ranks 1-4: kernel values, flat and row-major, in centred order
            // (zero offset at index points / 2 per axis), sized to the product of the point counts.
            DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const std::vector<Axis> & axes,
                             const std::vector<double> & resolution);

            ~DetectorLevelPDF();

            ///@name Factories (return a SignalPDFPtr; used by the Python bindings)
            ///@{
            static SignalPDFPtr make(const ObservableCache & cache, const QualifiedName & signal_name, const QualifiedName & resolution_name, const Options & options,
                                     const std::vector<Axis> & axes);

            // Sampled resolution, ranks 1-4 (see the constructor above).
            static SignalPDFPtr make_from_grid(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const std::vector<Axis> & axes,
                                               const std::vector<double> & resolution);

            // 1D convenience: forwards to make_from_grid() with a single axis.
            static SignalPDFPtr make_1d(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const Axis & axis,
                                        const std::vector<double> & resolution);
            ///@}

            ///@name SignalPDF interface
            ///@{
            virtual const QualifiedName & name() const;

            virtual double evaluate() const;

            virtual double evaluate_linear() const;

            virtual double normalization() const;

            virtual ObservablePtr unnormalized_pdf() const;

            virtual ObservablePtr normalization_observable() const;

            virtual Kinematics kinematics();

            virtual Parameters parameters();

            virtual Options options();

            virtual DensityPtr clone() const;

            virtual DensityPtr clone(const Parameters & parameters) const;

            // Clone bound to an existing cache rather than to a fresh private one. The clone
            // registers its own batch with `cache`; blocks that must share the LogLikelihood's
            // cache use this.
            std::shared_ptr<DetectorLevelPDF> clone(const ObservableCache & cache) const;

            virtual Density::Iterator begin() const;

            virtual Density::Iterator end() const;
            ///@}

            ///@name Grid access
            ///@{
            // The smeared grid. Stable for the lifetime of this PDF: the buffer is allocated once
            // and never reallocated.
            std::span<const double> grid() const;

            // Recompute the smeared grid in place from the current cache contents. A no-op when the
            // cache's update generation has not advanced since the last call.
            void update_grid() const;

            // The cache this PDF is tied to.
            const ObservableCache & cache() const;

            // The sampling axes, in axis order.
            const std::vector<Axis> & axes() const;

            // Query points bound to grid(), with the cell index and fractional offsets precomputed
            // once (see ResolutionConvolution::make_interpolator()).
            std::unique_ptr<ResolutionConvolution::Interpolator> make_interpolator(const std::vector<std::vector<double>> & points) const;
            ///@}

        private:
            struct Data;

            std::unique_ptr<Data> _data;

            // Shared construction: build the engine, the query kinematics, and the per-grid-point
            // truth PDFs. `resolution_name` is empty in the pre-computed-grid case.
            DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const QualifiedName & resolution_name, const Options & options,
                             const std::vector<Axis> & axes, const std::vector<double> & resolution);
    };
} // namespace eos

#endif
