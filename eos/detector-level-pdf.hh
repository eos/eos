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

#ifndef EOS_GUARD_EOS_DETECTOR_LEVEL_PDF_HH
#define EOS_GUARD_EOS_DETECTOR_LEVEL_PDF_HH 1

#include <eos/maths/resolution-convolution.hh>
#include <eos/signal-pdf.hh>
#include <eos/utils/observable_cache.hh>

#include <memory>
#include <string>
#include <vector>

namespace eos
{
    /*
     * A DetectorLevelPDF represents a truth-level SignalPDF after convolution with a detector
     * resolution function. It is itself a SignalPDF, so it can be evaluated, plotted, and inspected
     * point-by-point like any other SignalPDF -- the convolution that the unbinned likelihood
     * performs internally is materialised here for debugging.
     *
     * The truth PDF is sampled on a uniform tensor-product grid, convolved with the resolution via
     * eos::ResolutionConvolution, and the smeared grid is interpolated at the point given by
     * kinematics(). The convolution is circular, so the grid must be padded with a region in which
     * both the PDF and the resolution are negligible (see eos::ResolutionConvolution).
     *
     * The resolution can be supplied in two ways:
     *   - as a SignalPDF over the per-axis *offset* variables (the general, D-dimensional case), or
     *   - as a pre-computed grid of kernel values in centred (pre-ifftshift) order (1D only).
     *
     * An ObservableCache is required at construction. At present the cache is only the common source
     * of Parameters and the grid of truth-PDF values is computed serially; requiring the cache here
     * prepares the interface for a future switch to batched, cache-driven grid evaluation without a
     * change to the external API.
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

            // 1D convenience: the resolution is a pre-computed grid of `axis.points` kernel values in
            // centred (pre-ifftshift) order, i.e. the zero offset sits at index points / 2.
            DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const Axis & axis, const std::vector<double> & resolution);

            ~DetectorLevelPDF();

            ///@name Factories (return a SignalPDFPtr; used by the Python bindings)
            ///@{
            static SignalPDFPtr make(const ObservableCache & cache, const QualifiedName & signal_name, const QualifiedName & resolution_name, const Options & options,
                                     const std::vector<Axis> & axes);

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

            virtual Kinematics kinematics();

            virtual Parameters parameters();

            virtual Options options();

            virtual DensityPtr clone() const;

            virtual DensityPtr clone(const Parameters & parameters) const;

            virtual Density::Iterator begin() const;

            virtual Density::Iterator end() const;
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
