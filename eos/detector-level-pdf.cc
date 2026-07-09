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

#include <eos/detector-level-pdf.hh>
#include <eos/utils/density-impl.hh>

#include <cmath>
#include <format>
#include <limits>

namespace eos
{
    struct DetectorLevelPDF::Data
    {
            // Construction inputs, retained so that clone() can rebuild an independent copy.
            ObservableCache   cache;
            QualifiedName     name;
            QualifiedName     signal_name;
            QualifiedName     resolution_name;
            Options           options;
            std::vector<Axis> axes;
            bool              has_resolution_pdf;

            Parameters parameters;

            // Grid geometry (row-major), mirroring eos::ResolutionConvolution.
            std::vector<std::size_t> dimensions;
            std::vector<std::size_t> strides;
            std::vector<double>      origin;
            std::vector<double>      spacing;
            std::vector<std::string> offset_variables;

            // The convolution engine, the truth PDF at each grid point, and (in the general case) the
            // resolution PDF swept over the offset grid.
            std::unique_ptr<ResolutionConvolution> engine;
            std::vector<SignalPDFPtr>              signal_pdfs;
            SignalPDFPtr                           resolution_pdf;

            // The query point exposed via kinematics(); interpolation reads its sampling variables.
            Kinematics kinematics;

            // Scratch buffers: the truth-PDF grid values and the (re-)sampled resolution kernel.
            std::vector<double> signal_values;
            std::vector<double> resolution_grid;

            // Empty; present so that Density::begin()/end() have something to iterate over.
            std::vector<ParameterDescription> descriptions;

            Data(const ObservableCache & cache, const QualifiedName & signal_name, const QualifiedName & resolution_name) :
                cache(cache),
                name(signal_name),
                signal_name(signal_name),
                resolution_name(resolution_name),
                has_resolution_pdf(false),
                parameters(cache.parameters())
            {
            }
    };

    DetectorLevelPDF::DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const QualifiedName & resolution_name, const Options & options,
                                       const std::vector<Axis> & axes, const std::vector<double> & resolution) :
        _data(new Data(cache, signal_name, resolution_name))
    {
        Data & d = *_data;

        if (axes.empty())
        {
            throw InternalError("DetectorLevelPDF: at least one axis is required");
        }

        d.options            = options;
        d.axes               = axes;
        // The general case supplies the resolution as a SignalPDF and an empty grid; the 1D
        // convenience case supplies a pre-computed grid.
        d.has_resolution_pdf = resolution.empty();

        const std::size_t rank = axes.size();

        // Assemble the grid geometry and construct the convolution engine (which validates the
        // per-axis point counts, and the dimensionality D in [1, 4]).
        std::vector<ResolutionConvolution::AxisGeometry> geometry;
        geometry.reserve(rank);
        d.dimensions.resize(rank);
        d.origin.resize(rank);
        d.spacing.resize(rank);
        for (std::size_t i = 0; i < rank; ++i)
        {
            if (axes[i].max <= axes[i].min)
            {
                throw InternalError(std::format("DetectorLevelPDF: axis '{}' has empty range [{}, {}]", axes[i].variable, axes[i].min, axes[i].max));
            }

            const double spacing = (axes[i].max - axes[i].min) / static_cast<double>(axes[i].points - 1);

            d.dimensions[i] = axes[i].points;
            d.origin[i]     = axes[i].min;
            d.spacing[i]    = spacing;
            geometry.push_back(ResolutionConvolution::AxisGeometry{ axes[i].min, spacing, axes[i].points });
        }
        d.engine = ResolutionConvolution::make(geometry);

        // Row-major strides, mirroring eos::ResolutionConvolution.
        d.strides.resize(rank);
        d.strides[rank - 1] = 1;
        for (std::size_t i = rank - 1; i-- > 0;)
        {
            d.strides[i] = d.strides[i + 1] * d.dimensions[i + 1];
        }

        const std::size_t N = d.engine->size();

        // The query point: one entry per sampling variable, initialised to the grid centre. The
        // [min, max] bounds are declared alongside (named '<variable>_min/_max') so that the PDF
        // presents the same variable/bound structure as any other SignalPDF.
        for (std::size_t i = 0; i < rank; ++i)
        {
            d.kinematics.declare(axes[i].variable, 0.5 * (axes[i].min + axes[i].max));
            d.kinematics.declare(axes[i].variable + "_min", axes[i].min);
            d.kinematics.declare(axes[i].variable + "_max", axes[i].max);
        }

        // Build the truth PDF at each grid point. Its kinematics carry the sampling variable at the
        // grid coordinate together with the [min, max] bounds (named '<variable>_min/_max', the
        // convention used by the SignalPDF normalisation observable).
        d.signal_pdfs.reserve(N);
        for (std::size_t flat = 0; flat < N; ++flat)
        {
            Kinematics k;
            for (std::size_t i = 0; i < rank; ++i)
            {
                const std::size_t index = (flat / d.strides[i]) % d.dimensions[i];

                k.declare(axes[i].variable, d.origin[i] + static_cast<double>(index) * d.spacing[i]);
                k.declare(axes[i].variable + "_min", axes[i].min);
                k.declare(axes[i].variable + "_max", axes[i].max);
            }
            d.signal_pdfs.push_back(SignalPDF::make(signal_name, d.parameters, k, options));
        }
        d.signal_values.assign(N, 0.0);

        if (d.has_resolution_pdf)
        {
            // Build the resolution PDF over the per-axis offset variables. Only its shape enters the
            // kernel (the engine renormalises it to unit sum), so its own normalisation is unused;
            // the offset bounds are declared solely so the SignalPDF constructs.
            d.offset_variables.resize(rank);
            Kinematics rk;
            for (std::size_t i = 0; i < rank; ++i)
            {
                const std::string offset_variable = axes[i].offset_variable.empty() ? axes[i].variable : axes[i].offset_variable;
                d.offset_variables[i]             = offset_variable;

                const double lo = -static_cast<double>(d.dimensions[i] / 2) * d.spacing[i];
                const double hi = static_cast<double>(d.dimensions[i] / 2 - 1) * d.spacing[i];

                rk.declare(offset_variable, 0.0);
                rk.declare(offset_variable + "_min", lo);
                rk.declare(offset_variable + "_max", hi);
            }
            d.resolution_pdf = SignalPDF::make(resolution_name, d.parameters, rk, options);
            d.resolution_grid.assign(N, 0.0);
        }
        else
        {
            // Pre-computed grid (centred, pre-ifftshift). Fixed for the lifetime of the PDF, so its
            // spectrum is cached once here.
            if (resolution.size() != N)
            {
                throw InternalError(std::format("DetectorLevelPDF: pre-computed resolution has {} entries but the grid has {} points", resolution.size(), N));
            }

            d.resolution_grid = resolution;
            d.engine->set_resolution(d.resolution_grid);
        }
    }

    DetectorLevelPDF::DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const QualifiedName & resolution_name, const Options & options,
                                       const std::vector<Axis> & axes) :
        DetectorLevelPDF(cache, signal_name, resolution_name, options, axes, std::vector<double>{})
    {
    }

    DetectorLevelPDF::DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const Axis & axis,
                                       const std::vector<double> & resolution) :
        DetectorLevelPDF(cache, signal_name, signal_name /* unused resolution name */, options, std::vector<Axis>{ axis }, resolution)
    {
        if (resolution.empty())
        {
            throw InternalError("DetectorLevelPDF: pre-computed resolution grid must not be empty");
        }
    }

    DetectorLevelPDF::~DetectorLevelPDF() = default;

    SignalPDFPtr
    DetectorLevelPDF::make(const ObservableCache & cache, const QualifiedName & signal_name, const QualifiedName & resolution_name, const Options & options,
                           const std::vector<Axis> & axes)
    {
        return SignalPDFPtr(new DetectorLevelPDF(cache, signal_name, resolution_name, options, axes));
    }

    SignalPDFPtr
    DetectorLevelPDF::make_1d(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const Axis & axis, const std::vector<double> & resolution)
    {
        return SignalPDFPtr(new DetectorLevelPDF(cache, signal_name, options, axis, resolution));
    }

    const QualifiedName &
    DetectorLevelPDF::name() const
    {
        return _data->name;
    }

    double
    DetectorLevelPDF::evaluate_linear() const
    {
        Data & d = *_data;

        const std::size_t rank = d.axes.size();
        const std::size_t N    = d.engine->size();

        // Sample the truth PDF on the grid. The grid is computed serially, one point at a time.
        for (std::size_t flat = 0; flat < N; ++flat)
        {
            d.signal_values[flat] = d.signal_pdfs[flat]->evaluate_linear();
        }

        // In the general case the resolution may depend on the parameters, so its kernel is
        // re-sampled from the resolution PDF (on the centred offset grid) on every evaluation.
        if (d.has_resolution_pdf)
        {
            Kinematics rk = d.resolution_pdf->kinematics();
            for (std::size_t flat = 0; flat < N; ++flat)
            {
                for (std::size_t i = 0; i < rank; ++i)
                {
                    const std::size_t index  = (flat / d.strides[i]) % d.dimensions[i];
                    const double      offset = (static_cast<double>(index) - static_cast<double>(d.dimensions[i] / 2)) * d.spacing[i];
                    rk.set(d.offset_variables[i], offset);
                }
                d.resolution_grid[flat] = d.resolution_pdf->evaluate_linear();
            }
            d.engine->set_resolution(d.resolution_grid);
        }

        const std::vector<double> & convolved = d.engine->convolve(d.signal_values);

        // Interpolate the smeared grid at the query point.
        std::vector<double> point(rank);
        for (std::size_t i = 0; i < rank; ++i)
        {
            point[i] = d.kinematics[d.axes[i].variable].evaluate();
        }

        const double value = d.engine->interpolate(convolved, point);

        return value > 0.0 ? value : 0.0;
    }

    double
    DetectorLevelPDF::evaluate() const
    {
        const double value = this->evaluate_linear();

        if (value > 0.0) [[likely]]
        {
            return std::log(value);
        }

        return -std::numeric_limits<double>::infinity();
    }

    double
    DetectorLevelPDF::normalization() const
    {
        // The convolution preserves the truth PDF's normalisation (the kernel is renormalised to
        // unit sum), so the normalisation is that of the truth PDF. It depends only on the [min, max]
        // bounds, which are identical for every grid point, so any grid point's PDF may be used.
        return _data->signal_pdfs.front()->normalization();
    }

    ObservablePtr
    DetectorLevelPDF::unnormalized_pdf() const
    {
        // A detector-level PDF is a convolution over the whole grid of truth-PDF observables, not a
        // single observable, so it cannot expose one via this interface (which the ObservableCache
        // uses to batch a PDF's linear value). Throw for the time being.
        throw InternalError("DetectorLevelPDF::unnormalized_pdf: a detector-level PDF is not backed by a single observable");
    }

    Kinematics
    DetectorLevelPDF::kinematics()
    {
        return _data->kinematics;
    }

    Parameters
    DetectorLevelPDF::parameters()
    {
        return _data->parameters;
    }

    Options
    DetectorLevelPDF::options()
    {
        return _data->options;
    }

    DensityPtr
    DetectorLevelPDF::clone() const
    {
        return this->clone(_data->parameters.clone());
    }

    DensityPtr
    DetectorLevelPDF::clone(const Parameters & parameters) const
    {
        const ObservableCache cache(parameters);

        if (_data->has_resolution_pdf)
        {
            return DensityPtr(new DetectorLevelPDF(cache, _data->signal_name, _data->resolution_name, _data->options, _data->axes));
        }

        return DensityPtr(new DetectorLevelPDF(cache, _data->signal_name, _data->options, _data->axes.front(), _data->resolution_grid));
    }

    Density::Iterator
    DetectorLevelPDF::begin() const
    {
        return Density::Iterator(_data->descriptions.cbegin());
    }

    Density::Iterator
    DetectorLevelPDF::end() const
    {
        return Density::Iterator(_data->descriptions.cend());
    }
} // namespace eos
