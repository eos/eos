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

#include <eos/utils/density-impl.hh>
#include <eos/utils/detector-level-pdf.hh>
#include <eos/utils/expression-observable.hh>
#include <eos/utils/expression-parser.hh>

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <span>

namespace eos
{
    struct DetectorLevelPDF::Data
    {
            // Construction inputs, retained so that clone() can rebuild an independent copy.
            ObservableCache   cache;
            QualifiedName     name;
            QualifiedName     signal_name;
            Options           options;
            std::vector<Axis> axes;
            ResolutionKind    kind;
            std::string       resolution_expression;

            Parameters parameters;

            // Grid geometry (row-major), mirroring eos::ResolutionConvolution.
            std::vector<std::size_t> dimensions;
            std::vector<std::size_t> strides;
            std::vector<double>      origin;
            std::vector<double>      spacing;
            std::vector<std::string> offset_variables;

            // The convolution. Per-grid-point truth PDFs are not retained; the cache batch keeps
            // their observables alive.
            std::unique_ptr<ResolutionConvolution> convolution;

            // The batch of resolution samples on the offset grid, if the resolution is an expression.
            ObservableCache::BatchId resolution_batch_id;

            // The batch of per-grid-point truth-PDF observables, and the (identical at every grid
            // point) truth normalization, both registered with `cache` at construction.
            ObservableCache::BatchId      batch_id;
            ObservableCache::ObservableId normalization_id;

            // The cache generation at the last convolve(); update_grid() is a no-op while it is
            // still current. Default-constructed, it matches no generation of any cache.
            ObservableCache::Generation last_generation;

            // The query point exposed via kinematics(); interpolation reads its sampling variables.
            Kinematics kinematics;

            // Scratch buffers: the truth-PDF grid values (clamped, copied out of the cache), the
            // (re-)sampled resolution kernel, and the interpolation query point.
            std::vector<double> signal_values;
            std::vector<double> resolution_grid;
            std::vector<double> query_point;

            // The smeared grid exposed via grid(); allocated once, so spans and Interpolators over it stay valid.
            std::vector<double> smeared;

            // Empty; present so that Density::begin()/end() have something to iterate over.
            std::vector<ParameterDescription> descriptions;

            Data(const ObservableCache & cache, const QualifiedName & signal_name) :
                cache(cache),
                name(signal_name),
                signal_name(signal_name),
                kind(ResolutionKind::grid),
                parameters(cache.parameters())
            {
            }
    };

    DetectorLevelPDF::DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const std::vector<Axis> & axes,
                                       ResolutionKind kind, const std::vector<double> & resolution, const std::string & resolution_expression) :
        _data(new Data(cache, signal_name))
    {
        Data & d = *_data;

        if (axes.empty())
        {
            throw InternalError("DetectorLevelPDF: at least one axis is required");
        }

        d.options               = options;
        d.axes                  = axes;
        d.kind                  = kind;
        d.resolution_expression = resolution_expression;

        const std::size_t rank = axes.size();

        // Assemble the grid geometry and construct the convolution (which validates the
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
        d.convolution = ResolutionConvolution::make(geometry);

        // Row-major strides, mirroring eos::ResolutionConvolution.
        d.strides.resize(rank);
        d.strides[rank - 1] = 1;
        for (std::size_t i = rank - 1; i-- > 0;)
        {
            d.strides[i] = d.strides[i + 1] * d.dimensions[i + 1];
        }

        const std::size_t N = d.convolution->size();

        // The query point: one entry per sampling variable, initialised to the grid centre. The
        // [min, max] bounds are declared alongside (named '<variable>_min/_max') so that the PDF
        // presents the same variable/bound structure as any other SignalPDF.
        for (std::size_t i = 0; i < rank; ++i)
        {
            d.kinematics.declare(axes[i].variable, 0.5 * (axes[i].min + axes[i].max));
            d.kinematics.declare(axes[i].variable + "_min", axes[i].min);
            d.kinematics.declare(axes[i].variable + "_max", axes[i].max);
        }

        // The offset variables, declared with their bounds.
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

        // Validate the resolution before registering anything with the (possibly shared) cache,
        // so that a failed construction leaves the cache untouched.
        exp::ExpressionPtr         resolution_expression_ptr;
        std::vector<ObservablePtr> resolution_observables;
        if (ResolutionKind::expression == kind)
        {
            // Sample the resolution at every offset, as one batch registered with the cache below.
            resolution_expression_ptr = exp::parse_expression(resolution_expression);
            resolution_observables.reserve(N);
            for (std::size_t flat = 0; flat < N; ++flat)
            {
                Kinematics k = rk.clone();
                for (std::size_t i = 0; i < rank; ++i)
                {
                    const std::size_t index = (flat / d.strides[i]) % d.dimensions[i];
                    k.set(d.offset_variables[i], (static_cast<double>(index) - static_cast<double>(d.dimensions[i] / 2)) * d.spacing[i]);
                }
                resolution_observables.push_back(
                        ObservablePtr(new ExpressionObservable(QualifiedName("DetectorLevelPDF::resolution"), d.parameters, k, options, resolution_expression_ptr)));
            }
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
            d.convolution->set_resolution(d.resolution_grid);
        }

        // Register each grid point's truth PDF observable as one cache batch; the normalization is
        // identical at every grid point (same [min, max] bounds), so capture it once.
        std::vector<ObservablePtr> unnormalized_pdfs;
        unnormalized_pdfs.reserve(N);
        ObservablePtr normalization_observable;
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

            SignalPDFPtr signal_pdf = SignalPDF::make(signal_name, d.parameters, k, options);
            if (! signal_pdf.get())
            {
                throw InternalError("DetectorLevelPDF: '" + signal_name.str() + "' is not a valid signal PDF name");
            }
            unnormalized_pdfs.push_back(signal_pdf->unnormalized_pdf());

            if (! normalization_observable)
            {
                normalization_observable = signal_pdf->normalization_observable();
            }
        }
        d.batch_id         = d.cache.add_batch(std::move(unnormalized_pdfs));
        d.normalization_id = d.cache.add(normalization_observable);

        if (ResolutionKind::expression == kind)
        {
            d.resolution_batch_id = d.cache.add_batch(std::move(resolution_observables));
        }
        d.signal_values.assign(N, 0.0);
        d.query_point.assign(rank, 0.0);
        d.smeared.assign(N, 0.0);
    }

    DetectorLevelPDF::DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const std::vector<Axis> & axes,
                                       const std::vector<double> & resolution) :
        DetectorLevelPDF(cache, signal_name, options, axes, ResolutionKind::grid, resolution, std::string())
    {
    }

    DetectorLevelPDF::DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const std::vector<Axis> & axes,
                                       const std::string & resolution) :
        DetectorLevelPDF(cache, signal_name, options, axes, ResolutionKind::expression, std::vector<double>{}, resolution)
    {
    }

    DetectorLevelPDF::~DetectorLevelPDF() = default;

    SignalPDFPtr
    DetectorLevelPDF::make_from_grid(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const std::vector<Axis> & axes,
                                     const std::vector<double> & resolution)
    {
        return SignalPDFPtr(new DetectorLevelPDF(cache, signal_name, options, axes, resolution));
    }

    SignalPDFPtr
    DetectorLevelPDF::make_1d(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const Axis & axis, const std::vector<double> & resolution)
    {
        return make_from_grid(cache, signal_name, options, std::vector<Axis>{ axis }, resolution);
    }

    SignalPDFPtr
    DetectorLevelPDF::make_from_expression(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const std::vector<Axis> & axes,
                                           const std::string & resolution)
    {
        return SignalPDFPtr(new DetectorLevelPDF(cache, signal_name, options, axes, resolution));
    }

    const QualifiedName &
    DetectorLevelPDF::name() const
    {
        return _data->name;
    }

    void
    DetectorLevelPDF::update_grid() const
    {
        Data & d = *_data;

        if (d.cache.generation() == d.last_generation)
        {
            return;
        }

        const std::size_t N = d.convolution->size();

        const std::span<const double> batch = d.cache[d.batch_id];
        for (std::size_t flat = 0; flat < N; ++flat)
        {
            d.signal_values[flat] = std::max(0.0, batch[flat]);
        }

        // An expression resolution may depend on parameters, so its kernel is re-sampled on every
        // convolution; a pre-computed grid is fixed, so its spectrum stays cached from construction.
        if (ResolutionKind::expression == d.kind)
        {
            const std::span<const double> kernel = d.cache[d.resolution_batch_id];
            std::copy(kernel.begin(), kernel.end(), d.resolution_grid.begin());

            d.convolution->set_resolution(d.resolution_grid);
        }

        const std::vector<double> & result = d.convolution->convolve(d.signal_values);
        std::copy(result.cbegin(), result.cend(), d.smeared.begin());

        d.last_generation = d.cache.generation();
    }

    std::span<const double>
    DetectorLevelPDF::grid() const
    {
        return _data->smeared;
    }

    const ObservableCache &
    DetectorLevelPDF::cache() const
    {
        return _data->cache;
    }

    const std::vector<DetectorLevelPDF::Axis> &
    DetectorLevelPDF::axes() const
    {
        return _data->axes;
    }

    std::unique_ptr<ResolutionConvolution::Interpolator>
    DetectorLevelPDF::make_interpolator(std::span<const double> points) const
    {
        return _data->convolution->make_interpolator(this->grid(), points);
    }

    double
    DetectorLevelPDF::evaluate_linear() const
    {
        Data & d = *_data;

        this->update_grid();

        const std::size_t rank = d.axes.size();

        for (std::size_t i = 0; i < rank; ++i)
        {
            d.query_point[i] = d.kinematics[d.axes[i].variable].evaluate();
        }

        // interpolate() rather than a single-point Interpolator: this is the hot path (plotting,
        // sampling), and building an Interpolator per call would allocate.
        const double smeared = d.convolution->interpolate(this->grid(), d.query_point);

        return smeared > 0.0 ? smeared : 0.0;
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
        // As for every SignalPDF, this is the *logarithm* of the normalization. The convolution
        // preserves the truth PDF's normalization (the kernel has unit sum), and it is identical at
        // every grid point since it depends only on the [min, max] bounds.
        const double norm = _data->cache[_data->normalization_id];

        if (norm > 0.0) [[likely]]
        {
            return std::log(norm);
        }

        return -std::numeric_limits<double>::infinity();
    }

    ObservablePtr
    DetectorLevelPDF::unnormalized_pdf() const
    {
        // A detector-level PDF is not backed by a single observable, so it cannot expose one here;
        // this is also the only path into ObservableCache::add()/add_batch(), so the throw is deliberate.
        throw InternalError("DetectorLevelPDF::unnormalized_pdf: a detector-level PDF is not backed by a single observable");
    }

    ObservablePtr
    DetectorLevelPDF::normalization_observable() const
    {
        // The truth PDF's normalization observable, registered with the cache at construction; an
        // ordinary observable that legitimately belongs in a cache.
        return _data->cache.observable(_data->normalization_id);
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
        throw InternalError("DetectorLevelPDF::clone: use clone(const ObservableCache &) to bind the clone to a cache that is updated");
    }

    DensityPtr
    DetectorLevelPDF::clone(const Parameters &) const
    {
        throw InternalError("DetectorLevelPDF::clone: use clone(const ObservableCache &) to bind the clone to a cache that is updated");
    }

    std::shared_ptr<DetectorLevelPDF>
    DetectorLevelPDF::clone(const ObservableCache & cache) const
    {
        // The clone has its own result buffer, so any span or Interpolator held over the original
        // PDF's grid() is invalid for the clone.
        switch (_data->kind)
        {
            case ResolutionKind::expression: return std::make_shared<DetectorLevelPDF>(cache, _data->signal_name, _data->options, _data->axes, _data->resolution_expression);
            default:                         return std::make_shared<DetectorLevelPDF>(cache, _data->signal_name, _data->options, _data->axes, _data->resolution_grid);
        }
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
