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

#include <eos/maths/cauchy-transform.hh>
#include <eos/observable.hh>
#include <eos/utils/detector-level-pdf.hh>
#include <eos/utils/expression-kernel-decomposer.hh>
#include <eos/utils/expression-observable.hh>
#include <eos/utils/expression-parser.hh>
#include <eos/utils/log.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <functional>
#include <limits>
#include <span>

namespace eos
{
    namespace
    {
        // Largest relative difference between the grid sum of a kernel resolution and the integral of its kernels;
        // power-law tails, such as those of Kernel::CrystalBall, reach beyond any practical grid.
        constexpr double kernel_normalization_tolerance = 1.0e-4;

        // Number of grid points used to extract a subtraction coefficient from the grid.
        constexpr std::size_t pole_stencil_points = 7;

        /*
         * The subtraction coefficient of a pole at z = x_p - i gamma, from one line of the truth PDF p on a uniform grid.
         * The function g(x) = |x - z|^2 p(x) is free of the pole and varies only on scales much larger than gamma.
         * At x = conj(z), the analytic continuation of g yields c = (i / gamma) conj(g(conj(z))). We continue g by
         * Lagrange interpolation through the grid points nearest to x_p, which must not straddle an edge of the
         * phase space or another pole.
         */
        complex<double>
        coefficient_from_grid(std::span<const double> p, const std::size_t & base, const std::size_t & stride, const std::size_t & points, const double & origin,
                              const double & spacing, const complex<double> & z)
        {
            const double    gamma = -z.imag();
            const long long first = std::llround((z.real() - origin) / spacing) - static_cast<long long>(pole_stencil_points / 2);
            if ((first < 0) || (static_cast<std::size_t>(first) + pole_stencil_points > points))
            {
                throw InternalError(std::format("DetectorLevelPDF: the pole at {} lies too close to an end of the grid to extract its coefficient", z.real()));
            }

            std::array<double, pole_stencil_points> t;
            for (std::size_t a = 0; a < pole_stencil_points; ++a)
            {
                t[a] = origin + static_cast<double>(first + static_cast<long long>(a)) * spacing - z.real();
            }

            complex<double> g = 0.0;
            for (std::size_t a = 0; a < pole_stencil_points; ++a)
            {
                complex<double> weight = std::norm(t[a] + 1.0i * gamma) * p[base + (static_cast<std::size_t>(first) + a) * stride];
                for (std::size_t b = 0; b < pole_stencil_points; ++b)
                {
                    if (b != a)
                    {
                        weight *= (1.0i * gamma - t[b]) / (t[a] - t[b]);
                    }
                }
                g += weight;
            }

            return 1.0i / gamma * std::conj(g);
        }

        // The integral over [a, b] of the hat function of width 2 h centred on the grid node x, i.e. the weight of
        // that node in the integral of the linear interpolation.
        double
        hat_integral(const double & a, const double & b, const double & x, const double & h)
        {
            const auto antiderivative = [](const double & u) -> double
            {
                if (u <= -1.0)
                {
                    return 0.0;
                }
                if (u <= 0.0)
                {
                    return 0.5 * (1.0 + u) * (1.0 + u);
                }
                if (u <= 1.0)
                {
                    return 1.0 - 0.5 * (1.0 - u) * (1.0 - u);
                }
                return 1.0;
            };

            return h * (antiderivative((b - x) / h) - antiderivative((a - x) / h));
        }

        // Nodes and weights of the four-point Gauss-Legendre rule on [-1, 1].
        constexpr std::array<double, 4> gauss_legendre_nodes{ -0.8611363115940526, -0.3399810435848563, +0.3399810435848563, +0.8611363115940526 };
        constexpr std::array<double, 4> gauss_legendre_weights{ 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538 };
    } // namespace

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
            std::vector<Pole> poles;

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

            // Per pole: the index of its axis, and the cache entries of its position and, if given, of
            // its coefficient, one batch entry per grid line along that axis in row-major order.
            struct PoleIds
            {
                    std::size_t                             axis;
                    ObservableCache::ObservableId           position_re;
                    ObservableCache::ObservableId           position_im;
                    std::optional<ObservableCache::BatchId> coefficient_re{};
                    std::optional<ObservableCache::BatchId> coefficient_im{};

                    // The subtraction coefficients of the last update, one per grid line, given or extracted from the grid.
                    std::vector<complex<double>> coefficients{};

                    // The pole position of the last update, and per kernel term the coefficients of its Cauchy transform
                    // in the add-back, one per grid line: (c * w_t), rescaled like the grid resolution.
                    complex<double>                  position{};
                    std::vector<std::vector<double>> addback_re{}, addback_im{};

                    // The resolution as a linear combination of kernels in the offset variable of `axis`, the
                    // batches of their coefficients, sampled on the offsets of the other axes, and of their
                    // parameters (flat, in term order), and their Cauchy transforms; batches, since the cache would
                    // merge the equally named observables. For rank > 1, the coefficients are convolved over the other axes.
                    std::vector<exp::KernelTerm>                        kernel_terms{};
                    std::optional<ObservableCache::BatchId>             kernel_coefficients{};
                    std::optional<ObservableCache::BatchId>             kernel_parameters{};
                    std::vector<std::unique_ptr<CauchyTransform>>       cauchy_transforms{};
                    std::vector<std::unique_ptr<ResolutionConvolution>> line_convolutions{};

                    // whether the width of the pole has been reported as too large for a first-order Cauchy transform
                    bool width_reported{};
            };

            std::vector<PoleIds> pole_ids;

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

            // Scratch buffer for the add-back: one value per grid line along a pole axis.
            std::vector<double> line_values;

            // The add-back of the convolved singular terms at a point, which is evaluated there rather than
            // interpolated from the grid; the coefficients are interpolated over the axes other than the pole's.
            // Called concurrently for distinct events, hence the per-thread scratch.
            double
            addback(std::span<const double> point) const
            {
                thread_local std::vector<double> other_point;

                double result = 0.0;
                for (const auto & ids : pole_ids)
                {
                    other_point.clear();
                    for (std::size_t i = 0; i < point.size(); ++i)
                    {
                        if (i != ids.axis)
                        {
                            other_point.push_back(point[i]);
                        }
                    }

                    for (std::size_t t = 0; t < ids.cauchy_transforms.size(); ++t)
                    {
                        const complex<double> c  = ids.line_convolutions.empty() ? complex<double>(ids.addback_re[t][0], ids.addback_im[t][0])
                                                                                 : complex<double>(ids.line_convolutions[t]->interpolate(ids.addback_re[t], other_point),
                                                                                                  ids.line_convolutions[t]->interpolate(ids.addback_im[t], other_point));
                        result                  += (c * ids.cauchy_transforms[t]->evaluate(point[ids.axis] - ids.position)).real();
                    }
                }

                return result;
            }

            // The smeared grid exposed via grid(); allocated once, so spans and Interpolators over it stay valid.
            std::vector<double> smeared;

            // The first grid node along `axis` whose hat function overlaps [a, b], and the weights of the nodes from there on.
            std::pair<std::size_t, std::vector<double>>
            axis_weights(const std::size_t & axis, const double & a, const double & b) const
            {
                const double      h     = spacing[axis];
                const std::size_t first = static_cast<std::size_t>(std::max(0.0, std::floor((a - origin[axis]) / h)));
                const std::size_t last  = std::min(dimensions[axis] - 1, static_cast<std::size_t>(std::ceil((b - origin[axis]) / h)));

                std::vector<double> weights;
                weights.reserve(last - first + 1);
                for (std::size_t k = first; k <= last; ++k)
                {
                    weights.push_back(hat_integral(a, b, origin[axis] + static_cast<double>(k) * h, h));
                }

                return { first, std::move(weights) };
            }

            // The integral of the multilinear interpolation of the smeared grid over one window, given its axis_weights().
            double
            grid_integral(const std::vector<std::pair<std::size_t, std::vector<double>>> & window_weights) const
            {
                const std::size_t rank = axes.size();

                double                   result = 0.0;
                std::vector<std::size_t> index(rank, 0);
                while (true)
                {
                    std::size_t flat   = 0;
                    double      weight = 1.0;
                    for (std::size_t i = 0; i < rank; ++i)
                    {
                        const auto & [first, weights]  = window_weights[i];
                        flat                          += (first + index[i]) * strides[i];
                        weight                        *= weights[index[i]];
                    }
                    result += weight * smeared[flat];

                    std::size_t i = rank;
                    for (; i > 0; --i)
                    {
                        if (++index[i - 1] < window_weights[i - 1].second.size())
                        {
                            break;
                        }
                        index[i - 1] = 0;
                    }
                    if (0 == i)
                    {
                        break;
                    }
                }

                return result;
            }

            // The integral of the add-back over one window, given its axis_weights(). Along the pole axis, the Cauchy transforms vary on the scale of
            // the resolution, which the grid resolves, so a Gauss-Legendre rule on panels of one grid spacing suffices. Along
            // the other axes, the multilinear interpolation of the coefficients is integrated exactly.
            double
            addback_integral(const Window & window, const std::vector<std::pair<std::size_t, std::vector<double>>> & window_weights) const
            {
                const std::size_t rank = axes.size();

                double result = 0.0;
                for (const auto & ids : pole_ids)
                {
                    // the weight of each coefficient line, numbered row-major over the other axes
                    std::vector<std::vector<double>> weights(rank);
                    for (std::size_t i = 0; i < rank; ++i)
                    {
                        if (i != ids.axis)
                        {
                            const auto & [first, w] = window_weights[i];
                            weights[i].assign(dimensions[i], 0.0);
                            std::copy(w.cbegin(), w.cend(), weights[i].begin() + first);
                        }
                    }

                    std::vector<double> line_weights(ids.coefficients.size(), 1.0);
                    for (std::size_t line = 0; line < line_weights.size(); ++line)
                    {
                        std::size_t remainder = line;
                        for (std::size_t i = rank; i-- > 0;)
                        {
                            if (i == ids.axis)
                            {
                                continue;
                            }

                            line_weights[line] *= weights[i][remainder % dimensions[i]];
                            remainder          /= dimensions[i];
                        }
                    }

                    const auto [a, b]   = window[ids.axis];
                    const auto   panels = static_cast<std::size_t>(std::max(1.0, std::ceil((b - a) / spacing[ids.axis])));
                    const double width  = (b - a) / static_cast<double>(panels);

                    for (std::size_t t = 0; t < ids.cauchy_transforms.size(); ++t)
                    {
                        complex<double> c = 0.0;
                        if (ids.line_convolutions.empty())
                        {
                            c = complex<double>(ids.addback_re[t][0], ids.addback_im[t][0]);
                        }
                        else
                        {
                            for (std::size_t line = 0; line < line_weights.size(); ++line)
                            {
                                c += line_weights[line] * complex<double>(ids.addback_re[t][line], ids.addback_im[t][line]);
                            }
                        }

                        complex<double> S = 0.0;
                        for (std::size_t panel = 0; panel < panels; ++panel)
                        {
                            const double centre = a + (static_cast<double>(panel) + 0.5) * width;
                            for (std::size_t n = 0; n < gauss_legendre_nodes.size(); ++n)
                            {
                                S += gauss_legendre_weights[n] * ids.cauchy_transforms[t]->evaluate(centre + 0.5 * width * gauss_legendre_nodes[n] - ids.position);
                            }
                        }

                        result += (c * S).real() * 0.5 * width;
                    }
                }

                return result;
            }

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
                                       ResolutionKind kind, const std::vector<double> & resolution, const std::string & resolution_expression, const std::vector<Pole> & poles) :
        _data(new Data(cache, signal_name))
    {
        Data & d = *_data;

        if (axes.empty())
        {
            throw InternalError("DetectorLevelPDF: at least one axis is required");
        }

        // Kinematics::declare() overwrites a repeated name, which would collapse two axes into one.
        for (std::size_t i = 0; i < axes.size(); ++i)
        {
            for (std::size_t j = 0; j < i; ++j)
            {
                if (axes[i].variable == axes[j].variable)
                {
                    throw InternalError(std::format("DetectorLevelPDF: variable '{}' names more than one axis", axes[i].variable));
                }

                const std::string & offset_i = axes[i].offset_variable.empty() ? axes[i].variable : axes[i].offset_variable;
                const std::string & offset_j = axes[j].offset_variable.empty() ? axes[j].variable : axes[j].offset_variable;
                if (offset_i == offset_j)
                {
                    throw InternalError(std::format("DetectorLevelPDF: offset variable '{}' belongs to more than one axis", offset_i));
                }
            }
        }

        d.options               = options;
        d.axes                  = axes;
        d.kind                  = kind;
        d.resolution_expression = resolution_expression;
        d.poles                 = poles;

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

        // Construct the pole observables before registering anything with the cache.
        struct PoleObservables
        {
                std::size_t                                         axis;
                ObservablePtr                                       position_re;
                ObservablePtr                                       position_im;
                std::vector<ObservablePtr>                          coefficient_re{};
                std::vector<ObservablePtr>                          coefficient_im{};
                std::vector<exp::KernelTerm>                        kernel_terms{};
                std::vector<ObservablePtr>                          kernel_coefficients{};
                std::vector<ObservablePtr>                          kernel_parameters{};
                std::vector<std::unique_ptr<CauchyTransform>>       cauchy_transforms{};
                std::vector<std::unique_ptr<ResolutionConvolution>> line_convolutions{};
        };

        std::vector<PoleObservables> pole_observables;
        for (const auto & pole : poles)
        {
            // the coefficient of a pole in another variable would vary on the scale of the first pole's width along its axis
            if (pole.variable != poles.front().variable)
            {
                throw InternalError(
                        std::format("DetectorLevelPDF: all poles must lie in the same sampling variable, but lie in '{}' and '{}'", poles.front().variable, pole.variable));
            }

            const auto axis = std::find_if(axes.cbegin(), axes.cend(), [&](const Axis & a) { return a.variable == pole.variable; });
            if (axis == axes.cend())
            {
                throw InternalError(std::format("DetectorLevelPDF: pole variable '{}' is not a sampling variable", pole.variable));
            }

            if (pole.coefficient_re.has_value() != pole.coefficient_im.has_value())
            {
                throw InternalError(std::format("DetectorLevelPDF: pole in '{}' must provide both or neither of coefficient_re and coefficient_im", pole.variable));
            }

            auto make_observable = [&](const QualifiedName & name, const Kinematics & k) -> ObservablePtr
            {
                ObservablePtr result;
                try
                {
                    result = Observable::make(name, d.parameters, k, options);
                }
                catch (const UnknownObservableError &)
                {
                }

                if (! result)
                {
                    throw InternalError("DetectorLevelPDF: '" + name.full() + "' is neither a known observable nor a parameter");
                }

                return result;
            };

            PoleObservables o{ static_cast<std::size_t>(axis - axes.cbegin()), make_observable(pole.position_re, Kinematics()), make_observable(pole.position_im, Kinematics()) };

            if (pole.coefficient_re.has_value())
            {
                // One grid line along the pole axis per combination of the other axes' indices.
                const std::size_t lines = N / d.dimensions[o.axis];
                for (std::size_t line = 0; line < lines; ++line)
                {
                    Kinematics  k;
                    std::size_t remainder = line;
                    for (std::size_t i = rank; i-- > 0;)
                    {
                        if (i == o.axis)
                        {
                            continue;
                        }

                        const std::size_t index  = remainder % d.dimensions[i];
                        remainder               /= d.dimensions[i];
                        k.declare(axes[i].variable, d.origin[i] + static_cast<double>(index) * d.spacing[i]);
                    }

                    o.coefficient_re.push_back(make_observable(*pole.coefficient_re, k));
                    o.coefficient_im.push_back(make_observable(*pole.coefficient_im, k));
                }
            }

            // The add-back needs the resolution as a linear combination of kernels in this pole's offset variable.
            const std::string & offset_variable = d.offset_variables[o.axis];
            try
            {
                o.kernel_terms = exp::ExpressionKernelDecomposer(offset_variable).decompose(resolution_expression_ptr);
            }
            catch (const exp::ExpressionError & e)
            {
                throw InternalError(std::format("DetectorLevelPDF: the resolution is not a linear combination of kernels in '{}': {}", offset_variable, e.what()));
            }

            // The kernel coefficients are sampled on the centred offsets of the other axes, row-major; the kernel
            // parameters must not depend on them.
            std::vector<ResolutionConvolution::AxisGeometry> other_axes;
            for (std::size_t i = 0; i < rank; ++i)
            {
                if (i != o.axis)
                {
                    other_axes.push_back(ResolutionConvolution::AxisGeometry{ d.origin[i], d.spacing[i], d.dimensions[i] });
                }
            }
            const std::size_t lines  = N / d.dimensions[o.axis];
            const double      extent = 2.0 * (axes[o.axis].max - axes[o.axis].min);
            for (const auto & term : o.kernel_terms)
            {
                for (std::size_t line = 0; line < lines; ++line)
                {
                    Kinematics  k         = rk.clone();
                    std::size_t remainder = line;
                    for (std::size_t i = rank; i-- > 0;)
                    {
                        if (i == o.axis)
                        {
                            continue;
                        }

                        const std::size_t index  = remainder % d.dimensions[i];
                        remainder               /= d.dimensions[i];
                        k.set(d.offset_variables[i], (static_cast<double>(index) - static_cast<double>(d.dimensions[i] / 2)) * d.spacing[i]);
                    }

                    o.kernel_coefficients.push_back(
                            ObservablePtr(new ExpressionObservable(QualifiedName("DetectorLevelPDF::kernel_coefficient"), d.parameters, k, options, term.coefficient)));
                }
                if (1 < rank)
                {
                    o.line_convolutions.push_back(ResolutionConvolution::make(other_axes));
                }
                for (const auto & parameter : term.parameters)
                {
                    o.kernel_parameters.push_back(
                            ObservablePtr(new ExpressionObservable(QualifiedName("DetectorLevelPDF::kernel_parameter"), d.parameters, Kinematics(), options, parameter)));
                }
                o.cauchy_transforms.push_back(CauchyTransform::make(term.kernel, extent));
            }

            pole_observables.push_back(std::move(o));
        }

        // A single signal PDF provides the normalization, which is identical at every grid point (same [min, max]
        // bounds), and the name and options of its numerator, which is then made directly at every grid point
        // and registered with the cache as one batch.
        SignalPDFPtr signal_pdf;
        {
            Kinematics k;
            for (std::size_t i = 0; i < rank; ++i)
            {
                k.declare(axes[i].variable, d.origin[i]);
                k.declare(axes[i].variable + "_min", axes[i].min);
                k.declare(axes[i].variable + "_max", axes[i].max);
            }

            signal_pdf = SignalPDF::make(signal_name, d.parameters, k, options);
            if (! signal_pdf.get())
            {
                throw InternalError("DetectorLevelPDF: '" + signal_name.str() + "' is not a valid signal PDF name");
            }
        }
        const ObservablePtr normalization_observable = signal_pdf->normalization_observable();
        const ObservablePtr numerator                = signal_pdf->unnormalized_pdf();
        const QualifiedName numerator_name           = numerator->name();
        const Options       numerator_options        = numerator->options();

        std::vector<ObservablePtr> unnormalized_pdfs;
        unnormalized_pdfs.reserve(N);
        for (std::size_t flat = 0; flat < N; ++flat)
        {
            Kinematics k;
            for (std::size_t i = 0; i < rank; ++i)
            {
                const std::size_t index = (flat / d.strides[i]) % d.dimensions[i];

                k.declare(axes[i].variable, d.origin[i] + static_cast<double>(index) * d.spacing[i]);
            }

            unnormalized_pdfs.push_back(Observable::make(numerator_name, d.parameters, k, numerator_options));
        }
        d.batch_id         = d.cache.add_batch(std::move(unnormalized_pdfs));
        d.normalization_id = d.cache.add(normalization_observable);

        for (auto & o : pole_observables)
        {
            Data::PoleIds ids{ o.axis, d.cache.add(o.position_re), d.cache.add(o.position_im) };
            if (! o.coefficient_re.empty())
            {
                ids.coefficient_re = d.cache.add_batch(std::move(o.coefficient_re));
                ids.coefficient_im = d.cache.add_batch(std::move(o.coefficient_im));
            }
            ids.kernel_terms        = std::move(o.kernel_terms);
            ids.kernel_coefficients = d.cache.add_batch(std::move(o.kernel_coefficients));
            ids.kernel_parameters   = d.cache.add_batch(std::move(o.kernel_parameters));
            ids.cauchy_transforms   = std::move(o.cauchy_transforms);
            ids.line_convolutions   = std::move(o.line_convolutions);
            ids.coefficients.assign(N / d.dimensions[o.axis], 0.0);
            ids.addback_re.assign(ids.cauchy_transforms.size(), std::vector<double>(N / d.dimensions[o.axis], 0.0));
            ids.addback_im.assign(ids.cauchy_transforms.size(), std::vector<double>(N / d.dimensions[o.axis], 0.0));
            d.pole_ids.push_back(std::move(ids));
        }

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
        DetectorLevelPDF(cache, signal_name, options, axes, ResolutionKind::grid, resolution, std::string(), std::vector<Pole>{})
    {
    }

    DetectorLevelPDF::DetectorLevelPDF(const ObservableCache & cache, const QualifiedName & signal_name, const Options & options, const std::vector<Axis> & axes,
                                       const std::string & resolution, const std::vector<Pole> & poles) :
        DetectorLevelPDF(cache, signal_name, options, axes, ResolutionKind::expression, std::vector<double>{}, resolution, poles)
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
                                           const std::string & resolution, const std::vector<Pole> & poles)
    {
        return SignalPDFPtr(new DetectorLevelPDF(cache, signal_name, options, axes, resolution, poles));
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

        // The flat index of point k on a grid line along an axis, with the line numbered row-major over the other axes.
        const auto flat_index = [](const std::size_t & line, const std::size_t & k, const std::size_t & points, const std::size_t & stride) -> std::size_t
        { return (line / stride) * points * stride + k * stride + line % stride; };
        // The volume of a grid cell over all axes but one.
        const auto other_volume = [&d](const std::size_t & axis) -> double
        {
            double result = 1.0;
            for (std::size_t i = 0; i < d.spacing.size(); ++i)
            {
                result *= (i == axis) ? 1.0 : d.spacing[i];
            }
            return result;
        };

        for (const auto & ids : d.pole_ids)
        {
            // also rejects NaN
            if (! (d.cache[ids.position_im] < 0.0))
            {
                throw InternalError(
                        std::format("DetectorLevelPDF: pole in '{}' has Im z_p = {}, but must lie in the lower half plane", d.axes[ids.axis].variable, d.cache[ids.position_im]));
            }
        }

        const std::span<const double> batch = d.cache[d.batch_id];
        for (std::size_t flat = 0; flat < N; ++flat)
        {
            // argument order matters: std::max(NaN, 0.0) keeps the NaN, which then invalidates the whole grid
            d.signal_values[flat] = std::max(batch[flat], 0.0);
        }

        // The grid cannot resolve the pole terms Re[c / (x - z)]; they are subtracted here, one line along
        // the pole axis at a time, and their convolution is added back analytically below.
        for (auto & ids : d.pole_ids)
        {
            const complex<double> z(d.cache[ids.position_re], d.cache[ids.position_im]);
            const std::size_t     points = d.dimensions[ids.axis], stride = d.strides[ids.axis];
            for (std::size_t line = 0; line < ids.coefficients.size(); ++line)
            {
                const std::size_t base = flat_index(line, 0, points, stride);
                // the truth PDF is assumed non-negative near a pole, where the clipping would otherwise bias the extracted coefficient
                ids.coefficients[line] = ids.coefficient_re.has_value() ? complex<double>(d.cache[*ids.coefficient_re][line], d.cache[*ids.coefficient_im][line])
                                                                        : coefficient_from_grid(batch, base, stride, points, d.origin[ids.axis], d.spacing[ids.axis], z);
                for (std::size_t k = 0; k < points; ++k)
                {
                    const double x                      = d.origin[ids.axis] + static_cast<double>(k) * d.spacing[ids.axis];
                    d.signal_values[base + k * stride] -= (ids.coefficients[line] / (x - z)).real();
                }
            }
        }

        // An expression resolution may depend on parameters, so its kernel is re-sampled on every
        // convolution; a pre-computed grid is fixed, so its spectrum stays cached from construction.
        if (ResolutionKind::expression == d.kind)
        {
            const std::span<const double> kernel = d.cache[d.resolution_batch_id];
            std::copy(kernel.begin(), kernel.end(), d.resolution_grid.begin());

            // the grid sum must match the integral of the kernels, or the add-back would use a different kernel
            for (const auto & ids : d.pole_ids)
            {
                double norm = 0.0;
                for (const double coefficient : d.cache[*ids.kernel_coefficients])
                {
                    norm += coefficient;
                }
                norm *= other_volume(ids.axis);

                double grid_norm = 0.0;
                for (const double value : d.resolution_grid)
                {
                    grid_norm += value;
                }
                grid_norm *= other_volume(ids.axis) * d.spacing[ids.axis];

                if (! (std::abs(grid_norm - norm) <= kernel_normalization_tolerance * std::abs(norm)))
                {
                    throw InternalError(
                            std::format("DetectorLevelPDF: the resolution sums to {} on the grid, but its kernels integrate to {}; the grid is too narrow or too coarse",
                                        grid_norm,
                                        norm));
                }
            }

            d.convolution->set_resolution(d.resolution_grid);
        }

        const std::vector<double> & result = d.convolution->convolve(d.signal_values);
        std::copy(result.cbegin(), result.cend(), d.smeared.begin());

        // The add-back sum_t Re[(c * w_t)(y) S_t(x - z)] is evaluated at the query points, see Data::addback(). Here, the
        // coefficient lines are convolved with each term's coefficient w_t over the other axes, and rescaled like the grid resolution.
        for (auto & ids : d.pole_ids)
        {
            ids.position = complex<double>(d.cache[ids.position_re], d.cache[ids.position_im]);

            const std::size_t             lines      = ids.coefficients.size();
            const std::span<const double> weights    = d.cache[*ids.kernel_coefficients];
            const std::span<const double> parameters = d.cache[*ids.kernel_parameters];

            double norm = 0.0;
            for (const double weight : weights)
            {
                norm += weight;
            }

            d.line_values.resize(lines);
            std::size_t offset = 0;
            for (std::size_t t = 0; t < ids.cauchy_transforms.size(); ++t)
            {
                const std::size_t arity = ids.kernel_terms[t].parameters.size();
                ids.cauchy_transforms[t]->prepare(parameters.subspan(offset, arity));
                // the Cauchy transform of Kernel::CrystalBall is first order in gamma / sigma
                if (("Kernel::CrystalBall" == ids.kernel_terms[t].kernel) && (! ids.width_reported) && (-ids.position.imag() > 0.1 * parameters[offset + 1]))
                {
                    ids.width_reported = true;
                    Log::instance()->message("DetectorLevelPDF::update_grid", ll_warning)
                            << "The pole in '" << d.axes[ids.axis].variable << "' has the width " << -ids.position.imag() << ", more than 0.1 times the width "
                            << parameters[offset + 1] << " of the Crystal Ball resolution; its convolution neglects terms of second order in their ratio";
                }
                offset += arity;

                const std::span<const double> w     = weights.subspan(t * lines, lines);
                double                        w_sum = 0.0;
                for (const double value : w)
                {
                    w_sum += value;
                }

                const double scale = w_sum / norm;
                if (0.0 == w_sum)
                {
                    std::fill(ids.addback_re[t].begin(), ids.addback_re[t].end(), 0.0);
                    std::fill(ids.addback_im[t].begin(), ids.addback_im[t].end(), 0.0);
                }
                else if (ids.line_convolutions.empty())
                {
                    ids.addback_re[t][0] = scale * ids.coefficients[0].real();
                    ids.addback_im[t][0] = scale * ids.coefficients[0].imag();
                }
                else
                {
                    // the convolution renormalizes w_t to unit sum; its share w_sum / norm of the resolution enters as scale
                    auto & convolution = *ids.line_convolutions[t];
                    std::transform(w.begin(), w.end(), d.line_values.begin(), [&w_sum](const double & value) { return w_sum > 0.0 ? value : -value; });
                    convolution.set_resolution(d.line_values);

                    std::transform(ids.coefficients.begin(), ids.coefficients.end(), d.line_values.begin(), [](const complex<double> & c) { return c.real(); });
                    const std::vector<double> & re = convolution.convolve(d.line_values);
                    std::transform(re.begin(), re.end(), ids.addback_re[t].begin(), [&scale](const double & value) { return scale * value; });

                    std::transform(ids.coefficients.begin(), ids.coefficients.end(), d.line_values.begin(), [](const complex<double> & c) { return c.imag(); });
                    const std::vector<double> & im = convolution.convolve(d.line_values);
                    std::transform(im.begin(), im.end(), ids.addback_im[t].begin(), [&scale](const double & value) { return scale * value; });
                }
            }
        }

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

    const std::vector<DetectorLevelPDF::Pole> &
    DetectorLevelPDF::poles() const
    {
        return _data->poles;
    }

    namespace
    {
        // Interpolates the smeared remainder from the grid, and adds the convolved singular terms at each point.
        class PoleInterpolator : public ResolutionConvolution::Interpolator
        {
            private:
                std::unique_ptr<ResolutionConvolution::Interpolator> _remainder;
                std::function<double(std::span<const double>)>       _addback;
                std::vector<double>                                  _points;
                std::size_t                                          _rank;

            public:
                PoleInterpolator(std::unique_ptr<ResolutionConvolution::Interpolator> && remainder, const std::function<double(std::span<const double>)> & addback,
                                 std::span<const double> points, const std::size_t & rank) :
                    _remainder(std::move(remainder)),
                    _addback(addback),
                    _points(points.begin(), points.end()),
                    _rank(rank)
                {
                }

                virtual std::size_t
                size() const override
                {
                    return _remainder->size();
                }

                virtual double
                operator() (std::size_t i) const override
                {
                    return (*_remainder)(i) + _addback(std::span<const double>(_points.data() + i * _rank, _rank));
                }
        };
    } // namespace

    std::unique_ptr<ResolutionConvolution::Interpolator>
    DetectorLevelPDF::make_interpolator(std::span<const double> points) const
    {
        auto remainder = _data->convolution->make_interpolator(this->grid(), points);
        if (_data->pole_ids.empty())
        {
            return remainder;
        }

        Data * d = _data.get();
        return std::make_unique<PoleInterpolator>(std::move(remainder), [d](std::span<const double> point) { return d->addback(point); }, points, _data->axes.size());
    }

    void
    DetectorLevelPDF::validate(const std::vector<Window> & windows) const
    {
        const std::vector<Axis> & axes = _data->axes;

        if (windows.empty())
        {
            throw InternalError("DetectorLevelPDF: at least one window is required");
        }

        for (const auto & window : windows)
        {
            if (window.size() != axes.size())
            {
                throw InternalError(std::format("DetectorLevelPDF: a window has {} intervals, but the PDF has rank {}", window.size(), axes.size()));
            }

            for (std::size_t i = 0; i < axes.size(); ++i)
            {
                // also rejects NaN
                const auto & [a, b] = window[i];
                if (! ((axes[i].min <= a) && (a < b) && (b <= axes[i].max)))
                {
                    throw InternalError(
                            std::format("DetectorLevelPDF: window [{}, {}] in '{}' is empty or exceeds the grid [{}, {}]", a, b, axes[i].variable, axes[i].min, axes[i].max));
                }
            }
        }

        for (std::size_t w = 0; w < windows.size(); ++w)
        {
            for (std::size_t v = w + 1; v < windows.size(); ++v)
            {
                bool overlap = true;
                for (std::size_t i = 0; i < axes.size(); ++i)
                {
                    overlap = overlap && (std::max(windows[w][i].first, windows[v][i].first) < std::min(windows[w][i].second, windows[v][i].second));
                }

                if (overlap)
                {
                    throw InternalError(std::format("DetectorLevelPDF: windows {} and {} overlap", w, v));
                }
            }
        }
    }

    double
    DetectorLevelPDF::integral(const std::vector<Window> & windows) const
    {
        return (*this->make_window_integrator(windows))();
    }

    std::unique_ptr<DetectorLevelPDF::WindowIntegrator>
    DetectorLevelPDF::make_window_integrator(const std::vector<Window> & windows) const
    {
        this->validate(windows);

        std::vector<std::vector<WindowIntegrator::AxisWeights>> weights;
        weights.reserve(windows.size());
        for (const auto & window : windows)
        {
            std::vector<WindowIntegrator::AxisWeights> window_weights;
            window_weights.reserve(window.size());
            for (std::size_t i = 0; i < window.size(); ++i)
            {
                window_weights.push_back(_data->axis_weights(i, window[i].first, window[i].second));
            }
            weights.push_back(std::move(window_weights));
        }

        return std::unique_ptr<WindowIntegrator>(new WindowIntegrator(*this, windows, std::move(weights)));
    }

    DetectorLevelPDF::WindowIntegrator::WindowIntegrator(const DetectorLevelPDF & pdf, const std::vector<Window> & windows, std::vector<std::vector<AxisWeights>> && weights) :
        _pdf(pdf),
        _windows(windows),
        _weights(std::move(weights))
    {
    }

    double
    DetectorLevelPDF::WindowIntegrator::operator() () const
    {
        _pdf.update_grid();

        double result = 0.0;
        for (std::size_t w = 0; w < _weights.size(); ++w)
        {
            result += _pdf._data->grid_integral(_weights[w]) + _pdf._data->addback_integral(_windows[w], _weights[w]);
        }

        if (result > 0.0) [[likely]]
        {
            return std::log(result);
        }

        return -std::numeric_limits<double>::infinity();
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
        double smeared = d.convolution->interpolate(this->grid(), d.query_point);
        if (! d.pole_ids.empty())
        {
            smeared += d.addback(d.query_point);
        }

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

    SignalPDFPtr
    DetectorLevelPDF::clone(const ObservableCache & cache) const
    {
        // The clone has its own result buffer, so any span or Interpolator held over the original
        // PDF's grid() is invalid for the clone.
        switch (_data->kind)
        {
            case ResolutionKind::expression:
                return std::make_shared<DetectorLevelPDF>(cache, _data->signal_name, _data->options, _data->axes, _data->resolution_expression, _data->poles);
            default: return std::make_shared<DetectorLevelPDF>(cache, _data->signal_name, _data->options, _data->axes, _data->resolution_grid);
        }
    }
} // namespace eos
