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

#include <eos/maths/resolution-convolution.hh>
#include <eos/utils/exception.hh>

#include <test/test.hh>

#include <cmath>
#include <vector>

using namespace test;
using namespace eos;

// A grid axis with `points` points, spanning [0, points) with unit spacing (origin 0).
static ResolutionConvolution::AxisGeometry
unit_axis(std::size_t points)
{
    return ResolutionConvolution::AxisGeometry{ 0.0, 1.0, points };
}

// Discrete delta kernel at the centre index (points / 2); convolving with it is the identity.
static std::vector<double>
delta_kernel_1d(std::size_t points)
{
    std::vector<double> kernel(points, 0.0);
    kernel[points / 2] = 1.0;
    return kernel;
}

class ResolutionConvolutionIdentityTest : public TestCase
{
    public:
        ResolutionConvolutionIdentityTest() :
            TestCase("resolution_convolution_identity_test")
        {
        }

        virtual void
        run() const override
        {
            const std::size_t N = 8;

            auto engine = ResolutionConvolution::make({ unit_axis(N) });
            TEST_CHECK_EQUAL(engine->rank(), 1u);
            TEST_CHECK_EQUAL(engine->size(), N);

            engine->set_resolution(delta_kernel_1d(N));

            // Convolution with a delta at zero offset returns the signal unchanged.
            std::vector<double> signal{ 0.0, 1.0, 3.0, 7.0, 2.0, 0.5, 0.25, 0.0 };
            const auto &        convolved = engine->convolve(signal);

            for (std::size_t i = 0; i < N; ++i)
            {
                TEST_CHECK_NEARLY_EQUAL(convolved[i], signal[i], 1.0e-12);
            }
        }
} resolution_convolution_identity_test;

class ResolutionConvolutionNormalisationTest : public TestCase
{
    public:
        ResolutionConvolutionNormalisationTest() :
            TestCase("resolution_convolution_normalisation_test")
        {
        }

        virtual void
        run() const override
        {
            const std::size_t N = 8;

            auto engine = ResolutionConvolution::make({ unit_axis(N) });

            // Uniform, un-normalised kernel; renormalised to unit sum, it maps every point to the signal's mean.
            std::vector<double> kernel(N, 2.0);
            engine->set_resolution(kernel);

            std::vector<double> signal{ 0.0, 1.0, 3.0, 7.0, 2.0, 0.5, 0.25, 0.25 };
            double              mean = 0.0;
            for (const auto & s : signal)
            {
                mean += s;
            }
            mean /= static_cast<double>(N);

            const auto & convolved = engine->convolve(signal);
            for (std::size_t i = 0; i < N; ++i)
            {
                TEST_CHECK_NEARLY_EQUAL(convolved[i], mean, 1.0e-12);
            }
        }
} resolution_convolution_normalisation_test;

class ResolutionConvolutionConstantSignalTest : public TestCase
{
    public:
        ResolutionConvolutionConstantSignalTest() :
            TestCase("resolution_convolution_constant_signal_test")
        {
        }

        virtual void
        run() const override
        {
            const std::size_t N = 16;

            auto engine = ResolutionConvolution::make({ unit_axis(N) });

            // A symmetric three-point smoothing kernel, centred.
            std::vector<double> kernel(N, 0.0);
            kernel[N / 2 - 1] = 0.25;
            kernel[N / 2]     = 0.50;
            kernel[N / 2 + 1] = 0.25;
            engine->set_resolution(kernel);

            // A constant signal must be preserved exactly by any unit-sum kernel (up to roundoff).
            std::vector<double> signal(N, 3.5);
            const auto &        convolved = engine->convolve(signal);
            for (std::size_t i = 0; i < N; ++i)
            {
                TEST_CHECK_NEARLY_EQUAL(convolved[i], 3.5, 1.0e-12);
            }
        }
} resolution_convolution_constant_signal_test;

class ResolutionConvolutionShiftTest : public TestCase
{
    public:
        ResolutionConvolutionShiftTest() :
            TestCase("resolution_convolution_shift_test")
        {
        }

        virtual void
        run() const override
        {
            const std::size_t N = 8;

            auto engine = ResolutionConvolution::make({ unit_axis(N) });

            // Delta displaced by +1 from the centre; convolving circularly shifts the signal by +1.
            std::vector<double> kernel(N, 0.0);
            kernel[N / 2 + 1] = 1.0;
            engine->set_resolution(kernel);

            std::vector<double> signal{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
            const auto &        convolved = engine->convolve(signal);

            for (std::size_t i = 0; i < N; ++i)
            {
                TEST_CHECK_NEARLY_EQUAL(convolved[i], signal[(i + N - 1) % N], 1.0e-12);
            }
        }
} resolution_convolution_shift_test;

class ResolutionConvolutionInterpolateTest : public TestCase
{
    public:
        ResolutionConvolutionInterpolateTest() :
            TestCase("resolution_convolution_interpolate_test")
        {
        }

        virtual void
        run() const override
        {
            const std::size_t N = 8;

            // Axis spanning [10, 17] with unit spacing.
            auto engine = ResolutionConvolution::make({
                ResolutionConvolution::AxisGeometry{ 10.0, 1.0, N }
            });

            // Grid samples f(x) = 2 x - 5; multilinear interpolation is exact on affine data.
            std::vector<double> grid(N);
            for (std::size_t i = 0; i < N; ++i)
            {
                grid[i] = 2.0 * (10.0 + static_cast<double>(i)) - 5.0;
            }

            TEST_CHECK_NEARLY_EQUAL(engine->interpolate(grid, { 10.0 }), 15.0, 1.0e-12); // f(10)
            TEST_CHECK_NEARLY_EQUAL(engine->interpolate(grid, { 13.5 }), 22.0, 1.0e-12); // f(13.5)
            TEST_CHECK_NEARLY_EQUAL(engine->interpolate(grid, { 17.0 }), 29.0, 1.0e-12); // f(17), upper boundary

            // A point outside the grid must throw.
            TEST_CHECK_THROWS(InternalError, engine->interpolate(grid, { 9.0 }));
            TEST_CHECK_THROWS(InternalError, engine->interpolate(grid, { 17.5 }));
        }
} resolution_convolution_interpolate_test;

class ResolutionConvolution2DIdentityTest : public TestCase
{
    public:
        ResolutionConvolution2DIdentityTest() :
            TestCase("resolution_convolution_2d_identity_test")
        {
        }

        virtual void
        run() const override
        {
            const std::size_t N0 = 4;
            const std::size_t N1 = 8;

            auto engine = ResolutionConvolution::make({ unit_axis(N0), unit_axis(N1) });
            TEST_CHECK_EQUAL(engine->rank(), 2u);
            TEST_CHECK_EQUAL(engine->size(), N0 * N1);

            // 2D delta at the centre (row N0/2, column N1/2), row-major.
            std::vector<double> kernel(N0 * N1, 0.0);
            kernel[(N0 / 2) * N1 + (N1 / 2)] = 1.0;
            engine->set_resolution(kernel);

            std::vector<double> signal(N0 * N1);
            for (std::size_t i = 0; i < signal.size(); ++i)
            {
                signal[i] = static_cast<double>(i * i % 13);
            }

            const auto & convolved = engine->convolve(signal);
            for (std::size_t i = 0; i < signal.size(); ++i)
            {
                TEST_CHECK_NEARLY_EQUAL(convolved[i], signal[i], 1.0e-12);
            }
        }
} resolution_convolution_2d_identity_test;

class ResolutionConvolutionErrorsTest : public TestCase
{
    public:
        ResolutionConvolutionErrorsTest() :
            TestCase("resolution_convolution_errors_test")
        {
        }

        virtual void
        run() const override
        {
            // Unsupported dimensionality.
            TEST_CHECK_THROWS(InternalError, ResolutionConvolution::make({}));
            TEST_CHECK_THROWS(InternalError, ResolutionConvolution::make({ unit_axis(4), unit_axis(4), unit_axis(4), unit_axis(4), unit_axis(4) }));

            // Odd / too-small point counts.
            TEST_CHECK_THROWS(InternalError, ResolutionConvolution::make({ unit_axis(7) }));
            TEST_CHECK_THROWS(InternalError,
                              ResolutionConvolution::make({
                                  ResolutionConvolution::AxisGeometry{ 0.0, 1.0, 0u }
            }));

            // Non-positive spacing.
            TEST_CHECK_THROWS(InternalError,
                              ResolutionConvolution::make({
                                  ResolutionConvolution::AxisGeometry{ 0.0, 0.0, 4u }
            }));

            auto engine = ResolutionConvolution::make({ unit_axis(8) });

            // Kernel size mismatch.
            TEST_CHECK_THROWS(InternalError, engine->set_resolution(std::vector<double>(4, 1.0)));

            // Kernel with non-positive total weight.
            TEST_CHECK_THROWS(InternalError, engine->set_resolution(std::vector<double>(8, 0.0)));

            // convolve() before set_resolution().
            TEST_CHECK_THROWS(InternalError, engine->convolve(std::vector<double>(8, 1.0)));

            // Signal size mismatch after a valid resolution.
            engine->set_resolution(delta_kernel_1d(8));
            TEST_CHECK_THROWS(InternalError, engine->convolve(std::vector<double>(4, 1.0)));
        }
} resolution_convolution_errors_test;

// Unnormalised Gaussian density.
static double
gaussian_density(double x, double mu, double sigma)
{
    const double z = (x - mu) / sigma;
    return std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * M_PI));
}

class ResolutionConvolutionGaussianTest : public TestCase
{
    public:
        ResolutionConvolutionGaussianTest() :
            TestCase("resolution_convolution_gaussian_test")
        {
        }

        virtual void
        run() const override
        {
            const double mu1 = 0.0, s1 = 1.0, s2 = 0.5;
            const double sigma_total = std::sqrt(s1 * s1 + s2 * s2);

            // Pad by ~7 sigma_total on each side; N even.
            const double      spacing    = 0.06;
            const std::size_t N          = 256;
            const double      half_width = 0.5 * static_cast<double>(N) * spacing;
            const double      origin     = -half_width;

            auto engine = ResolutionConvolution::make({
                ResolutionConvolution::AxisGeometry{ origin, spacing, N }
            });

            std::vector<double> signal(N);
            for (std::size_t i = 0; i < N; ++i)
            {
                signal[i] = gaussian_density(origin + static_cast<double>(i) * spacing, mu1, s1);
            }

            std::vector<double> kernel(N);
            for (std::size_t i = 0; i < N; ++i)
            {
                const double offset = (static_cast<double>(i) - static_cast<double>(N / 2)) * spacing;
                kernel[i]           = gaussian_density(offset, 0.0, s2);
            }
            engine->set_resolution(kernel);

            const auto & convolved = engine->convolve(signal);

            // Reference: N(mu1, sigma_total), discretised and normalised to unit sum times the spacing.
            std::vector<double> reference(N);
            double              total = 0.0;
            for (std::size_t i = 0; i < N; ++i)
            {
                reference[i]  = gaussian_density(origin + static_cast<double>(i) * spacing, mu1, sigma_total);
                total        += reference[i];
            }
            total *= spacing;
            for (auto & r : reference)
            {
                r /= total;
            }

            // Compare well inside the padded region, at +-{0, 1, 2, 3} sigma_total.
            for (int k = -3; k <= 3; ++k)
            {
                const double x   = mu1 + static_cast<double>(k) * sigma_total;
                const auto   idx = static_cast<std::size_t>(std::llround((x - origin) / spacing));

                TEST_CHECK_RELATIVE_ERROR(convolved[idx], reference[idx], 1.0e-6);
            }
        }
} resolution_convolution_gaussian_test;

class ResolutionConvolutionBoxTest : public TestCase
{
    public:
        ResolutionConvolutionBoxTest() :
            TestCase("resolution_convolution_box_test")
        {
        }

        virtual void
        run() const override
        {
            // Box widths are odd multiples of the spacing, so edges fall on half-integer offsets.
            const double spacing = 0.1;
            const double w1      = 7.0 * spacing;
            const double w2      = 11.0 * spacing;
            const double a       = 0.5 * w1;
            const double b       = 0.5 * w2;

            const std::size_t N          = 64; // half-width 3.2, more than 3 box widths of padding
            const double      half_width = 0.5 * static_cast<double>(N) * spacing;
            const double      origin     = -half_width;

            auto engine = ResolutionConvolution::make({
                ResolutionConvolution::AxisGeometry{ origin, spacing, N }
            });

            std::vector<double> signal(N, 0.0);
            for (std::size_t i = 0; i < N; ++i)
            {
                const double x = origin + static_cast<double>(i) * spacing;
                if (std::abs(x) < a)
                {
                    signal[i] = 1.0 / w1;
                }
            }

            std::vector<double> kernel(N, 0.0);
            for (std::size_t i = 0; i < N; ++i)
            {
                const double offset = (static_cast<double>(i) - static_cast<double>(N / 2)) * spacing;
                if (std::abs(offset) < b)
                {
                    kernel[i] = 1.0 / w2;
                }
            }
            engine->set_resolution(kernel);

            const auto & convolved = engine->convolve(signal);

            // Trapezoid reference for the convolution of two boxes of half-widths a <= b.
            auto trapezoid = [a, b, w1, w2](double x) -> double
            {
                const double ax = std::abs(x);
                if (ax <= b - a)
                {
                    return 1.0 / w2;
                }
                if (ax <= a + b)
                {
                    return (a + b - ax) / (w1 * w2);
                }
                return 0.0;
            };

            for (double x : { 0.0, 0.1, -0.1, 0.2, -0.2, 0.5, -0.5, 0.9, -0.9, 1.2, -1.2 })
            {
                const auto   idx      = static_cast<std::size_t>(std::llround((x - origin) / spacing));
                const double expected = trapezoid(x);

                // Absolute comparison: expected rounds to a tiny nonzero at the a + b edge, where a relative check would spuriously fail.
                TEST_CHECK_NEARLY_EQUAL(convolved[idx], expected, 1.0e-6);
            }
        }
} resolution_convolution_box_test;

class ResolutionConvolutionInterpolatorAgreesTest : public TestCase
{
    public:
        ResolutionConvolutionInterpolatorAgreesTest() :
            TestCase("resolution_convolution_interpolator_agrees_test")
        {
        }

        virtual void
        run() const override
        {
            // 1D: reuse the affine grid f(x) = 2 x - 5 on [10, 17].
            {
                const std::size_t N      = 8;
                auto              engine = ResolutionConvolution::make({
                    ResolutionConvolution::AxisGeometry{ 10.0, 1.0, N }
                });

                std::vector<double> grid(N);
                for (std::size_t i = 0; i < N; ++i)
                {
                    grid[i] = 2.0 * (10.0 + static_cast<double>(i)) - 5.0;
                }

                // Points scattered across the grid, including the first and last cells.
                const std::vector<std::vector<double>> points{ { 10.0 }, { 10.3 }, { 13.5 }, { 16.7 }, { 17.0 } };

                auto interpolator = engine->make_interpolator(grid, points);
                TEST_CHECK_EQUAL(interpolator->size(), points.size());

                for (std::size_t i = 0; i < points.size(); ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL((*interpolator)(i), engine->interpolate(grid, points[i]), 1.0e-14);
                }
            }

            // 2D: an arbitrary grid, checking the interpolator against interpolate() point by point.
            {
                const std::size_t N0     = 4;
                const std::size_t N1     = 6;
                auto              engine = ResolutionConvolution::make({
                    ResolutionConvolution::AxisGeometry{ 0.0, 1.0, N0 },
                    ResolutionConvolution::AxisGeometry{ 0.0, 1.0, N1 }
                });

                std::vector<double> grid(N0 * N1);
                for (std::size_t i0 = 0; i0 < N0; ++i0)
                {
                    for (std::size_t i1 = 0; i1 < N1; ++i1)
                    {
                        grid[i0 * N1 + i1] = 3.0 * static_cast<double>(i0) - 2.0 * static_cast<double>(i1) + 1.0;
                    }
                }

                const std::vector<std::vector<double>> points{
                    { 0.0, 0.0 },
                    { 0.2, 0.4 },
                    { 1.5, 2.5 },
                    { 2.9, 4.9 },
                    { 3.0, 5.0 }
                };

                auto interpolator = engine->make_interpolator(grid, points);
                TEST_CHECK_EQUAL(interpolator->size(), points.size());

                for (std::size_t i = 0; i < points.size(); ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL((*interpolator)(i), engine->interpolate(grid, points[i]), 1.0e-14);
                }
            }
        }
} resolution_convolution_interpolator_agrees_test;

class ResolutionConvolutionInterpolatorErrorsTest : public TestCase
{
    public:
        ResolutionConvolutionInterpolatorErrorsTest() :
            TestCase("resolution_convolution_interpolator_errors_test")
        {
        }

        virtual void
        run() const override
        {
            const std::size_t N      = 8;
            auto              engine = ResolutionConvolution::make({
                ResolutionConvolution::AxisGeometry{ 10.0, 1.0, N }
            });

            std::vector<double> grid(N, 0.0);

            // Grid size mismatch.
            TEST_CHECK_THROWS(InternalError, engine->make_interpolator(std::vector<double>(4, 0.0), { { 10.0 } }));

            // Wrong coordinate count.
            TEST_CHECK_THROWS(InternalError,
                              engine->make_interpolator(grid,
                                                        {
                                                            { 10.0, 0.0 }
            }));

            // Point outside the grid.
            TEST_CHECK_THROWS(InternalError, engine->make_interpolator(grid, { { 9.0 } }));
            TEST_CHECK_THROWS(InternalError, engine->make_interpolator(grid, { { 17.5 } }));
        }
} resolution_convolution_interpolator_errors_test;
