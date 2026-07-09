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

#include <test/test.hh>
#include <eos/maths/resolution-convolution.hh>
#include <eos/utils/exception.hh>

#include <vector>

using namespace test;
using namespace eos;

// A grid axis with `points` points, spanning [0, points) with unit spacing (origin 0).
static ResolutionConvolution::AxisGeometry unit_axis(std::size_t points)
{
    return ResolutionConvolution::AxisGeometry{ 0.0, 1.0, points };
}

// A discrete delta kernel: 1 at the centre (zero-offset) index, 0 elsewhere. In centred order the
// centre sits at index points / 2 along each axis. Convolving with it is the identity.
static std::vector<double> delta_kernel_1d(std::size_t points)
{
    std::vector<double> kernel(points, 0.0);
    kernel[points / 2] = 1.0;
    return kernel;
}

class ResolutionConvolutionIdentityTest :
    public TestCase
{
    public:
        ResolutionConvolutionIdentityTest() :
            TestCase("resolution_convolution_identity_test")
        {
        }

        virtual void run() const override
        {
            const std::size_t N = 8;

            auto engine = ResolutionConvolution::make({ unit_axis(N) });
            TEST_CHECK_EQUAL(engine->rank(), 1u);
            TEST_CHECK_EQUAL(engine->size(), N);

            engine->set_resolution(delta_kernel_1d(N));

            // Convolution with a delta at zero offset returns the signal unchanged.
            std::vector<double> signal{ 0.0, 1.0, 3.0, 7.0, 2.0, 0.5, 0.25, 0.0 };
            const auto & convolved = engine->convolve(signal);

            for (std::size_t i = 0 ; i < N ; ++i)
                TEST_CHECK_NEARLY_EQUAL(convolved[i], signal[i], 1.0e-12);
        }
} resolution_convolution_identity_test;

class ResolutionConvolutionNormalisationTest :
    public TestCase
{
    public:
        ResolutionConvolutionNormalisationTest() :
            TestCase("resolution_convolution_normalisation_test")
        {
        }

        virtual void run() const override
        {
            const std::size_t N = 8;

            auto engine = ResolutionConvolution::make({ unit_axis(N) });

            // A uniform, deliberately un-normalised kernel. The engine renormalises it to unit sum,
            // so convolving with it must map every grid point to the mean of the signal.
            std::vector<double> kernel(N, 2.0);
            engine->set_resolution(kernel);

            std::vector<double> signal{ 0.0, 1.0, 3.0, 7.0, 2.0, 0.5, 0.25, 0.25 };
            double mean = 0.0;
            for (const auto & s : signal)
                mean += s;
            mean /= static_cast<double>(N);

            const auto & convolved = engine->convolve(signal);
            for (std::size_t i = 0 ; i < N ; ++i)
                TEST_CHECK_NEARLY_EQUAL(convolved[i], mean, 1.0e-12);
        }
} resolution_convolution_normalisation_test;

class ResolutionConvolutionConstantSignalTest :
    public TestCase
{
    public:
        ResolutionConvolutionConstantSignalTest() :
            TestCase("resolution_convolution_constant_signal_test")
        {
        }

        virtual void run() const override
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
            const auto & convolved = engine->convolve(signal);
            for (std::size_t i = 0 ; i < N ; ++i)
                TEST_CHECK_NEARLY_EQUAL(convolved[i], 3.5, 1.0e-12);
        }
} resolution_convolution_constant_signal_test;

class ResolutionConvolutionShiftTest :
    public TestCase
{
    public:
        ResolutionConvolutionShiftTest() :
            TestCase("resolution_convolution_shift_test")
        {
        }

        virtual void run() const override
        {
            const std::size_t N = 8;

            auto engine = ResolutionConvolution::make({ unit_axis(N) });

            // A delta displaced by +1 from the centre. Convolving shifts the signal by +1 (with the
            // circular wrap-around), which pins down the centred -> wrap-around ordering convention.
            std::vector<double> kernel(N, 0.0);
            kernel[N / 2 + 1] = 1.0;
            engine->set_resolution(kernel);

            std::vector<double> signal{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
            const auto & convolved = engine->convolve(signal);

            for (std::size_t i = 0 ; i < N ; ++i)
                TEST_CHECK_NEARLY_EQUAL(convolved[i], signal[(i + N - 1) % N], 1.0e-12);
        }
} resolution_convolution_shift_test;

class ResolutionConvolutionInterpolateTest :
    public TestCase
{
    public:
        ResolutionConvolutionInterpolateTest() :
            TestCase("resolution_convolution_interpolate_test")
        {
        }

        virtual void run() const override
        {
            const std::size_t N = 8;

            // Axis spanning [10, 17] with unit spacing.
            auto engine = ResolutionConvolution::make({ ResolutionConvolution::AxisGeometry{ 10.0, 1.0, N } });

            // A grid sampling the affine function f(x) = 2 x - 5. Multilinear interpolation is exact
            // on affine data.
            std::vector<double> grid(N);
            for (std::size_t i = 0 ; i < N ; ++i)
                grid[i] = 2.0 * (10.0 + static_cast<double>(i)) - 5.0;

            TEST_CHECK_NEARLY_EQUAL(engine->interpolate(grid, { 10.0  }),  15.0, 1.0e-12); // f(10)
            TEST_CHECK_NEARLY_EQUAL(engine->interpolate(grid, { 13.5  }),  22.0, 1.0e-12); // f(13.5)
            TEST_CHECK_NEARLY_EQUAL(engine->interpolate(grid, { 17.0  }),  29.0, 1.0e-12); // f(17), upper boundary

            // A point outside the grid must throw.
            TEST_CHECK_THROWS(InternalError, engine->interpolate(grid, {  9.0 }));
            TEST_CHECK_THROWS(InternalError, engine->interpolate(grid, { 17.5 }));
        }
} resolution_convolution_interpolate_test;

class ResolutionConvolution2DIdentityTest :
    public TestCase
{
    public:
        ResolutionConvolution2DIdentityTest() :
            TestCase("resolution_convolution_2d_identity_test")
        {
        }

        virtual void run() const override
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
            for (std::size_t i = 0 ; i < signal.size() ; ++i)
                signal[i] = static_cast<double>(i * i % 13);

            const auto & convolved = engine->convolve(signal);
            for (std::size_t i = 0 ; i < signal.size() ; ++i)
                TEST_CHECK_NEARLY_EQUAL(convolved[i], signal[i], 1.0e-12);
        }
} resolution_convolution_2d_identity_test;

class ResolutionConvolutionErrorsTest :
    public TestCase
{
    public:
        ResolutionConvolutionErrorsTest() :
            TestCase("resolution_convolution_errors_test")
        {
        }

        virtual void run() const override
        {
            // Unsupported dimensionality.
            TEST_CHECK_THROWS(InternalError, ResolutionConvolution::make({}));
            TEST_CHECK_THROWS(InternalError, ResolutionConvolution::make(
                { unit_axis(4), unit_axis(4), unit_axis(4), unit_axis(4), unit_axis(4) }));

            // Odd / too-small point counts.
            TEST_CHECK_THROWS(InternalError, ResolutionConvolution::make({ unit_axis(7) }));
            TEST_CHECK_THROWS(InternalError, ResolutionConvolution::make(
                { ResolutionConvolution::AxisGeometry{ 0.0, 1.0, 0u } }));

            // Non-positive spacing.
            TEST_CHECK_THROWS(InternalError, ResolutionConvolution::make(
                { ResolutionConvolution::AxisGeometry{ 0.0, 0.0, 4u } }));

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
