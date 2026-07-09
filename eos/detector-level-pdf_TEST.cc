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
#include <eos/signal-pdf.hh>

#include <test/test.hh>

#include <cmath>
#include <vector>

using namespace test;
using namespace eos;

namespace
{
    // Grid used throughout: the sampling variable 'z' on a uniform grid over [-1.5, 5.5].
    static const std::size_t N       = 256;
    static const double      z_min   = -1.5;
    static const double      z_max   = +5.5;
    static const double      spacing = (z_max - z_min) / static_cast<double>(N - 1);

    double
    grid_coordinate(std::size_t i)
    {
        return z_min + static_cast<double>(i) * spacing;
    }
} // namespace

// A delta resolution (centred, i.e. peak at index N/2) leaves the truth PDF unchanged, so the
// detector-level PDF must reproduce the (clamped) truth PDF sampled on the grid.
class DetectorLevelPDFIdentityTest : public TestCase
{
    public:
        DetectorLevelPDFIdentityTest() :
            TestCase("detector_level_pdf_identity_test")
        {
        }

        virtual void
        run() const override
        {
            Parameters      p = Parameters::Defaults();
            ObservableCache cache(p);

            std::vector<double> resolution(N, 0.0);
            resolution[N / 2] = 1.0;

            SignalPDFPtr pdf = DetectorLevelPDF::make_1d(cache, "TestLegendre1D::P(z)", Options{}, DetectorLevelPDF::Axis{ "z", z_min, z_max, N }, resolution);

            // A standalone truth PDF for reference.
            Kinematics k;
            k.declare("z", 0.0);
            k.declare("z_min", z_min);
            k.declare("z_max", z_max);
            SignalPDFPtr truth = SignalPDF::make("TestLegendre1D::P(z)", p, k, Options{});

            for (std::size_t i : { std::size_t(80), std::size_t(128), std::size_t(170) })
            {
                const double z = grid_coordinate(i);

                pdf->kinematics().set("z", z);
                truth->kinematics().set("z", z);

                TEST_CHECK_NEARLY_EQUAL(pdf->evaluate_linear(), truth->evaluate_linear(), 1.0e-10);
                TEST_CHECK_NEARLY_EQUAL(pdf->evaluate(), std::log(truth->evaluate_linear()), 1.0e-10);
            }
        }
} detector_level_pdf_identity_test;

// The general case (resolution supplied as a SignalPDF) must agree with the pre-computed-grid case
// when the grid is obtained by sampling that same SignalPDF on the centred offset grid.
class DetectorLevelPDFGeneralVsGridTest : public TestCase
{
    public:
        DetectorLevelPDFGeneralVsGridTest() :
            TestCase("detector_level_pdf_general_vs_grid_test")
        {
        }

        virtual void
        run() const override
        {
            Parameters      p = Parameters::Defaults();
            ObservableCache cache(p);

            // Sample the resolution SignalPDF on the centred offset grid, exactly as the general case
            // does internally.
            Kinematics rk;
            rk.declare("z", 0.0);
            rk.declare("z_min", -static_cast<double>(N / 2) * spacing);
            rk.declare("z_max", +static_cast<double>(N / 2 - 1) * spacing);
            SignalPDFPtr resolution_pdf = SignalPDF::make("TestLegendre1D::P(z)", p, rk, Options{});

            std::vector<double> resolution(N, 0.0);
            for (std::size_t i = 0; i < N; ++i)
            {
                const double offset = (static_cast<double>(i) - static_cast<double>(N / 2)) * spacing;
                resolution_pdf->kinematics().set("z", offset);
                resolution[i] = resolution_pdf->evaluate_linear();
            }

            SignalPDFPtr from_grid = DetectorLevelPDF::make_1d(cache, "TestLegendre1D::P(z)", Options{}, DetectorLevelPDF::Axis{ "z", z_min, z_max, N }, resolution);

            SignalPDFPtr general = DetectorLevelPDF::make(cache,
                                                          "TestLegendre1D::P(z)",
                                                          "TestLegendre1D::P(z)",
                                                          Options{
            },
                                                          std::vector<DetectorLevelPDF::Axis>{ DetectorLevelPDF::Axis{ "z", z_min, z_max, N, "z" } });

            for (std::size_t i : { std::size_t(90), std::size_t(128), std::size_t(160), std::size_t(200) })
            {
                const double z = grid_coordinate(i);

                from_grid->kinematics().set("z", z);
                general->kinematics().set("z", z);

                TEST_CHECK_NEARLY_EQUAL(general->evaluate_linear(), from_grid->evaluate_linear(), 1.0e-10);
            }
        }
} detector_level_pdf_general_vs_grid_test;

// A cloned detector-level PDF must evaluate identically to its origin.
class DetectorLevelPDFCloneTest : public TestCase
{
    public:
        DetectorLevelPDFCloneTest() :
            TestCase("detector_level_pdf_clone_test")
        {
        }

        virtual void
        run() const override
        {
            Parameters      p = Parameters::Defaults();
            ObservableCache cache(p);

            std::vector<double> resolution(N, 0.0);
            resolution[N / 2 - 1] = 0.25;
            resolution[N / 2]     = 0.50;
            resolution[N / 2 + 1] = 0.25;

            SignalPDFPtr pdf    = DetectorLevelPDF::make_1d(cache, "TestLegendre1D::P(z)", Options{}, DetectorLevelPDF::Axis{ "z", z_min, z_max, N }, resolution);
            SignalPDFPtr cloned = std::static_pointer_cast<SignalPDF>(pdf->clone());

            for (std::size_t i : { std::size_t(100), std::size_t(140) })
            {
                const double z = grid_coordinate(i);

                pdf->kinematics().set("z", z);
                cloned->kinematics().set("z", z);

                TEST_CHECK_NEARLY_EQUAL(cloned->evaluate_linear(), pdf->evaluate_linear(), 1.0e-12);
            }
        }
} detector_level_pdf_clone_test;

class DetectorLevelPDFErrorsTest : public TestCase
{
    public:
        DetectorLevelPDFErrorsTest() :
            TestCase("detector_level_pdf_errors_test")
        {
        }

        virtual void
        run() const override
        {
            Parameters      p = Parameters::Defaults();
            ObservableCache cache(p);

            // Pre-computed grid of the wrong size.
            TEST_CHECK_THROWS(InternalError,
                              DetectorLevelPDF::make_1d(cache, "TestLegendre1D::P(z)", Options{}, DetectorLevelPDF::Axis{ "z", z_min, z_max, N }, std::vector<double>(N / 2, 1.0)));

            // Empty pre-computed grid.
            TEST_CHECK_THROWS(InternalError,
                              DetectorLevelPDF::make_1d(cache, "TestLegendre1D::P(z)", Options{}, DetectorLevelPDF::Axis{ "z", z_min, z_max, N }, std::vector<double>{}));

            // Unsupported dimensionality D = 5 (the general case dispatches through ResolutionConvolution).
            std::vector<DetectorLevelPDF::Axis> axes;
            for (std::size_t i = 0; i < 5; ++i)
            {
                axes.push_back(DetectorLevelPDF::Axis{ "z", z_min, z_max, 4 });
            }
            TEST_CHECK_THROWS(InternalError, DetectorLevelPDF::make(cache, "TestLegendre1D::P(z)", "TestLegendre1D::P(z)", Options{}, axes));

            // Empty range.
            TEST_CHECK_THROWS(InternalError,
                              DetectorLevelPDF::make_1d(cache, "TestLegendre1D::P(z)", Options{}, DetectorLevelPDF::Axis{ "z", 1.0, 1.0, N }, std::vector<double>(N, 1.0)));
        }
} detector_level_pdf_errors_test;
