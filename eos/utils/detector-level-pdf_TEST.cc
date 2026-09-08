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

#include <eos/maths/power-of.hh>
#include <eos/observable.hh>
#include <eos/signal-pdf.hh>
#include <eos/utils/detector-level-pdf.hh>

#include <test/test-pdfs.hh>
#include <test/test.hh>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

using namespace test;
using namespace eos;

namespace
{
    // Grid used throughout the legacy tests: 'z' on a uniform grid over [-1.0, 5.0], chosen so that
    // TestLegendre1D::P(z) = z (4 - z)'s normalization (its integral over the range) stays positive.
    static const std::size_t N       = 256;
    static const double      z_min   = -1.0;
    static const double      z_max   = +5.0;
    static const double      spacing = (z_max - z_min) / static_cast<double>(N - 1);

    double
    grid_coordinate(std::size_t i)
    {
        return z_min + static_cast<double>(i) * spacing;
    }

    // Analytic references, deliberately independent of the closed forms the test PDFs use, so that
    // the tests do not verify the implementation against itself.
    double
    gaussian_norm(double mu, double sigma, double lo, double hi)
    {
        const double s = sigma * M_SQRT2;

        return sigma * std::sqrt(M_PI / 2.0) * (std::erf((hi - mu) / s) - std::erf((lo - mu) / s));
    }

    double
    gaussian_density(double x, double mu, double sigma)
    {
        return std::exp(-(x - mu) * (x - mu) / (2.0 * sigma * sigma)) / (sigma * std::sqrt(2.0 * M_PI));
    }

    double
    box_norm(double centre, double width, double lo, double hi)
    {
        const double a = std::max(lo, centre - 0.5 * width);
        const double b = std::min(hi, centre + 0.5 * width);

        return std::max(0.0, b - a);
    }

    // box(w1) (x) box(w2), both of unit height; a trapezoid, a triangle when w1 == w2.
    double
    box_convolution(double t, double w1, double w2)
    {
        const double a  = std::min(w1, w2);
        const double b  = std::max(w1, w2);
        const double at = std::abs(t);

        if (at <= 0.5 * (b - a))
        {
            return a;
        }

        if (at < 0.5 * (a + b))
        {
            return 0.5 * (a + b) - at;
        }

        return 0.0;
    }
} // namespace

// A delta resolution leaves the truth PDF's shape unchanged, so the result must equal the
// (clamped) truth PDF, normalized by the truth's own normalization.
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

            cache.update();

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

                const double reference = truth->evaluate_linear();

                TEST_CHECK_NEARLY_EQUAL(pdf->evaluate_linear(), reference, 1.0e-10);
                TEST_CHECK_NEARLY_EQUAL(pdf->evaluate(), std::log(reference), 1.0e-10);
                TEST_CHECK_NEARLY_EQUAL(pdf->normalization(), truth->normalization(), 1.0e-10);
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

            cache.update();

            for (std::size_t i : { std::size_t(90), std::size_t(128), std::size_t(160), std::size_t(200) })
            {
                const double z = grid_coordinate(i);

                from_grid->kinematics().set("z", z);
                general->kinematics().set("z", z);

                TEST_CHECK_NEARLY_EQUAL(general->evaluate_linear(), from_grid->evaluate_linear(), 1.0e-10);
            }
        }
} detector_level_pdf_general_vs_grid_test;

// A cloned detector-level PDF must evaluate identically to its origin. The clone is tied to its own,
// independent ObservableCache, which must be updated separately.
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

            cache.update();

            ObservableCache cloned_cache = std::dynamic_pointer_cast<DetectorLevelPDF>(cloned)->cache();
            cloned_cache.update();

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

            // A detector-level PDF is never backed by a single observable: it must never be able to
            // reach the ObservableCache as an evaluated observable.
            std::vector<double> resolution(N, 0.0);
            resolution[N / 2] = 1.0;
            SignalPDFPtr pdf  = DetectorLevelPDF::make_1d(cache, "TestLegendre1D::P(z)", Options{}, DetectorLevelPDF::Axis{ "z", z_min, z_max, N }, resolution);
            TEST_CHECK_THROWS(InternalError, pdf->unnormalized_pdf());
        }
} detector_level_pdf_errors_test;

// update_grid() must be a no-op until the cache's generation advances, then recompute in place
// (the buffer address never changes); cache() must identify the cache passed in at construction.
class DetectorLevelPDFGridAccessTest : public TestCase
{
    public:
        DetectorLevelPDFGridAccessTest() :
            TestCase("detector_level_pdf_grid_access_test")
        {
        }

        virtual void
        run() const override
        {
            Parameters      p = Parameters::Defaults();
            ObservableCache cache(p);

            static const std::size_t M  = 64;
            static const double      lo = -20.0;
            static const double      hi = +20.0;

            std::vector<double> resolution(M, 0.0);
            resolution[M / 2] = 1.0;

            SignalPDFPtr base = DetectorLevelPDF::make_1d(cache, "TestGaussian1D::P(x)", Options{}, DetectorLevelPDF::Axis{ "x", lo, hi, M }, resolution);
            auto         pdf  = std::dynamic_pointer_cast<DetectorLevelPDF>(base);
            TEST_CHECK(pdf != nullptr);

            TEST_CHECK(pdf->cache() == cache);

            cache.update();
            pdf->update_grid();

            const double         value_before   = pdf->grid()[M / 2];
            const double * const address_before = pdf->grid().data();

            // Moving a parameter without an intervening cache.update() must leave the grid untouched.
            p["TestGaussian1D::mu"] = 3.0;
            pdf->update_grid();
            TEST_CHECK_EQUAL(pdf->grid()[M / 2], value_before);
            TEST_CHECK(pdf->grid().data() == address_before);

            // Once the cache is updated, update_grid() must pick up the new parameter value, without
            // reallocating the grid buffer.
            cache.update();
            pdf->update_grid();
            TEST_CHECK(pdf->grid()[M / 2] != value_before);
            TEST_CHECK(pdf->grid().data() == address_before);
        }
} detector_level_pdf_grid_access_test;

// A Gaussian signal convolved with a Gaussian resolution: the result is N(mu, sqrt(s1^2 + s2^2)).
class DetectorLevelPDFGaussianTest : public TestCase
{
    public:
        DetectorLevelPDFGaussianTest() :
            TestCase("detector_level_pdf_gaussian_test")
        {
        }

        virtual void
        run() const override
        {
            Parameters      p = Parameters::Defaults();
            ObservableCache cache(p);

            static const std::size_t M       = 512;
            static const double      x_min   = -15.0;
            static const double      x_max   = +15.0;
            static const double      spacing = (x_max - x_min) / static_cast<double>(M - 1);

            SignalPDFPtr pdf = DetectorLevelPDF::make(cache,
                                                      "TestGaussian1D::P(x)",
                                                      "TestGaussianResolution1D::P(x)",
                                                      Options{
            },
                                                      std::vector<DetectorLevelPDF::Axis>{ DetectorLevelPDF::Axis{ "x", x_min, x_max, M } });

            cache.update();

            const double mu      = p["TestGaussian1D::mu"].evaluate();
            const double sigma_t = p["TestGaussian1D::sigma"].evaluate();
            const double sigma_r = p["TestGaussianResolution1D::sigma"].evaluate();
            const double sigma_c = std::sqrt(sigma_t * sigma_t + sigma_r * sigma_r);

            static const double tolerance = 1.0e-6;

            for (std::size_t i : { std::size_t(176),
                                   std::size_t(196),
                                   std::size_t(216),
                                   std::size_t(236),
                                   std::size_t(246),
                                   std::size_t(266),
                                   std::size_t(276),
                                   std::size_t(296),
                                   std::size_t(316),
                                   std::size_t(336) })
            {
                const double x = x_min + static_cast<double>(i) * spacing;

                pdf->kinematics().set("x", x);

                TEST_CHECK_RELATIVE_ERROR(std::exp(pdf->evaluate() - pdf->normalization()), gaussian_density(x, mu, sigma_c), tolerance);
            }

            TEST_CHECK_RELATIVE_ERROR(std::exp(pdf->normalization()), gaussian_norm(mu, sigma_t, x_min, x_max), tolerance);
        }
} detector_level_pdf_gaussian_test;

// A box-shaped signal convolved with a box-shaped resolution (equal widths give a triangle): an
// exact, compactly supported reference, supplied here as a pre-computed grid rather than a SignalPDF.
class DetectorLevelPDFBoxTest : public TestCase
{
    public:
        DetectorLevelPDFBoxTest() :
            TestCase("detector_level_pdf_box_test")
        {
        }

        virtual void
        run() const override
        {
            Parameters p = Parameters::Defaults();

            // The resolution's offset grid has a node at offset 0, so an even width would land its
            // edges exactly on a node; override the default (even) width with an odd one to avoid that.
            p["TestBoxResolution1D::width"] = 1.9;

            ObservableCache cache(p);

            // The truth's own grid is shifted by half a spacing (x = 0 sits exactly halfway between
            // two nodes), so an *even* multiple of the spacing keeps its edges off its grid nodes.
            static const std::size_t M       = 256;
            static const double      spacing = 0.1;
            static const double      x_min   = -(static_cast<double>(M / 2) - 0.5) * spacing;
            static const double      x_max   = x_min + static_cast<double>(M - 1) * spacing;

            const double centre    = p["TestBox1D::centre"].evaluate();
            const double width     = p["TestBox1D::width"].evaluate();
            const double width_res = p["TestBoxResolution1D::width"].evaluate();

            std::vector<double> resolution(M, 0.0);
            for (std::size_t i = 0; i < M; ++i)
            {
                const double offset = (static_cast<double>(i) - static_cast<double>(M / 2)) * spacing;
                resolution[i]       = (std::abs(offset) <= 0.5 * width_res) ? 1.0 : 0.0;
            }

            SignalPDFPtr pdf = DetectorLevelPDF::make_from_grid(cache,
                                                                "TestBox1D::P(x)",
                                                                Options{
            },
                                                                std::vector<DetectorLevelPDF::Axis>{ DetectorLevelPDF::Axis{ "x", x_min, x_max, M } },
                                                                resolution);

            cache.update();

            static const double tolerance = 1.0e-6;

            for (std::size_t i : { std::size_t(110),
                                   std::size_t(115),
                                   std::size_t(120),
                                   std::size_t(124),
                                   std::size_t(126),
                                   std::size_t(128),
                                   std::size_t(130),
                                   std::size_t(132),
                                   std::size_t(136),
                                   std::size_t(142) })
            {
                const double x = x_min + static_cast<double>(i) * spacing;

                pdf->kinematics().set("x", x);

                const double reference = box_convolution(x - centre, width, width_res) / (width * width_res);

                TEST_CHECK_RELATIVE_ERROR(std::exp(pdf->evaluate() - pdf->normalization()), reference, tolerance);
            }

            TEST_CHECK_RELATIVE_ERROR(std::exp(pdf->normalization()), box_norm(centre, width, x_min, x_max), tolerance);
        }
} detector_level_pdf_box_test;

// Separable 2D case: a Gaussian signal convolved with a Gaussian resolution along x and y, built
// both from the general (SignalPDF) resolution and an equivalent pre-computed grid, which must agree.
class DetectorLevelPDF2DTest : public TestCase
{
    public:
        DetectorLevelPDF2DTest() :
            TestCase("detector_level_pdf_2d_test")
        {
        }

        virtual void
        run() const override
        {
            Parameters      p = Parameters::Defaults();
            ObservableCache cache(p);

            // Both axes share the same range and point count; y_min/y_max are aliases of x_min/x_max
            // (not independent values) so that using the wrong one for either axis below is a no-op.
            static const std::size_t M       = 128;
            static const double      x_min   = -12.8;
            static const double      x_max   = +12.8;
            static const double      y_min   = x_min;
            static const double      y_max   = x_max;
            static const double      spacing = (x_max - x_min) / static_cast<double>(M - 1);

            const std::vector<DetectorLevelPDF::Axis> axes{
                DetectorLevelPDF::Axis{ "x", x_min, x_max, M },
                DetectorLevelPDF::Axis{ "y", y_min, y_max, M }
            };

            SignalPDFPtr general = DetectorLevelPDF::make(cache, "TestGaussian2D::P(x,y)", "TestGaussianResolution2D::P(x,y)", Options{}, axes);

            // Sample the same resolution SignalPDF onto the centred offset grid, exactly as the
            // general case does internally, to build an equivalent pre-computed-grid PDF.
            Kinematics rk;
            rk.declare("x", 0.0);
            rk.declare("x_min", -static_cast<double>(M / 2) * spacing);
            rk.declare("x_max", +static_cast<double>(M / 2 - 1) * spacing);
            rk.declare("y", 0.0);
            rk.declare("y_min", -static_cast<double>(M / 2) * spacing);
            rk.declare("y_max", +static_cast<double>(M / 2 - 1) * spacing);
            SignalPDFPtr resolution_pdf = SignalPDF::make("TestGaussianResolution2D::P(x,y)", p, rk, Options{});

            std::vector<double> resolution(M * M, 0.0);
            for (std::size_t ix = 0; ix < M; ++ix)
            {
                const double ox = (static_cast<double>(ix) - static_cast<double>(M / 2)) * spacing;
                for (std::size_t iy = 0; iy < M; ++iy)
                {
                    const double oy = (static_cast<double>(iy) - static_cast<double>(M / 2)) * spacing;

                    resolution_pdf->kinematics().set("x", ox);
                    resolution_pdf->kinematics().set("y", oy);
                    resolution[ix * M + iy] = resolution_pdf->evaluate_linear();
                }
            }

            SignalPDFPtr from_grid = DetectorLevelPDF::make_from_grid(cache, "TestGaussian2D::P(x,y)", Options{}, axes, resolution);

            cache.update();

            const double mu_x     = p["TestGaussian2D::mu_x"].evaluate();
            const double sigma_x  = p["TestGaussian2D::sigma_x"].evaluate();
            const double mu_y     = p["TestGaussian2D::mu_y"].evaluate();
            const double sigma_y  = p["TestGaussian2D::sigma_y"].evaluate();
            const double sigma_rx = p["TestGaussianResolution2D::sigma_x"].evaluate();
            const double sigma_ry = p["TestGaussianResolution2D::sigma_y"].evaluate();

            const double sigma_cx = std::sqrt(sigma_x * sigma_x + sigma_rx * sigma_rx);
            const double sigma_cy = std::sqrt(sigma_y * sigma_y + sigma_ry * sigma_ry);

            static const double tolerance = 1.0e-5;

            for (std::size_t ix : { std::size_t(56), std::size_t(64), std::size_t(72) })
            {
                for (std::size_t iy : { std::size_t(56), std::size_t(64), std::size_t(72) })
                {
                    const double x = x_min + static_cast<double>(ix) * spacing;
                    const double y = y_min + static_cast<double>(iy) * spacing;

                    general->kinematics().set("x", x);
                    general->kinematics().set("y", y);
                    from_grid->kinematics().set("x", x);
                    from_grid->kinematics().set("y", y);

                    const double reference = gaussian_density(x, mu_x, sigma_cx) * gaussian_density(y, mu_y, sigma_cy);

                    TEST_CHECK_RELATIVE_ERROR(std::exp(general->evaluate() - general->normalization()), reference, tolerance);
                    TEST_CHECK_RELATIVE_ERROR(std::exp(from_grid->evaluate() - from_grid->normalization()), reference, tolerance);
                }
            }

            const double reference_norm = gaussian_norm(mu_x, sigma_x, x_min, x_max) * gaussian_norm(mu_y, sigma_y, y_min, y_max);

            TEST_CHECK_RELATIVE_ERROR(std::exp(general->normalization()), reference_norm, tolerance);
            TEST_CHECK_RELATIVE_ERROR(std::exp(from_grid->normalization()), reference_norm, tolerance);
        }
} detector_level_pdf_2d_test;

// A truth PDF whose support is a strict sub-interval of the padded grid. Its normalization is
// asked for over the grid range, and since the PDF vanishes identically outside its support that
// must give the same number as the integral over the support alone. B->Dlnu::NormalizationPDF is
// ~5.7e-12, below the integrator's absolute tolerance of 1e-9, so hcubature accepts its first
// unrefined rule and misses the true integral by ~5% once the range extends past the support.
class DetectorLevelPDFPaddedSupportTest : public TestCase
{
    public:
        DetectorLevelPDFPaddedSupportTest() :
            TestCase("detector_level_pdf_padded_support_test")
        {
        }

        virtual void
        run() const override
        {
            Parameters      p = Parameters::Defaults();
            ObservableCache cache(p);

            const Options options{
                { "l"_ok, "mu"_ov },
                { "q"_ok,  "d"_ov }
            };

            const double q2_min = power_of<2>(p["mass::mu"].evaluate());
            const double q2_max = power_of<2>(p["mass::B_d"].evaluate() - p["mass::D_d"].evaluate());

            // The grid is padded by ~2 GeV^2 on either side of the physical support.
            static const std::size_t M        = 256;
            static const double      grid_min = -2.0;
            static const double      grid_max = 13.6;

            std::vector<double> resolution(M, 0.0);
            resolution[M / 2] = 1.0;

            SignalPDFPtr pdf = DetectorLevelPDF::make_1d(cache, "B->Dlnu::P(q2)", options, DetectorLevelPDF::Axis{ "q2", grid_min, grid_max, M }, resolution);

            cache.update();

            // Reference: the same normalization observable summed over 64 panels covering the
            // physical support, so that each panel gets its own cubature rule.
            Kinematics k;
            k.declare("q2_min", q2_min);
            k.declare("q2_max", q2_max);
            ObservablePtr normalization = Observable::make("B->Dlnu::NormalizationPDF(q2)", p, k, options);

            static const std::size_t panels = 64;

            double reference_norm = 0.0;
            for (std::size_t i = 0; i < panels; ++i)
            {
                const double lo = q2_min + static_cast<double>(i) * (q2_max - q2_min) / static_cast<double>(panels);
                const double hi = q2_min + static_cast<double>(i + 1) * (q2_max - q2_min) / static_cast<double>(panels);

                k["q2_min"] = lo;
                k["q2_max"] = hi;

                reference_norm += normalization->evaluate();
            }

            static const double tolerance = 1.0e-3;

            TEST_CHECK_RELATIVE_ERROR(std::exp(pdf->normalization()), reference_norm, tolerance);

            // The same excess in the normalized density: with a delta resolution it is the truth
            // PDF's own normalized density, whose integral over the support is 1 by construction.
            double integral = 0.0;
            for (std::size_t i = 0; i <= 2 * panels; ++i)
            {
                const double q2     = q2_min + static_cast<double>(i) * (q2_max - q2_min) / static_cast<double>(2 * panels);
                const double weight = ((0 == i) || (2 * panels == i)) ? 0.5 : 1.0;

                pdf->kinematics().set("q2", q2);

                integral += weight * std::exp(pdf->evaluate() - pdf->normalization());
            }
            integral *= (q2_max - q2_min) / static_cast<double>(2 * panels);

            TEST_CHECK_RELATIVE_ERROR(integral, 1.0, 1.0e-2);
        }
} detector_level_pdf_padded_support_test;
