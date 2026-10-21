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

#include <eos/maths/faddeeva.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/observable.hh>
#include <eos/signal-pdf.hh>
#include <eos/utils/detector-level-pdf.hh>

#include <test/test-pdfs.hh>
#include <test/test.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <functional>
#include <memory>
#include <set>
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

    // The integral of 1 + Re[c C_{R_G}(x - z)] over [a, b] for a Gaussian resolution, by Simpson's rule.
    double
    smeared_pole_integral(double a, double b, const complex<double> & c, const complex<double> & z, double sigma)
    {
        static const std::size_t n = 20000;
        const double             h = (b - a) / static_cast<double>(n);

        double result = 0.0;
        for (std::size_t k = 0; k <= n; ++k)
        {
            const double          x       = a + static_cast<double>(k) * h;
            const complex<double> C       = -1.0i * std::sqrt(M_PI / 2.0) / sigma * faddeeva_w((x - z) / (std::sqrt(2.0) * sigma));
            const double          weight  = ((0 == k) || (n == k)) ? 1.0 : ((k % 2) ? 4.0 : 2.0);
            result                       += weight * (1.0 + (c * C).real());
        }

        return result * h / 3.0;
    }

    // The Cauchy transform int du R(u) / (zeta - u), Im zeta > 0, of the kernel R normalized over [-reach, reach], by adaptive quadrature.
    complex<double>
    brute_force_cauchy_transform(const std::function<double(const double &)> & R, const complex<double> & zeta, const double & reach, std::vector<double> points)
    {
        const double s = zeta.real(), g = zeta.imag();
        for (const double k : { -100.0, -10.0, 0.0, 10.0, 100.0 })
        {
            points.push_back(s + k * g);
        }
        points.push_back(-reach);
        points.push_back(+reach);
        std::sort(points.begin(), points.end());

        const auto config = GSL::QAGS::Config().epsrel(1.0e-12).epsabs(0.0);
        double     norm = 0.0, re = 0.0, im = 0.0;
        for (std::size_t i = 1; i < points.size(); ++i)
        {
            norm += integrate<GSL::QAGS>(R, points[i - 1], points[i], config);
            re   += integrate<GSL::QAGS>([&](const double & u) { return R(u) * (s - u) / ((s - u) * (s - u) + g * g); }, points[i - 1], points[i], config);
            im   += integrate<GSL::QAGS>([&](const double & u) { return -R(u) * g / ((s - u) * (s - u) + g * g); }, points[i - 1], points[i], config);
        }

        return complex<double>(re, im) / norm;
    }

    // Kernel::CrystalBall(u, mu, sigma, alpha_L, n_L, alpha_R, n_R), up to its normalization.
    double
    crystal_ball_shape(const double & u, const double & mu, const double & sigma, const double & alpha_L, const double & n_L, const double & alpha_R, const double & n_R)
    {
        const double t = (u - mu) / sigma;
        if (t < -alpha_L)
        {
            return std::exp(-0.5 * alpha_L * alpha_L) * std::pow(1.0 + alpha_L * (-t - alpha_L) / n_L, -n_L);
        }
        if (t > alpha_R)
        {
            return std::exp(-0.5 * alpha_R * alpha_R) * std::pow(1.0 + alpha_R * (t - alpha_R) / n_R, -n_R);
        }

        return std::exp(-0.5 * t * t);
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

// A detector-level PDF cloned into another cache must evaluate identically to its origin once that
// cache is updated; the clones that would bind to a private cache must throw.
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

            SignalPDFPtr pdf = DetectorLevelPDF::make_1d(cache, "TestLegendre1D::P(z)", Options{}, DetectorLevelPDF::Axis{ "z", z_min, z_max, N }, resolution);

            ObservableCache cloned_cache(p.clone());
            SignalPDFPtr    cloned = std::dynamic_pointer_cast<DetectorLevelPDF>(pdf)->clone(cloned_cache);

            cache.update();
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

            // A failed construction must register neither the batch nor the normalization.
            {
                ObservableCache                           pristine(p);
                const std::vector<DetectorLevelPDF::Axis> axis{
                    DetectorLevelPDF::Axis{ "z", z_min, z_max, N }
                };
                TEST_CHECK_THROWS(InternalError, DetectorLevelPDF::make_from_grid(pristine, "TestLegendre1D::P(z)", Options{}, axis, std::vector<double>(N / 2, 1.0)));
                TEST_CHECK_THROWS(InternalError, DetectorLevelPDF::make_from_grid(pristine, "TestLegendre1D::P(z)", Options{}, axis, std::vector<double>{}));
                TEST_CHECK_THROWS(UnknownObservableError,
                                  DetectorLevelPDF::make_from_expression(pristine, "TestLegendre1D::P(z)", Options{}, axis, "<<TestLegendre1D::NOSUCH(z)>>"));
                TEST_CHECK_EQUAL(pristine.size(), 0u);
                TEST_CHECK_EQUAL(pristine.add_batch({}).value(), 0u);
            }

            // Pre-computed grid of the wrong size.
            TEST_CHECK_THROWS(InternalError,
                              DetectorLevelPDF::make_1d(cache, "TestLegendre1D::P(z)", Options{}, DetectorLevelPDF::Axis{ "z", z_min, z_max, N }, std::vector<double>(N / 2, 1.0)));

            // Empty pre-computed grid.
            TEST_CHECK_THROWS(InternalError,
                              DetectorLevelPDF::make_1d(cache, "TestLegendre1D::P(z)", Options{}, DetectorLevelPDF::Axis{ "z", z_min, z_max, N }, std::vector<double>{}));

            // Unsupported dimensionality D = 5 (dispatched through ResolutionConvolution).
            std::vector<DetectorLevelPDF::Axis> axes;
            for (std::size_t i = 0; i < 5; ++i)
            {
                axes.push_back(DetectorLevelPDF::Axis{ "z" + std::to_string(i), z_min, z_max, 4 });
            }
            TEST_CHECK_THROWS(InternalError, DetectorLevelPDF::make_from_expression(cache, "TestLegendre1D::P(z)", Options{}, axes, "1.0"));

            // A variable that names two axes, and two axes that share an offset variable.
            {
                const std::vector<DetectorLevelPDF::Axis> twice{
                    DetectorLevelPDF::Axis{ "z", z_min, z_max, 4 },
                    DetectorLevelPDF::Axis{ "z", z_min, z_max, 4 }
                };
                TEST_CHECK_THROWS(InternalError, DetectorLevelPDF::make_from_expression(cache, "TestLegendre1D::P(z)", Options{}, twice, "1.0"));

                const std::vector<DetectorLevelPDF::Axis> shared{
                    DetectorLevelPDF::Axis{ "z", z_min, z_max, 4, "u" },
                    DetectorLevelPDF::Axis{ "y", z_min, z_max, 4, "u" }
                };
                TEST_CHECK_THROWS(InternalError, DetectorLevelPDF::make_from_expression(cache, "TestLegendre1D::P(z)", Options{}, shared, "1.0"));
            }

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

            SignalPDFPtr pdf = DetectorLevelPDF::make_from_expression(cache,
                                                                      "TestGaussian1D::P(x)",
                                                                      Options{
            },
                                                                      std::vector<DetectorLevelPDF::Axis>{ DetectorLevelPDF::Axis{ "x", x_min, x_max, M } },
                                                                      "<<TestGaussianResolution1D::UnnormalizedPDF(x)>>");

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

            // the cache reports the resolution's parameters
            const std::set<Parameter::Id> used(cache.begin_used_parameter_ids(), cache.end_used_parameter_ids());
            TEST_CHECK(used.contains(p["TestGaussianResolution1D::sigma"].id()));
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
// both from a resolution expression and an equivalent pre-computed grid, which must agree.
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

            SignalPDFPtr from_expression =
                    DetectorLevelPDF::make_from_expression(cache, "TestGaussian2D::P(x,y)", Options{}, axes, "<<TestGaussianResolution2D::UnnormalizedPDF(x,y)>>");

            // Sample the same resolution density onto the centred offset grid, as the expression is
            // sampled internally, to build an equivalent pre-computed-grid PDF.
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

                    from_expression->kinematics().set("x", x);
                    from_expression->kinematics().set("y", y);
                    from_grid->kinematics().set("x", x);
                    from_grid->kinematics().set("y", y);

                    const double reference = gaussian_density(x, mu_x, sigma_cx) * gaussian_density(y, mu_y, sigma_cy);

                    TEST_CHECK_RELATIVE_ERROR(std::exp(from_expression->evaluate() - from_expression->normalization()), reference, tolerance);
                    TEST_CHECK_RELATIVE_ERROR(std::exp(from_grid->evaluate() - from_grid->normalization()), reference, tolerance);
                }
            }

            const double reference_norm = gaussian_norm(mu_x, sigma_x, x_min, x_max) * gaussian_norm(mu_y, sigma_y, y_min, y_max);

            TEST_CHECK_RELATIVE_ERROR(std::exp(from_expression->normalization()), reference_norm, tolerance);
            TEST_CHECK_RELATIVE_ERROR(std::exp(from_grid->normalization()), reference_norm, tolerance);
        }
} detector_level_pdf_2d_test;

// A truth PDF supported on a strict sub-interval of the padded grid: its normalization over the grid
// must equal that over the support. B->Dlnu::NormalizationPDF (~5.7e-12) is below hcubature's absolute
// tolerance of 1e-9, so it accepts its first, unrefined rule and misses the integral by ~5% once the
// range extends past the support.
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

// A resolution given as an expression is sampled through the cache, resampled when its parameters
// change, and agrees with the same kernel given as a pre-computed grid.
class DetectorLevelPDFExpressionTest : public TestCase
{
    public:
        DetectorLevelPDFExpressionTest() :
            TestCase("detector_level_pdf_expression_test")
        {
        }

        virtual void
        run() const override
        {
            static const std::size_t M       = 64;
            static const double      lo      = -20.0;
            static const double      hi      = +20.0;
            static const double      spacing = (hi - lo) / static_cast<double>(M - 1);

            const DetectorLevelPDF::Axis axis{ "x", lo, hi, M };

            // the kernel sampled on the centred offset grid
            auto sampled = [&](const std::function<double(const double &)> & kernel)
            {
                std::vector<double> result(M);
                for (std::size_t i = 0; i < M; ++i)
                {
                    result[i] = kernel((static_cast<double>(i) - static_cast<double>(M / 2)) * spacing);
                }
                return result;
            };

            auto check_grids = [&](const std::shared_ptr<DetectorLevelPDF> & a, const std::shared_ptr<DetectorLevelPDF> & b)
            {
                a->update_grid();
                b->update_grid();
                for (std::size_t i = 0; i < M; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(a->grid()[i], b->grid()[i], 1e-14);
                }
            };

            // kernels and arbitrary expressions agree with their pre-computed grids
            {
                Parameters      p = Parameters::Defaults();
                ObservableCache cache(p);

                auto kernel = std::dynamic_pointer_cast<DetectorLevelPDF>(
                        DetectorLevelPDF::make_from_expression(cache, "TestGaussian1D::P(x)", Options{}, { axis }, "Kernel::Gaussian({x}, 0, 2)"));
                auto kernel_grid = std::dynamic_pointer_cast<DetectorLevelPDF>(
                        DetectorLevelPDF::make_1d(cache, "TestGaussian1D::P(x)", Options{}, axis, sampled([](const double & u) { return gaussian_density(u, 0.0, 2.0); })));
                auto other = std::dynamic_pointer_cast<DetectorLevelPDF>(
                        DetectorLevelPDF::make_from_expression(cache, "TestGaussian1D::P(x)", Options{}, { axis }, "exp(0 - {x}^2 / 8)"));
                auto other_grid = std::dynamic_pointer_cast<DetectorLevelPDF>(
                        DetectorLevelPDF::make_1d(cache, "TestGaussian1D::P(x)", Options{}, axis, sampled([](const double & u) { return std::exp(-u * u / 8.0); })));

                cache.update();
                check_grids(kernel, kernel_grid);
                check_grids(other, other_grid);

                ObservableCache clone_cache(p);
                auto            clone = std::dynamic_pointer_cast<DetectorLevelPDF>(kernel->clone(clone_cache));
                clone_cache.update();
                check_grids(clone, kernel_grid);
            }

            // the kernel is resampled when its parameters change
            {
                Parameters      p = Parameters::Defaults();
                ObservableCache cache(p);

                auto pdf = std::dynamic_pointer_cast<DetectorLevelPDF>(
                        DetectorLevelPDF::make_from_expression(cache, "TestGaussian1D::P(x)", Options{}, { axis }, "Kernel::Gaussian({x}, [[TestBox1D::centre]], 2)"));
                auto shifted = std::dynamic_pointer_cast<DetectorLevelPDF>(
                        DetectorLevelPDF::make_1d(cache, "TestGaussian1D::P(x)", Options{}, axis, sampled([](const double & u) { return gaussian_density(u, 3.0, 2.0); })));

                p["TestBox1D::centre"] = 3.0;
                cache.update();
                check_grids(pdf, shifted);
            }

            // an invalid expression leaves the cache untouched
            {
                Parameters      p = Parameters::Defaults();
                ObservableCache pristine(p);

                TEST_CHECK_THROWS(ParsingError, DetectorLevelPDF::make_from_expression(pristine, "TestGaussian1D::P(x)", Options{}, { axis }, "Kernel::Gaussian({x}, 0"));
                TEST_CHECK_EQUAL(pristine.size(), 0u);
            }
        }
} detector_level_pdf_expression_test;

// Poles, with given or extracted coefficients, on either axis of grids of rank 1 and 2, subtracted and convolved analytically.
class DetectorLevelPDFPolesTest : public TestCase
{
    public:
        DetectorLevelPDFPolesTest() :
            TestCase("detector_level_pdf_poles_test")
        {
        }

        virtual void
        run() const override
        {
            Observables().insert("DetectorLevelPDFTest::Re_z_p", "", Unit::None(), Options(), "[[TestGaussian1D::sigma]]");
            Observables().insert("DetectorLevelPDFTest::Im_z_p", "", Unit::None(), Options(), "0 - [[TestGaussian1D::sigma]]");
            Observables().insert("DetectorLevelPDFTest::Re_c", "", Unit::None(), Options(), "0.5");
            Observables().insert("DetectorLevelPDFTest::Im_c", "", Unit::None(), Options(), "0.25");
            Observables().insert("DetectorLevelPDFTest::Re_c_x", "", Unit::None(), Options(), "{x}");

            static const std::size_t M  = 64;
            static const double      lo = -20.0;
            static const double      hi = +20.0;

            const std::vector<DetectorLevelPDF::Axis> axes{
                DetectorLevelPDF::Axis{ "x", lo, hi, M }
            };
            const std::string resolution = "Kernel::Gaussian({x}, 0, 2)";

            const DetectorLevelPDF::Pole valid{ "x", "DetectorLevelPDFTest::Re_z_p", "DetectorLevelPDFTest::Im_z_p", "DetectorLevelPDFTest::Re_c", "DetectorLevelPDFTest::Im_c" };

            // Invalid poles leave the cache untouched.
            {
                Parameters      p = Parameters::Defaults();
                ObservableCache pristine(p);

                const DetectorLevelPDF::Pole unknown_variable{ "y",
                                                               "DetectorLevelPDFTest::Re_z_p",
                                                               "DetectorLevelPDFTest::Im_z_p",
                                                               "DetectorLevelPDFTest::Re_c",
                                                               "DetectorLevelPDFTest::Im_c" };
                const DetectorLevelPDF::Pole unknown_name{ "x",
                                                           "DetectorLevelPDFTest::Re_z_p",
                                                           "DetectorLevelPDFTest::NOSUCH",
                                                           "DetectorLevelPDFTest::Re_c",
                                                           "DetectorLevelPDFTest::Im_c" };
                DetectorLevelPDF::Pole       one_coefficient{ "x", "DetectorLevelPDFTest::Re_z_p", "DetectorLevelPDFTest::Im_z_p" };
                one_coefficient.coefficient_re = QualifiedName("DetectorLevelPDFTest::Re_c");

                TEST_CHECK_THROWS(InternalError, DetectorLevelPDF::make_from_expression(pristine, "TestGaussian1D::P(x)", Options{}, axes, resolution, { unknown_variable }));
                TEST_CHECK_THROWS(InternalError, DetectorLevelPDF::make_from_expression(pristine, "TestGaussian1D::P(x)", Options{}, axes, resolution, { unknown_name }));
                TEST_CHECK_THROWS(InternalError, DetectorLevelPDF::make_from_expression(pristine, "TestGaussian1D::P(x)", Options{}, axes, resolution, { one_coefficient }));

                // all poles must lie in the same variable
                const std::vector<DetectorLevelPDF::Axis> axes_xy{
                    DetectorLevelPDF::Axis{ "x", lo, hi, M },
                    DetectorLevelPDF::Axis{ "y", lo, hi, M }
                };
                const DetectorLevelPDF::Pole pole_y{ "y",
                                                     "DetectorLevelPDFTest::Re_z_p",
                                                     "DetectorLevelPDFTest::Im_z_p",
                                                     "DetectorLevelPDFTest::Re_c",
                                                     "DetectorLevelPDFTest::Im_c" };
                TEST_CHECK_THROWS(InternalError,
                                  DetectorLevelPDF::make_from_expression(pristine,
                                                                         "TestGaussian2D::P(x,y)",
                                                                         Options{},
                                                                         axes_xy,
                                                                         "Kernel::Gaussian({x}, 0, 2) * Kernel::Gaussian({y}, 0, 2)",
                                                                         { valid, pole_y }));

                // the resolution must be a linear combination of kernels in the pole's offset variable
                TEST_CHECK_THROWS(InternalError, DetectorLevelPDF::make_from_expression(pristine, "TestGaussian1D::P(x)", Options{}, axes, "exp(0 - {x}^2)", { valid }));
                TEST_CHECK_THROWS(InternalError,
                                  DetectorLevelPDF::make_from_expression(pristine, "TestGaussian1D::P(x)", Options{}, axes, "Kernel::Gaussian({x}, 0, 2) + 0.1", { valid }));

                // coefficients do not receive the pole variable
                TEST_CHECK_THROWS(UnknownKinematicVariableError,
                                  DetectorLevelPDF::make_from_expression(pristine,
                                                                         "TestGaussian1D::P(x)",
                                                                         Options{
                },
                                                                         axes,
                                                                         resolution,
                                                                         { DetectorLevelPDF::Pole{ "x",
                                                                                                   "DetectorLevelPDFTest::Re_z_p",
                                                                                                   "DetectorLevelPDFTest::Im_z_p",
                                                                                                   "DetectorLevelPDFTest::Re_c_x",
                                                                                                   "DetectorLevelPDFTest::Im_c" } }));

                TEST_CHECK_EQUAL(pristine.size(), 0u);
            }

            // A pole is reported, and cloned along with the PDF.
            {
                Parameters      p = Parameters::Defaults();
                ObservableCache cache(p);

                auto pdf =
                        std::dynamic_pointer_cast<DetectorLevelPDF>(DetectorLevelPDF::make_from_expression(cache, "TestGaussian1D::P(x)", Options{}, axes, resolution, { valid }));

                TEST_CHECK_EQUAL(pdf->poles().size(), 1u);
                TEST_CHECK_EQUAL(pdf->poles()[0].variable, "x");
                TEST_CHECK_EQUAL(pdf->poles()[0].position_re.full(), "DetectorLevelPDFTest::Re_z_p");
                TEST_CHECK_EQUAL(pdf->poles()[0].position_im.full(), "DetectorLevelPDFTest::Im_z_p");
                TEST_CHECK_EQUAL(pdf->poles()[0].coefficient_re->full(), "DetectorLevelPDFTest::Re_c");
                TEST_CHECK_EQUAL(pdf->poles()[0].coefficient_im->full(), "DetectorLevelPDFTest::Im_c");

                ObservableCache other(p);
                TEST_CHECK_EQUAL(std::dynamic_pointer_cast<DetectorLevelPDF>(pdf->clone(other))->poles().size(), 1u);
            }

            // A parameter serves as a position; Im z_p >= 0 is rejected on update.
            {
                Parameters      p = Parameters::Defaults();
                ObservableCache cache(p);

                const DetectorLevelPDF::Pole pole{ "x", "DetectorLevelPDFTest::Re_z_p", "TestGaussian1D::mu", "DetectorLevelPDFTest::Re_c", "DetectorLevelPDFTest::Im_c" };

                auto pdf =
                        std::dynamic_pointer_cast<DetectorLevelPDF>(DetectorLevelPDF::make_from_expression(cache, "TestGaussian1D::P(x)", Options{}, axes, resolution, { pole }));

                p["TestGaussian1D::mu"] = 0.0;
                cache.update();
                TEST_CHECK_THROWS(InternalError, pdf->update_grid());

                p["TestGaussian1D::mu"] = -1.0;
                cache.update();
                pdf->update_grid();
            }

            // The kernels must integrate to the grid sum of the resolution.
            {
                Parameters      p = Parameters::Defaults();
                ObservableCache cache(p);

                auto weighted  = DetectorLevelPDF::make_from_expression(cache,
                                                                       "TestGaussian1D::P(x)",
                                                                       Options{},
                                                                       axes,
                                                                       "0.25 * Kernel::Gaussian({x}, 0, 2) + Kernel::Gaussian({x}, 1, 3) / 2",
                                                                        { valid });
                auto truncated = DetectorLevelPDF::make_from_expression(cache, "TestGaussian1D::P(x)", Options{}, axes, "Kernel::Gaussian({x}, 0, 15)", { valid });
                cache.update();
                TEST_CHECK_NO_THROW(weighted->evaluate());
                TEST_CHECK_THROWS(InternalError, truncated->evaluate());
            }

            // For p(x) = |1 + R / (x - z)|^2 = 1 + Re[c / (x - z)], c = 2 R + i |R|^2 / gamma, the subtracted
            // remainder is constant, so the smeared grid is exactly 1 + Re[c S(x - z)], with
            // S(zeta) = -i sqrt(pi / 2) / sigma w(zeta / (sqrt(2) sigma)) for a Gaussian resolution.
            {
                const double          x_p = 0.3, gamma = 1.0e-3, sigma = 0.5;
                const complex<double> R = 0.01 * std::polar(1.0, M_PI / 4.0);
                const complex<double> z(x_p, -gamma);
                const complex<double> c = 2.0 * R + 1.0i * std::norm(R) / gamma;

                const std::string pole    = std::format("(({} * ({{x}} - {}) + {} * {}) / (({{x}} - {})^2 + {}^2))", R.real(), x_p, R.imag(), gamma, x_p, gamma);
                const std::string pole_im = std::format("(({} * ({{x}} - {}) - {} * {}) / (({{x}} - {})^2 + {}^2))", R.imag(), x_p, R.real(), gamma, x_p, gamma);
                Observables().insert("DetectorLevelPDFTest::p(x)", "", Unit::None(), Options(), std::format("(1 + {})^2 + {}^2", pole, pole_im));
                Observables().insert("DetectorLevelPDFTest::N", "", Unit::None(), Options(), "{x_max} - {x_min}");
                Observables().insert("DetectorLevelPDFTest::Re_z", "", Unit::None(), Options(), std::format("{}", x_p));
                Observables().insert("DetectorLevelPDFTest::Im_z", "", Unit::None(), Options(), std::format("{}", -gamma));
                Observables().insert("DetectorLevelPDFTest::Re_c_exact", "", Unit::None(), Options(), std::format("{}", c.real()));
                Observables().insert("DetectorLevelPDFTest::Im_c_exact", "", Unit::None(), Options(), std::format("{}", c.imag()));
                SignalPDFs().insert("DetectorLevelPDFTest::P(x)", "", Options(), "DetectorLevelPDFTest::p(x)", { "x" }, "DetectorLevelPDFTest::N", { "x_min", "x_max" });

                const std::size_t                         points = 1024;
                const std::vector<DetectorLevelPDF::Axis> fine_axes{
                    DetectorLevelPDF::Axis{ "x", lo, hi, points }
                };
                const std::string            gaussian = std::format("Kernel::Gaussian({{x}}, 0, {})", sigma);
                const DetectorLevelPDF::Pole exact{ "x",
                                                    "DetectorLevelPDFTest::Re_z",
                                                    "DetectorLevelPDFTest::Im_z",
                                                    "DetectorLevelPDFTest::Re_c_exact",
                                                    "DetectorLevelPDFTest::Im_c_exact" };

                Parameters      p = Parameters::Defaults();
                ObservableCache cache(p);
                auto            subtracted = std::dynamic_pointer_cast<DetectorLevelPDF>(
                        DetectorLevelPDF::make_from_expression(cache, "DetectorLevelPDFTest::P(x)", Options{}, fine_axes, gaussian, { exact }));
                auto naive =
                        std::dynamic_pointer_cast<DetectorLevelPDF>(DetectorLevelPDF::make_from_expression(cache, "DetectorLevelPDFTest::P(x)", Options{}, fine_axes, gaussian));
                // without a coefficient, it is extracted from the grid; |x - z|^2 p(x) = |x - z + R|^2 is interpolated exactly
                auto extracted = std::dynamic_pointer_cast<DetectorLevelPDF>(
                        DetectorLevelPDF::make_from_expression(cache,
                                                               "DetectorLevelPDFTest::P(x)",
                                                               Options{
                },
                                                               fine_axes,
                                                               gaussian,
                                                               { DetectorLevelPDF::Pole{ "x", "DetectorLevelPDFTest::Re_z", "DetectorLevelPDFTest::Im_z" } }));
                cache.update();
                subtracted->update_grid();
                naive->update_grid();
                extracted->update_grid();
                // the remainder is constant and the add-back is evaluated at the query points, so even points
                // between the grid nodes are exact
                const double        spacing = (hi - lo) / static_cast<double>(points - 1);
                std::vector<double> query;
                for (std::size_t k = 0; k + 1 < points; ++k)
                {
                    const double x = lo + (static_cast<double>(k) + 0.37) * spacing;
                    if (std::abs(x - x_p) <= 5.0 * sigma)
                    {
                        query.push_back(x);
                    }
                }
                const auto values = [&query](const std::shared_ptr<DetectorLevelPDF> & pdf) -> std::vector<double>
                {
                    const auto          interpolator = pdf->make_interpolator(query);
                    std::vector<double> result;
                    for (std::size_t i = 0; i < interpolator->size(); ++i)
                    {
                        result.push_back((*interpolator)(i));
                    }
                    return result;
                };
                const auto values_subtracted = values(subtracted), values_naive = values(naive), values_extracted = values(extracted);

                double max_subtracted = 0.0, max_naive = 0.0;
                for (std::size_t i = 0; i < query.size(); ++i)
                {
                    const complex<double> S        = -1.0i * std::sqrt(M_PI / 2.0) / sigma * faddeeva_w((query[i] - z) / (std::sqrt(2.0) * sigma));
                    const double          expected = 1.0 + (c * S).real();
                    max_subtracted                 = std::max(max_subtracted, std::abs(values_subtracted[i] / expected - 1.0));
                    max_naive                      = std::max(max_naive, std::abs(values_naive[i] / expected - 1.0));
                    TEST_CHECK_RELATIVE_ERROR(values_extracted[i], values_subtracted[i], 1.0e-10);
                }

                // the same holds for a single point
                auto variable = subtracted->kinematics()["x"];
                variable.set(query.front());
                TEST_CHECK_RELATIVE_ERROR(subtracted->evaluate_linear(), values_subtracted.front(), 1.0e-14);
                TEST_CHECK(max_subtracted < 1.0e-10);
                // the grid spacing of 0.04 does not resolve gamma = 1e-3
                TEST_CHECK(max_naive > 1.0e-2);

                // the integral over windows includes the integral of the add-back
                const std::vector<DetectorLevelPDF::Window> windows{ { { x_p - 1.3, x_p + 0.71 } }, { { x_p + 2.0, x_p + 3.0 } } };
                const double expected_integral = smeared_pole_integral(x_p - 1.3, x_p + 0.71, c, z, sigma) + smeared_pole_integral(x_p + 2.0, x_p + 3.0, c, z, sigma);
                TEST_CHECK_NEARLY_EQUAL(subtracted->integral(windows), std::log(expected_integral), 1.0e-10);
                TEST_CHECK_NEARLY_EQUAL(extracted->integral(windows), std::log(expected_integral), 1.0e-10);

                // extracting a coefficient requires grid points on either side of the pole
                Observables().insert("DetectorLevelPDFTest::Re_z_edge", "", Unit::None(), Options(), std::format("{}", hi - 0.01));
                auto edge = DetectorLevelPDF::make_from_expression(cache,
                                                                   "DetectorLevelPDFTest::P(x)",
                                                                   Options{
                },
                                                                   fine_axes,
                                                                   gaussian,
                                                                   { DetectorLevelPDF::Pole{ "x", "DetectorLevelPDFTest::Re_z_edge", "DetectorLevelPDFTest::Im_z" } });
                cache.update();
                TEST_CHECK_THROWS(InternalError, edge->evaluate());
            }

            // For p(x, y) = q(y) |1 + R / (x - z)|^2 with q(y) = 1 + cos(2 pi (y - y_min) / L) / 2, periodic on the grid,
            // the coefficient is c q(y) and the remainder is q(y). With the discrete resolution in y, whose cosine moment is
            // kappa, the smeared grid is exactly (1 + kappa cos(...) / 2) (1 + Re[c S(x - z)]).
            {
                const double          x_p = 0.3, gamma = 1.0e-3, sigma_x = 0.5, sigma_y = 0.1;
                const complex<double> R = 0.01 * std::polar(1.0, M_PI / 4.0);
                const complex<double> z(x_p, -gamma);
                const complex<double> c = 2.0 * R + 1.0i * std::norm(R) / gamma;

                const std::size_t points_x = 1024, points_y = 32;
                const double      y_min = -1.0, y_max = 1.0, spacing_y = (y_max - y_min) / static_cast<double>(points_y - 1);
                const double      omega = 2.0 * M_PI / (static_cast<double>(points_y) * spacing_y);

                const std::string q       = std::format("(1 + 0.5 * cos({} * ({{y}} - {})))", omega, y_min);
                const std::string pole    = std::format("(({} * ({{x}} - {}) + {} * {}) / (({{x}} - {})^2 + {}^2))", R.real(), x_p, R.imag(), gamma, x_p, gamma);
                const std::string pole_im = std::format("(({} * ({{x}} - {}) - {} * {}) / (({{x}} - {})^2 + {}^2))", R.imag(), x_p, R.real(), gamma, x_p, gamma);
                Observables().insert("DetectorLevelPDFTest::p(x,y)", "", Unit::None(), Options(), std::format("{} * ((1 + {})^2 + {}^2)", q, pole, pole_im));
                Observables().insert("DetectorLevelPDFTest::N2", "", Unit::None(), Options(), "({x_max} - {x_min}) * ({y_max} - {y_min})");
                Observables().insert("DetectorLevelPDFTest::Re_c(y)", "", Unit::None(), Options(), std::format("{} * {}", c.real(), q));
                Observables().insert("DetectorLevelPDFTest::Im_c(y)", "", Unit::None(), Options(), std::format("{} * {}", c.imag(), q));
                SignalPDFs().insert("DetectorLevelPDFTest::P(x,y)",
                                    "",
                                    Options(),
                                    "DetectorLevelPDFTest::p(x,y)",
                                    { "x", "y" },
                                    "DetectorLevelPDFTest::N2",
                                    { "x_min", "x_max", "y_min", "y_max" });

                const DetectorLevelPDF::Axis axis_x{ "x", lo, hi, points_x }, axis_y{ "y", y_min, y_max, points_y };
                const std::string            resolution_2d = std::format("Kernel::Gaussian({{x}}, 0, {}) * Kernel::Gaussian({{y}}, 0, {})", sigma_x, sigma_y);

                // the cosine moment of the y resolution, sampled on the centred offsets and normalized to unit sum
                double kernel_sum = 0.0, kernel_moment = 0.0;
                for (std::size_t m = 0; m < points_y; ++m)
                {
                    const double v       = (static_cast<double>(m) - static_cast<double>(points_y / 2)) * spacing_y;
                    const double weight  = std::exp(-0.5 * v * v / (sigma_y * sigma_y));
                    kernel_sum          += weight;
                    kernel_moment       += weight * std::cos(omega * v);
                }
                const double kappa = kernel_moment / kernel_sum;

                // the pole lies on the leading axis, and on the trailing one
                for (const bool x_leads : { true, false })
                {
                    const std::vector<DetectorLevelPDF::Axis> axes_2d =
                            x_leads ? std::vector<DetectorLevelPDF::Axis>{ axis_x, axis_y } : std::vector<DetectorLevelPDF::Axis>{ axis_y, axis_x };
                    // a point or an interval, in the order of the axes
                    const auto in_order = [x_leads](const auto & x, const auto & y) { return x_leads ? std::array{ x, y } : std::array{ y, x }; };

                    Parameters      p = Parameters::Defaults();
                    ObservableCache cache(p);
                    auto            given     = std::dynamic_pointer_cast<DetectorLevelPDF>(DetectorLevelPDF::make_from_expression(cache,
                                                                                                                    "DetectorLevelPDFTest::P(x,y)",
                                                                                                                    Options{
                    },
                                                                                                                    axes_2d,
                                                                                                                    resolution_2d,
                                                                                                                                   { DetectorLevelPDF::Pole{ "x",
                                                                                                                                              "DetectorLevelPDFTest::Re_z",
                                                                                                                                              "DetectorLevelPDFTest::Im_z",
                                                                                                                                              "DetectorLevelPDFTest::Re_c(y)",
                                                                                                                                              "DetectorLevelPDFTest::Im_c(y)" } }));
                    auto            extracted = std::dynamic_pointer_cast<DetectorLevelPDF>(
                            DetectorLevelPDF::make_from_expression(cache,
                                                                   "DetectorLevelPDFTest::P(x,y)",
                                                                   Options{
                    },
                                                                   axes_2d,
                                                                   resolution_2d,
                                                                              { DetectorLevelPDF::Pole{ "x", "DetectorLevelPDFTest::Re_z", "DetectorLevelPDFTest::Im_z" } }));
                    cache.update();
                    given->update_grid();
                    extracted->update_grid();

                    // between the grid nodes in x, where the remainder is constant, and on the nodes in y
                    const double        spacing_x = (hi - lo) / static_cast<double>(points_x - 1);
                    std::vector<double> query;
                    for (std::size_t i = 0; i + 1 < points_x; ++i)
                    {
                        for (std::size_t n = 0; n < points_y; ++n)
                        {
                            for (const double value : in_order(lo + (static_cast<double>(i) + 0.37) * spacing_x, y_min + static_cast<double>(n) * spacing_y))
                            {
                                query.push_back(value);
                            }
                        }
                    }
                    const auto interpolator_given = given->make_interpolator(query), interpolator_extracted = extracted->make_interpolator(query);

                    double max_given = 0.0, max_extracted = 0.0;
                    for (std::size_t i = 0; i < query.size() / 2; ++i)
                    {
                        const double          x         = query[x_leads ? 2 * i : 2 * i + 1];
                        const std::size_t     n         = i % points_y;
                        const complex<double> S         = -1.0i * std::sqrt(M_PI / 2.0) / sigma_x * faddeeva_w((x - z) / (std::sqrt(2.0) * sigma_x));
                        const double          smeared_x = 1.0 + (c * S).real();
                        const double          expected  = (1.0 + 0.5 * kappa * std::cos(omega * static_cast<double>(n) * spacing_y)) * smeared_x;
                        max_given                       = std::max(max_given, std::abs((*interpolator_given)(i) / expected - 1.0));
                        max_extracted                   = std::max(max_extracted, std::abs((*interpolator_extracted)(i) / expected - 1.0));
                    }
                    TEST_CHECK(max_given < 1.0e-10);
                    TEST_CHECK(max_extracted < 1.0e-10);

                    // the integral over a window aligned with the nodes in y: the trapezoidal rule in y times the integral in x
                    const std::size_t n_min = 3, n_max = 20;
                    double            trapezoid_y = 0.0;
                    for (std::size_t n = n_min; n <= n_max; ++n)
                    {
                        const double weight  = ((n_min == n) || (n_max == n)) ? 0.5 : 1.0;
                        trapezoid_y         += weight * spacing_y * (1.0 + 0.5 * kappa * std::cos(omega * static_cast<double>(n) * spacing_y));
                    }
                    const auto                                  intervals = in_order(std::pair{ x_p - 1.3, x_p + 0.71 },
                                                    std::pair{ y_min + static_cast<double>(n_min) * spacing_y, y_min + static_cast<double>(n_max) * spacing_y });
                    const std::vector<DetectorLevelPDF::Window> windows{
                        { intervals[0], intervals[1] }
                    };
                    const double expected_integral = trapezoid_y * smeared_pole_integral(x_p - 1.3, x_p + 0.71, c, z, sigma_x);
                    TEST_CHECK_NEARLY_EQUAL(given->integral(windows), std::log(expected_integral), 1.0e-10);
                    TEST_CHECK_NEARLY_EQUAL(extracted->integral(windows), std::log(expected_integral), 1.0e-10);
                }
            }

            // As for p(x) = |1 + R / (x - z)|^2 above, the smeared grid is 1 + Re[c S(x - z)] for any resolution; here S is
            // the Cauchy transform of the normalized resolution, a sum over its kernels for several of them.
            const std::size_t                         fine_points  = 1024;
            const double                              fine_spacing = (hi - lo) / static_cast<double>(fine_points - 1);
            const std::vector<DetectorLevelPDF::Axis> fine_axes{
                DetectorLevelPDF::Axis{ "x", lo, hi, fine_points }
            };
            // the largest relative deviation from the expected smeared grid between the grid nodes within [a, b]
            const auto max_deviation =
                    [&](const std::shared_ptr<DetectorLevelPDF> & pdf, const double & a, const double & b, const std::function<double(const double &)> & expected) -> double
            {
                std::vector<double> query;
                for (std::size_t k = 0; k + 1 < fine_points; ++k)
                {
                    const double x = lo + (static_cast<double>(k) + 0.37) * fine_spacing;
                    if ((a <= x) && (x <= b))
                    {
                        query.push_back(x);
                    }
                }

                const auto interpolator = pdf->make_interpolator(query);
                double     result       = 0.0;
                for (std::size_t i = 0; i < query.size(); ++i)
                {
                    result = std::max(result, std::abs((*interpolator)(i) / expected(query[i]) - 1.0));
                }

                return result;
            };
            const auto gaussian_transform = [](const complex<double> & zeta, const double & mu, const double & sigma) -> complex<double>
            { return -1.0i * std::sqrt(M_PI / 2.0) / sigma * faddeeva_w((zeta - mu) / (std::sqrt(2.0) * sigma)); };

            // An asymmetric resolution, against a brute-force convolution: a shifted Gaussian, and a shifted Crystal Ball
            // with unequal tails, whose add-back is accurate to first order in gamma / sigma = 4e-3.
            {
                const double          x_p = 0.3, gamma = 1.0e-3;
                const complex<double> R = 0.01 * std::polar(1.0, M_PI / 4.0);
                const complex<double> z(x_p, -gamma);
                const complex<double> c = 2.0 * R + 1.0i * std::norm(R) / gamma;

                const DetectorLevelPDF::Pole exact{ "x",
                                                    "DetectorLevelPDFTest::Re_z",
                                                    "DetectorLevelPDFTest::Im_z",
                                                    "DetectorLevelPDFTest::Re_c_exact",
                                                    "DetectorLevelPDFTest::Im_c_exact" };

                const double mu = 0.1, sigma = 0.25, alpha_L = 1.5, n_L = 4.0, alpha_R = 2.0, n_R = 5.0;

                struct Resolution
                {
                        std::string                           expression;
                        std::function<double(const double &)> kernel;
                        std::vector<double>                   points;
                        double                                tolerance;
                };

                const std::vector<Resolution> resolutions{
                    { std::format("Kernel::Gaussian({{x}}, {}, {})", mu, sigma),
                     [=](const double & u) { return std::exp(-0.5 * power_of<2>((u - mu) / sigma)); },
                     { mu },
                     1.0e-12 },
                    { std::format("Kernel::CrystalBall({{x}}, {}, {}, {}, {}, {}, {})", mu, sigma, alpha_L, n_L, alpha_R, n_R),
                     [=](const double & u) { return crystal_ball_shape(u, mu, sigma, alpha_L, n_L, alpha_R, n_R); },
                     { mu - alpha_L * sigma, mu, mu + alpha_R * sigma },
                     2.0e-6 },
                };

                for (const auto & resolution : resolutions)
                {
                    Parameters      p = Parameters::Defaults();
                    ObservableCache cache(p);
                    auto            pdf = std::dynamic_pointer_cast<DetectorLevelPDF>(
                            DetectorLevelPDF::make_from_expression(cache, "DetectorLevelPDFTest::P(x)", Options{}, fine_axes, resolution.expression, { exact }));
                    cache.update();
                    pdf->update_grid();

                    const auto expected = [&](const double & x) -> double
                    { return 1.0 + (c * brute_force_cauchy_transform(resolution.kernel, x - z, 1.0e4 * sigma, resolution.points)).real(); };
                    TEST_CHECK(max_deviation(pdf, x_p - 5.0 * sigma, x_p + 5.0 * sigma, expected) < resolution.tolerance);
                }
            }

            // A resolution of several kernels in the pole variable: S is the weighted sum of their Cauchy transforms.
            {
                const double          x_p = 0.3, gamma = 1.0e-3;
                const complex<double> R = 0.01 * std::polar(1.0, M_PI / 4.0);
                const complex<double> z(x_p, -gamma);
                const complex<double> c = 2.0 * R + 1.0i * std::norm(R) / gamma;

                const DetectorLevelPDF::Pole exact{ "x",
                                                    "DetectorLevelPDFTest::Re_z",
                                                    "DetectorLevelPDFTest::Im_z",
                                                    "DetectorLevelPDFTest::Re_c_exact",
                                                    "DetectorLevelPDFTest::Im_c_exact" };

                Parameters      p = Parameters::Defaults();
                ObservableCache cache(p);
                auto            pdf = std::dynamic_pointer_cast<DetectorLevelPDF>(
                        DetectorLevelPDF::make_from_expression(cache,
                                                               "DetectorLevelPDFTest::P(x)",
                                                               Options{},
                                                               fine_axes,
                                                               "0.25 * Kernel::Gaussian({x}, 0, 0.5) + Kernel::Gaussian({x}, 0.2, 0.8) / 2",
                                                                          { exact }));
                cache.update();
                pdf->update_grid();

                const auto expected = [&](const double & x) -> double
                { return 1.0 + (c * (0.25 * gaussian_transform(x - z, 0.0, 0.5) + 0.5 * gaussian_transform(x - z, 0.2, 0.8)) / 0.75).real(); };
                TEST_CHECK(max_deviation(pdf, x_p - 4.0, x_p + 4.0, expected) < 1.0e-10);
            }

            // Two poles in the same variable: p(x) = |1 + R_1 / (x - z_1) + R_2 / (x - z_2)|^2 = 1 + sum_j Re[c_j / (x - z_j)],
            // with c_j = 2 R_j (1 + sum_k conj(R_k) / (z_j - conj(z_k))), so that the smeared grid is 1 + sum_j Re[c_j S(x - z_j)].
            {
                const double                         sigma = 0.5;
                const std::array<double, 2>          x_p{ 0.3, 1.1 }, gamma{ 1.0e-3, 2.0e-3 };
                const std::array<complex<double>, 2> R{ 0.01 * std::polar(1.0, M_PI / 4.0), 0.008 * std::polar(1.0, -M_PI / 3.0) };
                const std::array<complex<double>, 2> z{ complex<double>(x_p[0], -gamma[0]), complex<double>(x_p[1], -gamma[1]) };
                std::array<complex<double>, 2>       c;
                for (std::size_t j = 0; j < 2; ++j)
                {
                    c[j] = 2.0 * R[j] * (1.0 + std::conj(R[0]) / (z[j] - std::conj(z[0])) + std::conj(R[1]) / (z[j] - std::conj(z[1])));
                }

                std::string re_sum = "1", im_sum = "0";
                for (std::size_t j = 0; j < 2; ++j)
                {
                    re_sum += std::format(" + (({} * ({{x}} - {}) + {} * {}) / (({{x}} - {})^2 + {}^2))", R[j].real(), x_p[j], R[j].imag(), gamma[j], x_p[j], gamma[j]);
                    im_sum += std::format(" + (({} * ({{x}} - {}) - {} * {}) / (({{x}} - {})^2 + {}^2))", R[j].imag(), x_p[j], R[j].real(), gamma[j], x_p[j], gamma[j]);

                    Observables().insert(std::format("DetectorLevelPDFTest::Re_z_{}", j + 1), "", Unit::None(), Options(), std::format("{}", x_p[j]));
                    Observables().insert(std::format("DetectorLevelPDFTest::Im_z_{}", j + 1), "", Unit::None(), Options(), std::format("{}", -gamma[j]));
                    Observables().insert(std::format("DetectorLevelPDFTest::Re_c_{}", j + 1), "", Unit::None(), Options(), std::format("{}", c[j].real()));
                    Observables().insert(std::format("DetectorLevelPDFTest::Im_c_{}", j + 1), "", Unit::None(), Options(), std::format("{}", c[j].imag()));
                }
                Observables().insert("DetectorLevelPDFTest::p_two(x)", "", Unit::None(), Options(), std::format("({})^2 + ({})^2", re_sum, im_sum));
                Observables().insert("DetectorLevelPDFTest::N_two", "", Unit::None(), Options(), "{x_max} - {x_min}");
                SignalPDFs()
                        .insert("DetectorLevelPDFTest::P_two(x)", "", Options(), "DetectorLevelPDFTest::p_two(x)", { "x" }, "DetectorLevelPDFTest::N_two", { "x_min", "x_max" });

                std::vector<DetectorLevelPDF::Pole> poles;
                for (std::size_t j = 1; j <= 2; ++j)
                {
                    poles.push_back(DetectorLevelPDF::Pole{ "x",
                                                            std::format("DetectorLevelPDFTest::Re_z_{}", j),
                                                            std::format("DetectorLevelPDFTest::Im_z_{}", j),
                                                            std::format("DetectorLevelPDFTest::Re_c_{}", j),
                                                            std::format("DetectorLevelPDFTest::Im_c_{}", j) });
                }

                Parameters      p = Parameters::Defaults();
                ObservableCache cache(p);
                auto            pdf = std::dynamic_pointer_cast<DetectorLevelPDF>(DetectorLevelPDF::make_from_expression(cache,
                                                                                                              "DetectorLevelPDFTest::P_two(x)",
                                                                                                              Options{},
                                                                                                              fine_axes,
                                                                                                              std::format("Kernel::Gaussian({{x}}, 0, {})", sigma),
                                                                                                              poles));
                cache.update();
                pdf->update_grid();

                const auto expected = [&](const double & x) -> double
                { return 1.0 + (c[0] * gaussian_transform(x - z[0], 0.0, sigma)).real() + (c[1] * gaussian_transform(x - z[1], 0.0, sigma)).real(); };
                TEST_CHECK(max_deviation(pdf, x_p[0] - 5.0 * sigma, x_p[1] + 5.0 * sigma, expected) < 1.0e-10);
            }
        }
} detector_level_pdf_poles_test;
