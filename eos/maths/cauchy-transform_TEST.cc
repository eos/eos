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
#include <eos/utils/exception.hh>

#include <test/test.hh>

#include <array>
#include <vector>

using namespace test;
using namespace eos;

class CauchyTransformTest : public TestCase
{
    public:
        CauchyTransformTest() :
            TestCase("cauchy_transform_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_THROWS(InternalError, CauchyTransform::make("Kernel::NoSuchKernel", 1.0));

            // Gaussian, sigma = 5e-3 and gamma = 1.5e-4, against adaptive quadrature
            {
                auto transform = CauchyTransform::make("Kernel::Gaussian", 1.0);

                struct Reference
                {
                        double          mu, y;
                        complex<double> value;
                };

                const std::vector<Reference> references{
                    { 0.0,    0.0,         { 0.0, -244.773851 } },
                    { 0.0,  0.004,  { 125.586174, -179.164649 } },
                    { 0.0, -0.009, { -136.012493, -51.0448018 } },
                    // the location shifts the argument
                    { 0.3,  0.304,  { 125.586174, -179.164649 } },
                };

                for (const auto & r : references)
                {
                    const std::array<double, 2> parameters{ r.mu, 5.0e-3 };
                    transform->prepare(parameters);

                    const complex<double> value = transform->evaluate(complex<double>(r.y, 1.5e-4));
                    TEST_CHECK_NEARLY_EQUAL(value.real(), r.value.real(), 1.0e-6 * std::abs(r.value));
                    TEST_CHECK_NEARLY_EQUAL(value.imag(), r.value.imag(), 1.0e-6 * std::abs(r.value));
                }

                // far from the kernel, S(zeta) approaches 1 / zeta
                const std::array<double, 2> parameters{ 0.0, 5.0e-3 };
                transform->prepare(parameters);
                const complex<double> zeta(1.0, 1.0e-4);
                TEST_CHECK_RELATIVE_ERROR(std::abs(transform->evaluate(zeta) * zeta), 1.0, 1.0e-4);
            }

            // Crystal Ball, sigma = 5e-3 and gamma = 1.25e-5, against adaptive quadrature; the first order in
            // gamma / sigma is kept, and the table reaches 200 and 10 sigma, respectively
            for (const double extent : { 1.0, 0.05 })
            {
                auto transform = CauchyTransform::make("Kernel::CrystalBall", extent);

                struct Reference
                {
                        double          mu, y;
                        complex<double> value;
                };

                const std::vector<Reference> references{
                    { 0.0,    0.0,  { 2.695766568986e+00, -2.329943298856e+02 } },
                    { 0.0,  0.004,  { 1.224152880653e+02, -1.693028353130e+02 } },
                    { 0.0, -0.009, { -1.182138633527e+02, -4.995865998241e+01 } },
                    { 0.0,   0.02,  { 5.315033078765e+01, -1.710099772471e+00 } },
                    { 0.0,  -0.05, { -2.083558777097e+01, -5.292126267826e-01 } },
                    { 0.3,  0.304,  { 1.224152880653e+02, -1.693028353130e+02 } },
                    { 0.0,    0.5,  { 1.995598200143e+00, -5.009775430005e-05 } },
                };

                for (const auto & r : references)
                {
                    const std::array<double, 6> parameters{ r.mu, 5.0e-3, 1.5, 3.0, 2.0, 5.0 };
                    transform->prepare(parameters);

                    const complex<double> value = transform->evaluate(complex<double>(r.y, 1.25e-5));
                    TEST_CHECK_NEARLY_EQUAL(value.real(), r.value.real(), 2.0e-6 * std::abs(r.value));
                    TEST_CHECK_NEARLY_EQUAL(value.imag(), r.value.imag(), 2.0e-6 * std::abs(r.value));
                }
            }

            // Crystal Ball beyond its table, two ulps below a breakpoint of the principal value at s = 10
            {
                auto                        transform = CauchyTransform::make("Kernel::CrystalBall", 5.0);
                const std::array<double, 6> parameters{ 0.0, 1.0, 1.5, 3.0, 2.0, 5.0 };
                transform->prepare(parameters);

                const complex<double> below = transform->evaluate(complex<double>(9.999999999999996, 1.0e-3));
                const complex<double> at    = transform->evaluate(complex<double>(10.0, 1.0e-3));
                TEST_CHECK_RELATIVE_ERROR(below.real(), at.real(), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(below.imag(), at.imag(), 1.0e-10);
            }

            // Crystal Ball with heavy tails, n_L = 1.2, which reach far beyond the table, and with steep tails,
            // n = 200, for which (n / alpha)^n overflows; as above, but heavy tails enlarge the second order in gamma / sigma
            {
                struct Reference
                {
                        double          y;
                        complex<double> value;
                };

                const std::vector<std::pair<std::array<double, 6>, std::vector<Reference>>> cases{
                    {     { 0.0, 5.0e-3, 1.5, 1.2, 2.0, 5.0 },
                     {
                     { 0.0, { 4.974401021918e+00, -1.710844382362e+02 } },
                     { 0.004, { 9.245897254465e+01, -1.243165670510e+02 } },
                     { -0.009, { -8.029189331370e+01, -3.807197608631e+01 } },
                     { 0.02, { 4.077076550540e+01, -1.256090710662e+00 } },
                     { -0.05, { -1.620610118543e+01, -2.934829464808e+00 } },
                     { 0.5, { 1.712955044919e+00, -4.165156000222e-05 } },
                     } },
                    { { 0.0, 5.0e-3, 1.5, 200.0, 2.0, 200.0 },
                     {
                     { 0.0, { 1.132300682109e+00, -2.442178599667e+02 } },
                     { 0.004, { 1.271200383455e+02, -1.774581715560e+02 } },
                     { -0.009, { -1.295913539339e+02, -5.080562747822e+01 } },
                     { 0.02, { 5.422228119609e+01, -6.715885037435e-01 } },
                     { -0.05, { -2.035498135065e+01, -5.597572133275e-03 } },
                     { 0.5, { 1.999359473974e+00, -4.997396313636e-05 } },
                     } },
                };

                for (const auto & [parameters, references] : cases)
                {
                    auto transform = CauchyTransform::make("Kernel::CrystalBall", 1.0);
                    transform->prepare(parameters);

                    for (const auto & r : references)
                    {
                        const complex<double> value = transform->evaluate(complex<double>(r.y, 1.25e-5));
                        TEST_CHECK_NEARLY_EQUAL(value.real(), r.value.real(), 5.0e-6 * std::abs(r.value));
                        TEST_CHECK_NEARLY_EQUAL(value.imag(), r.value.imag(), 5.0e-6 * std::abs(r.value));
                    }
                }
            }
        }
} cauchy_transform_test;
