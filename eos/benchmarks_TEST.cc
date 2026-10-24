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

#include <eos/benchmarks.hh>
#include <eos/maths/complex.hh>
#include <eos/observable.hh>
#include <eos/signal-pdf.hh>

#include <test/test.hh>

#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class BenchmarkB1Test : public TestCase
{
    public:
        BenchmarkB1Test() :
            TestCase("benchmark_b1_test")
        {
        }

        virtual void
        run() const
        {
            register_benchmarks();
            // registering twice has no further effect
            register_benchmarks();

            Parameters p = Parameters::Defaults();

            // reference values from an independent implementation in Python, for the inputs as doubles
            {
                Kinematics k{
                    { "m_mumu", 0.0 }
                };
                auto pdf = Observable::make("B1::P(m_mumu)@RvDV:2026A", p, k, Options{});

                const std::vector<std::array<double, 2>> references{
                    {     0.5, 9.5296627272080705e-01 },
                    {     3.0, 8.3577157861915513e-02 },
                    {     3.1, 3.9252483381915365e+06 },
                    { 3.10003, 2.8880108292038264e+06 },
                    {     3.7, 1.0696123350256778e+05 },
                    {     4.5, 1.0356561162982925e+00 },
                };
                for (const auto & [x, value] : references)
                {
                    k.set("m_mumu", x);
                    TEST_CHECK_RELATIVE_ERROR(pdf->evaluate(), value, 1.0e-12);
                }

                // the PDF vanishes outside the box
                for (const double x : { -0.1, 4.7 })
                {
                    k.set("m_mumu", x);
                    TEST_CHECK_EQUAL(pdf->evaluate(), 0.0);
                }
            }

            // the normalization agrees with a numerical integration in Python
            {
                Kinematics k{
                    { "m_mumu_min", 0.0 },
                    { "m_mumu_max", 0.0 }
                };
                auto norm = Observable::make("B1::N@RvDV:2026A", p, k, Options{});

                const std::vector<std::array<double, 3>> references{
                    { -0.25, 4.849, 6.7155328263836616e+02 },
                    {   1.0,   3.5, 6.1896748250322108e+02 },
                    { 3.099, 3.101, 5.9696880890434650e+02 },
                    {  3.69,  3.71, 4.9949706645975568e+01 },
                };
                for (const auto & [lo, hi, value] : references)
                {
                    k.set("m_mumu_min", lo);
                    k.set("m_mumu_max", hi);
                    TEST_CHECK_RELATIVE_ERROR(norm->evaluate(), value, 1.0e-10);
                }
            }

            // the poles are z_j = x_j - i gamma_j
            {
                const std::vector<std::array<double, 3>> references{
                    { 1.0, 3.1, -5.0e-5 },
                    { 2.0, 3.7, -1.5e-4 },
                };
                for (const auto & [j, re, im] : references)
                {
                    const std::string index = std::to_string(static_cast<int>(j));
                    TEST_CHECK_RELATIVE_ERROR(Observable::make("B1::Re{mu_" + index + "}@RvDV:2026A", p, Kinematics(), Options{})->evaluate(), re, 1.0e-14);
                    TEST_CHECK_RELATIVE_ERROR(Observable::make("B1::Im{mu_" + index + "}@RvDV:2026A", p, Kinematics(), Options{})->evaluate(), im, 1.0e-14);
                }
            }

            // with the coefficients, P - Re[c_j / (x - z_j)] is smooth across the pole
            for (const std::string index : { "1", "2" })
            {
                const double re_mu = Observable::make("B1::Re{mu_" + index + "}@RvDV:2026A", p, Kinematics(), Options{})->evaluate();
                const double im_mu = Observable::make("B1::Im{mu_" + index + "}@RvDV:2026A", p, Kinematics(), Options{})->evaluate();
                const double re_c  = Observable::make("B1::Re{c_" + index + "}@RvDV:2026A", p, Kinematics(), Options{})->evaluate();
                const double im_c  = Observable::make("B1::Im{c_" + index + "}@RvDV:2026A", p, Kinematics(), Options{})->evaluate();

                Kinematics k{
                    { "m_mumu", re_mu }
                };
                auto pdf = Observable::make("B1::P(m_mumu)@RvDV:2026A", p, k, Options{});

                const auto remainder = [&](const double & x) -> double
                {
                    k.set("m_mumu", x);
                    return pdf->evaluate() - (complex<double>(re_c, im_c) / (x - complex<double>(re_mu, im_mu))).real();
                };

                // the PDF varies by more than a factor of 100 over +-10 gamma; the remainder, which keeps the other pole, by less than 1e-2
                const double gamma = -im_mu, at_pole = remainder(re_mu);
                for (const double offset : { -10.0, -1.0, 1.0, 10.0 })
                {
                    TEST_CHECK_RELATIVE_ERROR(remainder(re_mu + offset * gamma), at_pole, 1.0e-2);
                }
            }

            // the signal PDF is registered and normalized
            {
                Kinematics k{
                    {     "m_mumu",   0.5 },
                    { "m_mumu_min", -0.25 },
                    { "m_mumu_max", 4.849 }
                };
                auto pdf = SignalPDF::make("benchmark::B1@RvDV:2026A", p, k, Options{});
                TEST_CHECK(pdf.get() != nullptr);
                TEST_CHECK_RELATIVE_ERROR(std::exp(pdf->evaluate() - pdf->normalization()), 9.5296627272080705e-01 / 6.7155328263836616e+02, 1.0e-10);
            }
        }
} benchmark_b1_test;
