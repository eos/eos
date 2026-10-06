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

class BenchmarkB2Test : public TestCase
{
    public:
        BenchmarkB2Test() :
            TestCase("benchmark_b2_test")
        {
        }

        virtual void
        run() const
        {
            register_benchmarks();
            // registering twice has no further effect
            register_benchmarks();

            Parameters p = Parameters::Defaults();

            // reference values from an independent implementation in Python
            {
                Kinematics k{
                    {       "m_mumu", 0.0 },
                    { "cos(theta_l)", 0.0 }
                };
                auto pdf = Observable::make("B2::P(m_mumu,cos(theta_l))@RvDV:2026A", p, k, Options{});

                const std::vector<std::array<double, 3>> references{
                    { 1.0,  0.3, 2.670137119275402e+04 },
                    { 3.0, -0.7, 5.040137275152078e+03 },
                    { 3.1,  0.0, 1.660319093725659e+11 },
                    { 4.5,  0.9, 3.734501255305918e+03 },
                };
                for (const auto & [m, c, value] : references)
                {
                    k.set("m_mumu", m);
                    k.set("cos(theta_l)", c);
                    TEST_CHECK_RELATIVE_ERROR(pdf->evaluate(), value, 1.0e-12);
                }

                // the PDF vanishes outside the phase space and outside |cos(theta_l)| <= 1
                for (const auto & [m, c] : std::vector<std::array<double, 2>>{
                         { 0.2,   0.0 },
                         { 4.8,   0.0 },
                         { 2.0,  1.01 },
                         { 2.0, -1.01 }
                })
                {
                    k.set("m_mumu", m);
                    k.set("cos(theta_l)", c);
                    TEST_CHECK_EQUAL(pdf->evaluate(), 0.0);
                }
            }

            // the normalization agrees with a brute-force integration in Python
            {
                Kinematics k{
                    {       "m_mumu_min",  0.0 },
                    {       "m_mumu_max",  5.0 },
                    { "cos(theta_l)_min", -1.2 },
                    { "cos(theta_l)_max",  1.2 }
                };
                auto norm = Observable::make("B2::N@RvDV:2026A", p, k, Options{});
                TEST_CHECK_RELATIVE_ERROR(norm->evaluate(), 3.787096935763e+07, 1.0e-8);

                k.set("m_mumu_min", 2.0);
                k.set("m_mumu_max", 4.0);
                k.set("cos(theta_l)_min", -0.5);
                k.set("cos(theta_l)_max", 0.8);
                TEST_CHECK_RELATIVE_ERROR(norm->evaluate(), 3.081454711403e+07, 1.0e-8);

                // ranges of a few widths around a pole, where the singular terms dominate
                k.set("m_mumu_min", 3.099);
                k.set("m_mumu_max", 3.101);
                k.set("cos(theta_l)_min", -1.0);
                k.set("cos(theta_l)_max", 1.0);
                TEST_CHECK_RELATIVE_ERROR(norm->evaluate(), 3.366785110072e+07, 1.0e-8);

                k.set("m_mumu_min", 3.0999);
                k.set("m_mumu_max", 3.1001);
                k.set("cos(theta_l)_min", -0.3);
                k.set("cos(theta_l)_max", 0.6);
                TEST_CHECK_RELATIVE_ERROR(norm->evaluate(), 1.505502810363e+07, 1.0e-8);
            }

            // the poles in m_mumu are the exact square roots of s_j = m_j^2 - i m_j Gamma_j
            {
                const std::vector<std::array<double, 3>> references{
                    { 1.0, 3.100000000403226, -4.9999999993496355e-05 },
                    { 2.0, 3.700000003040541, -1.4999999987673482e-04 },
                };
                for (const auto & [j, re, im] : references)
                {
                    const std::string index = std::to_string(static_cast<int>(j));
                    TEST_CHECK_RELATIVE_ERROR(Observable::make("B2::Re{mu_" + index + "}@RvDV:2026A", p, Kinematics(), Options{})->evaluate(), re, 1.0e-14);
                    TEST_CHECK_RELATIVE_ERROR(Observable::make("B2::Im{mu_" + index + "}@RvDV:2026A", p, Kinematics(), Options{})->evaluate(), im, 1.0e-10);
                }
            }

            // with the coefficients, P - Re[c_j / (m - mu_j)] is smooth across the pole
            for (const std::string index : { "1", "2" })
            {
                const double re_mu = Observable::make("B2::Re{mu_" + index + "}@RvDV:2026A", p, Kinematics(), Options{})->evaluate();
                const double im_mu = Observable::make("B2::Im{mu_" + index + "}@RvDV:2026A", p, Kinematics(), Options{})->evaluate();

                Kinematics k{
                    {       "m_mumu", re_mu },
                    { "cos(theta_l)",   0.4 }
                };
                auto pdf  = Observable::make("B2::P(m_mumu,cos(theta_l))@RvDV:2026A", p, k, Options{});
                auto re_c = Observable::make("B2::Re{c_" + index + "}(cos(theta_l))@RvDV:2026A", p, k, Options{});
                auto im_c = Observable::make("B2::Im{c_" + index + "}(cos(theta_l))@RvDV:2026A", p, k, Options{});

                const auto remainder = [&](const double & m) -> double
                {
                    k.set("m_mumu", m);
                    const complex<double> c(re_c->evaluate(), im_c->evaluate());
                    return pdf->evaluate() - (c / (m - complex<double>(re_mu, im_mu))).real();
                };

                // the PDF varies by more than a factor of 100 over +-10 gamma; the remainder, which keeps the other pole, by less than 5e-3
                const double gamma = -im_mu, at_pole = remainder(re_mu);
                for (const double offset : { -10.0, -1.0, 1.0, 10.0 })
                {
                    TEST_CHECK_RELATIVE_ERROR(remainder(re_mu + offset * gamma), at_pole, 5.0e-3);
                }
            }

            // the signal PDF is registered and normalized
            {
                Kinematics k{
                    {           "m_mumu",  1.0 },
                    {     "cos(theta_l)",  0.3 },
                    {       "m_mumu_min",  0.0 },
                    {       "m_mumu_max",  5.0 },
                    { "cos(theta_l)_min", -1.2 },
                    { "cos(theta_l)_max",  1.2 }
                };
                auto pdf = SignalPDF::make("benchmark::B2@RvDV:2026A", p, k, Options{});
                TEST_CHECK(pdf.get() != nullptr);
                TEST_CHECK_RELATIVE_ERROR(std::exp(pdf->evaluate() - pdf->normalization()), 2.670137119275402e+04 / 3.787096935763e+07, 1.0e-8);
            }
        }
} benchmark_b2_test;
