
/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Danny van Dyk
 * Copyright (c) 2026 Dominik Suelmann
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

#include <eos/utils/rge-impl.hh>

#include <test/test.hh>

#include <array>

using namespace test;
using namespace eos;

class MultiplicativeRGELLTest : public TestCase
{
    public:
        MultiplicativeRGELLTest() :
            TestCase("multiplicative_rge_ll_test")
        {
        }

        virtual void
        run() const
        {
            // trivial test case (nf = 5, dim = 10)
            {
                const std::array<double, 10u>                  gamma_0_ev{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
                const std::array<std::array<double, 10u>, 10u> V{
                    { { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } }
                };
                const MultiplicativeRenormalizationGroupEvolution<accuracy::LL, 5u, 10u> rge(gamma_0_ev, V);
                const std::array<double, 10u>                                            c_0_0{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

                const double            alpha_s_mu = 0.218017;
                const double            alpha_s_0  = 0.121864;
                std::array<double, 10u> result     = rge.evolve(alpha_s_mu, alpha_s_0, c_0_0);

                static const double eps = 1.0e-15;
                TEST_CHECK_NEARLY_EQUAL(result[0], 0.1, eps);
                TEST_CHECK_NEARLY_EQUAL(result[1], 0.2, eps);
                TEST_CHECK_NEARLY_EQUAL(result[2], 0.3, eps);
                TEST_CHECK_NEARLY_EQUAL(result[3], 0.4, eps);
                TEST_CHECK_NEARLY_EQUAL(result[4], 0.5, eps);
                TEST_CHECK_NEARLY_EQUAL(result[5], 0.6, eps);
                TEST_CHECK_NEARLY_EQUAL(result[6], 0.7, eps);
                TEST_CHECK_NEARLY_EQUAL(result[7], 0.8, eps);
                TEST_CHECK_NEARLY_EQUAL(result[8], 0.9, eps);
                TEST_CHECK_NEARLY_EQUAL(result[9], 1.0, eps);
            }

            // current-current test case (nf = 5, dim = 2), checked by S. Meiser 2023/07/10
            {
                const double                                 sq2 = std::sqrt(2.0);
                const std::array<double, 2u>                 gamma_0_ev{ -8.0, +4.0 };
                const std::array<std::array<double, 2u>, 2u> V{
                    {
                     { -1.0 / sq2, +1.0 / sq2 },
                     { +1.0 / sq2, +1.0 / sq2 },
                     }
                };
                const MultiplicativeRenormalizationGroupEvolution<accuracy::LL, 5u, 2u> rge(gamma_0_ev, V);
                const std::array<double, 2u>                                            c_0_0{ 0.0, 1.0 };

                const double           alpha_s_mu = 0.218017;
                const double           alpha_s_0  = 0.121864;
                std::array<double, 2u> result     = rge.evolve(alpha_s_mu, alpha_s_0, c_0_0);

                static const double eps = 1.0e-6;
                TEST_CHECK_NEARLY_EQUAL(result[0], -0.247675, eps);
                TEST_CHECK_NEARLY_EQUAL(result[1], +1.106887, eps);
            }

            // non-symmetric test case (nf = 5, dim = 2), checked by S. Meiser 2023/07/10
            {
                const std::array<double, 2u>                 gamma_0_ev{ -16.0, +2.0 };
                const std::array<std::array<double, 2u>, 2u> V{
                    {
                     { 1.0 / 6.0, -4.0 / 3.0 },
                     { 1.0, 1.0 },
                     }
                };
                const MultiplicativeRenormalizationGroupEvolution<accuracy::LL, 5u, 2u> rge(gamma_0_ev, V);
                const std::array<double, 2u>                                            c_0_0{ 0.0, 1.0 };

                const double           alpha_s_mu = 0.218017;
                const double           alpha_s_0  = 0.121864;
                std::array<double, 2u> result     = rge.evolve(alpha_s_mu, alpha_s_0, c_0_0);

                static const double eps = 1.0e-6;
                TEST_CHECK_NEARLY_EQUAL(result[0], +0.134504, eps);
                TEST_CHECK_NEARLY_EQUAL(result[1], +1.733962, eps);
            }
        }
} multiplicative_rge_ll_test;

class MultiplicativeRGENLLTest : public TestCase
{
    public:
        MultiplicativeRGENLLTest() :
            TestCase("multiplicative_rge_nll_test")
        {
        }

        virtual void
        run() const
        {
            // trivial test case (nf = 5, dim = 10)
            {
                const std::array<double, 10u>                  gamma_0_ev{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
                const std::array<std::array<double, 10u>, 10u> V{
                    { { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } }
                };
                const std::array<std::array<double, 10u>, 10u> gamma_1{
                    { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } }
                };
                const MultiplicativeRenormalizationGroupEvolution<accuracy::NLL, 5u, 10u> rge(gamma_0_ev, V, gamma_1);
                const std::array<double, 10u>                                             c_0_0{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
                const std::array<double, 10u>                                             c_0_1{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

                const double            alpha_s_mu = 0.3 * 4.0 * M_PI;
                const double            alpha_s_0  = 0.1 * 4.0 * M_PI;
                std::array<double, 10u> result     = rge.evolve(alpha_s_mu, alpha_s_0, c_0_0, c_0_1);

                static const double eps = 1.0e-15;
                TEST_CHECK_NEARLY_EQUAL(result[0], 0.11, eps);
                TEST_CHECK_NEARLY_EQUAL(result[1], 0.22, eps);
                TEST_CHECK_NEARLY_EQUAL(result[2], 0.33, eps);
                TEST_CHECK_NEARLY_EQUAL(result[3], 0.44, eps);
                TEST_CHECK_NEARLY_EQUAL(result[4], 0.55, eps);
                TEST_CHECK_NEARLY_EQUAL(result[5], 0.66, eps);
                TEST_CHECK_NEARLY_EQUAL(result[6], 0.77, eps);
                TEST_CHECK_NEARLY_EQUAL(result[7], 0.88, eps);
                TEST_CHECK_NEARLY_EQUAL(result[8], 0.99, eps);
                TEST_CHECK_NEARLY_EQUAL(result[9], 1.10, eps);
            }

            // current-current test case (nf = 5, dim = 2)
            {
                const double                                 sq2 = std::sqrt(2.0);
                const std::array<double, 2u>                 gamma_0_ev{ -8.0, +4.0 };
                const std::array<std::array<double, 2u>, 2u> V{
                    {
                     { -1.0 / sq2, +1.0 / sq2 },
                     { +1.0 / sq2, +1.0 / sq2 },
                     }
                };
                const std::array<std::array<double, 2u>, 2u> gamma_1{
                    { { -209.0 / 18.0, +41.0 / 6.0 }, { +41.0 / 6.0, -209.0 / 18.0 } }
                };
                const MultiplicativeRenormalizationGroupEvolution<accuracy::NLL, 5u, 2u> rge(gamma_0_ev, V, gamma_1);
                const std::array<double, 2u>                                             c_0_0{ 0.0, 1.0 };
                const std::array<double, 2u>                                             c_0_1{ 11.0 / 2.0, -11.0 / 6.0 };

                const double           alpha_s_mu = 0.218017;
                const double           alpha_s_0  = 0.121864;
                std::array<double, 2u> result     = rge.evolve(alpha_s_mu, alpha_s_0, c_0_0, c_0_1);

                static const double eps = 1.0e-6;
                TEST_CHECK_NEARLY_EQUAL(result[0], -0.172203, eps);
                TEST_CHECK_NEARLY_EQUAL(result[1], +1.073145, eps);
            }

            // non-symmetric test case (nf = 5, dim = 2)
            {
                const std::array<double, 2u>                 gamma_0_ev{ -16.0, +2.0 };
                const std::array<std::array<double, 2u>, 2u> V{
                    {
                     { 1.0 / 6.0, -4.0 / 3.0 },
                     { 1.0, 1.0 },
                     }
                };
                const std::array<std::array<double, 2u>, 2u> gamma_1{
                    {
                     { -28.0 / 3.0, -374.0 / 3.0 },
                     { -2044.0 / 27.0, -2975.0 / 18.0 },
                     }
                };
                const MultiplicativeRenormalizationGroupEvolution<accuracy::NLL, 5u, 2u> rge(gamma_0_ev, V, gamma_1);
                const std::array<double, 2u>                                             c_0_0{ 0.0, 1.0 };
                const std::array<double, 2u>                                             c_0_1{ 11.0 / 2.0, -11.0 / 6.0 };

                const double           alpha_s_mu = 0.218017;
                const double           alpha_s_0  = 0.121864;
                std::array<double, 2u> result     = rge.evolve(alpha_s_mu, alpha_s_0, c_0_0, c_0_1);

                static const double eps = 1.0e-6;
                TEST_CHECK_NEARLY_EQUAL(result[0], +0.229589, eps);
                TEST_CHECK_NEARLY_EQUAL(result[1], +1.826340, eps);
            }

            // test case evolve_intermediate (nf = 5, dim = 2)
            {
                const std::array<double, 2u>                 gamma_0_ev{ -16.0, +2.0 };
                const std::array<std::array<double, 2u>, 2u> V{
                    {
                     { 1.0 / 6.0, -4.0 / 3.0 },
                     { 1.0, 1.0 },
                     }
                };
                const std::array<std::array<double, 2u>, 2u> gamma_1{
                    {
                     { -28.0 / 3.0, -374.0 / 3.0 },
                     { -2044.0 / 27.0, -2975.0 / 18.0 },
                     }
                };
                const MultiplicativeRenormalizationGroupEvolution<accuracy::NLL, 5u, 2u> rge(gamma_0_ev, V, gamma_1);
                const std::array<double, 2u>                                             c_0_0{ 0.0, 1.0 };
                const std::array<double, 2u>                                             c_0_1{ 11.0 / 2.0, -11.0 / 6.0 };

                const double                                                     alpha_s_mu = 0.121864;
                const double                                                     alpha_s_0  = 0.121864;
                const std::tuple<std::array<double, 2u>, std::array<double, 2u>> result_int = rge.evolve_intermediate(alpha_s_mu, alpha_s_0, c_0_0, c_0_1);
                const std::array<double, 2u>                                     result     = rge.evolve(alpha_s_mu, alpha_s_0, c_0_0, c_0_1);

                static const double eps = 1.0e-6;
                static const double as  = (alpha_s_mu / (4 * M_PI));
                TEST_CHECK_NEARLY_EQUAL(std::get<0>(result_int)[0], 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::get<0>(result_int)[1], 1.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::get<1>(result_int)[0], 11.0 / 2.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::get<1>(result_int)[1], -11.0 / 6.0, eps);
                TEST_CHECK_NEARLY_EQUAL(result[0], std::get<0>(result_int)[0] + as * std::get<1>(result_int)[0], std::abs(10 * as * as * std::get<1>(result_int)[0]));
                TEST_CHECK_NEARLY_EQUAL(result[1], std::get<0>(result_int)[1] + as * std::get<1>(result_int)[1], std::abs(10 * as * as * std::get<1>(result_int)[1]));

                const double                                                     alpha_s_mu2 = 0.218017;
                const std::tuple<std::array<double, 2u>, std::array<double, 2u>> result_int2 = rge.evolve_intermediate(alpha_s_mu2, alpha_s_0, c_0_0, c_0_1);
                const std::array<double, 2u>                                     result2     = rge.evolve(alpha_s_mu2, alpha_s_0, c_0_0, c_0_1);
                static const double                                              as2         = (alpha_s_mu2 / (4 * M_PI));

                TEST_CHECK_NEARLY_EQUAL(result2[0], std::get<0>(result_int2)[0] + as2 * std::get<1>(result_int2)[0], std::abs(10 * as2 * as2 * std::get<1>(result_int2)[0]));
                TEST_CHECK_NEARLY_EQUAL(result2[1], std::get<0>(result_int2)[1] + as2 * std::get<1>(result_int2)[1], std::abs(10 * as2 * as2 * std::get<1>(result_int2)[1]));
            }
        }
} multiplicative_rge_nll_test;

class MultiplicativeRGENNLLTest : public TestCase
{
    public:
        MultiplicativeRGENNLLTest() :
            TestCase("multiplicative_rge_nnll_test")
        {
        }

        virtual void
        run() const
        {
            // trivial test case (nf = 5, dim = 10)
            {
                const std::array<double, 10u>                  gamma_0_ev{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
                const std::array<std::array<double, 10u>, 10u> V{
                    { { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } }
                };
                const std::array<std::array<double, 10u>, 10u> gamma_1{
                    { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } }
                };
                const std::array<std::array<double, 10u>, 10u> gamma_2{
                    { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                     { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } }
                };
                const MultiplicativeRenormalizationGroupEvolution<accuracy::NNLL, 5u, 10u> rge(gamma_0_ev, V, gamma_1, gamma_2);
                const std::array<double, 10u>                                              c_0_0{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
                const std::array<double, 10u>                                              c_0_1{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
                const std::array<double, 10u>                                              c_0_2{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

                const double            alpha_s_mu = 0.3 * 4.0 * M_PI;
                const double            alpha_s_0  = 0.1 * 4.0 * M_PI;
                std::array<double, 10u> result     = rge.evolve(alpha_s_mu, alpha_s_0, c_0_0, c_0_1, c_0_2);

                static const double eps = 1.0e-15;
                TEST_CHECK_NEARLY_EQUAL(result[0], 0.111, eps);
                TEST_CHECK_NEARLY_EQUAL(result[1], 0.222, eps);
                TEST_CHECK_NEARLY_EQUAL(result[2], 0.333, eps);
                TEST_CHECK_NEARLY_EQUAL(result[3], 0.444, eps);
                TEST_CHECK_NEARLY_EQUAL(result[4], 0.555, eps);
                TEST_CHECK_NEARLY_EQUAL(result[5], 0.666, eps);
                TEST_CHECK_NEARLY_EQUAL(result[6], 0.777, eps);
                TEST_CHECK_NEARLY_EQUAL(result[7], 0.888, eps);
                TEST_CHECK_NEARLY_EQUAL(result[8], 0.999, eps);
                TEST_CHECK_NEARLY_EQUAL(result[9], 1.11, eps);
            }

            // current-current test case (nf = 5, dim = 2)
            {
                const double                                 sq2 = std::sqrt(2.0);
                const std::array<double, 2u>                 gamma_0_ev{ -8.0, +4.0 };
                const std::array<std::array<double, 2u>, 2u> V{
                    {
                     { -1.0 / sq2, +1.0 / sq2 },
                     { +1.0 / sq2, +1.0 / sq2 },
                     }
                };
                const std::array<std::array<double, 2u>, 2u> gamma_1{
                    { { -209.0 / 18.0, +41.0 / 6.0 }, { +41.0 / 6.0, -209.0 / 18.0 } }
                };
                const std::array<std::array<double, 2u>, 2u> gamma_2{
                    { { +239.0 / 18.0, -43.0 / 6.0 }, { -43.0 / 6.0, +239.0 / 18.0 } }
                };
                const MultiplicativeRenormalizationGroupEvolution<accuracy::NNLL, 5u, 2u> rge(gamma_0_ev, V, gamma_1, gamma_2);
                const std::array<double, 2u>                                              c_0_0{ 0.0, 1.0 };
                const std::array<double, 2u>                                              c_0_1{ 11.0 / 2.0, -11.0 / 6.0 };
                const std::array<double, 2u>                                              c_0_2{ -9.0 / 2.0, -13.0 / 6.0 };

                const double           alpha_s_mu = 0.218017;
                const double           alpha_s_0  = 0.121864;
                std::array<double, 2u> result     = rge.evolve(alpha_s_mu, alpha_s_0, c_0_0, c_0_1, c_0_2);

                static const double eps = 1.0e-8;
                TEST_CHECK_NEARLY_EQUAL(result[0], -0.171944007908276, eps);
                TEST_CHECK_NEARLY_EQUAL(result[1], +1.072074409306490, eps);
            }

            // non-symmetric test case (nf = 5, dim = 2)
            {
                const std::array<double, 2u>                 gamma_0_ev{ -16.0, +2.0 };
                const std::array<std::array<double, 2u>, 2u> V{
                    {
                     { 1.0 / 6.0, -4.0 / 3.0 },
                     { 1.0, 1.0 },
                     }
                };
                const std::array<std::array<double, 2u>, 2u> gamma_1{
                    {
                     { -28.0 / 3.0, -374.0 / 3.0 },
                     { -2044.0 / 27.0, -2975.0 / 18.0 },
                     }
                };
                const std::array<std::array<double, 2u>, 2u> gamma_2{
                    {
                     { -25.0 / 3.0, -377.0 / 3.0 },
                     { -2047.0 / 27.0, -2978.0 / 18.0 },
                     }
                };
                const MultiplicativeRenormalizationGroupEvolution<accuracy::NNLL, 5u, 2u> rge(gamma_0_ev, V, gamma_1, gamma_2);
                const std::array<double, 2u>                                              c_0_0{ 0.0, 1.0 };
                const std::array<double, 2u>                                              c_0_1{ 11.0 / 2.0, -11.0 / 6.0 };
                const std::array<double, 2u>                                              c_0_2{ -9.0 / 2.0, -13.0 / 6.0 };

                const double           alpha_s_mu = 0.218017;
                const double           alpha_s_0  = 0.121864;
                std::array<double, 2u> result     = rge.evolve(alpha_s_mu, alpha_s_0, c_0_0, c_0_1, c_0_2);

                static const double eps = 1.0e-8;
                TEST_CHECK_NEARLY_EQUAL(result[0], +0.23134566937815798, eps);
                TEST_CHECK_NEARLY_EQUAL(result[1], +1.82581433200249621, eps);
            }
        }
} multiplicative_rge_nnll_test;
