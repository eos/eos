/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Philip Lüghausen
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
#include <eos/form-factors/b-lcdas-flvd2022.hh>

#include <numeric>

using namespace test;
using namespace eos;
using namespace b_lcdas;

class FLvD2022Test :
    public TestCase
{
    public:
        FLvD2022Test() :
            TestCase("b_lcdas_flvd2022_test")
        {
        }

        virtual void run() const
        {
            // basic test
            {
                Parameters p = Parameters::Defaults();

                p["B_u::mu_0@FLvD2022"] = 1.0;
                p["B_u::omega_0@FLvD2022"] = 1.0;

                std::array<double, 9> ref = { 1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0 };
                for (size_t k = 0; k < ref.size(); k++)
                {
                    p["B_u::a^phi+_" + std::to_string(k) + "@FLvD2022"] = ref[k];
                }

                FLvD2022 blcdas(p, Options());

                auto [c, c_end] = blcdas.coefficient_range(1.0);
                for (auto it = c ; it != c_end ; ++it)
                {
                    TEST_CHECK_NEARLY_EQUAL(*it, ref[std::distance(c, it)], 1e-15);
                }
            }

            // RG evolution
            {
                Parameters p = Parameters::Defaults();

                p["B_u::mu_0@FLvD2022"] = 1.0;
                p["B_u::omega_0@FLvD2022"] = 0.55;
                const double mu = 1.5;

                // Evolved coefficients below are obtained using numerical evaluation of the exact 1-loop evolution
                // The test checks the implemented fast *approximation*
                std::array<double, 9> ref = { 1.6, 0.8, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0 };
                std::array<double, 9> res_evolved = { 1.4134123399029324, 0.5565559530817631, 0.2504760836874369, 0.12325426840513773, 0.09162847424850495, -0.019314494689795878, 0.015306298689121587, -0.013140355927689062, 0.011614395045637823 };
                for (size_t k = 0; k < ref.size(); k++)
                {
                    p["B_u::a^phi+_" + std::to_string(k) + "@FLvD2022"] = ref[k];
                }

                FLvD2022 blcdas(p, Options());

                auto [c, c_end] = blcdas.coefficient_range(mu);
                for (auto it = c ; it != c_end ; ++it)
                {
                    TEST_CHECK_NEARLY_EQUAL(*it, res_evolved[std::distance(c, it)], 1e-7);
                }
            }
        }
} b_lcdas_flvd2022_test;
