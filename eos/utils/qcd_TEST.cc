/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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
#include <eos/utils/qcd.hh>

#include <cmath>
#include <array>
#include <vector>

using namespace test;
using namespace eos;

class QCDMassesTest :
    public TestCase
{
    public:
        QCDMassesTest() :
            TestCase("qcd_masses_test")
        {
        }

        virtual void run() const
        {
            // m_q_msbar(m_q_pole) at 3-loop
            {
                static const double eps = 1e-7;
                static const std::vector<std::array<double, 5>> tests
                {
                                        // pole,  nf, alpha_s,   MSbar-ref
                    std::array<double, 5>{{ 172.0, 5.0,    0.10,  162.6620051 }},
                    std::array<double, 5>{{ 170.0, 5.0,    0.10,  160.7705865 }},
                    std::array<double, 5>{{ 168.0, 5.0,    0.10,  158.8791678 }},
                    std::array<double, 5>{{   4.9, 4.0,    0.22,    4.0271606 }},
                    std::array<double, 5>{{   4.8, 4.0,    0.22,    3.9449736 }},
                    std::array<double, 5>{{   4.7, 4.0,    0.22,    3.8627867 }},
                };

                for (const auto & test : tests)
                {
                    double m_q_pole_input = test[0], nf = test[1], alpha_s = test[2], m_q_msbar_ref = test[3];
                    double m_q_msbar = QCD::m_q_msbar(m_q_pole_input, alpha_s, nf);

                    TEST_CHECK_RELATIVE_ERROR(m_q_msbar, m_q_msbar_ref, eps);
                }
            }

            // m_q_pole(m_q_msbar) at orders 1 to 3
            {
                static const double eps = 1e-15; // exact to float precision
                TEST_CHECK_NEARLY_EQUAL(QCD::m_q_pole(4.0, 0.22, 5.0, 0), 4.0,               eps);
                TEST_CHECK_NEARLY_EQUAL(QCD::m_q_pole(4.0, 0.22, 5.0, 1), 4.373483599788981, eps);
                TEST_CHECK_NEARLY_EQUAL(QCD::m_q_pole(4.0, 0.22, 5.0, 2), 4.535117636490992, eps);
                TEST_CHECK_NEARLY_EQUAL(QCD::m_q_pole(4.0, 0.22, 5.0, 3), 4.636150134920265, eps);
            }

            // Safety check: higher-than-implemented loop order throws error
            {
                const constexpr unsigned int higher_loop_order = 4;
                TEST_CHECK_THROWS(InternalError, QCD::m_q_pole(4.0, 0.22, 5.0, higher_loop_order));
            }
        }
} qcd_masses_test;
