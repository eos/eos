/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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
#include <src/utils/qcd.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace eos;

class QCDTest :
    public TestCase
{
    public:
        QCDTest() :
            TestCase("qcd_test")
        {
        }

        virtual void run() const
        {
#if 0
            static const double eps = 1e-6;
            TEST_CHECK(std::abs(0.224302 - QCD::alpha_s( 4.2)) <= eps);
            TEST_CHECK(std::abs(0.215643 - QCD::alpha_s( 4.8)) <= eps);
            TEST_CHECK(std::abs(0.213133 - QCD::alpha_s( 5.0)) <= eps);
            TEST_CHECK(std::abs(0.120246 - QCD::alpha_s(80.0)) <= eps);
#endif
        }
} qcd_test;

class QCDBMassesTest :
    public TestCase
{
    public:
        QCDBMassesTest() :
            TestCase("qcd_b_masses_test")
        {
        }

        virtual void run() const
        {
#if 0
            static const double eps = 1e-6;

            /* b quark pole mass */
            TEST_CHECK(std::abs(4.780082 - QCD::mb_pole(4.1)) <= eps);
            TEST_CHECK(std::abs(4.888672 - QCD::mb_pole(4.2)) <= eps);
            TEST_CHECK(std::abs(4.997244 - QCD::mb_pole(4.3)) <= eps);
#endif
        }
} qcd_b_masses_test;
