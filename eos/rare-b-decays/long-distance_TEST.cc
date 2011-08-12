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
#include <eos/rare-b-decays/charm-loops.hh>
#include <eos/rare-b-decays/long-distance.hh>

#include <cmath>
#include <vector>

#include <iostream>

using namespace test;
using namespace eos;

class LongDistanceTest :
    public TestCase
{
    public:
        LongDistanceTest() :
            TestCase("long_distance_test")
        {
        }

        virtual void run() const
        {
            static const double m_c = 1.2;
            static const double eps = 1e-5;
            static const std::vector<double> inputs{ 14.00, 15.00, 16.00, 19.21 };
            static const std::vector<complex<double>> results{
                complex<double>(1.13014, 0.381747),
                complex<double>(1.60574, 0.766175),
                complex<double>(0.57689, 2.113815),
                complex<double>(1.02058, 1.780107),
            };
            auto r = results.cbegin();
            for (auto i = inputs.cbegin() ; inputs.cend() != i ; ++i, ++r)
            {
                complex<double> g = LongDistance::g_had_ccbar(*i, m_c);

                TEST_CHECK_NEARLY_EQUAL(real(*r), real(g), eps);
                TEST_CHECK_NEARLY_EQUAL(imag(*r), imag(g), eps);
            }
        }
} long_distance_test;
