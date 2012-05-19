/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
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
#include <eos/utils/welford.hh>

#include <vector>

using namespace test;
using namespace eos;

class WelfordTest :
    public TestCase
{
    public:
        WelfordTest() :
            TestCase("welford_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-14;

            // compare with numpy.average and numpy.var results
            {
                std::vector<double> samples { 1.23, 413.132, 213.12 };

                Welford w;
                for (auto s = samples.cbegin() ; s != samples.cend() ; ++s)
                {
                    w.add(*s);
                }

                TEST_CHECK_EQUAL(w.number_of_elements(), 3u);
                TEST_CHECK_RELATIVE_ERROR(w.mean(), 209.16066666666665697, eps);
                TEST_CHECK_RELATIVE_ERROR(w.variance(), 42427.57164133333571954, eps);
            }
        }
} welford_test;
