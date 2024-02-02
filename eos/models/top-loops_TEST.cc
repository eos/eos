/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
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
#include <eos/models/top-loops.hh>

#include <cmath>

using namespace test;
using namespace eos;

class TopLoopsTest :
    public TestCase
{
    public:
        TopLoopsTest() :
            TestCase("top_loops_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-9;

            /* One Loop */

            /* TopLoops::A0 */
            TEST_CHECK_NEARLY_EQUAL(TopLoops::A0(1.01), -1.068365257524, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::A0(1.2),  -1.049298669426, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::A0(2.0),  -0.989621888924, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::A0(4.0),  -0.906216652958, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::A0(5.0),  -0.880122806829, eps);

            /* TopLoops::B0 */
            TEST_CHECK_NEARLY_EQUAL(TopLoops::B0(1.01), 0.1245854042490, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::B0(1.2),  0.1174116759550, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::B0(2.0),  0.0965735902800, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::B0(4.0),  0.0706993734578, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::B0(5.0),  0.0632373369089, eps);

            /* TopLoops::C0 */
            TEST_CHECK_NEARLY_EQUAL(TopLoops::C0(1.01), 0.189582291687,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::C0(1.2),  0.228752692673,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::C0(2.0),  0.386294361120,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::C0(4.0),  0.744895614204,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::C0(5.0),  0.912517363726,  eps);

            /* TopLoops::D0 */
            TEST_CHECK_NEARLY_EQUAL(TopLoops::D0(1.01), 1.009965964272,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::D0(1.2),  0.942637752306,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::D0(2.0),  0.750780172847,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::D0(4.0),  0.511391014634,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::D0(5.0),  0.439870399317,  eps);

            /* TopLoops::E0 */
            TEST_CHECK_NEARLY_EQUAL(TopLoops::E0(1.01), 0.705188995053,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::E0(1.2),  0.652807696367,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::E0(2.0),  0.520248203698,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::E0(4.0),  0.390023196843,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::E0(5.0),  0.358633084449,  eps);

            /* TopLoops::F0 */
            TEST_CHECK_NEARLY_EQUAL(TopLoops::F0(1.01), -0.541169152595, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::F0(1.2),  -0.532564994828, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::F0(2.0),  -0.507783583307, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::F0(4.0),  -0.478135004113, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::F0(5.0),  -0.469971659962, eps);

            /* Two Loop */

            static const double mu_t = 120, m_t = 170.067834276559, log_t = log(mu_t / m_t);

            /* TopLoops::A1 */
            TEST_CHECK_NEARLY_EQUAL(TopLoops::A1(1.2, log_t), +9.1998495895399, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::A1(2.0, log_t), +8.8159805356153, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::A1(4.0, log_t), +8.3546426020653, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::A1(5.0, log_t), +8.2252286129223, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::A1(6.0, log_t), +8.1270411885495, eps);

            TEST_CHECK_NEARLY_EQUAL(TopLoops::B1(1.2, log_t), +0.7468955630420, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::B1(2.0, log_t), +0.6278464641750, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::B1(4.0, log_t), +0.4697994022820, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::B1(5.0, log_t), +0.4224942388700, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::B1(6.0, log_t), +0.3857046769660, eps);

            TEST_CHECK_NEARLY_EQUAL(TopLoops::C1(1.2, log_t), -0.3947275884940, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::C1(2.0, log_t), -1.2195006654023, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::C1(4.0, log_t), -3.2748716039000, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::C1(5.0, log_t), -4.2970357867133, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::C1(6.0, log_t), -5.3167141273327, eps);

            TEST_CHECK_NEARLY_EQUAL(TopLoops::D1(1.2, log_t), +1.6736393150700, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::D1(2.0, log_t), +1.2260072732383, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::D1(4.0, log_t), +0.7675637354830, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::D1(5.0, log_t), +0.6580815813680, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::D1(6.0, log_t), +0.5814206574330, eps);

            TEST_CHECK_NEARLY_EQUAL(TopLoops::E1(1.2, log_t), +1.0785040184548, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::E1(2.0, log_t), +0.8428693243740, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::E1(4.0, log_t), +0.6051834608610, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::E1(5.0, log_t), +0.5472498941840, eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::E1(6.0, log_t), +0.5058619620550, eps);

            TEST_CHECK_NEARLY_EQUAL(TopLoops::F1(1.2, log_t), +3.541936481881,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::F1(2.0, log_t), +3.439077304370,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::F1(4.0, log_t), +3.335560345606,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::F1(5.0, log_t), +3.310086286383,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::F1(6.0, log_t), +3.291755326435,  eps);

            TEST_CHECK_NEARLY_EQUAL(TopLoops::G1(1.2, log_t), +2.119194210716,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::G1(2.0, log_t), +1.714648791598,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::G1(4.0, log_t), +1.372579451930,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::G1(5.0, log_t), +1.299043502121,  eps);
            TEST_CHECK_NEARLY_EQUAL(TopLoops::G1(6.0, log_t), +1.248922424830,  eps);
        }
} top_loops_test;
