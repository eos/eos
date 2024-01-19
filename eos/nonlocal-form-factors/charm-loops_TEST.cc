/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2017 Danny van Dyk
 * Copyright (c) 2010, 2011 Christian Wacker
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
#include <eos/nonlocal-form-factors/charm-loops.hh>

#include <cmath>

using namespace test;
using namespace eos;

class OneLoopTest :
    public TestCase
{
    public:
        OneLoopTest() :
            TestCase("one_loop_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Christoph Bobeth's result from May 2010 */
            /* One-Loop h */
            {
                static const double mu = 4.2, s = 1.0, m_c = 1.4, m_b = 4.8, eps = 0.00001;
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, s)), +1.57192, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, s)), +1.39626, eps);

                TEST_CHECK_NEARLY_EQUAL(CharmLoops::h(mu, s, m_c), +0.58013, eps);

                TEST_CHECK_NEARLY_EQUAL(CharmLoops::h(mu, s, m_b), -0.55926, eps);
            }

            /* One-Loop h, m = m_c, as a function of s, in interval around s = 0 */
            {
                static const double mu = 4.2, m_c = 1.4, eps = 0.00001;

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, -9.0, m_c)), +0.24674, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, -8.0, m_c)), +0.27049, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, -7.0, m_c)), +0.29564, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, -6.0, m_c)), +0.32236, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, -5.0, m_c)), +0.35086, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, -4.0, m_c)), +0.38141, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, -3.0, m_c)), +0.41433, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, -2.0, m_c)), +0.45003, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, -1.0, m_c)), +0.48906, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu,  0.0, m_c)), +0.53210, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, +1.0, m_c)), +0.58014, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, +2.0, m_c)), +0.63454, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, +3.0, m_c)), +0.69738, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, +4.0, m_c)), +0.77198, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, +7.0, m_c)), +1.17626, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, +8.0, m_c)), +1.68173, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::h(mu, +9.0, m_c)), +1.48786, eps);

                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, -4.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, -3.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, -2.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, -1.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu,  0.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, +1.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, +2.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, +3.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, +4.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, +8.0, m_c)), +0.29422, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, +9.0, m_c)), +0.71961, eps);
            }

            /* One-Loop h, massless */
            {
                /* The imaginary part for massless h is always 4/9 Pi */
                static const double mu = 4.2, eps = 0.00001;

                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu,  1.00)), 4.0 / 9.0 * M_PI, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu,  6.00)), 4.0 / 9.0 * M_PI, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, 14.00)), 4.0 / 9.0 * M_PI, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::h(mu, 19.21)), 4.0 / 9.0 * M_PI, eps);
            }

            /* One-Loop B0, m = m_c */
            {
                static const double m_c = 1.4, eps = 0.00001;

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(-9.0, m_c)), -2.54698, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(-8.0, m_c)), -2.50056, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(-7.0, m_c)), -2.45159, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(-6.0, m_c)), -2.39975, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(-5.0, m_c)), -2.34468, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(-4.0, m_c)), -2.28592, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(-3.0, m_c)), -2.22288, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(-2.0, m_c)), -2.15488, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(-1.0, m_c)), -2.08099, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0( 0.0, m_c)), -2.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(+1.0, m_c)), -1.91028, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(+2.0, m_c)), -1.80952, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(+3.0, m_c)), -1.69427, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(+4.0, m_c)), -1.55906, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(+8.0, m_c)), -0.04026, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B0(+9.0, m_c)), -0.26980, eps);

                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0(-4.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0(-3.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0(-2.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0(-1.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0( 0.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0(+1.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0(+2.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0(+3.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0(+4.0, m_c)), +0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0(+8.0, m_c)), +0.44429, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B0(+9.0, m_c)), +1.12787, eps);
            }
        }
} one_loop_test;

class HelperTest :
    public TestCase
{
    public:
        HelperTest() :
            TestCase("helper_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Mathematica results from July 2010 */

            /* C0 */
            {
                static const double m_b = 4.45, s_one = m_b * m_b, eps = 0.000001;
                // real parts
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0( 0.0,  m_b)), -1.64493406685, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0( 0.5,  m_b)), -1.648607,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0( 1.0,  m_b)), -1.652304,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0( 2.0,  m_b)), -1.659779,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0( 3.0,  m_b)), -1.667360,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0( 6.0,  m_b)), -1.690774,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0( 9.0,  m_b)), -1.715257,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0(12.0,  m_b)), -1.740899,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0(15.0,  m_b)), -1.767803,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0(18.0,  m_b)), -1.796088,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0(19.21, m_b)), -1.807916,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0(19.80, m_b)), -1.813774,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0(s_one, m_b)), -1.813799364,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0(19.81, m_b)), -1.813874,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0(20.0,  m_b)), -1.815774,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0(21.0,  m_b)), -1.825884,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::C0(22.0,  m_b)), -1.836178,      eps);

                // imag parts
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::C0( 1.0, m_b)), +0.000000000,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::C0( 6.0, m_b)), +0.000000000,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::C0(11.0, m_b)), +0.000000000,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::C0(16.0, m_b)), +0.000000000,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::C0(19.0, m_b)), +0.000000000,  eps);
            }
        }
} helper_test;

class SeidelFunctionTest :
    public TestCase
{
    public:
        SeidelFunctionTest() :
            TestCase("seidel_function_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Mathematica results */
            {
                static const double eps = 1e-7;
                static const double mu  = 4.2;
                static const double m_b = 4.6;

                // CharmLoops::A
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.1   * m_b * m_b, m_b)), 0.9773295847097823, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.2   * m_b * m_b, m_b)), 0.9214363985136449, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.3   * m_b * m_b, m_b)), 0.8763703269301939, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.4   * m_b * m_b, m_b)), 0.8378238727298053, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.5   * m_b * m_b, m_b)), 0.8038113953761623, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.6   * m_b * m_b, m_b)), 0.7731991760882022, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.7   * m_b * m_b, m_b)), 0.7452609937428958, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.8   * m_b * m_b, m_b)), 0.7194961842495493, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.9   * m_b * m_b, m_b)), 0.6955413360449034, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.989 * m_b * m_b, m_b)), 0.6755202491234607, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.991 * m_b * m_b, m_b)), 0.6750830687385123, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 0.999 * m_b * m_b, m_b)), 0.6733396774972732, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 1.0   * m_b * m_b, m_b)), 0.6731223501151583, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 1.001 * m_b * m_b, m_b)), 0.6729051546645827, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 1.009 * m_b * m_b, m_b)), 0.6711723181761371, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 1.011 * m_b * m_b, m_b)), 0.6707404153377090, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 1.1   * m_b * m_b, m_b)), 0.6520262388101052, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 1.2   * m_b * m_b, m_b)), 0.6320834983833805, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 1.3   * m_b * m_b, m_b)), 0.6131565526952305, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 1.4   * m_b * m_b, m_b)), 0.595131879233407,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::A(mu, 1.5   * m_b * m_b, m_b)), 0.5779144700352841, eps);

                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.1   * m_b * m_b, m_b)), 0.6146283866916408, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.2   * m_b * m_b, m_b)), 0.6636082140569280, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.3   * m_b * m_b, m_b)), 0.6989293588461240, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.4   * m_b * m_b, m_b)), 0.7267965974570413, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.5   * m_b * m_b, m_b)), 0.7498451601160822, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.6   * m_b * m_b, m_b)), 0.7694868420677202, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.7   * m_b * m_b, m_b)), 0.7865792246366115, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.8   * m_b * m_b, m_b)), 0.8016869872542616, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.9   * m_b * m_b, m_b)), 0.8152036571681478, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.989 * m_b * m_b, m_b)), 0.8261288291878165, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.991 * m_b * m_b, m_b)), 0.8263637305615362, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 0.999 * m_b * m_b, m_b)), 0.8272989394799615, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 1.0   * m_b * m_b, m_b)), 0.8274153490936081, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 1.001 * m_b * m_b, m_b)), 0.8275316501089502, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 1.009 * m_b * m_b, m_b)), 0.8284581709044996, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 1.011 * m_b * m_b, m_b)), 0.8286887280252740, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 1.1   * m_b * m_b, m_b)), 0.8385370486103412, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 1.2   * m_b * m_b, m_b)), 0.8487346187930029, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 1.3   * m_b * m_b, m_b)), 0.8581388308167810, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 1.4   * m_b * m_b, m_b)), 0.8668546781765816, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, 1.5   * m_b * m_b, m_b)), 0.8749677701646030, eps);
            }

            {
                static const double eps = 5e-7;
                static const double mu  = 4.2;
                static const double m_b = 4.6;

                // CharmLoops::B
                // compare with results from Mathematica
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.1    * m_b * m_b, m_b)),  1.774836109825061, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.2    * m_b * m_b, m_b)),  0.660841461662527, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.3    * m_b * m_b, m_b)),  0.054907084216395, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.4    * m_b * m_b, m_b)), -0.349267491647551, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.5    * m_b * m_b, m_b)), -0.645893534420443, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.6    * m_b * m_b, m_b)), -0.876221241512994, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.7    * m_b * m_b, m_b)), -1.061898918819812, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.8    * m_b * m_b, m_b)), -1.215647742767820, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.9    * m_b * m_b, m_b)), -1.345549747764342, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.989  * m_b * m_b, m_b)), -1.445551305876865, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.991  * m_b * m_b, m_b)), -1.447653911549639, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 0.999  * m_b * m_b, m_b)), -1.456005423324896, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 1.0    * m_b * m_b, m_b)), -1.457042788460703, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 1.001  * m_b * m_b, m_b)), -1.458078703478798, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 1.009  * m_b * m_b, m_b)), -1.466314200418869, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 1.011  * m_b * m_b, m_b)), -1.468358794660238, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 1.1    * m_b * m_b, m_b)), -1.553950465199060, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 1.2    * m_b * m_b, m_b)), -1.639057020304387, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 1.3    * m_b * m_b, m_b)), -1.714448422091565, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 1.4    * m_b * m_b, m_b)), -1.781724905891100, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::B(mu, 1.5    * m_b * m_b, m_b)), -1.842138876544684, eps);

                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.1    * m_b * m_b, m_b)),  5.553840794326391, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.2    * m_b * m_b, m_b)),  4.820249145955325, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.3    * m_b * m_b, m_b)),  4.360108135584585, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.4    * m_b * m_b, m_b)),  4.017944023522331, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.5    * m_b * m_b, m_b)),  3.743067928416354, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.6    * m_b * m_b, m_b)),  3.512167767255114, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.7    * m_b * m_b, m_b)),  3.312463519105912, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.8    * m_b * m_b, m_b)),  3.136141749913582, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.9    * m_b * m_b, m_b)),  2.978054117771269, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.989  * m_b * m_b, m_b)),  2.849764927794023, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.991  * m_b * m_b, m_b)),  2.846999895193758, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 0.999  * m_b * m_b, m_b)),  2.835988381444127, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 1.0    * m_b * m_b, m_b)),  2.834617374526911, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 1.001  * m_b * m_b, m_b)),  2.833247567362475, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 1.009  * m_b * m_b, m_b)),  2.822332037234379, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 1.011  * m_b * m_b, m_b)),  2.819614998962907, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 1.1    * m_b * m_b, m_b)),  2.703230935969235, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 1.2    * m_b * m_b, m_b)),  2.581943911053763, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 1.3    * m_b * m_b, m_b)),  2.469253015913550, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 1.4    * m_b * m_b, m_b)),  2.363973892039418, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::B(mu, 1.5    * m_b * m_b, m_b)),  2.265155869018138, eps);
            }
        }
} seidel_function_test;


class FormFactorsTest :
    public TestCase
{
    public:
        FormFactorsTest() :
            TestCase("form_factors_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Christoph Bobeth's result from May 2010 */

            /* Formfactors, massless loops */
            {
                static const double mu = 4.2, s = 6.0, m_b = 4.6, eps = 0.0000001;
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F17_massless(mu, s, m_b)), - 0.8832611, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F17_massless(mu, s, m_b)), - 0.6937322, eps);

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F27_massless(mu, s, m_b)), + 5.2995666, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F27_massless(mu, s, m_b)), + 4.1623936, eps);

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F19_massless(mu, s, m_b)), + 3.3632062, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F19_massless(mu, s, m_b)), - 6.9078480, eps);

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F29_massless(mu, s, m_b)), + 3.4455298, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F29_massless(mu, s, m_b)), +24.6919276, eps);

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F87_massless(mu, s, m_b)), - 1.2486221, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F87_massless(mu, s, m_b)), - 2.7925269, eps);

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F89_massless(s, m_b)), - 3.2730189, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F89_massless(s, m_b)),   0.0000000, eps);
            }

            /* Formfactors for O_8 are problematic near the zero recoil point */
            {
                static const double mu = 4.2, s = 19.2, m_b = 4.6, eps = 0.0000001;

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F87_massless(mu, s, m_b)), - 0.9708796,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F87_massless(mu, s, m_b)), - 2.7925268,  eps);

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F89_massless(s, m_b)), - 2.0208146,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F89_massless(s, m_b)),   0.0000000,  eps);
            }

            /* Check imaginary parts of the form factors at the boundaries of the Low Recoil region */
            {
                static const double mu = 4.2, s1 = 14.0, s2 = 19.2, m_b = 4.6, eps = 0.0000001;

                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, s1, m_b)), +0.7802809, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::A(mu, s2, m_b)), +0.8161455, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F87_massless(mu, s1, m_b)), -2.7925268, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F87_massless(mu, s2, m_b)), -2.7925268, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F19_massless(mu, s1, m_b)), -5.8682288, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F19_massless(mu, s2, m_b)), -5.4492617, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F29_massless(mu, s1, m_b)), 18.4542117, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F29_massless(mu, s2, m_b)), 15.9404096, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F89_massless(s1, m_b)),  0.0000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F89_massless(s2, m_b)),  0.0000000, eps);
            }

            /* Check F87_massless and F89_massless near s_hat = 1 */
            {
                static const double mu = 4.2, m_b = 4.6, eps = 6e-7;

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F87_massless(mu, 0.989 * m_b * m_b, m_b)), -0.951276751910152, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F87_massless(mu, 0.991 * m_b * m_b, m_b)), -0.950828012015162, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F87_massless(mu, 0.999 * m_b * m_b, m_b)), -0.949047612614735, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F87_massless(mu, 1.0   * m_b * m_b, m_b)), -0.948826693927197, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F87_massless(mu, 1.001 * m_b * m_b, m_b)), -0.948606136377817, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F87_massless(mu, 1.009 * m_b * m_b, m_b)), -0.946854628625707, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F87_massless(mu, 1.011 * m_b * m_b, m_b)), -0.946420334996110, eps);

                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F87_massless(mu, 0.989 * m_b * m_b, m_b)), -2.792526803190927, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F87_massless(mu, 0.991 * m_b * m_b, m_b)), -2.792526803190927, eps);

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F89_massless(0.989 * m_b * m_b, m_b)), -1.944265886425612, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F89_massless(0.991 * m_b * m_b, m_b)), -1.942498970870095, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F89_massless(0.999 * m_b * m_b, m_b)), -1.935479699592446, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F89_massless(1.0   * m_b * m_b, m_b)), -1.934607689969977, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F89_massless(1.001 * m_b * m_b, m_b)), -1.933736871247196, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F89_massless(1.009 * m_b * m_b, m_b)), -1.926812876524436, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F89_massless(1.011 * m_b * m_b, m_b)), -1.925093602661525, eps);
            }

            /* Formfactors, massive loops at timelike s/q^2 */
            {
                static const double mu = 4.2, s = 6.0, m_b = 4.6, m_c = 1.2, eps = 1e-7;
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F17_massive(mu, s, m_b, m_c)), - 0.73093991, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F17_massive(mu, s, m_b, m_c)), - 0.17771334, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F27_massive(mu, s, m_b, m_c)), + 4.38563254, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F27_massive(mu, s, m_b, m_c)), + 1.06627403, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F19_massive(mu, s, m_b, m_c)), -34.40870331, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F19_massive(mu, s, m_b, m_c)), - 0.25864665, eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F29_massive(mu, s, m_b, m_c)), + 6.27364439, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F29_massive(mu, s, m_b, m_c)), + 1.55195807, eps);
            }

            /* Formfactors, massive loops at spacelike s/q^2 */
            {
                static const double mu = 4.2, m_b = 4.6, m_c = 1.2, eps = 1e-5;

                TEST_CHECK_RELATIVE_ERROR(real(CharmLoops::F17_massive(mu, -6.0, m_b, m_c)), - 0.5851990,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(CharmLoops::F17_massive(mu, -6.0, m_b, m_c)), - 0.0622661,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(CharmLoops::F17_massive(mu, -1.0, m_b, m_c)), - 0.6507510,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(CharmLoops::F17_massive(mu, -1.0, m_b, m_c)), - 0.0921000,  eps);

                TEST_CHECK_RELATIVE_ERROR(real(CharmLoops::F27_massive(mu, -6.0, m_b, m_c)), + 3.5112500,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(CharmLoops::F27_massive(mu, -6.0, m_b, m_c)), + 0.3736050,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(CharmLoops::F27_massive(mu, -1.0, m_b, m_c)), + 3.9045200,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(CharmLoops::F27_massive(mu, -1.0, m_b, m_c)), + 0.5526040,  eps);

                TEST_CHECK_RELATIVE_ERROR(real(CharmLoops::F19_massive(mu, -6.0, m_b, m_c)), - 3.2450800,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(CharmLoops::F19_massive(mu, -6.0, m_b, m_c)), + 0.1208170,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(CharmLoops::F19_massive(mu, -1.0, m_b, m_c)), -10.1066000,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(CharmLoops::F19_massive(mu, -1.0, m_b, m_c)), + 0.1100320,  eps);

                TEST_CHECK_RELATIVE_ERROR(real(CharmLoops::F29_massive(mu, -6.0, m_b, m_c)), + 4.4729700,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(CharmLoops::F29_massive(mu, -6.0, m_b, m_c)), - 0.7247960,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(CharmLoops::F29_massive(mu, -1.0, m_b, m_c)), + 4.0282600,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(CharmLoops::F29_massive(mu, -1.0, m_b, m_c)), - 0.6601020,  eps);
            }
        }
} two_loop_test;

class LowRecoilTest :
    public TestCase
{
    public:
        LowRecoilTest() :
            TestCase("low_recoil_test")
        {
        }

        virtual void run() const
        {
            using std::imag;
            using std::real;

            /* Comparison with Stefan Meinel from September 2019 */

            Parameters p  = Parameters::Defaults();
            auto       m  = Model::make("WET", p, Options());

            {
                static const double mu = 4.2, s = 15.0, eps = 1e-2;
                const double        alpha_s = m->alpha_s(mu);
                const double        m_b_PS = m->m_b_ps(2.0);
                const auto          m_c = m->m_c_msbar(mu);
                const auto          wc = m->wilson_coefficients_b_to_s(mu, LeptonFlavor::muon, false);

                TEST_CHECK_NEARLY_EQUAL(m_b_PS,          4.46, eps);

                const auto c7eff = ShortDistanceLowRecoil::c7eff(s, mu, alpha_s, m_b_PS, true, wc);
                TEST_CHECK_NEARLY_EQUAL(real(c7eff),    -0.39, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(c7eff),    -0.10, eps);

                const auto c9eff = ShortDistanceLowRecoil::c9eff(s, mu, alpha_s, m_b_PS, m_c, true, false, 0.0, wc);
                TEST_CHECK_NEARLY_EQUAL(real(c9eff),    +4.66, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(c9eff),     0.55, eps);
            }
        }
} low_recoil_test;


class CharmlessTest :
    public TestCase
{
    public:
        CharmlessTest() :
            TestCase("charmless_test")
        {
        }

        virtual void run() const
        {
            using std::imag;
            using std::real;

            /* Test cubic spline interpolation by comparison with Javier Virto's notebook */
            {

                static const double eps = 1e-3;

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F17_massive_Qsb(2.0)), -0.0715022,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F17_massive_Qsb(2.0)), -0.00894292,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F19_massive_Qsb(2.0)),  0.195383,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F19_massive_Qsb(2.0)),  0.25814,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F27_massive_Qsb(2.0)),  0.429013,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F27_massive_Qsb(2.0)),  0.0536575,    eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F29_massive_Qsb(2.0)), -1.1723,       eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F29_massive_Qsb(2.0)), -1.54884,      eps);

                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F17_massive_Qsb(8.0)), -0.0683473,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F17_massive_Qsb(8.0)), -0.0153357,    eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F19_massive_Qsb(8.0)),  0.102639,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F19_massive_Qsb(8.0)),  0.191314,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F27_massive_Qsb(8.0)),  0.410084,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F27_massive_Qsb(8.0)),  0.0920144,    eps);
                TEST_CHECK_NEARLY_EQUAL(real(CharmLoops::F29_massive_Qsb(8.0)), -0.615832,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(CharmLoops::F29_massive_Qsb(8.0)), -1.14788,      eps);

            }
        }
} charmless_test;
