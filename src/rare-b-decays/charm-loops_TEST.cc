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
#include <src/rare-b-decays/charm-loops.hh>

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

            /* One-Loop */
            {
                static const double mu = 4.2, s = 1.0, m_c = 1.4, m_b = 4.8, eps = 0.00001;
                TEST_CHECK_NEARLY_EQUAL(+1.57192, real(CharmLoops::h(mu, s)), eps);
                TEST_CHECK_NEARLY_EQUAL(+1.39626, imag(CharmLoops::h(mu, s)), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.58013, CharmLoops::h(mu, s, m_c), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.55926, CharmLoops::h(mu, s, m_b), eps);
            }

            /* One-Loop, massless */
            {
                /* The imaginary part for massless h is always 4/9 Pi */
                static const double mu = 4.2, eps = 0.00001;

                TEST_CHECK_NEARLY_EQUAL(4.0 / 9.0 * M_PI, imag(CharmLoops::h(mu,  1.00)), eps);
                TEST_CHECK_NEARLY_EQUAL(4.0 / 9.0 * M_PI, imag(CharmLoops::h(mu,  6.00)), eps);
                TEST_CHECK_NEARLY_EQUAL(4.0 / 9.0 * M_PI, imag(CharmLoops::h(mu, 14.00)), eps);
                TEST_CHECK_NEARLY_EQUAL(4.0 / 9.0 * M_PI, imag(CharmLoops::h(mu, 19.21)), eps);
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
                static const double m_b = 4.45, eps = 0.000001;
                // real parts
                TEST_CHECK_NEARLY_EQUAL(-1.64493406685, real(CharmLoops::C0( 0.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.648607,      real(CharmLoops::C0( 0.5,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.652304,      real(CharmLoops::C0( 1.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.659779,      real(CharmLoops::C0( 2.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.667360,      real(CharmLoops::C0( 3.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.690774,      real(CharmLoops::C0( 6.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.715257,      real(CharmLoops::C0( 9.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.740899,      real(CharmLoops::C0(12.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.767803,      real(CharmLoops::C0(15.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.796088,      real(CharmLoops::C0(18.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.807916,      real(CharmLoops::C0(19.21, m_b)), eps);

                // imag parts
                TEST_CHECK_NEARLY_EQUAL(+0.000000000,  imag(CharmLoops::C0( 1.0, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000,  imag(CharmLoops::C0( 6.0, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000,  imag(CharmLoops::C0(11.0, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000,  imag(CharmLoops::C0(16.0, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000,  imag(CharmLoops::C0(19.0, m_b)), eps);
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
                static const double mu = 4.2;
                static const double m_b = 4.6;

                // CharmLoops::A
                TEST_CHECK_NEARLY_EQUAL(0.9773295847097823, real(CharmLoops::A(mu, 0.1 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.9214363985136449, real(CharmLoops::A(mu, 0.2 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8763703269301939, real(CharmLoops::A(mu, 0.3 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8378238727298053, real(CharmLoops::A(mu, 0.4 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8038113953761623, real(CharmLoops::A(mu, 0.5 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.7731991760882022, real(CharmLoops::A(mu, 0.6 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.7452609937428958, real(CharmLoops::A(mu, 0.7 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.7194961842495493, real(CharmLoops::A(mu, 0.8 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.6955413360449034, real(CharmLoops::A(mu, 0.9 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.6731223501151583, real(CharmLoops::A(mu, 1.0 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.6520262388101052, real(CharmLoops::A(mu, 1.1 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.6320834983833805, real(CharmLoops::A(mu, 1.2 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.6131565526952305, real(CharmLoops::A(mu, 1.3 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.595131879233407,  real(CharmLoops::A(mu, 1.4 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.5779144700352841, real(CharmLoops::A(mu, 1.5 * m_b * m_b, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(0.6146283866916408, imag(CharmLoops::A(mu, 0.1 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.6636082140569280, imag(CharmLoops::A(mu, 0.2 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.6989293588461240, imag(CharmLoops::A(mu, 0.3 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.7267965974570413, imag(CharmLoops::A(mu, 0.4 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.7498451601160822, imag(CharmLoops::A(mu, 0.5 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.7694868420677202, imag(CharmLoops::A(mu, 0.6 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.7865792246366115, imag(CharmLoops::A(mu, 0.7 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8016869872542616, imag(CharmLoops::A(mu, 0.8 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8152036571681478, imag(CharmLoops::A(mu, 0.9 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8274153490936081, imag(CharmLoops::A(mu, 1.0 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8385370486103412, imag(CharmLoops::A(mu, 1.1 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8487346187930029, imag(CharmLoops::A(mu, 1.2 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8581388308167810, imag(CharmLoops::A(mu, 1.3 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8668546781765816, imag(CharmLoops::A(mu, 1.4 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.8749677701646030, imag(CharmLoops::A(mu, 1.5 * m_b * m_b, m_b)), eps);
            }

            {
                static const double eps = 1e-7;
                static const double mu = 4.2;
                static const double m_b = 4.6;

                // CharmLoops::B
                TEST_CHECK_NEARLY_EQUAL( 1.7748361098250610, real(CharmLoops::B(mu, 0.1 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.6608414616625271, real(CharmLoops::B(mu, 0.2 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0549070842163954, real(CharmLoops::B(mu, 0.3 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.3492674916475516, real(CharmLoops::B(mu, 0.4 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.6458935344204433, real(CharmLoops::B(mu, 0.5 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.8762212415129942, real(CharmLoops::B(mu, 0.6 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.0618989188198120, real(CharmLoops::B(mu, 0.7 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.2156477427678207, real(CharmLoops::B(mu, 0.8 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.3455497477643420, real(CharmLoops::B(mu, 0.9 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.4570427884607030, real(CharmLoops::B(mu, 1.0 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.5539504651990608, real(CharmLoops::B(mu, 1.1 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.6390570203043875, real(CharmLoops::B(mu, 1.2 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.7144484220915650, real(CharmLoops::B(mu, 1.3 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.7817249058911009, real(CharmLoops::B(mu, 1.4 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.8421388765446847, real(CharmLoops::B(mu, 1.5 * m_b * m_b, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(5.553840794326391, imag(CharmLoops::B(mu, 0.1 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(4.820249145955325, imag(CharmLoops::B(mu, 0.2 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(4.360108135584585, imag(CharmLoops::B(mu, 0.3 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(4.017944023522331, imag(CharmLoops::B(mu, 0.4 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(3.743067928416354, imag(CharmLoops::B(mu, 0.5 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(3.512167767255114, imag(CharmLoops::B(mu, 0.6 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(3.312463519105912, imag(CharmLoops::B(mu, 0.7 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(3.136141749913582, imag(CharmLoops::B(mu, 0.8 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(2.978054117771269, imag(CharmLoops::B(mu, 0.9 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(2.834617374526911, imag(CharmLoops::B(mu, 1.0 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(2.703230935969235, imag(CharmLoops::B(mu, 1.1 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(2.581943911053763, imag(CharmLoops::B(mu, 1.2 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(2.469253015913550, imag(CharmLoops::B(mu, 1.3 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(2.363973892039418, imag(CharmLoops::B(mu, 1.4 * m_b * m_b, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(2.265155869018138, imag(CharmLoops::B(mu, 1.5 * m_b * m_b, m_b)), eps);
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
                TEST_CHECK_NEARLY_EQUAL(- 0.8832611, real(CharmLoops::F17_massless(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 0.6937322, imag(CharmLoops::F17_massless(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 5.2995666, real(CharmLoops::F27_massless(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 4.1623936, imag(CharmLoops::F27_massless(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 3.3632062, real(CharmLoops::F19_massless(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 6.9078480, imag(CharmLoops::F19_massless(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 3.4455298, real(CharmLoops::F29_massless(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+24.6919276, imag(CharmLoops::F29_massless(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(- 1.2486221, real(CharmLoops::F87_massless(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 2.7925269, imag(CharmLoops::F87_massless(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(- 3.2730189, real(CharmLoops::F89_massless(s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(  0.0000000, imag(CharmLoops::F89_massless(s, m_b)), eps);
            }

            /* Formfactors for O_8 are problematic near the zero recoil point */
            {
                static const double mu = 4.2, s = 19.2, m_b = 4.6, eps = 0.0000001;

                TEST_CHECK_NEARLY_EQUAL(- 0.9708796,  real(CharmLoops::F87_massless(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 2.7925268,  imag(CharmLoops::F87_massless(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(- 2.0208146,  real(CharmLoops::F89_massless(s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(  0.0000000,  imag(CharmLoops::F89_massless(s, m_b)), eps);
            }

            /* Check imaginary parts of the form factors at the boundaries of the Low Recoil region */
            {
                static const double mu = 4.2, s1 = 14.0, s2 = 19.2, m_b = 4.6, eps = 0.0000001;

                TEST_CHECK_NEARLY_EQUAL(+0.7802809, imag(CharmLoops::A(mu, s1, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.8161455, imag(CharmLoops::A(mu, s2, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-2.7925268, imag(CharmLoops::F87_massless(mu, s1, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-2.7925268, imag(CharmLoops::F87_massless(mu, s2, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-5.8682288, imag(CharmLoops::F19_massless(mu, s1, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-5.4492617, imag(CharmLoops::F19_massless(mu, s2, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(18.4542117, imag(CharmLoops::F29_massless(mu, s1, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(15.9404096, imag(CharmLoops::F29_massless(mu, s2, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0000000, imag(CharmLoops::F89_massless(s1, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0000000, imag(CharmLoops::F89_massless(s2, m_b)), eps);
            }

            /* Formfactors, massive loops */
            {
                static const double mu = 4.2, s = 6.0, m_b = 4.6, m_c = 1.2, eps = 0.0000001;
                TEST_CHECK_NEARLY_EQUAL(- 0.73093991, real(CharmLoops::F17_massive(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 0.17771334, imag(CharmLoops::F17_massive(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 4.38563254, real(CharmLoops::F27_massive(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 1.06627403, imag(CharmLoops::F27_massive(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(-34.40870331, real(CharmLoops::F19_massive(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 0.25864665, imag(CharmLoops::F19_massive(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 6.27364439, real(CharmLoops::F29_massive(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 1.55195807, imag(CharmLoops::F29_massive(mu, s, m_b, m_c)), eps);
            }
        }
} two_loop_test;
