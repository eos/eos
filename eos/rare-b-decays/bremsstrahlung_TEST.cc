/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2013 Danny van Dyk
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
#include <eos/rare-b-decays/bremsstrahlung.hh>

#include <cmath>
#include <iostream>
#include <limits>

using namespace test;
using namespace eos;

class IntegralsTest :
    public TestCase
{
    public:
        IntegralsTest() :
            TestCase("integrals_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Mathematica results from May 2010 */

            static const double eps = 1.0e-11;

            /* G_{-1}, t < 4 */
            TEST_CHECK_NEARLY_EQUAL(-0.261172396481, real(Bremsstrahlung::G_m1(0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,            imag(Bremsstrahlung::G_m1(0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL(-0.548311355616, real(Bremsstrahlung::G_m1(1.0)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,            imag(Bremsstrahlung::G_m1(1.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-1.23370055014,  real(Bremsstrahlung::G_m1(2.0)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,            imag(Bremsstrahlung::G_m1(2.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-2.92543799360,  real(Bremsstrahlung::G_m1(3.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,            imag(Bremsstrahlung::G_m1(3.5)), eps);

            /* G_{-1}, t > 4 */
            TEST_CHECK_NEARLY_EQUAL(-0.76663543016,  real(Bremsstrahlung::G_m1(20.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-9.07062920657,  imag(Bremsstrahlung::G_m1(20.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-1.46604599600,  real(Bremsstrahlung::G_m1(16.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-8.27469050813,  imag(Bremsstrahlung::G_m1(16.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-3.38116340075,  real(Bremsstrahlung::G_m1(8.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-5.53783357210,  imag(Bremsstrahlung::G_m1(8.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-4.69457569359,  real(Bremsstrahlung::G_m1(4.5)), eps);
            TEST_CHECK_NEARLY_EQUAL(-2.17758609030,  imag(Bremsstrahlung::G_m1(4.5)), eps);

            /* G_0, t < 4 */
            TEST_CHECK_NEARLY_EQUAL(-0.0878249162964, real(Bremsstrahlung::G_0(0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,             imag(Bremsstrahlung::G_0(0.5)), eps);
            TEST_CHECK_NEARLY_EQUAL(-0.186200635766,  real(Bremsstrahlung::G_0(1.0)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,             imag(Bremsstrahlung::G_0(1.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-0.429203673205,  real(Bremsstrahlung::G_0(2.0)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,             imag(Bremsstrahlung::G_0(2.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-1.08575745738,   real(Bremsstrahlung::G_0(3.5)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,             imag(Bremsstrahlung::G_0(3.5)), eps);

            /* G_0, t > 4 */
            TEST_CHECK_NEARLY_EQUAL(+0.58245364578, real(Bremsstrahlung::G_0(20.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-2.80992589242, imag(Bremsstrahlung::G_0(20.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(+0.28103798890, real(Bremsstrahlung::G_0(16.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-2.72069904635, imag(Bremsstrahlung::G_0(16.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-0.75354951972, real(Bremsstrahlung::G_0(8.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-2.22144146908, imag(Bremsstrahlung::G_0(8.0)), eps);
            TEST_CHECK_NEARLY_EQUAL(-1.76895093981, real(Bremsstrahlung::G_0(4.5)), eps);
            TEST_CHECK_NEARLY_EQUAL(-1.04719755120, imag(Bremsstrahlung::G_0(4.5)), eps);
        }
} integrals_test;

class BuildingBlocksTest :
    public TestCase
{
    public:
        BuildingBlocksTest() :
            TestCase("building_blocks_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Mathematica results from May 2010 */

            static const double z = pow(1.4 / 4.8, 2);

            {
                static const double eps = 1.0e-11;

                /* Deltai_23 */
                TEST_CHECK_NEARLY_EQUAL(+0.0220301686064, real(Bremsstrahlung::Deltai_23(0.04, 0.05, z)), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,             imag(Bremsstrahlung::Deltai_23(0.04, 0.05, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.65483498715,   real(Bremsstrahlung::Deltai_23(0.04, 0.50, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(+2.66785798352,   imag(Bremsstrahlung::Deltai_23(0.04, 0.50, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.67280659329,   real(Bremsstrahlung::Deltai_23(0.04, 0.75, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(+2.59469236324,   imag(Bremsstrahlung::Deltai_23(0.04, 0.75, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.25492085502,   real(Bremsstrahlung::Deltai_23(0.04, 1.00, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(+2.31191553533,   imag(Bremsstrahlung::Deltai_23(0.04, 1.00, z)), eps);

                /* Deltai_27 */
                TEST_CHECK_NEARLY_EQUAL(0.0438646709279, real(Bremsstrahlung::Deltai_27(0.04, 0.05, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.0,             imag(Bremsstrahlung::Deltai_27(0.04, 0.05, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(2.38748258093,   real(Bremsstrahlung::Deltai_27(0.04, 0.50, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(3.55121967487,   imag(Bremsstrahlung::Deltai_27(0.04, 0.50, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(1.03104194045,   real(Bremsstrahlung::Deltai_27(0.04, 0.75, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(4.64401909136,   imag(Bremsstrahlung::Deltai_27(0.04, 0.75, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(0.15251167173,   real(Bremsstrahlung::Deltai_27(0.04, 1.00, z)), eps);
                TEST_CHECK_NEARLY_EQUAL(5.10340958495,   imag(Bremsstrahlung::Deltai_27(0.04, 1.00, z)), eps);
            }

            {
                static const double eps = 5.0e-11;

                /* tau_22 */
                TEST_CHECK_NEARLY_EQUAL(40.8001578505,   real(Bremsstrahlung::tau_22(0.04, 1.0/3.0, z)), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,            imag(Bremsstrahlung::tau_22(0.04, 1.0/3.0, z)), eps);

                /* tau_27 */
                TEST_CHECK_NEARLY_EQUAL(36.2911473431,   real(Bremsstrahlung::tau_27(0.04, 1.0/3.0, z)), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,            imag(Bremsstrahlung::tau_27(0.04, 1.0/3.0, z)), eps);

                /* tau_29 */
                TEST_CHECK_NEARLY_EQUAL(-0.646314459835, real(Bremsstrahlung::tau_28(0.04, 1.0/3.0, z)), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,            imag(Bremsstrahlung::tau_28(0.04, 1.0/3.0, z)), eps);

                /* tau_29 */
                TEST_CHECK_NEARLY_EQUAL( 1.44526998657,  real(Bremsstrahlung::tau_29(0.04, 1.0/3.0, z)), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,            imag(Bremsstrahlung::tau_29(0.04, 1.0/3.0, z)), eps);
            }
        }
} building_blocks_test;

class IntegrationTest :
    public TestCase
{
    public:
        IntegrationTest() :
            TestCase("integration_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Mathematica results from May 2010 */

            static const double z = pow(1.4 / 4.8, 2);

            static const double eps = 1.0e-2; // all of the tau_2x are logarithmically divergent at w = s_hat.

            /* tau_22 */
            TEST_CHECK_NEARLY_EQUAL(18.6097,   real(Bremsstrahlung::itau_22(0.04, z)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,      imag(Bremsstrahlung::itau_22(0.04, z)), eps);

            /* tau_27 */
            TEST_CHECK_NEARLY_EQUAL( 9.4495,   real(Bremsstrahlung::itau_27(0.04, z)), eps);
            TEST_CHECK_NEARLY_EQUAL(21.2929,   imag(Bremsstrahlung::itau_27(0.04, z)), eps);

            /* tau_28 */
            TEST_CHECK_NEARLY_EQUAL( 1.0778,   real(Bremsstrahlung::itau_28(0.04, z)), eps);
            TEST_CHECK_NEARLY_EQUAL(-2.48846,  imag(Bremsstrahlung::itau_28(0.04, z)), eps);

            /* tau_29 */
            TEST_CHECK_NEARLY_EQUAL( 0.594135, real(Bremsstrahlung::itau_29(0.04, z)), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.762988, imag(Bremsstrahlung::itau_29(0.04, z)), eps);

        }
} integration_test;

class TypeATest :
    public TestCase
{
    public:
        TypeATest() :
            TestCase("type_A_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-4;

            TEST_CHECK_NEARLY_EQUAL( 79.656100,      Bremsstrahlung::tau_78(0.04), eps);
            TEST_CHECK_NEARLY_EQUAL(  5.311830,      Bremsstrahlung::tau_78(0.20), eps);
            TEST_CHECK_NEARLY_EQUAL( 46.770500,      Bremsstrahlung::tau_88(0.04), eps);
            TEST_CHECK_NEARLY_EQUAL(  0.622935,      Bremsstrahlung::tau_88(0.20), eps);
            TEST_CHECK_NEARLY_EQUAL(  3.517410,      Bremsstrahlung::tau_89(0.04), eps);
            TEST_CHECK_NEARLY_EQUAL(  0.762566,      Bremsstrahlung::tau_89(0.20), eps);
        }
} type_a_test;
