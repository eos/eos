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
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_m1(0.5)), -0.261172396481, eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_m1(0.5)),  0.0,            eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_m1(1.0)), -0.548311355616, eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_m1(1.0)),  0.0,            eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_m1(2.0)), -1.23370055014,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_m1(2.0)),  0.0,            eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_m1(3.5)), -2.92543799360,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_m1(3.5)),  0.0,            eps);

            /* G_{-1}, t > 4 */
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_m1(20.0)), -0.76663543016,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_m1(20.0)), -9.07062920657,  eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_m1(16.0)), -1.46604599600,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_m1(16.0)), -8.27469050813,  eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_m1(8.0)),  -3.38116340075,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_m1(8.0)),  -5.53783357210,  eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_m1(4.5)),  -4.69457569359,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_m1(4.5)),  -2.17758609030,  eps);

            /* G_0, t < 4 */
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_0(0.5)), -0.0878249162964, eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_0(0.5)),  0.0,             eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_0(1.0)), -0.186200635766,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_0(1.0)),  0.0,             eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_0(2.0)), -0.429203673205,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_0(2.0)),  0.0,             eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_0(3.5)), -1.08575745738,   eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_0(3.5)),  0.0,             eps);

            /* G_0, t > 4 */
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_0(20.0)), +0.58245364578, eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_0(20.0)), -2.80992589242, eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_0(16.0)), +0.28103798890, eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_0(16.0)), -2.72069904635, eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_0(8.0)),  -0.75354951972, eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_0(8.0)),  -2.22144146908, eps);
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::G_0(4.5)),  -1.76895093981, eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::G_0(4.5)),  -1.04719755120, eps);
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
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::Deltai_23(0.04, 0.05, z)), +0.0220301686064, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::Deltai_23(0.04, 0.05, z)),  0.0,             eps);
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::Deltai_23(0.04, 0.50, z)), +0.65483498715,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::Deltai_23(0.04, 0.50, z)), +2.66785798352,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::Deltai_23(0.04, 0.75, z)), -0.67280659329,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::Deltai_23(0.04, 0.75, z)), +2.59469236324,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::Deltai_23(0.04, 1.00, z)), -1.25492085502,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::Deltai_23(0.04, 1.00, z)), +2.31191553533,   eps);

                /* Deltai_27 */
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::Deltai_27(0.04, 0.05, z)), 0.0438646709279, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::Deltai_27(0.04, 0.05, z)), 0.0,             eps);
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::Deltai_27(0.04, 0.50, z)), 2.38748258093,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::Deltai_27(0.04, 0.50, z)), 3.55121967487,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::Deltai_27(0.04, 0.75, z)), 1.03104194045,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::Deltai_27(0.04, 0.75, z)), 4.64401909136,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::Deltai_27(0.04, 1.00, z)), 0.15251167173,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::Deltai_27(0.04, 1.00, z)), 5.10340958495,   eps);
            }

            {
                static const double eps = 5.0e-11;

                /* tau_22 */
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::tau_22(0.04, 1.0/3.0, z)), 40.8001578505,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::tau_22(0.04, 1.0/3.0, z)),  0.0,            eps);

                /* tau_27 */
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::tau_27(0.04, 1.0/3.0, z)), 36.2911473431,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::tau_27(0.04, 1.0/3.0, z)),  0.0,            eps);

                /* tau_29 */
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::tau_28(0.04, 1.0/3.0, z)), -0.646314459835, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::tau_28(0.04, 1.0/3.0, z)),  0.0,            eps);

                /* tau_29 */
                TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::tau_29(0.04, 1.0/3.0, z)),  1.44526998657,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::tau_29(0.04, 1.0/3.0, z)),  0.0,            eps);
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
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::itau_22(0.04, z)), 18.6097,   eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::itau_22(0.04, z)),  0.0,      eps);

            /* tau_27 */
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::itau_27(0.04, z)),  9.4495,   eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::itau_27(0.04, z)), 21.2929,   eps);

            /* tau_28 */
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::itau_28(0.04, z)),  1.0778,   eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::itau_28(0.04, z)), -2.48846,  eps);

            /* tau_29 */
            TEST_CHECK_NEARLY_EQUAL(real(Bremsstrahlung::itau_29(0.04, z)),  0.594135, eps);
            TEST_CHECK_NEARLY_EQUAL(imag(Bremsstrahlung::itau_29(0.04, z)),  0.762988, eps);

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

            TEST_CHECK_NEARLY_EQUAL(Bremsstrahlung::tau_78(0.04),  79.656100,  eps);
            TEST_CHECK_NEARLY_EQUAL(Bremsstrahlung::tau_78(0.20),   5.311830,  eps);
            TEST_CHECK_NEARLY_EQUAL(Bremsstrahlung::tau_88(0.04),  46.770500,  eps);
            TEST_CHECK_NEARLY_EQUAL(Bremsstrahlung::tau_88(0.20),   0.622935,  eps);
            TEST_CHECK_NEARLY_EQUAL(Bremsstrahlung::tau_89(0.04),   3.517410,  eps);
            TEST_CHECK_NEARLY_EQUAL(Bremsstrahlung::tau_89(0.20),   0.762566,  eps);
        }
} type_a_test;
