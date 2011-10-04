/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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
#include <eos/rare-b-decays/qcdf_integrals.hh>

#include <iostream>

using namespace test;
using namespace eos;

class QCDFIntegralsSToZeroTest :
    public TestCase
{
    public:
        QCDFIntegralsSToZeroTest() :
            TestCase("qcdf_integrals_s_to_zero_test")
        {
        }

        virtual void run() const
        {
            static const double m_B = 5.279, m_Kstar = 0.892, mu = 4.2;
            static const double m_b = 4.8, m_c = 1.6;
            static const double eps = 1e-10;

            // Asymptotic LCDA, shat = 0, bottom
            {
                QCDFIntegrals::Results results = QCDFIntegrals::photon_bottom_case(m_b, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+3.0,                    real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+3.0,                    real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-0.44013571266972111917, real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.44013571266972111917, real(results.j1_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(-0.50461362909944732565, real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.50461362909944732565, real(results.j4_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j4_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(-1.5740636314593566396,  real(results.j5_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(-0.53472500117995465696, real(results.j6_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+8.7095926471677816409,  results.j7_perp,           eps);
            }

            // Asymptotic LCDA, shat = 0, charm
            {
                QCDFIntegrals::Results results = QCDFIntegrals::photon_charm_case(m_c, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+3.0,                    real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+3.0,                    real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-3.9822376720419326578,  real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-3.9822376720419326578,  real(results.j1_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(-8.9150230855055024194,  imag(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-8.9150230855055024194,  imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(+0.8983587092605624785,  real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.8983587092605624785,  real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.7311284420639714829,  imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.7311284420639714829,  imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+2.4384255761102856504,  real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.1922235139632131500,  imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+0.7700334334248615859,  real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.2305475359496208336,  imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+8.7095926471677816409,  results.j7_perp,           eps);
            }

            // Asymptotic LCDA, shat = 0, massless
            {
                QCDFIntegrals::Results results = QCDFIntegrals::photon_massless_case(m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+3.0,                    real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+3.0,                    real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(+3.0,                    real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+3.0,                    real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(+0.4634203017314163927,  real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615,  imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.4634203017314163927,  real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615,  imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+2.2791497940831380669,  real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+4.1887902047863909846,  imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+0.9078647461758608371,  real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615,  imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+8.7095926471677816409,  results.j7_perp,           eps);
            }


            // Full LCDA, shat = 0, bottom
            {
                QCDFIntegrals::Results results = QCDFIntegrals::photon_bottom_case(m_b, m_B, m_Kstar, mu, +1.0, +2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+12.0,                   real(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-0.41907060371815448608, real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.41287205075314462248, real(results.j1_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(-0.53860597545218056924, real(results.j4_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.54504606267212781403, real(results.j4_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(-6.6468035561425338498,  real(results.j5_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(-0.5534251381472675701,  real(results.j6_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+35.571542061706523982,  results.j7_perp,           eps);
            }

            // Full LCDA, shat = 0, charm
            {
                QCDFIntegrals::Results results = QCDFIntegrals::photon_charm_case(m_c, m_B, m_Kstar, mu, +1.0, +2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+12.0,                   real(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-2.9828183887166422735,  real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+7.9494977585951323552,  imag(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-13.358084235247275442,  real(results.j1_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(-19.157182706492162991,  imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(+0.26453100712750606627, real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.18745275844653243462, imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.7237114167923608615,  real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.3523683673934455985,  imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+5.7162706573906113261,  real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.7583413149864139741,  imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+0.84364449498494953339, real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.20124338212004477811, imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+35.571542061706523982,  results.j7_perp,           eps);
            }

            // Full LCDA, shat = 0, massless
            {
                QCDFIntegrals::Results results = QCDFIntegrals::photon_massless_case(m_B, m_Kstar, mu, +1.0, +2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+12.0,                   real(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(+12.0,                   real(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(+1.0634203017314163927,  real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615,  imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.5300869683980830593,  real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615,  imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+16.449932509665885601,  real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+16.755160819145563938,  imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+0.9745314128425276445,  real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615,  imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+35.571542061706523982,  results.j7_perp,           eps);
            }
        }
} qcdf_integrals_s_to_zero_test;

class QCDFIntegralsDileptonBottomTest :
    public TestCase
{
    public:
        QCDFIntegralsDileptonBottomTest() :
            TestCase("qcdf_dilepton_bottom_test")
        {
        }

        virtual void run() const
        {
            static const double m_b = 4.8, m_B = 5.279, m_Kstar = 0.892, mu = 4.2;
            static const double eps = 1e-10;

            // Asymptotic LCDA, s = 1 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_bottom_case(1.0, m_b, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+2.5438661014683875729,   real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+2.5438661014683875729,   real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-0.129581012569702833082, real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.129581012569702833082, real(results.j1_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(+0.61348263455239097749,  real(results.j2_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(+0.22346879971623810801,  real(results.j3_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(-0.50244603005062536430,  real(results.j4_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.50244603005062536430,  real(results.j4_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(-1.3169030952519950247,   real(results.j5_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(-0.53165427751500952580,  real(results.j6_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+6.8087856446925005,      results.j7_perp,           eps);
            }

            // Full LCDA, s = 1 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_bottom_case(1.0, m_b, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+7.5060595661946593437,  real(results.j0_perp),      eps);
                TEST_CHECK_RELATIVE_ERROR(+1.5344070918431294599,  real(results.j0_parallel),  eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),      eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel),  eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-0.19022979798387140377, real(results.j1_perp),      eps);
                TEST_CHECK_RELATIVE_ERROR(-0.10123452954925091472, real(results.j1_parallel),  eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j1_perp),      eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j1_parallel),  eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(+2.3767573440564869425,  real(results.j2_perp),      eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j2_perp),      eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(+0.19337963084186422126, real(results.j3_parallel),  eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j3_parallel),  eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(-0.53548524560201187502,  real(results.j4_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.54153642537272670884,  real(results.j4_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(-4.1156844701148025195,   real(results.j5_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(-0.5497760398280187033,   real(results.j6_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                     imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+23.649061913234279,      results.j7_perp,           eps);
            }

            // Asymptotic LCDA, s = 6 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_bottom_case(6.0, m_b, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+1.8152337379049806314,  real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.8152337379049806314,  real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-0.23208465387143242229, real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.23208465387143242229, real(results.j1_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(+0.51894979960588456348, real(results.j2_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(+0.24996635642094325187, real(results.j3_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(-0.49140482018284548016, real(results.j4_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.49140482018284548016, real(results.j4_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(-0.90507537649861638238, real(results.j5_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(-0.51593046142938424122, real(results.j6_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+3.3218659138027814,     results.j7_perp,           eps);
            }

            // Full LCDA, s = 6 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_bottom_case(6.0, m_b, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+3.2577075222152920859,  real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.9483244402401200377,  real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-0.64044720707780513634, real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.08998732072510043348, real(results.j1_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(+2.0203980709249948114,  real(results.j2_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(+0.13687651166349362156, real(results.j3_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(-0.51943282456822609501, real(results.j4_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.52366886593109383238, real(results.j4_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(-1.7211609067159169416,  real(results.j5_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(-0.53107001510168672098, real(results.j6_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+6.5286206952110030,    results.j7_perp,           eps);
            }
        }
} qcdf_integrals_dilepton_bottom_test;

class QCDFIntegralsDileptonCharmTest :
    public TestCase
{
    public:
        QCDFIntegralsDileptonCharmTest() :
            TestCase("qcdf_dilepton_charm_test")
        {
        }

        virtual void run() const
        {
            static const double m_c = 1.6, m_B = 5.279, m_Kstar = 0.892, mu = 4.2;
            static const double eps = 1e-10;

            // Asymptotic LCDA, s = 1 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_charm_case(1.0, m_c, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+2.5438661014683875729,  real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+2.5438661014683875729,  real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-1.2084795500217594518,  real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-2.5731752833018131227,  imag(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-1.2084795500217594518,  real(results.j1_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(-2.5731752833018131227,  imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(+7.5499378384070829504,  real(results.j2_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+3.8575235933762310573,  imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(+2.4927197132985721754,  real(results.j3_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+2.2579515246103501287,  imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(+0.90662000663822865676, real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.76810245512688750554, imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.90662000663822865676, real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.76810245512688750554, imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+2.2357284385466520804,  real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.2480087814547703380,  imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+0.80930949867046960465, real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.25725993069133215304, imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+6.8087856446925005,    results.j7_perp,           eps);
            }

            // Full LCDA, s = 1 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_charm_case(1.0, m_c, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+7.5060595661946593437,  real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.5344070918431294599,  real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-2.0259892377536594235,  real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+2.0901758815693526036,  imag(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-3.3603065663597425813,  real(results.j1_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(-5.9062192258577044449,  imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(+25.633224198327682961,  real(results.j2_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-3.2821702413938853550,  imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(+5.3764386950330671286,  real(results.j3_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+3.0051332554454714949,  imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(+0.3707067450311682311,  real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.2143763905157959554,  imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.7177556396681478383,  real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.4878614401751314190,  imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+4.0690777118525885441,  real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.7738438357039390889,  imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+0.9217204757571799689,  real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.1614905141057737525,  imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(23.649061913234279,      results.j7_perp,           eps);
            }

            // Asymptotic LCDA, s = 6 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_charm_case(6.0, m_c, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+1.8152337379049806314,  real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.8152337379049806314,  real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-2.5668375066590532237,  real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-6.2820634450983950390,  imag(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-2.5668375066590532237,  real(results.j1_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(-6.2820634450983950390,  imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(+8.3202230206970707438,  real(results.j2_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+7.3987151431737863524,  imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(+3.1801081101889446211,  real(results.j3_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+4.1971096225878986166,  imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(+0.90380432437933200684, real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.98191507574626079235, imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.90380432437933200684, real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.98191507574626079235, imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+1.7655917043498167032,  real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.5510177821620955223,  imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+1.02630625774849761007, real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.47714758980293173711, imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+3.3218659138027814,     results.j7_perp,           eps);
            }

            // Full LCDA, s = 6 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_charm_case(6.0, m_c, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+3.2577075222152920859,  real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.9483244402401200377,  real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                    imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(-16.318405730056516574,  real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-5.4853614861278124770,  imag(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.3732425912900265450,  real(results.j1_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(-13.668520917564193531,  imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(+47.588404767679882834,  real(results.j2_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+2.3701030060135062110,  imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(+1.7642242844228323230,  real(results.j3_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+8.3285909048480953900,  imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(+1.0761522299504679615,  real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.0275263680199818626,  imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3896136112058304143,  real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.1061607651804707851,  imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+3.5779030782052170195,  real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(-0.0983563145623310670,  imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+1.3033931133205887143,  real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.2500736812081865157,  imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+6.5286206952110030,    results.j7_perp,           eps);
            }
        }
} qcdf_integrals_dilepton_charm_test;

class QCDFIntegralsDileptonMasslessTest :
    public TestCase
{
    public:
        QCDFIntegralsDileptonMasslessTest() :
            TestCase("qcdf_dilepton_massless_test")
        {
        }

        virtual void run() const
        {
            static const double m_B = 5.279, m_Kstar = 0.892, mu = 4.2;
            static const double eps = 1e-13;

            // Asymptotic LCDA, s = 1 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_massless_case(1.0, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+2.5438661014683875729, real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+2.5438661014683875729, real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(+3.0,                   real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+3.0,                   real(results.j1_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(-27.431079556392130558, real(results.j2_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(-6.8921108596679530233, real(results.j3_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(+0.4354036325983407521, real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.4354036325983407521, real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+1.5127703537653745058, real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+3.5519071360396417498, imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+0.7784849884802046572, real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+6.8087856446925005,    results.j7_perp,           eps);
            }

            // Full LCDA, s = 1 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_massless_case(1.0, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+7.5060595661946593437, real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.5344070918431294599, real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(+12.0,                  real(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   real(results.j1_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(-187.90193454111942863, real(results.j2_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(-4.0012672343594354200, real(results.j3_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(+0.9317103385752230940, real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.5291482770336041636, real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+7.7196342667907491874, real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+10.480436262473125396, imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+0.8790914572324735055, real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+23.649061913234279,    results.j7_perp,           eps);
            }

            // Asymptotic LCDA, s = 6 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_massless_case(6.0, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+1.8152337379049806314, real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.8152337379049806314, real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(+3.0,                   real(results.j1_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+3.0,                   real(results.j1_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(-7.5994040642306060309, real(results.j2_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(-3.3303478297578050525, real(results.j3_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(+0.3349587633592996386, real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.3349587633592996386, real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+0.6938301178137627689, real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+2.5345444335780565994, imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+0.5161218162270731122, real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+3.3218659138027814,    results.j7_perp,           eps);
            }

            // Full LCDA, s = 6 GeV
            {
                QCDFIntegrals::Results results = QCDFIntegrals::dilepton_massless_case(6.0, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(+3.2577075222152920859, real(results.j0_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.9483244402401200377, real(results.j0_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j0_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j0_parallel), eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(+12.0,                  real(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   real(results.j1_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j1_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j1_parallel), eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(-38.580065526099502420, real(results.j2_perp),     eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j2_perp),     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(-1.0701497634748736200, real(results.j3_parallel), eps);
                TEST_CHECK_NEARLY_EQUAL(   0.0,                   imag(results.j3_parallel), eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(+0.5937069566838658442, real(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j4_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+0.4446370078996704249, real(results.j4_parallel), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j4_parallel), eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(+2.0577789868031868146, real(results.j5_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+4.5486177863714532314, imag(results.j5_perp),     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(+0.5992101322006142539, real(results.j6_perp),     eps);
                TEST_CHECK_RELATIVE_ERROR(+1.3962634015954636615, imag(results.j6_perp),     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(+6.5286206952110030,    results.j7_perp,           eps);
            }
        }
} qcdf_integrals_dilepton_massless_test;
