/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2016 Danny van Dyk
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
#include <eos/rare-b-decays/qcdf-integrals.hh>

#include <iostream>
#include <string>

using namespace test;
using namespace eos;

template <typename Tag_>
class QCDFIntegralsPhotonTest :
    public TestCase
{
    public:
        using Calculator = QCDFIntegralCalculator<BToKstarDilepton, Tag_>;

        QCDFIntegralsPhotonTest() :
            TestCase("qcdf_integrals_photon_test" + Tag_::name + ">")
        {
        }

        virtual void run() const
        {
            static const double m_B = 5.279, m_Kstar = 0.892, mu = 4.2;
            static const double m_b = 4.8, m_c = 1.6;
            static const double eps = 1e-10;

            // Asymptotic LCDA, shat = 0, bottom
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::photon_bottom_case(m_b, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +3.0,                    eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +3.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     -0.44013571266972111917, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), -0.44013571266972111917, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_parallel),  0.0,                    eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     -0.50461362909944732565, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), -0.50461362909944732565, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_parallel),  0.0,                    eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     -1.5740636314593566396,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j5_perp),      0.0,                    eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     -0.53472500117995465696, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j6_perp),      0.0,                    eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +8.7095926471677816409,  eps);
            }

            // Asymptotic LCDA, shat = 0, charm
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::photon_charm_case(m_c, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +3.0,                    eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +3.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     -3.9822376720419326578,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), -3.9822376720419326578,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j1_perp),     -8.9150230855055024194,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j1_parallel), -8.9150230855055024194,  eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +0.8983587092605624785,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +0.8983587092605624785,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     +0.7311284420639714829,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +0.7311284420639714829,  eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +2.4384255761102856504,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     +1.1922235139632131500,  eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +0.7700334334248615859,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     +0.2305475359496208336,  eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +8.7095926471677816409,  eps);
            }

            // Asymptotic LCDA, shat = 0, massless
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::photon_massless_case(m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +3.0,                    eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +3.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +3.0,                    eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +3.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +0.4634203017314163927,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     +1.3962634015954636615,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +0.4634203017314163927,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +1.3962634015954636615,  eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +2.2791497940831380669,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     +4.1887902047863909846,  eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +0.9078647461758608371,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     +1.3962634015954636615,  eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +8.7095926471677816409,  eps);
            }


            // Full LCDA, shat = 0, bottom
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::photon_bottom_case(m_b, m_B, m_Kstar, mu, +1.0, +2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +12.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  real(results.j0_parallel),  0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     -0.41907060371815448608, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), -0.41287205075314462248, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_parallel),  0.0,                    eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     -0.53860597545218056924, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_perp),      0.0,                    eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), -0.54504606267212781403, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_parallel),  0.0,                    eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     -6.6468035561425338498,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j5_perp),      0.0,                    eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     -0.5534251381472675701,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j6_perp),      0.0,                    eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +35.571542061706523982,  eps);
            }

            // Full LCDA, shat = 0, charm
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::photon_charm_case(m_c, m_B, m_Kstar, mu, +1.0, +2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +12.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  real(results.j0_parallel),  0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     -2.9828183887166422735,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j1_perp),     +7.9494977585951323552,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), -13.358084235247275442,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j1_parallel), -19.157182706492162991,  eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +0.26453100712750606627, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     -0.18745275844653243462, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +1.7237114167923608615,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +0.3523683673934455985,  eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +5.7162706573906113261,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     -0.7583413149864139741,  eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +0.84364449498494953339, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     -0.20124338212004477811, eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +35.571542061706523982,  eps);
            }

            // Full LCDA, shat = 0, massless
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::photon_massless_case(m_B, m_Kstar, mu, +1.0, +2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +12.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  real(results.j0_parallel),  0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +12.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  real(results.j0_parallel),  0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J2
                TEST_CHECK(std::isnan(real(results.j2_perp)));
                TEST_CHECK(std::isnan(imag(results.j2_perp)));

                // J3
                TEST_CHECK(std::isnan(real(results.j3_parallel)));
                TEST_CHECK(std::isnan(imag(results.j3_parallel)));

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +1.0634203017314163927,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     +1.3962634015954636615,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +0.5300869683980830593,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +1.3962634015954636615,  eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +16.449932509665885601,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     +16.755160819145563938,  eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +0.9745314128425276445,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     +1.3962634015954636615,  eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +35.571542061706523982,  eps);
            }
        }
};
QCDFIntegralsPhotonTest<tag::Analytical> qcdf_integrals_photon_test_analytical;
QCDFIntegralsPhotonTest<tag::Mixed> qcdf_integrals_photon_test_mixed;

template <typename Tag_>
class QCDFIntegralsDileptonBottomTest :
    public TestCase
{
    public:
        using Calculator = QCDFIntegralCalculator<BToKstarDilepton, Tag_>;

        QCDFIntegralsDileptonBottomTest() :
            TestCase("qcdf_dilepton_bottom_test<" + Tag_::name + ">")
        {
        }

        virtual void run() const
        {
            static const double m_b = 4.8, m_B = 5.279, m_Kstar = 0.892, mu = 4.2;
            static const double eps = 1e-3;

            // Asymptotic LCDA, s = 1 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_bottom_case(1.0, m_b, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +2.5438661014683875729,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +2.5438661014683875729,   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                     eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                     eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     -0.129581012569702833082, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), -0.129581012569702833082, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_perp),      0.0,                     eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_parallel),  0.0,                     eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),     +0.61348263455239097749,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j2_perp),      0.0,                     eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel), +0.22346879971623810801,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j3_parallel),  0.0,                     eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     -0.50244603005062536430,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_perp),      0.0,                     eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), -0.50244603005062536430,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_parallel),  0.0,                     eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     -1.3169030952519950247,   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j5_perp),      0.0,                     eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     -0.53165427751500952580,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j6_perp),      0.0,                     eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +6.8087856446925005,      eps);
            }

            // Full LCDA, s = 1 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_bottom_case(1.0, m_b, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),      +7.5060595661946593437,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel),  +1.5344070918431294599,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),       0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),   0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),      -0.19022979798387140377, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel),  -0.10123452954925091472, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_perp),       0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_parallel),   0.0,                    eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),      +2.3767573440564869425,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j2_perp),       0.0,                    eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel),  +0.19337963084186422126, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j3_parallel),   0.0,                    eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),      -0.53548524560201187502, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_perp),       0.0,                    eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel),  -0.54153642537272670884, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_parallel),   0.0,                    eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),      -4.1156844701148025195,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j5_perp),       0.0,                    eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),      -0.5497760398280187033,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j6_perp),       0.0,                    eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,            +23.649061913234279,     eps);
            }

            // Asymptotic LCDA, s = 6 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_bottom_case(6.0, m_b, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),      +1.8152337379049806314,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel),  +1.8152337379049806314,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),       0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),   0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),      -0.23208465387143242229, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel),  -0.23208465387143242229, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_perp),       0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_parallel),   0.0,                    eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),      +0.51894979960588456348, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j2_perp),       0.0,                    eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel),  +0.24996635642094325187, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j3_parallel),   0.0,                    eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),      -0.49140482018284548016, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_perp),       0.0,                    eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel),  -0.49140482018284548016, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_parallel),   0.0,                    eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),      -0.90507537649861638238, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j5_perp),       0.0,                    eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),      -0.51593046142938424122, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j6_perp),       0.0,                    eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,            +3.3218659138027814,     eps);
            }

            // Full LCDA, s = 6 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_bottom_case(6.0, m_b, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),      +3.2577075222152920859,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel),  +1.9483244402401200377,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),       0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),   0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),      -0.64044720707780513634, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel),  -0.08998732072510043348, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_perp),       0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_parallel),   0.0,                    eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),      +2.0203980709249948114,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j2_perp),       0.0,                    eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel),  +0.13687651166349362156, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j3_parallel),   0.0,                    eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),      -0.51943282456822609501, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_perp),       0.0,                    eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel),  -0.52366886593109383238, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j4_parallel),   0.0,                    eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),      -1.7211609067159169416,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j5_perp),       0.0,                    eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),      -0.53107001510168672098, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j6_perp),       0.0,                    eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,            +6.5286206952110030,     eps);
            }
        }
};
QCDFIntegralsDileptonBottomTest<tag::Analytical> qcdf_integrals_dilepton_bottom_test_analytical;
QCDFIntegralsDileptonBottomTest<tag::Mixed> qcdf_integrals_dilepton_bottom_test_mixed;
QCDFIntegralsDileptonBottomTest<tag::Numerical> qcdf_integrals_dilepton_bottom_test_numerical;

template <typename Tag_>
class QCDFIntegralsDileptonCharmTest :
    public TestCase
{
    public:
        using Calculator = QCDFIntegralCalculator<BToKstarDilepton, Tag_>;

        QCDFIntegralsDileptonCharmTest() :
            TestCase("qcdf_dilepton_charm_test<" + Tag_::name + ">")
        {
        }

        virtual void run() const
        {
            static const double m_c = 1.6, m_B = 5.279, m_Kstar = 0.892, mu = 4.2;
            static const double eps = 5.5e-3;

            // Asymptotic LCDA, s = 1 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_charm_case(1.0, m_c, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +2.5438661014683875729,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +2.5438661014683875729,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     -1.2084795500217594518,  eps);
                if (("numerical" == Tag_::name) || ("mixed" == Tag_::name))
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_perp),     -1.6953885768872480000,  eps); // correct
                }
                else
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_perp),     -2.5731752833018131227,  eps); // wrong!
                }
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), -1.2084795500217594518,  eps);
                if (("numerical" == Tag_::name) || ("mixed" == Tag_::name))
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_parallel), -1.6953885768872480000,  eps); // correct
                }
                else
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_parallel), -2.5731752833018131227,  eps); // wrong!
                }

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),     +7.5499378384070829504,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j2_perp),     +3.8575235933762310573,  eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel), +2.4927197132985721754,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j3_parallel), +2.2579515246103501287,  eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +0.90662000663822865676, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     +0.76810245512688750554, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +0.90662000663822865676, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +0.76810245512688750554, eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +2.2357284385466520804,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     +1.2480087814547703380,  eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +0.80930949867046960465, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     +0.25725993069133215304, eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +6.8087856446925005,     eps);
            }

            // Full LCDA, s = 1 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_charm_case(1.0, m_c, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +7.5060595661946593437,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +1.5344070918431294599,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     -2.0259892377536594235,  eps);
                if (("numerical" == Tag_::name) || ("mixed" == Tag_::name))
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_perp),     +1.3250338631629440000,  eps); // correct
                }
                else
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_perp),     +2.0901758815693526036,  eps); // wrong!
                }
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), -3.3603065663597425813,  eps);
                if (("numerical" == Tag_::name) || ("mixed" == Tag_::name))
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_parallel), -3.0735616834264360000,  eps); // correct
                }
                else
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_parallel), -5.9062192258577044449,  eps); // wrong!
                }

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),     +25.633224198327682961,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j2_perp),     -3.2821702413938853550,  eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel), +5.3764386950330671286,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j3_parallel), +3.0051332554454714949,  eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +0.3707067450311682311,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     -0.2143763905157959554,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +1.7177556396681478383,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +0.4878614401751314190,  eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +4.0690777118525885441,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     -0.7738438357039390889,  eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +0.9217204757571799689,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     -0.1614905141057737525,  eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           23.649061913234279,      eps);
            }

            // Asymptotic LCDA, s = 6 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_charm_case(6.0, m_c, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +1.8152337379049806314,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +1.8152337379049806314,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     -2.5668375066590532237,  eps);
                if (("numerical" == Tag_::name) || ("mixed" == Tag_::name))
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_perp),     -3.9575320372001630000,  eps); // correct
                }
                else
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_perp),     -6.2820634450983950390,  eps); // wrong!
                }
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), -2.5668375066590532237,  eps);
                if (("numerical" == Tag_::name) || ("mixed" == Tag_::name))
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_parallel), -3.9575320372001630000,  eps); // correct
                }
                else
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_parallel), -6.2820634450983950390,  eps); // wrong!
                }

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),     +8.3202230206970707438,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j2_perp),     +7.3987151431737863524,  eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel), +3.1801081101889446211,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j3_parallel), +4.1971096225878986166,  eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +0.90380432437933200684, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     +0.98191507574626079235, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +0.90380432437933200684, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +0.98191507574626079235, eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +1.7655917043498167032,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     +1.5510177821620955223,  eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +1.02630625774849761007, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     +0.47714758980293173711, eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +3.3218659138027814,     eps);
            }

            // Full LCDA, s = 6 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_charm_case(6.0, m_c, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +3.2577075222152920859,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +1.9483244402401200377,  eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                    eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                    eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     -16.318405730056516574,  eps);
                if (("numerical" == Tag_::name) || ("mixed" == Tag_::name))
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_perp),     -0.9626512455812828000,  eps); // correct
                }
                else
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_perp),     -5.4853614861278124770,  eps); // wrong!
                }
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), +0.3732425912900265450,  eps);
                if (("numerical" == Tag_::name) || ("mixed" == Tag_::name))
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_parallel), -8.7757310767704430000,  eps); // correct
                }
                else
                {
                    TEST_CHECK_RELATIVE_ERROR(imag(results.j1_parallel), -13.668520917564193531,  eps); // wrong!
                }

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),     +47.588404767679882834,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j2_perp),     +2.3701030060135062110,  eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel), +1.7642242844228323230,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j3_parallel), +8.3285909048480953900,  eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +1.0761522299504679615,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     -0.0275263680199818626,  10 * eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +1.3896136112058304143,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +1.1061607651804707851,  eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +3.5779030782052170195,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     -0.0983563145623310670,  eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +1.3033931133205887143,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     +0.2500736812081865157,  eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +6.5286206952110030,     eps);
            }
        }
};
QCDFIntegralsDileptonCharmTest<tag::Analytical> qcdf_integrals_dilepton_charm_test_analytical;
QCDFIntegralsDileptonCharmTest<tag::Mixed> qcdf_integrals_dilepton_charm_test_mixed;
QCDFIntegralsDileptonCharmTest<tag::Numerical> qcdf_integrals_dilepton_charm_test_numerical;

template <typename Tag_>
class QCDFIntegralsDileptonMasslessTest :
    public TestCase
{
    public:
        using Calculator = QCDFIntegralCalculator<BToKstarDilepton, Tag_>;

        QCDFIntegralsDileptonMasslessTest() :
            TestCase("qcdf_dilepton_massless_test<" + Tag_::name + ">")
        {
        }

        virtual void run() const
        {
            static const double m_B = 5.279, m_Kstar = 0.892, mu = 4.2;
            static const double eps = 5e-3;

            // Asymptotic LCDA, s = 1 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_massless_case(1.0, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +2.5438661014683875729, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +2.5438661014683875729, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                   eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     +3.0,                   eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), +3.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_perp),      0.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_parallel),  0.0,                   eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),     -27.431079556392130558, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j2_perp),      0.0,                   eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel), -6.8921108596679530233, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j3_parallel),  0.0,                   eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +0.4354036325983407521, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     +1.3962634015954636615, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +0.4354036325983407521, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +1.3962634015954636615, eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +1.5127703537653745058, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     +3.5519071360396417498, eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +0.7784849884802046572, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     +1.3962634015954636615, eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +6.8087856446925005,    eps);
            }

            // Full LCDA, s = 1 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_massless_case(1.0, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +7.5060595661946593437, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +1.5344070918431294599, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                   eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     +12.0,                  eps);
                TEST_CHECK_NEARLY_EQUAL(  real(results.j1_parallel),  0.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_perp),      0.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_parallel),  0.0,                   eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),     -187.90193454111942863, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j2_perp),      0.0,                   eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel), -4.0012672343594354200, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j3_parallel),  0.0,                   eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +0.9317103385752230940, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     +1.3962634015954636615, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +0.5291482770336041636, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +1.3962634015954636615, eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +7.7196342667907491874, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     +10.480436262473125396, eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +0.8790914572324735055, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     +1.3962634015954636615, eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +23.649061913234279,    eps);
            }

            // Asymptotic LCDA, s = 6 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_massless_case(6.0, m_B, m_Kstar, mu, 0.0, 0.0, 0.0, 0.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +1.8152337379049806314, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +1.8152337379049806314, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                   eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     +3.0,                   eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_parallel), +3.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_perp),      0.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_parallel),  0.0,                   eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),     -7.5994040642306060309, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j2_perp),      0.0,                   eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel), -3.3303478297578050525, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j3_parallel),  0.0,                   eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +0.3349587633592996386, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     +1.3962634015954636615, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +0.3349587633592996386, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +1.3962634015954636615, eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +0.6938301178137627689, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     +2.5345444335780565994, eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +0.5161218162270731122, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     +1.3962634015954636615, eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +3.3218659138027814,    eps);
            }

            // Full LCDA, s = 6 GeV
            {
                QCDFIntegrals<BToKstarDilepton> results = Calculator::dilepton_massless_case(6.0, m_B, m_Kstar, mu, 1.0, 2.0, 1.0, -2.0);

                // J0
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_perp),     +3.2577075222152920859, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j0_parallel), +1.9483244402401200377, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_perp),      0.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j0_parallel),  0.0,                   eps);

                // J1
                TEST_CHECK_RELATIVE_ERROR(real(results.j1_perp),     +12.0,                  eps);
                TEST_CHECK_NEARLY_EQUAL(  real(results.j1_parallel),  0.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_perp),      0.0,                   eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j1_parallel),  0.0,                   eps);

                // J2
                TEST_CHECK_RELATIVE_ERROR(real(results.j2_perp),     -38.580065526099502420, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j2_perp),      0.0,                   eps);

                // J3
                TEST_CHECK_RELATIVE_ERROR(real(results.j3_parallel), -1.0701497634748736200, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(results.j3_parallel),  0.0,                   eps);

                // J4
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_perp),     +0.5937069566838658442, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_perp),     +1.3962634015954636615, eps);
                TEST_CHECK_RELATIVE_ERROR(real(results.j4_parallel), +0.4446370078996704249, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j4_parallel), +1.3962634015954636615, eps);

                // J5
                TEST_CHECK_RELATIVE_ERROR(real(results.j5_perp),     +2.0577789868031868146, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j5_perp),     +4.5486177863714532314, eps);

                // J6
                TEST_CHECK_RELATIVE_ERROR(real(results.j6_perp),     +0.5992101322006142539, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(results.j6_perp),     +1.3962634015954636615, eps);

                // J7
                TEST_CHECK_RELATIVE_ERROR(results.j7_perp,           +6.5286206952110030,    eps);
            }
        }
};
QCDFIntegralsDileptonMasslessTest<tag::Analytical> qcdf_integrals_dilepton_massless_test_analytical;
QCDFIntegralsDileptonMasslessTest<tag::Mixed> qcdf_integrals_dilepton_massless_test_mixed;
QCDFIntegralsDileptonMasslessTest<tag::Numerical>  qcdf_integrals_dilepton_massless_test_numerical;
