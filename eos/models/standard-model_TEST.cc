/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015 Danny van Dyk
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
#include <eos/models/model.hh>
#include <eos/models/standard-model.hh>

#include <cmath>

using namespace test;
using namespace eos;

Parameters
reference_parameters()
{
    Parameters result = Parameters::Defaults();
    result["QCD::alpha_s(MZ)"] = 0.117620;
    result["QCD::mu_t"] = 170.0;
    result["QCD::mu_b"] = 4.2;
    result["QCD::mu_c"] = 1.2;
    result["mass::t(pole)"] = 173.3;
    result["mass::b(MSbar)"] = 4.2;
    result["mass::c"] = 1.27;
    result["mass::s(2GeV)"] = 0.101;
    result["mass::d(2GeV)"] = 0.0048;
    result["mass::u(2GeV)"] = 0.0032;
    result["mass::W"] = 80.398;
    result["mass::Z"] = 91.1876;
    result["GSW::sin^2(theta)"] = 0.23116;
    result["CKM::A"] = 0.812;
    result["CKM::lambda"] = 0.22543;
    result["CKM::rhobar"] = 0.144;
    result["CKM::etabar"] = 0.342;
    // WET sectors
    result["sbsb::mu"]   =  4.2;
    result["sbcu::mu_0"] = 80.0;
    result["sbcu::mu"  ] =  4.2;
    result["dbcu::mu_0"] = 80.0;
    result["dbcu::mu"  ] =  4.2;

    return result;
}

class MakeTest :
    public TestCase
{
    public:
        MakeTest() :
            TestCase("sm_make_test")
        {
        }

        virtual void run() const
        {
            try
            {
                std::shared_ptr<Model> m = Model::make("SM", reference_parameters(), Options());
            }
            catch (NoSuchModelError &)
            {
                TEST_CHECK_FAILED("Model::make does not know the model 'SM'");
            }
            catch (...)
            {
                TEST_CHECK_FAILED("Unknown Exception while making 'SM'");
            }
        }
} sm_make_test;

class AlphaSTest :
    public TestCase
{
    public:
        AlphaSTest() :
            TestCase("sm_alpha_s_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-6;

            StandardModel model(reference_parameters());

            // Calculation of alpha_s is not self-consistent:
            //   alpha_s(mu) != alpha_s_0
            // So check for relative error
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s(91.1876), 0.117620, 5e-5);

            // Data in agreement with RunDec, cf. [CKS2000]
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s(120.0), 0.112968, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s(80.403),0.119918, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s(80.0 ), 0.120011, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s(40.0 ), 0.134400, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s(20.0 ), 0.152867, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s(10.0 ), 0.177507, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 9.6 ), 0.179220, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 4.8 ), 0.214716, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 4.45), 0.219518, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 4.2 ), 0.223342, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 3.0 ), 0.252878, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 2.4 ), 0.277227, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 2.0 ), 0.301404, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 1.2 ), 0.405724, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 1.0 ), 0.490620, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 0.7 ), 0.883896, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 0.6 ), 1.524938, eps);
            TEST_CHECK_NEARLY_EQUAL(model.alpha_s( 0.5 ), 5.709652, eps);
        }
} sm_alpha_s_test;

class TMassesTest :
    public TestCase
{
    public:
        TMassesTest() :
            TestCase("sm_t_masses_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            StandardModel model(reference_parameters());

            TEST_CHECK_RELATIVE_ERROR(model.m_t_msbar(120.0), 167.794, eps);
            TEST_CHECK_RELATIVE_ERROR(model.m_t_msbar( 80.0), 173.647, eps);
        }
} sm_t_masses_test;

class BMassesTest :
    public TestCase
{
    public:
        BMassesTest() :
            TestCase("sm_b_masses_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = reference_parameters();
            StandardModel model(p);

            TEST_CHECK_NEARLY_EQUAL(model.m_b_msbar(9.6), 3.67956, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_b_msbar(4.8), 4.10051, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_b_msbar(4.2), 4.20000, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_b_msbar(2.4), 4.75221, eps);

            TEST_CHECK_NEARLY_EQUAL(model.m_b_pole(),     4.74167, 1e-3); // Precision is hard-limited in fixed-point routine

            TEST_CHECK_NEARLY_EQUAL(model.m_b_ps(1.0),    4.60728, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_b_ps(1.5),    4.54012, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_b_ps(2.0),    4.47735, eps);

            TEST_CHECK_NEARLY_EQUAL(model.m_b_kin(0.75),  4.63362, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_b_kin(1.00),  4.56114, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_b_kin(1.25),  4.49203, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_b_kin(1.50),  4.42520, eps);
        }
} sm_b_masses_test;

class CMassesTest :
    public TestCase
{
    public:
        CMassesTest() :
            TestCase("sm_c_masses_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-6;

            Parameters p = reference_parameters();
            StandardModel model(p);

            TEST_CHECK_NEARLY_EQUAL(model.m_c_msbar(4.8),  0.891000, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_c_msbar(4.2),  0.912618, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_c_msbar(1.27), 1.270000, eps);

            TEST_CHECK_NEARLY_EQUAL(model.m_c_pole(),      1.595301, eps);

            TEST_CHECK_NEARLY_EQUAL(model.m_c_kin(0.75),   1.203723, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_c_kin(1.00),   1.060682, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_c_kin(1.25),   0.931772, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_c_kin(1.50),   0.813366, eps);
        }
} sm_c_masses_test;

class SMassesTest :
    public TestCase
{
    public:
        SMassesTest() :
            TestCase("sm_s_masses_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-6;

            Parameters p = reference_parameters();
            StandardModel model(p);

            TEST_CHECK_NEARLY_EQUAL(model.m_s_msbar(1.0), 0.136682, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_s_msbar(1.7), 0.106128, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_s_msbar(2.0), 0.101000, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_s_msbar(4.2), 0.084980, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_s_msbar(4.8), 0.082967, eps);
        }
} sm_s_masses_test;

class UDMassesTest :
    public TestCase
{
    public:
        UDMassesTest() :
            TestCase("sm_ud_masses_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-6;

            Parameters p = reference_parameters();
            StandardModel model(p);

            TEST_CHECK_NEARLY_EQUAL(model.m_ud_msbar(1.0), 0.010826, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_ud_msbar(2.0), 0.008000, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_ud_msbar(3.0), 0.007223, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_ud_msbar(4.0), 0.006803, eps);
            TEST_CHECK_NEARLY_EQUAL(model.m_ud_msbar(5.0), 0.006525, eps);
        }
} sm_ud_masses_test;

class CKMElementsTest :
    public TestCase
{
    public:
        CKMElementsTest() :
            TestCase("ckm_elements_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-8;

            // central values
            {
                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                // values
                TEST_CHECK_NEARLY_EQUAL(real(model.ckm_ud()), +0.974253267, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(model.ckm_ud()), +0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_ud()),  +0.97425,     1e-5);

                TEST_CHECK_NEARLY_EQUAL(real(model.ckm_us()), +0.225428590, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(model.ckm_us()), +0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_us()),  +0.22543,     1e-5);

                TEST_CHECK_NEARLY_EQUAL(real(model.ckm_ub()), +0.001372189, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(model.ckm_ub()), -0.003264270, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_ub()),  +0.00354,     1e-5);

                TEST_CHECK_NEARLY_EQUAL(real(model.ckm_cd()), +0.225296132, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(model.ckm_cd()), +0.000138121, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_cd()),  +0.22529,     1e-5);

                TEST_CHECK_NEARLY_EQUAL(real(model.ckm_cs()), +0.973416767, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(model.ckm_cs()), -0.000030365, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_cs()),  +0.97342,     1e-5);

                TEST_CHECK_NEARLY_EQUAL(real(model.ckm_cb()), +0.041264513, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(model.ckm_cb()), +0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_cb()),  +0.04126,     1e-5);

                TEST_CHECK_NEARLY_EQUAL(real(model.ckm_td()), +0.007966605, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(model.ckm_td()), -0.003177489, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_td()),  +0.00858,     1e-5);

                TEST_CHECK_NEARLY_EQUAL(real(model.ckm_ts()), -0.040511671, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(model.ckm_ts()), -0.000735237, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_ts()),  +0.04052,     1e-5);

                TEST_CHECK_NEARLY_EQUAL(real(model.ckm_tb()), +0.999141977, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(model.ckm_tb()), +0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(abs(model.ckm_tb()),  +0.999141,    1e-6);

                // angles
                double alpha = arg(-1.0 * model.ckm_td() * conj(model.ckm_tb()) / model.ckm_ud() / conj(model.ckm_ub()));
                TEST_CHECK_NEARLY_EQUAL(alpha,                                  +1.589220699, eps);
                TEST_CHECK_NEARLY_EQUAL(std::sin(2.0 * alpha),                  -0.036840406, eps);

                double beta  = arg(-1.0 * model.ckm_cd() * conj(model.ckm_cb()) / model.ckm_td() / conj(model.ckm_tb()));
                TEST_CHECK_NEARLY_EQUAL(beta,                                   -2.761464006, eps);
                TEST_CHECK_NEARLY_EQUAL(std::sin(2.0 * beta),                   +0.689107918, eps);

                double gamma = arg(-1.0 * model.ckm_ud() * conj(model.ckm_ub()) / model.ckm_cd() / conj(model.ckm_cb()));
                TEST_CHECK_NEARLY_EQUAL(gamma,                                  -1.969349346, eps);
                TEST_CHECK_NEARLY_EQUAL(std::abs(std::sin(2.0 * beta + gamma)), +0.935295092, eps);

                complex<double> lambda_t = model.ckm_tb() * conj(model.ckm_ts());
                TEST_CHECK_NEARLY_EQUAL(std::abs(lambda_t), +0.040483577, eps);
                complex<double> lambda_c = model.ckm_cb() * conj(model.ckm_cs());
                TEST_CHECK_NEARLY_EQUAL(std::abs(lambda_c), +0.040167570, eps);
                complex<double> lambda_u = model.ckm_ub() * conj(model.ckm_us());
                TEST_CHECK_NEARLY_EQUAL(std::abs(lambda_u), +0.000798232, eps);

                // unitarity
                TEST_CHECK_NEARLY_EQUAL(real(lambda_t + lambda_c + lambda_u), -1.131956683e-8, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(lambda_t + lambda_c + lambda_u), +0.0,            eps);
            }
        }
} ckm_elements_test;

class WilsonCoefficientsBToSTest :
    public TestCase
{
    public:
        WilsonCoefficientsBToSTest() :
            TestCase("wilson_coefficients_b_to_s_test")
        {
        }

        virtual void run() const
        {
            /* Test for 5 active flavors, evolving from mu_0c = 80, mu_0t = 120 to mu = 4.350516515 */
            {
                static const double eps = 1e-4;
                static const double mu = 4.350516515; // Stems from older, lower-order calculations of alpha_s

                Parameters parameters = reference_parameters();
                parameters["sb::mu"] = mu;
                StandardModel model(parameters);

                TEST_CHECK_NEARLY_EQUAL(model.alpha_s(mu), +0.2209967815, eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(mu, "mu", false);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c1()), -0.279801085, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c2()), +1.009683640, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c3()), -0.005775920, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c4()), -0.083977609, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c5()), +0.000401406, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c6()), +0.001072008, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c7()), -0.334390556, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c8()), -0.180952245, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c9()), +4.256827890, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c10()),-4.160202020, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10()), +0.0, eps);
            }

            /* Test for 5 active flavors, evolving from mu_0c = 80, mu_0t = 120 to mu = 4.2 */
            {
                static const double eps = 1e-4;
                static const double mu = 4.2; // approximate m_b(m_b) MSbar mass

                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                TEST_CHECK_NEARLY_EQUAL(model.alpha_s(mu), +0.2233419372, eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(mu, "mu", false);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c1()), -0.28768333, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c2()), +1.01013250, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c3()), -0.00600697, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c4()), -0.08597076, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c5()), +0.00041824, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c6()), +0.00112410, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c7()), -0.33613067, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c8()), -0.18205267, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c9()), +4.27450580, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.c10()),-4.16020202, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9()),  +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10()), +0.0, eps);
            }

            /* Test for equality between SM Wilson coefficients and default parameter values */
            {
                // Do NOT use the reference parameters here!
                static const double eps = 1e-4;
                static const double mu = 4.2; // approximate m_b(m_b) MSbar mass

                Parameters parameters = Parameters::Defaults();
                parameters["sb::mu"] = mu;
                StandardModel model(parameters);

                TEST_CHECK_NEARLY_EQUAL(model.alpha_s(mu), +0.2263282172, eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(mu, "mu", false);
                TEST_CHECK_RELATIVE_ERROR(parameters["b->s::c1"],           real(wc.c1()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["b->s::c2"],           real(wc.c2()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["b->s::c3"],           real(wc.c3()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["b->s::c4"],           real(wc.c4()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["b->s::c5"],           real(wc.c5()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["b->s::c6"],           real(wc.c6()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["b->s::Re{c7}"],       real(wc.c7()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["b->s::c8"],           real(wc.c8()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["b->smumu::Re{c9}"],   real(wc.c9()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["b->smumu::Re{c10}"],  real(wc.c10()), eps);

                TEST_CHECK_NEARLY_EQUAL(parameters["b->s::Im{c7}"],         imag(wc.c7()),  eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["b->smumu::Im{c9}"],     imag(wc.c9()),  eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["b->smumu::Im{c10}"],    imag(wc.c10()), eps);
            }
        }
} wilson_coefficients_b_to_s_test;

class WilsonCoefficientsSBSBTest :
    public TestCase
{
    public:
        WilsonCoefficientsSBSBTest() :
            TestCase("wilson_coefficients_sbsb_test")
        {
        }

        virtual void run() const
        {
            /* Test for 5 active flavors, evolving from mu_0 = 120 GeV to mu = 4.2 GeV */
            {
                static const double eps = 1e-8;

                Parameters parameters = reference_parameters(); // set scale sbsb::mu
                StandardModel model(parameters);

                WilsonCoefficients<wc::SBSB> wc = model.wet_sbsb();
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1()),  +0.001313228, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1()),   0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2()),   0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2()),   0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3()),   0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3()),   0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c4()),   0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4()),   0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c5()),   0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5()),   0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1p()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1p()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2p()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2p()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3p()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3p()),  0.000000000, eps);
            }
        }
} wilson_coefficients_sbsb_test;

class WilsonCoefficientsSBNuNuTest :
    public TestCase
{
    public:
        WilsonCoefficientsSBNuNuTest() :
            TestCase("wilson_coefficients_sbnunu_test")
        {
        }

        virtual void run() const
        {
            /* Test for 5 active flavors, evolving from mu_0 = 120 GeV to mu = 4.2 GeV */
            {
                static const double eps = 1e-8;

                Parameters parameters = reference_parameters(); // set scale sbnunu::mu
                StandardModel model(parameters);

                WilsonCoefficients<wc::SBNuNu> wc = model.wet_sbnunu(false);
                TEST_CHECK_NEARLY_EQUAL(real(wc.cVL()),  6.605426281, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.cVL()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.cVR()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.cVR()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.cSL()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.cSL()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.cSR()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.cSR()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.cTL()),  0.000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.cTL()),  0.000000000, eps);
            }
        }
} wilson_coefficients_sbnunu_test;

class WilsonCoefficientsSBCUTest :
    public TestCase
{
    public:
        WilsonCoefficientsSBCUTest() :
            TestCase("wilson_coefficients_sbcu_test")
        {
        }

        virtual void run() const
        {
            /* Test for 5 active flavors, evolving from mu_0 = 80 GeV to mu = 4.2 GeV */
            {
                static const double eps = 1e-8;

                Parameters parameters = reference_parameters(); // set scale sbcu::mu
                StandardModel model(parameters);

                TEST_CHECK_NEARLY_EQUAL( model.alpha_s(80.0), 0.12001051, 1e-6);
                TEST_CHECK_NEARLY_EQUAL( model.alpha_s( 4.2), 0.22334194, 1e-6);

                WilsonCoefficients<wc::SBCU> wc = model.wet_sbcu(false);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1()),   -0.041858794, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2()),   -0.896743838, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3()),    0.011274504, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c4()),    0.194524251, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c5()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c6()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c7()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c8()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c9()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c10()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c4p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c5p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c6p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c7p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c8p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c9p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c10p()),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10p()),  0.0,         eps);
            }
        }
} wilson_coefficients_sbcu_test;

class WilsonCoefficientsDBCUTest :
    public TestCase
{
    public:
        WilsonCoefficientsDBCUTest() :
            TestCase("wilson_coefficients_dbcu_test")
        {
        }

        virtual void run() const
        {
            /* Test for 5 active flavors, evolving from mu_0 = 80 GeV to mu = 4.2 GeV */
            {
                static const double eps = 1e-8;

                Parameters parameters = reference_parameters(); // set scale sbcu::mu
                StandardModel model(parameters);

                TEST_CHECK_NEARLY_EQUAL( model.alpha_s(80.0), 0.12001051, 1e-6);
                TEST_CHECK_NEARLY_EQUAL( model.alpha_s( 4.2), 0.22334194, 1e-6);

                WilsonCoefficients<wc::DBCU> wc = model.wet_dbcu(false);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1()),   -0.041858794, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2()),   -0.896743838, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3()),    0.011274504, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c4()),    0.194524251, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c5()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c6()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c7()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c8()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c9()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9()),    0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c10()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c4p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c5p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c6p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c7p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c8p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c9p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9p()),   0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c10p()),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10p()),  0.0,         eps);
            }
        }
} wilson_coefficients_dbcu_test;
