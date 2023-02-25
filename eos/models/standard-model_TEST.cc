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
            TEST_CHECK_NEARLY_EQUAL(0.117620, model.alpha_s(91.1876), 5e-5);

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

            TEST_CHECK_RELATIVE_ERROR(167.794, model.m_t_msbar(120.0), eps);
            TEST_CHECK_RELATIVE_ERROR(173.647, model.m_t_msbar( 80.0), eps);
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

            TEST_CHECK_NEARLY_EQUAL(3.67956, model.m_b_msbar(9.6), eps);
            TEST_CHECK_NEARLY_EQUAL(4.10051, model.m_b_msbar(4.8), eps);
            TEST_CHECK_NEARLY_EQUAL(4.20000, model.m_b_msbar(4.2), eps);
            TEST_CHECK_NEARLY_EQUAL(4.75221, model.m_b_msbar(2.4), eps);

            TEST_CHECK_NEARLY_EQUAL(4.74167, model.m_b_pole(), 1e-3); // Precision is hard-limited in fixed-point routine

            TEST_CHECK_NEARLY_EQUAL(4.60728, model.m_b_ps(1.0), eps);
            TEST_CHECK_NEARLY_EQUAL(4.54012, model.m_b_ps(1.5), eps);
            TEST_CHECK_NEARLY_EQUAL(4.47735, model.m_b_ps(2.0), eps);

            TEST_CHECK_NEARLY_EQUAL(4.63362, model.m_b_kin(0.75), eps);
            TEST_CHECK_NEARLY_EQUAL(4.56114, model.m_b_kin(1.00), eps);
            TEST_CHECK_NEARLY_EQUAL(4.49203, model.m_b_kin(1.25), eps);
            TEST_CHECK_NEARLY_EQUAL(4.42520, model.m_b_kin(1.50), eps);
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

            TEST_CHECK_NEARLY_EQUAL(0.891000, model.m_c_msbar(4.8),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.912618, model.m_c_msbar(4.2),  eps);
            TEST_CHECK_NEARLY_EQUAL(1.270000, model.m_c_msbar(1.27), eps);

            TEST_CHECK_NEARLY_EQUAL(1.595301, model.m_c_pole(), eps);

            TEST_CHECK_NEARLY_EQUAL(1.203723, model.m_c_kin(0.75), eps);
            TEST_CHECK_NEARLY_EQUAL(1.060682, model.m_c_kin(1.00), eps);
            TEST_CHECK_NEARLY_EQUAL(0.931772, model.m_c_kin(1.25), eps);
            TEST_CHECK_NEARLY_EQUAL(0.813366, model.m_c_kin(1.50), eps);
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

            TEST_CHECK_NEARLY_EQUAL(0.136682, model.m_s_msbar(1.0),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.106128, model.m_s_msbar(1.7),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.101000, model.m_s_msbar(2.0),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.084980, model.m_s_msbar(4.2),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.082967, model.m_s_msbar(4.8),  eps);
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

            TEST_CHECK_NEARLY_EQUAL(0.010826, model.m_ud_msbar(1.0),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.008000, model.m_ud_msbar(2.0),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.007223, model.m_ud_msbar(3.0),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.006803, model.m_ud_msbar(4.0),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.006525, model.m_ud_msbar(5.0),  eps);
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
                TEST_CHECK_NEARLY_EQUAL(+0.974253267, real(model.ckm_ud()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_ud()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.97425,     abs(model.ckm_ud()), 1e-5);

                TEST_CHECK_NEARLY_EQUAL(+0.225428590, real(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.22543,     abs(model.ckm_us()), 1e-5);

                TEST_CHECK_NEARLY_EQUAL(+0.001372189, real(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.003264270, imag(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.00354,     abs(model.ckm_ub()), 1e-5);

                TEST_CHECK_NEARLY_EQUAL(+0.225296132, real(model.ckm_cd()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000138121, imag(model.ckm_cd()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.22529,     abs(model.ckm_cd()), 1e-5);

                TEST_CHECK_NEARLY_EQUAL(+0.973416767, real(model.ckm_cs()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.000030365, imag(model.ckm_cs()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.97342,     abs(model.ckm_cs()), 1e-5);

                TEST_CHECK_NEARLY_EQUAL(+0.041264513, real(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.04126,     abs(model.ckm_cb()), 1e-5);

                TEST_CHECK_NEARLY_EQUAL(+0.007966605, real(model.ckm_td()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.003177489, imag(model.ckm_td()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.00858,     abs(model.ckm_td()), 1e-5);

                TEST_CHECK_NEARLY_EQUAL(-0.040511671, real(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.000735237, imag(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.04052,     abs(model.ckm_ts()), 1e-5);

                TEST_CHECK_NEARLY_EQUAL(+0.999141977, real(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.999141,    abs(model.ckm_tb()), 1e-6);

                // angles
                double alpha = arg(-1.0 * model.ckm_td() * conj(model.ckm_tb()) / model.ckm_ud() / conj(model.ckm_ub()));
                TEST_CHECK_NEARLY_EQUAL(+1.589220699,    alpha, eps);
                TEST_CHECK_NEARLY_EQUAL(-0.036840406,    std::sin(2.0 * alpha), eps);

                double beta  = arg(-1.0 * model.ckm_cd() * conj(model.ckm_cb()) / model.ckm_td() / conj(model.ckm_tb()));
                TEST_CHECK_NEARLY_EQUAL(-2.761464006,    beta, eps);
                TEST_CHECK_NEARLY_EQUAL(+0.689107918,    std::sin(2.0 * beta), eps);

                double gamma = arg(-1.0 * model.ckm_ud() * conj(model.ckm_ub()) / model.ckm_cd() / conj(model.ckm_cb()));
                TEST_CHECK_NEARLY_EQUAL(-1.969349346,    gamma, eps);
                TEST_CHECK_NEARLY_EQUAL(+0.935295092,    std::abs(std::sin(2.0 * beta + gamma)), eps);

                complex<double> lambda_t = model.ckm_tb() * conj(model.ckm_ts());
                TEST_CHECK_NEARLY_EQUAL(+0.040483577,    std::abs(lambda_t), eps);
                complex<double> lambda_c = model.ckm_cb() * conj(model.ckm_cs());
                TEST_CHECK_NEARLY_EQUAL(+0.040167570,    std::abs(lambda_c), eps);
                complex<double> lambda_u = model.ckm_ub() * conj(model.ckm_us());
                TEST_CHECK_NEARLY_EQUAL(+0.000798232,    std::abs(lambda_u), eps);

                // unitarity
                TEST_CHECK_NEARLY_EQUAL(-1.131956683e-8, real(lambda_t + lambda_c + lambda_u), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,            imag(lambda_t + lambda_c + lambda_u), eps);
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

                TEST_CHECK_NEARLY_EQUAL(+0.2209967815, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(mu, "mu", false);
                TEST_CHECK_RELATIVE_ERROR(-0.279801085, real(wc.c1()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+1.009683640, real(wc.c2()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.005775920, real(wc.c3()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.083977609, real(wc.c4()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+0.000401406, real(wc.c5()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+0.001072008, real(wc.c6()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.334390556, real(wc.c7()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.180952245, real(wc.c8()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+4.256827890, real(wc.c9()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-4.160202020, real(wc.c10()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c1()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c2()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c3()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c4()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c5()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c6()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c7()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c8()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c9()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c10()),  eps);
            }

            /* Test for 5 active flavors, evolving from mu_0c = 80, mu_0t = 120 to mu = 4.2 */
            {
                static const double eps = 1e-4;
                static const double mu = 4.2; // approximate m_b(m_b) MSbar mass

                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(mu, "mu", false);
                TEST_CHECK_RELATIVE_ERROR(-0.28768333, real(wc.c1()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+1.01013250, real(wc.c2()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.00600697, real(wc.c3()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.08597076, real(wc.c4()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+0.00041824, real(wc.c5()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+0.00112410, real(wc.c6()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.33613067, real(wc.c7()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.18205267, real(wc.c8()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+4.27450580, real(wc.c9()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-4.16020202, real(wc.c10()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c1()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c2()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c3()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c4()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c5()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c6()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c7()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c8()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c9()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c10()),  eps);
            }

            /* Test for equality between SM Wilson coefficients and default parameter values */
            {
                // Do NOT use the reference parameters here!
                static const double eps = 1e-4;
                static const double mu = 4.2; // approximate m_b(m_b) MSbar mass

                Parameters parameters = Parameters::Defaults();
                parameters["sb::mu"] = mu;
                StandardModel model(parameters);

                TEST_CHECK_NEARLY_EQUAL(+0.2263282172, model.alpha_s(mu), eps);

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
                TEST_CHECK_NEARLY_EQUAL(+0.001313228, real(wc.c1()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.c1()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.c2()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.c2()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.c3()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.c3()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.c4()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.c4()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.c5()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.c5()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.c1p()), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.c1p()), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.c2p()), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.c2p()), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.c3p()), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.c3p()), eps);
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
                TEST_CHECK_NEARLY_EQUAL( 6.605426281, real(wc.cVL()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.cVL()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.cVR()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.cVR()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.cSL()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.cSL()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.cSR()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.cSR()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, real(wc.cTL()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.000000000, imag(wc.cTL()),  eps);
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
