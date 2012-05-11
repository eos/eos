/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013 Danny van Dyk
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
#include <eos/utils/model.hh>
#include <eos/utils/standard-model.hh>

#include <cmath>

#include <iostream> // <-- Remove!

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
    result["mass::b(MSbar)"] = 4.2;
    result["mass::W"] = 80.398;
    result["mass::Z"] = 91.1876;
    result["mass::t(pole)"] = 173.3;

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

            TEST_CHECK_NEARLY_EQUAL(4.74167, model.m_b_pole(), eps);

            TEST_CHECK_NEARLY_EQUAL(4.60728, model.m_b_ps(1.0), eps);
            TEST_CHECK_NEARLY_EQUAL(4.54012, model.m_b_ps(1.5), eps);
            TEST_CHECK_NEARLY_EQUAL(4.47735, model.m_b_ps(2.0), eps);

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

            TEST_CHECK_NEARLY_EQUAL(0.101000, model.m_s_msbar(2.0),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.084980, model.m_s_msbar(4.2),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.082967, model.m_s_msbar(4.8),  eps);
        }
} sm_s_masses_test;

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

            // A raised
            {
                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                parameters["CKM::A"] = parameters["CKM::A"].max();

                TEST_CHECK_NEARLY_EQUAL(+0.225428541, real(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_us()), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.041160228, real(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.000746979, imag(model.ckm_ts()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.001394070, real(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.003316490, imag(model.ckm_ub()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.041925143, real(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_cb()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.999114272, real(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_tb()), eps);
            }

            // A lowered
            {
                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                parameters["CKM::A"] = parameters["CKM::A"].min();

                TEST_CHECK_NEARLY_EQUAL(+0.225428679, real(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_us()), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.039164663, real(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.000710846, imag(model.ckm_ts()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.001326740, real(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.003155800, imag(model.ckm_ub()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.039892433, real(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_cb()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.999198111, real(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_tb()), eps);
            }

            // lambda raised
            {
                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                parameters["CKM::lambda"] = parameters["CKM::lambda"].max();

                TEST_CHECK_NEARLY_EQUAL(+0.226198552, real(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_us()), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.040783654, real(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.000745458, imag(model.ckm_ts()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.001386510, real(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.003298420, imag(model.ckm_ub()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.041546883, real(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_cb()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.999130143, real(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_tb()), eps);
            }

            // lambda lowered
            {
                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                parameters["CKM::lambda"] = parameters["CKM::lambda"].min();

                TEST_CHECK_NEARLY_EQUAL(+0.224658620, real(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_us()), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.040240544, real(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.000725122, imag(model.ckm_ts()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.001357970, real(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.003230360, imag(model.ckm_ub()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.040983100, real(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_cb()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.999153689, real(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_tb()), eps);
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
                StandardModel model(parameters);

                parameters["mu"] = mu;
                TEST_CHECK_NEARLY_EQUAL(+0.2209967815, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(false);
                TEST_CHECK_RELATIVE_ERROR(-0.28058190, real(wc.c1()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+1.00972828, real(wc.c2()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.00581526, real(wc.c3()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.08408026, real(wc.c4()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+0.00040465, real(wc.c5()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+0.00108870, real(wc.c6()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.32561211, real(wc.c7()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.17593283, real(wc.c8()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+4.25821127, real(wc.c9()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-4.15077942, real(wc.c10()), eps);
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

                parameters["mu"] = mu;
                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(false);
                TEST_CHECK_RELATIVE_ERROR(-0.28846675, real(wc.c1()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+1.01017822, real(wc.c2()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.00604628, real(wc.c3()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.08607506, real(wc.c4()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+0.00042146, real(wc.c5()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+0.00114089, real(wc.c6()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.32741917, real(wc.c7()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-0.17707354, real(wc.c8()),  eps);
                TEST_CHECK_RELATIVE_ERROR(+4.27584793, real(wc.c9()),  eps);
                TEST_CHECK_RELATIVE_ERROR(-4.15077943, real(wc.c10()), eps);
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
                StandardModel model(parameters);

                parameters["mu"] = mu;
                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(false);
                TEST_CHECK_RELATIVE_ERROR(parameters["c1"],       real(wc.c1()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["c2"],       real(wc.c2()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["c3"],       real(wc.c3()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["c4"],       real(wc.c4()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["c5"],       real(wc.c5()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["c6"],       real(wc.c6()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["Abs{c7}"],  abs(wc.c7()),   eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["Arg{c7}"],  arg(wc.c7()),   eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["Re{c7}"],   real(wc.c7()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["c8"],       real(wc.c8()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["Abs{c9}"],  abs(wc.c9()),   eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["Re{c9}"],   real(wc.c9()),  eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["Abs{c10}"], abs(wc.c10()),   eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["Arg{c10}"], arg(wc.c10()),   eps);
                TEST_CHECK_RELATIVE_ERROR(parameters["Re{c10}"],  real(wc.c10()),  eps);

                TEST_CHECK_NEARLY_EQUAL(parameters["Im{c7}"],   imag(wc.c7()),  eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["Arg{c9}"],  arg(wc.c9()),   eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["Im{c9}"],   imag(wc.c9()),  eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["Im{c10}"],  imag(wc.c10()),  eps);
            }
        }
} wilson_coefficients_b_to_s_test;
