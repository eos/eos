/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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
#include <src/utils/model.hh>
#include <src/utils/standard-model.hh>

#include <cmath>
#include <iostream>

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
    result["mass::Z"] = 91.1876;

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
                std::shared_ptr<Model> m = Model::make("SM", reference_parameters());
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

            // Data in agreement with RanDec, cf. [CKS2000]
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

            StandardModel model(reference_parameters());

            TEST_CHECK_NEARLY_EQUAL(3.67956, model.m_b_msbar(9.6), eps);
            TEST_CHECK_NEARLY_EQUAL(4.10051, model.m_b_msbar(4.8), eps);
            TEST_CHECK_NEARLY_EQUAL(4.20000, model.m_b_msbar(4.2), eps);
            TEST_CHECK_NEARLY_EQUAL(4.75221, model.m_b_msbar(2.4), eps);

            TEST_CHECK_NEARLY_EQUAL(4.88402, model.m_b_pole(), eps);

            TEST_CHECK_NEARLY_EQUAL(4.60728, model.m_b_ps(1.0), eps);
            TEST_CHECK_NEARLY_EQUAL(4.54012, model.m_b_ps(1.5), eps);
            TEST_CHECK_NEARLY_EQUAL(4.47735, model.m_b_ps(2.0), eps);
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

            StandardModel model(reference_parameters());

            TEST_CHECK_NEARLY_EQUAL(0.891000, model.m_c_msbar(4.8),  eps);
            TEST_CHECK_NEARLY_EQUAL(0.912618, model.m_c_msbar(4.2),  eps);
            TEST_CHECK_NEARLY_EQUAL(1.270000, model.m_c_msbar(1.27), eps);

            TEST_CHECK_NEARLY_EQUAL(1.891359, model.m_c_pole(), eps);
        }
} sm_c_masses_test;

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
            static const double eps = 1e-6;

            // central values
            {
                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                TEST_CHECK_NEARLY_EQUAL(+0.225429000, real(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_us()), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.040511700, real(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000734841, imag(model.ckm_ts()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.001371630, real(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.003264270, imag(model.ckm_ub()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.041264500, real(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_cb()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.999142000, real(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_tb()), eps);
            }

            // A raised
            {
                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                parameters["CKM::A"] = parameters["CKM::A"].max();

                TEST_CHECK_NEARLY_EQUAL(+0.225429000, real(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_us()), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.041160200, real(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000746979, imag(model.ckm_ts()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.001394070, real(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.003316490, imag(model.ckm_ub()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.041925100, real(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_cb()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.999114000, real(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_tb()), eps);
            }

            // A lowered
            {
                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                parameters["CKM::A"] = parameters["CKM::A"].min();

                TEST_CHECK_NEARLY_EQUAL(+0.225429000, real(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_us()), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.039164700, real(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000710846, imag(model.ckm_ts()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.001326740, real(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.003155800, imag(model.ckm_ub()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.039892400, real(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_cb()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.999198000, real(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_tb()), eps);
            }

            // lambda raised
            {
                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                parameters["CKM::lambda"] = parameters["CKM::lambda"].max();

                TEST_CHECK_NEARLY_EQUAL(+0.226199000, real(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_us()), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.040783700, real(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000745458, imag(model.ckm_ts()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.001386510, real(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.003298420, imag(model.ckm_ub()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.041546900, real(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_cb()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.999130000, real(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_tb()), eps);
            }

            // lambda lowered
            {
                Parameters parameters = reference_parameters();
                StandardModel model(parameters);

                parameters["CKM::lambda"] = parameters["CKM::lambda"].min();

                TEST_CHECK_NEARLY_EQUAL(+0.224659000, real(model.ckm_us()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_us()), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.040240500, real(model.ckm_ts()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000725122, imag(model.ckm_ts()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.001357970, real(model.ckm_ub()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.003230360, imag(model.ckm_ub()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.040983100, real(model.ckm_cb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_cb()), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.999154000, real(model.ckm_tb()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000, imag(model.ckm_tb()), eps);
            }
        }
} ckm_elements_test;
