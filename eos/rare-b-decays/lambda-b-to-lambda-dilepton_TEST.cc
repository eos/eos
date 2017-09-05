/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Thomas Blake
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
#include <eos/observable.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-dilepton.hh>
#include <eos/utils/complex.hh>

using namespace test;
using namespace eos;

class LambdaBToLambdaDileptonLowRecoilTest :
    public TestCase
{
    public:
        LambdaBToLambdaDileptonLowRecoilTest() :
            TestCase("lambda_b_to_lambda_dilepton_test")
        {
        }

        virtual void run() const
        {
            // Standard Model
            {
                static const double eps = 1e-4;

                // unpolarised SM
                {
                    Options oo;
                    oo.set("model",                   "WilsonScan");
                    oo.set("q",                       "d");
                    oo.set("l",                       "mu");
                    oo.set("production-polarisation", "unpolarised");

                    Parameters p = Parameters::Defaults();

                    LambdaBToLambdaDilepton<LowRecoil> d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(16.0), +1.0923e-7, eps);
		   
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m1(15.0,19.0),   0.3536, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m2(15.0,19.0),   0.2928, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m3(15.0,19.0),  -0.2451, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m4(15.0,19.0),  -0.2055, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m5(15.0,19.0),  -0.1604, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m6(15.0,19.0),   0.1842, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m7(15.0,19.0),  -0.0228, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m8(15.0,19.0),  -0.0888, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m9(15.0,19.0),   0.0004, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m10(15.0,19.0), -0.0006, eps);

		    TEST_CHECK_EQUAL(d.integrated_m11(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m12(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m13(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m14(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m15(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m16(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m17(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m18(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m19(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m20(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m21(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m22(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m23(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m24(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m25(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m26(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m27(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m28(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m29(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m30(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m31(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m32(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m33(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m34(15.0,19.0), 0);
 
                }

                // LHCb-polarised SM
                {
                    Options oo;
                    oo.set("model",                   "WilsonScan");
                    oo.set("q",                       "d");
                    oo.set("l",                       "mu");
                    oo.set("production-polarisation", "LHCb");

                    Parameters p = Parameters::Defaults();

                    LambdaBToLambdaDilepton<LowRecoil> d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(16.0), +1.0923e-7, eps);

		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m1(15.0,19.0),   0.3536, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m2(15.0,19.0),   0.2928, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m3(15.0,19.0),  -0.2451, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m4(15.0,19.0),  -0.2055, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m5(15.0,19.0),  -0.1604, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m6(15.0,19.0),   0.1842, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m7(15.0,19.0),  -0.0228, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m8(15.0,19.0),  -0.0888, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m9(15.0,19.0),   0.0004, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m10(15.0,19.0), -0.0006, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m11(15.0,19.0), -0.0042, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m12(15.0,19.0),  0.0150, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m13(15.0,19.0), -0.0172, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m14(15.0,19.0),  0.0023, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m15(15.0,19.0), -0.0113, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m16(15.0,19.0),  0.0094, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m17(15.0,19.0),  0.0054, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m18(15.0,19.0),  0.0013, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m19(15.0,19.0),  0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m20(15.0,19.0), -0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m21(15.0,19.0),  0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m22(15.0,19.0), -0.0001, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m23(15.0,19.0), -0.0188, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m24(15.0,19.0),  0.0203, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m25(15.0,19.0), -0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m26(15.0,19.0),  0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m27(15.0,19.0),  0.0133, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m28(15.0,19.0), -0.0118, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m29(15.0,19.0),  0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m30(15.0,19.0), -0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m31(15.0,19.0),  0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m32(15.0,19.0), -0.0024, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m33(15.0,19.0), -0.0028, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m34(15.0,19.0),  0.0000, eps);
		}

                // unpolarised BMP
                {
                    Options oo;
                    oo.set("model",                   "WilsonScan");
                    oo.set("q",                       "d");
                    oo.set("l",                       "mu");
                    oo.set("production-polarisation", "unpolarised");

                    Parameters p = Parameters::Defaults();
                    p["b->smumu::Re{c9}"]  = +3.2734;
                    p["b->smumu::Re{c9'}"] = +1.0000;

                    LambdaBToLambdaDilepton<LowRecoil> d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(16.0), +0.8251e-7, eps);

		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m1(15.0,19.0), 0.3567, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m2(15.0,19.0), 0.2867, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m3(15.0,19.0), -0.2639, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m4(15.0,19.0), -0.2131, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m5(15.0,19.0), -0.1674, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m6(15.0,19.0), 0.1730, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m7(15.0,19.0), -0.0225, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m8(15.0,19.0), -0.0351, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m9(15.0,19.0), 0.0006, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m10(15.0,19.0), -0.0008, eps);
		    
		    TEST_CHECK_EQUAL(d.integrated_m11(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m12(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m13(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m14(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m15(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m16(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m17(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m18(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m19(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m20(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m21(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m22(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m23(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m24(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m25(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m26(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m27(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m28(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m29(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m30(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m31(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m32(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m33(15.0,19.0), 0);
		    TEST_CHECK_EQUAL(d.integrated_m34(15.0,19.0), 0);
                }

                // LHCb-polarised BMP
                {
                    Options oo;
                    oo.set("model",                   "WilsonScan");
                    oo.set("q",                       "d");
                    oo.set("l",                       "mu");
                    oo.set("production-polarisation", "LHCb");

                    Parameters p = Parameters::Defaults();
                    p["b->smumu::Re{c9}"]  = +3.2734;
                    p["b->smumu::Re{c9'}"] = +1.0000;

                    LambdaBToLambdaDilepton<LowRecoil> d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(16.0), +0.8251e-7, eps);

		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m1(15.0,19.0),   0.3567, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m2(15.0,19.0),   0.2867, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m3(15.0,19.0),  -0.2639, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m4(15.0,19.0),  -0.2131, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m5(15.0,19.0),  -0.1674, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m6(15.0,19.0),   0.1730, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m7(15.0,19.0),  -0.0225, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m8(15.0,19.0),  -0.0351, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m9(15.0,19.0),   0.0006, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m10(15.0,19.0), -0.0008, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m11(15.0,19.0), -0.0043, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m12(15.0,19.0),  0.0156, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m13(15.0,19.0), -0.0162, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m14(15.0,19.0),  0.0027, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m15(15.0,19.0), -0.0110, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m16(15.0,19.0),  0.0102, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m17(15.0,19.0),  0.0026, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m18(15.0,19.0),  0.0015, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m19(15.0,19.0), -0.0005, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m20(15.0,19.0), -0.0001, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m21(15.0,19.0),  0.0001, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m22(15.0,19.0), -0.0002, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m23(15.0,19.0), -0.0195, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m24(15.0,19.0),  0.0196, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m25(15.0,19.0),  0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m26(15.0,19.0),  0.0001, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m27(15.0,19.0),  0.0134, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m28(15.0,19.0), -0.0128, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m29(15.0,19.0),  0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m30(15.0,19.0),  0.0003, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m31(15.0,19.0),  0.0000, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m32(15.0,19.0), -0.0006, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m33(15.0,19.0), -0.0017, eps);
		    TEST_CHECK_NEARLY_EQUAL(d.integrated_m34(15.0,19.0),  0.0002, eps);   
                }
            }
        }
} lambda_b_to_lambda_dilepton_low_recoil_test;
