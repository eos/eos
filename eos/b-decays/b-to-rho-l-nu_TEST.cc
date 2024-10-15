/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
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
#include <eos/b-decays/b-to-vec-l-nu.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

using namespace test;
using namespace eos;

class BToRhoLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToRhoLeptonNeutrinoTest() :
            TestCase("b_to_rho_l_nu_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();
                p["CKM::abs(V_ub)"]              =  3.32e-3;
                p["B->rho::alpha^A0_0@BSZ2015" ] = +0.36;
                p["B->rho::alpha^A0_1@BSZ2015" ] = -0.83;
                p["B->rho::alpha^A0_2@BSZ2015" ] = +1.33;
                p["B->rho::alpha^A1_0@BSZ2015" ] = +0.26;
                p["B->rho::alpha^A1_1@BSZ2015" ] = +0.39;
                p["B->rho::alpha^A1_2@BSZ2015" ] = +0.16;
                p["B->rho::alpha^A12_1@BSZ2015"] = +0.76;
                p["B->rho::alpha^A12_2@BSZ2015"] = +0.46;
                p["B->rho::alpha^V_0@BSZ2015"  ] = +0.33;
                p["B->rho::alpha^V_1@BSZ2015"  ] = -0.86;
                p["B->rho::alpha^V_2@BSZ2015"  ] = +1.80;
                p["B->rho::alpha^T1_0@BSZ2015" ] = +0.27;
                p["B->rho::alpha^T1_1@BSZ2015" ] = -0.74;
                p["B->rho::alpha^T1_2@BSZ2015" ] = +1.45;
                p["B->rho::alpha^T2_1@BSZ2015" ] = +0.47;
                p["B->rho::alpha^T2_2@BSZ2015" ] = +0.58;
                p["B->rho::alpha^T23_0@BSZ2015"] = +0.75;
                p["B->rho::alpha^T23_1@BSZ2015"] = +1.90;
                p["B->rho::alpha^T23_2@BSZ2015"] = +2.93;

                Options ooplus
                {
                    { "model",        "CKM"     },
                    { "form-factors", "BSZ2015" },
                    { "V",            "rho"     },
                    { "q",            "d"       },
                    { "l",            "e"       },
                };

                BToVectorLeptonNeutrino dplus(p, ooplus);

                const double eps = 1e-5;

                TEST_CHECK_RELATIVE_ERROR(dplus.integrated_branching_ratio( 0.01,  2.00), 2.27193e-05, eps);
                TEST_CHECK_RELATIVE_ERROR(dplus.integrated_branching_ratio( 2.00,  4.00), 2.62836e-05, eps);
                TEST_CHECK_RELATIVE_ERROR(dplus.integrated_branching_ratio( 4.00,  6.00), 2.95276e-05, eps);
                TEST_CHECK_RELATIVE_ERROR(dplus.integrated_branching_ratio( 6.00,  8.00), 3.24542e-05, eps);
                TEST_CHECK_RELATIVE_ERROR(dplus.integrated_branching_ratio( 8.00, 10.00), 3.49032e-05, eps);
                TEST_CHECK_RELATIVE_ERROR(dplus.integrated_branching_ratio(10.00, 12.00), 3.66198e-05, eps);

                Options oozero
                {
                    { "model",        "CKM"     },
                    { "form-factors", "BSZ2015" },
                    { "V",            "rho"     },
                    { "q",            "u"       },
                    { "l",            "e"       },
                };

                BToVectorLeptonNeutrino dzero(p, oozero);

                TEST_CHECK_RELATIVE_ERROR(dzero.integrated_branching_ratio( 0.01,  2.00), 1.22754e-05, eps);
                TEST_CHECK_RELATIVE_ERROR(dzero.integrated_branching_ratio( 2.00,  4.00), 1.42002e-05, eps);
                TEST_CHECK_RELATIVE_ERROR(dzero.integrated_branching_ratio( 4.00,  6.00), 1.59517e-05, eps);
                TEST_CHECK_RELATIVE_ERROR(dzero.integrated_branching_ratio( 6.00,  8.00), 1.75316e-05, eps);
                TEST_CHECK_RELATIVE_ERROR(dzero.integrated_branching_ratio( 8.00, 10.00), 1.88532e-05, eps);
                TEST_CHECK_RELATIVE_ERROR(dzero.integrated_branching_ratio(10.00, 12.00), 1.97787e-05, eps);

            }
        }
} b_to_rho_l_nu_test;
