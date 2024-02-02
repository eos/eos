/*
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2022 Méril Reboud
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
#include <eos/form-factors/parametric-dkmr2017.hh>
#include <eos/form-factors/parametric-dkmr2017-impl.hh>

using namespace test;
using namespace eos;

class DKMR2017FormFactorsTest :
public TestCase
{
public:
    DKMR2017FormFactorsTest() :
    TestCase("dkmr2017_form_factors_test")
    {
    }

    virtual void run() const
    {
        static const double eps = 1e-3;

        Parameters p = Parameters::Defaults();
        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> ff = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda_c::DKMR2017", p);
        // fix z-expansion parameters @DKMR2017 - all a_2's are zero by default
        // Note that a0_perp_A = a0_long_A ve a0_perp_T5 = a0_long_T5
        p["Lambda_b->Lambda_c::a_0_time^V@DKMR2017"]    =  0.33;
        p["Lambda_b->Lambda_c::a_1_time^V@DKMR2017"]    = -0.43;
        p["Lambda_b->Lambda_c::a_2_time^V@DKMR2017"]    = +0.00;

        p["Lambda_b->Lambda_c::a_0_long^V@DKMR2017"]    =  0.44;
        p["Lambda_b->Lambda_c::a_1_long^V@DKMR2017"]    =  0.27;
        p["Lambda_b->Lambda_c::a_2_long^V@DKMR2017"]    =  0.00;

        p["Lambda_b->Lambda_c::a_0_perp^V@DKMR2017"]    = -0.52;
        p["Lambda_b->Lambda_c::a_1_perp^V@DKMR2017"]    = +0.42;
        p["Lambda_b->Lambda_c::a_2_perp^V@DKMR2017"]    = +0.00;

        p["Lambda_b->Lambda_c::a_0_time^A@DKMR2017"]    = +1.75;
        p["Lambda_b->Lambda_c::a_1_time^A@DKMR2017"]    = -1.1;
        p["Lambda_b->Lambda_c::a_2_time^A@DKMR2017"]    = -0.00;

        p["Lambda_b->Lambda_c::a_0_long^A@DKMR2017"]    =  2.14;
        p["Lambda_b->Lambda_c::a_1_long^A@DKMR2017"]    = -1.58;
        p["Lambda_b->Lambda_c::a_2_long^A@DKMR2017"]    = -0.00;
        //
        p["Lambda_b->Lambda_c::a_1_perp^A@DKMR2017"]    = -1.30;
        p["Lambda_b->Lambda_c::a_2_perp^A@DKMR2017"]    = -0.00;

        p["Lambda_b->Lambda_c::a_0_long^T@DKMR2017"]    =  0.71;
        p["Lambda_b->Lambda_c::a_1_long^T@DKMR2017"]    =  0.21;
        p["Lambda_b->Lambda_c::a_2_long^T@DKMR2017"]    =  0.00;

        p["Lambda_b->Lambda_c::a_0_perp^T@DKMR2017"]    = -0.24;
        p["Lambda_b->Lambda_c::a_1_perp^T@DKMR2017"]    = -0.24;
        p["Lambda_b->Lambda_c::a_2_perp^T@DKMR2017"]    = -0.00;

        p["Lambda_b->Lambda_c::a_0_long^T5@DKMR2017"]   =  0.14;
        p["Lambda_b->Lambda_c::a_1_long^T5@DKMR2017"]   = -0.14;
        p["Lambda_b->Lambda_c::a_2_long^T5@DKMR2017"]   = -0.00;
        //
        p["Lambda_b->Lambda_c::a_1_perp^T5@DKMR2017"]   =  3.14;
        p["Lambda_b->Lambda_c::a_2_perp^T5@DKMR2017"]   =  0.00;

        p["mass::Lambda_b"]                             = 5.61951;
        p["mass::Lambda_c"]                             = 2.2865;

        // test of the ten baryonic FFs
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_v( 0.0),  0.299748,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_v( 5.0),  0.351122,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_v(10.0),  0.419267,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_v(15.0),  0.513241,  eps);

        TEST_CHECK_NEARLY_EQUAL(ff->f_long_v( 0.0),  0.461852,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_v( 5.0),  0.517426,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_v(10.0),  0.589584,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_v(15.0),  0.687469,  eps);

        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_v( 0.0),  -0.486008, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_v( 5.0),  -0.571162, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_v(10.0),  -0.687539, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_v(15.0),  -0.855001, eps);

        TEST_CHECK_NEARLY_EQUAL(ff->f_time_a( 0.0),  1.65909,   eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_a( 5.0),  1.94289,   eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_a(10.0),  2.3313,    eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_a(15.0),  2.89206,   eps);

        TEST_CHECK_NEARLY_EQUAL(ff->f_long_a( 0.0),  2.03046,   eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_a( 5.0),  2.33035,   eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_a(10.0),  2.72177,   eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_a(15.0),  3.25185,   eps);

        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_a( 0.0),  2.04987,   eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_a( 5.0),  2.34308,   eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_a(10.0),  2.72459,   eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_a(15.0),  3.23947,   eps);

        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t( 0.0),  0.726996,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t( 5.0),  0.822619,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t(10.0),  0.948552,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t(15.0),  1.122315,  eps);

        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t( 0.0),  -0.259424, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t( 5.0),  -0.287293, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t(10.0),  -0.322751, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t(15.0),  -0.369644, eps);

        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t5( 0.0), 0.130294,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t5( 5.0), 0.150786,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t5(10.0), 0.177691,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t5(15.0), 0.214357,  eps);

        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t5( 0.0), 0.357693,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t5( 5.0), 0.299983,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t5(10.0), 0.210694,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t5(15.0), 0.069375,  eps);
    }
} dkmr2017_form_factors_test;
