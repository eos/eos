/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Stefan Meiser
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
#include <eos/b-decays/bq-to-dstarq-vec.hh>
#include <eos/maths/complex.hh>

using namespace test;
using namespace eos;

class BqToDstarqVecTest :
    public TestCase
{
    public:
        BqToDstarqVecTest() :
            TestCase("bq_to_dstarq_vec_test")
        {
        }

        virtual void run() const
        {
            // Standard Model NLO
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_d"] = 1.519e-12;
                // Meson masses
                p["mass::B_d"] = 5.27966;
                p["mass::D_d^*"] = 2.01026;
                p["mass::K_u^*"] = 0.89166;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                p["mass::c"] = 1.2;
                // decay constant
                p["decay-constant::K_u^*"] = 0.204;
                // fermi constant
                p["WET::G_Fermi"] = 1.16637e-05;
                // CKM matrix elements
                p["CKM::abs(V_cb)"] = 0.0408;
                p["CKM::arg(V_cb)"] = 0.0;
                p["CKM::abs(V_us)"] = 0.2243;
                p["CKM::arg(V_us)"] = 0.0;
                // form factor
                p["B->D^*K^*::a_0(MKstar2)"] = 0.70023;
                // WC Values
                p["sbcu::Re{c1}" ] = -0.04235657776117585;
                p["sbcu::Im{c1}" ] = 0.0;
                p["sbcu::Re{c2}" ] = -0.8948941708221622;
                p["sbcu::Im{c2}" ] = 0.0;
                p["sbcu::Re{c3}" ] = 0.011381250932999982;
                p["sbcu::Im{c3}" ] = 0.0;
                p["sbcu::Re{c4}" ] = 0.19426386543613433;
                p["sbcu::Im{c4}" ] = 0.0;
                Options oo
                {
                    { "accuracy",     "LO+NLO" },
                    { "q",            "d"      },
                    { "model",        "WET"    }
                };
                BqToDstarqVector d(p, oo);
                {
                    const double eps = 1.0e-6;
                    TEST_CHECK_RELATIVE_ERROR(d.re_a_1(), 1.0522133268488623,   eps);
                    TEST_CHECK_RELATIVE_ERROR(d.im_a_1(), 0.013052982861163655, eps);

                    TEST_CHECK_RELATIVE_ERROR(d.decay_width(),     2.2493931485120835e-16, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(), 0.0005191075508483368 , eps);
                }
            }

            // BSM NLO test case
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_d"] = 1.519e-12;
                // Meson masses
                p["mass::B_d"] = 5.27966;
                p["mass::D_d^*"] = 2.01026;
                p["mass::K_u^*"] = 0.89166;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                p["mass::c"] = 1.2;
                // decay constant
                p["decay-constant::K_u^*"] = 0.204;
                // fermi constant
                p["WET::G_Fermi"] = 1.16637e-05;
                // CKM matrix elements
                p["CKM::abs(V_cb)"] = 0.0408;
                p["CKM::arg(V_cb)"] = 0.0;
                p["CKM::abs(V_us)"] = 0.2243;
                p["CKM::arg(V_us)"] = 0.0;
                // form factor
                p["B->D^*K^*::a_0(MKstar2)"] = 0.70023;
                // BSM test case: WET parameter point
                p["sbcu::Re{c1}"  ] = -1.72424;
                p["sbcu::Im{c1}"  ] = -1.56379;
                p["sbcu::Re{c1'}" ] = -1.05356;
                p["sbcu::Im{c1'}" ] = -0.791464;
                p["sbcu::Re{c2}"  ] = -2.84324;
                p["sbcu::Im{c2}"  ] = -1.10401;
                p["sbcu::Re{c2'}" ] = +1.10235;
                p["sbcu::Im{c2'}" ] = +2.0774;
                p["sbcu::Re{c3}"  ] = +1.61473;
                p["sbcu::Im{c3}"  ] = +1.23153;
                p["sbcu::Re{c3'}" ] = -2.95587;
                p["sbcu::Im{c3'}" ] = -2.28859;
                p["sbcu::Re{c4}"  ] = +2.72844;
                p["sbcu::Im{c4}"  ] = +2.4199;
                p["sbcu::Re{c4'}" ] = +1.42602;
                p["sbcu::Im{c4'}" ] = +2.15745;
                p["sbcu::Re{c5}"  ] = +2.1994;
                p["sbcu::Im{c5}"  ] = -1.4183;
                p["sbcu::Re{c5'}" ] = +1.28771;
                p["sbcu::Im{c5'}" ] = -2.51855;
                p["sbcu::Re{c6}"  ] = -1.148;
                p["sbcu::Im{c6}"  ] = +2.69186;
                p["sbcu::Re{c6'}" ] = -0.857562;
                p["sbcu::Im{c6'}" ] = -1.25387;
                p["sbcu::Re{c7}"  ] = -0.0232947;
                p["sbcu::Im{c7}"  ] = +0.746233;
                p["sbcu::Re{c7'}" ] = +0.925099;
                p["sbcu::Im{c7'}" ] = +2.16794;
                p["sbcu::Re{c8}"  ] = -0.787739;
                p["sbcu::Im{c8}"  ] = +2.30108;
                p["sbcu::Re{c8'}" ] = -2.67008;
                p["sbcu::Im{c8'}" ] = -0.331634;
                p["sbcu::Re{c9}"  ] = -1.60631;
                p["sbcu::Im{c9}"  ] = -1.09823;
                p["sbcu::Re{c9'}" ] = +0.601768;
                p["sbcu::Im{c9'}" ] = -0.224144;
                p["sbcu::Re{c10}" ] = +0.25629;
                p["sbcu::Im{c10}" ] = -2.96255;
                p["sbcu::Re{c10'}"] = +2.03425;
                p["sbcu::Im{c10'}"] = +1.24073;

                Options oo
                {
                    { "accuracy",     "LO+NLO" },
                    { "q",            "d"      },
                    { "model",        "WET"    }
                };
                BqToDstarqVector d(p, oo);

                {
                    const double eps = 1.0e-5;
                    TEST_CHECK_RELATIVE_ERROR(d.re_a_1(), 15.054768650304595, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.im_a_1(), -12.30566901066119, eps);

                    TEST_CHECK_RELATIVE_ERROR(d.decay_width(),     7.680144381680253e-14, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(), 0.17723984545666605,   eps);
                }
            }
        }
} bq_to_dstarq_vec_test;
