/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Stefan Meiser
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
#include <eos/b-decays/bq-to-dq-vec.hh>
#include <eos/maths/complex.hh>

using namespace test;
using namespace eos;

class BqToDqVecTest :
    public TestCase
{
    public:
        BqToDqVecTest() :
            TestCase("bq_to_dq_vec_test")
        {
        }

        virtual void run() const
        {
            // Standard Model
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_s"] = 1.521e-12;
                // Meson masses
                p["mass::B_s"] = 5.36692;
                p["mass::D_s"] = 1.96835;
                p["mass::rho^+"] = 0.77526;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                p["mass::c"] = 1.2;
                // decay constant
                p["decay-constant::rho"] = 0.216;
                // fermi constant
                p["WET::G_Fermi"] = 1.16637e-05;
                // CKM matrix elements
                p["CKM::abs(V_cb)"] = 0.0408;
                p["CKM::arg(V_cb)"] = 0;
                p["CKM::abs(V_ud)"] = 0.97373;
                p["CKM::arg(V_ud)"] = 0;
                // form factor
                p["B_s->D_srho::f_p(Mrho2)"] = 0.685068;
                // WC Values
                p["dbcu::Re{c1}" ] = -0.04235657776117585;
                p["dbcu::Im{c1}" ] = 0;
                p["dbcu::Re{c2}" ] = -0.8948941708221622;
                p["dbcu::Im{c2}" ] = 0;
                p["dbcu::Re{c3}" ] = 0.011381250932999982;
                p["dbcu::Im{c3}" ] = 0;
                p["dbcu::Re{c4}" ] = 0.19426386543613433;
                p["dbcu::Im{c4}" ] = 0;
                Options oo
                {
                    { "accuracy",     "LO+NLO" },
                    { "q",            "s"      },
                    { "model",        "WET"    }
                };
                BqToDqVector d(p, oo);

                {
                    const double eps = 1.0e-05;
                    TEST_CHECK_RELATIVE_ERROR(d.re_a_1(), -1.051961367517211,    eps);
                    TEST_CHECK_RELATIVE_ERROR(d.im_a_1(), -0.016829537324562847, eps);

                    TEST_CHECK_RELATIVE_ERROR(d.decay_width(),     4.70856e-15, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(), 0.0108806  , eps);
                }
            }

            // BSM test case
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_s"] = 1.521e-12;
                // Meson masses
                p["mass::B_s"] = 5.36692;
                p["mass::D_s"] = 1.96835;
                p["mass::rho^+"] = 0.77526;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                p["mass::c"] = 1.2;
                // decay constant
                p["decay-constant::rho"] = 0.216;
                // fermi constant
                p["WET::G_Fermi"] = 1.16637e-05;
                // CKM matrix elements
                p["CKM::abs(V_cb)"] = 0.0408;
                p["CKM::arg(V_cb)"] = 0;
                p["CKM::abs(V_ud)"] = 0.97373;
                p["CKM::arg(V_ud)"] = 0;
                // form factor
                p["B_s->D_srho::f_p(Mrho2)"] = 0.685068;
                // BSM test case: WET parameter point
                p["dbcu::Re{c1}"  ] = -1.72424;
                p["dbcu::Im{c1}"  ] = -1.56379;
                p["dbcu::Re{c1'}" ] = -1.05356;
                p["dbcu::Im{c1'}" ] = -0.791464;
                p["dbcu::Re{c2}"  ] = -2.84324;
                p["dbcu::Im{c2}"  ] = -1.10401;
                p["dbcu::Re{c2'}" ] = +1.10235;
                p["dbcu::Im{c2'}" ] = +2.0774;
                p["dbcu::Re{c3}"  ] = +1.61473;
                p["dbcu::Im{c3}"  ] = +1.23153;
                p["dbcu::Re{c3'}" ] = -2.95587;
                p["dbcu::Im{c3'}" ] = -2.28859;
                p["dbcu::Re{c4}"  ] = +2.72844;
                p["dbcu::Im{c4}"  ] = +2.4199;
                p["dbcu::Re{c4'}" ] = +1.42602;
                p["dbcu::Im{c4'}" ] = +2.15745;
                p["dbcu::Re{c5}"  ] = +2.1994;
                p["dbcu::Im{c5}"  ] = -1.4183;
                p["dbcu::Re{c5'}" ] = +1.28771;
                p["dbcu::Im{c5'}" ] = -2.51855;
                p["dbcu::Re{c6}"  ] = -1.148;
                p["dbcu::Im{c6}"  ] = +2.69186;
                p["dbcu::Re{c6'}" ] = -0.857562;
                p["dbcu::Im{c6'}" ] = -1.25387;
                p["dbcu::Re{c7}"  ] = -0.0232947;
                p["dbcu::Im{c7}"  ] = +0.746233;
                p["dbcu::Re{c7'}" ] = +0.925099;
                p["dbcu::Im{c7'}" ] = +2.16794;
                p["dbcu::Re{c8}"  ] = -0.787739;
                p["dbcu::Im{c8}"  ] = +2.30108;
                p["dbcu::Re{c8'}" ] = -2.67008;
                p["dbcu::Im{c8'}" ] = -0.331634;
                p["dbcu::Re{c9}"  ] = -1.60631;
                p["dbcu::Im{c9}"  ] = -1.09823;
                p["dbcu::Re{c9'}" ] = +0.601768;
                p["dbcu::Im{c9'}" ] = -0.224144;
                p["dbcu::Re{c10}" ] = +0.25629;
                p["dbcu::Im{c10}" ] = -2.96255;
                p["dbcu::Re{c10'}"] = +2.03425;
                p["dbcu::Im{c10'}"] = +1.24073;

                Options oo
                {
                    { "accuracy",     "LO+NLO" },
                    { "q",            "s"      },
                    { "model",        "WET"    }
                };
                BqToDqVector d(p, oo);

                {
                    const double eps = 1.0e-5;
                    TEST_CHECK_RELATIVE_ERROR(d.re_a_1(), -13.343504533329863, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.im_a_1(), -39.050977290302455,  eps);

                    TEST_CHECK_RELATIVE_ERROR(d.decay_width(),     7.244337799221826e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(), 16.740256623422844,    eps);
                }
            }
        }
} bq_to_dq_vec_test;
