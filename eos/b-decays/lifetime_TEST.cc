/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Danny van Dyk
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
#include <eos/b-decays/lifetime.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <vector>


using namespace test;
using namespace eos;

class LifetimeTest :
    public TestCase
{
    public:
        LifetimeTest() :
            TestCase("lifetime_test")
        {
        }

        virtual void run() const
        {
            // SM tests
            {
                Parameters p = Parameters::Defaults();
                p["WET::G_Fermi"]        = 1.1664e-05;
                p["CKM::abs(V_ud)"]      = 1.000;
                p["CKM::abs(V_us)"]      = 0.2254;
                p["CKM::abs(V_cb)"]      = 4.2e-2;
                p["mass::B_u"]           = 5.279;
                p["mass::B_d"]           = 5.279;
                p["mass::B_s"]           = 5.366;
                p["decay-constant::B_u"] = 0.1905;
                p["decay-constant::B_d"] = 0.8940;
                p["decay-constant::B_s"] = 0.2307;
                // SM test for the B^-
                {
                    static const double eps = 1.0e-6;

                    Options oo
                    {
                        { "model", "CKM" },
                        { "q",     "u"   }
                    };

                    Lifetime d(p, oo);

                    // compare to known decay width of
                    // (1.638 +/- 0.004 ps)^-1 = (0.6105 +/- 0.0015) ps^-1
                    TEST_CHECK_NEARLY_EQUAL(d.decay_width_dbcu(), 0.0153181,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.decay_width_sbcu(), 0.000778238, eps);
                }
                // SM test for the B_d
                {
                    static const double eps = 1.0e-7;

                    Options oo
                    {
                        { "model", "CKM" },
                        { "q",     "d"   }
                    };

                    Lifetime d(p, oo);

                    // compare to known decay width of
                    // (1.519 +/- 0.004 ps)^-1 = (0.6583 +/- 0.0017) ps^-1
                    TEST_CHECK_NEARLY_EQUAL(d.decay_width_dbcu(), 0.0138916, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.decay_width_sbcu(), 0.0,       eps);
                }
                // SM test for the B_s
                {
                    static const double eps = 1.0e-7;

                    Options oo
                    {
                        { "model", "CKM" },
                        { "q",     "s"   }
                    };

                    Lifetime d(p, oo);

                    // compare to known decay width of
                    // (1.521 +/- 0.005 ps)^-1 = (0.6575 +/- 0.0022) ps^-1
                    TEST_CHECK_NEARLY_EQUAL(d.decay_width_dbcu(), 0.0,          eps);
                    TEST_CHECK_NEARLY_EQUAL(d.decay_width_sbcu(), 0.0000477726, eps);
                }
            }
            // BSM tests
            {
                Parameters p = Parameters::Defaults();
                p["WET::G_Fermi"]        = 1.1664e-05;
                p["CKM::abs(V_ud)"]      = 1.000;
                p["CKM::abs(V_us)"]      = 0.2254;
                p["CKM::abs(V_cb)"]      = 4.2e-2;
                p["mass::B_u"]           = 5.279;
                p["mass::B_d"]           = 5.279;
                p["mass::B_s"]           = 5.366;
                p["decay-constant::B_u"] = 0.1905;
                p["decay-constant::B_d"] = 0.8940;
                p["decay-constant::B_s"] = 0.2307;
                // dbcu WC
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
                // sbcu WC
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
                // BSM test for the B^-
                {
                    static const double eps = 1.0e-5;

                    Options oo
                    {
                        { "model", "WET" },
                        { "q",     "u"   }
                    };

                    Lifetime d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.decay_width_dbcu(), 3683.55940755512,   eps);
                    TEST_CHECK_RELATIVE_ERROR(d.decay_width_sbcu(), 187.1438250703431, eps);
                }
                // BSM test for the B_d
                {
                    static const double eps = 1.0e-6;

                    Options oo
                    {
                        { "model", "WET" },
                        { "q",     "d"   }
                    };

                    Lifetime d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.decay_width_dbcu(), 49570.26112482211, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.decay_width_sbcu(),   0.0,               eps);
                }
                // BSM test for the B_s
                {
                    static const double eps = 1.0e-6;

                    Options oo
                    {
                        { "model", "WET" },
                        { "q",     "s"   }
                    };

                    Lifetime d(p, oo);

                    TEST_CHECK_NEARLY_EQUAL(d.decay_width_dbcu(),   0.0,                eps);
                    TEST_CHECK_RELATIVE_ERROR(d.decay_width_sbcu(), 170.47008899539506, eps);
                }
            }
        }
} lifetime_test;
