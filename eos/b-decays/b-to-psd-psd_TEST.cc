/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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
#include <eos/b-decays/b-to-psd-psd.hh>
#include <eos/maths/complex.hh>

using namespace test;
using namespace eos;

class BToPseudoscalarPseudoscalarTest :
    public TestCase
{
    public:
        BToPseudoscalarPseudoscalarTest() :
            TestCase("b_to_psd_psd_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();
            p["nonleptonic::Re{AT3}@SU3F"]  =  0.1;
            p["nonleptonic::Re{CT3}@SU3F"]  = -0.2;
            p["nonleptonic::Re{AT6}@SU3F"]  =  0.3;
            p["nonleptonic::Re{CT6}@SU3F"]  = -0.4;
            p["nonleptonic::Re{AT15}@SU3F"] =  0.5;
            p["nonleptonic::Re{CT15}@SU3F"] = -0.6;
            p["nonleptonic::Re{BT3}@SU3F"]  =  0.7;
            p["nonleptonic::Re{BT6}@SU3F"]  = -0.8;
            p["nonleptonic::Re{BT15}@SU3F"] =  0.9;
            p["nonleptonic::Re{DT3}@SU3F"]  = -1.0;
            p["nonleptonic::Re{AP3}@SU3F"]  =  1.1;
            p["nonleptonic::Re{CP3}@SU3F"]  = -1.2;
            p["nonleptonic::Re{AP6}@SU3F"]  =  1.3;
            p["nonleptonic::Re{CP6}@SU3F"]  = -1.4;
            p["nonleptonic::Re{AP15}@SU3F"] =  1.5;
            p["nonleptonic::Re{CP15}@SU3F"] = -1.6;
            p["nonleptonic::Re{BP3}@SU3F"]  =  1.7;
            p["nonleptonic::Re{BP6}@SU3F"]  = -1.8;
            p["nonleptonic::Re{BP15}@SU3F"] =  1.9;
            p["nonleptonic::Re{DP3}@SU3F"]  = -2.0;
            p["eta::theta_18"]              =  0.0;

            static const double eps = 1.0e-6;

            Options o
            {
                { "representation", "SU3F" },
                { "q", "u" },
                { "P1", "pi^0" },
                { "P2", "pi^+" },
                { "model", "CKM" }
            };

            BToPseudoscalarPseudoscalar d(p, o);

            TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(), 0.004160580429772407, eps);

            Options oo
            {
                { "representation", "SU3F" },
                { "q", "s" },
                { "P1", "eta" },
                { "P2", "Kbar_d" },
                { "model", "CKM" }
            };

            BToPseudoscalarPseudoscalar dd(p, oo);

            TEST_CHECK_RELATIVE_ERROR(dd.branching_ratio(), 5.276489800357915e-05, eps);
        }

} b_to_psd_psd_test;
