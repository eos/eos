/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021-2025 Danny van Dyk
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
#include <eos/meson-mixing/bq-mixing.hh>

using namespace test;
using namespace eos;

class BsMixingTest :
    public TestCase
{
    public:
        BsMixingTest() :
            TestCase("b_s_mixing_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();
                p["CKM::lambda"]           =  0.22535;
                p["CKM::A"]                =  0.827;
                p["CKM::rhobar"]           =  0.132;
                p["CKM::etabar"]           =  0.350;
                // Using [DDHLMSW:2019A] inputs for the reduced matrix elements.
                p["B_s<->Bbar_s::R^1"]     =  0.54200;
                p["B_s<->Bbar_s::R^2"]     = -0.54500;
                p["B_s<->Bbar_s::R^3"]     =  0.10900;
                p["B_s<->Bbar_s::R^4"]     =  0.91250;
                p["B_s<->Bbar_s::R^5"]     =  0.48625;

                Options oo
                {
                    { "model"_ok,        "SM"         },
                    { "q"_ok,            "s"          },
                };

                BMixing process(p, oo);

                const double eps = 1.0e-5;

                TEST_CHECK_RELATIVE_ERROR(process.delta_m(), 17.26561, eps); // in units of ps^-1
            }
        }
} b_s_mixing_test;
