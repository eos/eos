/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Danny van Dyk
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
#include <eos/b-decays/properties.hh>

using namespace test;
using namespace eos;

class BMesonPropertiesTest :
    public TestCase
{
    public:
        BMesonPropertiesTest() :
            TestCase("b_meson_properties_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();
            p["mass::b(MSbar)"] = 4.19;
            p["B->B::mu_G^2@1GeV"] = 0.35;

            Options oo;

            BMesonProperties d(p, oo);

            const double eps = 1e-5;
            TEST_CHECK_NEARLY_EQUAL(d.mass_splitting_j1_j0(), 4.57252e-2, eps);
        }
} b_meson_properties_test;
