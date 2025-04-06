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
#include <eos/utils/destringify.hh>
#include <eos/utils/quantum-numbers.hh>
#include <eos/utils/stringify.hh>

using namespace test;
using namespace eos;

class LeptonFlavorTest :
    public TestCase
{
    public:
        LeptonFlavorTest() :
            TestCase("lepton_flavor_test")
        {
        }

        virtual void run() const
        {
            TEST_CHECK_EQUAL_STR("e",   stringify(LeptonFlavor::electron));
            TEST_CHECK_EQUAL_STR("mu",  stringify(LeptonFlavor::muon));
            TEST_CHECK_EQUAL_STR("tau", stringify(LeptonFlavor::tauon));
        }
} lepton_flavor_test;

class QuarkFlavorTest :
    public TestCase
{
    public:
        QuarkFlavorTest() :
            TestCase("quark_flavor_test")
        {
        }

        virtual void run() const
        {
            TEST_CHECK_EQUAL_STR("u", stringify(QuarkFlavor::up));
            TEST_CHECK_EQUAL_STR("d", stringify(QuarkFlavor::down));
            TEST_CHECK_EQUAL_STR("s", stringify(QuarkFlavor::strange));
            TEST_CHECK_EQUAL_STR("c", stringify(QuarkFlavor::charm));
            TEST_CHECK_EQUAL_STR("b", stringify(QuarkFlavor::bottom));
            TEST_CHECK_EQUAL_STR("t", stringify(QuarkFlavor::top));
        }
} quark_flavor_test;

class IsospinTest :
    public TestCase
{
    public:
        IsospinTest() :
            TestCase("isospin_test")
        {
        }

        virtual void run() const
        {
            TEST_CHECK_EQUAL_STR("",      stringify(Isospin::none));
            TEST_CHECK_EQUAL_STR("0",     stringify(Isospin::zero));
            TEST_CHECK_EQUAL_STR("1",     stringify(Isospin::one));
            TEST_CHECK_EQUAL_STR("1/2",   stringify(Isospin::onehalf));
            TEST_CHECK_EQUAL_STR("2",     stringify(Isospin::two));
            TEST_CHECK_EQUAL_STR("3/2",   stringify(Isospin::threehalves));

            TEST_CHECK_EQUAL_STR("0|1",   stringify(Isospin::zero | Isospin::one));
            TEST_CHECK_EQUAL_STR("0|3/2", stringify(Isospin::zero | Isospin::threehalves));
            TEST_CHECK_EQUAL_STR("1|2",   stringify(Isospin::one | Isospin::two));

            TEST_CHECK_EQUAL(destringify<Isospin>("0|1"),   Isospin::zero | Isospin::one);
            TEST_CHECK_EQUAL(destringify<Isospin>("0|3/2"), Isospin::zero | Isospin::threehalves);
            TEST_CHECK_EQUAL(destringify<Isospin>("1|2"),   Isospin::one | Isospin::two);
        }
} isospin_test;

class IsospinRepresentationTest :
    public TestCase
{
    public:
        IsospinRepresentationTest() :
            TestCase("isospin_representation_test")
        {
        }

        virtual void run() const
        {
            TEST_CHECK_EQUAL_STR("0",   stringify(IsospinRepresentation::zero));
            TEST_CHECK_EQUAL_STR("1",   stringify(IsospinRepresentation::one));
            TEST_CHECK_EQUAL_STR("2",   stringify(IsospinRepresentation::two));
            TEST_CHECK_EQUAL_STR("1/2", stringify(IsospinRepresentation::onehalf));
            TEST_CHECK_EQUAL_STR("3/2", stringify(IsospinRepresentation::threehalves));
        }
} isospin_representation_test;

class LightMesonTest :
    public TestCase
{
    public:
        LightMesonTest() :
            TestCase("light_meson_test")
        {
        }

        virtual void run() const
        {
            TEST_CHECK_EQUAL_STR("pi^0",        stringify(LightMeson::pi0));
            TEST_CHECK_EQUAL_STR("pi^+",        stringify(LightMeson::piplus));
            TEST_CHECK_EQUAL_STR("pi^-",        stringify(LightMeson::piminus));
            TEST_CHECK_EQUAL_STR("K_d",         stringify(LightMeson::K0));
            TEST_CHECK_EQUAL_STR("Kbar_d",      stringify(LightMeson::K0bar));
            TEST_CHECK_EQUAL_STR("K_S",         stringify(LightMeson::KS));
            TEST_CHECK_EQUAL_STR("K_u",         stringify(LightMeson::Kplus));
            TEST_CHECK_EQUAL_STR("Kbar_u",      stringify(LightMeson::Kminus));
            TEST_CHECK_EQUAL_STR("eta",         stringify(LightMeson::eta));
            TEST_CHECK_EQUAL_STR("eta_prime",   stringify(LightMeson::etap));

        }
} light_meson_test;
