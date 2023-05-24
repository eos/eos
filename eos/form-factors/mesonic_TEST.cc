/*
 * Copyright (c) 2023 Danny van Dyk
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
#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/mesonic.hh>

#include <eos/models/model.hh>
#include <eos/maths/power-of.hh>

using namespace test;
using namespace eos;

class PToPFormFactorsTest :
    public TestCase
{
    public:
        PToPFormFactorsTest() :
            TestCase("p_to_p_form_factor_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToP>::create("Foo->Bar::BSZ2015", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToP>::create("B->pi::FooBar",     parameter, options));
            }
        }
} p_to_p_form_factor_test;

class PToPPFormFactorsTest :
    public TestCase
{
    public:
        PToPPFormFactorsTest() :
            TestCase("p_to_pp_form_factor_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToPP>::create("Foo->BarBaz::FvDV2018", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToPP>::create("B->pipi::BazBar",       parameter, options));
            }
        }
} p_to_pp_form_factor_test;

class PToVFormFactorsTest :
    public TestCase
{
    public:
        PToVFormFactorsTest() :
            TestCase("p_to_v_form_factor_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToV>::create("Foo->Baz::BSZ2015", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToV>::create("B->rho::FooBaz",    parameter, options));
            }
        }
} p_to_v_form_factor_test;

class PToGammaFormFactorsTest :
    public TestCase
{
    public:
        PToGammaFormFactorsTest() :
            TestCase("p_to_gamma_form_factor_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToGamma>::create("Foo->gluon::FLvD2022QCDF", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToGamma>::create("B->gamma::FooBaz",         parameter, options));
            }
        }
} p_to_gamma_form_factor_test;

class PToGammaOffShellFormFactorsTest :
    public TestCase
{
    public:
        PToGammaOffShellFormFactorsTest() :
            TestCase("p_to_gamma_off_shell_form_factor_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToGammaOffShell>::create("Foo->gluon^*::KKvDZ2022", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToGammaOffShell>::create("B->gamma^*::FooBaz",      parameter, options));
            }
        }
} p_to_gamma_off_shell_form_factor_test;

class VToPFormFactorsTest :
    public TestCase
{
    public:
        VToPFormFactorsTest() :
            TestCase("v_to_p_form_factor_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<VToP>::create("Foo->Baz::BGJvD2019", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<VToP>::create("B^*->D::FooBaz",      parameter, options));
            }
        }
} v_to_p_form_factor_test;

class VToVFormFactorsTest :
    public TestCase
{
    public:
        VToVFormFactorsTest() :
            TestCase("v_to_v_form_factor_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<VToV>::create("Foo->Baz::BGJvD2019", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<VToV>::create("B^*->D^*::FooBaz",    parameter, options));
            }
        }
} v_to_v_form_factor_test;
