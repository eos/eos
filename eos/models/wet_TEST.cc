/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2023 Danny van Dyk
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
#include <eos/maths/complex.hh>
#include <eos/models/model.hh>
#include <eos/models/wet.hh>
#include <eos/models/standard-model.hh>


#include <algorithm>
#include <cmath>
#include <list>

using namespace test;
using namespace eos;

class MakeTest :
    public TestCase
{
    public:
        MakeTest() :
            TestCase("wcm_make_test")
        {
        }

        virtual void run() const
        {
            std::list<std::string> models = {"WET", "WET-SMEFT"};
            for (const auto & name : models)
            {
                try
                {
                    std::shared_ptr<Model> m = Model::make(name, Parameters::Defaults(), Options());
                }
                catch (NoSuchModelError &)
                {
                    TEST_CHECK_FAILED("Model::make does not know the model '" + name + "'");
                }
                catch (...)
                {
                    throw;
                    TEST_CHECK_FAILED("Unknown Exception while making '" + name + "'");
                }
            }
        }
} wsm_make_test;

class WilsonCoefficientsBToSTest :
    public TestCase
{
    public:
        WilsonCoefficientsBToSTest() :
            TestCase("wilson_coefficients_b_to_s_test")
        {
        }

        virtual void run() const
        {
            /* Test passing of SM parameters via cartesian parametrisations */
            {
                static const double eps = 1e-8;
                static const double mu = 4.2; // approximate m_b(m_b) MSbar mass

                Parameters p = Parameters::Defaults();
                p["QCD::alpha_s(MZ)"] = 0.117620;
                p["QCD::mu_t"] = 170.0;
                p["QCD::mu_b"] = 4.2;
                p["QCD::mu_c"] = 1.2;
                p["mass::W"] = 80.398;
                p["mass::Z"] = 91.1876;
                p["mass::t(pole)"] = 173.3;
                p["sb::mu"] = mu;

                Options o;
                o.declare("scan-mode", "cartesian");

                WilsonScanModel model(p, o);

                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(mu, "mu", false);
                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, wc._alpha_s, eps);
                TEST_CHECK_NEARLY_EQUAL(-0.29063621, real(wc.c1()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+1.01029623, real(wc.c2()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.00616220, real(wc.c3()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.08730376, real(wc.c4()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.00042854, real(wc.c5()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.00115807, real(wc.c6()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.33726473, real(wc.c7()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.18288898, real(wc.c8()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+4.27342842, real(wc.c9()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-4.16611761, real(wc.c10()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c1()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c2()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c3()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c4()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c5()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c6()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c7()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c8()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c9()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c10()),  eps);
            }

            /* Test passing of non-SM parameters via cartesian parametrisations */
            {
                static const double eps = 1e-8;
                static const double mu = 4.2; // approximate m_b(m_b) MSbar mass

                Parameters p = Parameters::Defaults();
                p["QCD::alpha_s(MZ)"] = 0.117620;
                p["QCD::mu_t"] = 170.0;
                p["QCD::mu_b"] = 4.2;
                p["QCD::mu_c"] = 1.2;
                p["mass::W"] = 80.398;
                p["mass::Z"] = 91.1876;
                p["mass::t(pole)"] = 173.3;
                p["sb::mu"] = mu;
                p["b->s::Re{c7'}"] = 0.008;
                p["b->s::Im{c7'}"] = M_PI;
                p["b->s::c8'"] = 0.012;
                p["b->see::Re{c9}"] = 3.27;
                p["b->see::Re{c9'}"] = 0.007;
                p["b->see::Im{c9'}"] = 0.01;
                p["b->see::Re{c10'}"] = 0.006;
                p["b->see::Im{c10'}"] = -M_PI+0.01;
                p["b->smumu::Re{c9'}"] = 0.006;
                p["b->smumu::Im{c9'}"] = 0.0;
                p["b->smumu::Re{c10'}"] = 0.005;
                p["b->smumu::Im{c10'}"] = -M_PI;

                Options o;
                o.declare("scan-mode", "cartesian");

                WilsonScanModel model(p, o);

                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(mu, "mu", false);
                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, wc._alpha_s, eps);
                TEST_CHECK_NEARLY_EQUAL(-0.29063621, real(wc.c1()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+1.01029623, real(wc.c2()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.00616220, real(wc.c3()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.08730376, real(wc.c4()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.00042854, real(wc.c5()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.00115807, real(wc.c6()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.33726473, real(wc.c7()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.18288898, real(wc.c8()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+4.27342842, real(wc.c9()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-4.16611761, real(wc.c10()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c1()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c2()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c3()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c4()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c5()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c6()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c7()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c8()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c9()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c10()),  eps);

                TEST_CHECK_NEARLY_EQUAL(+0.008,      real(wc.c7prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.012,      real(wc.c8prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.006,      real(wc.c9prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.005,      real(wc.c10prime()), eps);
                TEST_CHECK_NEARLY_EQUAL(+M_PI,       imag(wc.c7prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c9prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-M_PI,       imag(wc.c10prime()), eps);

                wc = model.wilson_coefficients_b_to_s(mu, "e", false);
                TEST_CHECK_NEARLY_EQUAL(+3.27,       real(wc.c9()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.007,      real(wc.c9prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.006,      real(wc.c10prime()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.01,       imag(wc.c9prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-M_PI+0.01,  imag(wc.c10prime()), eps);
            }
        }
} wilson_coefficients_b_to_s_test;

class WilsonCoefficientsSBSBTest :
    public TestCase
{
    public:
        WilsonCoefficientsSBSBTest() :
            TestCase("wilson_coefficients_sbsb_test")
        {
        }

        virtual void run() const
        {
            /* Test passing of WC via cartesian parametrisations */
            {
                static const double eps = 1e-6;

                Parameters p = Parameters::Defaults();
                p["sbsb::Re{c1}" ] =  0.123456;
                p["sbsb::Im{c1}" ] = -0.234567;
                p["sbsb::Re{c1'}"] = -0.345678;
                p["sbsb::Im{c1'}"] =  0.456789;
                p["sbsb::Re{c2}" ] =  0.567890;
                p["sbsb::Im{c2}" ] = -0.678901;
                p["sbsb::Re{c2'}"] = -0.789012;
                p["sbsb::Im{c2'}"] =  0.890123;
                p["sbsb::Re{c3}" ] =  0.901234;
                p["sbsb::Im{c3}" ] = -0.012345;
                p["sbsb::Re{c3'}"] = -0.123456;
                p["sbsb::Im{c3'}"] =  0.234567;
                p["sbsb::Re{c4}" ] =  0.345678;
                p["sbsb::Im{c4}" ] = -0.456789;
                p["sbsb::Re{c5}" ] = -0.567890;
                p["sbsb::Im{c5}" ] =  0.678901;
                p["sbsb::mu"]      = 4.2;

                Options o{};

                WilsonScanModel model(p, o);

                const auto wc = model.wet_sbsb();
                TEST_CHECK_NEARLY_EQUAL( 0.123456, real(wc.c1()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.234567, imag(wc.c1()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.345678, real(wc.c1p()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.456789, imag(wc.c1p()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.567890, real(wc.c2()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.678901, imag(wc.c2()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.789012, real(wc.c2p()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.890123, imag(wc.c2p()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.901234, real(wc.c3()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.012345, imag(wc.c3()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.123456, real(wc.c3p()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.234567, imag(wc.c3p()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.345678, real(wc.c4()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.456789, imag(wc.c4()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.567890, real(wc.c5()),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.678901, imag(wc.c5()),  eps);
            }
        }
} wilson_coefficients_sbsb_test;

class WilsonCoefficientsDBCUTest :
    public TestCase
{
    public:
        WilsonCoefficientsDBCUTest() :
            TestCase("wilson_coefficients_dbcu_test")
        {
        }

        virtual void run() const
        {
            /* Test default values against the SM */
            {
                static const double eps = 1e-6;

                Parameters p = Parameters::Defaults();
                Options o{};
                StandardModel sm(p);
                WilsonScanModel wet(p, o);

                const auto wc_sm  = sm.wet_dbcu(false);
                const auto wc_wet = wet.wet_dbcu(false);

                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c1()),   real(wc_sm.c1()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c2()),   real(wc_sm.c2()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c3()),   real(wc_sm.c3()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c4()),   real(wc_sm.c4()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c5()),   real(wc_sm.c5()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c6()),   real(wc_sm.c6()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c7()),   real(wc_sm.c7()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c8()),   real(wc_sm.c8()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c9()),   real(wc_sm.c9()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c10()),  real(wc_sm.c10()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c1p()),  real(wc_sm.c1p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c2p()),  real(wc_sm.c2p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c3p()),  real(wc_sm.c3p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c4p()),  real(wc_sm.c4p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c5p()),  real(wc_sm.c5p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c6p()),  real(wc_sm.c6p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c7p()),  real(wc_sm.c7p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c8p()),  real(wc_sm.c8p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c9p()),  real(wc_sm.c9p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c10p()), real(wc_sm.c10p()), eps);
            }

            /* Test passing of WC via cartesian parametrisations */
            {
                static const double eps = 1e-6;

                Parameters p = Parameters::Defaults();
                p["dbcu::Re{c1}"  ] =  0.123456;
                p["dbcu::Im{c1}"  ] = -0.234567;
                p["dbcu::Re{c1'}" ] = -0.345678;
                p["dbcu::Im{c1'}" ] =  0.456789;
                p["dbcu::Re{c2}"  ] =  0.567890;
                p["dbcu::Im{c2}"  ] = -0.678901;
                p["dbcu::Re{c2'}" ] = -0.789012;
                p["dbcu::Im{c2'}" ] =  0.890123;
                p["dbcu::Re{c3}"  ] =  0.901234;
                p["dbcu::Im{c3}"  ] = -0.012345;
                p["dbcu::Re{c3'}" ] = -0.123456;
                p["dbcu::Im{c3'}" ] =  0.234567;
                p["dbcu::Re{c4}"  ] =  0.345678;
                p["dbcu::Im{c4}"  ] = -0.456789;
                p["dbcu::Re{c4'}" ] = -0.567890;
                p["dbcu::Im{c4'}" ] =  0.678901;
                p["dbcu::Re{c5}"  ] =  0.789012;
                p["dbcu::Im{c5}"  ] = -0.890123;
                p["dbcu::Re{c5'}" ] = -0.901234;
                p["dbcu::Im{c5'}" ] =  0.012345;
                p["dbcu::Re{c6}"  ] =  0.123456;
                p["dbcu::Im{c6}"  ] = -0.234567;
                p["dbcu::Re{c6'}" ] = -0.345678;
                p["dbcu::Im{c6'}" ] =  0.456789;
                p["dbcu::Re{c7}"  ] =  0.901234;
                p["dbcu::Im{c7}"  ] = -0.012345;
                p["dbcu::Re{c7'}" ] = -0.123456;
                p["dbcu::Im{c7'}" ] =  0.234567;
                p["dbcu::Re{c8}"  ] =  0.345678;
                p["dbcu::Im{c8}"  ] = -0.456789;
                p["dbcu::Re{c8'}" ] = -0.567890;
                p["dbcu::Im{c8'}" ] =  0.678901;
                p["dbcu::Re{c9}"  ] =  0.901234;
                p["dbcu::Im{c9}"  ] = -0.012345;
                p["dbcu::Re{c9'}" ] = -0.123456;
                p["dbcu::Im{c9'}" ] =  0.234567;
                p["dbcu::Re{c10}" ] =  0.345678;
                p["dbcu::Im{c10}" ] = -0.456789;
                p["dbcu::Re{c10'}"] = -0.567890;
                p["dbcu::Im{c10'}"] =  0.678901;
                p["dbcu::mu"]      = 4.2;

                Options o{};

                WilsonScanModel model(p, o);

                const auto wc = model.wet_dbcu(false);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1()),    0.123456, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1()),   -0.234567, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1p()),  -0.345678, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1p()),   0.456789, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2()),    0.56789, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2()),   -0.678901, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2p()),  -0.789012, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2p()),   0.890123, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3()),    0.901234, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3()),   -0.012345, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3p()),  -0.123456, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3p()),   0.234567, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c4()),    0.345678, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4()),   -0.456789, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c4p()),  -0.56789, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4p()),   0.678901, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c5()),    0.789012, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5()),   -0.890123, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c5p()),  -0.901234, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5p()),   0.012345, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c6()),    0.123456, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6()),   -0.234567, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c6p()),  -0.345678, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6p()),   0.456789, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c7()),    0.901234, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7()),   -0.012345, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c7p()),  -0.123456, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7p()),   0.234567, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c8()),    0.345678, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8()),   -0.456789, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c8p()),  -0.56789, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8p()),   0.678901, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c9()),    0.901234, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9()),   -0.012345, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c9p()),  -0.123456, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9p()),   0.234567, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c10()),   0.345678, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10()),  -0.456789, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c10p()), -0.56789, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10p()),  0.678901, eps);
            }
        }
} wilson_coefficients_dbcu_test;

class WilsonCoefficientsSBCUTest :
    public TestCase
{
    public:
        WilsonCoefficientsSBCUTest() :
            TestCase("wilson_coefficients_sbcu_test")
        {
        }

        virtual void run() const
        {
            /* Test default values against the SM */
            {
                static const double eps = 1e-6;

                Parameters p = Parameters::Defaults();
                Options o{};
                StandardModel sm(p);
                WilsonScanModel wet(p, o);

                const auto wc_sm  = sm.wet_sbcu(false);
                const auto wc_wet = wet.wet_sbcu(false);

                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c1()),   real(wc_sm.c1()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c2()),   real(wc_sm.c2()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c3()),   real(wc_sm.c3()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c4()),   real(wc_sm.c4()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c5()),   real(wc_sm.c5()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c6()),   real(wc_sm.c6()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c7()),   real(wc_sm.c7()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c8()),   real(wc_sm.c8()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c9()),   real(wc_sm.c9()),   eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c10()),  real(wc_sm.c10()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c1p()),  real(wc_sm.c1p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c2p()),  real(wc_sm.c2p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c3p()),  real(wc_sm.c3p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c4p()),  real(wc_sm.c4p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c5p()),  real(wc_sm.c5p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c6p()),  real(wc_sm.c6p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c7p()),  real(wc_sm.c7p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c8p()),  real(wc_sm.c8p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c9p()),  real(wc_sm.c9p()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc_wet.c10p()), real(wc_sm.c10p()), eps);
            }

            /* Test passing of WC via cartesian parametrisations */
            {
                static const double eps = 1e-6;

                Parameters p = Parameters::Defaults();
                p["sbcu::Re{c1}"  ] =  0.123456;
                p["sbcu::Im{c1}"  ] = -0.234567;
                p["sbcu::Re{c1'}" ] = -0.345678;
                p["sbcu::Im{c1'}" ] =  0.456789;
                p["sbcu::Re{c2}"  ] =  0.567890;
                p["sbcu::Im{c2}"  ] = -0.678901;
                p["sbcu::Re{c2'}" ] = -0.789012;
                p["sbcu::Im{c2'}" ] =  0.890123;
                p["sbcu::Re{c3}"  ] =  0.901234;
                p["sbcu::Im{c3}"  ] = -0.012345;
                p["sbcu::Re{c3'}" ] = -0.123456;
                p["sbcu::Im{c3'}" ] =  0.234567;
                p["sbcu::Re{c4}"  ] =  0.345678;
                p["sbcu::Im{c4}"  ] = -0.456789;
                p["sbcu::Re{c4'}" ] = -0.567890;
                p["sbcu::Im{c4'}" ] =  0.678901;
                p["sbcu::Re{c5}"  ] =  0.789012;
                p["sbcu::Im{c5}"  ] = -0.890123;
                p["sbcu::Re{c5'}" ] = -0.901234;
                p["sbcu::Im{c5'}" ] =  0.012345;
                p["sbcu::Re{c6}"  ] =  0.123456;
                p["sbcu::Im{c6}"  ] = -0.234567;
                p["sbcu::Re{c6'}" ] = -0.345678;
                p["sbcu::Im{c6'}" ] =  0.456789;
                p["sbcu::Re{c7}"  ] =  0.901234;
                p["sbcu::Im{c7}"  ] = -0.012345;
                p["sbcu::Re{c7'}" ] = -0.123456;
                p["sbcu::Im{c7'}" ] =  0.234567;
                p["sbcu::Re{c8}"  ] =  0.345678;
                p["sbcu::Im{c8}"  ] = -0.456789;
                p["sbcu::Re{c8'}" ] = -0.567890;
                p["sbcu::Im{c8'}" ] =  0.678901;
                p["sbcu::Re{c9}"  ] =  0.901234;
                p["sbcu::Im{c9}"  ] = -0.012345;
                p["sbcu::Re{c9'}" ] = -0.123456;
                p["sbcu::Im{c9'}" ] =  0.234567;
                p["sbcu::Re{c10}" ] =  0.345678;
                p["sbcu::Im{c10}" ] = -0.456789;
                p["sbcu::Re{c10'}"] = -0.567890;
                p["sbcu::Im{c10'}"] =  0.678901;
                p["sbcu::mu"]      = 4.2;

                Options o{};

                WilsonScanModel model(p, o);

                const auto wc = model.wet_sbcu(false);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1()),    0.123456, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1()),   -0.234567, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1p()),  -0.345678, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1p()),   0.456789, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2()),    0.56789, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2()),   -0.678901, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2p()),  -0.789012, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2p()),   0.890123, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3()),    0.901234, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3()),   -0.012345, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3p()),  -0.123456, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3p()),   0.234567, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c4()),    0.345678, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4()),   -0.456789, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c4p()),  -0.56789, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4p()),   0.678901, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c5()),    0.789012, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5()),   -0.890123, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c5p()),  -0.901234, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5p()),   0.012345, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c6()),    0.123456, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6()),   -0.234567, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c6p()),  -0.345678, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6p()),   0.456789, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c7()),    0.901234, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7()),   -0.012345, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c7p()),  -0.123456, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7p()),   0.234567, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c8()),    0.345678, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8()),   -0.456789, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c8p()),  -0.56789, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8p()),   0.678901, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c9()),    0.901234, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9()),   -0.012345, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c9p()),  -0.123456, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9p()),   0.234567, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c10()),   0.345678, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10()),  -0.456789, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c10p()), -0.56789, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10p()),  0.678901, eps);
            }
        }
} wilson_coefficients_sbcu_test;

class WilsonCoefficientsSBNuNuTest :
    public TestCase
{
    public:
        WilsonCoefficientsSBNuNuTest() :
            TestCase("wilson_coefficients_sbnunu_test")
        {
        }

        virtual void run() const
        {
            /* Test comparing WC of WET and SM */
            {
                static const double eps = 1e-8;

                Parameters p = Parameters::Defaults();

                Options o{};

                StandardModel sm(p);

                const auto sm_wc = sm.wet_sbnunu(false);

                WilsonScanModel wsm(p, o);

                const auto wsm_wc = wsm.wet_sbnunu(false);

                TEST_CHECK_NEARLY_EQUAL(real(sm_wc.cVL()), real(wsm_wc.cVL()),  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(sm_wc.cVL()), imag(wsm_wc.cVL()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(sm_wc.cVR()), real(wsm_wc.cVR()),  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(sm_wc.cVR()), imag(wsm_wc.cVR()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(sm_wc.cSL()), real(wsm_wc.cSL()),  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(sm_wc.cSL()), imag(wsm_wc.cSL()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(sm_wc.cSR()), real(wsm_wc.cSR()),  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(sm_wc.cSR()), imag(wsm_wc.cSR()),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(sm_wc.cTL()), real(wsm_wc.cTL()),  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(sm_wc.cTL()), imag(wsm_wc.cTL()),  eps);
            }
        }
} wilson_coefficients_sbnunu_test;

class ConstrainedWilsonScanModelTest:
    public TestCase
{
    public:
        ConstrainedWilsonScanModelTest() :
            TestCase("constrained_wilson_scan_model_test")
        {
        }

        virtual void run() const
        {
            static const double mu  = 4.2; // approximate value of the b quark mass in the MSbar scheme
            static const double eps = 1e-15;

            /* Vary parameters that should be ignored */
            {
                Parameters p = Parameters::Defaults();
                Options o;
                ConstrainedWilsonScanModel model(p, o);

                p["b->s::Re{c7}"] = 1.008;
                p["b->smumu::Re{cS}"] = 42;
                p["b->smumu::Re{cP}"] = 100;
                p["b->smumu::Im{cS'}"] = -12;
                p["b->smumu::Im{cP'}"] = -135;
                p["b->smumu::Re{cT}"] = 2.0;
                p["b->smumu::Re{cT5}"] = -43.0;

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(mu, "mu", false);

                TEST_CHECK_RELATIVE_ERROR(std::real(wc.c7()),  1.008, eps);

                /* C_P should be ignored, and always equal -C_S */
                TEST_CHECK_RELATIVE_ERROR(std::real(wc.cS()),  42,    eps);
                TEST_CHECK_RELATIVE_ERROR(std::real(wc.cP()), -42,    eps);

                TEST_CHECK_RELATIVE_ERROR(imag(wc.cSprime()), -12,    eps);
                TEST_CHECK_RELATIVE_ERROR(imag(wc.cPprime()), -12,    eps);

                /* C_T and C_T5 vanish */
                TEST_CHECK_NEARLY_EQUAL(std::real(wc.cT()),     0.0,  eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(wc.cT()),     0.0,  eps);

                TEST_CHECK_NEARLY_EQUAL(std::real(wc.cT5()),    0.0,  eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(wc.cT5()),    0.0,  eps);

                /* Used parameters registered */
                TEST_CHECK(std::find(std::begin(model), std::end(model), p["b->smumu::Re{cS}"].id()) != std::end(model));
                TEST_CHECK(std::find(std::begin(model), std::end(model), p["b->smumu::Im{cS}"].id()) != std::end(model));

                std::list<Parameter::Id> unused_ids = {
                        p["b->smumu::Re{cP}"].id(),
                        p["b->smumu::Im{cP}"].id(),
                        p["b->smumu::Re{cP'}"].id(),
                        p["b->smumu::Im{cP'}"].id(),
                        p["b->smumu::Re{cT}"].id(),
                        p["b->smumu::Im{cT}"].id(),
                        p["b->smumu::Re{cT5}"].id(),
                        p["b->smumu::Im{cT5}"].id(),
                };
                for (auto & id : model)
                {
                    TEST_CHECK(std::find(unused_ids.begin(), unused_ids.end(), id) == unused_ids.end());
                }
            }

            /* cartesian parametrisation */
            {
                Parameters p = Parameters::Defaults();
                Options o;
                ConstrainedWilsonScanModel model(p, o);

                p["b->s::Re{c7}"] = 1.008;
                p["b->smumu::Re{cS}"] = 42;
                p["b->smumu::Im{cS}"] = 0.5;
                p["b->smumu::Re{cS'}"] = 3.2;
                p["b->smumu::Im{cS'}"] = 1.2;
                p["b->smumu::Re{cP}"] = 100;
                p["b->smumu::Im{cP'}"] = 35;
                p["b->smumu::Re{cT}"] = 2.0;
                p["b->smumu::Im{cT}"] = 9.0;
                p["b->smumu::Re{cT5}"] = -43.0;
                p["b->smumu::Im{cT5}"] = M_PI;

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(mu, "mu", false);

                TEST_CHECK_RELATIVE_ERROR(real(wc.c7()),      1.008, eps);

                /* C_P should be ignored, and always equal -C_S */
                TEST_CHECK_RELATIVE_ERROR(real(wc.cS()),     42.0,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(wc.cS()),      0.5,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.cP()),    -42.0,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(wc.cP()),     -0.5,   eps);

                TEST_CHECK_RELATIVE_ERROR(real(wc.cSprime()), 3.2,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(wc.cSprime()), 1.2,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.cPprime()), 3.2,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(wc.cPprime()), 1.2,   eps);

                /* C_T and C_T5 vanish */
                TEST_CHECK_NEARLY_EQUAL(std::real(wc.cT()),   0.0,   eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(wc.cT()),   0.0,   eps);

                TEST_CHECK_NEARLY_EQUAL(std::real(wc.cT5()),  0.0,   eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(wc.cT5()),  0.0,   eps);

                /* Used parameters registered */
                TEST_CHECK(std::find(std::begin(model), std::end(model), p["b->smumu::Re{cS}"].id()) != std::end(model));
                TEST_CHECK(std::find(std::begin(model), std::end(model), p["b->smumu::Re{cS}"].id()) != std::end(model));

                std::list<Parameter::Id> unused_ids = {
                    p["b->smumu::Re{cP}"].id(),
                    p["b->smumu::Im{cP}"].id(),
                    p["b->smumu::Re{cP'}"].id(),
                    p["b->smumu::Im{cP'}"].id(),
                    p["b->smumu::Re{cT}"].id(),
                    p["b->smumu::Im{cT}"].id(),
                    p["b->smumu::Re{cT5}"].id(),
                    p["b->smumu::Im{cT5}"].id(),
                };
                for (auto & id : model)
                {
                    TEST_CHECK(std::find(unused_ids.begin(), unused_ids.end(), id) == unused_ids.end());
                }
            }

            /* most parameters identical to the usual WilsonScanModel */
            {
                Parameters p = Parameters::Defaults();
                Options o;
                o.declare("scan-mode", "cartesian");

                p["b->s::Re{c7}"] = 1.008;
                p["b->smumu::Re{cS}"] = 42;   p["b->smumu::Re{cP}"] = -1.0 * p["b->smumu::Re{cS}"]();
                p["b->smumu::Im{cS'}"] = -12; p["b->smumu::Im{cP'}"] = p["b->smumu::Im{cS'}"]();
                p["b->smumu::Re{cT}"] = 0.0;  p["b->smumu::Im{cT}"] = 0.0;
                p["b->smumu::Re{cT5}"] = 0.0; p["b->smumu::Im{cT5}"] = 0.0;

                ConstrainedWilsonScanModel constrained_model(p, o);
                WilsonScanModel unconstrained_model(p, o);

                WilsonCoefficients<BToS> constrained_wc = constrained_model.wilson_coefficients_b_to_s(mu, "mu", false);
                WilsonCoefficients<BToS> unconstrained_wc = constrained_model.wilson_coefficients_b_to_s(mu, "mu", false);

                auto ux = unconstrained_wc._sm_like_coefficients.begin();
                for (auto & x : constrained_wc._sm_like_coefficients)
                {
                    TEST_CHECK_EQUAL(x, *ux);
                    ++ux;
                }

                ux = unconstrained_wc._primed_coefficients.begin();
                for (auto & x : constrained_wc._primed_coefficients)
                {
                    TEST_CHECK_EQUAL(x, *ux);
                    ++ux;
                }

                ux = unconstrained_wc._scalar_tensor_coefficients.begin();
                for (auto & x : constrained_wc._scalar_tensor_coefficients)
                {
                    TEST_CHECK_EQUAL(x, *ux);
                    ++ux;
                }
            }
        }
} constrained_wilson_scan_model_test;
