/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013 Danny van Dyk
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
#include <eos/utils/model.hh>
#include <eos/utils/wilson_scan_model.hh>

#include <cmath>

using namespace test;
using namespace eos;

Parameters
reference_parameters()
{
    Parameters result = Parameters::Defaults();
    result["QCD::alpha_s(MZ)"] = 0.117620;
    result["QCD::mu_t"] = 170.0;
    result["QCD::mu_b"] = 4.2;
    result["QCD::mu_c"] = 1.2;
    result["mass::W"] = 80.398;
    result["mass::Z"] = 91.1876;
    result["mass::t(pole)"] = 173.3;

    return result;
}

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
            try
            {
                std::shared_ptr<Model> m = Model::make("WilsonScan", reference_parameters(), Options());
            }
            catch (NoSuchModelError &)
            {
                TEST_CHECK_FAILED("Model::make does not know the model 'WilsonScan'");
            }
            catch (...)
            {
                throw;
                TEST_CHECK_FAILED("Unknown Exception while making 'WilsonScan'");
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
            /* Test passing of SM parameters via polar parametrisations */
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
                p["mu"] = mu;

                Options o;
                o.set("scan-mode", "polar");

                WilsonScanModel model(p, o);

                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(false);
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
                p["mu"] = mu;

                Options o;
                o.set("scan-mode", "cartesian");

                WilsonScanModel model(p, o);

                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(false);
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

            /* Test passing of non-SM parameters via polar parametrisations */
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
                p["mu"] = mu;
                p["Abs{c7'}"] = 0.008;
                p["Arg{c7'}"] = M_PI;
                p["c8'"] = 0.012;
                p["Abs{c9'}"] = 0.006;
                p["Arg{c9'}"] = 0.0;
                p["Abs{c10'}"] = 0.005;
                p["Arg{c10'}"] = M_PI;

                Options o;
                o.set("scan-mode", "polar");

                WilsonScanModel model(p, o);

                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(false);
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
                TEST_CHECK_NEARLY_EQUAL(-0.0,        imag(wc.c10()), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.008,      real(wc.c7prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.012,      real(wc.c8prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.006,      real(wc.c9prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.005,      real(wc.c10prime()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c7prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c9prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0,        imag(wc.c10prime()), eps);
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
                p["mu"] = mu;
                p["Re{c7'}"] = 0.008;
                p["Im{c7'}"] = M_PI;
                p["c8'"] = 0.012;
                p["Re{c9'}"] = 0.006;
                p["Im{c9'}"] = 0.0;
                p["Re{c10'}"] = 0.005;
                p["Im{c10'}"] = -M_PI;

                Options o;
                o.set("scan-mode", "cartesian");

                WilsonScanModel model(p, o);

                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(false);
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
            }
        }
} wilson_coefficients_b_to_s_test;
