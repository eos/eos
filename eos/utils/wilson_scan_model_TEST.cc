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
#include <eos/utils/complex.hh>
#include <eos/utils/model.hh>
#include <eos/utils/wilson_scan_model.hh>

#include <algorithm>
#include <cmath>
#include <list>

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
            std::list<std::string> models = {"WilsonScan", "ConstrainedWilsonScan"};
            for (const auto & name : models)
            {
                try
                {
                    std::shared_ptr<Model> m = Model::make(name, reference_parameters(), Options());
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
            static const double eps = 1e-15;

            /* Vary parameters that should be ignored */
            {
                Parameters p = Parameters::Defaults();
                Options o;
                o.set("scan-mode", "cartesian");
                ConstrainedWilsonScanModel model(p, o);

                p["Re{c7}"] = 1.008;
                p["Re{cS}"] = 42;
                p["Re{cP}"] = 100;
                p["Im{cS'}"] = -12;
                p["Im{cP'}"] = -135;
                p["Re{cT}"] = 2.0;
                p["Re{cT5}"] = -43.0;

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(false);

                TEST_CHECK_RELATIVE_ERROR(std::real(wc.c7()), 1.008, eps);

                /* C_P should be ignored, and always equal -C_S */
                TEST_CHECK_RELATIVE_ERROR(std::real(wc.cS()), 42, eps);
                TEST_CHECK_RELATIVE_ERROR(std::real(wc.cP()), -42, eps);

                TEST_CHECK_RELATIVE_ERROR(imag(wc.cSprime()), -12, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(wc.cPprime()), -12, eps);

                /* C_T and C_T5 vanish */
                TEST_CHECK_NEARLY_EQUAL(std::real(wc.cT()), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(wc.cT()), 0.0, eps);

                TEST_CHECK_NEARLY_EQUAL(std::real(wc.cT5()), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(wc.cT5()), 0.0, eps);

                /* Used parameters registered */
                TEST_CHECK(std::find(std::begin(model), std::end(model), p["Re{cS}"].id()) != std::end(model));
                TEST_CHECK(std::find(std::begin(model), std::end(model), p["Im{cS}"].id()) != std::end(model));

                std::list<Parameter::Id> unused_ids = {
                        p["Re{cP}"].id(),
                        p["Im{cP}"].id(),
                        p["Re{cP'}"].id(),
                        p["Im{cP'}"].id(),
                        p["Re{cT}"].id(),
                        p["Im{cT}"].id(),
                        p["Re{cT5}"].id(),
                        p["Im{cT5}"].id(),
                };
                for (auto & id : model)
                {
                    TEST_CHECK(std::find(unused_ids.begin(), unused_ids.end(), id) == unused_ids.end());
                }
            }

            /* polar parametrisation */
            {
                Parameters p = Parameters::Defaults();
                Options o;
                o.set("scan-mode", "polar");
                ConstrainedWilsonScanModel model(p, o);

                p["Abs{c7}"] = 1.008;
                p["Abs{cS}"] = 42;
                p["Arg{cS}"] = 0.5;
                p["Abs{cP}"] = 100;
                p["Abs{cS'}"] = 3.2;
                p["Arg{cS'}"] = 1.2;
                p["Abs{cP'}"] = 35;
                p["Arg{cP'}"] = -0.2;
                p["Abs{cT}"] = 2.0;
                p["Abs{cT5}"] = -43.0;

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s(false);

                TEST_CHECK_RELATIVE_ERROR(abs(wc.c7()), 1.008, eps);

                /* C_P should be ignored, and always equal -C_S */
                TEST_CHECK_RELATIVE_ERROR(abs(wc.cS()), 42, eps);
                TEST_CHECK_RELATIVE_ERROR(abs(wc.cP()), 42, eps);
                TEST_CHECK_RELATIVE_ERROR(real(wc.cP()), -real(wc.cS()), eps);
                TEST_CHECK_RELATIVE_ERROR(imag(wc.cP()), -imag(wc.cS()), eps);

                TEST_CHECK_RELATIVE_ERROR(abs(wc.cSprime()), 3.2, eps);
                TEST_CHECK_RELATIVE_ERROR(arg(wc.cSprime()), 1.2, eps);
                TEST_CHECK_RELATIVE_ERROR(abs(wc.cPprime()), 3.2, eps);
                TEST_CHECK_RELATIVE_ERROR(arg(wc.cPprime()), 1.2, eps);

                /* C_T and C_T5 vanish */
                TEST_CHECK_NEARLY_EQUAL(std::real(wc.cT()), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(wc.cT()), 0.0, eps);

                TEST_CHECK_NEARLY_EQUAL(std::real(wc.cT5()), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(wc.cT5()), 0.0, eps);

                /* Used parameters registered */
                TEST_CHECK(std::find(std::begin(model), std::end(model), p["Abs{cS}"].id()) != std::end(model));
                TEST_CHECK(std::find(std::begin(model), std::end(model), p["Arg{cS}"].id()) != std::end(model));
                TEST_CHECK(std::find(std::begin(model), std::end(model), p["Re{cS}"].id())  == std::end(model));
                TEST_CHECK(std::find(std::begin(model), std::end(model), p["Im{cS}"].id())  == std::end(model));

                std::list<Parameter::Id> unused_ids = {
                    p["Abs{cP}"].id(),
                    p["Arg{cP}"].id(),
                    p["Abs{cP'}"].id(),
                    p["Arg{cP'}"].id(),
                    p["Abs{cT}"].id(),
                    p["Arg{cT}"].id(),
                    p["Abs{cT5}"].id(),
                    p["Arg{cT5}"].id(),
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
                o.set("scan-mode", "cartesian");

                p["Re{c7}"] = 1.008;
                p["Re{cS}"] = 42;   p["Re{cP}"] = -1.0 * p["Re{cS}"]();
                p["Im{cS'}"] = -12; p["Im{cP'}"] = p["Im{cS'}"]();
                p["Re{cT}"] = 0.0;  p["Im{cT}"] = 0.0;
                p["Re{cT5}"] = 0.0; p["Im{cT5}"] = 0.0;

                ConstrainedWilsonScanModel constrained_model(p, o);
                WilsonScanModel unconstrained_model(p, o);

                WilsonCoefficients<BToS> constrained_wc = constrained_model.wilson_coefficients_b_to_s(false);
                WilsonCoefficients<BToS> unconstrained_wc = constrained_model.wilson_coefficients_b_to_s(false);

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
