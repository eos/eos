/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2015 Danny van Dyk
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

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s("mu", false);
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
                p["mu"] = mu;
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
                o.set("scan-mode", "cartesian");

                WilsonScanModel model(p, o);

                TEST_CHECK_NEARLY_EQUAL(+0.2233419372, model.alpha_s(mu), eps);

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s("mu", false);
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

                wc = model.wilson_coefficients_b_to_s("e", false);
                TEST_CHECK_NEARLY_EQUAL(+3.27,       real(wc.c9()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.007,      real(wc.c9prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(+0.006,      real(wc.c10prime()), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.01,       imag(wc.c9prime()),  eps);
                TEST_CHECK_NEARLY_EQUAL(-M_PI+0.01,  imag(wc.c10prime()), eps);
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
                ConstrainedWilsonScanModel model(p, o);

                p["b->s::Re{c7}"] = 1.008;
                p["b->smumu::Re{cS}"] = 42;
                p["b->smumu::Re{cP}"] = 100;
                p["b->smumu::Im{cS'}"] = -12;
                p["b->smumu::Im{cP'}"] = -135;
                p["b->smumu::Re{cT}"] = 2.0;
                p["b->smumu::Re{cT5}"] = -43.0;

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s("mu", false);

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

                WilsonCoefficients<BToS> wc = model.wilson_coefficients_b_to_s("mu", false);

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
                o.set("scan-mode", "cartesian");

                p["b->s::Re{c7}"] = 1.008;
                p["b->smumu::Re{cS}"] = 42;   p["b->smumu::Re{cP}"] = -1.0 * p["b->smumu::Re{cS}"]();
                p["b->smumu::Im{cS'}"] = -12; p["b->smumu::Im{cP'}"] = p["b->smumu::Im{cS'}"]();
                p["b->smumu::Re{cT}"] = 0.0;  p["b->smumu::Im{cT}"] = 0.0;
                p["b->smumu::Re{cT5}"] = 0.0; p["b->smumu::Im{cT5}"] = 0.0;

                ConstrainedWilsonScanModel constrained_model(p, o);
                WilsonScanModel unconstrained_model(p, o);

                WilsonCoefficients<BToS> constrained_wc = constrained_model.wilson_coefficients_b_to_s("mu", false);
                WilsonCoefficients<BToS> unconstrained_wc = constrained_model.wilson_coefficients_b_to_s("mu", false);

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
