/*
 * Copyright (c) 2023-2026 Danny van Dyk
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

#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/observable.hh>
#include <eos/utils/observable_cache.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class PToPFormFactorsTest : public TestCase
{
    public:
        PToPFormFactorsTest() :
            TestCase("p_to_p_form_factor_test")
        {
        }

        virtual void
        run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToP>::create("Foo->Bar::BSZ2015", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToP>::create("B->pi::FooBar", parameter, options));
            }
        }
} p_to_p_form_factor_test;

class PToPPFormFactorsTest : public TestCase
{
    public:
        PToPPFormFactorsTest() :
            TestCase("p_to_pp_form_factor_test")
        {
        }

        virtual void
        run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToPP>::create("Foo->BarBaz::FvDV2018", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToPP>::create("B->pipi::BazBar", parameter, options));
            }
        }
} p_to_pp_form_factor_test;

class PToVFormFactorsTest : public TestCase
{
    public:
        PToVFormFactorsTest() :
            TestCase("p_to_v_form_factor_test")
        {
        }

        virtual void
        run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToV>::create("Foo->Baz::BSZ2015", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToV>::create("B->rho::FooBaz", parameter, options));
            }
        }
} p_to_v_form_factor_test;

class PToGammaFormFactorsTest : public TestCase
{
    public:
        PToGammaFormFactorsTest() :
            TestCase("p_to_gamma_form_factor_test")
        {
        }

        virtual void
        run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToGamma>::create("Foo->gluon::FLvD2022QCDF", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToGamma>::create("B->gamma::FooBaz", parameter, options));
            }
        }
} p_to_gamma_form_factor_test;

class PToGammaOffShellFormFactorsTest : public TestCase
{
    public:
        PToGammaOffShellFormFactorsTest() :
            TestCase("p_to_gamma_off_shell_form_factor_test")
        {
        }

        virtual void
        run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToGammaOffShell>::create("Foo->gluon^*::KKvDZ2022", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<PToGammaOffShell>::create("B->gamma^*::FooBaz", parameter, options));
            }
        }
} p_to_gamma_off_shell_form_factor_test;

class VToPFormFactorsTest : public TestCase
{
    public:
        VToPFormFactorsTest() :
            TestCase("v_to_p_form_factor_test")
        {
        }

        virtual void
        run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<VToP>::create("Foo->Baz::BGJvD2019", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<VToP>::create("B^*->D::FooBaz", parameter, options));
            }
        }
} v_to_p_form_factor_test;

class VToVFormFactorsTest : public TestCase
{
    public:
        VToVFormFactorsTest() :
            TestCase("v_to_v_form_factor_test")
        {
        }

        virtual void
        run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<VToV>::create("Foo->Baz::BGJvD2019", parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<VToV>::create("B^*->D^*::FooBaz", parameter, options));
            }
        }
} v_to_v_form_factor_test;

class VacuumToPPCacheableObservableTest : public TestCase
{
    public:
        VacuumToPPCacheableObservableTest() :
            TestCase("vacuum_to_pp_cacheable_observable_test")
        {
        }

        virtual void
        run() const
        {
            // observables that share an intermediate result agree with observables that do not
            {
                static const double eps = 1.0e-14;

                Parameters p                        = Parameters::Defaults();
                p["mass::pi^+"]                     = 0.13957;
                p["mass::K_d"]                      = 0.497611;
                p["0->pipi::t_0@KKRvD2024"]         = -1.0;
                p["0->pipi::b_(+,1)^2@KKRvD2024"]   = -0.0182238;
                p["0->pipi::b_(+,1)^3@KKRvD2024"]   = -0.0225337;
                p["0->pipi::M_(+,1)@KKRvD2024"]     = 0.760895;
                p["0->pipi::Gamma_(+,1)@KKRvD2024"] = 0.146155;
                p["0->Kpi::t_0@KSvD2025"]           = -1.0;
                p["0->Kpi::M_(+,0)@KSvD2025"]       = 0.890;
                p["0->Kpi::Gamma_(+,0)@KSvD2025"]   = 0.026;
                p["0->Kpi::M_(+,1)@KSvD2025"]       = 1.368;
                p["0->Kpi::Gamma_(+,1)@KSvD2025"]   = 0.106;
                p["0->Kpi::b_+^1@KSvD2025"]         = 0.1;
                p["0->Kpi::b_+^2@KSvD2025"]         = 0.05;
                p["0->Kpi::b_+^3@KSvD2025"]         = 0.01;

                Options o_pipi{
                    { "form-factors"_ok, "KKRvD2024"_ov }
                };
                Options o_Kpi{
                    { "form-factors"_ok, "KSvD2025"_ov }
                };

                Kinematics k_1{
                    { "q2", 0.1 }
                };
                Kinematics k_2{
                    { "q2", 0.5 }
                };
                Kinematics k_3{
                    { "Re{q2}", 0.5 },
                    { "Im{q2}", 0.1 }
                };

                const std::array<std::tuple<QualifiedName, Kinematics, Options>, 5> observables{
                    std::make_tuple(QualifiedName("0->pipi::Abs{f_+}^2(q2)"), k_1, o_pipi), std::make_tuple(QualifiedName("0->pipi::Arg{f_+}(q2)"), k_1, o_pipi),
                    std::make_tuple(QualifiedName("0->pipi::Abs{f_+}^2(q2)"), k_2, o_pipi), std::make_tuple(QualifiedName("0->pipi::Re{f_+}(Re{q2},Im{q2})"), k_3, o_pipi),
                    std::make_tuple(QualifiedName("0->Kpi::Abs{f_+}^2(q2)"), k_1, o_Kpi),
                };

                ObservableCache                            cache(p);
                std::vector<ObservableCache::ObservableId> ids;
                std::vector<double>                        expectations;

                for (const auto & [name, kinematics, options] : observables)
                {
                    ObservablePtr observable = Observable::make(name, p, kinematics, options);

                    TEST_CHECK(nullptr != observable);

                    ids.push_back(cache.add(observable));
                    expectations.push_back(Observable::make(name, p, kinematics, options)->evaluate());
                }

                TEST_CHECK_NO_THROW(cache.update());

                for (unsigned i = 0; i < ids.size(); ++i)
                {
                    TEST_CHECK_RELATIVE_ERROR(cache[ids[i]], expectations[i], eps);
                }
            }
        }
} vacuum_to_pp_cacheable_observable_test;
