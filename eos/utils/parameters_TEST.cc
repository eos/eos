/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2024 Danny van Dyk
 * Copyright (c) 2021 Philip Lüghausen
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

#include <eos/utils/parameters.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class ParametersTest : public TestCase
{
    public:
        ParametersTest() :
            TestCase("parameters_test")
        {
        }

        virtual void
        run() const
        {
            // Setting and retrieval
            {
                Parameters original = Parameters::Defaults();
                Parameter  m_c      = original["mass::c"];

                TEST_CHECK_EQUAL(m_c(), m_c.central());

                m_c = 0.0;
                TEST_CHECK_EQUAL(m_c(), 0.0);

                m_c = m_c.central();
                TEST_CHECK_EQUAL(m_c(), m_c.central());
            }

            // Declaring a new parameter
            {
                Parameters::declare("mass::boeing747", R"(\text{Boeing 747})", Unit::Undefined(), 100000.0, 90000.0, 110000.0);
                Parameters parameters    = Parameters::Defaults();
                Parameter  new_parameter = parameters["mass::boeing747"];

                TEST_CHECK_EQUAL(new_parameter.name(), "mass::boeing747");
                TEST_CHECK_EQUAL(new_parameter.latex(), R"(\text{Boeing 747})");
                TEST_CHECK_EQUAL(new_parameter.unit(), Unit::Undefined());
                TEST_CHECK_EQUAL(new_parameter.evaluate(), 100000.0);
                TEST_CHECK_EQUAL(new_parameter.min(), 90000.0);
                TEST_CHECK_EQUAL(new_parameter.max(), 110000.0);
            }

            // Cloning
            {
                Parameters original = Parameters::Defaults();
                Parameters clone    = original.clone();

                Parameter m_c_original = original["mass::c"];
                Parameter m_c_clone    = clone["mass::c"];

                TEST_CHECK_EQUAL(m_c_original(), m_c_original.central());
                TEST_CHECK_EQUAL(m_c_clone(), m_c_clone.central());

                m_c_clone = 0.0;
                TEST_CHECK_EQUAL(m_c_original(), m_c_original.central());
                TEST_CHECK_EQUAL(m_c_clone(), 0.0);

                m_c_clone = m_c_clone.central();
                TEST_CHECK_EQUAL(m_c_original(), m_c_original.central());
                TEST_CHECK_EQUAL(m_c_clone(), m_c_clone.central());

                m_c_original = 0.0;
                TEST_CHECK_EQUAL(m_c_original(), 0.0);
                TEST_CHECK_EQUAL(m_c_clone(), m_c_clone.central());
            }

            // Parameters::has
            {
                Parameters p = Parameters::Defaults();

                TEST_CHECK_EQUAL(p.has("mass::tau"), true);
                TEST_CHECK_EQUAL(p.has("mass::boing747"), false);
            }

            // Parameters::declare_and_insert
            {
                Parameters p = Parameters::Defaults();

                TEST_CHECK_EQUAL(p.has("mass::boing747"), false);

                p.declare_and_insert("mass::boing747", R"(\text{Boeing 747})", Unit::Undefined(), 100000.0, 90000.0, 110000.0);

                TEST_CHECK_EQUAL(p.has("mass::boing747"), true);
            }

            // Parameters::redirect
            {
                Parameters p = Parameters::Defaults();

                Parameter     p_tau  = p["ubtaunutau::Re{cVL}"];
                Parameter::Id id_tau = p_tau.id();
                p_tau.set(-9.87);

                TEST_CHECK_NEARLY_EQUAL(p_tau.evaluate(), -9.87, 1e-12);

                Parameter     p_ell  = p.declare_and_insert("ublnul::Re{cVL}", R"(\text{Re} C_{V_L}^{ub\ell\nu_\ell})", Unit::None(), 1.23, -1.0, 1.0);
                Parameter::Id id_ell = p_ell.id();

                TEST_CHECK_NEARLY_EQUAL(p_ell.evaluate(), +1.23, 1e-12);

                // redirect the tau parameter name to the lepton-flavor universal parameter
                p.redirect_and_apply("ubtaunutau::Re{cVL}", id_ell);

                // check that the old tau parameter still has the old parameter id
                TEST_CHECK_EQUAL(p_tau.id(), id_tau);

                // re-access the tau parameter
                p_tau = p["ubtaunutau::Re{cVL}"];

                // check that the tau parameter now has the value and id of the lepton-flavor universal parameter
                TEST_CHECK_NEARLY_EQUAL(p_tau.evaluate(), +1.23, 1e-12);
                TEST_CHECK_EQUAL(p_tau.id(), id_ell);

                // check a new Parameters object
                {
                    Parameters p2     = Parameters::Defaults();
                    Parameter  p2_tau = p2["ubtaunutau::Re{cVL}"];
                    TEST_CHECK_NEARLY_EQUAL(p2_tau.evaluate(), +1.23, 1e-12);
                }

                // undo the redirect
                p.redirect_and_apply("ubtaunutau::Re{cVL}", id_tau);

                // re-access the tau parameter
                p_tau = p["ubtaunutau::Re{cVL}"];

                // check that the tau parameter now has the value and id as at the beginning
                TEST_CHECK_NEARLY_EQUAL(p_tau.evaluate(), -9.87, 1e-12);
                TEST_CHECK_EQUAL(p_tau.id(), id_tau);

                // check a new Parameters object
                {
                    Parameters p2     = Parameters::Defaults();
                    Parameter  p2_tau = p2["ubtaunutau::Re{cVL}"];
                    TEST_CHECK_NEARLY_EQUAL(p2_tau.evaluate(), +1.00, 1e-12); // default value
                }
            }

            // Loading of templated parameters
            {
                Parameters p = Parameters::Defaults();

                // 'B->pi::a^%1%_%2%@G2026' is templated over the cartesian product of
                //   - %1% in { "f+", "f0", "fT" }
                //   - %2% in { 0, 1, 2, 3, 4 }
                // with the latex template '$a_%2%^{%1%,B \to \pi,\mathrm{G2026}}$' and the
                // latex substitutions { "f+" -> "f_+", "f0" -> "f_0", "fT" -> "f_T" }.

                // all instances of the cartesian product must be present, ...
                for (const std::string & ff : { "f+", "f0", "fT" })
                {
                    for (unsigned i = 0; i <= 4; ++i)
                    {
                        const std::string name = "B->pi::a^" + ff + "_" + std::to_string(i) + "@G2026";
                        TEST_CHECK_MSG(p.has(name), "templated parameter '" + name + "' is missing");
                    }
                }

                // ... while instances outside of the cartesian product must be absent.
                // (the raw template name 'a^%1%_%2%' cannot be queried directly, since '%'
                // is not a valid character in a QualifiedName.)
                TEST_CHECK(! p.has("B->pi::a^f+_5@G2026"));
                TEST_CHECK(! p.has("B->pi::a^fX_0@G2026"));

                // the positional substitutions (note: %2% precedes %1% in the latex template)
                // and the latex map must be applied correctly.
                Parameter a_fp_0 = p["B->pi::a^f+_0@G2026"];
                TEST_CHECK_EQUAL(a_fp_0.name(), "B->pi::a^f+_0@G2026");
                TEST_CHECK_EQUAL(a_fp_0.latex(), R"($a_0^{f_+,B \to \pi,\mathrm{G2026}}$)");
                TEST_CHECK_NEARLY_EQUAL(a_fp_0.central(), 0.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(a_fp_0.min(), -1.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(a_fp_0.max(), +1.0, 1e-12);

                Parameter a_fT_2 = p["B->pi::a^fT_2@G2026"];
                TEST_CHECK_EQUAL(a_fT_2.name(), "B->pi::a^fT_2@G2026");
                TEST_CHECK_EQUAL(a_fT_2.latex(), R"($a_2^{f_T,B \to \pi,\mathrm{G2026}}$)");
            }
        }
} parameters_test;
