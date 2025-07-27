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
        }
} parameters_test;
