/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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
#include <eos/utils/parameters.hh>

using namespace test;
using namespace eos;

class ParametersTest :
    public TestCase
{
    public:
        ParametersTest() :
            TestCase("parameters_test")
        {
        }

        virtual void run() const
        {
            // Setting and retrieval
            {
                Parameters original = Parameters::Defaults();
                Parameter m_c = original["mass::c"];

                TEST_CHECK_EQUAL(m_c(), m_c.central());

                m_c = 0.0;
                TEST_CHECK_EQUAL(m_c(), 0.0);

                m_c = m_c.central();
                TEST_CHECK_EQUAL(m_c(), m_c.central());
            }

            // Cloning
            {
                Parameters original = Parameters::Defaults();
                Parameters clone = original.clone();

                Parameter m_c_original = original["mass::c"];
                Parameter m_c_clone = clone["mass::c"];

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
        }
} parameters_test;
