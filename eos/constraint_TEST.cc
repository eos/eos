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
#include <eos/constraint.hh>

#include <iostream>
#include <vector>

using namespace test;
using namespace eos;

class ConstraintTest :
    public TestCase
{
    public:
        ConstraintTest() :
            TestCase("constraint_test")
        {
        }

        virtual void run() const
        {
            /* Test making constraints */
            {
                std::vector<std::string> constraint_names
                {
                    /* 2000 */
                    // CLEO
                    "B^0->K^*0gamma::BR@CLEO-2000",
                    "B^+->K^*+gamma::BR@CLEO-2000",
                    /* 2004 */
                    // Belle
                    "B^0->K^*0gamma::BR@Belle-2004",
                    "B^+->K^*+gamma::BR@Belle-2004",
                    /* 2006 */
                    // Belle
                    "B^0->K^*0gamma::S_K@Belle-2006",
                    "B^0->K^*0gamma::C_K@Belle-2006",
                    /* 2008 */
                    // BaBar
                    "B^0->K^*0gamma::S_K@BaBar-2008",
                    "B^0->K^*0gamma::C_K@BaBar-2008",
                    /* 2009 */
                    // BaBar
                    "B^0->K^*0gamma::BR@BaBar-2009",
                    "B^+->K^*+gamma::BR@BaBar-2009",
                    // Belle
                    // B^+ -> K^+ mu^+ mu^-
                    "B^+->K^+mu^+mu^-::BR[1.00,6.00]@Belle-2009",
                    "B^+->K^+mu^+mu^-::BR[14.18,16.00]@Belle-2009",
                    "B^+->K^+mu^+mu^-::BR[16.00,22.86]@Belle-2009",
                    /* The following commented observables have not yet been implemented! */
                    //"B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@Belle-2009",
                    //"B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@Belle-2009",
                    //"B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@Belle-2009",
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@Belle-2009",
                    "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@Belle-2009",
                    "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@Belle-2009",
                    "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@Belle-2009",
                    "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@Belle-2009",
                    "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@Belle-2009",
                    "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@Belle-2009",
                    "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@Belle-2009",
                    "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@Belle-2009",
                    /* 2011 */
                    // CDF
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@CDF-2011",
                    // B^+ -> K^*+ mu^+ mu^-
                    "B^+->K^*+mu^+mu^-::BR[1.00,6.00]@CDF-2011",
                    "B^+->K^*+mu^+mu^-::BR[14.18,16.00]@CDF-2011",
                    "B^+->K^*+mu^+mu^-::BR[16.00,19.21]@CDF-2011",
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::A_T_2[1.00,6.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::A_T_2[14.18,16.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::A_T_2[16.00,19.21]@CDF-2011",
                    // B^0 -> K^0 mu^+ mu^-
                    "B^0->K^0mu^+mu^-::BR[1.00,6.00]@CDF-2011",
                    "B^0->K^0mu^+mu^-::BR[14.18,16.00]@CDF-2011",
                    "B^0->K^0mu^+mu^-::BR[16.00,22.86]@CDF-2011",
                    // B^+ -> K^+ mu^+ mu^-
                    "B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2011",
                    "B^+->K^+mu^+mu^-::BR[14.18,16.00]@CDF-2011",
                    "B^+->K^+mu^+mu^-::BR[16.00,22.86]@CDF-2011",
                    /* The following commented observables have not yet been implemented! */
                    //"B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2011",
                    //"B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2011",
                    //"B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2011",
                    // LHCb
                    "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@LHCb-2011",
                };

                for (auto n = constraint_names.cbegin(), n_end = constraint_names.cend() ; n != n_end ; ++n)
                {
                    Constraint c = Constraint::make(*n, Options());
                    TEST_CHECK_EQUAL(c.name(), *n);
                    TEST_CHECK(std::distance(c.begin_observables(), c.end_observables()) > 0);
                    TEST_CHECK(std::distance(c.begin_blocks(), c.end_blocks()) > 0);
                }
            }
        }
} constraint_test;
