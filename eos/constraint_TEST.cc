/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2014 Danny van Dyk
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
#include <eos/utils/log_likelihood.hh>

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

                    // BaBar
                    "B->X_sll::BR[1.0,6.0]@BaBar-2004A",

                    // Belle
                    "B^0->K^*0gamma::BR@Belle-2004",
                    "B^+->K^*+gamma::BR@Belle-2004",


                    /* 2005 */

                    // Belle
                    "B->X_sll::BR[1.0,6.0]@Belle-2005A",


                    /* 2006 */

                    // Belle
                    "B^0->K^*0gamma::S_K@Belle-2006",
                    "B^0->K^*0gamma::C_K@Belle-2006",


                    /* 2008 */

                    // BaBar
                    "B^0->K^*0gamma::S_K@BaBar-2008",
                    "B^0->K^*0gamma::C_K@BaBar-2008",

                    // Belle
                    "B->X_sgamma::E_1[1.8]@Belle-2008",
                    "B->X_sgamma::E_2[1.8]@Belle-2008",
                    "B->X_sgamma::E_1[1.8]+E_2[1.8]@Belle-2008",


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

                    // Belle
                    // B -> X_s gamma
                    "B->X_sgamma::BR[1.8]@Belle-2009B",
                    "B->X_sgamma::BR[1.8]+E_1[1.8]+E_2[1.8]@Belle-2009B",


                    /* 2010 */
                    // BaBar
                    "B^0->pi^+lnu::BR[0.0,4.0]@BaBar-2010A",
                    "B^0->pi^+lnu::BR[4.0,8.0]@BaBar-2010A",
                    "B^0->pi^+lnu::BR[8.0,12.0]@BaBar-2010A",

                    "B^0->pi^+lnu::BR@BaBar-2010B",

                    // Belle
                    "B^0->pi^+lnu::BR@Belle-2010A",

                    /* 2011 */

                    // HFAG
                    "B^0->K^*0gamma::S_K+C_K@HFAG-2011",

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
                    "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.21]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@CDF-2011",
                    "B^0->K^*0mu^+mu^-::A_im[16.00,19.21]@CDF-2011",
                    // B^0 -> K^0 mu^+ mu^-
                    "B^0->K^0mu^+mu^-::BR[1.00,6.00]@CDF-2011",
                    "B^0->K^0mu^+mu^-::BR[14.18,16.00]@CDF-2011",
                    "B^0->K^0mu^+mu^-::BR[16.00,22.86]@CDF-2011",
                    // B^+ -> K^+ mu^+ mu^-
                    "B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2011",
                    "B^+->K^+mu^+mu^-::BR[14.18,16.00]@CDF-2011",
                    "B^+->K^+mu^+mu^-::BR[16.00,22.86]@CDF-2011",

                    "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2011",
                    "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2011",
                    "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2011",

                    "B^0_s->mu^+mu^-::BR_limit@CDF-2011",

                    // LHCb
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2011",
                    "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@LHCb-2011",
                    // B^0_s -> mu^+ mu^-
                    "B^0_s->mu^+mu^-::BR_limit@LHCb-CMS-2011",
                    "B^0_s->mu^+mu^-::BR_limit@LHCb-CMS-2011-Bayes",


                    /* 2012 */

                    // BaBar
                    // B^+ -> K^+ mu^+ mu^-
                    "B^+->K^+mu^+mu^-::BR[1.00,6.00]@BaBar-2012",
                    "B^+->K^+mu^+mu^-::BR[14.21,16.00]@BaBar-2012",
                    "B^+->K^+mu^+mu^-::BR[16.00,22.86]@BaBar-2012",
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@BaBar-2012",
                    "B^0->K^*0mu^+mu^-::BR[14.21,16.00]@BaBar-2012",
                    "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@BaBar-2012",
                    "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@BaBar-2012",
                    "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@BaBar-2012",
                    "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@BaBar-2012",
                    "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@BaBar-2012",
                    "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@BaBar-2012",
                    "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@BaBar-2012",
                    // B -> X_s gamma
                    "B->X_sgamma::BR[1.8]@BaBar-2012",
                    "B->X_sgamma::E_1[1.8]@BaBar-2012",
                    "B->X_sgamma::E_2[1.8]@BaBar-2012",
                    // B^0 -> pi l nu
                    "B^0->pi^+lnu::BR@BaBar-2012D",

                    // CDF
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@CDF-2012",
                    // B^+ -> K^*+ mu^+ mu^-
                    "B^+->K^*+mu^+mu^-::BR[1.00,6.00]@CDF-2012",
                    "B^+->K^*+mu^+mu^-::BR[14.18,16.00]@CDF-2012",
                    "B^+->K^*+mu^+mu^-::BR[16.00,19.21]@CDF-2012",
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.21]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@CDF-2012",
                    "B^0->K^*0mu^+mu^-::A_im[16.00,19.21]@CDF-2012",
                    // B^0 -> K^0 mu^+ mu^-
                    "B^0->K^0mu^+mu^-::BR[1.00,6.00]@CDF-2012",
                    "B^0->K^0mu^+mu^-::BR[14.18,16.00]@CDF-2012",
                    "B^0->K^0mu^+mu^-::BR[16.00,22.86]@CDF-2012",
                    // B^+ -> K^+ mu^+ mu^-
                    "B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2012",
                    "B^+->K^+mu^+mu^-::BR[14.18,16.00]@CDF-2012",
                    "B^+->K^+mu^+mu^-::BR[16.00,22.86]@CDF-2012",

                    "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2012",
                    "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2012",
                    "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2012",

                    // LHCb
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::A_im[16.00,19.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::S_3[1.00,6.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::S_3[14.18,16.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::S_3[16.00,19.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::A_CP[1.00,6.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::A_CP[14.18,16.00]@LHCb-2012",
                    "B^0->K^*0mu^+mu^-::A_CP[16.00,20.00]@LHCb-2012",

                    // B^0_s -> mu^+ mu^-
                    "B^0_s->mu^+mu^-::BR_limit@LHCb-2012",
                    "B^0_s->mu^+mu^-::BR_limit@LHCb-Nov-2012",

                    // B^+ -> K^+ mu^+ mu^-
                    "B^+->K^+mu^+mu^-::BR[1.00,6.00]@LHCb-2012",
                    "B^+->K^+mu^+mu^-::BR[14.18,16.00]@LHCb-2012",
                    "B^+->K^+mu^+mu^-::BR[16.00,18.00]@LHCb-2012",
                    "B^+->K^+mu^+mu^-::BR[18.00,22.00]@LHCb-2012",

                    "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@LHCb-2012",
                    "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@LHCb-2012",
                    "B^+->K^+mu^+mu^-::A_FB[16.00,18.00]@LHCb-2012",
                    "B^+->K^+mu^+mu^-::A_FB[18.00,22.00]@LHCb-2012",

                    "B^+->K^+mu^+mu^-::F_H[1.00,6.00]@LHCb-2012",
                    "B^+->K^+mu^+mu^-::F_H[14.18,16.00]@LHCb-2012",
                    "B^+->K^+mu^+mu^-::F_H[16.00,18.00]@LHCb-2012",
                    "B^+->K^+mu^+mu^-::F_H[18.00,22.00]@LHCb-2012",

                    /* 2013 */

                    // ATLAS
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@ATLAS-2013A",
                    "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@ATLAS-2013A",
                    "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@ATLAS-2013A",
                    "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@ATLAS-2013A",
                    "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@ATLAS-2013A",
                    "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@ATLAS-2013A",

                    // Belle
                    "B^0->pi^+lnu::BR@Belle-2013A",

                    // CMS
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CMS-2013A",
                    "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CMS-2013A",
                    "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@CMS-2013A",
                    "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CMS-2013A",
                    "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CMS-2013A",
                    "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@CMS-2013A",
                    "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CMS-2013A",
                    "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CMS-2013A",
                    "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@CMS-2013A",

                    // B^0_s -> mu^+ mu^-
                    "B^0_s->mu^+mu^-::BR@CMS-2013B",

                    // LHCb
                    // B^0 -> K^*0 mu^+ mu^-
                    "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::S_3[1.00,6.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::S_3[14.18,16.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::S_3[16.00,19.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::S_9[1.00,6.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::S_9[14.18,16.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::S_9[16.00,19.00]@LHCb-2013",
                    // The following commented observables have not yet been implemented!
#if 0
                    "B^0->K^*0mu^+mu^-::A_9[1.00,6.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::A_9[14.18,16.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::A_9[16.00,19.00]@LHCb-2013",
#endif
                    "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::A_T^re[1.00,6.00]@LHCb-2013",
                    // The following constraint does not conform to a gaussian likelihood.
#if 0
                    "B^0->K^*0mu^+mu^-::A_T^re[14.18,16.00]@LHCb-2013",
#endif
                    "B^0->K^*0mu^+mu^-::A_T^re[16.00,19.00]@LHCb-2013",

                    "B^0->K^*0mu^+mu^-::P'_4[1.00,6.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::P'_4[14.18,16.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::P'_4[16.00,19.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::P'_5[1.00,6.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::P'_5[14.18,16.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::P'_5[16.00,19.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::P'_6[1.00,6.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::P'_6[14.18,16.00]@LHCb-2013",
                    "B^0->K^*0mu^+mu^-::P'_6[16.00,19.00]@LHCb-2013",

                    // B^0_s -> mu^+ mu^-
                    "B^0_s->mu^+mu^-::BR@LHCb-2013D",

                    /* Theory Constraints */
                    /* 2013 */
                    // disabled, since it needs the option "form-factors" specified.
                    //"B->K::f_+@HPQCD-2013A",

                    //"B->K^*::V@HPQCD-2013B",
                    //"B->K^*::A_1@HPQCD-2013B",
                    //"B->K^*::A_12V@HPQCD-2013B",


                    // LHCb
                    // B^+ -> K^+ mu^+ mu^-
                    "B^+->K^+mu^+mu^-::BR[1.10,2.00]@LHCb-2014",
                    "B^+->K^+mu^+mu^-::BR[2.00,3.00]@LHCb-2014",
                    "B^+->K^+mu^+mu^-::BR[3.00,4.00]@LHCb-2014",
                    "B^+->K^+mu^+mu^-::BR[4.00,5.00]@LHCb-2014",
                    "B^+->K^+mu^+mu^-::BR[5.00,6.00]@LHCb-2014",
                    "B^+->K^+mu^+mu^-::BR[1.10,6.00]@LHCb-2014",
                    "B^+->K^+mu^+mu^-::BR[15.00,22.00]@LHCb-2014",

                    "B^+->K^+mu^+mu^-::A_FB[1.10,6.00]@LHCb-2014",
                    "B^+->K^+mu^+mu^-::A_FB[15.00,22.00]@LHCb-2014",
                    "B^+->K^+mu^+mu^-::F_H[1.10,6.00]@LHCb-2014",
                    "B^+->K^+mu^+mu^-::F_H[15.00,22.00]@LHCb-2014",
                };

                std::cout << "# Constraints :" << std::endl;

                for (auto n = constraint_names.cbegin(), n_end = constraint_names.cend() ; n != n_end ; ++n)
                {
                    std::cout << "#  " << *n << ": ";

                    Constraint c = Constraint::make(*n, Options());
                    TEST_CHECK_EQUAL(c.name(), *n);
                    TEST_CHECK(std::distance(c.begin_observables(), c.end_observables()) > 0);
                    TEST_CHECK(std::distance(c.begin_blocks(), c.end_blocks()) > 0);

                    for (auto o = c.begin_observables(), o_end = c.end_observables(); o != o_end ; ++o)
                    {
                        std::cout << (**o).name() << '['
                                << (**o).kinematics().as_string() << ']'
                                << " with options: " << (**o).options().as_string();
                    }
                    for (auto b = c.begin_blocks(), b_end = c.end_blocks(); b != b_end ; ++b)
                    {
                        std::cout << ", " << (**b).as_string();
                    }
                    std::cout << std::endl;
                }
            }
        }
} constraint_test;
