/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Ahmet Kokulu
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
#include <eos/observable.hh>
#include <eos/b-decays/b-to-d-l-nu.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <vector>


using namespace test;
using namespace eos;

class BToDLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToDLeptonNeutrinoTest() :
            TestCase("b_to_d_l_nu_test")
        {
        }

        virtual void run() const
        {

            // tests for normalized SM Br (|V_cb|=1)
            {
                Parameters p1 = Parameters::Defaults();
                p1["B->D::f_+(0)@BCL2008"]  = +0.660;
                p1["B->D::b_+^1@BCL2008"]   = -4.00;
                p1["B->D::b_+^2@BCL2008"]   = -0.80;
                p1["B->D::b_0^1@BCL2008"]   = +0.40;
                p1["B->D::b_0^2@BCL2008"]   = -1.20;
                p1["mass::B_d"]             = 5.279;
                p1["mass::D_d"]             = 1.870;
                p1["b->cmunumu::Re{cVL}"]   = 1.0;
                
                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BCL2008");
                oo.set("q", "d");
                
                BToDLeptonNeutrino d(p1, oo);
                
                const double eps = 1e-5;
                
                // test for different s-bins - the default lepton is muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(4.0, 8.0), 4.413884800954554, eps);
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(8.0, 11.62), 1.0206050186357505, eps);
                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(0.011164, 11.62), 13.02629683128414, eps);
            }

            // tests for normalized NP Br (|V_cb|=1)
            {
                Parameters p2 = Parameters::Defaults();
                p2["B->D::f_+(0)@BCL2008"]  = +0.660;
                p2["B->D::b_+^1@BCL2008"]   = -4.00;
                p2["B->D::b_+^2@BCL2008"]   = -0.800;
                p2["B->D::b_0^1@BCL2008"]   = +0.400;
                p2["B->D::b_0^2@BCL2008"]   = -1.20;
                p2["mass::B_d"]             = 5.279;
                p2["mass::D_d"]             = 1.870;
                // fix the scale
                p2["mu"]                    = 4.18;
                p2["mass::b(MSbar)"]        = 4.18;
                p2["mass::c"]               = 1.275;
                p2["b->cmunumu::Re{cVL}"]   =  1.0;
                p2["b->cmunumu::Im{cVL}"]   = -1.0;
                p2["b->cmunumu::Re{cVR}"]   =  2.0;
                p2["b->cmunumu::Im{cVR}"]   = -2.0;
                p2["b->cmunumu::Re{cSL}"]   =  3.0;
                p2["b->cmunumu::Im{cSL}"]   = -3.0;
                p2["b->cmunumu::Re{cSR}"]   =  4.0;
                p2["b->cmunumu::Im{cSR}"]   = -4.0;
                
                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BCL2008");
                
                BToDLeptonNeutrino d(p2, oo);
                
                const double eps = 1e-5;
                
                // test for different s-bins - the default lepton is muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(4.0, 8.0), 749.0934190079766, eps);
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(8.0, 11.62), 616.8445050995084, eps);
                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(0.011164, 11.62), 1774.909629523804, eps);
            }
        }
} b_to_d_l_nu_test;
