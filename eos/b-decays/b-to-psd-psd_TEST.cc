/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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
#include <eos/b-decays/b-to-psd-psd.hh>
#include <eos/maths/complex.hh>

using namespace test;
using namespace eos;

class BToPseudoscalarPseudoscalarTest :
    public TestCase
{
    public:
        BToPseudoscalarPseudoscalarTest() :
            TestCase("b_to_psd_psd_test")
        {
        }

        virtual void run() const
        {

            /*Tests with Huber's paper values*/

            {
                Parameters p = Parameters::Defaults();
                p["nonleptonic::Re{AT3}@SU3F"]  =  0.029 * std::cos( -3.083);
                p["nonleptonic::Im{AT3}@SU3F"]  =  0.029 * std::sin( -3.083);
                p["nonleptonic::Re{CT3}@SU3F"]  =  0.258 * std::cos(-0.105);
                p["nonleptonic::Im{CT3}@SU3F"]  =  0.258 * std::sin(-0.105);
                p["nonleptonic::Re{AT6}@SU3F"]  =  0. * std::cos( 0.3);
                p["nonleptonic::Im{AT6}@SU3F"]  =  0. * std::sin( 0.3);
                p["nonleptonic::Re{CT6}@SU3F"]  =  0.235 * std::cos(-0.079);
                p["nonleptonic::Im{CT6}@SU3F"]  =  0.235 * std::sin(-0.079);
                p["nonleptonic::Re{AT15}@SU3F"] =  0.029 * std::cos(-3.083);
                p["nonleptonic::Im{AT15}@SU3F"] =  0.029 * std::sin(-3.083);
                p["nonleptonic::Re{CT15}@SU3F"] =  0.151 * std::cos(0.061);
                p["nonleptonic::Im{CT15}@SU3F"] =  0.151 * std::sin(0.061);
                p["nonleptonic::Re{BT3}@SU3F"]  =  0.034 * std::cos(3.087);
                p["nonleptonic::Im{BT3}@SU3F"]  =  0.034 * std::sin(3.087);
                p["nonleptonic::Re{BT6}@SU3F"]  =  0.033 * std::cos(-0.286);
                p["nonleptonic::Im{BT6}@SU3F"]  =  0.033 * std::sin(-0.286);
                p["nonleptonic::Re{BT15}@SU3F"] =  0.008 * std::cos(-1.892);
                p["nonleptonic::Im{BT15}@SU3F"] =  0.008 * std::sin(-1.892);
                p["nonleptonic::Re{DT3}@SU3F"]  =  0.055 * std::cos(2.942);
                p["nonleptonic::Im{DT3}@SU3F"]  =  0.055 * std::sin(2.942);
                p["nonleptonic::Re{AP3}@SU3F"]  =  0.014 * std::cos(-1.328);
                p["nonleptonic::Im{AP3}@SU3F"]  =  0.014 * std::sin(-1.328);
                p["nonleptonic::Re{CP3}@SU3F"]  =  0.008 * std::cos(0.);
                p["nonleptonic::Im{CP3}@SU3F"]  =  0.008 * std::sin(0.);
                p["nonleptonic::Re{AP6}@SU3F"]  =  0. * std::cos( 1.3);
                p["nonleptonic::Im{AP6}@SU3F"]  =  0. * std::sin( 1.3);
                p["nonleptonic::Re{CP6}@SU3F"]  =  0.145 * std::cos(-2.881);
                p["nonleptonic::Im{CP6}@SU3F"]  =  0.145 * std::sin(-2.881);
                p["nonleptonic::Re{AP15}@SU3F"] =  0.003 * std::cos( 2.234);
                p["nonleptonic::Im{AP15}@SU3F"] =  0.003 * std::sin( 2.234);
                p["nonleptonic::Re{CP15}@SU3F"] =  0.003 * std::cos(-0.608);
                p["nonleptonic::Im{CP15}@SU3F"] =  0.003 * std::sin(-0.608);
                p["nonleptonic::Re{BP3}@SU3F"]  =  0.043 * std::cos( 2.367);
                p["nonleptonic::Im{BP3}@SU3F"]  =  0.043 * std::sin( 2.367);
                p["nonleptonic::Re{BP6}@SU3F"]  =  0.099 * std::cos(0.353);
                p["nonleptonic::Im{BP6}@SU3F"]  =  0.099 * std::sin(0.353);
                p["nonleptonic::Re{BP15}@SU3F"] =  0.031 * std::cos( -0.690);
                p["nonleptonic::Im{BP15}@SU3F"] =  0.031 * std::sin( -0.690);
                p["nonleptonic::Re{DP3}@SU3F"]  =  0.030 * std::cos( 0.477);
                p["nonleptonic::Im{DP3}@SU3F"]  =  0.030 * std::sin( 0.477);
                p["eta::theta_18"]              = -0.3273166181245092;



                static const double eps = 1.0e-6;


                Options o
                {
                    { "representation", "SU3F" },
                    { "q", "d" },
                    { "P1", "eta" },
                    { "P2", "K_d" },
                    { "model", "CKM" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar d(p, o);


                Options oo
                {
                    { "representation", "SU3F" },
                    { "q", "d" },
                    { "P1", "eta" },
                    { "P2", "K_d" },
                    { "model", "CKM" },
                    { "cp-conjugate", "false"}
                };

                BToPseudoscalarPseudoscalar dd(p, oo);



                TEST_CHECK_NEARLY_EQUAL(0.5 * (d.branching_ratio() + dd.branching_ratio()), 1.4164690361731448e-6, eps);

                TEST_CHECK_NEARLY_EQUAL((d.branching_ratio() - dd.branching_ratio()) / (d.branching_ratio() + dd.branching_ratio()), -0.0008872201130703394, eps);
            }

            /*Tests with Topological Amplitudes*/
            {
                Parameters p = Parameters::Defaults();
                p["nonleptonic::Re{T}@Topological"]    =  0.1 * std::cos( 0.1);
                p["nonleptonic::Im{T}@Topological"]    =  0.1 * std::sin( 0.1);
                p["nonleptonic::Re{C}@Topological"]    = -0.2 * std::cos(-0.2);
                p["nonleptonic::Im{C}@Topological"]    = -0.2 * std::sin(-0.2);
                p["nonleptonic::Re{A}@Topological"]    =  0.3 * std::cos( 0.3);
                p["nonleptonic::Im{A}@Topological"]    =  0.3 * std::sin( 0.3);
                p["nonleptonic::Re{E}@Topological"]    = -0.4 * std::cos(-0.4);
                p["nonleptonic::Im{E}@Topological"]    = -0.4 * std::sin(-0.4);
                p["nonleptonic::Re{TES}@Topological"]  =  0.5 * std::cos( 0.5);
                p["nonleptonic::Im{TES}@Topological"]  =  0.5 * std::sin( 0.5);
                p["nonleptonic::Re{TAS}@Topological"]  = -0.6 * std::cos(-0.6);
                p["nonleptonic::Im{TAS}@Topological"]  = -0.6 * std::sin(-0.6);
                p["nonleptonic::Re{TS}@Topological"]   =  0.7 * std::cos( 0.7);
                p["nonleptonic::Im{TS}@Topological"]   =  0.7 * std::sin( 0.7);
                p["nonleptonic::Re{TPA}@Topological"]  = -0.8 * std::cos(-0.8);
                p["nonleptonic::Im{TPA}@Topological"]  = -0.8 * std::sin(-0.8);
                p["nonleptonic::Re{TP}@Topological"]   =  0.9 * std::cos( 0.9);
                p["nonleptonic::Im{TP}@Topological"]   =  0.9 * std::sin( 0.9);
                p["nonleptonic::Re{TSS}@Topological"]  = -1.0 * std::cos(-1.0);
                p["nonleptonic::Im{TSS}@Topological"]  = -1.0 * std::sin(-1.0);
                p["nonleptonic::Re{P}@Topological"]    =  1.1 * std::cos( 1.1);
                p["nonleptonic::Im{P}@Topological"]    =  1.1 * std::sin( 1.1);
                p["nonleptonic::Re{PT}@Topological"]   = -1.2 * std::cos(-1.2);
                p["nonleptonic::Im{PT}@Topological"]   = -1.2 * std::sin(-1.2);
                p["nonleptonic::Re{S}@Topological"]    =  1.3 * std::cos( 1.3);
                p["nonleptonic::Im{S}@Topological"]    =  1.3 * std::sin( 1.3);
                p["nonleptonic::Re{PC}@Topological"]   = -1.4 * std::cos(-1.4);
                p["nonleptonic::Im{PC}@Topological"]   = -1.4 * std::sin(-1.4);
                p["nonleptonic::Re{PTA}@Topological"]  =  1.5 * std::cos( 1.5);
                p["nonleptonic::Im{PTA}@Topological"]  =  1.5 * std::sin( 1.5);
                p["nonleptonic::Re{PA}@Topological"]   = -1.6 * std::cos(-1.6);
                p["nonleptonic::Im{PA}@Topological"]   = -1.6 * std::sin(-1.6);
                p["nonleptonic::Re{PTE}@Topological"]  =  1.7 * std::cos( 1.7);
                p["nonleptonic::Im{PTE}@Topological"]  =  1.7 * std::sin( 1.7);
                p["nonleptonic::Re{PAS}@Topological"]  = -1.8 * std::cos(-1.8);
                p["nonleptonic::Im{PAS}@Topological"]  = -1.8 * std::sin(-1.8);
                p["nonleptonic::Re{PSS}@Topological"]  =  1.9 * std::cos( 1.9);
                p["nonleptonic::Im{PSS}@Topological"]  =  1.9 * std::sin( 1.9);
                p["nonleptonic::Re{PES}@Topological"]  = -2.0 * std::cos(-2.0);
                p["nonleptonic::Im{PES}@Topological"]  = -2.0 * std::sin(-2.0);
                p["eta::theta_18"]                     =  1.0;

                static const double eps = 1.0e-6;

                Options o
                {
                    { "representation", "topological" },
                    { "q", "d" },
                    { "P1", "pi^+" },
                    { "P2", "pi^-" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar d(p, o);

                Options oo
                {
                    { "representation", "topological" },
                    { "q", "d" },
                    { "P1", "pi^+" },
                    { "P2", "pi^-" },
                    { "cp-conjugate", "false"}
                };

                BToPseudoscalarPseudoscalar dd(p, oo);

                TEST_CHECK_NEARLY_EQUAL(0.5 * (d.branching_ratio() + dd.branching_ratio()), 0.0022804835159666867, eps);
                TEST_CHECK_NEARLY_EQUAL((d.branching_ratio() - dd.branching_ratio()) / (d.branching_ratio() + dd.branching_ratio()), 0.09133792735291475, eps);

                Options ooo
                {
                    { "representation", "topological" },
                    { "q", "s" },
                    { "P1", "pi^0" },
                    { "P2", "Kbar_d" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar ddd(p, ooo);

                Options o4
                {
                    { "representation", "topological" },
                     { "q", "s" },
                    { "P1", "pi^0" },
                    { "P2", "Kbar_d" },
                    { "cp-conjugate", "false"}
                };

                BToPseudoscalarPseudoscalar d4(p, o4);

                TEST_CHECK_NEARLY_EQUAL(0.5 * (ddd.branching_ratio() + d4.branching_ratio()), 0.00001936527065018223, eps);
                TEST_CHECK_NEARLY_EQUAL((ddd.branching_ratio() - d4.branching_ratio()) / (ddd.branching_ratio() + d4.branching_ratio()), 0.7366560008503581, eps);

                Options o5
                {
                    { "representation", "topological" },
                    { "q", "u" },
                    { "P1", "eta" },
                    { "P2", "pi^+" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar d5(p, o5);

                Options o6
                {
                    { "representation", "topological" },
                    { "q", "u" },
                    { "P1", "eta" },
                    { "P2", "pi^+" },
                    { "cp-conjugate", "false"}
                };

                BToPseudoscalarPseudoscalar d6(p, o6);

                TEST_CHECK_NEARLY_EQUAL(0.5 * (d5.branching_ratio() + d6.branching_ratio()), 0.0022644592272822654, eps);
                TEST_CHECK_NEARLY_EQUAL((d5.branching_ratio() - d6.branching_ratio()) / (d5.branching_ratio() + d6.branching_ratio()), -0.16086877559501292, eps);

                Options o7
                {
                    { "representation", "topological" },
                    { "q", "d" },
                    { "P1", "eta'" },
                    { "P2", "K_d" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar d7(p, o7);

                Options o8
                {
                    { "representation", "topological" },
                    { "q", "d" },
                    { "P1", "eta'" },
                    { "P2", "K_d" },
                    { "cp-conjugate", "false"}
                };

                BToPseudoscalarPseudoscalar d8(p, o8);

                TEST_CHECK_NEARLY_EQUAL(0.5 * (d7.branching_ratio() + d8.branching_ratio()), 0.00537825569161075, eps);
                TEST_CHECK_NEARLY_EQUAL((d7.branching_ratio() - d8.branching_ratio()) / (d7.branching_ratio() + d8.branching_ratio()), 0.006891941917256833, eps);


                Options o9
                {
                    { "representation", "topological" },
                    { "q", "d" },
                    { "P1", "eta" },
                    { "P2", "eta" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar d9(p, o9);

                Options o10
                {
                    { "representation", "topological" },
                    { "q", "d" },
                    { "P1", "eta" },
                    { "P2", "eta" },
                    { "cp-conjugate", "false"}
                };

                BToPseudoscalarPseudoscalar d10(p, o10);

                TEST_CHECK_NEARLY_EQUAL(0.5 * (d9.branching_ratio() + d10.branching_ratio()), 0.004473515098754955, eps);
                TEST_CHECK_NEARLY_EQUAL((d9.branching_ratio() - d10.branching_ratio()) / (d9.branching_ratio() + d10.branching_ratio()), 0.14535763351380676, eps);
            }

            /*Tests with su3 parameterization */
            {
                Parameters p = Parameters::Defaults();
                p["nonleptonic::Re{AT3}@SU3F"]  =  0.1 * std::cos( 0.1);
                p["nonleptonic::Im{AT3}@SU3F"]  =  0.1 * std::sin( 0.1);
                p["nonleptonic::Re{CT3}@SU3F"]  = -0.2 * std::cos(-0.2);
                p["nonleptonic::Im{CT3}@SU3F"]  = -0.2 * std::sin(-0.2);
                p["nonleptonic::Re{AT6}@SU3F"]  =  0.3 * std::cos( 0.3);
                p["nonleptonic::Im{AT6}@SU3F"]  =  0.3 * std::sin( 0.3);
                p["nonleptonic::Re{CT6}@SU3F"]  = -0.4 * std::cos(-0.4);
                p["nonleptonic::Im{CT6}@SU3F"]  = -0.4 * std::sin(-0.4);
                p["nonleptonic::Re{AT15}@SU3F"] =  0.5 * std::cos( 0.5);
                p["nonleptonic::Im{AT15}@SU3F"] =  0.5 * std::sin( 0.5);
                p["nonleptonic::Re{CT15}@SU3F"] = -0.6 * std::cos(-0.6);
                p["nonleptonic::Im{CT15}@SU3F"] = -0.6 * std::sin(-0.6);
                p["nonleptonic::Re{BT3}@SU3F"]  =  0.7 * std::cos( 0.7);
                p["nonleptonic::Im{BT3}@SU3F"]  =  0.7 * std::sin( 0.7);
                p["nonleptonic::Re{BT6}@SU3F"]  = -0.8 * std::cos(-0.8);
                p["nonleptonic::Im{BT6}@SU3F"]  = -0.8 * std::sin(-0.8);
                p["nonleptonic::Re{BT15}@SU3F"] =  0.9 * std::cos( 0.9);
                p["nonleptonic::Im{BT15}@SU3F"] =  0.9 * std::sin( 0.9);
                p["nonleptonic::Re{DT3}@SU3F"]  = -1.0 * std::cos(-1.0);
                p["nonleptonic::Im{DT3}@SU3F"]  = -1.0 * std::sin(-1.0);
                p["nonleptonic::Re{AP3}@SU3F"]  =  1.1 * std::cos( 1.1);
                p["nonleptonic::Im{AP3}@SU3F"]  =  1.1 * std::sin( 1.1);
                p["nonleptonic::Re{CP3}@SU3F"]  = -1.2 * std::cos(-1.2);
                p["nonleptonic::Im{CP3}@SU3F"]  = -1.2 * std::sin(-1.2);
                p["nonleptonic::Re{AP6}@SU3F"]  =  1.3 * std::cos( 1.3);
                p["nonleptonic::Im{AP6}@SU3F"]  =  1.3 * std::sin( 1.3);
                p["nonleptonic::Re{CP6}@SU3F"]  = -1.4 * std::cos(-1.4);
                p["nonleptonic::Im{CP6}@SU3F"]  = -1.4 * std::sin(-1.4);
                p["nonleptonic::Re{AP15}@SU3F"] =  1.5 * std::cos( 1.5);
                p["nonleptonic::Im{AP15}@SU3F"] =  1.5 * std::sin( 1.5);
                p["nonleptonic::Re{CP15}@SU3F"] = -1.6 * std::cos(-1.6);
                p["nonleptonic::Im{CP15}@SU3F"] = -1.6 * std::sin(-1.6);
                p["nonleptonic::Re{BP3}@SU3F"]  =  1.7 * std::cos( 1.7);
                p["nonleptonic::Im{BP3}@SU3F"]  =  1.7 * std::sin( 1.7);
                p["nonleptonic::Re{BP6}@SU3F"]  = -1.8 * std::cos(-1.8);
                p["nonleptonic::Im{BP6}@SU3F"]  = -1.8 * std::sin(-1.8);
                p["nonleptonic::Re{BP15}@SU3F"] =  1.9 * std::cos( 1.9);
                p["nonleptonic::Im{BP15}@SU3F"] =  1.9 * std::sin( 1.9);
                p["nonleptonic::Re{DP3}@SU3F"]  = -2.0 * std::cos(-2.0);
                p["nonleptonic::Im{DP3}@SU3F"]  = -2.0 * std::sin(-2.0);
                p["eta::theta_18"]              =  1.0;

                static const double eps = 1.0e-6;


                Options o
                {
                    { "representation", "SU3F" },
                    { "q", "d" },
                    { "P1", "pi^+" },
                    { "P2", "pi^-" },
                    { "model", "CKM" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar d(p, o);


                Options oo
                {
                    { "representation", "SU3F" },
                    { "q", "d" },
                    { "P1", "pi^+" },
                    { "P2", "pi^-" },
                    { "model", "CKM" },
                    { "cp-conjugate", "false"}
                 };

                BToPseudoscalarPseudoscalar dd(p, oo);

                TEST_CHECK_NEARLY_EQUAL(0.5 * (d.branching_ratio() + dd.branching_ratio()),0.00420439042550259, eps);
                TEST_CHECK_NEARLY_EQUAL((d.branching_ratio() - dd.branching_ratio()) / (d.branching_ratio() + dd.branching_ratio()), 0.14470492226141188, eps);

                Options o3
                {
                    { "representation", "SU3F" },
                    { "q", "s" },
                    { "P1", "pi^0" },
                    { "P2", "Kbar_d" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar d3(p, o3);

                Options o4
                {
                    { "representation", "SU3F" },
                     { "q", "s" },
                    { "P1", "pi^0" },
                    { "P2", "Kbar_d" },
                    { "cp-conjugate", "false"}
                };

                BToPseudoscalarPseudoscalar d4(p, o4);

                TEST_CHECK_NEARLY_EQUAL(0.5 * (d3.branching_ratio() + d4.branching_ratio()), 0.0015865861860772175, eps);
                TEST_CHECK_NEARLY_EQUAL((d3.branching_ratio() - d4.branching_ratio()) / (d3.branching_ratio() + d4.branching_ratio()), 0.13929047286150023, eps);


                Options o5
                {
                    { "representation", "SU3F" },
                    { "q", "u" },
                    { "P1", "eta" },
                    { "P2", "pi^+" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar d5(p, o5);

                Options o6
                {
                    { "representation", "SU3F" },
                    { "q", "u" },
                    { "P1", "eta" },
                    { "P2", "pi^+" },
                    { "cp-conjugate", "false"}
                };

                BToPseudoscalarPseudoscalar d6(p, o6);

                TEST_CHECK_NEARLY_EQUAL(0.5 * (d5.branching_ratio() + d6.branching_ratio()), 0.016622890508282546, eps);
                TEST_CHECK_NEARLY_EQUAL((d5.branching_ratio() - d6.branching_ratio()) / (d5.branching_ratio() + d6.branching_ratio()), -0.06464070737112038, eps);

                Options o7
                {
                    { "representation", "SU3F" },
                    { "q", "d" },
                    { "P1", "eta'" },
                    { "P2", "K_d" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar d7(p, o7);

                Options o8
                {
                    { "representation", "SU3F" },
                    { "q", "d" },
                    { "P1", "eta'" },
                    { "P2", "K_d" },
                    { "cp-conjugate", "false"}
                };

                BToPseudoscalarPseudoscalar d8(p, o8);

                TEST_CHECK_NEARLY_EQUAL(0.5 * (d7.branching_ratio() + d8.branching_ratio()), 0.0013884379780450893, eps);
                TEST_CHECK_NEARLY_EQUAL((d7.branching_ratio() - d8.branching_ratio()) / (d7.branching_ratio() + d8.branching_ratio()), 0.025695932529644958, eps);

                Options o9
                {
                    { "representation", "SU3F" },
                    { "q", "d" },
                    { "P1", "eta" },
                    { "P2", "eta" },
                    { "cp-conjugate", "true"}
                };

                BToPseudoscalarPseudoscalar d9(p, o9);

                Options o10
                {
                    { "representation", "SU3F" },
                    { "q", "d" },
                    { "P1", "eta" },
                    { "P2", "eta" },
                    { "cp-conjugate", "false"}
                };

                BToPseudoscalarPseudoscalar d10(p, o10);

                TEST_CHECK_NEARLY_EQUAL(0.5 * (d9.branching_ratio() + d10.branching_ratio()), 0.0025284115535334, eps);
                TEST_CHECK_NEARLY_EQUAL((d9.branching_ratio() - d10.branching_ratio()) / (d9.branching_ratio() + d10.branching_ratio()), 0.06574658066886276, eps);
            }
        }
} b_to_psd_psd_test;
