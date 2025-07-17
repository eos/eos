/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
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
                p["life_time::Delta_B_s"]       =  0.089050132e+12;


                static const double eps = 1.0e-4;


                Options o
                {
                    { "representation"_ok, "SU3F" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "eta" },
                    { "P2"_ok, "K_d" },
                    { "model"_ok, "CKM" },
                };

                BToPseudoscalarPseudoscalar d(p, o);

                TEST_CHECK_RELATIVE_ERROR(d.avg_branching_ratio(), 1.4164770202846702e-6, eps);
                TEST_CHECK_RELATIVE_ERROR(d.cp_asymmetry(), 0.0008872201130703394, eps);

                Options o2
                {
                    { "representation"_ok, "SU3F" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "pi^+" },
                    { "P2"_ok, "pi^-" },
                    { "model"_ok, "CKM" },
                };

                BToPseudoscalarPseudoscalar d2(p, o2);

                TEST_CHECK_RELATIVE_ERROR(d2.avg_branching_ratio(), 6.382209068937381e-6, eps);
                TEST_CHECK_RELATIVE_ERROR(d2.cp_asymmetry(), -0.35001688891865884, eps);

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
                p["life_time::Delta_B_s"]       =  0.089050132e+12;

                static const double eps = 1.0e-6;

                Options o
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "pi^+" },
                    { "P2"_ok, "pi^-" },
                };

                BToPseudoscalarPseudoscalar d(p, o);

                TEST_CHECK_NEARLY_EQUAL( d.avg_branching_ratio(), 0.0022804835159666867, eps);
                TEST_CHECK_NEARLY_EQUAL(d.cp_asymmetry(), -0.09133792735291475, eps);

                Options oo
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "pi^0" },
                    { "P2"_ok, "Kbar_d" },
                };

                BToPseudoscalarPseudoscalar dd(p, oo);

                TEST_CHECK_NEARLY_EQUAL( dd.avg_branching_ratio(), 0.00001936527065018223, eps);
                TEST_CHECK_NEARLY_EQUAL( dd.cp_asymmetry(), -0.7366560008503581, eps);

                Options o3
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "u" },
                    { "P1"_ok, "eta" },
                    { "P2"_ok, "pi^+" },
                };

                BToPseudoscalarPseudoscalar d3(p, o3);


                TEST_CHECK_NEARLY_EQUAL(d3.avg_branching_ratio(), 0.0022644592272822654, eps);
                TEST_CHECK_NEARLY_EQUAL(d3.cp_asymmetry(), 0.16086877559501292, eps);

                Options o4
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "eta_prime" },
                    { "P2"_ok, "K_d" },
                };

                BToPseudoscalarPseudoscalar d4(p, o4);

                TEST_CHECK_NEARLY_EQUAL(d4.avg_branching_ratio(), 0.00537825569161075, eps);
                TEST_CHECK_NEARLY_EQUAL(d4.cp_asymmetry(), -0.006891941917256833, eps);


                Options o5
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "eta" },
                    { "P2"_ok, "eta" },
                };

                BToPseudoscalarPseudoscalar d5(p, o5);

                TEST_CHECK_NEARLY_EQUAL(d5.avg_branching_ratio(), 0.004473515098754955, eps);
                TEST_CHECK_NEARLY_EQUAL(d5.cp_asymmetry(), -0.14535763351380676, eps);

                Options o6
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "u" },
                    { "P1"_ok, "K_u" },
                    { "P2"_ok, "Kbar_d" },
                };

                BToPseudoscalarPseudoscalar d6(p, o6);

                TEST_CHECK_NEARLY_EQUAL(d6.avg_branching_ratio(), 0.00032802869044737436, eps);

                Options o7
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "u" },
                    { "P1"_ok, "pi^0" },
                    { "P2"_ok, "pi^+" },
                };

                BToPseudoscalarPseudoscalar d7(p, o7);

                TEST_CHECK_RELATIVE_ERROR(d7.exp_branching_ratio(), 0.0001673012642752443, eps);
                TEST_CHECK_RELATIVE_ERROR(d7.cp_asymmetry(), -0.024799122608299445, eps);

                Options o8
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "K_d" },
                    { "P2"_ok, "pi^0" },
                };

                BToPseudoscalarPseudoscalar d8(p, o8);

                TEST_CHECK_RELATIVE_ERROR(d8.exp_branching_ratio(), 0.00034695779913672127, eps);
                TEST_CHECK_RELATIVE_ERROR(d8.cp_asymmetry(), 0.04186432399144578, eps);

                Options o9
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "K_S" },
                    { "P2"_ok, "pi^0" },
                };

                BToPseudoscalarPseudoscalar d9(p, o9);

                TEST_CHECK_RELATIVE_ERROR(d9.mixing_induced_cp_asymmetry(), -0.7079345748007424, eps);
                TEST_CHECK_RELATIVE_ERROR(d9.a_Delta_Gamma(), 0.7050361807584282, eps);


                Options o10
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "Kbar_d" },
                    { "P2"_ok, "pi^0" },
                };

                BToPseudoscalarPseudoscalar d10(p, o10);

                TEST_CHECK_RELATIVE_ERROR(d10.exp_branching_ratio(), 0.000019453554976036367, eps);
                TEST_CHECK_RELATIVE_ERROR(d10.cp_asymmetry(), -0.7366560008503581, eps);

                Options o11
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "K_S" },
                    { "P2"_ok, "pi^0" },
                };

                BToPseudoscalarPseudoscalar d11(p, o11);

                TEST_CHECK_RELATIVE_ERROR(d11.cp_asymmetry(), -0.7366560008503581, eps);
                TEST_CHECK_RELATIVE_ERROR(d11.mixing_induced_cp_asymmetry(), 0.26599025244408186, eps);
                TEST_CHECK_RELATIVE_ERROR(d11.a_Delta_Gamma(), 0.6217613063032237, eps);

                Options o12
                {
                    { "representation"_ok, "topological" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "Kbar_u" },
                    { "P2"_ok, "pi^+" },
                };

                BToPseudoscalarPseudoscalar d12(p, o12);

                TEST_CHECK_RELATIVE_ERROR(d12.exp_branching_ratio(), 0.00020527347065702994, eps);
                TEST_CHECK_RELATIVE_ERROR(d12.cp_asymmetry(), 0.23578152390612866, eps);
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
                p["life_time::Delta_B_s"]       =  0.089050132e+12;

                static const double eps = 1.0e-6;


                Options o
                {
                    { "representation"_ok, "SU3F" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "pi^+" },
                    { "P2"_ok, "pi^-" },
                    { "model"_ok, "CKM" },
                };

                BToPseudoscalarPseudoscalar d(p, o);

                TEST_CHECK_NEARLY_EQUAL(d.avg_branching_ratio(),0.00420439042550259, eps);
                TEST_CHECK_NEARLY_EQUAL(d.cp_asymmetry(), -0.14470492226141188, eps);

                Options o2
                {
                    { "representation"_ok, "SU3F" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "pi^0" },
                    { "P2"_ok, "Kbar_d" },
                };

                BToPseudoscalarPseudoscalar d2(p, o2);

                TEST_CHECK_NEARLY_EQUAL(d2.avg_branching_ratio(), 0.0015865861860772175, eps);
                TEST_CHECK_NEARLY_EQUAL(d2.cp_asymmetry(), -0.13929047286150023, eps);


                Options o3
                {
                    { "representation"_ok, "SU3F" },
                    { "q"_ok, "u" },
                    { "P1"_ok, "eta" },
                    { "P2"_ok, "pi^+" },
                };

                BToPseudoscalarPseudoscalar d3(p, o3);

                TEST_CHECK_NEARLY_EQUAL(d3.avg_branching_ratio(), 0.016622890508282546, eps);
                TEST_CHECK_NEARLY_EQUAL(d3.cp_asymmetry(), 0.06464070737112038, eps);

                Options o4
                {
                    { "representation"_ok, "SU3F" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "eta_prime" },
                    { "P2"_ok, "K_d" },
                };

                BToPseudoscalarPseudoscalar d4(p, o4);

                TEST_CHECK_NEARLY_EQUAL(d4.avg_branching_ratio(), 0.0013884379780450893, eps);
                TEST_CHECK_NEARLY_EQUAL(d4.cp_asymmetry(), -0.025695932529644958, eps);

                Options o5
                {
                    { "representation"_ok, "SU3F" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "eta" },
                    { "P2"_ok, "eta" },
                };

                BToPseudoscalarPseudoscalar d5(p, o5);

                TEST_CHECK_NEARLY_EQUAL(d5.avg_branching_ratio(), 0.0025284115535334, eps);
                TEST_CHECK_NEARLY_EQUAL(d5.cp_asymmetry(), -0.06574658066886276, eps);
            }

            /*Test with QCDF amplitudes*/

            {
                Parameters p = Parameters::Defaults();

                p["nonleptonic::Re{alpha1}@QCDF"]     =  0.1 * std::cos( 0.1);
                p["nonleptonic::Im{alpha1}@QCDF"]     =  0.1 * std::sin( 0.1);
                p["nonleptonic::Re{alpha2}@QCDF"]     = -0.2 * std::cos(-0.2);
                p["nonleptonic::Im{alpha2}@QCDF"]     = -0.2 * std::sin(-0.2);
                p["nonleptonic::Re{b2}@QCDF"]         =  0.3 * std::cos( 0.3);
                p["nonleptonic::Im{b2}@QCDF"]         =  0.3 * std::sin( 0.3);
                p["nonleptonic::Re{b1}@QCDF"]         = -0.4 * std::cos(-0.4);
                p["nonleptonic::Im{b1}@QCDF"]         = -0.4 * std::sin(-0.4);
                p["nonleptonic::Re{bS2}@QCDF"]        =  0.5 * std::cos( 0.5);
                p["nonleptonic::Im{bS2}@QCDF"]        =  0.5 * std::sin( 0.5);
                p["nonleptonic::Re{bS1}@QCDF"]        = -0.6 * std::cos(-0.6);
                p["nonleptonic::Im{bS1}@QCDF"]        = -0.6 * std::sin(-0.6);
                p["nonleptonic::Re{alpha4_u}@QCDF"]   =  0.7 * std::cos( 0.7);
                p["nonleptonic::Im{alpha4_u}@QCDF"]   =  0.7 * std::sin( 0.7);
                p["nonleptonic::Re{alpha3_u}@QCDF"]   = -0.8 * std::cos(-0.8);
                p["nonleptonic::Im{alpha3_u}@QCDF"]   = -0.8 * std::sin(-0.8);
                p["nonleptonic::Re{b4_u}@QCDF"]       =  0.9 * std::cos( 0.9);
                p["nonleptonic::Im{b4_u}@QCDF"]       =  0.9 * std::sin( 0.9);
                p["nonleptonic::Re{bS4_u}@QCDF"]      = -1.0 * std::cos(-1.0);
                p["nonleptonic::Im{bS4_u}@QCDF"]      = -1.0 * std::sin(-1.0);
                p["nonleptonic::Re{alpha4EW_c}@QCDF"] =  1.1 * std::cos( 1.1);
                p["nonleptonic::Im{alpha4EW_c}@QCDF"] =  1.1 * std::sin( 1.1);
                p["nonleptonic::Re{alpha3EW_c}@QCDF"] = -1.2 * std::cos(-1.2);
                p["nonleptonic::Im{alpha3EW_c}@QCDF"] = -1.2 * std::sin(-1.2);
                p["nonleptonic::Re{b3EW_c}@QCDF"]     =  1.3 * std::cos( 1.3);
                p["nonleptonic::Im{b3EW_c}@QCDF"]     =  1.3 * std::sin( 1.3);
                p["nonleptonic::Re{b4EW_c}@QCDF"]     = -1.4 * std::cos(-1.4);
                p["nonleptonic::Im{b4EW_c}@QCDF"]     = -1.4 * std::sin(-1.4);
                p["nonleptonic::Re{bS3EW_c}@QCDF"]    =  1.5 * std::cos( 1.5);
                p["nonleptonic::Im{bS3EW_c}@QCDF"]    =  1.5 * std::sin( 1.5);
                p["nonleptonic::Re{bS4EW_c}@QCDF"]    = -1.6 * std::cos(-1.6);
                p["nonleptonic::Im{bS4EW_c}@QCDF"]    = -1.6 * std::sin(-1.6);
                p["nonleptonic::Re{alpha4_c}@QCDF"]   =  1.7 * std::cos( 1.7);
                p["nonleptonic::Im{alpha4_c}@QCDF"]   =  1.7 * std::sin( 1.7);
                p["nonleptonic::Re{alpha3_c}@QCDF"]   = -1.8 * std::cos(-1.8);
                p["nonleptonic::Im{alpha3_c}@QCDF"]   = -1.8 * std::sin(-1.8);
                p["nonleptonic::Re{b4_c}@QCDF"]       =  1.9 * std::cos( 1.9);
                p["nonleptonic::Im{b4_c}@QCDF"]       =  1.9 * std::sin( 1.9);
                p["nonleptonic::Re{bS4_c}@QCDF"]      = -2.0 * std::cos(-2.0);
                p["nonleptonic::Im{bS4_c}@QCDF"]      = -2.0 * std::sin(-2.0);
                p["eta::theta_18"]                    =  0.0;
                p["eta::theta_FKS"]                   =  0.5;
                p["life_time::Delta_B_s"]       =  0.089050132e+12;

                static const double eps = 1.0e-6;


                Options o
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "pi^+" },
                    { "P2"_ok, "pi^-" },
                };

                BToPseudoscalarPseudoscalar d(p, o);

                TEST_CHECK_RELATIVE_ERROR(d.avg_branching_ratio(),0.0006535297808412619, eps);
                TEST_CHECK_RELATIVE_ERROR(d.cp_asymmetry(), -0.12500385563385133, eps);

                Options o2
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "pi^0" },
                    { "P2"_ok, "Kbar_d" },
                };

                BToPseudoscalarPseudoscalar d2(p, o2);

                TEST_CHECK_RELATIVE_ERROR(d2.avg_branching_ratio(), 6.578184442118153e-6, eps);
                TEST_CHECK_RELATIVE_ERROR(d2.cp_asymmetry(), 0.5974882599801686, eps);

                Options o7
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "u" },
                    { "P1"_ok, "pi^0" },
                    { "P2"_ok, "pi^+" },
                };

                BToPseudoscalarPseudoscalar d7(p, o7);

                TEST_CHECK_RELATIVE_ERROR(d7.exp_branching_ratio(), 0.00035056299610878797, eps);
                TEST_CHECK_RELATIVE_ERROR(d7.cp_asymmetry(), 0.021852766998392625, eps);

                Options o8
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "K_d" },
                    { "P2"_ok, "pi^0" },
                };

                BToPseudoscalarPseudoscalar d8(p, o8);

                TEST_CHECK_RELATIVE_ERROR(d8.exp_branching_ratio(), 0.00015267686389138193, eps);
                TEST_CHECK_RELATIVE_ERROR(d8.cp_asymmetry(), -0.00454281035369326, eps);

                Options o9
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "K_S" },
                    { "P2"_ok, "pi^0" },
                };

                BToPseudoscalarPseudoscalar d9(p, o9);

                TEST_CHECK_RELATIVE_ERROR(d9.exp_branching_ratio(), 0.00015267686389138193/2, eps);
                TEST_CHECK_RELATIVE_ERROR(d9.cp_asymmetry(), -0.00454281035369326, eps);
                TEST_CHECK_RELATIVE_ERROR(d9.mixing_induced_cp_asymmetry(), -0.7443181962311635, eps);
                TEST_CHECK_RELATIVE_ERROR(d9.a_Delta_Gamma(), 0.6678096926769465, eps);

                Options o10
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "Kbar_d" },
                    { "P2"_ok, "pi^0" },
                };

                BToPseudoscalarPseudoscalar d10(p, o10);

                TEST_CHECK_RELATIVE_ERROR(d10.exp_branching_ratio(), 6.608293479283136e-6, eps);
                TEST_CHECK_RELATIVE_ERROR(d10.cp_asymmetry(), 0.5974882599801683, eps);

                Options o11
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "K_S" },
                    { "P2"_ok, "pi^0" },
                };

                BToPseudoscalarPseudoscalar d11(p, o11);

                TEST_CHECK_RELATIVE_ERROR(d11.exp_branching_ratio(), 6.660246957763282e-6 / 2, eps);
                TEST_CHECK_RELATIVE_ERROR(d11.cp_asymmetry(), 0.5974882599801683, eps);
                TEST_CHECK_RELATIVE_ERROR(d11.mixing_induced_cp_asymmetry(), 0.7933738414149348, eps);
                TEST_CHECK_RELATIVE_ERROR(d11.a_Delta_Gamma(), 0.11647200068849438, eps);

                Options o12
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "Kbar_u" },
                    { "P2"_ok, "pi^+" },
                };

                BToPseudoscalarPseudoscalar d12(p, o12);

                TEST_CHECK_RELATIVE_ERROR(d12.exp_branching_ratio(), 0.0006458820066104115, eps);
                TEST_CHECK_RELATIVE_ERROR(d12.cp_asymmetry(), -0.12499343458253753, eps);
                /*TEST_CHECK_RELATIVE_ERROR(d12.mixing_induced_cp_asymmetry(), 0.0, eps);
                TEST_CHECK_RELATIVE_ERROR(d12.a_Delta_Gamma(), 0.0, eps);*/

                Options o13
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "pi^+" },
                    { "P2"_ok, "pi^-" },
                };

                BToPseudoscalarPseudoscalar d13(p, o13);

                TEST_CHECK_RELATIVE_ERROR(d13.exp_branching_ratio(), 4.821560613831973e-7, eps);
                TEST_CHECK_RELATIVE_ERROR(d13.cp_asymmetry(), 0.0072671672100170165, eps);
                TEST_CHECK_RELATIVE_ERROR(d13.mixing_induced_cp_asymmetry(), -0.028884541837291628, eps);
                TEST_CHECK_RELATIVE_ERROR(d13.a_Delta_Gamma(), -0.9995563373435195, eps);


                Options o14
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "d" },
                    { "P1"_ok, "pi^0" },
                    { "P2"_ok, "eta" },
                };

                BToPseudoscalarPseudoscalar d14(p, o14);

                TEST_CHECK_RELATIVE_ERROR(d14.exp_branching_ratio(), 0.00036489812391703603, eps);
                TEST_CHECK_RELATIVE_ERROR(d14.cp_asymmetry(), -0.028354067521044263, eps);
                TEST_CHECK_RELATIVE_ERROR(d14.mixing_induced_cp_asymmetry(), 0.5539879807384019, eps);
                TEST_CHECK_RELATIVE_ERROR(d14.a_Delta_Gamma(), -0.8320416840834355, eps);

                Options o15
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "u" },
                    { "P1"_ok, "K_u" },
                    { "P2"_ok, "eta_prime" },
                };

                BToPseudoscalarPseudoscalar d15(p, o15);

                TEST_CHECK_RELATIVE_ERROR(d15.exp_branching_ratio(), 0.0662594630144483, eps);
                TEST_CHECK_RELATIVE_ERROR(d15.cp_asymmetry(), -0.0026511491186821377, eps);

                Options o16
                {
                    { "representation"_ok, "QCDF" },
                    { "q"_ok, "s" },
                    { "P1"_ok, "eta_prime" },
                    { "P2"_ok, "eta_prime" },
                };

                BToPseudoscalarPseudoscalar d16(p, o16);

                TEST_CHECK_RELATIVE_ERROR(d16.exp_branching_ratio(), 0.03037051193347059, eps);
                TEST_CHECK_RELATIVE_ERROR(d16.cp_asymmetry(), -0.0036695733788966248, eps);
                TEST_CHECK_RELATIVE_ERROR(d16.mixing_induced_cp_asymmetry(), -0.0270033482282267, eps);
                TEST_CHECK_RELATIVE_ERROR(d16.a_Delta_Gamma(), -0.9996286077417366, eps);

            }

        }
} b_to_psd_psd_test;
