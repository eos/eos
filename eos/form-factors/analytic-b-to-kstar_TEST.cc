/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
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
#include <eos/form-factors/analytic-b-to-kstar.hh>
#include <eos/form-factors/mesonic.hh>

#include <vector>
#include <utility>

using namespace test;
using namespace eos;

class KMO2006FormFactorsTest :
    public TestCase
{
    public:
        KMO2006FormFactorsTest() :
            TestCase("kmo2006_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> K^* diagnostic values */
            {
                Parameters p = Parameters::Defaults();
                p["mu"]                  = 4.2;
                p["lambda_B_p"]          = 0.460;
                p["mass::b(MSbar)"]      = 4.17;
                p["mass::B_d"]           = 5.2795;
                p["mass::K^*_d"]         = 0.896;
                p["decay-constant::B_d"] = 0.180;
                p["K^*::f_para"]         = 0.217;
                p["B->K^*::s_0@LCSR"]    = 1.7;
                p["B->K^*::M^2@LCSR"]    = 1.0;

                AnalyticFormFactorBToKstarKMO2006 ff{ p, Options{} };
                auto diagnostics = ff.diagnostics();

                //std::cout << "Diagnostics:" << std::endl;
                //for (auto & d : diagnostics)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                //std::cout << "Diagnostics ended" << std::endl;

                // references taken from Javier's Mathematica notebook
                static const std::vector<std::pair<double, double>> reference
                {
                    // V
                    // skipping Iota_1 = 0
                    std::make_pair(-0.0139373,  1.0e-7), // Iota_2V(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0139373,  1.0e-7), // Iota_2V(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0139373,  1.0e-7), // Iota_2V(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0174407,  1.0e-7), // Iota_2V(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0174407,  1.0e-7), // Iota_2V(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0174407,  1.0e-7), // Iota_2V(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.00719578, 1.0e-7), // Iota_3V(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.00625785, 1.0e-7), // Iota_3V(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.00602337, 1.0e-7), // Iota_3V(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0122255,  1.0e-7), // Iota_3V(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.010575,   1.0e-7), // Iota_3V(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0101624,  1.0e-7), // Iota_3V(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.25717,    1.0e-6), // DIota_3V(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.222014,   1.0e-6), // DIota_3V(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.213224,   1.0e-6), // DIota_3V(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.238114,   1.0e-6), // DIota_3V(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.203049,   1.0e-6), // DIota_3V(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.194282,   1.0e-6), // DIota_3V(q^2 = -0, sigma = 0.06)

                    // A_0
                    // skipping Iota_1 = 0
                    std::make_pair(-0.0104636,  1.0e-7), // Iota_2A0(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0104636,  1.0e-7), // Iota_2A0(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0104636,  1.0e-7), // Iota_2A0(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0190116,  1.0e-7), // Iota_2A0(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0190116,  1.0e-7), // Iota_2A0(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0190116,  1.0e-7), // Iota_2A0(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.00436787, 1.0e-7), // Iota_3A0(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.00379854, 1.0e-7), // Iota_3A0(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.00365621, 1.0e-7), // Iota_3A0(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0100129,  1.0e-7), // Iota_3A0(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.00866109, 1.0e-7), // Iota_3A0(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.00832314, 1.0e-7), // Iota_3A0(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.241322,   1.0e-6), // DIota_3A0(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.208874,   1.0e-6), // DIota_3A0(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.200762,   1.0e-6), // DIota_3A0(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.309905,   1.0e-6), // DIota_3A0(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.265675,   1.0e-6), // DIota_3A0(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.254618,   1.0e-6), // DIota_3A0(q^2 = -0, sigma = 0.06)

                    // A_1
                    std::make_pair(-0.0026399,  1.0e-7), // Iota_1A1(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0026399,  1.0e-7), // Iota_1A1(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0026399,  1.0e-7), // Iota_1A1(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.00330348, 1.0e-7), // Iota_1A1(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.00330348, 1.0e-7), // Iota_1A1(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.00330348, 1.0e-7), // Iota_1A1(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.0801903,  1.0e-7), // Iota_2A1(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0690130,  1.0e-7), // Iota_2A1(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0662187,  1.0e-7), // Iota_2A1(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0976144,  1.0e-7), // Iota_2A1(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0832444,  1.0e-7), // Iota_2A1(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0796519,  1.0e-7), // Iota_2A1(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.0435693,  1.0e-7), // Iota_3A1(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0329515,  1.0e-7), // Iota_3A1(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0305284,  1.0e-7), // Iota_3A1(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.072989,   1.0e-7), // Iota_3A1(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0546117,  1.0e-7), // Iota_3A1(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.050433,   1.0e-7), // Iota_3A1(q^2 = -0, sigma = 0.06)

                    std::make_pair(-1.52653,    3.0e-6), // DIota_3A1(q^2 = -5, sigma = 0.04)
                    std::make_pair(-1.13729,    3.0e-6), // DIota_3A1(q^2 = -1, sigma = 0.04)
                    std::make_pair(-1.04889,    3.0e-6), // DIota_3A1(q^2 = -0, sigma = 0.04)

                    std::make_pair(-1.37016,    3.0e-6), // DIota_3A1(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.995026,   3.0e-6), // DIota_3A1(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.910518,   3.0e-6), // DIota_3A1(q^2 = -0, sigma = 0.06)

                    // A_2
                    // skipping Iota_1 = 0
                    std::make_pair(-0.0781378,  1.0e-7), // Iota_2A2(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0781378,  1.0e-7), // Iota_2A2(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0781378,  1.0e-7), // Iota_2A2(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.100711,   1.0e-6), // Iota_2A2(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.100711,   1.0e-6), // Iota_2A2(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.100711,   1.0e-6), // Iota_2A2(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.0448074,  1.0e-7), // Iota_3A2(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0402517,  1.0e-7), // Iota_3A2(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0391128,  1.0e-7), // Iota_3A2(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0798836,  1.0e-6), // Iota_3A2(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0722155,  1.0e-6), // Iota_3A2(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0702985,  1.0e-6), // Iota_3A2(q^2 = -0, sigma = 0.06)

                    std::make_pair(-1.72831,    5.0e-6), // DIota_3A2(q^2 = -5, sigma = 0.04)
                    std::make_pair(-1.56746,    5.0e-6), // DIota_3A2(q^2 = -1, sigma = 0.04)
                    std::make_pair(-1.52724,    5.0e-6), // DIota_3A2(q^2 = -0, sigma = 0.04)

                    std::make_pair(-1.71681,    5.0e-6), // DIota_3A2(q^2 = -5, sigma = 0.06)
                    std::make_pair(-1.57132,    5.0e-6), // DIota_3A2(q^2 = -1, sigma = 0.06)
                    std::make_pair(-1.53495,    5.0e-6), // DIota_3A2(q^2 = -0, sigma = 0.06)

                    // T_1
                    // skipping Iota_1 = 0
                    std::make_pair(-0.0678584,  1.0e-7), // Iota_2T1(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0678584,  1.0e-7), // Iota_2T1(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0678584,  1.0e-7), // Iota_2T1(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0821403,  1.0e-7), // Iota_2T1(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0821403,  1.0e-7), // Iota_2T1(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0821403,  1.0e-7), // Iota_2T1(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.0364705,  1.0e-7), // Iota_3T1(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0317168,  1.0e-7), // Iota_3T1(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0305284,  1.0e-7), // Iota_3T1(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0606717,  1.0e-7), // Iota_3T1(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0524808,  1.0e-7), // Iota_3T1(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0504330,  1.0e-7), // Iota_3T1(q^2 = -0, sigma = 0.06)

                    std::make_pair(-1.26543,    5.0e-6), // DIota_3T1(q^2 = -5, sigma = 0.04)
                    std::make_pair(-1.0922,     5.0e-6), // DIota_3T1(q^2 = -1, sigma = 0.04)
                    std::make_pair(-1.04889,    5.0e-6), // DIota_3T1(q^2 = -0, sigma = 0.04)

                    std::make_pair(-1.11715,    5.0e-6), // DIota_3T1(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.951845,   5.0e-6), // DIota_3T1(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.910518,   5.0e-6), // DIota_3T1(q^2 = -0, sigma = 0.06)

                    // T_23A
                    // skipping Iota_1 = 0
                    std::make_pair(-0.0678584,  1.0e-7), // Iota_2T23A(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0678584,  1.0e-7), // Iota_2T23A(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0678584,  1.0e-7), // Iota_2T23A(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0821403,  1.0e-7), // Iota_2T23A(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0821403,  1.0e-7), // Iota_2T23A(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0821403,  1.0e-7), // Iota_2T23A(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.0346288,  1.0e-7), // Iota_3T23A(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0313484,  1.0e-7), // Iota_3T23A(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0305284,  1.0e-7), // Iota_3T23A(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0562664,  1.0e-7), // Iota_3T23A(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0515997,  1.0e-7), // Iota_3T23A(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0504330,  1.0e-7), // Iota_3T23A(q^2 = -0, sigma = 0.06)

                    std::make_pair(-1.15977,    5.0e-6), // DIota_3T23A(q^2 = -5, sigma = 0.04)
                    std::make_pair(-1.07107,    5.0e-6), // DIota_3T23A(q^2 = -1, sigma = 0.04)
                    std::make_pair(-1.04889,    5.0e-6), // DIota_3T23A(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.971319,   5.0e-6), // DIota_3T23A(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.922678,   5.0e-6), // DIota_3T23A(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.910518,   5.0e-6), // DIota_3T23A(q^2 = -0, sigma = 0.06)

                    std::make_pair(+0.20191,    1.0e-3), // T_23A(q^2 = -5)
                    std::make_pair(+0.264059,   1.0e-3), // T_23A(q^2 = -1)
                    std::make_pair(+0.284136,   1.0e-3), // T_23A(q^2 = -0)

                    // T_23B
                    // skipping Iota_1 = 0
                    std::make_pair(+0.00614296, 1.0e-7), // Iota_2T23B(q^2 = -5, sigma = 0.04)
                    std::make_pair(+0.00614296, 1.0e-7), // Iota_2T23B(q^2 = -1, sigma = 0.04)
                    std::make_pair(+0.00614296, 1.0e-7), // Iota_2T23B(q^2 = -0, sigma = 0.04)

                    std::make_pair(+0.0112214,  1.0e-7), // Iota_2T23B(q^2 = -5, sigma = 0.06)
                    std::make_pair(+0.0112214,  1.0e-7), // Iota_2T23B(q^2 = -1, sigma = 0.06)
                    std::make_pair(+0.0112214,  1.0e-7), // Iota_2T23B(q^2 = -0, sigma = 0.06)

                    std::make_pair(+0.00719938, 1.0e-7), // Iota_3T23B(q^2 = -5, sigma = 0.04)
                    std::make_pair(+0.0062953,  1.0e-7), // Iota_3T23B(q^2 = -1, sigma = 0.04)
                    std::make_pair(+0.0060928,  3.0e-5), // Iota_3T23B(q^2 = -0, sigma = 0.04)

                    std::make_pair(+0.0167722,  1.0e-7), // Iota_3T23B(q^2 = -5, sigma = 0.06)
                    std::make_pair(+0.0145997,  1.0e-7), // Iota_3T23B(q^2 = -1, sigma = 0.06)
                    std::make_pair(+0.0140566,  1.0e-7), // Iota_3T23B(q^2 = -0, sigma = 0.06)

                    std::make_pair(+0.405612,   5.0e-6), // DIota_3T23B(q^2 = -5, sigma = 0.04)
                    std::make_pair(+0.353532,   5.0e-6), // DIota_3T23B(q^2 = -1, sigma = 0.04)
                    std::make_pair(+0.340512,   5.0e-6), // DIota_3T23B(q^2 = -0, sigma = 0.04)

                    std::make_pair(+0.528061,   5.0e-6), // DIota_3T23B(q^2 = -5, sigma = 0.06)
                    std::make_pair(+0.455624,   5.0e-6), // DIota_3T23B(q^2 = -1, sigma = 0.06)
                    std::make_pair(+0.437514,   5.0e-6), // DIota_3T23B(q^2 = -0, sigma = 0.06)

                    std::make_pair(+0.0608932,  5.0e-4), // T_23B(q^2 = -5)
                    std::make_pair(+0.0740864,  5.0e-4), // T_23B(q^2 = -1)
                    std::make_pair(+0.0781295,  5.0e-4), // T_23B(q^2 = -0)

                    // T_2
                    // skipping Iota_1 = 0
                    std::make_pair(-0.0826571, 5.0e-6), // Iota_2T2(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0708181, 5.0e-6), // Iota_2T2(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0678584, 5.0e-6), // Iota_2T2(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.1014520, 6.0e-6), // Iota_2T2(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0860026, 6.0e-6), // Iota_2T2(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0821403, 6.0e-6), // Iota_2T2(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.0436818, 5.0e-6), // Iota_3T2(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0329711, 5.0e-6), // Iota_3T2(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0305284, 5.0e-6), // Iota_3T2(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0728501, 5.0e-6), // Iota_3T2(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0545836, 5.0e-6), // Iota_3T2(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0504330, 5.0e-6), // Iota_3T2(q^2 = -0, sigma = 0.06)

                    std::make_pair(-1.52372,   5.0e-4), // DIota_3T2(q^2 = -5, sigma = 0.04)
                    std::make_pair(-1.13673,   5.0e-4), // DIota_3T2(q^2 = -1, sigma = 0.04)
                    std::make_pair(-1.04889,   5.0e-4), // DIota_3T2(q^2 = -0, sigma = 0.04)

                    std::make_pair(-1.34569,   2.0e-4), // DIota_3T2(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.990406,  2.0e-4), // DIota_3T2(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.910518,  2.0e-4), // DIota_3T2(q^2 = -0, sigma = 0.06)

                    // T_3
                    // skipping Iota_1 = 0
                    std::make_pair(-0.0801443, 5.0e-6), // Iota_2T3(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0801443, 5.0e-6), // Iota_2T3(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0801443, 5.0e-6), // Iota_2T3(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.1045830, 6.0e-6), // Iota_2T3(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.1045830, 6.0e-6), // Iota_2T3(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.1045830, 6.0e-6), // Iota_2T3(q^2 = -0, sigma = 0.06)

                    std::make_pair(-0.0490275, 5.0e-6), // Iota_3T3(q^2 = -5, sigma = 0.04)
                    std::make_pair(-0.0439390, 5.0e-6), // Iota_3T3(q^2 = -1, sigma = 0.04)
                    std::make_pair(-0.0426669, 5.0e-6), // Iota_3T3(q^2 = -0, sigma = 0.04)

                    std::make_pair(-0.0898107, 5.0e-6), // Iota_3T3(q^2 = -5, sigma = 0.06)
                    std::make_pair(-0.0807991, 5.0e-6), // Iota_3T3(q^2 = -1, sigma = 0.06)
                    std::make_pair(-0.0785462, 5.0e-6), // Iota_3T3(q^2 = -0, sigma = 0.06)

                    std::make_pair(-1.97099,   5.0e-4), // DIota_3T3(q^2 = -5, sigma = 0.04)
                    std::make_pair(-1.77813,   5.0e-4), // DIota_3T3(q^2 = -1, sigma = 0.04)
                    std::make_pair(-1.72991,   5.0e-4), // DIota_3T3(q^2 = -0, sigma = 0.04)

                    std::make_pair(-2.02744,   2.0e-4), // DIota_3T3(q^2 = -5, sigma = 0.06)
                    std::make_pair(-1.83393,   2.0e-4), // DIota_3T3(q^2 = -1, sigma = 0.06)
                    std::make_pair(-1.78555,   2.0e-4), // DIota_3T3(q^2 = -0, sigma = 0.06)
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

            /* B -> K^* form factor values */
            {
                static const double eps = 3e-3; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*@KMO2006", p);

                p["mu"]                  = 4.2;
                p["lambda_B_p"]          = 0.460;
                p["mass::b(MSbar)"]      = 4.17;
                p["mass::B_d"]           = 5.2795;
                p["mass::K^*_d"]         = 0.896;
                p["decay-constant::B_d"] = 0.180;
                p["K^*::f_para"]         = 0.217;
                p["B->K^*::s_0@LCSR"]    = 1.7;
                p["B->K^*::M^2@LCSR"]    = 1.0;

                TEST_CHECK_RELATIVE_ERROR( 0.256730, ff->v(-5.0),          eps);
                TEST_CHECK_RELATIVE_ERROR( 0.321635, ff->v(-1.0),          eps);
                TEST_CHECK_RELATIVE_ERROR( 0.341850, ff->v( 0.0),          eps);

                TEST_CHECK_RELATIVE_ERROR( 0.185435, ff->a_0( -5.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.229611, ff->a_0( -1.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.243154, ff->a_0(  0.0),       eps);

                TEST_CHECK_RELATIVE_ERROR( 0.219723, ff->a_1( -5.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.240066, ff->a_1( -1.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.245726, ff->a_1(  0.0),       eps);

                TEST_CHECK_RELATIVE_ERROR( 0.168974, ff->a_2( -5.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.197871, ff->a_2( -1.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.205965, ff->a_2(  0.0),       eps);

                TEST_CHECK_RELATIVE_ERROR( 0.221752, ff->a_12(-5.0), 3.0 * eps);
                TEST_CHECK_RELATIVE_ERROR( 0.241467, ff->a_12(-1.0), 3.0 * eps);
                TEST_CHECK_RELATIVE_ERROR( 0.246990, ff->a_12( 0.0), 3.0 * eps);

                TEST_CHECK_RELATIVE_ERROR( 0.214249, ff->t_1( -5.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.267582, ff->t_1( -1.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.284136, ff->t_1(  0.0),       eps);

                TEST_CHECK_RELATIVE_ERROR( 0.227965, ff->t_2( -5.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.271092, ff->t_2( -1.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.284136, ff->t_2(  0.0),       eps);

                TEST_CHECK_RELATIVE_ERROR( 0.141060, ff->t_3( -5.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.189972, ff->t_3( -1.0),       eps);
                TEST_CHECK_RELATIVE_ERROR( 0.206006, ff->t_3(  0.0),       eps);

                TEST_CHECK_RELATIVE_ERROR( 0.476978, ff->t_23(-5.0), 2.0 * eps);
                TEST_CHECK_RELATIVE_ERROR( 0.536617, ff->t_23(-1.0), 2.0 * eps);
                TEST_CHECK_RELATIVE_ERROR( 0.554854, ff->t_23( 0.0), 2.0 * eps);
            }
        }
} kmo2006_form_factors_test;
