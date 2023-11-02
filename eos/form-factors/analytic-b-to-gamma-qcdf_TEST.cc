/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022-2024 Philip Lüghausen
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
#include <eos/form-factors/analytic-b-to-gamma-qcdf.hh>
#include <eos/observable.hh>
#include <eos/utils/diagnostics.hh>

using namespace test;
using namespace eos;

class AnalyticFormFactorBToGammaQCDFTest :
    public TestCase
{
    public:
        AnalyticFormFactorBToGammaQCDFTest():
            TestCase("analytic_b_to_gamma_qcdf_test")
        {
        }

        virtual void run() const
        {
            // Independent numerical reference
            {
                Parameters p = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"]  = 0.7;
                p["B_u::mu_0@FLvD2022"]  = 1.3;
                p["B->gamma::mu@FLvD2022QCDF"]  = 1.3;
                p["B->gamma::mu_h1@FLvD2022QCDF"]  = 4.7;
                p["B->gamma::mu_h2@FLvD2022QCDF"]  = 4.5;
                p["B->gamma::s_0@FLvD2022QCDF"] = 1.59;
                p["B->gamma::M^2@FLvD2022QCDF"] = 1.35;
                p["decay-constant::B_u"] = 0.129;
                p["mass::B_u"] = 5.27929;
                p["mass::b(MSbar)"] = 4.45432371854873; // fix m_b_pole@1-loop to 4.8
                p["mass::rho^+"] = 0.7;
                p["B::lambda_E^2"] = 0.0625;
                p["B::lambda_H^2"] = 0.125;
                p["B::LambdaBar"] = 0.5;
                p["B_u::a^phi+_0@FLvD2022"] =  1.868119356054707;
                p["B_u::a^phi+_1@FLvD2022"] =  0.151143197362311;
                p["B_u::a^phi+_2@FLvD2022"] =  1.203196552637887;
                p["B_u::a^phi+_3@FLvD2022"] =  0.429631987348729;
                p["B_u::a^phi+_4@FLvD2022"] =  0.304198191688109;
                p["B_u::a^phi+_5@FLvD2022"] = -0.324469147908141;
                p["B_u::a^phi+_6@FLvD2022"] =  0.381019563820993;
                p["B_u::a^phi+_7@FLvD2022"] = -0.246884872397705;
                p["B_u::a^phi+_8@FLvD2022"] = -0.058121797086248;


                // Diagnostics: check pieces against Mathematica implementation

                AnalyticFormFactorBToGammaQCDF ff(p, Options());

                static const std::vector<std::pair<double, double>> reference = {
                    std::make_pair(3.39713985820215,    1e-9), // L0()
                    std::make_pair(3.32067923218836,    1e-9), // L0_incomplete(8.0)
                    std::make_pair(2.9528207810186,     1e-9), // norm_incomplete(8.0)
                    std::make_pair(0.156908479594529,   1e-9), // lapltr_incomplete(8.0, 4.0)
                    std::make_pair(-0.186641425933295,  1e-9), // lapltr_incomplete_dsigma(8.0, 4.0)
                    std::make_pair(0.272354067021644,   1e-3), // L0_effective(3.0); numerical reference is imprecise
                    std::make_pair(0.0527171469619207,  1e-3), // L0_effective(2.16)
                    std::make_pair(10.4492075178413 + -10.4781709714967 + 6.58190087562423 + -8.92720937287174,     1e-3), // L0_incomplete_effective(3.0, 8.0); numerical reference is imprecise
                    std::make_pair(0.101195623867872 + -0.188854545332271 + 0.334217768087141 + -0.315849218854024, 1e-6), // lapltr_effective_incomplete(3.0, 8.0, 4.0)
                    std::make_pair(0.88941106522119,   1e-8), // C at Egamma=2.16
                    std::make_pair(0.92569114368575,   1e-8), // K_inv at Egamma=2.16
                    std::make_pair(0.882916019547774,  1e-8), // U at Egamma=2.16
                    std::make_pair(0.2597033704013,    1e-8), // F_leading_power(2.16)
                    std::make_pair(-0.0465690290256039 + 0.00392626152208998, 1e-5), // xi(2.16)
                    std::make_pair(0.0133804021454904 + 0.0, 1e-6), // delta_xi(2.16)
                };

                Diagnostics diagnostics = ff.diagnostics();
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);


                // Integration test: observable evaluation

                Kinematics k = Kinematics({ { "E_gamma", 2.16 } });
                Options o { {"form-factors", "FLvD2022QCDF"} };
                auto obs_F_V = Observable::make("B->gamma::F_V(E_gamma)", p, k, o);
                auto obs_F_A = Observable::make("B->gamma::F_A(E_gamma)", p, k, o);

                TEST_CHECK_NEARLY_EQUAL( ff.F_V(k["E_gamma"]),      0.2304410231690338,     1e-8 );
                TEST_CHECK_NEARLY_EQUAL( obs_F_V->evaluate(),       0.2304410231690338,     1e-8 );

                TEST_CHECK_NEARLY_EQUAL( ff.F_A(k["E_gamma"]),      0.2036809682121717,     1e-8 );
                TEST_CHECK_NEARLY_EQUAL( obs_F_A->evaluate(),       0.2036809682121717,     1e-8 );


                // Math integrity test: cross-check complete form factors against Mathematica implementation

                TEST_CHECK_NEARLY_EQUAL( ff.F_V(4.0),               0.1120729602946285,     1e-8 );
                TEST_CHECK_NEARLY_EQUAL( ff.F_V(12.0),              0.02616520856223276,    1e-8 );

                TEST_CHECK_NEARLY_EQUAL( ff.F_A(4.0),               0.109708682955886,      1e-8 );
                TEST_CHECK_NEARLY_EQUAL( ff.F_A(12.0),              0.02853002542154079,    1e-8 );
            }
        }
} analytic_b_to_gamma_qcdf_test;
