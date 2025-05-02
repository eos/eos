/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Philip Lüghausen
 * Copyright (c) 2023-2025 Danny van Dyk
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
#include <eos/b-decays/b-to-gamma-l-nu.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/observable.hh>

using namespace test;
using namespace eos;

class BToGammaLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToGammaLeptonNeutrinoTest() :
            TestCase("b_to_gamma-l_nu_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-6;

            Parameters p = Parameters::Defaults();

            // Values from [BBJW:2018A]
            p["B_u::mu_0@FLvD2022"]  = 1.5;
            p["B->gamma::mu@FLvD2022QCDF"]  = 1.5;
            p["B->gamma::mu_h1@FLvD2022QCDF"]  = 4.7;
            p["B->gamma::mu_h2@FLvD2022QCDF"]  = 4.5;
            p["B->gamma::s_0@FLvD2022QCDF"] = 1.59;
            p["B->gamma::M^2@FLvD2022QCDF"] = 1.35;
            p["decay-constant::B_u"] = 0.192;
            p["mass::B_u"] = 5.27929;
            p["mass::b(MSbar)"] = 4.453796188717916; // fix m_b_pole@1-loop to 4.8
            p["mass::rho^+"] = 0.77526;
            p["B::lambda_E^2"] = 0.0625;
            p["B::lambda_H^2"] = 0.125;
            p["B::LambdaBar"] = 1.033232013955; // m_B - m_b
            p["CKM::abs(V_ub)"] = 3.7e-3;
            p["life_time::B_u"] = 1.638e-12;
            p["WET::G_Fermi"] = 1.166378e-5;

            // Compare to the exponential model
            p["B_u::omega_0@FLvD2022"]  = 0.2; // Corresponds to lambda_B
            p["B_u::a^phi+_0@FLvD2022"] = 1.0;
            p["B_u::a^phi+_1@FLvD2022"] = 0.0;
            p["B_u::a^phi+_2@FLvD2022"] = 0.0;
            p["B_u::a^phi+_3@FLvD2022"] = 0.0;
            p["B_u::a^phi+_4@FLvD2022"] = 0.0;
            p["B_u::a^phi+_5@FLvD2022"] = 0.0;
            p["B_u::a^phi+_6@FLvD2022"] = 0.0;
            p["B_u::a^phi+_7@FLvD2022"] = 0.0;
            p["B_u::a^phi+_8@FLvD2022"] = 0.0;

            Options oo
            {
                { "form-factors"_ok, "FLvD2022QCDF" },
            };

            // Branching ratio tests
            {
                auto BR = [p, oo](const double & E_gamma_min, const double & lambda_B) -> double {
                    auto k = Kinematics({ { "E_gamma_min", E_gamma_min } });
                    auto parameters = p; parameters["B_u::omega_0@FLvD2022"] = lambda_B; // a_k fixed to exp. model => omega_0 = lambda_B

                    return Observable::make("B_u->gammalnu::BR(E_gamma_min)", p, k, oo)->evaluate();
                };

                // Values taken from Ref. [BBJW:2018A], Fig. (9)
                TEST_CHECK_NEARLY_EQUAL(BR(1.0, 0.2), 5.0e-6, 0.5e-6);
                TEST_CHECK_NEARLY_EQUAL(BR(1.0, 0.3), 2.4e-6, 0.4e-6);
                TEST_CHECK_NEARLY_EQUAL(BR(1.0, 0.4), 1.2e-6, 0.2e-6);
                TEST_CHECK_NEARLY_EQUAL(BR(1.0, 0.6), 0.2e-6, 0.2e-6);

                TEST_CHECK_NEARLY_EQUAL(BR(1.5, 0.2), 3.2e-6,  0.5e-6);
                TEST_CHECK_NEARLY_EQUAL(BR(2.0, 0.2), 1.25e-6, 0.5e-6);
            }

            // Consistency check of the angular-differential decay witdh and the integrated decay width
            {
                BToGammaLeptonNeutrino obs(p, oo);
                const double E_gamma_min = 2.0;

                auto integrand = [&](const std::array<double, 2> & x) -> double {
                    const double E_gamma = std::get<0>(x);
                    const double costheta = std::get<1>(x);
                    return obs.fully_differential_decay_width(E_gamma, costheta);
                };

                const double decay_width_analytical = obs.integrated_decay_width(E_gamma_min);
                const double A_FB_analytical = obs.forward_backward_asymmetry(E_gamma_min);

                // Check analytical decay width numerically
                {
                    const double integrated_angular_decay_width = integrate<2>(integrand, { E_gamma_min, 1.0 }, { p["mass::B_u"] / 2.0, -1.0 });
                    TEST_CHECK_NEARLY_EQUAL(decay_width_analytical, integrated_angular_decay_width, eps);
                }

                // Check analytical implementation of A_FB numerically
                {
                    const double Gamma_F = integrate<2>(integrand, { E_gamma_min, 1.0 }, { p["mass::B_u"] / 2.0,  0.0 });
                    const double Gamma_B = integrate<2>(integrand, { E_gamma_min, 0.0 }, { p["mass::B_u"] / 2.0, -1.0 });
                    TEST_CHECK_NEARLY_EQUAL((Gamma_F - Gamma_B)/(Gamma_F + Gamma_B), A_FB_analytical, eps);
                }

            }

            // Diagnostics
            {
                static const std::vector<std::pair<double, double>> reference = {
                    std::make_pair(0.0, 1e-9), // Gamma_F + Gamma_B - Gamma
                };

                BToGammaLeptonNeutrino obs(p, oo);
                Diagnostics diagnostics = obs.diagnostics();
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }
        }
} b_to_gamma_l_nu_test;
