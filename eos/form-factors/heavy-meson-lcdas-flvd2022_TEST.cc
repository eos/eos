/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022-2026 Danny van Dyk
 * Copyright (c) 2022-2023 Philip Lüghausen
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

#include <eos/form-factors/heavy-meson-lcdas-flvd2022.hh>
#include <eos/observable.hh>

#include <test/test.hh>

#include <numeric>

using namespace test;
using namespace eos;
using namespace heavy_meson_lcdas;

class FLvD2022Test : public TestCase
{
    public:
        FLvD2022Test() :
            TestCase("b_lcdas_flvd2022_test")
        {
        }

        virtual void
        run() const
        {
            // basic test
            {
                Parameters p = Parameters::Defaults();
                Options    o{
                       { "alpha-s"_ok, "naive"_ov }
                };

                p["B_u::mu_0@FLvD2022"]    = 1.0;
                p["B_u::omega_0@FLvD2022"] = 1.0;

                std::array<double, 9> ref = { 1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0 };
                for (size_t k = 0; k < ref.size(); k++)
                {
                    p["B_u::a^phi+_" + std::to_string(k) + "@FLvD2022"] = ref[k];
                }

                FLvD2022 blcdas(p, o);

                auto [c, c_end] = blcdas.coefficient_range(1.0);
                for (auto it = c; it != c_end; ++it)
                {
                    TEST_CHECK_NEARLY_EQUAL(*it, ref[std::distance(c, it)], 1e-15);
                }
            }

            // RG evolution
            {
                Parameters p = Parameters::Defaults();
                Options    o{
                       { "alpha-s"_ok, "naive"_ov }
                };

                p["B_u::mu_0@FLvD2022"]    = 1.0;
                p["B_u::omega_0@FLvD2022"] = 0.55;
                const double mu            = 1.5;

                // Regression reference for the coefficient evolution implemented in coefficient_range().
                // The inverse-moment test below provides the independent physics check.
                std::array<double, 9> ref         = { 1.6, 0.8, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0 };
                std::array<double, 9> res_evolved = { 1.3684832595861576,  0.5055458502529527, 0.2415628835292180,  0.1069513981135003, 0.0908610730053892,
                                                      -0.0238057216522631, 0.0190678282263298, -0.0164456740588282, 0.0145822817501768 };
                for (size_t k = 0; k < ref.size(); k++)
                {
                    p["B_u::a^phi+_" + std::to_string(k) + "@FLvD2022"] = ref[k];
                }

                FLvD2022 blcdas(p, o);

                auto [c, c_end] = blcdas.coefficient_range(mu);
                for (auto it = c; it != c_end; ++it)
                {
                    TEST_CHECK_NEARLY_EQUAL(*it, res_evolved[std::distance(c, it)], 1e-7);
                }
            }

            // Evolution of the inverse moment for the exponential model. The corresponding
            // O(alpha_s) values from [BR:2011A], Eq. (3.1), are 0.8492, 0.8938, and 0.9223.
            // These pure inverse-moment ratios differ from breakdown-FV.py's leading-power
            // ratios, which additionally include the evolution-induced change in the jet function J.
            {
                Parameters p = Parameters::Defaults();
                Options    o;

                p["B_u::mu_0@FLvD2022"]              = 1.0;
                p["QCD::alpha_s(MZ)"]                = 0.11739977629025072; // reproduces the [BBJW:2018A] Table 1 coupling
                const double                mu       = 1.5;
                const std::array<double, 3> omega_0  = { 0.20, 0.35, 0.50 };
                const std::array<double, 3> expected = { 0.8428320729585771, 0.8846413359198045, 0.9123645542739210 };

                for (unsigned i = 0; i < omega_0.size(); ++i)
                {
                    p["B_u::omega_0@FLvD2022"] = omega_0[i];
                    FLvD2022     blcdas(p, o);
                    const double evolution_factor = blcdas.inverse_moment(mu) / blcdas.inverse_moment(1.0);
                    TEST_CHECK_NEARLY_EQUAL(evolution_factor, expected[i], 1e-4);
                }
            }

            // pseudo-observables tildephi and the derivative
            {
                Parameters            p          = Parameters::Defaults();
                std::array<double, 9> parameters = { 1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0 };

                p["B_u::mu_0@FLvD2022"]    = 1.0;
                p["B_u::omega_0@FLvD2022"] = 0.3;
                for (size_t k = 0; k < parameters.size(); k++)
                {
                    p["B_u::a^phi+_" + std::to_string(k) + "@FLvD2022"] = parameters[k];
                }

                Kinematics k = Kinematics({
                    { "tau",                     0.4 },
                    {  "mu", p["B_u::mu_0@FLvD2022"] }
                });
                Options    o{
                       { "alpha-s"_ok, "naive"_ov }
                };

                auto phitilde = Observable::make("B::phitilde_+(-i*tau,mu)@FLvD2022", p, k, o);
                TEST_CHECK_NEARLY_EQUAL(phitilde->evaluate(), 11.558610659381285, 1e-12);

                auto t_d_dt_phitilde = Observable::make("B::tau*d_dtau_phitilde_+(-i*tau,mu)@FLvD2022", p, k, o);
                TEST_CHECK_NEARLY_EQUAL(t_d_dt_phitilde->evaluate(), -13.812451430614884, 1e-12);

                auto t2_d2_d2t_phitilde = Observable::make("B::tau^2*d2_d2tau_phitilde_+(-i*tau,mu)@FLvD2022", p, k, o);
                TEST_CHECK_NEARLY_EQUAL(t2_d2_d2t_phitilde->evaluate(), 20.320360292444562, 1e-12);
            }

            // inverse moment and first and second logarithmic moments
            {
                Parameters            p          = Parameters::Defaults();
                std::array<double, 9> parameters = { 1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0 };

                p["B_u::mu_0@FLvD2022"]    = 1.0;
                p["B_u::omega_0@FLvD2022"] = 0.3;
                for (size_t k = 0; k < parameters.size(); k++)
                {
                    p["B_u::a^phi+_" + std::to_string(k) + "@FLvD2022"] = parameters[k];
                }

                Kinematics k = Kinematics({
                    { "mu", 1.0 }
                });
                Options    o{
                       { "alpha-s"_ok, "naive"_ov }
                };

                auto L0 = Observable::make("B::L0@FLvD2022", p, k, o);
                TEST_CHECK_NEARLY_EQUAL(L0->evaluate(), 16.66666666667032, 1e-10);

                auto L1 = Observable::make("B::L1@FLvD2022", p, k, o);
                TEST_CHECK_NEARLY_EQUAL(L1->evaluate(), 36.95238095239132, 1e-10);

                auto L2 = Observable::make("B::L2@FLvD2022", p, k, o);
                TEST_CHECK_NEARLY_EQUAL(L2->evaluate(), 126.63249899776702, 1e-10);
            }
        }
} b_lcdas_flvd2022_test;
