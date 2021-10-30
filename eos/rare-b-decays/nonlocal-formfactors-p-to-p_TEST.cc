/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Nico Gubernari
 * Copyright (c) 2021 Méril Reboud
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
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

class NonlocalFormFactorGvDV2020Test :
    public TestCase
{
    public:
        NonlocalFormFactorGvDV2020Test() :
            TestCase("nonlocal_formfactor_GvDV2020_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                               = 5.27942;
                p["mass::K_d"]                               = 0.4936;
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 9.0;
                p["b->sccbar::t_s"]                          = -17.4724;
                p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;
                p["B->Kccbar::Re{alpha_0^plus}@GvDV2020"]    = 2.0;
                p["B->Kccbar::Im{alpha_0^plus}@GvDV2020"]    = 3.0;
                p["B->Kccbar::Re{alpha_1^plus}@GvDV2020"]    = 4.0;
                p["B->Kccbar::Im{alpha_1^plus}@GvDV2020"]    = 5.0;
                p["B->Kccbar::Re{alpha_2^plus}@GvDV2020"]    = 6.0;
                p["B->Kccbar::Im{alpha_2^plus}@GvDV2020"]    = 7.0;

                Options o = { { "model", "WET" } };

                auto nff = NonlocalFormFactor<nff::PToP>::make("B->K::GvDV2020", p, o);


                auto diagnostics = nff->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    /* outer functions */
                    std::make_pair( 0.0,      eps),            // Re{1/phi_+(q2 = 0.0)}
                    std::make_pair( 0.0,      eps),            // Im{1/phi_+(q2 = 0.0)}
                    std::make_pair( -16.3399, eps),            // Re{phi_+(q2 = 16.0)}
                    std::make_pair( 2.06922,  eps),            // Im{phi_+(q2 = 16.0)}

                    std::make_pair( 5.73999,  eps),            // Re{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, 2.0, 3.0, 4.0)}
                    std::make_pair( 0.0,      eps),            // Im{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, 2.0, 3.0, 4.0)}


                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(-1.0)),  0.131447,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(-1.0)),  0.156871,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(0.0)),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(0.0)),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(4.0)),  -1.09013,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(4.0)),  -1.2906,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(12.0)), 14.5673,  10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(12.0)), 16.9699,  10*eps);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus_residue_jpsi()),   19.0589,  10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus_residue_jpsi()),   22.3204,  10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus_residue_psi2s()), -6.68208,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus_residue_psi2s()), -7.75525,     eps);
            }
        }
} nonlocal_formfactor_gvdv2020_test;


class NonlocalFormFactorGvDV2021Test :
    public TestCase
{
    public:
        NonlocalFormFactorGvDV2021Test() :
            TestCase("nonlocal_formfactor_GvDV2021_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                                = 5.27942;
                p["mass::K_d"]                                = 0.4936;
                p["mass::J/psi"]                              = 3.0969;
                p["mass::psi(2S)"]                            = 3.6860;
                p["mass::B_s^*"]                              = 5.4154;
                p["mass::D^0"]                                = 1.86723;
                p["b->sccbar::t_0"]                           = 9.0;
                p["b->sccbar::t_s"]                           = -17.4724;
                p["b->sccbar::chiOPE@GRvDV2021"]              = 1.81e-4;
                p["B->Kccbar::Re{alpha_0^plus}@GRvDV2021"]    = 2.0;
                p["B->Kccbar::Im{alpha_0^plus}@GRvDV2021"]    = 3.0;
                p["B->Kccbar::Re{alpha_1^plus}@GRvDV2021"]    = 4.0;
                p["B->Kccbar::Im{alpha_1^plus}@GRvDV2021"]    = 5.0;
                p["B->Kccbar::Re{alpha_2^plus}@GRvDV2021"]    = 6.0;
                p["B->Kccbar::Im{alpha_2^plus}@GRvDV2021"]    = 7.0;

                Options o = { { "model", "WET" } };

                auto nff = NonlocalFormFactor<nff::PToP>::make("B->K::GRvDV2021", p, o);


                auto diagnostics = nff->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair( 5.31175,  eps),            // Re{phi_+(q2 = 16.0)}
                    std::make_pair( -3.8452,  eps),            // Im{phi_+(q2 = 16.0)}

                    std::make_pair( 1.16924,  eps),            // Re{P(q2 = 1.0, 2.0, 3.0, 4.0)}
                    std::make_pair( 0.0,      eps),            // Im{P(q2 = 1.0, 2.0, 3.0, 4.0)}
                    std::make_pair( 1.16924,  eps),            // Re{P(q2 = 1.0, (2.0,5.0), (3.0,6.0), (4.0,7.0))}
                    std::make_pair( 2.71519,  eps),            // Im{P(q2 = 1.0, (2.0,5.0), (3.0,6.0), (4.0,7.0))}

                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);


                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(-1.0)), -0.071174,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(-1.0)), -0.0983567,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(0.0)),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(0.0)),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(4.0)),   0.410774,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(4.0)),   0.58287,    eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(12.0)), -1.18903,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(12.0)), -1.88908,    eps);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus_residue_jpsi()),  -3.14489,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus_residue_jpsi()),  -4.76717,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus_residue_psi2s()),  0.334561,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus_residue_psi2s()),  0.484645,  eps);
            }
        }
} nonlocal_formfactor_grvdv2021_test;

