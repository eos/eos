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
                p["mass::B_d"]                               = 5.279;
                p["mass::K_d"]                               = 0.492;
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 4.0;
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
                    std::make_pair(   0.0,      eps),            // Re{1/phi_+(q2 = 0.0)}
                    std::make_pair(   0.0,      eps),            // Im{1/phi_+(q2 = 0.0)}
                    std::make_pair( -17.44509,  eps),            // Re{phi_+(q2 = 16.0)}
                    std::make_pair(   4.863096, eps),            // Im{phi_+(q2 = 16.0)}

                    std::make_pair( -18.00857,  eps),            // Re{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, {2.0, 3.0, 4.0, 5.0})}
                    std::make_pair(   0.0,      eps),            // Im{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, {2.0, 3.0, 4.0, 5.0})}


                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(-1.0)),  0.09205107389108,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(-1.0)),  0.11107379720400,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(0.0)),   0.,                 eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(0.0)),   0.,                 eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(4.0)),  -0.726740909982928,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(4.0)),  -0.868844878978099,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(12.0)),  7.94707073360654,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(12.0)),  9.306172848800037,  eps);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus_residue_jpsi()),   11.46205588287294,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus_residue_jpsi()),   13.52065260822002,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus_residue_psi2s()), -3.089134313454883,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus_residue_psi2s()), -3.595356292756863,  eps);
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
                p["mass::B_d"]                                = 5.279;
                p["mass::K_d"]                                = 0.492;
                p["mass::J/psi"]                              = 3.0969;
                p["mass::psi(2S)"]                            = 3.6860;
                p["mass::B_s^*"]                              = 5.4154;
                p["mass::D^0"]                                = 1.86723;
                p["b->sccbar::t_0"]                           = 4.0;
                p["b->sccbar::t_s"]                           = -17.4724;
                p["b->sccbar::chiOPE@GRvDV2021"]              = 1.81e-4;
                p["B->Kccbar::Re{alpha_0^plus}@GRvDV2021"]    = 2.0;
                p["B->Kccbar::Im{alpha_0^plus}@GRvDV2021"]    = 3.0;
                p["B->Kccbar::Re{alpha_1^plus}@GRvDV2021"]    = 4.0;
                p["B->Kccbar::Im{alpha_1^plus}@GRvDV2021"]    = 5.0;
                p["B->Kccbar::Re{alpha_2^plus}@GRvDV2021"]    = 6.0;
                p["B->Kccbar::Im{alpha_2^plus}@GRvDV2021"]    = 7.0;
                p["B->Kccbar::Re{alpha_3^plus}@GRvDV2021"]    = 0.0;
                p["B->Kccbar::Im{alpha_3^plus}@GRvDV2021"]    = 0.0;
                p["B->Kccbar::Re{alpha_4^plus}@GRvDV2021"]    = 0.0;
                p["B->Kccbar::Im{alpha_4^plus}@GRvDV2021"]    = 0.0;
                p["B->Kccbar::Re{alpha_5^plus}@GRvDV2021"]    = 0.0;
                p["B->Kccbar::Im{alpha_5^plus}@GRvDV2021"]    = 0.0;

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
                    std::make_pair( 0.8835558,  eps),  // Re{P(q2 = 1.0, 2.0, 3.0, 4.0)}
                    std::make_pair( 0.0,        eps),  // Im{P(q2 = 1.0, 2.0, 3.0, 4.0)}
                    std::make_pair( 0.8835558,  eps),  // Re{P(q2 = 1.0, (2.0,5.0), (3.0,6.0), (4.0,7.0))}
                    std::make_pair( 2.165236,   eps),  // Im{P(q2 = 1.0, (2.0,5.0), (3.0,6.0), (4.0,7.0))}
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);


                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(-1.0)),  0.0174583192607786,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(-1.0)),  0.0253239353342969,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(0.0)),   0.,                  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(0.0)),   0.,                  eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(4.0)),  -0.0959414569689401,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(4.0)),  -0.1439121854534102,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus(12.0)),  0.3846562670363007,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(12.0)),  0.6019561272742963,  eps);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus_residue_jpsi()),   0.7846262367768398,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus_residue_jpsi()),   1.242956592032815,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_plus_residue_psi2s()), -0.1495996738226507,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus_residue_psi2s()), -0.2065313797574414,  eps);
            }
        }
} nonlocal_formfactor_grvdv2021_test;

