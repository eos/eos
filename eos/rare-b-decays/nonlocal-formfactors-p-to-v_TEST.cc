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
                p["mass::K_d^*"]                             = 0.89555;
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 9.0;
                p["b->sccbar::t_s"]                          = -17.4724;
                p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;
                p["B->K^*ccbar::Re{alpha_0^perp}@GvDV2020"]  = 2.0;
                p["B->K^*ccbar::Im{alpha_0^perp}@GvDV2020"]  = 3.0;
                p["B->K^*ccbar::Re{alpha_1^perp}@GvDV2020"]  = 4.0;
                p["B->K^*ccbar::Im{alpha_1^perp}@GvDV2020"]  = 5.0;
                p["B->K^*ccbar::Re{alpha_2^perp}@GvDV2020"]  = 6.0;
                p["B->K^*ccbar::Im{alpha_2^perp}@GvDV2020"]  = 7.0;
                p["B->K^*ccbar::Re{alpha_0^para}@GvDV2020"]  = 8.0;
                p["B->K^*ccbar::Im{alpha_0^para}@GvDV2020"]  = 9.0;
                p["B->K^*ccbar::Re{alpha_1^para}@GvDV2020"]  = 10.0;
                p["B->K^*ccbar::Im{alpha_1^para}@GvDV2020"]  = 11.0;
                p["B->K^*ccbar::Re{alpha_2^para}@GvDV2020"]  = 12.0;
                p["B->K^*ccbar::Im{alpha_2^para}@GvDV2020"]  = 13.0;
                p["B->K^*ccbar::Re{alpha_0^long}@GvDV2020"]  = 14.0;
                p["B->K^*ccbar::Im{alpha_0^long}@GvDV2020"]  = 15.0;
                p["B->K^*ccbar::Re{alpha_1^long}@GvDV2020"]  = 16.0;
                p["B->K^*ccbar::Im{alpha_1^long}@GvDV2020"]  = 17.0;
                p["B->K^*ccbar::Re{alpha_2^long}@GvDV2020"]  = 18.0;
                p["B->K^*ccbar::Im{alpha_2^long}@GvDV2020"]  = 19.0;

                Options o = {
		    { "model", "WET" },
		    { "q", "d" }
		};

                auto nff = NonlocalFormFactor<nff::PToV>::make("B->K^*::GvDV2020", p, o);


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
                    std::make_pair( 0.0,      eps),            // Re{1/phi_long(q2 = 0.0)}
                    std::make_pair( 0.0,      eps),            // Im{1/phi_long(q2 = 0.0)}

                    std::make_pair( -36.5755, 10*eps),         // Re{phi_long(q2 = 16.0)}
                    std::make_pair( 4.63177,  10*eps),         // Im{phi_long(q2 = 16.0)}

                    std::make_pair( 24.6148,  10*eps),         // Re{phi_perp(q2 = 16.0)}
                    std::make_pair( -13.2048, 10*eps)          // Im{phi_perp(q2 = 16.0)}
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp(16.0)),  -2.36353,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp(16.0)),  -1.27642,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_para(16.0)),  -4.48563,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para(16.0)),  -2.2198,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_long(16.0)),   5.53271,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long(16.0)),   0.443831,    eps);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp_residue_jpsi()),   -52.3353,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp_residue_jpsi()),   -61.0889,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp_residue_psi2s()),    7.67603,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp_residue_psi2s()),    8.90398,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_para_residue_jpsi()),  -104.857,   100*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para_residue_jpsi()),  -113.610,   100*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_para_residue_psi2s()),   15.0437,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para_residue_psi2s()),   16.2717,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_long_residue_jpsi()),    49.1116,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long_residue_jpsi()),    51.8432,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_long_residue_psi2s()),  -13.3074,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long_residue_psi2s()),  -14.0365,   10*eps);
            }
        }
} nonlocal_formfactor_gvdv2020_test;


class NonlocalFormFactorGRvDV2021Test :
    public TestCase
{
    public:
        NonlocalFormFactorGRvDV2021Test() :
            TestCase("nonlocal_formfactor_GRvDV2021_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                                = 5.27942;
                p["mass::K_d^*"]                              = 0.89555;
                p["mass::J/psi"]                              = 3.0969;
                p["mass::psi(2S)"]                            = 3.6860;
                p["mass::B_s^*"]                              = 5.4154;
                p["mass::D^0"]                                = 1.86723;
                p["b->sccbar::t_0"]                           = 9.0;
                p["b->sccbar::t_s"]                           = -17.4724;
                p["b->sccbar::chiOPE@GRvDV2021"]              = 1.81e-4;
                p["B->K^*ccbar::Re{alpha_0^perp}@GRvDV2021"]  = 2.0;
                p["B->K^*ccbar::Im{alpha_0^perp}@GRvDV2021"]  = 3.0;
                p["B->K^*ccbar::Re{alpha_1^perp}@GRvDV2021"]  = 4.0;
                p["B->K^*ccbar::Im{alpha_1^perp}@GRvDV2021"]  = 5.0;
                p["B->K^*ccbar::Re{alpha_2^perp}@GRvDV2021"]  = 6.0;
                p["B->K^*ccbar::Im{alpha_2^perp}@GRvDV2021"]  = 7.0;
                p["B->K^*ccbar::Re{alpha_0^para}@GRvDV2021"]  = 8.0;
                p["B->K^*ccbar::Im{alpha_0^para}@GRvDV2021"]  = 9.0;
                p["B->K^*ccbar::Re{alpha_1^para}@GRvDV2021"]  = 10.0;
                p["B->K^*ccbar::Im{alpha_1^para}@GRvDV2021"]  = 11.0;
                p["B->K^*ccbar::Re{alpha_2^para}@GRvDV2021"]  = 12.0;
                p["B->K^*ccbar::Im{alpha_2^para}@GRvDV2021"]  = 13.0;
                p["B->K^*ccbar::Re{alpha_0^long}@GRvDV2021"]  = 14.0;
                p["B->K^*ccbar::Im{alpha_0^long}@GRvDV2021"]  = 15.0;
                p["B->K^*ccbar::Re{alpha_1^long}@GRvDV2021"]  = 16.0;
                p["B->K^*ccbar::Im{alpha_1^long}@GRvDV2021"]  = 17.0;
                p["B->K^*ccbar::Re{alpha_2^long}@GRvDV2021"]  = 18.0;
                p["B->K^*ccbar::Im{alpha_2^long}@GRvDV2021"]  = 19.0;

                Options o = {
                    { "model", "WET" },
                    { "q", "d"}
                };
                auto nff = NonlocalFormFactor<nff::PToV>::make("B->K^*::GRvDV2021", p, o);


                auto diagnostics = nff->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair( 11.8899,  eps),         // Re{phi_long(q2 = 16.0)}
                    std::make_pair( -8.60714, eps),         // Im{phi_long(q2 = 16.0)}

                    std::make_pair( -6.07403, eps),         // Re{phi_perp(q2 = 16.0)}
                    std::make_pair( 9.3159,   eps)          // Im{phi_perp(q2 = 16.0)}
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp(-1.)), -2.9621,      eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp(-1.)), -4.09339,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_para(-1.)), -9.74979,  10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para(-1.)), -10.8811,  10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_long(-1.)), -0.412137,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long(-1.)), -0.44033,     eps);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp(0.)), -2.96771,      eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp(0.)), -4.11808,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_para(0.)), -9.8699,    10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para(0.)), -11.0203,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_long(0.)),   0.,          eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long(0.)),   0.,          eps);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp(4.)), -3.20673,      eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp(4.)), -4.55021,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_para(4.)), -11.2676,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para(4.)), -12.6111,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_long(4.)),   2.126,       eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long(4.)),   2.27378,     eps);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp(12.)),   1.54707,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp(12.)),   2.4579,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_para(12.)),  7.01207,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para(12.)),  7.9229,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_long(12.)), -5.5288,      eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long(12.)), -5.93241,     eps);

                TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp_residue_jpsi()),    6.65637,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp_residue_jpsi()),   10.09,        eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp_residue_psi2s()),  -0.294468,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp_residue_psi2s()),  -0.426566,    eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_para_residue_jpsi()),   27.2582,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para_residue_jpsi()),   30.6919,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_para_residue_psi2s()),  -1.08705,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para_residue_psi2s()),  -1.21915,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_long_residue_jpsi()),  -14.9353,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long_residue_jpsi()),  -16.0068,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nff->H_long_residue_psi2s()),   1.11609,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long_residue_psi2s()),   1.19452,     eps);
            }
        }
} nonlocal_formfactor_grvdv2021_test;
