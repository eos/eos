/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Nico Gubernari
 * Copyright (c) 2021 Méril Reboud
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
#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>

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
                p["mass::K_d^*"]                             = 0.896;
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 4.0;
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
                    { "model"_ok, "WET" },
                    { "q"_ok, "d" }
                };

                auto nff = NonlocalFormFactor<PToV>::make("B->K^*::GvDV2020", p, o);


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
                    std::make_pair( 0.0,       eps),        // Re{1/phi_long(q2 = 0.0)}
                    std::make_pair( 0.0,       eps),        // Im{1/phi_long(q2 = 0.0)}

                    std::make_pair(-39.01168,  eps),        // Re{phi_long(q2 = 16.0)}
                    std::make_pair( 10.87513,  eps),        // Im{phi_long(q2 = 16.0)}

                    std::make_pair( 24.64525,  eps),        // Re{phi_perp(q2 = 16.0)}
                    std::make_pair(-18.28392,  eps)         // Im{phi_perp(q2 = 16.0)}
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_RELATIVE_ERROR(real(nff->H_perp(16.0)), -0.765303,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_perp(16.0)), -1.031884,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_para(16.0)), -1.475739,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_para(16.0)), -1.888891,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_long(16.0)),  2.291979,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_long(16.0)),  1.349007,   eps);

                TEST_CHECK_RELATIVE_ERROR(real(nff->H_perp_residue_jpsi()),  -31.142,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_perp_residue_jpsi()),  -36.5503,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_perp_residue_psi2s()), 3.50665,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_perp_residue_psi2s()), 4.0751,    eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_para_residue_jpsi()),  -63.5918,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_para_residue_jpsi()),  -69.0001,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_para_residue_psi2s()), 6.91734,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_para_residue_psi2s()), 7.48579,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_long_residue_jpsi()),  29.9732,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_long_residue_jpsi()),  31.6611,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_long_residue_psi2s()), -6.13303,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_long_residue_psi2s()), -6.47059,  eps);
            }
        }
} nonlocal_formfactor_gvdv2020_test;
