/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013, 2015 Danny van Dyk
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
#include <eos/b-decays/bs-to-kstar-l-nu.hh>
#include <eos/maths/complex.hh>

using namespace test;
using namespace eos;

class BsToKstarLeptonNeutrinoTest :
    public TestCase
{
    public:
        BsToKstarLeptonNeutrinoTest() :
            TestCase("bs_to_kstar_l_nu_test")
        {
        }

        virtual void run() const
        {
            /* Low Recoil (SM) */
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_d"] = 1.516e-12;
                // PDG 2012 CKM parameters
                p["CKM::A"] = 0.827;
                p["CKM::lambda"] = 0.22535;
                p["CKM::rhobar"] = 0.132;
                p["CKM::etabar"] = 0.340;
                // CKM matrix elements corresponding to the above Wolfenstein parameters
                p["CKM::abs(V_ub)"] =  0.003540609803917236;
                p["CKM::arg(V_ub)"] = -1.2010727175261147;

                // Kaon mass
                p["mass::K_u^*"] = 0.89166;
                // B mass
                p["mass::B_s"] = 5.3668;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                // mu lepton mass
                p["mass::mu"] = 0.1056583715;

                Options oo;
                oo.declare("model", "WET");
                oo.declare("form-factors", "BSZ2015");

                BsToKstarLeptonNeutrino d(p, oo);

                /* q^2 = [14.00, 19.21] */
                {
                    const double eps = 1e-4;

                    TEST_CHECK_NEARLY_EQUAL(-0.4125863683, d.integrated_forward_backward_asymmetry(14.00, 19.21), eps);
                    TEST_CHECK_NEARLY_EQUAL( 0.3482936714, d.integrated_longitudinal_polarisation(14.00, 19.21),  eps);
                    TEST_CHECK_NEARLY_EQUAL(-0.5132406718, d.integrated_transverse_asymmetry_2(14.00, 19.21),     eps);
                    TEST_CHECK_NEARLY_EQUAL( 1.7577913835, d.integrated_transverse_asymmetry_3(14.00, 19.21),     eps);
                    TEST_CHECK_NEARLY_EQUAL( 0.5655315082, d.integrated_transverse_asymmetry_4(14.00, 19.21),     eps);
                    TEST_CHECK_NEARLY_EQUAL( 0.0775462723, d.integrated_transverse_asymmetry_5(14.00, 19.21),     eps);
                    TEST_CHECK_NEARLY_EQUAL(-0.8441151078, d.integrated_transverse_asymmetry_re(14.00, 19.21),    eps);
                    TEST_CHECK_NEARLY_EQUAL( 0,            d.integrated_transverse_asymmetry_im(14.00, 19.21),    eps);
                    TEST_CHECK_NEARLY_EQUAL( 0.9969214819, d.integrated_h_1(14.00, 19.21),                        eps);
                    TEST_CHECK_NEARLY_EQUAL(-0.9940864123, d.integrated_h_2(14.00, 19.21),                        eps);
                    TEST_CHECK_NEARLY_EQUAL(-0.9835366074, d.integrated_h_3(14.00, 19.21),                        eps);
                    TEST_CHECK_NEARLY_EQUAL( 0,            d.integrated_h_4(14.00, 19.21),                        eps);
                    TEST_CHECK_NEARLY_EQUAL(-0,            d.integrated_h_5(14.00, 19.21),                        eps);
                }

                /* q^2 = [16.00, 19.21] */
                {
                    const double eps = 1e-4;
                    TEST_CHECK_NEARLY_EQUAL(-0.3959778457, d.integrated_forward_backward_asymmetry(16.00, 19.21), eps);
                    TEST_CHECK_NEARLY_EQUAL( 0.3354919677, d.integrated_longitudinal_polarisation(16.00, 19.21),  eps);
                    TEST_CHECK_NEARLY_EQUAL(-0.5932022373, d.integrated_transverse_asymmetry_2(16.00, 19.21),     eps);
                    TEST_CHECK_NEARLY_EQUAL( 1.9770010813, d.integrated_transverse_asymmetry_3(16.00, 19.21),     eps);
                    TEST_CHECK_NEARLY_EQUAL( 0.5022878134, d.integrated_transverse_asymmetry_4(16.00, 19.21),     eps);
                    TEST_CHECK_NEARLY_EQUAL( 0.0648764771, d.integrated_transverse_asymmetry_5(16.00, 19.21),     eps);
                    TEST_CHECK_NEARLY_EQUAL(-0.7945283357, d.integrated_transverse_asymmetry_re(16.00, 19.21),    eps);
                    TEST_CHECK_NEARLY_EQUAL( 0,            d.integrated_transverse_asymmetry_im(16.00, 19.21),    eps);
                    TEST_CHECK_NEARLY_EQUAL( 0.9989890567, d.integrated_h_1(16.00, 19.21),                        eps);
                    TEST_CHECK_NEARLY_EQUAL(-0.9930235503, d.integrated_h_2(16.00, 19.21),                        eps);
                    TEST_CHECK_NEARLY_EQUAL(-0.9869261553, d.integrated_h_3(16.00, 19.21),                        eps);
                    TEST_CHECK_NEARLY_EQUAL( 0,            d.integrated_h_4(16.00, 19.21),                        eps);
                    TEST_CHECK_NEARLY_EQUAL(-0,            d.integrated_h_5(16.00, 19.21),                        eps);
                }
            }
        }
} bs_to_kstar_lepton_neutrino_low_recoil_test;
