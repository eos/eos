/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013, 2015 Danny van Dyk
 * Copyright (c) 2024 Carolina Bolognani
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
#include <eos/b-decays/b-to-vec-l-nu.hh>
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

                // Resonance masses for the form-factors
                p["mass::B_d,1@BSZ2015"] = 5.723;

                Options oo
                {
                    { "model",        "WET"     },
                    { "form-factors", "BSZ2015" },
                    { "U",            "u"       },
                    { "I",            "1/2"     },
                    { "q",            "s"       },
                    { "l",            "mu"      }
                };

                BToVectorLeptonNeutrino d(p, oo);

                /* q^2 = [14.00, 19.21] */
                {
                    const double eps = 1e-4;
                    const auto ir = d.prepare(14.00,19.21);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_a_fb_leptonic(ir),             0.4120073563,       eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_amplitude_polarization_L(ir),  2.031721313e-17,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_amplitude_polarization_T(ir),  3.798507563e-17,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_f_L(ir),                       0.348480541,        eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(14.00,19.21),  0.0001007117913,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J1c(ir),                       9.122178885e-13,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J1s(ir),                       1.278178116e-12,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J2c(ir),                      -9.099664166e-13,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J2s(ir),                       4.25672643e-13,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J3(ir),                       -4.369577063e-13,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J4(ir),                        7.63268221e-13,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J5(ir),                       -2.21800317e-15,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J6c(ir),                       1.438237834e-12,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J6s(ir),                       0,                  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J7(ir),                        0,                  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J8(ir),                        0,                  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J9(ir),                        0,                  eps);
                }

                /* q^2 = [16.00, 19.21] */
                {
                    const double eps = 1e-4;
                    const auto ir = d.prepare(16.00,19.21);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_a_fb_leptonic(ir),             0.3955070687,       eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_amplitude_polarization_L(ir),  1.156941644e-17,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_amplitude_polarization_T(ir),  2.290102635e-17,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_f_L(ir),                       0.335632951,        eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(16.00,19.21),  5.954448983e-05,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J1c(ir),                       5.19407435e-13,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J1s(ir),                       7.706130679e-13,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J2c(ir),                      -5.183059103e-13,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J2s(ir),                       2.566522528e-13,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J3(ir),                       -3.044967454e-13,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J4(ir),                        4.598992481e-13,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J5(ir),                        4.627663896e-13,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J6c(ir),                      -1.068492527e-15,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J6s(ir),                       8.161887531e-13,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J7(ir),                        0,                  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J8(ir),                        0,                  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_J9(ir),                        0,                  eps);
                }

            }
        }
} bs_to_kstar_lepton_neutrino_low_recoil_test;
