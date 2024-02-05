/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018, 2019 Ahmet Kokulu
 * Copyright (c) 2019 Danny van Dyk
 * Copyright (c) 2021 Christoph Bobeth
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
#include <eos/b-decays/b-to-psd-l-nu.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <vector>


using namespace test;
using namespace eos;

class BToDLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToDLeptonNeutrinoTest() :
            TestCase("b_to_d_l_nu_test")
        {
        }

        virtual void run() const
        {
            // comparison with Martin Jung in 3/2/1 model
            {
                Parameters p = Parameters::Defaults();

                p["B(*)->D(*)::xi'(1)@HQET"].set(-1.06919);
                p["B(*)->D(*)::xi''(1)@HQET"].set(1.66581);
                p["B(*)->D(*)::xi'''(1)@HQET"].set(-2.91356);
                p["B(*)->D(*)::chi_2(1)@HQET"].set(-0.0600518);
                p["B(*)->D(*)::chi_2'(1)@HQET"].set(-0.0000101998);
                p["B(*)->D(*)::chi_2''(1)@HQET"].set(-0.085385);
                p["B(*)->D(*)::chi_3'(1)@HQET"].set(0.0400584);
                p["B(*)->D(*)::chi_3''(1)@HQET"].set(-0.0522346);
                p["B(*)->D(*)::eta(1)@HQET"].set(0.586099);
                p["B(*)->D(*)::eta'(1)@HQET"].set(-0.0233426);
                p["B(*)->D(*)::eta''(1)@HQET"].set(-0.0288193);
                p["B(*)->D(*)::l_1(1)@HQET"].set(0.113962);
                p["B(*)->D(*)::l_1'(1)@HQET"].set(-8.15957);
                p["B(*)->D(*)::l_2(1)@HQET"].set(-1.90706);
                p["B(*)->D(*)::l_2'(1)@HQET"].set(-3.16227);
                p["B(*)->D(*)::l_3(1)@HQET"].set(-3.41861);
                p["B(*)->D(*)::l_3'(1)@HQET"].set(5.6966);
                p["B(*)->D(*)::l_4(1)@HQET"].set(-1.89465);
                p["B(*)->D(*)::l_4'(1)@HQET"].set(0.220492);
                p["B(*)->D(*)::l_5(1)@HQET"].set(4.97017);
                p["B(*)->D(*)::l_5'(1)@HQET"].set(-2.34767);
                p["B(*)->D(*)::l_6(1)@HQET"].set(1.98608);
                p["B(*)->D(*)::l_6'(1)@HQET"].set(1.42747);

                p["CKM::abs(V_cb)"].set(1.0);
                p["mass::e"].set(0.000001);
                p["B(*)->D(*)::a@HQET"].set(1.000);
                p["mass::B_d"].set(5.27942);
                p["mass::D_u"].set(1.86723);
                p["mass::D_d"].set(1.86723);
                p["life_time::B_d"].set(1.520e-12);

                Options o{
                    { "l",             "e"         },
                    { "model",         "CKM"       },
                    { "U",             "c"         },
                    { "q",             "d"         },
                    { "I",             "1/2"       },
                    { "z-order-lp",    "3"         },
                    { "z-order-slp",   "2"         },
                    { "z-order-sslp",  "1"         },
                    { "form-factors",  "BGJvD2019" }
                };

                BToPseudoscalarLeptonNeutrino d(p, o);

                const double eps = 1e-3;
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(0.001, 11.643), 13.462, eps);
            }

            // comparison with Martin Jung in 3/2/1 model
            {
                Parameters p = Parameters::Defaults();

                p["B(*)->D(*)::xi'(1)@HQET"].set(-1.06919);
                p["B(*)->D(*)::xi''(1)@HQET"].set(1.66581);
                p["B(*)->D(*)::xi'''(1)@HQET"].set(-2.91356);
                p["B(*)->D(*)::chi_2(1)@HQET"].set(-0.0600518);
                p["B(*)->D(*)::chi_2'(1)@HQET"].set(-0.0000101998);
                p["B(*)->D(*)::chi_2''(1)@HQET"].set(-0.085385);
                p["B(*)->D(*)::chi_3'(1)@HQET"].set(0.0400584);
                p["B(*)->D(*)::chi_3''(1)@HQET"].set(-0.0522346);
                p["B(*)->D(*)::eta(1)@HQET"].set(0.586099);
                p["B(*)->D(*)::eta'(1)@HQET"].set(-0.0233426);
                p["B(*)->D(*)::eta''(1)@HQET"].set(-0.0288193);
                p["B(*)->D(*)::l_1(1)@HQET"].set(0.113962);
                p["B(*)->D(*)::l_1'(1)@HQET"].set(-8.15957);
                p["B(*)->D(*)::l_2(1)@HQET"].set(-1.90706);
                p["B(*)->D(*)::l_2'(1)@HQET"].set(-3.16227);
                p["B(*)->D(*)::l_3(1)@HQET"].set(-3.41861);
                p["B(*)->D(*)::l_3'(1)@HQET"].set(5.6966);
                p["B(*)->D(*)::l_4(1)@HQET"].set(-1.89465);
                p["B(*)->D(*)::l_4'(1)@HQET"].set(0.220492);
                p["B(*)->D(*)::l_5(1)@HQET"].set(4.97017);
                p["B(*)->D(*)::l_5'(1)@HQET"].set(-2.34767);
                p["B(*)->D(*)::l_6(1)@HQET"].set(1.98608);
                p["B(*)->D(*)::l_6'(1)@HQET"].set(1.42747);

                p["CKM::abs(V_cb)"].set(1.0);
                p["mass::e"].set(0.000001);
                p["B(*)->D(*)::a@HQET"].set(1.000);
                p["mass::B_d"].set(5.27942);
                p["mass::D_u"].set(1.86723);
                p["mass::D_d"].set(1.86723);
                p["life_time::B_d"].set(1.520e-12);

                Options o{
                    { "l",             "tau"       },
                    { "model",         "CKM"       },
                    { "U",             "c"         },
                    { "q",             "d"         },
                    { "I",             "1/2"       },
                    { "z-order-lp",    "3"         },
                    { "z-order-slp",   "2"         },
                    { "z-order-sslp",  "1"         },
                    { "form-factors",  "BGJvD2019" }
                };

                BToPseudoscalarLeptonNeutrino d(p, o);

                const double eps = 1e-3;
                TEST_CHECK_NEARLY_EQUAL(d.integrated_lepton_polarization(3.157, 11.643), +0.320914, eps);
            }

            // SM tests
            {
                Parameters p = Parameters::Defaults();
                p["B->D::f_+(0)@BCL2008"]  = +0.660;
                p["B->D::f_T(0)@BCL2008"]  = +0.00;
                p["B->D::b_+^1@BCL2008"]   = -4.00;
                p["B->D::b_+^2@BCL2008"]   = -0.80;
                p["B->D::b_0^1@BCL2008"]   = +0.40;
                p["B->D::b_0^2@BCL2008"]   = -1.20;
                p["B->D::b_T^1@BCL2008"]   = +0.00;
                p["B->D::b_T^2@BCL2008"]   = +0.00;
                p["mass::B_d"]             =  5.279;
                p["mass::D_d"]             =  1.870;
                // by default, all other couplings are zero in eos
                p["CKM::abs(V_cb)"]      =  0.041996951916414726;
                p["cbmunumu::Re{cVL}"]   =  1.0066;  // include Sirlin correction
                p["cbtaunutau::Re{cVL}"] =  1.0066;  // include Sirlin correction

                Options oo
                {
                    { "model",        "WET" },
                    { "form-factors", "BCL2008"    },
                    { "U",            "c"          },
                    { "q",            "d"          },
                    { "I",            "1/2"        },
                    { "l",            "mu"         }
                };


                const double eps = 1e-3;
                {
                    BToPseudoscalarLeptonNeutrino d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(0.011164, 11.62), 13.1988, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(0.011164, 11.62), -0.0138762, eps);

                    oo.declare("l", "tau");
                    auto k_tau = Kinematics{
                        { "q2_min",  3.15702  },
                        { "q2_max", 11.62     }
                    };
                    auto obs_BRtau = Observable::make("B->Dlnu::BR", p, k_tau, oo);
                    TEST_CHECK(obs_BRtau.get() != nullptr);
                    TEST_CHECK_RELATIVE_ERROR(obs_BRtau->evaluate(), 0.0069634, eps);

                    oo.declare("l", "mu");
                    auto k_mu = Kinematics{
                        { "q2_min",   0.011164 },
                        { "q2_max",  11.62     },
                    };
                    auto obs_BRmu = Observable::make("B->Dlnu::BR", p, k_mu, oo);
                    TEST_CHECK(obs_BRmu.get() != nullptr);
                    TEST_CHECK_RELATIVE_ERROR(obs_BRmu->evaluate(), 0.0232794, eps);

                    oo =
                    {
                        { "model",        "WET" },
                        { "form-factors", "BCL2008"    },
                        { "U",            "c"          },
                        { "q",            "d"          },
                        { "I",            "1/2"        }
                    };
                    auto k = Kinematics{
                        { "q2_mu_min",   0.011164 },
                        { "q2_mu_max",  11.62     },
                        { "q2_tau_min",  3.15702  },
                        { "q2_tau_max", 11.62     }
                    };
                    auto obs_RD = Observable::make("B->Dlnu::R_D", p, k, oo);
                    TEST_CHECK(obs_RD.get() != nullptr);
                    TEST_CHECK_RELATIVE_ERROR(obs_RD->evaluate(), 0.299132, eps);
                }
            }

            // NP tests
            {
                const double etaEW = 1.0066;

                Parameters p = Parameters::Defaults();
                p["B->D::f_+(0)@BCL2008"]  = +0.660;
                p["B->D::f_T(0)@BCL2008"]  = +1.00;
                p["B->D::b_+^1@BCL2008"]   = -4.00;
                p["B->D::b_+^2@BCL2008"]   = -0.800;
                p["B->D::b_0^1@BCL2008"]   = +0.400;
                p["B->D::b_0^2@BCL2008"]   = -1.20;
                p["B->D::b_T^1@BCL2008"]   = +3.00;
                p["B->D::b_T^2@BCL2008"]   = -0.60;
                p["mass::B_d"]             =  5.279;
                p["mass::D_d"]             =  1.870;
                // fix the scale
                p["cbmunumu::mu"]          =  4.18;
                p["cbtaunutau::mu"]        =  4.18;
                p["mass::b(MSbar)"]        =  4.18;
                p["mass::c"]               =  1.275;
                // CKM
                p["CKM::abs(V_cb)"]        =  0.041996951916414726;
                // mu mode
                p["cbmunumu::Re{cVL}"]         = +1.0 * etaEW;
                p["cbmunumu::Im{cVL}"]         = -2.0 * etaEW;
                p["cbmunumu::Re{cVR}"]         = +2.0 * etaEW;
                p["cbmunumu::Im{cVR}"]         = -2.0 * etaEW;
                p["cbmunumu::Re{cSL}"]         = +3.0 * etaEW;
                p["cbmunumu::Im{cSL}"]         = -3.0 * etaEW;
                p["cbmunumu::Re{cSR}"]         = +4.0 * etaEW;
                p["cbmunumu::Im{cSR}"]         = -4.0 * etaEW;
                p["cbmunumu::Re{cT}"]          = +5.0 * etaEW;
                p["cbmunumu::Im{cT}"]          = -5.0 * etaEW;
                // tau mode
                p["cbtaunutau::Re{cVL}"]       = +1.0 * etaEW;
                p["cbtaunutau::Im{cVL}"]       = -5.0 * etaEW;
                p["cbtaunutau::Re{cVR}"]       = +2.1 * etaEW;
                p["cbtaunutau::Im{cVR}"]       = -6.0 * etaEW;
                p["cbtaunutau::Re{cSL}"]       = +3.1 * etaEW;
                p["cbtaunutau::Im{cSL}"]       = -7.0 * etaEW;
                p["cbtaunutau::Re{cSR}"]       = +4.1 * etaEW;
                p["cbtaunutau::Im{cSR}"]       = -8.0 * etaEW;
                p["cbtaunutau::Re{cT}"]        = +5.1 * etaEW;
                p["cbtaunutau::Im{cT}"]        = -9.0 * etaEW;

                Options oo
                {
                    { "model",        "WET" },
                    { "form-factors", "BCL2008"    },
                    { "U",            "c"          },
                    { "q",            "d"          },
                    { "I",            "1/2"        },
                    { "l",            "mu"         }
                };

                const double eps = 1e-3;
                {
                    BToPseudoscalarLeptonNeutrino d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(0.011164, 11.62), 2615.77, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(0.011164, 11.62), -0.621944, eps);

                    auto k      = Kinematics{
                        { "q2_mu_min",   0.011164 },
                        { "q2_mu_max",  11.62     },
                        { "q2_tau_min",  3.15702  },
                        { "q2_tau_max", 11.62     }
                    };
                    oo =
                    {
                        { "model",        "WET" },
                        { "form-factors", "BCL2008"    },
                        { "U",            "c"          },
                        { "q",            "d"          },
                        { "I",            "1/2"        }
                    };
                    auto obs_RD = Observable::make("B->Dlnu::R_D", p, k, oo);

                    TEST_CHECK(obs_RD.get() != nullptr);
                    TEST_CHECK_RELATIVE_ERROR(obs_RD->evaluate(), 1.43554, eps);
                }
            }
        }
} b_to_d_l_nu_test;
