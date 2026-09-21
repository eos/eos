/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#include <eos/maths/complex.hh>
#include <eos/observable.hh>
#include <eos/rare-b-decays/b-to-kstar-gamma.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class BToKstarGammaNaiveTest : public TestCase
{
    public:
        BToKstarGammaNaiveTest() :
            TestCase("b_to_kstar_gamma_naive_test")
        {
        }

        virtual void
        run() const
        {
            // Standard Model
            {
                Parameters p                  = Parameters::Defaults();
                p["b->s::c1"]                 = -0.32300000;
                p["b->s::c2"]                 = +1.00931000;
                p["b->s::c3"]                 = -0.00522869;
                p["b->s::c4"]                 = -0.08794730;
                p["b->s::c5"]                 = +0.00037476;
                p["b->s::c6"]                 = +0.00105859;
                p["b->s::Re{c7}"]             = -0.331;
                p["b->s::Re{c7'}"]            = -0.00659; // m_s(m_b) / m_b(m_b) * Abs{c7} = 85 / 4200 * Abs{c7}
                p["b->s::c8"]                 = -0.181;
                // PDG 2010 CKM parameters
                p["CKM::A"]                   = 0.812;
                p["CKM::lambda"]              = 0.22543;
                p["CKM::rhobar"]              = 0.144;
                p["CKM::etabar"]              = 0.342;
                p["CKM::abs(V_ub)"]           = 0.003540950873054711;
                p["CKM::arg(V_ub)"]           = -1.1728563751359748;
                p["CKM::abs(V_cb)"]           = 0.04126451344307112;
                p["CKM::arg(V_cb)"]           = 0.0;
                p["CKM::abs(V_tb)"]           = 0.9991419776905534;
                p["CKM::arg(V_tb)"]           = 0.0;
                p["CKM::abs(V_td)"]           = 0.008576901910577167;
                p["CKM::arg(V_td)"]           = -0.37951557931964897;
                p["CKM::abs(V_us)"]           = 0.22542858674178629;
                p["CKM::arg(V_us)"]           = 0.0;
                p["CKM::abs(V_cs)"]           = 0.9734167680132911;
                p["CKM::arg(V_cs)"]           = -3.119448393424795e-05;
                p["CKM::abs(V_ts)"]           = 0.04051834255894421;
                p["CKM::arg(V_ts)"]           = -3.123445879630718;
                p["mass::b(MSbar)"]           = 4.2;
                p["mass::B_d"]                = 5.27958;
                p["mass::K_d^*"]              = 0.89594;
                // the KMPW2010 parameters default to zero; use the central values of [KMPW:2010A], Table 4, p. 31
                p["B->K^*::F^T1(0)@KMPW2010"] = +0.31;
                p["B->K^*::b^T1_1@KMPW2010"]  = -4.6;

                Options oo{
                    {        "model"_ok,      "WET"_ov },
                    {          "tag"_ok,    "Naive"_ov },
                    { "form-factors"_ok, "KMPW2010"_ov }
                };

                BToKstarGamma d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(), +2.58953e-5, eps);

                TEST_CHECK_RELATIVE_ERROR(Observable::make("B->K^*gamma::BR", p, Kinematics(), oo)->evaluate(), +2.58953e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(Observable::make("B->K^*gamma::S_K^*gamma", p, Kinematics(), oo)->evaluate(), -3.00135e-2, eps);
                // the amplitudes carry no strong phase at this order
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->K^*gamma::C_K^*gamma", p, Kinematics(), oo)->evaluate(), 0.0, 1.0e-15);

                // the right-handed amplitude is driven by C_7' alone
                p["b->s::Re{c7'}"] = 0.0;

                BToKstarGamma e(p, oo);

                TEST_CHECK_NEARLY_EQUAL(e.real_a_right(), 0.0, 1.0e-15);
                TEST_CHECK_NEARLY_EQUAL(e.imag_a_right(), 0.0, 1.0e-15);
            }
        }
} b_to_kstar_gamma_naive_test;
