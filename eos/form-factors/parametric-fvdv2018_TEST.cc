/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2022 Danny van Dyk
 * Copyright (c) 2018 Keri Vos
 * Copyright (c) 2025 Florian Herren
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
#include <eos/form-factors/parametric-fvdv2018.hh>

#include <vector>

using namespace test;
using namespace eos;

class BToPiPiFvDV2018FormFactorsTest :
    public TestCase
{
    public:
        BToPiPiFvDV2018FormFactorsTest() :
        TestCase("b_to_pi_pi_fvdv2018_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 5.1e-3;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToPP>> ff = FormFactorFactory<PToPP>::create("B->pipi::FvDV2018", p, Options{ });
            TEST_CHECK(ff.get() != nullptr);

            // time
            p["B->pipi::a^Ftime_0_0@FvDV2018"] = -0.36816929771;
            p["B->pipi::a^Ftime_0_1@FvDV2018"] =  2.24240299735;
            p["B->pipi::a^Ftime_0_2@FvDV2018"] = -4.24917695998;
            p["B->pipi::a^Ftime_0_3@FvDV2018"] =  2.35541703192;
            p["B->pipi::a^Ftime_1_0@FvDV2018"] =  5.24115374042;
            p["B->pipi::a^Ftime_1_1@FvDV2018"] = -22.0437225418;
            p["B->pipi::a^Ftime_1_2@FvDV2018"] =  20.2542056471;
            p["B->pipi::b^Ftime_0_0@FvDV2018"] = -0.78770063503;
            p["B->pipi::b^Ftime_0_1@FvDV2018"] =  10.2021697105;
            p["B->pipi::b^Ftime_0_2@FvDV2018"] = -36.9473507702;
            p["B->pipi::b^Ftime_0_3@FvDV2018"] =  41.0909785695;
            p["B->pipi::b^Ftime_1_0@FvDV2018"] =  12.1068819568;
            p["B->pipi::b^Ftime_1_1@FvDV2018"] = -86.2325621190;
            p["B->pipi::b^Ftime_1_2@FvDV2018"] =  147.922534873;
            p["B->pipi::c^Ftime_0_0@FvDV2018"] =  1.56132165878;
            p["B->pipi::c^Ftime_0_1@FvDV2018"] = -15.8641613048;
            p["B->pipi::c^Ftime_0_2@FvDV2018"] =  51.0904259777;
            p["B->pipi::c^Ftime_0_3@FvDV2018"] = -52.9767513762;
            p["B->pipi::c^Ftime_1_0@FvDV2018"] = -23.6236882937;
            p["B->pipi::c^Ftime_1_1@FvDV2018"] =  141.591407772;
            p["B->pipi::c^Ftime_1_2@FvDV2018"] = -209.769464670;
            // long
            p["B->pipi::a^Flong_0_0@FvDV2018"] = -0.16139192740;
            p["B->pipi::a^Flong_0_1@FvDV2018"] =  1.09519092754;
            p["B->pipi::a^Flong_0_2@FvDV2018"] = -2.41027721486;
            p["B->pipi::a^Flong_0_3@FvDV2018"] =  1.70033143492;
            p["B->pipi::a^Flong_1_0@FvDV2018"] =  1.99824594658;
            p["B->pipi::a^Flong_1_1@FvDV2018"] = -8.62878627878;
            p["B->pipi::a^Flong_1_2@FvDV2018"] =  9.35509122572;
            p["B->pipi::b^Flong_0_0@FvDV2018"] = -0.44278949412;
            p["B->pipi::b^Flong_0_1@FvDV2018"] =  4.74822500206;
            p["B->pipi::b^Flong_0_2@FvDV2018"] = -15.9081401854;
            p["B->pipi::b^Flong_0_3@FvDV2018"] =  17.0309192822;
            p["B->pipi::b^Flong_1_0@FvDV2018"] = -5.92715117703;
            p["B->pipi::b^Flong_1_1@FvDV2018"] =  32.8787517243;
            p["B->pipi::b^Flong_1_2@FvDV2018"] = -44.5974162489;
            p["B->pipi::c^Flong_0_0@FvDV2018"] =  0.77102061709;
            p["B->pipi::c^Flong_0_1@FvDV2018"] = -6.98963227925;
            p["B->pipi::c^Flong_0_2@FvDV2018"] =  21.1956324818;
            p["B->pipi::c^Flong_0_3@FvDV2018"] = -21.2213423789;
            p["B->pipi::c^Flong_1_0@FvDV2018"] =  4.11524737745;
            p["B->pipi::c^Flong_1_1@FvDV2018"] = -23.9234305792;
            p["B->pipi::c^Flong_1_2@FvDV2018"] =  32.8488722373;
            // para
            p["B->pipi::a^Fpara_0_0@FvDV2018"] = -0.73839711162;
            p["B->pipi::a^Fpara_0_1@FvDV2018"] =  5.33340671033;
            p["B->pipi::a^Fpara_0_2@FvDV2018"] = -13.0548108152;
            p["B->pipi::a^Fpara_0_3@FvDV2018"] =  10.7823778312;
            p["B->pipi::a^Fpara_1_0@FvDV2018"] =  4.90990913468;
            p["B->pipi::a^Fpara_1_1@FvDV2018"] =  3.07689657262;
            p["B->pipi::a^Fpara_1_2@FvDV2018"] = -51.8342261186;
            p["B->pipi::b^Fpara_0_0@FvDV2018"] = -0.77364251324;
            p["B->pipi::b^Fpara_0_1@FvDV2018"] =  8.65261164441;
            p["B->pipi::b^Fpara_0_2@FvDV2018"] = -29.0610911551;
            p["B->pipi::b^Fpara_0_3@FvDV2018"] =  30.4513866623;
            p["B->pipi::b^Fpara_1_0@FvDV2018"] =  16.8847387994;
            p["B->pipi::b^Fpara_1_1@FvDV2018"] = -203.754603759;
            p["B->pipi::b^Fpara_1_2@FvDV2018"] =  456.700391055;
            p["B->pipi::c^Fpara_0_0@FvDV2018"] =  1.91667289560;
            p["B->pipi::c^Fpara_0_1@FvDV2018"] = -16.8682304272;
            p["B->pipi::c^Fpara_0_2@FvDV2018"] =  48.4938486882;
            p["B->pipi::c^Fpara_0_3@FvDV2018"] = -45.0828744220;
            p["B->pipi::c^Fpara_1_0@FvDV2018"] = -32.1123316945;
            p["B->pipi::c^Fpara_1_1@FvDV2018"] =  273.753395427;
            p["B->pipi::c^Fpara_1_2@FvDV2018"] = -531.002613990;
            // perp
            p["B->pipi::a^Fperp_0_0@FvDV2018"] = -2.05535479854;
            p["B->pipi::a^Fperp_0_1@FvDV2018"] =  16.9988604895;
            p["B->pipi::a^Fperp_0_2@FvDV2018"] = -46.6407789359;
            p["B->pipi::a^Fperp_0_3@FvDV2018"] =  42.5157927594;
            p["B->pipi::a^Fperp_1_0@FvDV2018"] =  25.4079749879;
            p["B->pipi::a^Fperp_1_1@FvDV2018"] = -139.728111016;
            p["B->pipi::a^Fperp_1_2@FvDV2018"] =  193.460283599;
            p["B->pipi::b^Fperp_0_0@FvDV2018"] =  8.23644248367;
            p["B->pipi::b^Fperp_0_1@FvDV2018"] = -67.3709309476;
            p["B->pipi::b^Fperp_0_2@FvDV2018"] =  183.591398737;
            p["B->pipi::b^Fperp_0_3@FvDV2018"] = -166.684702751;
            p["B->pipi::b^Fperp_1_0@FvDV2018"] = -112.935412057;
            p["B->pipi::b^Fperp_1_1@FvDV2018"] =  624.097066875;
            p["B->pipi::b^Fperp_1_2@FvDV2018"] = -867.134281277;
            p["B->pipi::c^Fperp_0_0@FvDV2018"] = -6.96163577513;
            p["B->pipi::c^Fperp_0_1@FvDV2018"] =  56.0519566307;
            p["B->pipi::c^Fperp_0_2@FvDV2018"] = -151.061706725;
            p["B->pipi::c^Fperp_0_3@FvDV2018"] =  136.078387363;
            p["B->pipi::c^Fperp_1_0@FvDV2018"] =  106.107678887;
            p["B->pipi::c^Fperp_1_1@FvDV2018"] = -588.669971295;
            p["B->pipi::c^Fperp_1_2@FvDV2018"] =  819.519985258;

            // time
            TEST_CHECK_NEARLY_EQUAL(real(ff->f_time(0.60, 19.0, +0.5)),  0.0,      eps);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_time(0.60, 19.0, +0.5)),  0.068642, eps);
            TEST_CHECK_NEARLY_EQUAL(real(ff->f_time(0.05, 16.0, -0.5)),  0.0,      eps);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_time(0.05, 16.0, -0.5)),  0.433295, eps);
            // long
            TEST_CHECK_NEARLY_EQUAL(real(ff->f_long(0.60, 19.0, +0.5)),  0.0,      eps);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_long(0.60, 19.0, +0.5)),  0.032153, eps);
            TEST_CHECK_NEARLY_EQUAL(real(ff->f_long(0.05, 16.0, -0.5)),  0.0,      eps);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_long(0.05, 16.0, -0.5)),  0.413832, eps);
            // para
            TEST_CHECK_NEARLY_EQUAL(real(ff->f_para(0.60, 19.0, +0.5)),  0.0,      eps);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_para(0.60, 19.0, +0.5)), -0.021664, eps);
            TEST_CHECK_NEARLY_EQUAL(real(ff->f_para(0.05, 16.0, -0.5)),  0.0,      eps);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_para(0.05, 16.0, -0.5)), -0.013560, eps);
            // perp
            TEST_CHECK_NEARLY_EQUAL(real(ff->f_perp(0.60, 19.0, +0.5)), 0.0,      eps);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_perp(0.60, 19.0, +0.5)), 0.001048, eps);
            TEST_CHECK_NEARLY_EQUAL(real(ff->f_perp(0.05, 16.0, -0.5)), 0.0,      eps);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_perp(0.05, 16.0, -0.5)), 0.001215, eps);

            // partial waves
            auto time_pw = ff->f_time(0.60, 19.0);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_time(0.60, 19.0, 1.0)) - imag(time_pw[0] + std::sqrt(3.0) * time_pw[1] + std::sqrt(5.0) * time_pw[2] + std::sqrt(7.0) * time_pw[3]), 0.0, eps);

            auto long_pw = ff->f_long(0.60, 19.0);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_long(0.60, 19.0, 1.0)) - imag(long_pw[0] + std::sqrt(3.0) * long_pw[1] + std::sqrt(5.0) * long_pw[2] + std::sqrt(7.0) * long_pw[3]), 0.0, eps);

            auto perp_pw = ff->f_perp(0.60, 19.0);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_perp(0.60, 19.0, 0.0)) - imag(std::sqrt(3.0) * perp_pw[1] - 1.5 * std::sqrt(7.0) * perp_pw[3]), 0.0, eps);

            auto para_pw = ff->f_para(0.60, 19.0);
            TEST_CHECK_NEARLY_EQUAL(imag(ff->f_para(0.60, 19.0, 0.0)) - imag(std::sqrt(3.0) * para_pw[1] - 1.5 * std::sqrt(7.0) * para_pw[3]), 0.0, eps);
        }
} b_to_pi_pi_fvdv2018_form_factors_test;
