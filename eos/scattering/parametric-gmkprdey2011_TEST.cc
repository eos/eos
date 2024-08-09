/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Florian Herren
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
#include <eos/scattering/parametric-gmkprdey2011.hh>

using namespace test;
using namespace eos;

class GMKPRDEY2011ScatteringAmplitudesTest :
    public TestCase
{
    public:
        GMKPRDEY2011ScatteringAmplitudesTest() :
            TestCase("pi_pi_to_pi_pi_gmkprdey2011_scattering_amplitudes_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;
            /* Factory */
            {
                Parameters p = Parameters::Defaults();
                std::shared_ptr<ScatteringAmplitudes<PPToPP>> amp = ScatteringAmplitudeFactory<PPToPP>::create("pipi->pipi::GMKPRDEY2011", p, Options{ });

                TEST_CHECK(0 != amp.get());
            }

            {
                Parameters p = Parameters::Defaults();
                p["pipi->pipi::P1_B_0@GMKPRDEY2011"]    = 1.055;
                p["pipi->pipi::P1_B_1@GMKPRDEY2011"]    = 0.15;
                p["pipi->pipi::P1_lam_1@GMKPRDEY2011"]  = 1.57;
                p["pipi->pipi::P1_lam_2@GMKPRDEY2011"]  = -1.96;
                
                p["mass::pi^+@GMKPRDEY2011"]            = 0.13957;
                p["mass::eta@GMKPRDEY2011"]             = 0.54751;
                p["mass::K_u@GMKPRDEY2011"]             = 0.496;
                p["mass::rho^0@GMKPRDEY2011"]           = 0.7736;
                p["mass::f_2@GMKPRDEY2011"]             = 1.2754;

                p["pipi->pipi::P1_s0@GMKPRDEY2011"]     = 1.1025;
                p["pipi->pipi::D0_s0@GMKPRDEY2011"]     = 1.1025;
                p["pipi->pipi::D0_sh@GMKPRDEY2011"]     = 2.1025;

                p["pipi->pipi::P1_n@GMKPRDEY2011"]      = 0.75;
                p["pipi->pipi::D0_n@GMKPRDEY2011"]      = 10.0;

                GMKPRDEY2011ScatteringAmplitudes amp(p, Options{ });
                Diagnostics diagnostics = amp.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  -1.0,      eps),       // w(s =  0)
                    std::make_pair(  0.514972,  eps),       // w(s =  1)

                    std::make_pair(  0.690453,  eps),       // del_S0(s =  0.25)
                    std::make_pair(  1.605099,  eps),       // del_S0(s =  0.72)
                    std::make_pair(  2.087624,  eps),       // del_S0(s =  0.9)
                    std::make_pair(  4.570380,  eps),       // del_S0(s =  1.44)
                    std::make_pair(  6.276073,  eps),       // del_S0(s =  4.0)

                    std::make_pair(  0.099257,  eps),       // del_P1(s =  0.25)
                    std::make_pair(  2.614059,  eps),       // del_P1(s =  0.9)
                    std::make_pair(  2.693794,  eps),       // del_P1(s =  1.0)
                    std::make_pair(  3.030938,  eps),       // del_P1(s =  4.0)

                    std::make_pair(  0.006052,  eps),       // del_D0(s =  0.25)
                    std::make_pair(  0.175553,  eps),       // del_D0(s =  0.9)
                    std::make_pair(  0.803720,  eps),       // del_D0(s =  1.44)
                    std::make_pair(  3.140443,  eps),       // del_D0(s =  4.0)
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 0 , "0").real() , -0.449485 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 0 , "0").imag() ,  0.790851 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 1 , "1").real() , -0.498251 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 1 , "1").imag() ,  0.290285 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 2 , "0").real() ,  0.215659 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 2 , "0").imag() ,  0.038253 ,  eps);

                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.55  , 1 , "1").real()       ,  2.542132 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.55  , 1 , "1").imag()       ,  5.375608 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.9   , 1 , "1").real()       , -1.760142 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.9   , 1 , "1").imag()       ,  1.025474 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(9.0   , 1 , "1").real()       , -0.081634 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(9.0   , 1 , "1").imag()       ,  0.005939 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(100.0 , 1 , "1").real()       , -0.006934 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(100.0 , 1 , "1").imag()       ,  0.000104 ,  eps);
                
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.55  , 2 , "0").real()       ,  1.611722 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.55  , 2 , "0").imag()       ,  0.079599 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(1.44  , 2 , "0").real()       ,  4.823117 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(1.44  , 2 , "0").imag()       ,  5.003168 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(1.65   , 2 , "0").real()      , -0.917640 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(1.65   , 2 , "0").imag()      ,  8.643233 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(100.0  , 2 , "0").real()      , -0.016038 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(100.0  , 2 , "0").imag()      ,  0.0      ,  eps);
            }
        }
} pi_pi_to_pi_pi_gmkprdey2011_scattering_amplitudes_test;
