/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024-2025 Florian Herren
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

                p["pipi->pipi::S0_B_0@GMKPRDEY2011"]    = 7.26;
                p["pipi->pipi::S0_B_1@GMKPRDEY2011"]    = -25.3;
                p["pipi->pipi::S0_B_2@GMKPRDEY2011"]    = -33.1;
                p["pipi->pipi::S0_B_3@GMKPRDEY2011"]    = -26.6;
                p["pipi->pipi::S0_d_0@GMKPRDEY2011"]    = 3.96364;
                p["pipi->pipi::S0_c@GMKPRDEY2011"]      = -11.5192;
                p["pipi->pipi::S0_B@GMKPRDEY2011"]      = 1.64061;
                p["pipi->pipi::S0_C@GMKPRDEY2011"]      = 0.705113;
                p["pipi->pipi::S0_D@GMKPRDEY2011"]      = -1.51669;

                p["pipi->pipi::P1_B_0@GMKPRDEY2011"]    = 1.055;
                p["pipi->pipi::P1_B_1@GMKPRDEY2011"]    = 0.15;
                p["pipi->pipi::P1_lam_1@GMKPRDEY2011"]  = 1.57;
                p["pipi->pipi::P1_lam_2@GMKPRDEY2011"]  = -1.96;

                p["pipi->pipi::D0_B_0@GMKPRDEY2011"]    = 12.47;
                p["pipi->pipi::D0_B_1@GMKPRDEY2011"]    = 10.12;
                p["pipi->pipi::D0_Bh_1@GMKPRDEY2011"]   = 43.7;

                p["mass::pi^+@GMKPRDEY2011"]            = 0.13957;
                p["mass::eta@GMKPRDEY2011"]             = 0.54751;
                p["mass::K_u@GMKPRDEY2011"]             = 0.496;
                p["mass::rho^0@GMKPRDEY2011"]           = 0.7736;
                p["mass::f_2@GMKPRDEY2011"]             = 1.2754;

                p["pipi->pipi::S0_sM@GMKPRDEY2011"]     = 0.7225;
                p["pipi->pipi::P1_s0@GMKPRDEY2011"]     = 1.1025;
                p["pipi->pipi::D0_s0@GMKPRDEY2011"]     = 1.1025;
                p["pipi->pipi::D0_sh@GMKPRDEY2011"]     = 2.1025;

                p["pipi->pipi::S0_n@GMKPRDEY2011"]      = 8.0;
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

                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 0 , IsospinRepresentation::zero).real() , -0.449485 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 0 , IsospinRepresentation::zero).imag() ,  0.790851 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 1 , IsospinRepresentation::one).real()  , -0.498251 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 1 , IsospinRepresentation::one).imag()  ,  0.290285 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 2 , IsospinRepresentation::zero).real() ,  0.215659 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 2 , IsospinRepresentation::zero).imag() ,  0.038253 ,  eps);

                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.55  , 1 , IsospinRepresentation::one).real()    ,  2.542132 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.55  , 1 , IsospinRepresentation::one).imag()    ,  5.375608 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.9   , 1 , IsospinRepresentation::one).real()    , -1.760142 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.9   , 1 , IsospinRepresentation::one).imag()    ,  1.025474 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(9.0   , 1 , IsospinRepresentation::one).real()    , -0.081634 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(9.0   , 1 , IsospinRepresentation::one).imag()    ,  0.005939 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(100.0 , 1 , IsospinRepresentation::one).real()    , -0.006934 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(100.0 , 1 , IsospinRepresentation::one).imag()    ,  0.000104 ,  eps);

                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.55  , 2 , IsospinRepresentation::zero).real()   ,  1.611722 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.55  , 2 , IsospinRepresentation::zero).imag()   ,  0.079599 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(1.44  , 2 , IsospinRepresentation::zero).real()   ,  4.823093 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(1.44  , 2 , IsospinRepresentation::zero).imag()   ,  5.003143 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(1.65  , 2 , IsospinRepresentation::zero).real()   , -0.917640 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(1.65  , 2 , IsospinRepresentation::zero).imag()   ,  8.643233 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(100.0 , 2 , IsospinRepresentation::zero).real()   , -0.016038 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(100.0 , 2 , IsospinRepresentation::zero).imag()   ,  0.0      ,  eps);

                double sp = 4.0 * 0.496 * 0.496;
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_outer_function(-100000.0 , sp , 0.0 , 1e-5, 1 , IsospinRepresentation::one).real(),  0.0      ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_outer_function(      0.0 , sp , 0.0 , 1e-5, 1 , IsospinRepresentation::one).real(),  0.26187  ,  eps);

                p["pipi->pipi::P1_B_0@GMKPRDEY2011"]    = 1.066;
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 1 , IsospinRepresentation::one).real() , -0.495833 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 1 , IsospinRepresentation::one).imag() ,  0.286062 ,  eps);
            }
        }
} pi_pi_to_pi_pi_gmkprdey2011_scattering_amplitudes_test;
