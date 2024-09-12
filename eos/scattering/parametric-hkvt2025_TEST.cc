/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025, Florian Herren
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
#include <eos/scattering/parametric-hkvt2025.hh>

using namespace test;
using namespace eos;

class HKVT2025ScatteringAmplitudesTest :
    public TestCase
{
    public:
        HKVT2025ScatteringAmplitudesTest() :
            TestCase("pi_pi_to_pi_pi_hkvt2025_scattering_amplitudes_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;
            /* Factory */
            {
                Parameters p = Parameters::Defaults();
                std::shared_ptr<ScatteringAmplitudes<PPToPP>> amp = ScatteringAmplitudeFactory<PPToPP>::create("pipi->pipi::HKvT2025", p, Options{ });

                TEST_CHECK(0 != amp.get());
            }

            {
                Parameters p = Parameters::Defaults();

                p["pipi->pipi::D0_B_0@GMKPRDEY2011"]    = 12.47;
                p["pipi->pipi::D0_B_1@GMKPRDEY2011"]    = 10.12;
                p["pipi->pipi::D0_Bh_1@GMKPRDEY2011"]   = 43.7;
                p["mass::pi^+@GMKPRDEY2011"]            = 0.13957;
                p["mass::f_2@GMKPRDEY2011"]             = 1.2754;
                p["pipi->pipi::D0_s0@GMKPRDEY2011"]     = 1.1025;
                p["pipi->pipi::D0_sh@GMKPRDEY2011"]     = 2.1025;

                p["pipi->pipi::P1_n@GMKPRDEY2011"]      = 0.75;
                p["pipi->pipi::D0_n@GMKPRDEY2011"]      = 10.0;

                p["pipi->pipi::Gamman0_pi@HKvT2025"]    = 0.984;
                p["pipi->pipi::Gamman0_K@HKvT2025"]     = 0.5;

                HKVT2025ScatteringAmplitudes amp(p, Options{ });
                Diagnostics diagnostics = amp.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  -1.0,      eps),       // w(s =  0)
                    std::make_pair(  0.514972,  eps),       // w(s =  1)

                    std::make_pair(  0.103826,  eps),       // del_P1(s =  0.25)
                    std::make_pair(  2.652486,  eps),       // del_P1(s =  0.9)
                    std::make_pair(  2.926285,  eps),       // del_P1(s =  1.44)
                    std::make_pair(  3.035405,  eps),       // del_P1(s =  4.0)

                    std::make_pair(  0.006052,  eps),       // del_D0(s =  0.25)
                    std::make_pair(  0.175553,  eps),       // del_D0(s =  0.9)
                    std::make_pair(  0.803720,  eps),       // del_D0(s =  1.44)
                    std::make_pair(  3.140443,  eps),       // del_D0(s =  4.0)
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 2 , IsospinRepresentation::zero).real() ,  0.215659 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.scattering_amplitude(0.9 , 2 , IsospinRepresentation::zero).imag() ,  0.038253 ,  eps);

                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.700332 , 0 , IsospinRepresentation::zero).real()   ,  0.0646509 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.700332 , 0 , IsospinRepresentation::zero).imag()   ,  1.2225023 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(25.0     , 0 , IsospinRepresentation::zero).real()   , -0.0375782 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(25.0     , 0 , IsospinRepresentation::zero).imag()   ,  0.0017212 ,  eps);

                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.700332 , 1 , IsospinRepresentation::one).real()   , -2.5211430 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.700332 , 1 , IsospinRepresentation::one).imag()   ,  3.1724329 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(25.0     , 1 , IsospinRepresentation::one).real()   , -0.0277992 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(25.0     , 1 , IsospinRepresentation::one).imag()   ,  0.0010054 ,  eps);

                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.55  , 2 , IsospinRepresentation::zero).real()   ,  1.611722 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(0.55  , 2 , IsospinRepresentation::zero).imag()   ,  0.079599 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(1.44  , 2 , IsospinRepresentation::zero).real()   ,  4.823093 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(1.44  , 2 , IsospinRepresentation::zero).imag()   ,  5.003143 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(100.0 , 2 , IsospinRepresentation::zero).real()   , -0.016038 ,  eps);
                TEST_CHECK_NEARLY_EQUAL( amp.omnes_factor(100.0 , 2 , IsospinRepresentation::zero).imag()   ,  0.0      ,  eps);

            }
        }
} pi_pi_to_pi_pi_hkvt2025_scattering_amplitudes_test;
