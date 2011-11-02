/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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
#include <eos/utils/random_number_generator.hh>

#include <cmath>
#include <random>

#include <iostream>

using namespace test;
using namespace eos;

class RandomNumberGeneratorTest :
    public TestCase
{
    public:
        RandomNumberGeneratorTest() :
            TestCase("rng_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 3e-15;

            // Uniform [0, 1)
            {
                RandomNumberGenerator rng(1723);

                TEST_CHECK_RELATIVE_ERROR(0.755526696098968390, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.640330279245972630, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.212258085142821070, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.956574363866820930, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.512321577174589040, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.137894445098936560, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.733577476348727940, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.648340581450611350, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.512517530936747790, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.710519498679786920, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.100748437456786630, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.036182452691718936, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.793245769571512940, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.316090840846300130, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.910442729713395240, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.137844955082982780, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.863410061690956350, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.640690742991864680, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.414283346850425000, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.541501202620565890, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.354803816881030800, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.084285020828247070, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.098871880210936069, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.709438384976238010, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.273271531565114860, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.461453695315867660, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.750975034898146990, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.589485029224306340, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.351696515223011370, rng(), eps);
                TEST_CHECK_RELATIVE_ERROR(0.993107097456231710, rng(), eps);
            }

            // Uniform [-0.321, 1.234)
            {
                static const double min = -0.321, max = 1.234;
                RandomNumberGenerator rng(1723);
                #if __GNUC__ >= 4 && __GNUC_MINOR__ < 5
                std::uniform_real<double> dist(min, max);
                #elif __GNUC__ >= 4 && __GNUC_MINOR__ >= 5
                std::uniform_real_distribution<double> dist(min, max);
                #endif
                TEST_CHECK_RELATIVE_ERROR(+0.8538440124338957400, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.6747135842274874000, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.0090613223970867351, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.1664731358129066000, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.4756600525064859600, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(-0.1065741378711536600, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.8197129757222718600, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.6871696041557007200, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.4759647606066427800, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.7838578204470687500, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(-0.1643361797546968000, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(-0.2647362860643770500, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.9124971716837027000, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.1705212575159966700, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.0947384447043296000, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(-0.1066510948459618000, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.0216026459294372000, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.6752741053523494400, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.3232106043524108400, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.5210343700749799900, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.2307199352500028400, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(-0.1899367926120758300, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(-0.1672542262719944300, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.7821766886380501200, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.1039372315837535600, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.3965604962161742300, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.8467661792666185700, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.5956492204437964200, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+0.2258880811717826600, dist(rng), eps);
                TEST_CHECK_RELATIVE_ERROR(+1.2232815365444403000, dist(rng), eps);
            }
        }
} rng_test;
