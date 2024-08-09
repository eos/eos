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
#include <eos/scattering/scattering-amplitudes.hh>
#include <eos/scattering/single-channel.hh>

#include <eos/models/model.hh>
#include <eos/maths/power-of.hh>

using namespace test;
using namespace eos;

class PPToPPScatteringAmplitudesTest :
    public TestCase
{
    public:
        PPToPPScatteringAmplitudesTest() :
            TestCase("pp_to_pp_scattering_amplitude_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchScatteringAmplitudeError, ScatteringAmplitudeFactory<PPToPP>::create("Foo->Bar::ABC2024", parameter, options));
            }
        }
} pp_to_pp_scattering_amplitude_test;
