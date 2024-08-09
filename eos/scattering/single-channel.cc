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

#include <eos/scattering/scattering-amplitudes.hh>
#include <eos/scattering/parametric-gmkprdey2011.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/qualified-name.hh>

#include <cmath>
#include <limits>
#include <map>

namespace eos
{
    /* PP -> PP Processes */

    ScatteringAmplitudes<PPToPP>::~ScatteringAmplitudes()
    {
    }

    const std::map<ScatteringAmplitudeFactory<PPToPP>::KeyType, ScatteringAmplitudeFactory<PPToPP>::ValueType>
    ScatteringAmplitudeFactory<PPToPP>::scattering_amplitudes
    {
        { "pipi->pipi::GMKPRDEY2011",    &GMKPRDEY2011ScatteringAmplitudes::make    }
    };

    std::shared_ptr<ScatteringAmplitudes<PPToPP>>
    ScatteringAmplitudeFactory<PPToPP>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When creating a PP->PP scattering amplitude");

        std::shared_ptr<ScatteringAmplitudes<PPToPP>> result;

        auto i = scattering_amplitudes.find(name);
        if (scattering_amplitudes.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
            return result;
        }

        throw NoSuchScatteringAmplitudeError(name.prefix_part().str(), name.name_part().str());
        return result;
    }

    OptionSpecification
    ScatteringAmplitudeFactory<PPToPP>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "scattering-amplitudes", {}, "" };
        for (const auto & t : ScatteringAmplitudeFactory<PPToPP>::scattering_amplitudes)
        {
            if (process == std::get<0>(t).prefix_part())
                result.allowed_values.push_back(std::get<0>(t).name_part().str());
        }

        return result;
    }

    OptionSpecification
    ScatteringAmplitudeFactory<PPToPP>::option_specification()
    {
        std::set<std::string> allowed_values;
        for (const auto & t : ScatteringAmplitudeFactory<PPToPP>::scattering_amplitudes)
        {
            allowed_values.insert(std::get<0>(t).name_part().str());
        }

        OptionSpecification result { "scattering-amplitudes", { allowed_values.cbegin(), allowed_values.cend() }, "" };
        return result;
    }

}
