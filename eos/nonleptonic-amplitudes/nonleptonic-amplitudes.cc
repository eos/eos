/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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

#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes.hh>
#include <eos/nonleptonic-amplitudes/topological-amplitudes.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/options-impl.hh>

#include <map>

namespace eos
{
    NoSuchNonleptonicAmplitudeError::NoSuchNonleptonicAmplitudeError(const std::string & process, const std::string & tag) :
        Exception("No amplitude found for process '" + process + "' and tag '" + tag + "'!")
    {
    }

    NonleptonicAmplitudes<PToPP>::~NonleptonicAmplitudes() {};

    const std::map<NonleptonicAmplitudeFactory<PToPP>::KeyType, NonleptonicAmplitudeFactory<PToPP>::ValueType>
    NonleptonicAmplitudeFactory<PToPP>::amplitudes
    {
        { "B->PP::topological",         &TopologicalRepresentation<PToPP>::make},
        //{ "B->PP::su3_amplitudes",      &SU3Representation<PToPP>::make},
    };

    std::shared_ptr<NonleptonicAmplitudes<PToPP>>
    NonleptonicAmplitudeFactory<PToPP>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When creating a P->PP nonleptonic amplitude");

        std::shared_ptr<NonleptonicAmplitudes<PToPP>> result;

        auto i = NonleptonicAmplitudeFactory<PToPP>::amplitudes.find(name);
        if (NonleptonicAmplitudeFactory<PToPP>::amplitudes.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
            return result;
        }

        throw NoSuchNonleptonicAmplitudeError(name.prefix_part().str(), name.name_part().str());
        return result;
    }

    OptionSpecification
    NonleptonicAmplitudeFactory<PToPP>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "representation", {}, "" };
        for (const auto & ff : NonleptonicAmplitudeFactory<PToPP>::amplitudes)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    OptionSpecification
    NonleptonicAmplitudeFactory<PToPP>::option_specification()
    {
        std::set<std::string> allowed_values;
        for (const auto & ff : NonleptonicAmplitudeFactory<PToPP>::amplitudes)
        {
            allowed_values.insert(std::get<0>(ff).name_part().str());
        }

        OptionSpecification result { "representation", { allowed_values.cbegin(), allowed_values.cend() }, "" };
        return result;
    }

}
