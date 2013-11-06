/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Danny van Dyk
 * Copyright (c) 2013 Frederik Beaujean
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

#include <eos/statistics/density.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

namespace eos
{
    template class WrappedForwardIterator<Density::IteratorTag, const ParameterDescription>;

    Density::~Density()
    {
    }

    void
    Density::dump_descriptions(hdf5::File & file, const std::string & data_set_base) const
    {
        auto data_set = file.create_data_set(data_set_base + "/parameters", Density::Output::description_type());

        auto record = Density::Output::description_record();

        // loop over descriptions
        for (auto & d : *this)
        {
            std::get<0>(record) = d.parameter->name().c_str();
            std::get<1>(record) = d.min;
            std::get<2>(record) = d.max;
            std::get<3>(record) = int(d.nuisance);

            data_set << record;
        }
    }

    Density::Output::DescriptionType
    Density::Output::description_type()
    {
        return
        DescriptionType
        {
         "parameter description",
         hdf5::Scalar<const char *>("name"),
         hdf5::Scalar<double>("min"),
         hdf5::Scalar<double>("max"),
         hdf5::Scalar<int>("nuisance"),
        };
    }

    std::tuple<const char *, double, double, int>
    Density::Output::description_record()
    {
        return std::make_tuple("name", 1.0, 2.0, 3);
    }
}
