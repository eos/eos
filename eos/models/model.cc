/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2014 Danny van Dyk
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

#include <eos/models/ckm.hh>
#include <eos/models/model.hh>
#include <eos/models/standard-model.hh>
#include <eos/models/wet.hh>

#include <map>

namespace eos
{
    Model::~Model() {}

    const std::map<Model::KeyType, Model::ValueType> Model::models{
        std::make_pair("CKM", &CKMScanModel::make),
        std::make_pair("SM", &StandardModel::make),
        std::make_pair("WET", &WilsonScanModel::make),
        std::make_pair("WET-SMEFT", &ConstrainedWilsonScanModel::make),
    };

    std::shared_ptr<Model>
    Model::make(const std::string & name, const Parameters & parameters, const Options & options)
    {
        auto i = Model::models.find(name);

        if (Model::models.cend() == i)
        {
            throw NoSuchModelError(name);
        }

        return i->second(parameters, options);
    }

    OptionSpecification
    Model::option_specification()
    {
        OptionSpecification result{ "model", {}, "SM" };

        for (const auto & m : Model::models)
        {
            result.allowed_values.push_back(std::get<0>(m));
        }

        return result;
    }

    NoSuchModelError::NoSuchModelError(const std::string & name) :
        Exception("No such model: '" + name + "'")
    {
    }
} // namespace eos
