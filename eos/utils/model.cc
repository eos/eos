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

#include <eos/utils/ckm_scan_model.hh>
#include <eos/utils/model.hh>
#include <eos/utils/standard-model.hh>
#include <eos/utils/wilson_scan_model.hh>

#include <map>

namespace eos
{
    Model::~Model()
    {
    }

    std::shared_ptr<Model>
    Model::make(const std::string & name, const Parameters & parameters, const Options & options)
    {
        typedef std::function<std::shared_ptr<Model> (const Parameters &, const Options &)> ModelMaker;
        static const std::map<std::string, ModelMaker> model_makers
        {
            std::make_pair("CKMScan", &CKMScanModel::make),
            std::make_pair("SM", &StandardModel::make),
            std::make_pair("WilsonScan", &WilsonScanModel::make),
        };

        auto i = model_makers.find(name);

        if (model_makers.cend() == i)
            throw NoSuchModelError(name);

        return i->second(parameters, options);
    }

    NoSuchModelError::NoSuchModelError(const std::string & name) :
        Exception("No such model: '" + name + "'")
    {
    }
}
