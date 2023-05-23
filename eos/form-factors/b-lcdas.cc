/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Philip Lüghausen
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

#include <eos/form-factors/b-lcdas-exponential.hh>
#include <eos/form-factors/b-lcdas-flvd2022.hh>
#include <eos/utils/exception.hh>

#include <functional>
#include <map>
#include <string>

namespace eos
{
    std::shared_ptr<BMesonLCDAs>
    BMesonLCDAs::make(const std::string & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When making an object for pseudoscalar LCDAs");

        std::shared_ptr<BMesonLCDAs> result;

        using type = std::function<BMesonLCDAs * (const Parameters &, const Options &)>;
        const std::map<std::string, type> models = {
            { "exponential",      &b_lcdas::Exponential::make      },
            { "FLvD2022",         &b_lcdas::FLvD2022::make         }
        };

        auto i = models.find(name);
        if (i != models.end())
        {
            result.reset( i->second(parameters, options) );
            return result;
        }

        throw InternalError("Unknown B-meson LCDA model: " + name);
        return result;
    }
}
