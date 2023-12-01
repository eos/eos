/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Stefan Meiser
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

#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/vec-lcdas.hh>

#include <map>

namespace eos
{
    VectorLCDAs::~VectorLCDAs()
    {
    }

    std::shared_ptr<VectorLCDAs>
    VectorLCDAs::make(const std::string & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When making an object for vector LCDAs");

        using KeyType = std::string;
        using ValueType = std::function<VectorLCDAs * (const Parameters &, const Options &)>;
        static const std::map<KeyType, ValueType> lcdas
        {
        };

        std::shared_ptr<VectorLCDAs> result;
        auto i = lcdas.find(name);
        if (lcdas.cend() != i)
        {
            result.reset(i->second(parameters, options));
            return result;
        }

        throw InternalError("Unknown vector LCDAs for state: " + name);
        return result;
    }
}
