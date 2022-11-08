/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
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

#include <eos/form-factors/psd-lcdas.hh>
#include <eos/form-factors/pi-lcdas.hh>
#include <eos/form-factors/k-lcdas.hh>

#include <map>

namespace eos
{
    PseudoscalarLCDAs::~PseudoscalarLCDAs()
    {
    }

    std::shared_ptr<PseudoscalarLCDAs>
    PseudoscalarLCDAs::make(const std::string & name, const Parameters & parameters, const Options & options)
    {
        using KeyType = std::string;
        using ValueType = std::function<PseudoscalarLCDAs * (const Parameters &, const Options &)>;
        static const std::map<KeyType, ValueType> lcdas
        {
            { "pi",   &PionLCDAs::make     },
            { "K",    &KaonLCDAs::make     },
            { "Kbar", &AntiKaonLCDAs::make }
        };

        std::shared_ptr<PseudoscalarLCDAs> result;
        auto i = lcdas.find(name);
        if (lcdas.cend() != i)
        {
            result.reset(i->second(parameters, options));
        }

        return result;
    }
}
