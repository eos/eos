/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021-2025 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_UNITS_HH
#define EOS_GUARD_EOS_UTILS_UNITS_HH 1

#include <eos/utils/instantiation_policy.hh>

#include <string>

namespace eos
{
    class Unit
    {
        public:
            enum class Id;

        private:
            Id _id;

            static const std::vector<std::string> _latex_representations;

            static const std::vector<std::string> _internal_representations;

            Unit(const Id & id) :
                _id(id)
            {
            }

        public:
            enum class Id : int
            {
                undefined = 0,
                none,
                gev,
                gev2,
                gev3,
                inverse_gev,
                inverse_gev2,
                inverse_gev4,
                s,
                inverse_s,
                inverse_ps,
                gev_s,
                fm2
            };

            Unit(const std::string &);
            Unit(const Unit &) = default;
            Unit(Unit &&)      = default;
            ~Unit()            = default;

            const std::string & latex() const;
            const std::string & string() const;

            static Unit Undefined();
            static Unit None();
            static Unit GeV();
            static Unit GeV2();
            static Unit GeV3();
            static Unit InverseGeV();
            static Unit InverseGeV2();
            static Unit InverseGeV4();
            static Unit Second();
            static Unit InverseSecond();
            static Unit InversePicoSecond();
            static Unit GeVSecond();
            static Unit Femtometer2();

            bool   operator== (const Unit &) const;
            Unit & operator= (const Unit &) = default;
    };
} // namespace eos

#endif
