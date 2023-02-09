/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Danny van Dyk
 * Copyright (c) 2020 Nico Gubernari
 * Copyright (c) 2020 Christoph Bobeth
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

#include <eos/form-factors/parametric-bgl1997-impl.hh>

namespace eos
{
    const std::set<ReferenceName>
    BGL1997FormFactors<BToDstar>::references
    {
        "BGL:1997A"_rn
    };

    const std::vector<OptionSpecification>
    BGL1997FormFactors<BToDstar>::_options
    {
    };

    template class BGL1997FormFactors<BToDstar>;

    const std::set<ReferenceName>
    BGL1997FormFactors<BToD>::references
    {
        "BGL:1997A"_rn
    };

    const std::vector<OptionSpecification>
    BGL1997FormFactors<BToD>::_options
    {
    };

    template class BGL1997FormFactors<BToD>;
}
