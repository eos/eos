/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Nico Gubernari
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

#include <eos/form-factors/parametric-g2026-impl.hh>

namespace eos
{
    // b -> u,d
    template class G2026FormFactors<BToPi, PToP>;
    template class G2026FormFactors<BsToK, PToP>;
    template class G2026FormFactors<BsToKstar, PToV>;

    // b -> s
    template class G2026FormFactors<BToK, PToP>;
    template class G2026FormFactors<BToKstar, PToV>;
    template class G2026FormFactors<BsToPhi, PToV>;

    // b -> c
    template class G2026FormFactors<BToD, PToP>;
    template class G2026FormFactors<BsToDs, PToP>;
    template class G2026FormFactors<BToDstar, PToV>;
    template class G2026FormFactors<BsToDsstar, PToV>;
}
