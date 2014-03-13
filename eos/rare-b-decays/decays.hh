/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2014 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_DECAYS_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_DECAYS_HH 1

namespace eos
{
    /*
     * Phase space tags
     */
    struct LargeRecoil { };
    struct LowRecoil { };

    template <typename T_> class BToKstarDilepton;

    template <typename T_> class BToKDilepton;

    template <typename T_> class BToXsDilepton;

    template <typename T_> class BToXsGamma;

    template <typename T_> class LambdaBToLambdaDilepton;

    enum Helicity
    {
        left_handed = -1,
        right_handed = +1
    };
}

#endif
