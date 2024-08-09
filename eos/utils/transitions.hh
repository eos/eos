/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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

#ifndef EOS_GUARD_EOS_UTILS_TRANSITIONS_HH
#define EOS_GUARD_EOS_UTILS_TRANSITIONS_HH 1

namespace eos
{
    /* Mesonic Tags */

    struct PToV { };
    struct PToGamma { };
    struct PToGammaOffShell { };
    struct PToP { };
    struct PToPP { };
    struct VToP { };
    struct VToV { };

    /* Baryonic Tags */

    // J=1/2^+ -> J=1/2^+ transitions
    struct OneHalfPlusToOneHalfPlus { };
    // J=1/2^+ -> J=1/2^- transitions
    struct OneHalfPlusToOneHalfMinus { };
    // J=1/2^+ -> J=3/2^- transitions
    struct OneHalfPlusToThreeHalfMinus { };

    /* Scattering Tags */
    struct PPToPP {};
}
#endif
