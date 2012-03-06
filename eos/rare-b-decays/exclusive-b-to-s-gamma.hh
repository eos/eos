/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_GAMMA_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_GAMMA_HH 1

#include <eos/rare-b-decays/decays.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*!
     * Calculates observables in B->K^*gamma decays
     */
    class BToKstarGamma :
        public ParameterUser,
        public PrivateImplementationPattern<BToKstarGamma>
    {
        public:
            BToKstarGamma(const Parameters & parameters, const Options & options);
            ~BToKstarGamma();

            /// Branching Ratio
            double branching_ratio() const;

            /// Branching Ratio (CP averaged)
            double branching_ratio_cp_averaged() const;

            /// Direct CP asymmetry A_CP
            double cp_asymmetry() const;

            /// Time dependent CP asymmetry S_K^*gamma
            double s_kstar_gamma() const;

            /// Time dependent CP asymmetry C_K^*gamma
            double c_kstar_gamma() const;

            /// Isospin asymmetry
            double isospin_asymmetry() const;
    };
}

#endif
