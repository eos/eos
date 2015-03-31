/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013, 2015 Danny van Dyk
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

#ifndef BSTOKSTARLNU_GUARD_EOS_B_DECAYS_B_TO_L_NU_HH
#define BSTOKSTARLNU_GUARD_EOS_B_DECAYS_B_TO_L_NU_HH 1

#include <eos/decays.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    class BToLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToLeptonNeutrino>
    {
        public:
            BToLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToLeptonNeutrino();

            // Observables
            double branching_ratio() const;
            double decay_width() const;
    };
}

#endif

