/* vim: set sw=4 sts=4 tw=120 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Danny van Dyk
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

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /* Unitarity bound implemented as discussed in [BJvD2019] */
    class HQETUnitarityBounds :
        public virtual ParameterUser,
        public PrivateImplementationPattern<HQETUnitarityBounds>
    {
        public:
            HQETUnitarityBounds(const Parameters &, const Options &);
            ~HQETUnitarityBounds();

            double bound_0p() const;

            double bound_0m() const;

            double bound_1p() const;

            double bound_1m() const;
    };
}
