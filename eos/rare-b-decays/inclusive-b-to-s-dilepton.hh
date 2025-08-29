/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_INCLUSIVE_B_TO_S_DILEPTON_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_INCLUSIVE_B_TO_S_DILEPTON_HH 1

#include <eos/maths/complex.hh>
#include <eos/rare-b-decays/decays.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /*
     * Decay: B -> X_s l lbar
     */


    // As given in [HLMW:2005]
    struct HLMW2005
    {
    };

    template <>
    class BToXsDilepton<HLMW2005> :
        public ParameterUser,
        public PrivateImplementationPattern<BToXsDilepton<HLMW2005>>
    {
        public:
            BToXsDilepton(const Parameters & parameters, const Options & options);
            ~BToXsDilepton();

            // Differential Observables
            double differential_branching_ratio(const double & s) const;

            // Integrated Observables
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;

            // Diagnostic Values
            Diagnostics diagnostics() const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };
}

#endif
