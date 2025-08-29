/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_DILEPTON_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_DILEPTON_HH 1

#include <eos/observable.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /*!
     * Calculation according to [BEKU:2002] and [F:2012].
     */
    class BToDilepton :
        public ParameterUser,
        public PrivateImplementationPattern<BToDilepton>
    {
        public:
            ///@name Basic Functions
            ///@{
            /// Constructor.
            BToDilepton(const Parameters & parameters, const Options & options);

            /// Destructor.
            ~BToDilepton();
            ///@}

            ///@name Observables
            ///@{
            /// Branching ratio at time t = 0, considering B-Bbar mixing effects.
            double branching_ratio_time_zero() const;

            /// Time-integrated untagged branching ratio, considering B-Bbar mixing effects.
            double branching_ratio_untagged_integrated() const;

            /// CP asymmetry due to finite decay-width difference
            double cp_asymmetry_del_gamma() const;

            /// CP asymmetry due to mixing
            double cp_asymmetry_mixing_S() const;

            /// Effective lifetime
            double effective_lifetime() const;
            ///@}

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
