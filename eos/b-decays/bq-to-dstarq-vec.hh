/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Stefan Meiser
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

#ifndef EOS_GUARD_EOS_B_DECAYS_BQ_TO_DSTARQ_VEC_HH
#define EOS_GUARD_EOS_B_DECAYS_BQ_TO_DSTARQ_VEC_HH 1

#include <eos/maths/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    class BqToDstarqVector :
        public ParameterUser,
        public PrivateImplementationPattern<BqToDstarqVector>
    {
        public:
            BqToDstarqVector(const Parameters & parameters, const Options & options);
            ~BqToDstarqVector();

            // Observables
            double branching_ratio() const;
            double decay_width() const;

            // Reduced matrix elements
            double re_a_1() const;
            double im_a_1() const;

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
