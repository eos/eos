/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#ifndef EOS_GUARD_EOS_STATISTICS_GOODNESS_OF_FIT_HH
#define EOS_GUARD_EOS_STATISTICS_GOODNESS_OF_FIT_HH 1

#include <eos/statistics/log-posterior.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/wrapped_forward_iterator-fwd.hh>

namespace eos
{
    class GoodnessOfFit :
        public PrivateImplementationPattern<GoodnessOfFit>
    {
        public:
            GoodnessOfFit(const LogPosterior &);
            ~GoodnessOfFit();

            ///@name ChiSquare test statistics
            ///@{
            double total_chi_square() const;
            int total_degrees_of_freedom() const;

            struct ChiSquareIteratorTag;
            typedef WrappedForwardIterator<ChiSquareIteratorTag, const std::pair<const QualifiedName, test_statistics::ChiSquare>> ChiSquareIterator;

            ChiSquareIterator begin_chi_square() const;
            ChiSquareIterator end_chi_square() const;
            ///q}
    };
}

#endif
