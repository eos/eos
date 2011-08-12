/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
 * Copyright (c) 2011 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_RVALUE_HH
#define EOS_GUARD_SRC_UTILS_RVALUE_HH 1

#include <eos/utils/instantiation_policy.hh>

#include <vector>

namespace eos
{
    struct RValue :
        public InstantiationPolicy<RValue, NonInstantiable>
    {
        /*!
         * Calculate the R-value (actually the sqrt(R)) for a given quantity x (param, log(posterior), ...)
         * according to [GR1992], Eqs. (3),(4), p. 461. Using their notation.
         * Included DoF estimation for the t-distribution.
         *
         * @param chain_means the mean of x in each chain
         * @param chain_variances the sample variance  of x in each chain
         * @param chain_length The number of iterations used to calculate the
         * means and variances within each chain. Usually the length of the prerun until now.
         * Note that if only chain_length is increased and everything else is kept constant, the R-value should increase
         */
        static double gelman_rubin(const std::vector<double> & chain_means, const std::vector<double> & chain_variances,
                const unsigned & chain_length);

        /*!
         * Approximate the R-value (actually sqrt(R)) for a given quantitiy x.
         *
         * Here we use the approximation @$R \approx \sigma^2 / W@$, just as in
         * BAT v0.4.
         *
         * @param chain_means the mean of x in each chain
         * @param chain_variances the sample variance  of x in each chain
         * @param chain_length The number of iterations used to calculate the
         * means and variances within each chain. Usually the length of the prerun until now.
         * Note that if only chain_length is increased and everything else is kept constant, the R-value should increase
         */
        static double approximation(const std::vector<double> & chain_means, const std::vector<double> & chain_variances,
                const unsigned & chain_length);
    };
}

#endif
