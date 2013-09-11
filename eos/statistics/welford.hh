/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
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

#ifndef EOS_GUARD_EOS_UTILS_WELFORD_HH
#define EOS_GUARD_EOS_UTILS_WELFORD_HH 1

namespace eos
{
    /*!
     * Calculate running mean and variance according to Welford's method,
     * cf. Knuth, D. TAOCP, vol. 2, 3rd edition, p. 232,
     * or http://www.johndcook.com/standard%5Fdeviation.html
     */
    struct Welford
    {
        private:
            double new_mean, old_mean;

            double new_sum, old_sum;

            unsigned size;

        public:
            Welford();

            void add(const double & value);

            double mean() const;

            unsigned number_of_elements() const;

            double std_deviation() const;

            double variance() const;
    };
}

#endif
