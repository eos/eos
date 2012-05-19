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

#include <eos/utils/welford.hh>

#include <cmath>

namespace eos
{
    Welford::Welford() :
        new_mean(0),
        old_mean(0),
        new_sum(0),
        old_sum(0),
        size(0)
    {
    }

    void
    Welford::add(const double & value)
    {
        ++size;

        if (size == 1)
        {
            old_mean = new_mean = value;
            old_sum = 0;
        }
        else
        {
            new_mean = old_mean + (value - old_mean) / size;
            new_sum = old_sum + (value - old_mean) * (value - new_mean);

            // setup for next iteration
            old_mean = new_mean;
            old_sum = new_sum;
        }
    }

    double
    Welford::mean() const
    {
        return (size > 0 ) ? new_mean: 0;
    }

    unsigned
    Welford::number_of_elements() const
    {
        return size;
    }

    double
    Welford::std_deviation() const
    {
        return (size > 1) ? std::sqrt(new_sum / (size - 1)) : 0;
    }

    double
    Welford::variance() const
    {
        return (size > 1) ? new_sum / (size - 1) : 0;
    }
}
