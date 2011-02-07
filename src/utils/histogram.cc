/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <src/utils/exception.hh>
#include <src/utils/histogram.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/stringify.hh>
#include <src/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <algorithm>
#include <list>

namespace eos
{
    template class WrappedForwardIterator<Histogram<1>::ConstIteratorTag, const Histogram<1>::Bin>;

    template <> struct Implementation<Histogram<1>>
    {
        std::list<Histogram<1>::Bin> bins;

        unsigned entries;

        unsigned underflow, overflow;

        Implementation() :
            entries(0),
            underflow(0),
            overflow(0)
        {
        }

        static bool falls_into(const double & value, const Histogram<1>::Bin & bin)
        {
            return (bin.lower <= value) && (value < bin.upper);
        }

        void insert(const double & value)
        {
            auto b = std::find_if(bins.begin(), bins.end(), std::bind(&Implementation<Histogram<1>>::falls_into, value, std::placeholders::_1));

            if (bins.end() == b)
            {
                if (value < bins.front().lower)
                    ++underflow;

                if (value >= bins.back().upper)
                    ++overflow;

                throw InternalError("Histogram<1>::insert(): No bin found to cover value '" + stringify(value) + "'");
            }

            b->value += 1.0;
            ++entries;
        }

        static bool is_right_of(const Histogram<1>::Bin & new_bin, const Histogram<1>::Bin & bin)
        {
            return ! (bin < new_bin);
        }

        void insert(const Histogram<1>::Bin & bin)
        {
            auto b = std::find_if(bins.begin(), bins.end(), std::bind(&Implementation<Histogram<1>>::is_right_of, bin, std::placeholders::_1));

            bins.insert(b, bin);
        }
    };

    Histogram<1>::Histogram() :
        PrivateImplementationPattern<Histogram<1>>(new Implementation<Histogram<1>>)
    {
    }

    Histogram<1>::~Histogram()
    {
    }

    Histogram<1>
    Histogram<1>::WithEqualBinning(const double & lower, const double & upper, const unsigned & count)
    {
        Histogram<1> result;

        double bin_width = std::abs(upper - lower) / count;
        for (unsigned i = 0 ; i < count ; ++i)
        {
            result.insert(Bin(lower + bin_width * i, lower + bin_width * (i + 1), 0.0));
        }

        return result;
    }

    void
    Histogram<1>::insert(const Bin & bin)
    {
        _imp->insert(bin);
    }

    void
    Histogram<1>::insert(const double & value)
    {
        _imp->insert(value);
    }

    unsigned
    Histogram<1>::entries() const
    {
        return _imp->entries;
    }

    Histogram<1>::ConstIterator
    Histogram<1>::begin() const
    {
        return ConstIterator(_imp->bins.cbegin());
    }

    Histogram<1>::ConstIterator
    Histogram<1>::end() const
    {
        return ConstIterator(_imp->bins.cend());
    }

    Histogram<1>
    estimate_cumulative_distribution(const Histogram<1> & distribution)
    {
        Histogram<1> result;
        unsigned entries = distribution.entries();
        double previous = 0.0;

        for (auto b = distribution.begin(), b_end = distribution.end() ; b != b_end ; ++b)
        {
            double current = b->value / entries;

            result.insert(Histogram<1>::Bin(b->lower, b->upper, previous + current));
            previous += current;
        }

        return result;
    }
}
