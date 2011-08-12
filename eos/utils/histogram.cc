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

#include <eos/utils/exception.hh>
#include <eos/utils/histogram.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <algorithm>
#include <list>
#include <vector>

namespace eos
{
    /* Histogram<1> */

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

    /* Histogram<2> */

    template class WrappedForwardIterator<Histogram<2>::ConstIteratorTag, const Histogram<2>::Bin>;

    template class WrappedForwardIterator<Histogram<2>::IteratorTag, Histogram<2>::Bin>;

    template <> struct Implementation<Histogram<2>>
    {
        std::vector<Histogram<2>::Bin> bins;

        Implementation()
        {
        }

        static bool falls_into(const std::array<double, 2> & coordinates, const Histogram<2>::Bin & bin)
        {
            bool result = (((bin.lower[0] <= coordinates[0]) && (coordinates[0] < bin.upper[0]))
                && (bin.lower[1] <= coordinates[1]) && (coordinates[1] < bin.upper[1]));

            return result;
        }

        std::vector<Histogram<2>::Bin>::iterator find(const std::array<double, 2> & coordinates)
        {
            // perform a binary search for a matching bin
            std::vector<Histogram<2>::Bin>::iterator l = bins.begin(), r = bins.end(), x;

            while (r - l > 2)
            {
                x = l + (r - l) / 2;

                if (coordinates[0] < x->lower[0])
                {
                    r = x;
                    continue;
                }

                if (x->upper[0] < coordinates[0])
                {
                    l = x;
                    continue;
                }

                break;
            };

            auto result = std::find_if(l, r, std::bind(&Implementation<Histogram<2>>::falls_into, coordinates, std::placeholders::_1));
            if (result == r)
                result = bins.end();

            return result;
        }

        static bool is_next_of(const Histogram<2>::Bin & new_bin, const Histogram<2>::Bin & bin)
        {
            return ! (bin < new_bin);
        }

        void insert(const Histogram<2>::Bin & bin)
        {
            auto b = std::find_if(bins.begin(), bins.end(), std::bind(&Implementation<Histogram<2>>::is_next_of, bin, std::placeholders::_1));

            bins.insert(b, bin);
        }
    };

    Histogram<2>::Histogram() :
        PrivateImplementationPattern<Histogram<2>>(new Implementation<Histogram<2>>)
    {
    }

    Histogram<2>::~Histogram()
    {
    }

    Histogram<2>
    Histogram<2>::WithEqualBinning(const std::array<double, 2> & lower, const std::array<double, 2> & upper, const std::array<unsigned, 2> & count)
    {
        Histogram<2> result;

        double x_bin_width = std::abs(upper[0] - lower[0]) / count[0];
        double y_bin_width = std::abs(upper[1] - lower[1]) / count[1];
        for (unsigned i = 0 ; i < count[0] ; ++i)
        {
            for (unsigned j = 0 ; j < count[1] ; ++j)
            {
                result.insert(Bin(std::array<double, 2>{{lower[0] + x_bin_width * i, lower[1] + y_bin_width * j}},
                            std::array<double, 2>{{lower[0] + x_bin_width * (i + 1), lower[1] + y_bin_width * (j + 1)}},
                            0.0));
            }
        }

        return result;
    }

    void
    Histogram<2>::insert(const Bin & bin)
    {
        _imp->insert(bin);
    }

    Histogram<2>::Iterator
    Histogram<2>::begin()
    {
        return Iterator(_imp->bins.begin());
    }

    Histogram<2>::Iterator
    Histogram<2>::end()
    {
        return Iterator(_imp->bins.end());
    }

    Histogram<2>::Iterator
    Histogram<2>::find(const std::array<double, 2> & coordinates)
    {
        return Iterator(_imp->find(coordinates));
    }

    Histogram<2>::ConstIterator
    Histogram<2>::cbegin() const
    {
        return ConstIterator(_imp->bins.cbegin());
    }

    Histogram<2>::ConstIterator
    Histogram<2>::cend() const
    {
        return ConstIterator(_imp->bins.cend());
    }
}
