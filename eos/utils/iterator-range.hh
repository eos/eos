/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_ITERATOR_RANGE_HH
#define EOS_GUARD_EOS_UTILS_ITERATOR_RANGE_HH 1

namespace eos
{
    /**
     * Convenience wrapper around a pair of iterators, such
     * that C++11 range-based for loops can easily be used.
     */
    template <typename Iterator_> class IteratorRange
    {
        private:
            Iterator_ _begin;
            Iterator_ _end;

        public:
            IteratorRange(const Iterator_ & begin, const Iterator_ & end) :
                _begin(begin),
                _end(end)
            {
            }

            ~IteratorRange() = default;

            Iterator_
            begin()
            {
                return _begin;
            }

            Iterator_
            end()
            {
                return _end;
            }
    };
} // namespace eos

#endif
