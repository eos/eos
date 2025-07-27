/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
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

#ifndef EOS_GUARD_EOS_UTILS_CARTESIAN_PRODUCT_HH
#define EOS_GUARD_EOS_UTILS_CARTESIAN_PRODUCT_HH 1

#include <cstddef>
#include <vector>

namespace eos
{
    /*!
     * @brief CartesianProduct is a template class to represent a Cartesian Product of an
     * arbitrary number of T_. T_ should be a container providing a random access iterator.
     */
    template <typename T_> class CartesianProduct
    {
        private:
            using IteratorState = std::vector<typename T_::const_iterator>;

            // All stored containers.
            std::vector<T_> _data;

            // Iterators for the begin and end of each stored container.
            IteratorState _begin, _end;

            // The overall number of elements stored in CartesianProduct.
            size_t _size;

            using Sizes = std::vector<size_t>;
            Sizes _sizes;

            template <typename U_> class _Iterator
            {
                private:
                    using IteratorState = std::vector<typename U_::const_iterator>;

                    IteratorState     _state;
                    Sizes             _sizes;
                    bool              _at_end;
                    std::vector<long> _values;

                    _Iterator(const IteratorState & state, const Sizes & sizes, bool atEnd) :
                        _state(state),
                        _sizes(sizes),
                        _at_end(atEnd),
                        _values(state.size(), 0)
                    {
                    }

                public:
                    friend class CartesianProduct<U_>;

                    bool
                    operator== (const _Iterator & other) const
                    {
                        if (_at_end && other._at_end)
                        {
                            return true;
                        }

                        // compare each element
                        auto i = _state.cbegin();
                        auto j = other._state.cbegin();
                        for (; i != _state.end(); ++i, ++j)
                        {
                            if (*i != *j)
                            {
                                return false;
                            }
                        }

                        return true;
                    }

                    bool
                    operator!= (const _Iterator & other) const
                    {
                        return ! (*this == other);
                    }

                    _Iterator &
                    operator+= (long increment)
                    {
                        // Already at the last element?
                        if (_at_end)
                        {
                            return *this;
                        }

                        auto size_it  = _sizes.rbegin();
                        auto state_it = _state.rbegin();
                        auto value_it = _values.rbegin();
                        for (; size_it != _sizes.rend(); ++size_it, ++state_it, ++value_it)
                        {
                            long new_value = (*value_it + increment) % *size_it;
                            ;
                            long difference  = new_value - *value_it;
                            *value_it        = new_value;
                            increment       -= difference;
                            increment       /= *size_it;
                            *state_it       += difference;
                        }

                        if (increment != 0)
                        {
                            _at_end = true;
                        }

                        return *this;
                    }

                    _Iterator &
                    operator++ ()
                    {
                        return (*this += 1);
                    }

                    std::vector<typename U_::value_type>
                    operator* () const
                    {
                        std::vector<typename U_::value_type> result;
                        for (auto i = _state.begin(); _state.end() != i; ++i)
                        {
                            result.push_back(*(*i));
                        }

                        return result;
                    }
            };

        public:
            /*!
             * @brief This is a random access iterator for CartesianPoduct.
             */
            using Iterator = _Iterator<T_>;

            /*!
             * Initializes an empty CartesianProduct.
             */
            CartesianProduct() :
                _size(0)
            {
            }

            /*!
             * Adds another container of type T_ to the CartesianProduct. So in some sense a multiplication
             * with new_container.
             * @param new_container The new container to be added.
             */
            void
            over(const T_ & new_container)
            {
                // Store the new container and its iterators
                _data.push_back(new_container);
                _begin.push_back(_data.back().cbegin());
                _end.push_back(_data.back().cend());

                // Update the size of CartesianProduct
                if (0 == _size)
                {
                    _size = new_container.size();
                }
                else
                {
                    _size *= new_container.size();
                }

                _sizes.push_back(new_container.size());
            }

            /*!
             * Returns an iterator for the first element.
             */
            Iterator
            begin() const
            {
                return Iterator(_begin, _sizes, false);
            }

            /*!
             * Returns an iterator for the position after the last element.
             */
            Iterator
            end() const
            {
                return Iterator(_end, _sizes, true);
            }

            /*!
             * Returns the number of elements stored in CartesianProduct. This is the product of the
             * sizes of each stored container.
             */
            size_t
            size() const
            {
                return _size;
            }
    };
} // namespace eos

#endif
