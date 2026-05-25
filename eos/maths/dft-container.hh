/*
 * Copyright (c) 2026 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_MATHS_DFT_CONTAINER_HH
#define EOS_GUARD_EOS_MATHS_DFT_CONTAINER_HH 1

#include <array>
#include <complex>

namespace eos
{
    namespace dft
    {
        template <typename T_, std::size_t rank_>
        class Container;

        template <std::size_t rank_>
        class Container<double, rank_>
        {
            private:
                std::size_t _size;
                double * _data;
                std::array<std::size_t, rank_> _dimensions;

            public:
                Container(const std::array<std::size_t, rank_> & dimensions);
                Container(const Container & other);
                Container(Container && other);
                Container & operator=(const Container & other);
                Container & operator=(Container && other);
                ~Container();

                const std::array<std::size_t, rank_> & dimensions() const;

                double & operator()(const std::array<std::size_t, rank_> & index);
                const double & operator()(const std::array<std::size_t, rank_> & index) const;

                inline double * data() { return this->_data; };
                inline const double * data() const { return this->_data; };
        };

        extern template class Container<double, 1>;
        extern template class Container<double, 2>;
        extern template class Container<double, 3>;
        extern template class Container<double, 4>;

        template <std::size_t rank_>
        class Container<std::complex<double>, rank_>
        {
            private:
                std::size_t _size;
                std::complex<double> * _data;
                std::array<std::size_t, rank_> _dimensions;

            public:
                Container(const std::array<std::size_t, rank_> & dimensions);
                Container(const Container & other);
                Container(Container && other);
                Container & operator=(const Container & other);
                Container & operator=(Container && other);
                ~Container();

                const std::array<std::size_t, rank_> & dimensions() const;

                std::complex<double> & operator()(const std::array<std::size_t, rank_> & index);
                const std::complex<double> & operator()(const std::array<std::size_t, rank_> & index) const;

                inline std::complex<double> * data() { return this->_data; };
                inline const std::complex<double> * data() const { return this->_data; };

                Container & operator*=(const Container & rhs);
        };

        template <std::size_t rank_>
        Container<std::complex<double>, rank_> operator*(Container<std::complex<double>, rank_> lhs, const Container<std::complex<double>, rank_> & rhs)
        {
            lhs *= rhs;
            return lhs;
        }

        extern template class Container<std::complex<double>, 1>;
        extern template class Container<std::complex<double>, 2>;
        extern template class Container<std::complex<double>, 3>;
        extern template class Container<std::complex<double>, 4>;
    } // namespace eos::dft
} // namespace eos

#endif
