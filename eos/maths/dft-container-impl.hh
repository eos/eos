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

#ifndef EOS_GUARD_EOS_MATHS_DFT_CONTAINER_IMPL_HH
#define EOS_GUARD_EOS_MATHS_DFT_CONTAINER_IMPL_HH 1

#include <eos/maths/dft-container.hh>

#include <eos/utils/exception.hh>

#include <algorithm>
#include <cstring>
#include <format>
#include <functional>
#include <new>
#include <stdexcept>
#include <type_traits>

#include <fftw3.h>

namespace eos
{
    namespace dft
    {
        namespace impl
        {
            template <std::size_t rank_>
            std::size_t size_from_dimensions(const std::array<std::size_t, rank_> & dimensions)
            {
                std::size_t result = 1;
                for (std::size_t d : dimensions)
                {
                    result *= d;
                }

                return result;
            }

            template <std::size_t rank_>
            std::size_t index_from_multi_index(const std::array<std::size_t, rank_> & multi_index, const std::array<std::size_t, rank_> & dimensions)
            {
                // computing the index in row-major order, i.e. the last index is contiguous in memory
                std::size_t result = 0, stride = 1;

                for (std::size_t i = rank_ - 1; i > 0 ; --i)
                {
                    result += multi_index[i] * stride;
                    stride *= dimensions[i];
                }

                result += multi_index[0] * stride;

                return result;
            }
        } // namespace eos::dft::impl

        // partial specialization for T_ == double

        template <std::size_t rank_>
        Container<double, rank_>::Container(const std::array<std::size_t, rank_> & dimensions) :
            _size(impl::size_from_dimensions(dimensions)),
            _data(fftw_alloc_real(this->_size)),
            _dimensions(dimensions)
        {
            if (nullptr == this->_data)
            {
                throw std::bad_alloc();
            }
        }

        template <std::size_t rank_>
        Container<double, rank_>::Container(const Container & other) :
            _size(other._size),
            _data(fftw_alloc_real(this->_size)),
            _dimensions(other._dimensions)
        {
            if (nullptr == this->_data)
            {
                throw std::bad_alloc();
            }

            memcpy(this->_data, other._data, this->_size * sizeof(double));
        }

        template <std::size_t rank_>
        Container<double, rank_>::Container(Container && other) :
            _size(other._size),
            _data(other._data),
            _dimensions(other._dimensions)
        {
            other._data = nullptr;
            other._size = 0;
        }

        template <std::size_t rank_>
        Container<double, rank_> & Container<double, rank_>::operator=(const Container & other)
        {
            if (this == &other)
            {
                return *this;
            }

            if (nullptr != this->_data)
            {
                fftw_free(this->_data);
                this->_data = nullptr;
                this->_size = 0;
            }

            this->_size = other._size;
            this->_data = fftw_alloc_real(this->_size);

            if (nullptr == this->_data)
            {
                throw std::bad_alloc();
            }

            this->_dimensions = other._dimensions;
            std::memcpy(this->_data, other._data, this->_size * sizeof(double));

            return *this;
        }

        template <std::size_t rank_>
        Container<double, rank_> & Container<double, rank_>::operator=(Container && other)
        {
            if (this == &other)
            {
                return *this;
            }

            if (nullptr != this->_data)
            {
                fftw_free(this->_data);
                this->_data = nullptr;
                this->_size = 0;
            }

            this->_size = other._size;
            this->_data = other._data;
            this->_dimensions = other._dimensions;

            other._data = nullptr;
            other._size = 0;

            return *this;
        }

        template <std::size_t rank_>
        Container<double, rank_>::~Container()
        {
            if (nullptr != this->_data)
            {
                fftw_free(this->_data);
                this->_data = nullptr;
                this->_size = 0;
            }
        }

        template <std::size_t rank_>
        const std::array<std::size_t, rank_> & Container<double, rank_>::dimensions() const
        {
            return this->_dimensions;
        }

        template <std::size_t rank_>
        double & Container<double, rank_>::operator()(const std::array<std::size_t, rank_> & index)
        {
            std::size_t _index = impl::index_from_multi_index(index, this->_dimensions);
            if (_index >= this->_size)
            {
                throw std::out_of_range(std::format("dft::Container<double, {}>: index out of bounds", rank_));
            }

            return this->_data[_index];
        }

        template <std::size_t rank_>
        const double & Container<double, rank_>::operator()(const std::array<std::size_t, rank_> & index) const
        {
            std::size_t _index = impl::index_from_multi_index(index, this->_dimensions);

            if (_index >= this->_size)
            {
                throw std::out_of_range(std::format("dft::Container<double, {}>: index out of bounds", rank_));
            }

            return this->_data[_index];
        }

        // partial specialization for T_ == std::complex<double>

        template <std::size_t rank_>
        Container<std::complex<double>, rank_>::Container(const std::array<std::size_t, rank_> & dimensions) :
            _size(impl::size_from_dimensions(dimensions)),
            _data(reinterpret_cast<std::complex<double> *>(fftw_alloc_complex(this->_size))),
            _dimensions(dimensions)
        {
            static_assert(sizeof(std::complex<double>) == 2 * sizeof(double), "std::complex<double> is not compatible with fftw_complex (size mismatch)");
            static_assert(alignof(std::complex<double>) == alignof(fftw_complex), "std::complex<double> is not compatible with fftw_complex (alignment mismatch)");
            static_assert(std::is_trivially_copyable<std::complex<double>>::value, "std::complex<double> is not compatible with fftw_complex (not trivially copyable)");

            if (nullptr == this->_data)
            {
                throw std::bad_alloc();
            }
        }

        template <std::size_t rank_>
        Container<std::complex<double>, rank_>::Container(const Container & other) :
            _size(other._size),
            _data(reinterpret_cast<std::complex<double> *>(fftw_alloc_complex(this->_size))),
            _dimensions(other._dimensions)
        {
            if (nullptr == this->_data)
            {
                throw std::bad_alloc();
            }

            std::memcpy(this->_data, other._data, this->_size * 2 * sizeof(double));
        }

        template <std::size_t rank_>
        Container<std::complex<double>, rank_>::Container(Container && other) :
            _size(other._size),
            _data(other._data),
            _dimensions(other._dimensions)
        {
            other._data = nullptr;
            other._size = 0;
        }

        template <std::size_t rank_>
        Container<std::complex<double>, rank_> & Container<std::complex<double>, rank_>::operator=(const Container & other)
        {
            if (this == &other)
            {
                return *this;
            }

            if (nullptr != this->_data)
            {
                fftw_free(this->_data);
            }

            this->_size = other._size;
            this->_data = reinterpret_cast<std::complex<double> *>(fftw_alloc_complex(this->_size));

            if (nullptr == this->_data)
            {
                throw std::bad_alloc();
            }

            this->_dimensions = other._dimensions;
            std::memcpy(this->_data, other._data, this->_size * 2 * sizeof(double));

            return *this;
        }

        template <std::size_t rank_>
        Container<std::complex<double>, rank_> & Container<std::complex<double>, rank_>::operator=(Container && other)
        {
            if (this == &other)
            {
                return *this;
            }

            if (nullptr != this->_data)
            {
                fftw_free(this->_data);
                this->_data = nullptr;
            }

            this->_size = other._size;
            this->_data = other._data;
            this->_dimensions = other._dimensions;

            other._data = nullptr;
            other._size = 0;

            return *this;
        }

        template <std::size_t rank_>
        Container<std::complex<double>, rank_>::~Container()
        {
            if (nullptr != this->_data)
            {
                fftw_free(this->_data);
                this->_data = nullptr;
                this->_size = 0;
            }
        }

        template <std::size_t rank_>
        const std::array<std::size_t, rank_> & Container<std::complex<double>, rank_>::dimensions() const
        {
            return this->_dimensions;
        }

        template <std::size_t rank_>
        std::complex<double> & Container<std::complex<double>, rank_>::operator()(const std::array<std::size_t, rank_> & index)
        {
            std::size_t _index = impl::index_from_multi_index(index, this->_dimensions);
            if (_index >= this->_size)
            {
                throw std::out_of_range(std::format("dft::Container<std::complex<double>, {}>: index out of bounds", rank_));
            }

            return this->_data[_index];
        }

        template <std::size_t rank_>
        const std::complex<double> & Container<std::complex<double>, rank_>::operator()(const std::array<std::size_t, rank_> & index) const
        {
            std::size_t _index = impl::index_from_multi_index(index, this->_dimensions);

            if (_index >= this->_size)
            {
                throw std::out_of_range(std::format("dft::Container<std::complex<double>, {}>: index out of bounds", rank_));
            }

            return this->_data[_index];
        }

        template <std::size_t rank_>
        Container<std::complex<double>, rank_> & Container<std::complex<double>, rank_>::operator*=(const Container & rhs)
        {
            if (this->_dimensions != rhs._dimensions)
            {
                throw InternalError(std::format("dft::Container<std::complex<double>, {}>: dimension mismatch in operator*=", rank_));
            }

            std::transform(this->_data, this->_data + this->_size, rhs._data, this->_data,
                std::multiplies<>{});

            return *this;
        }
    }
}

#endif
