/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WILSON_FITTER_GUARD_SRC_UTILS_COMPLEX_IMPL_HH
#define WILSON_FITTER_GUARD_SRC_UTILS_COMPLEX_IMPL_HH 1

#include <src/utils/complex.hh>
#include <src/utils/exception.hh>

#include <cmath>

namespace wf
{
    template <typename T_>
    Complex<T_>::Complex(const T_ & x, const T_ & y) :
        _real_part(x),
        _imaginary_part(y)
    {
    }

    template <typename T_>
    Complex<T_>::Complex() :
        _real_part(0.0),
        _imaginary_part(0.0)
    {
    }

    template <typename T_>
    Complex<T_>::Complex(const Complex<T_> & other) :
        _real_part(other._real_part),
        _imaginary_part(other._imaginary_part)
    {
    }

    template <typename T_>
    Complex<T_>
    Complex<T_>::Cartesian(const T_ & real, const T_ & imaginary)
    {
        return Complex(real, imaginary);
    }

    template <typename T_>
    Complex<T_>
    Complex<T_>::Polar(const T_ & modulus, const T_ & argument)
    {
        return modulus * Complex(std::cos(argument), std::sin(argument));
    }

    template <typename T_>
    Complex<T_> &
    Complex<T_>::operator= (const Complex<T_> & rhs)
    {
        this->_real_part = rhs._real_part;
        this->_imaginary_part = rhs._imaginary_part;

        return *this;
    }

    template <typename T_>
    T_
    Complex<T_>::absolute() const
    {
        return std::sqrt(this->_real_part * this->_real_part + this->_imaginary_part * this->_imaginary_part);
    }

    template <typename T_>
    T_
    Complex<T_>::absolute_squared() const
    {
        return this->_real_part * this->_real_part + this->_imaginary_part * this->_imaginary_part;
    }

    template <typename T_>
    Complex<T_>
    Complex<T_>::conjugate() const
    {
        Complex<T_> result(*this);

        result._imaginary_part *= -1.0;

        return result;
    }

    template <typename T_>
    T_ &
    Complex<T_>::real()
    {
        return this->_real_part;
    }

    template <typename T_>
    T_ &
    Complex<T_>::imaginary()
    {
        return this->_imaginary_part;
    }

    template <typename T_>
    T_
    Complex<T_>::phase()
    {
        T_ r(this->_real_part), i(this->_imaginary_part);

        if ((0 == r) && (0 == i))
            throw InternalError("Calculating phase of 0 + 0i");

        return std::atan2(r, i);
    }

    template <typename T_>
    Complex<T_>
    operator* (const Complex<T_> & lhs, const Complex<T_> & rhs)
    {
        Complex<T_> result;

        result._real_part = lhs._real_part * rhs._real_part - lhs._imaginary_part * rhs._imaginary_part;
        result._imaginary_part = lhs._real_part * rhs._imaginary_part + lhs._imaginary_part * rhs._real_part;

        return result;
    }

    template <typename T_>
    Complex<T_>
    operator+ (const Complex<T_> & lhs, const Complex<T_> & rhs)
    {
        Complex<T_> result;

        result._real_part = lhs._real_part + rhs._real_part;
        result._imaginary_part = lhs._imaginary_part + rhs._imaginary_part;

        return result;
    }

    template <typename T_>
    Complex<T_>
    operator- (const Complex<T_> & lhs, const Complex<T_> & rhs)
    {
        Complex<T_> result;

        result._real_part = lhs._real_part - rhs._real_part;
        result._imaginary_part = lhs._imaginary_part - rhs._imaginary_part;

        return result;
    }
}

#endif
