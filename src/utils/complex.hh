/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WILSON_FITTER_GUARD_SRC_UTILS_COMPLEX_HH
#define WILSON_FITTER_GUARD_SRC_UTILS_COMPLEX_HH 1

namespace wf
{
    template <typename T_>
    class Complex;

    template <typename T_>
    Complex<T_> operator* (const Complex<T_> & lhs, const Complex<T_> & rhs);

    template <typename T_>
    Complex<T_> operator+ (const Complex<T_> & lhs, const Complex<T_> & rhs);

    template <typename T_>
    Complex<T_> operator- (const Complex<T_> & lhs, const Complex<T_> & rhs);

    template <typename T_>
    class Complex
    {
        private:
            T_ _real_part;

            T_ _imaginary_part;

            Complex(const T_ &, const T_ &);

        public:
            friend Complex<T_> wf::operator*<>(const Complex<T_> & lhs, const Complex<T_> & rhs);
            friend Complex<T_> wf::operator+<>(const Complex<T_> & lhs, const Complex<T_> & rhs);
            friend Complex<T_> wf::operator-<>(const Complex<T_> & lhs, const Complex<T_> & rhs);

            /* Constructors */

            Complex();

            Complex(const Complex<T_> &);

            static Complex<T_> Cartesian(const T_ & real, const T_ & imaginary);

            static Complex<T_> Polar(const T_ & modulus, const T_ & argument);

            /* Operators */
            Complex<T_> & operator= (const Complex<T_> & rhs);

            T_ absolute() const;

            T_ absolute_squared() const;

            Complex<T_> conjugate() const;

            T_ & real();

            T_ & imaginary();
    };

    template <typename T_>
    Complex<T_> operator* (const Complex<T_> & lhs, const T_ & rhs)
    {
        return lhs * Complex<T_>::Cartesian(rhs, 0);
    }

    template <typename T_>
    Complex<T_> operator* (const T_ & lhs, const Complex<T_> & rhs)
    {
        return rhs * lhs;
    }

    template <typename T_>
    Complex<T_> operator+ (const Complex<T_> & lhs, const T_ & rhs)
    {
        return lhs + Complex<T_>::Cartesian(rhs, 0);
    }

    template <typename T_>
    Complex<T_> operator+ (const T_ & lhs, const Complex<T_> & rhs)
    {
        return rhs + lhs;
    }

    template <typename T_>
    Complex<T_> operator- (const Complex<T_> & lhs, const T_ & rhs)
    {
        return lhs - Complex<T_>::Cartesian(rhs, 0);
    }

    template <typename T_>
    Complex<T_> operator- (const T_ & lhs, const Complex<T_> & rhs)
    {
        return Complex<T_>::Cartesian(lhs, 0) - rhs;
    }

    extern template class Complex<double>;
}

#endif
