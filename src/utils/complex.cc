/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/complex.hh>
#include <src/utils/complex-impl.hh>

namespace wf
{
    template class Complex<double>;

    template Complex<double> operator*<double> (const Complex<double> &, const Complex<double> &);

    template Complex<double> operator+<double> (const Complex<double> &, const Complex<double> &);

    template Complex<double> operator-<double> (const Complex<double> &, const Complex<double> &);
}

