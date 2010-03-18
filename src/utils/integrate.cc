/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/complex-impl.hh>
#include <src/utils/integrate.hh>

namespace wf
{
    double integrate(const std::tr1::function<double (const double &)> & f, unsigned n, const double & a, const double & b)
    {
        if (n & 0x1)
            n += 1;

        double h = (b - a) / n;

        double result = f(a) + f(b);
        for (unsigned i = 1 ; i <= n / 2 - 1 ; ++i)
        {
            result += 2.0 * f(a + h * (2.0 * i));
        }
        for (unsigned i = 1 ; i <= n / 2 ; ++i)
        {
            result += 4.0 * f(a + h * (2.0 * i - 1.0));
        }

        result *= h / 3;

        return result;
    }

    Complex<double> integrate(const std::tr1::function<Complex<double> (const double &)> & f, unsigned n, const double & a, const double & b)
    {
        if (n & 0x1)
            n += 1;

        double h = (b - a) / n;

        Complex<double> result = f(a) + f(b);
        for (unsigned i = 1 ; i <= n / 2 - 1 ; ++i)
        {
            result = result + 2.0 * f(a + h * (2.0 * i));
        }
        for (unsigned i = 1 ; i <= n / 2 ; ++i)
        {
            result = result + 4.0 * f(a + h * (2.0 * i - 1.0));
        }

        result = result * (h / 3.0);

        return result;
    }
}
