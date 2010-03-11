/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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
            result += 2.0 * f(a + h * (2 * i));
        }
        for (unsigned i = 1 ; i <= n / 2 ; ++i)
        {
            result += 4.0 * f(a + h * (2 * i - 1));
        }

        result *= h / 3;

        return result;
    }
}
