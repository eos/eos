/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/integrate.hh>

#include <limits>
#include <vector>

namespace wf
{
    using std::abs;
    using std::real;
    using std::imag;

    double integrate(const std::function<double (const double &)> & f, unsigned n, const double & a, const double & b)
    {
        if (n & 0x1)
            n += 1;

        if (n < 16)
            n = 16;

        double h = (b - a) / n;
        std::vector<double> y;

        for (unsigned k(0) ; k < n + 1 ; ++k)
        {
            y.push_back(f(a + k * h));
        }

        double Q0 = 0.0, Q1 = 0.0, Q2 = 0.0;
        for (unsigned k(0) ; k < n / 8 ; ++k)
        {
            Q0 += y[8 * k] + 4.0 * y[8 * k + 4] + y[8 * k + 4];
        }
        for (unsigned k(0) ; k < n / 4 ; ++k)
        {
            Q1 += y[4 * k] + 4.0 * y[4 * k + 2] + y[4 * k + 4];
        }
        for (unsigned k(0) ; k < n / 2 ; ++k)
        {
            Q2 += y[2 * k] + 4.0 * y[2 * k + 1] + y[2 * k + 2];
        }

        Q0 = Q0 * h / 3.0 * 4.0;
        Q1 = Q1 * h / 3.0 * 2.0;
        Q2 = Q2 * h / 3.0;

        double denom = (Q0 + Q2 - 2.0 * Q1);
        double num = Q2 - Q1;
        double correction = num * num / denom;
        double result;

        if (isnan(correction))
        {
            result = Q2;
        }
        else if (abs(correction / Q2) < 1.0)
        {
            result = Q2 - correction;
        }
        else
        {
#if 0
            std::cerr << "Q0 = " << Q0 << std::endl;
            std::cerr << "Q1 = " << Q1 << std::endl;
            std::cerr << "Q2 = " << Q2 << std::endl;
            std::cerr << "Reintegrating with twice the number of data points" << std::endl;
#endif
            result = integrate(f, 2 * n, a, b);
        }

        return result;
    }

    complex<double> integrate(const std::function<complex<double> (const double &)> & f, unsigned n, const double & a, const double & b)
    {
        if (n & 0x1)
            n += 1;

        if (n < 16)
            n = 16;

        double h = (b - a) / n;
        std::vector<complex<double>> y;

        for (unsigned k(0) ; k < n + 1 ; ++k)
        {
            y.push_back(f(a + k * h));
        }

        complex<double> Q0 = 0.0, Q1 = 0.0, Q2 = 0.0;
        for (unsigned k(0) ; k < n / 8 ; ++k)
        {
            Q0 += y[8 * k] + 4.0 * y[8 * k + 4] + y[8 * k + 4];
        }
        for (unsigned k(0) ; k < n / 4 ; ++k)
        {
            Q1 += y[4 * k] + 4.0 * y[4 * k + 2] + y[4 * k + 4];
        }
        for (unsigned k(0) ; k < n / 2 ; ++k)
        {
            Q2 += y[2 * k] + 4.0 * y[2 * k + 1] + y[2 * k + 2];
        }

        Q0 = Q0 * h / 3.0 * 4.0;
        Q1 = Q1 * h / 3.0 * 2.0;
        Q2 = Q2 * h / 3.0;

        double denom_r = real(Q0 + Q2 - 2.0 * Q1), denom_i = imag(Q0 + Q2 - 2.0 * Q1);
        double num_r = real(Q2 - Q1), num_i = imag(Q2 - Q1);
        double correction_r = num_r * num_r / denom_r, correction_i = num_i * num_i / denom_i;
        complex<double> result;

        if (isnan(correction_r) || isnan(correction_i))
        {
            result = Q2;
        }
        else if ((abs(correction_r / real(Q2)) < 1.0) && (abs(correction_i / imag(Q2) < 1.0)))
        {
            result = Q2 - complex<double>(correction_r, correction_i);
        }
        else
        {
#if 0
            std::cerr << "Q0 = " << Q0 << std::endl;
            std::cerr << "Q1 = " << Q1 << std::endl;
            std::cerr << "Q2 = " << Q2 << std::endl;
            std::cerr << "Reintegrating with twice the number of data points" << std::endl;
#endif
            result = integrate(f, 2 * n, a, b);
        }

        return result;
    }
}
