/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <eos/utils/integrate.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/matrix.hh>

#include <limits>
#include <vector>

namespace eos
{
    static unsigned default_number_of_integration_points = 64;

    using std::abs;
    using std::real;
    using std::imag;

    using integration::validate;

    double
    integrate(const std::function<double (const double &)> & f, const double & a, const double & b, unsigned n)
    {
        validate(n);

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

        if (std::isnan(correction))
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
            result = integrate(f, a, b, 2 * n);
        }

        return result;
    }

    double
    integrate(const std::function<double (const double &)> & f, const double & a, const double & b)
    {
        return integrate(f, a, b, default_number_of_integration_points);
    }

    complex<double>
    integrate(const std::function<complex<double> (const double &)> & f, const double & a, const double & b, unsigned n)
    {
        validate(n);

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

        if (std::isnan(correction_r) || std::isnan(correction_i))
        {
            result = Q2;
        }
        else if ((abs(correction_r / real(Q2)) < 1.0) && (abs(correction_i / imag(Q2)) < 1.0))
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
            result = integrate(f, a, b, 2 * n);
        }

        return result;
    }

    complex<double>
    integrate(const std::function<complex<double> (const double &)> & f, const double & a, const double & b)
    {
        return integrate(f, a, b, default_number_of_integration_points);
    }

    namespace integration
    {
        InvalidNumberOfEvaluations::InvalidNumberOfEvaluations(unsigned n) :
            Exception("Parameter 'n' must be a power of two and at least 16 (got " + stringify(n) + ")")
        {
        }

        void
        set_n(unsigned n)
        {
            validate(n);
            default_number_of_integration_points = n;
        }

        unsigned
        get_n()
        {
            return default_number_of_integration_points;
        }

        void validate(unsigned n)
        {
            if ( n < 16 || (n & (n - 1)) ) // n & (n - 1) is 0 (false) only if n is a power of two.
                throw InvalidNumberOfEvaluations(n);
        }
    }
}
