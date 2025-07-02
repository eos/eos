/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2018 Frederik Beaujean
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

#ifndef EOS_GUARD_EOS_MATHS_INTEGRATE_IMPL_HH
#define EOS_GUARD_EOS_MATHS_INTEGRATE_IMPL_HH 1

#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-cubature.hh>
#include <eos/maths/matrix.hh>

#include <cassert>
#include <vector>

namespace eos
{
    template <std::size_t k> std::array<double, k> integrate1D(const std::function<std::array<double, k> (const double &)> & f, unsigned n, const double & a, const double & b)
    {
        if (n & 0x1)
            n += 1;

        if (n < 16)
            n = 16;

        // step width
        double h = (b - a) / n;

        // evaluate function for every sampling point
        std::vector<std::array<double, k>> y;
        for (unsigned i = 0 ; i < n + 1 ; ++i)
        {
            y.push_back(f(a + i * h));
        }

        std::array<double, k> Q0; Q0.fill(0.0);
        std::array<double, k> Q1; Q1.fill(0.0);
        std::array<double, k> Q2; Q2.fill(0.0);

        for (unsigned i = 0 ; i < n / 8 ; ++i)
        {
            Q0 = Q0 + y[8 * i] + 4.0 * y[8 * i + 4] + y[8 * i + 4];
        }
        for (unsigned i = 0 ; i < n / 4 ; ++i)
        {
            Q1 = Q1 + y[4 * i] + 4.0 * y[4 * i + 2] + y[4 * i + 4];
        }
        for (unsigned i = 0 ; i < n / 2 ; ++i)
        {
            Q2 = Q2 + y[2 * i] + 4.0 * y[2 * i + 1] + y[2 * i + 2];
        }

        Q0 = (h / 3.0 * 4.0) * Q0;
        Q1 = (h / 3.0 * 2.0) * Q1;
        Q2 = (h / 3.0) * Q2;

        std::array<double, k> denom = Q0 + Q2 - 2.0 * Q1;
        std::array<double, k> num = Q2 - Q1;
        std::array<double, k> correction = divide(mult(num, num), denom);

        bool correction_valid = true;
        for (unsigned i = 0 ; i < k ; ++i)
        {
            if (std::isnan(correction[i]))
            {
                correction_valid = false;
                break;
            }
        }

        if (!correction_valid)
        {
            return Q2;
        }
        else
        {
            bool correction_small = true;

            for (unsigned i = 0 ; i < k ; ++i)
            {
                if ((abs(correction[i] / Q2[i])) > 1.0)
                {
                    correction_small = false;
                    break;
                }
            }

            if (correction_small || n >= 1 << 10)
            {
                return Q2 - correction;
            }
            else
            {
#if 0
                std::cerr << "Q0 = " << Q0 << std::endl;
                std::cerr << "Q1 = " << Q1 << std::endl;
                std::cerr << "Q2 = " << Q2 << std::endl;
                std::cerr << "Reintegrating with twice the number of data points" << std::endl;
#endif
                return integrate1D(f, 2 * n, a, b);
            }
        }
    }

    namespace cubature
    {

        template <size_t ndim_>
        int scalar_integrand(unsigned ndim , const double *x, void *data,
                      unsigned fdim , double *fval)
        {
            assert(ndim == ndim_);
            assert(fdim == 1);

            auto& f = *static_cast<cubature::fdd<ndim_> *>(data);
            // TODO use std::array_view once available
            std::array<double, ndim_> args;
            std::copy(x, x + ndim_, args.data());
            *fval = f(args);

            return 0;
        }

    }

    template <size_t ndim_>
    double integrate(const cubature::fdd<ndim_> & f,
                     const std::array<double, ndim_> &a,
                     const std::array<double, ndim_> &b,
                     const cubature::Config &config)
    {
        // TODO Support infinite intervals by param trafo? Not for now.
        constexpr unsigned nintegrands = 1;
        double res;
        double err;
        if (hcubature(nintegrands, &cubature::scalar_integrand<ndim_>,
                      &const_cast<cubature::fdd<ndim_>&>(f), ndim_, a.data(), b.data(),
                      config.maxeval(), config.epsabs(), config.epsrel(), ERROR_L2, &res, &err))
        {
            throw IntegrationError("hcubature failed");
        }

        return res;
    }

}

#endif
