/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2018 Frederik Beaujean
 * Copyright (c) 2025 Florian Herren
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
    namespace cubature
    {
        template <size_t ndim_, size_t fdim_, typename T_>
        int integrand_wrapper(unsigned ndim , const double *x, void *data,
                      unsigned fdim , double *fval)
        {
            assert(ndim == ndim_);
            auto buffer_size = cubature::integrand_traits<ndim_, fdim_, T_>::buffer_size;
            assert(fdim == buffer_size);

            auto & f = *static_cast<cubature::integrand<ndim_, fdim_, T_> *>(data);
            typename cubature::integrand_traits<ndim_, fdim_, T_>::argument_type arguments;
            cubature::integrand_traits<ndim_, fdim_, T_>::copy_arguments(x, arguments);
            typename cubature::integrand_traits<ndim_, fdim_, T_>::result_type res = f(arguments);
            cubature::integrand_traits<ndim_, fdim_, T_>::copy_result(res, fval);

            return 0;
        }
    }

    template <size_t ndim_, size_t fdim_, typename T_>
    typename cubature::integrand_traits<ndim_, fdim_, T_>::result_type integrate(const cubature::integrand<ndim_, fdim_, T_> & f,
                                                                                 const typename cubature::integrand_traits<ndim_, fdim_, T_>::argument_type & a,
                                                                                 const typename cubature::integrand_traits<ndim_, fdim_, T_>::argument_type & b,
                                                                                 const cubature::Config &config)
    {
        using integrand = cubature::integrand<ndim_, fdim_, T_>;
        using integrand_traits = cubature::integrand_traits<ndim_, fdim_, T_>;
        using cubature::integrand_wrapper;

        constexpr unsigned nintegrands = integrand_traits::buffer_size;
        typename integrand_traits::buffer_type result_buffer;
        typename integrand_traits::buffer_type error_buffer;
        if (hcubature(nintegrands, &integrand_wrapper<ndim_, fdim_, T_>,
                      &const_cast<integrand&>(f), ndim_, integrand_traits::pointer_from_arguments(a),
                      integrand_traits::pointer_from_arguments(b), config.maxeval(), config.epsabs(), config.epsrel(),
                      ERROR_L2, integrand_traits::pointer_from_buffer(result_buffer), integrand_traits::pointer_from_buffer(error_buffer)))
        {
            throw IntegrationError("hcubature failed");
        }

        return integrand_traits::contruct_result(result_buffer);
    }
}

#endif
