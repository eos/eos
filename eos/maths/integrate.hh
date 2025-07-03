/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
 * Copyright (c) 2018 Danny van Dyk and Frederik Beaujean
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

#ifndef EOS_GUARD_EOS_MATHS_INTEGRATE_HH
#define EOS_GUARD_EOS_MATHS_INTEGRATE_HH 1

#include <eos/maths/complex.hh>
#include <eos/utils/exception.hh>

// TODO Didn't manage to forward declare C struct
// struct gsl_integration_workspace;
// using gsl_integration_workspace = gsl_integration_workspace;
// struct gsl_integration_workspace;
#include <gsl/gsl_integration.h>

#include <array>
#include <functional>

namespace eos
{
    /// @{
    /*!
     * Numerically integrate functions of one real-valued parameter.
     *
     * Uses the Delta^2-Rule by Aitkin to refine the result.
     *
     * @param f      Integrand.
     * @param n      Number of evaluations, must be a power of 2.
     * @param a      Lower limit of the domain of integration.
     * @param b      Upper limit of the domain of integration.
     */
    double integrate1D(const std::function<double (const double &)> & f, unsigned n, const double & a, const double & b);
    complex<double> integrate1D(const std::function<complex<double> (const double &)> & f, unsigned n, const double & a, const double & b);

    template <std::size_t k> std::array<double, k> integrate1D(const std::function<std::array<double, k> (const double &)> & f, unsigned n, const double & a, const double & b);
    /// @}

namespace GSL
{
    using fdd = std::function<double(const double &)>;
    struct QNG
    {
        class Config
        {
            public:
                Config();

                double epsabs() const;
                Config& epsabs(const double& x);

                double epsrel() const;
                Config& epsrel(const double& x);
            private:
              double _epsabs, _epsrel;
        };
    };

    struct QAGS
    {
        class Workspace
        {
        public:
            Workspace(int limit = 5000);
            Workspace(const Workspace &) = delete;
            Workspace(Workspace &&) = delete;
            ~Workspace();
            Workspace &operator=(const Workspace &) = delete;
            Workspace &operator=(Workspace &&) = delete;
            operator gsl_integration_workspace *() const { return _work_space; }
            int limit() const { return _work_space->limit; }
        private:
            gsl_integration_workspace *_work_space;
        };

        class Config
        {
            public:
                Config();

                double epsabs() const;
                Config& epsabs(const double& x);

                double epsrel() const;
                Config& epsrel(const double& x);

                int key() const;
                Config& key(const int&);
            private:
                QNG::Config _qng;
                int _key;
        };
    };

    static thread_local QAGS::Workspace work_space;
}

    /*!
     * Numerically integrate functions of one real-valued parameter.
     *
     * Two methods from the
     * GNU scientific library are wrapped:
     * 1) `QNG`: the non-adaptive Gauss-Kronrod rule
     * 2) `QAGS`: the adaptive Clenshaw-Kurtis rule
     */
    template <typename Method_>
    double integrate(const std::function<double(const double &)> & f,
                     const double &a, const double &b,
                     const typename Method_::Config &config = typename Method_::Config());

namespace cubature
{
    template <size_t ndim_, size_t fdim_> struct integrand_traits
    {
        // argument type for the integrand function
        using argument_type = std::array<double, ndim_>;

        // result type for the integrand function
        using result_type = std::array<double, fdim_>;

        // helpers to copy data between C-style arrays and C++ types
        static void copy_arguments(const double * from, argument_type & to) { std::copy(from, from + ndim_, to.begin()); };
        static void copy_result(result_type & from, double * to) { std::copy(from.cbegin(), from.cend(), to); };
        static const double * pointer_from_arguments(const argument_type & arguments) { return arguments.data(); };
        static double * pointer_from_result(result_type & result) { return result.data(); };

        // type of the integrand function
        using function_type = std::function<result_type (const argument_type &)>;
    };

    template <size_t ndim_> struct integrand_traits<ndim_, 1>
    {
        // argument type for the integrand function
        using argument_type = std::array<double, ndim_>;

        // result type for the integrand function
        using result_type = double;

        // helpers to copy data between C-style arrays and C++ types
        static void copy_arguments(const double * from, argument_type & to) { std::copy(from, from + ndim_, to.data()); };
        static void copy_result(result_type & from, double * to) { *to = from; };
        static double * pointer_from_result(result_type & result) { return &result; };
        static const double * pointer_from_arguments(const argument_type & arguments) { return arguments.data(); };

        // function type of the integrand function
        using function_type = std::function<result_type (const argument_type &)>;
    };

    template <size_t fdim_> struct integrand_traits<1, fdim_>
    {
        // argument type for the integrand function
        using argument_type = double;

        // result type for the integrand function
        using result_type = std::array<double, fdim_>;

        // helpers to copy data between C-style arrays and C++ types
        static void copy_arguments(const double * from, argument_type & to) { to = *from; };
        static void copy_result(result_type & from, double * to)
        {
            for (unsigned i = 0 ; i < fdim_ ; i++)
            {
                to[i] = from[i];
            }
        };
        static const double * pointer_from_arguments(const argument_type & arguments) { return &arguments; };
        static double * pointer_from_result(result_type & result) { return result.data(); };

        // function type of the integrand function
        using function_type = std::function<result_type (const argument_type &)>;
    };

    template <> struct integrand_traits<1, 1>
    {
        // argument type for the integrand function
        using argument_type = double;

        // result type for the integrand function
        using result_type = double;

        // helpers to copy data between C-style arrays and C++ types
        static void copy_arguments(const double * from, argument_type & to) { to = *from; };
        static void copy_result(result_type & from, double * to) { *to = from; };
        static const double * pointer_from_arguments(const argument_type & arguments) { return &arguments; };
        static double * pointer_from_result(result_type & result) { return &result; };

        // function type of the integrand function
        using function_type = std::function<result_type (const argument_type &)>;
    };

    template <size_t ndim_, size_t fdim_ = 1> using integrand = typename integrand_traits<ndim_, fdim_>::function_type;


    class Config
    {
    public:
        Config();

        double epsabs() const;
        Config& epsabs(const double& x);

        double epsrel() const;
        Config& epsrel(const double& x);

        size_t maxeval() const;
        Config& maxeval(const size_t& x);
    private:
        GSL::QNG::Config _qng;
        size_t _maxeval;
    };
}

    /*!
     * Numerically integrate functions of one or more than one variable with
     * cubature methods.
     */
    template <size_t ndim_, size_t fdim_ = 1>
    typename cubature::integrand_traits<ndim_, fdim_>::result_type integrate(const cubature::integrand<ndim_, fdim_> & f,
                                                                             const typename cubature::integrand_traits<ndim_, fdim_>::argument_type & a,
                                                                             const typename cubature::integrand_traits<ndim_, fdim_>::argument_type & b,
                                                                             const cubature::Config &config = cubature::Config());


    class IntegrationError :
        public Exception
    {
        public:
            IntegrationError(const std::string & message) throw ();
    };

}

#endif
