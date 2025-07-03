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
    template <size_t ndim_>
    using fdd = std::function<double(const std::array<double, ndim_> &)>;

    template <size_t ndim_, size_t fdim_>
    using fdd_v = std::function<std::array<double, fdim_>(const std::array<double, ndim_> &)>;

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
    template <size_t ndim_>
    double integrate(const std::function<double(const std::array<double, ndim_> &)> & f,
                     const std::array<double, ndim_> &a,
                     const std::array<double, ndim_> &b,
                     const cubature::Config &config = cubature::Config());

    template <size_t ndim_, size_t fdim_>
    std::array<double, fdim_> integrate(const std::function<std::array<double, fdim_>(const std::array<double, ndim_> &)> & f,
                                        const std::array<double, ndim_> &a,
                                        const std::array<double, ndim_> &b,
                                        const cubature::Config &config = cubature::Config());

    class IntegrationError :
        public Exception
    {
        public:
            IntegrationError(const std::string & message) throw ();
    };

}

#endif
