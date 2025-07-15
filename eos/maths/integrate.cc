/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/matrix.hh>

#include <gsl/gsl_errno.h>

#include <limits>
#include <vector>

namespace
{
    double gsl_function_adapter(double x, void *params)
    {
        const auto& f = *static_cast<eos::GSL::fdd*>(params);
        return f(x);
    }
}

namespace eos
{
    using std::abs;
    using std::real;
    using std::imag;

    double integrate1D(const std::function<double (const double &)> & f, unsigned n, const double & a, const double & b)
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
            result = integrate1D(f, 2 * n, a, b);
        }

        return result;
    }

    complex<double> integrate1D(const std::function<complex<double> (const double &)> & f, unsigned n, const double & a, const double & b)
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
            result = integrate1D(f, 2 * n, a, b);
        }

        return result;
    }

    namespace GSL
    {
        QNG::Config::Config() :
            _epsabs(0),
            _epsrel(1e-4)
        {
        }

        double QNG::Config::epsabs() const
        {
            return _epsabs;
        }

        QNG::Config & QNG::Config::epsabs(const double &x)
        {
            _epsabs = x;
            return *this;
        }

        double QNG::Config::epsrel() const
        {
            return _epsrel;
        }

        QNG::Config & QNG::Config::epsrel(const double &x)
        {
            _epsrel = x;
            return *this;
        }

        QAGS::Workspace::Workspace(int size) :
            _work_space(gsl_integration_workspace_alloc(size))
        {
        }

        QAGS::Workspace::~Workspace()
        {
            gsl_integration_workspace_free(_work_space);
        }

        QAGS::Config::Config() :
            _qng(),
            _key(2)
        {
        }

        double QAGS::Config::epsabs() const
        {
            return _qng.epsabs();
        }

        QAGS::Config & QAGS::Config::epsabs(const double &x)
        {
            _qng.epsabs(x);
            return *this;
        }

        double QAGS::Config::epsrel() const
        {
            return _qng.epsrel();
        }

        QAGS::Config & QAGS::Config::epsrel(const double &x)
        {
            _qng.epsrel(x);
            return *this;
        }

        int QAGS::Config::key() const
        {
            return _key;
        }

        QAGS::Config & QAGS::Config::key(const int & x)
        {
            _key = x;
            return *this;
        }
    }

    template <>
    double integrate<GSL::QNG>(const GSL::fdd &f, const double &a, const double &b, const GSL::QNG::Config &config)
    {
        double result, abserr;
        size_t neval;
        gsl_function F;
        F.function = &gsl_function_adapter;
        F.params = (void*)&f;

        auto status = gsl_integration_qng(&F, a, b, config.epsabs(), config.epsrel(),
                                          &result, &abserr, &neval);

        if (status)
        {
            throw IntegrationError(gsl_strerror(status));
        }
        return result;
    }

    template <>
    double integrate<GSL::QAGS>(const GSL::fdd &f, const double &a, const double &b, const GSL::QAGS::Config &config)
    {
        double result, abserr;
        gsl_function F;
        F.function = &gsl_function_adapter;
        F.params = (void*)&f;

        auto status = gsl_integration_qag(&F, a, b, config.epsabs(), config.epsrel(),
                                          GSL::work_space.limit(), config.key(),
                                          GSL::work_space,
                                          &result, &abserr);

        if (status)
        {
            throw IntegrationError(gsl_strerror(status));
        }
        return result;
    }

    namespace cubature
    {
        Config::Config() :
            _qng(),
            _maxeval(50000)
        {
        }

        double Config::epsabs() const
        {
            return _qng.epsabs();
        }

        Config & Config::epsabs(const double &x)
        {
            _qng.epsabs(x);
            return *this;
        }

        double Config::epsrel() const
        {
            return _qng.epsrel();
        }

        Config & Config::epsrel(const double &x)
        {
            _qng.epsrel(x);
            return *this;
        }

        size_t Config::maxeval() const
        {
            return _maxeval;
        }

        Config & Config::maxeval(const size_t& x)
        {
            _maxeval = x;
            return *this;
        }
    }

    double integrate(const cubature::fdd_s_s & f,
                     const double &a,
                     const double &b,
                     const cubature::Config &config)
    {
        constexpr unsigned nintegrands = 1;
        double res;
        double err;
        if (hcubature(nintegrands, &cubature::scalar_integrand<1>,
                      &const_cast<cubature::fdd_s_s&>(f), 1, &a, &b,
                      config.maxeval(), config.epsabs(), config.epsrel(), ERROR_L2, &res, &err))
        {
            throw IntegrationError("hcubature failed");
        }

        return res;
    }

    complex<double> integrate(const cubature::fdd_s_s_c & f,
                              const double &a,
                              const double &b,
                              const cubature::Config &config)
    {
        constexpr unsigned nintegrands = 2;
        std::array<double, 2> res;
        double err;
        if (hcubature(nintegrands, &cubature::complex_integrand<1>,
                      &const_cast<cubature::fdd_s_s_c&>(f), 1, &a, &b,
                      config.maxeval(), config.epsabs(), config.epsrel(), ERROR_L2, res.data(), &err))
        {
            throw IntegrationError("hcubature failed");
        }

        return complex<double>(res[0], res[1]);
    }

    IntegrationError::IntegrationError(const std::string & message) throw () :
        Exception(message)
    {
    }
}
