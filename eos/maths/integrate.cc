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

    IntegrationError::IntegrationError(const std::string & message) throw () :
        Exception(message)
    {
    }
}
