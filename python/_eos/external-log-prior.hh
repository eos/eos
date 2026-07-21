/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#include "eos/statistics/log-prior.hh"

#include <boost/python.hpp>

#ifndef EOS_PYTHON__EOS_EXTERNAL_LOG_PRIOR_HH
#  define EOS_PYTHON__EOS_EXTERNAL_LOG_PRIOR_HH 1

namespace eos
{
    class ExternalLogPrior : public LogPrior
    {
        private:
            boost::python::object _factory;
            boost::python::object _python_prior;
            boost::python::object _evaluate;
            boost::python::object _sample;
            boost::python::object _compute_cdf;
            bool                  _informative;

        public:
            ExternalLogPrior(const Parameters & parameters, boost::python::object factory);

            ~ExternalLogPrior();

            static LogPriorPtr make(const Parameters & parameters, boost::python::object factory);

            virtual std::string as_string() const override;

            virtual LogPriorPtr clone(const Parameters & parameters) const override;

            virtual double operator() () const override;

            virtual void sample() override;

            virtual void compute_cdf() override;

            virtual bool informative() const override;
    };
} // namespace eos

#endif // EOS_PYTHON__EOS_EXTERNAL_LOG_PRIOR_HH
