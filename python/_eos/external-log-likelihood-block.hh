/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2023-2025 Danny van Dyk
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

#include "eos/statistics/log-likelihood.hh"

#include <boost/python.hpp>

#ifndef EOS_PYTHON__EOS_EXTERNAL_LOG_LIKELIHOOD_BLOCK_HH
#  define EOS_PYTHON__EOS_EXTERNAL_LOG_LIKELIHOOD_BLOCK_HH 1

namespace eos
{
    class ExternalLogLikelihoodBlock : public LogLikelihoodBlock
    {
        private:
            ObservableCache       _cache;
            boost::python::object _factory;
            boost::python::object _python_llh_block;
            boost::python::object _evaluate;
            unsigned              _number_of_observations;

        public:
            ExternalLogLikelihoodBlock(const ObservableCache & cache, boost::python::object factory);

            ~ExternalLogLikelihoodBlock();

            static LogLikelihoodBlockPtr make(const ObservableCache & cache, boost::python::object factory);

            virtual std::string as_string() const override;

            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const override;

            virtual double evaluate() const override;

            virtual unsigned int number_of_observations() const override;

            virtual double sample(gsl_rng *) const override;

            virtual double significance() const override;

            virtual TestStatistic primary_test_statistic() const override;
    };
} // namespace eos

#endif // EOS_PYTHON__EOS_EXTERNAL_LOG_LIKELIHOOD_BLOCK_HH
