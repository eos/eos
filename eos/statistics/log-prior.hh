/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
 * Copyright (c) 2019-2022 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_STATISTICS_PRIOR_HH
#define EOS_GUARD_SRC_STATISTICS_PRIOR_HH 1

#include <eos/statistics/log-prior-fwd.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <vector>
#include <memory>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

namespace eos
{
    /*!
     * Base class for log(prior) distributions.
     *
     * Has a container of subclasses describing independent 1..k dimensional prior distributions.
     * Taken together, they specify the full N dimensional prior. Any actual calculation is done
     * by the subclasses.
     */
    class LogPrior
    {
        protected:
            /// Our associated Parameters object.
            Parameters _parameters;

            /// All parameters for which this prior provides information.
            std::vector<Parameter> _varied_parameters;

        public:
            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param parameters The Parameters objects from which we evaluate.
             */
            LogPrior(const Parameters & parameters);

            virtual ~LogPrior()
            {
            }

            virtual std::string as_string() const = 0;

            /*!
             * Create a clone (independent object) of this LogPrior object,
             *
             * @param parameters The Parameters object from which the clone shall evaluate.
             * @return The clone
             */
            virtual LogPriorPtr clone(const Parameters & parameters) const = 0;
            ///@}


            ///@name Iteration over descriptions
            ///@{
            struct IteratorTag;
            using Iterator = WrappedForwardIterator<IteratorTag, Parameter>;

            Iterator begin();
            Iterator end();
            ///@}

            ///@name Accessors
            ///@{

            /*!
             * Evaluate the natural logarithm of the prior.
             */
            virtual double operator() () const = 0;

            /*!
             * Generate a prior sample from the inverse CDF and a set of generator values.
             *
             * The generator values must have been passed to all Parameter objects via their respective
             * set_generator() method.
             */
            virtual void sample() = 0;

            /*!
             * Compute the vector of cummulative probabilities.
             *
             * The results are stored it in the generator values of the Parameter objects.
             */
            virtual void compute_cdf() = 0;

            /*!
             * Return whether or not this prior is informative.
             */
            virtual bool informative() const = 0;
            ///@}

            ///@name Named constructors for 1D prior distributions with finite support
            ///@{
            static LogPriorPtr Flat(const Parameters & parameters, const std::string & name, const double & min, const double & max);
            static LogPriorPtr CurtailedGauss(const Parameters & parameters, const std::string & name, const double & min, const double & max,
                    const double & lower, const double & central, const double & upper);
            static LogPriorPtr Scale(const Parameters & parameter, const std::string & name, const double & min, const double & max,
                    const double & mu_0, const double & lambda);
            ///@}

            ///@name Named constructors for prior distributions with infinite support
            ///@{
            static LogPriorPtr Gaussian(const Parameters & parameters, const QualifiedName & name, const double & mu, const double & sigma);
            static LogPriorPtr MultivariateGaussian(const Parameters & parameters, const std::vector<QualifiedName> & names,
                    gsl_vector * mean, gsl_matrix * covariance);
            static LogPriorPtr Poisson(const Parameters & parameters, const std::string & name, const double & k);
            static LogPriorPtr Transform(const Parameters & parameters, const std::vector<QualifiedName> & names, const std::vector<double> & shift, const std::vector<std::vector<double>> & transform,
                        const std::vector<double> &  min, const std::vector<double> & max);
            ///@}

            /*!
             * Construct a prior from its string representation.
             *
             * @param parameters The object from which the value of the parameter is retrieved.
             * @param serialization The prior in string form.
             * @return Pointer to the ready made prior.
             */
            static LogPriorPtr Make(const Parameters & parameters, const std::string & serialization);
    };

    extern template class WrappedForwardIterator<LogPrior::IteratorTag, Parameter>;
}

#endif
