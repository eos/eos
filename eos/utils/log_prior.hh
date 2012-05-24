/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
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

#ifndef EOS_GUARD_SRC_UTILS_PRIOR_HH
#define EOS_GUARD_SRC_UTILS_PRIOR_HH 1

#include <eos/utils/log_prior-fwd.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <vector>
#include <memory>

#include  <gsl/gsl_rng.h>

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
            std::vector<ParameterDescription> _parameter_descriptions;

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
            typedef WrappedForwardIterator<IteratorTag, ParameterDescription> Iterator;

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
             * @param rng the random number engine
             * @return a sample according to this prior distribution
             */
            virtual double sample(gsl_rng * rng) const = 0;

            /*!
             * Return the mean of the distribution.
             */
            virtual double mean() const = 0;

            /*!
             * Return the variance of this distribution.
             */
            virtual double variance() const = 0;
            ///@}

            ///@name Named constructors for 1D prior distributions
            ///@{
            static LogPriorPtr Discrete(const Parameters & parameters, const std::string & name, const std::set<double> & values);
            static LogPriorPtr Flat(const Parameters & parameters, const std::string & name, const ParameterRange & range);
            static LogPriorPtr Gauss(const Parameters & parameters, const std::string & name, const ParameterRange & range,
                    const double & lower, const double & central, const double & upper);

            /*!
             * The LogGamma distribution is useful to model continuous, unimodal, asymmetric priors in one dimension.
             *
             * @note The construction will typically fail if the asymmetry is less than ca. 5 %. Use a Gaussian instead.
             *
             * By construction, it is constructed to behave similarly to a Gaussian:
             * - the mode is at the central value
             * - the interval [lower, upper] contains 68% probability
             * - the density at lower is the same as at upper
             *
             * @param parameters The object from which the value of the parameter is retrieved.
             * @param name The parameter name.
             * @param range Total allowed range. The value of the pdf is rescaled such that the prior integrates to one over the range.
             *              For efficiency reasons, no checking if input parameter is in range. But conceptually, the prior is zero outside.
             * @param lower (central - lower) roughly corresponds to a Gaussian standard deviation for values left of the mode.
             * @param central The mode of the prior.
             * @param upper (upper - central) roughly corresponds to a Gaussian standard deviation for values right of the mode.
             * @return
             */
            static LogPriorPtr LogGamma(const Parameters & parameters, const std::string & name, const ParameterRange & range,
                    const double & lower, const double & central, const double & upper);

            /*!
             * Construct a prior from its string representation.
             *
             * @param parameters The object from which the value of the parameter is retrieved.
             * @param serialization The prior in string form.
             * @return Pointer to the ready made prior.
             */
            static LogPriorPtr Make(const Parameters & parameters, const std::string & serialization);
            ///@}
    };
}

#endif
