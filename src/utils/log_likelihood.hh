/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_LIKELIHOOD_HH
#define EOS_GUARD_SRC_UTILS_LIKELIHOOD_HH 1

#include <src/observable.hh>
#include <src/utils/observable_cache.hh>
#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern.hh>

#include <gsl/gsl_rng.h>

namespace eos
{
    // Forward declarations.
    class LogLikelihoodBlock;

    typedef std::shared_ptr<LogLikelihoodBlock> LogLikelihoodBlockPtr;

    /*!
     * LogLikelihoodBlock models the logarithm of the likelihood for a given
     * number of correlated observables which are independent of all other
     * observables such that total likelihood is just the product of
     * the independent blocks.
     *
     * Access to any LogLikelihoodBlock is coherent, i.e., changes to one object will propagate
     * to every other object copy. To create an independent instance, use clone().
     */
    class LogLikelihoodBlock
    {
        public:
            /// Destructor.
            virtual ~LogLikelihoodBlock() = 0;

            /*!
             * Fix the predictions for fixed parameters within given model.
             * Calculate any normalization constants to fix the sampling distribution
             * in order to speed up repeated calls of sample().
             */
            virtual void prepare_sampling();

            /// Compute the logarithm of the likelihood for this block.
            virtual double evaluate() const = 0;

            /*!
             * Sample from the logarithm of the likelihood for this block.
             * @warning Call prepare_sampling() before a call to sample() to
             * ensure that one really gets a log likelihood value from the correct distribution,
             * i.e.  llh ~ llh(model, fixed parameters)
             */
            virtual double sample(gsl_rng * rng) const = 0;

            /// Clone this block.
            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const = 0;

            /*!
             * Create a new LogLikelihoodBlock for one normally distributed observable.
             *
             * @param cache      The Observable cache from which we draw the predictions.
             * @param observable The Observable whose distribution we model.
             * @param min        The value one sigma below the mean of the experimental distribution.
             * @param central    The mean value of the experimental distribution.
             * @param max        The value one sigma above the mean of the experimental distribution.
             */
            static LogLikelihoodBlockPtr Gaussian(ObservableCache cache, const ObservablePtr & observable,
                    const double & min, const double & central, const double & max);
    };

    /*!
     * LogLikelihood handles a set of ObservablePtr with associated measurement data.
     *
     * Access to any LogLikelihood is coherent, i.e., changes to one object will propagate
     * to every other object copy. To create an independent instance, use clone().
     */
    class LogLikelihood :
        public PrivateImplementationPattern<LogLikelihood>
    {
        public:

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param parameters  The Parameters object to which all further ObservablePtr objects must be bound.
             */
            LogLikelihood(const Parameters & parameters);

            /*!
             * Destructor.
             */
            ~LogLikelihood();
            ///@}

            ///@name Access and evaluation
            ///@{
            /*!
             * Add an observable and its associated measurement.
             *
             * @param observable  The Observable that shall be added to the calculation of the likelihood.
             * @param min         The lower bound on the measurement.
             * @param central     The central value of the measurement.
             * @param max         The upper bound on the measurement.
             */
            void add(const ObservablePtr & observable, const double & min,
                    const double & central, const double & max);

            /*!
             * Calculate a p-value based on the \chi^2
             * test statistic for the current setting of the parameters.
             * @note   The p-value is _not_ corrected for degrees of freedom.
             * @param  datasets The number of simulated data sets
             * @return <p-value, uncertainty>, where the uncertainty is
             * estimated from the standard posterior for a Bernoulli experiment.
             */
            std::pair<double, double>
            bootstrap_p_value(const unsigned & datasets);

            /*!
             * Create an independent instance of this LogLikelihood that uses the same set of observables and measurements.
             */
            LogLikelihood clone() const;

            /*!
             * The number of independent observations used in the likelihood.
             * @note This may differ from the number of observables in case
             * two experiments reported their results on the same observable.
             */
            unsigned number_of_observations() const;

            /*!
             * Retrieve the underlying Parameters object.
             */
            Parameters parameters() const;

            /*!
             * Retrieve the cache of observables associated with this LogLikelihood.
             */
            ObservableCache observable_cache() const;

            /*!
             * Evaluate the log likelihood, i.e., return @f[ \log \mathcal{L} = \log P(D | \vec{\theta}, M)=  - \frac{\chi^2}{2} + C@f].
             * @note: all observables are recalculated
             */
            double operator()();

            /*!
             * Evaluate the likelihood, but recalculate only those observables which depend on
             * the parameter  with the given id
             *
             * @warning Do not use this call before __all__ observables have been added to the likelihood!
             * Since the dependence is checked only the first time a parameter id is encountered, adding
             * an observable later will not establish the link between the id and the observable, leading
             * to incorrect values of the likelihood
             */
            double operator()(const Parameter::Id& id);

            /*!
             * Reload previous observable values.
             *
             * @note this is only useful in conjunction with operator()(const Parameter& par)
             * and the Markov Chain sampler
             */
            void reset();
            ///@}
    };
}

#endif
