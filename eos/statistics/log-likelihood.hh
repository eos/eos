/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2014, 2017 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_STATISTICS_LIKELIHOOD_HH
#define EOS_GUARD_SRC_STATISTICS_LIKELIHOOD_HH 1

#include <eos/constraint.hh>
#include <eos/observable.hh>
#include <eos/statistics/log-likelihood-fwd.hh>
#include <eos/statistics/test-statistic.hh>
#include <eos/maths/matrix.hh>
#include <eos/utils/observable_cache.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <cmath>

namespace eos
{
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

            virtual std::string as_string() const = 0;

            /// Clone this block.
            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const = 0;

            /// Compute the logarithm of the likelihood for this block.
            virtual double evaluate() const = 0;

            /// The number of experimental observations (not observables!) used in this block.
            virtual unsigned number_of_observations() const = 0;

            /*!
             * Sample from the logarithm of the likelihood for this block.
             * @warning Call prepare_sampling() before a call to sample() to
             * ensure that one really gets a log likelihood value from the correct distribution,
             * i.e.  llh ~ llh(model, fixed parameters)
             *
             * @param rng The random number generator.
             */
            virtual double sample(gsl_rng * rng) const = 0;

            /*!
             * Calculate the significance of the deviation between
             * the observables' current value and the mode in
             * units of the standard Gaussian distribution.
             *
             * @example For a Gaussian around x=1 with sigma = 0.5,
             * a current value of x=2 will yield a significance of 2.
             *
             * @note The significance is >= 0.
             */
            virtual double significance() const  = 0;

            /*!
             * Calculate the primary test statistic (e.g. a chi^2 or -2 log L) for
             * the deviation between the observables' current values and
             * the mode.
             */
            virtual TestStatistic primary_test_statistic() const = 0;

            /*!
             * Create a new LogLikelihoodBlock for one normally distributed observable.
             *
             * @note For every dimension, this template and the corresponding implementation
             *       have to be instantiated explicitly.
             *
             * @param cache      The Observable cache from which we draw the predictions.
             * @param observable The Observable whose distribution we model.
             * @param min        The value one sigma below the mean of the experimental distribution.
             * @param central    The mean value of the experimental distribution.
             * @param max        The value one sigma above the mean of the experimental distribution.
             */
            static LogLikelihoodBlockPtr Gaussian(ObservableCache cache, const ObservablePtr & observable,
                    const double & min, const double & central, const double & max,
                    const unsigned & number_of_observations = 1u);

            /*!
             * Create a new LogLikelihoodBlock for one a single observable with asymmetric uncertainties.
             *
             * By construction, it is constructed to behave similarly to a Gaussian:
             * - the mode is at the central value
             * - the interval [min, max] contains 68% probability
             * - the density at min is the same as at max
             * For details see [C:2010A].
             *
             * @note Finding the correct parameter values needs to be accomplished by the user, by solving
             * a set of two equations corresponding to the condition above, numerically. This is unstable
             * if max and min uncertainty differ by less than 5%. In that case, just a Gaussian instead.
             * However, the consistency of (alpha, lambda) and (min, mode, max) is checked.
             *
             * @param cache      The Observable cache from which we draw the predictions.
             * @param observable The Observable whose distribution we model.
             * @param min        The value one sigma below the mean of the experimental distribution.
             * @param central    The mean value of the experimental distribution.
             * @param max        The value one sigma above the mean of the experimental distribution.
             * @param alpha      The shape parameter of a LogGamma distribution.
             * @param lambda     The scale parameter of a LogGamma distribution.
             */
            static LogLikelihoodBlockPtr LogGamma(ObservableCache cache, const ObservablePtr & observable,
                    const double & min, const double & central, const double & max,
                    const double & alpha, const double & lambda,
                    const unsigned & number_of_observations = 1u);

            /*!
             * A likelihood contribution representing an upper limit on a quantity x.
             *
             * Internally, it is represented by an Amoroso distribution [C:2010A] with
             * location parameter a set to the physical limit, scale parameter \theta and
             * shape parameters \alpha, \beta supplied by the user.
             *
             * @note The mode of the distribution is typically not at the physical limit.
             *
             * The limit values are required in order to check consistency with the parameter values.
             *
             * For example, consider a yet unobserved branching ratio, which has to be non-negative, x>= 0.
             *
             * @param cache The Observable cache from which we draw the predictions.
             * @param observable The Observable whose distribution we model.
             * @param mode The lower, physical limit. A branching ratio has to be >= 0.
             * @param upper_limit_10  With  10% probability, x < upper_limit_10.
             * @param upper_limit_50 With 50% probability, x < upper_limit_50.
             * @param upper_limit_90 With 90% probability, x < upper_limit_90.
             * @param theta scale parameter
             * @param alpha 1st shape parameter
             * @param beta  2nd shape parameter
             * @return
             */
            static LogLikelihoodBlockPtr Amoroso(ObservableCache cache, const ObservablePtr & observable,
                    const double & physical_limit, const double & upper_limit_10,
                    const double & upper_limit_50, const double & upper_limit_90,
                    const double & theta, const double & alpha, const double & beta,
                    const unsigned & number_of_observations = 1u);

            /*!
             * A likelihood contribution representing an upper limit on a quantity x.
             *
             * Internally, it is represented by an Amoroso distribution [C:2010A] with
             * location parameter a set to the physical limit, scale parameter \theta and
             * shape parameters \alpha, \beta supplied by the user.
             *
             * @note The mode of the distribution is typically not at the physical limit.
             *
             * For example, consider a yet unobserved branching ratio, which has to be non-negative, x>= 0.
             *
             * @param cache The Observable cache from which we draw the predictions.
             * @param observable The Observable whose distribution we model.
             * @param mode The lower, physical limit. A branching ratio has to be >= 0.
             * @param theta scale parameter
             * @param alpha 1st shape parameter
             * @param beta  2nd shape parameter
             * @return
             */
            static LogLikelihoodBlockPtr Amoroso(ObservableCache cache, const ObservablePtr & observable,
                    const double & physical_limit, const double & theta, const double & alpha, const double & beta,
                    const unsigned & number_of_observations = 1u);

            // todo document
            static LogLikelihoodBlockPtr Mixture(const std::vector<LogLikelihoodBlockPtr> & components,
                                                 const std::vector<double> & weights,
                                                 const std::vector<std::array<double, 2>> & test_stat);

            /*!
             * Create a new LogLikelihoodBlock for n observables distributed
             * according to a multivariate normal distribution.
             *
             * @param cache         The Observable cache from which we draw the predictions.
             * @param observables   The Observables whose distribution we model.
             * @param mean          The vector of means.
             * @param covariance    The covariance matrix
             */
            static LogLikelihoodBlockPtr MultivariateGaussian(ObservableCache cache, const std::vector<ObservablePtr> & observables,
                    gsl_vector * mean, gsl_matrix * covariance, gsl_matrix * response, const unsigned & number_of_observations);

            /*!
             * Templated helper to create a new LogLikelihoodBlock for n observables distributed
             * according to a multivariate normal distribution.
             *
             * @param cache         The Observable cache from which we draw the predictions.
             * @param observables   The Observables whose distribution we model.
             * @param mean          The vector of means.
             * @param covariance    The covariance matrix
             */
            template <std::size_t n_>
            static LogLikelihoodBlockPtr MultivariateGaussian(ObservableCache cache, const std::array<ObservablePtr, n_> & observables,
                    const std::array<double, n_> & mean, const std::array<std::array<double, n_>, n_> & covariance,
                    const unsigned & number_of_observations = n_)
            {
                std::vector<ObservablePtr> _observables(observables.begin(), observables.end());

                // create GSL vector for the mean
                gsl_vector * _mean = gsl_vector_alloc(n_);
                for (auto i = 0u ; i < n_ ; ++i)
                {
                    gsl_vector_set(_mean, i, mean[i]);
                }

                // create GSL matrix for the covariance
                gsl_matrix * _covariance = gsl_matrix_alloc(n_, n_);
                for (auto i = 0u ; i < n_ ; ++i)
                {
                    for (auto j = 0u ; j < n_ ; ++j)
                    {
                        gsl_matrix_set(_covariance, i, j, covariance[i][j]);
                    }
                }

                // create GSL matrix for the response matrix
                gsl_matrix * _response = gsl_matrix_calloc(n_, n_);
                gsl_matrix_set_identity(_response);

                return LogLikelihoodBlock::MultivariateGaussian(cache, _observables, _mean, _covariance, _response, number_of_observations);
            }

            /*!
             * Templated helper to create a new LogLikelihoodBlock for n observables distributed
             * according to a multivariate normal distribution.
             *
             * @param cache         The Observable cache from which we draw the predictions.
             * @param observables   The Observables whose distribution we model.
             * @param mean          The vector of means.
             * @param variances     The vector of variances.
             * @param correlation   The correlation matrix. Diagonal is assumed to be one.
             */
            template <std::size_t n_>
            static LogLikelihoodBlockPtr MultivariateGaussian(ObservableCache cache, const std::array<ObservablePtr, n_> & observables,
                    const std::array<double, n_> & mean, const std::array<double, n_> & variances,
                    const std::array<std::array<double, n_>, n_> & correlation,
                    const unsigned & number_of_observations = n_)
            {
                std::vector<ObservablePtr> _observables(observables.begin(), observables.end());

                // create GSL vector for the mean
                gsl_vector * _mean = gsl_vector_alloc(n_);
                for (auto i = 0u ; i < n_ ; ++i)
                {
                    gsl_vector_set(_mean, i, mean[i]);
                }

                // create GSL matrix for the covariance
                gsl_matrix * covariance = gsl_matrix_alloc(n_, n_);
                for (auto i = 0u ; i < n_ ; ++i)
                {
                    for (auto j = 0u ; j < n_ ; ++j)
                    {
                        double value = std::sqrt(variances[i] * variances[j]) * correlation[i][j];
                        gsl_matrix_set(covariance, i, j, value);
                    }
                }

                // create GSL matrix for the response matrix
                gsl_matrix * _response = gsl_matrix_calloc(n_, n_);
                gsl_matrix_set_identity(_response);

                return LogLikelihoodBlock::MultivariateGaussian(cache, _observables, _mean, covariance, _response, number_of_observations);
            }

            /*!
             * Create a new LogLikelihoodBlock for a uniform bound by summing individual contributions and applying a gaussian-like penalty
             * beyond the bound.
             *
             * @param observables The pseudo observables implementing the contributions to the bound.
             * @param bound The value of the bound.
             * @param uncetainty The uncertainty entering the gaussian-like penalty.
             */
            static LogLikelihoodBlockPtr UniformBound(ObservableCache cache, const std::vector<ObservablePtr> & observables,
                                                      const double & bound, const double & uncertainty);
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
             * @param observable                The Observable that shall be added to the calculation of the likelihood.
             * @param min                       The lower bound on the measurement.
             * @param central                   The central value of the measurement.
             * @param max                       The upper bound on the measurement.
             * @param number_of_observations    The number of observations associated with the measurement. Defaults to 1.
             */
            void add(const ObservablePtr & observable, const double & min,
                    const double & central, const double & max, const unsigned & number_of_observations = 1u);

            /*!
             * Add one of the libraries experimental constraints.
             *
             * @param constraint The experimental constraint, cf. Constraint::make
             */
            void add(const Constraint & constraint);

            /*!
             * Add an externally-defined log likelihood block.
             *
             * @param block The log likelihood block.
             */
            void add(const LogLikelihoodBlockPtr & block);

            ///@name Iteration and Access
            ///@{
            struct ConstraintIteratorTag;
            using ConstraintIterator = WrappedForwardIterator<ConstraintIteratorTag, Constraint>;

            /// Iterator to the first constraint.
            ConstraintIterator begin() const;

            /// Iterator pointing past the last constraint.
            ConstraintIterator end() const;
            ///@}

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
            double operator()() const;
            ///@}
    };

    extern template class WrappedForwardIterator<LogLikelihood::ConstraintIteratorTag, Constraint>;
}

#endif
