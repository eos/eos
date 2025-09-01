/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013-2019 Danny van Dyk
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

#include <eos/statistics/log-likelihood.hh>
#include <eos/statistics/test-statistic-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/observable_cache.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/verify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_vector.h>

#include <config.h>

#ifdef EOS_USE_GSL_LINALG_CHOLESKY_DECOMP
#  if (EOS_USE_GSL_LINALG_CHOLESKY_DECOMP == 1)
#    define GSL_LINALG_CHOLESKY_DECOMP gsl_linalg_cholesky_decomp
#  else
#    define GSL_LINALG_CHOLESKY_DECOMP gsl_linalg_cholesky_decomp1
#  endif
#else
#  error EOS_USE_GSL_LINALG_CHOLESKY_DECOMP not defined.
#endif

namespace eos
{
    namespace implementation
    {
        struct GaussianBlock :
            public LogLikelihoodBlock
        {
            ObservableCache cache;

            ObservableCache::Id id;

            const double mode, sigma_lower, sigma_upper;

            // coefficients needed for asymmetric Gaussian x^{+a}_{-b}
            // the pdf/cumulative is a piecewise function
            // CDF(x) = CDF_lower(x, sigma_lower) if x < central, else CDF_upper(x, sigma_upper)
            // To ensure that cumulative is
            // a) continuous at the central value
            // b) normalized to one
            // need to fix the coefficients in
            // P(y|x, a,b) = c_a \mathcal{N}(y|x,a) \theta(y-x) + c_b \mathcal{N}(y|x,b) \theta(x-y)
            // to give
            // c_a = 2 \frac{a}{a+b}, c_b = 2 \frac{b}{a+b}
            const double c_upper, c_lower;
            const double norm;

            const unsigned _number_of_observations;

            GaussianBlock(const ObservableCache & cache, ObservableCache::Id id, const double & min, const double & central, const double & max,
                    const unsigned & number_of_observations) :
                cache(cache),
                id(id),
                mode(central),
                sigma_lower(central - min),
                sigma_upper(max - central),
                c_upper(2.0 * sigma_upper / (sigma_upper + sigma_lower)),
                c_lower(sigma_lower / sigma_upper * c_upper),
                norm(std::log(std::sqrt(2.0/M_PI) / (sigma_upper + sigma_lower))),
                _number_of_observations(number_of_observations)
            {
            }

            virtual ~GaussianBlock()
            {
            }

            virtual std::string as_string() const
            {
                std::string result = "Gaussian: ";
                result += stringify(mode);
                if (sigma_upper == sigma_lower)
                {
                    result += " +- " + stringify(sigma_upper);
                }
                else
                {
                    result += " + " + stringify(sigma_upper) + " - " + stringify(sigma_lower);
                }

                if (0 == _number_of_observations)
                    result += "; no observation";

                return result;
            }

            virtual double evaluate() const
            {
                double value = cache[id];
                double sigma = 0.0;

                // allow for asymmetric Gaussian uncertainty
                if (value > mode)
                    sigma = sigma_upper;
                else
                    sigma = sigma_lower;

                double chi = (value - mode) / sigma;

                return norm - power_of<2>(chi) / 2.0;
            }

            virtual unsigned number_of_observations() const
            {
                return _number_of_observations;
            }

            /*!
             * Mirror and shift the experimental distribution.
             *
             * Why shift? we want to generate toy data for fixed
             * theory. We don't have a full forward model, so we have
             * to make an ad-hoc assumption. We choose the theory
             * prediction as the new most likely value, and take over
             * the uncertainties from experiment.
             *
             * Why mirror? If sigma_upper >> sigma_lower, and theory
             * value > mode, then the theory is in the slowly falling
             * tail.  If you flip the role of theory and exp., then a
             * theory value that is likely under exp. should yield
             * likely value of exp. assuming theory.
             *
             * This procedure is used in both sample() and significance()
             */

            virtual double sample(gsl_rng * rng) const
            {
                // find out if sample in upper or lower part
                double u = gsl_rng_uniform(rng);

                // mirror and shift the distribution
                const double & c_b = c_upper;
                const double & a = sigma_lower, & b = sigma_upper;

                // fixed theory prediction
                const double & theory = cache[id];

                // get a sample observable using the inverse-transform method
                double obs, sigma;
                if (u < b / (a + b))
                {
                    obs = gsl_cdf_gaussian_Pinv(u / c_b, b) + theory;
                    sigma = b;
                }
                else
                {
                    obs = gsl_cdf_gaussian_Pinv(u - 0.5 * c_b, a) + theory;
                    sigma = a;
                }

                // calculate the properly normalized log likelihood
                // note that we generate from theory,
                const double chi = (theory - obs) / sigma;
                return norm - power_of<2>(chi) / 2.0;
            }

            virtual double significance() const
            {
                const double value = cache[id];
                double sigma = 0.0;

                // flip and shift the exp. distribution!
                if (value > mode)
                    sigma = sigma_upper;
                else
                    sigma = sigma_lower;

                // Return positive significance if measured value exceeds predictions.
                // For the Gaussian, there still is 68% probability in [x-b, x+a], even if a != b
                return (mode - value) / sigma;
            }

            virtual TestStatistic primary_test_statistic() const
            {
                double chi = significance();
                return test_statistics::ChiSquare(chi * chi, 1, chi);
            }

            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                ObservablePtr observable = this->cache.observable(id)->clone(cache.parameters());

                return LogLikelihoodBlockPtr(new GaussianBlock(cache, cache.add(observable), mode - sigma_lower, mode, mode + sigma_upper, _number_of_observations));
            }
        };

        // For more details on the LogGamma distribution, see [C:2010A].
        struct LogGammaBlock :
            public LogLikelihoodBlock
        {
            ObservableCache cache;

            ObservableCache::Id id;

            double central, sigma_lower, sigma_upper;

            double nu, lambda, alpha;

            double norm;

            unsigned _number_of_observations;

            LogGammaBlock(const ObservableCache & cache, ObservableCache::Id id,
                          const double & min, const double & central, const double & max,
                          const double & alpha, const double & lambda,
                          const unsigned & number_of_observations) :
                cache(cache),
                id(id),
                central(central),
                sigma_lower(central - min),
                sigma_upper(max - central),
                nu(central - lambda * std::log(alpha)),
                lambda(lambda),
                alpha(alpha),
                _number_of_observations(number_of_observations)
            {
                const double sigma_plus( (sigma_upper > sigma_lower)? sigma_upper / sigma_lower : sigma_lower / sigma_upper);

                if (sigma_plus < 1 + 5e-2)
                {
                    Log::instance()->message("LogLikelihoodBlock::LogGamma.ctor", ll_warning)
                        << "For nearly symmetric uncertainties (" << sigma_lower << " vs " << sigma_upper
                        << "), this procedure may fail to find the correct parameter values. "
                           "Please use a Gaussian block instead.";
                }

                // check consistency
                static const double eps_cdf = 1.0e-4;
                if (std::abs(cdf(central + sigma_upper) - cdf(central - sigma_lower) - 0.68268949213708585 ) > eps_cdf)
                {
                    throw InternalError("LogLikelihoodBlock::LogGamma.ctor: For the current parameter values, "
                                        "the interval [lower, upper] doesn't contain approx. 68%; contents is " + stringify(cdf(central + sigma_upper) - cdf(central - sigma_lower)));

                }
                const double z_plus = (central + sigma_upper - nu) / lambda;
                const double z_minus = (central - sigma_lower - nu) / lambda;
                static const double eps_pdf = 2.5e-2;
                if (std::fabs(alpha * z_plus - std::exp(z_plus) - alpha * z_minus + std::exp(z_minus)) > eps_pdf)
                {
                    throw InternalError("LogLikelihoodBlock::LogGamma.ctor: For the current parameter values, "
                                        "the probability density at lower is not equal to the probability density at upper" + stringify(std::fabs(alpha * z_plus - std::exp(z_plus) - alpha * z_minus + std::exp(z_minus))));
                }

                // calculate normalization factors that are independent of x
                norm = -1.0 * gsl_sf_lngamma(alpha) - std::log(std::fabs(lambda));
            }

            virtual ~LogGammaBlock()
            {
            }

            virtual std::string as_string() const
            {
                std::string result = "LogGamma: ";
                result += stringify(central) + " + " + stringify(sigma_upper) + " - " + stringify(sigma_lower);
                result += " (nu = " + stringify(nu) + ", lambda = " + stringify(lambda) + ", alpha = " + stringify(alpha) + ")";

                if (0 == _number_of_observations)
                    result += "; no observation";

                return result;
            }

            double cdf(const double & x) const
            {
                // transform exp of standardized coordinates
                double z = std::exp((x - nu) / lambda);

                if (lambda < 0)
                    return gsl_sf_gamma_inc_Q(alpha, z);

                else
                    return 1.0 - gsl_sf_gamma_inc_Q(alpha, z);
            }

            virtual double evaluate() const
            {
                double value = (cache[id] - nu) / lambda;

                return norm + alpha * value - std::exp(value);
            }

            virtual unsigned number_of_observations() const
            {
                return _number_of_observations;
            }

            // todo remove 3 sigma limits, or remove whole block altogether
            // draw from standard gamma, apply log, then shift and rescale
            virtual double sample(gsl_rng * rng) const
            {
                double x = 0.0;

                // allow difference of three standard observations in either direction
                double range_min = central - 3.0 * sigma_lower;
                double range_max = central + 3.0 * sigma_upper;

                while(true)
                {
                    x = lambda * std::log(gsl_ran_gamma(rng, alpha, 1.0)) + nu;

                    if (range_min < x && x < range_max)
                        break;
                }

                // now x is a pseudo measurement
                // pretend it were the mode of the pdf
                double nu_pseudo = x - lambda * std::log(alpha);

                // compare with central value, not the prediction!
                // seems strange, but we just want the distribution of the test statistic
                // and hopefully it is independent of the best fit parameters chosen
                double value = (central - nu_pseudo) / lambda;

                return norm + alpha * value - std::exp(value);
            }

            /*
             * To find the significance, it is necessary to determine the smallest interval
             * around the mode. This is achieved by finding the mirror point
             * on the other side of the mode which has the same probability density.
             * The solution is found numerically by root finding.
             */
            virtual double significance() const
            {
                double value = cache[id];

                // find point on the opposite side of the mode as starting value
                // if value right of the central value, the mirror has to be to the left
                double mirror =  2.0 * central - value;

                // assign the functions whose root we want to find
                gsl_function_fdf f;
                f.f = &LogGammaBlock::significance_function_f;
                f.df = &LogGammaBlock::significance_function_df;
                f.fdf = &LogGammaBlock::significance_function_fdf;
                f.params = (void *)this;

                int status;
                int iter = 0, max_iter = 400;

                const gsl_root_fdfsolver_type * solver_type = gsl_root_fdfsolver_steffenson;
                gsl_root_fdfsolver * solver = gsl_root_fdfsolver_alloc (solver_type);
                gsl_root_fdfsolver_set(solver, &f, mirror);

                // store last root estimate to check convergence
                double previous_mirror = mirror;
                do
                {
                    iter++;
                    status = gsl_root_fdfsolver_iterate (solver);
                    mirror = gsl_root_fdfsolver_root (solver);

                    status = gsl_root_test_delta(previous_mirror, mirror, 0.0, 1e-7);
                    previous_mirror = mirror;
                }
                while (status == GSL_CONTINUE && iter < max_iter);

                if (GSL_SUCCESS != status)
                {
                    Log::instance()->message("LogGammaBlock::significance", ll_error)
                        << "Could not find the mirror point, stopped after "
                        << iter << " iterations with f(" << mirror << ") = " << f.f(mirror, f.params);
                }

                gsl_root_fdfsolver_free(solver);

                // find probability
                const double p = std::fabs(cdf(value) - cdf(mirror));

                // transform to Gaussian sigmas
                const double abs_significance =  gsl_cdf_ugaussian_Pinv((p + 1) / 2.0);

                // determine sign: + if measured value (i.e. *mode()*) exceeds predicted value (i.e. *value*)
                return (central > value ? +1.0 : -1.0) * abs_significance;
            }

            /*
             *
             * For standardized coordinates z = (x - \nu) / \lambda,
             * f(z_{-}) = \alpha (z_{+} - z_{-}) - (e^{z_{+}} - e^{z_{-}}).
             * This is the log of the pdf, up to constants.
             * Need this as interface for GSL root solving
             */
            static double significance_function_f(double x, void * data)
            {
                LogGammaBlock * log_gamma = static_cast<LogGammaBlock *>(data);

                // the standardized points
                double zp = (log_gamma->cache[log_gamma->id] - log_gamma->nu) / log_gamma->lambda;
                double zm = (x - log_gamma->nu) / log_gamma->lambda;

                // the function value
                return log_gamma->alpha * (zp - zm) - std::exp(zp) + std::exp(zm);
            }

            static double significance_function_df(double x, void * data)
            {
                LogGammaBlock * log_gamma = static_cast<LogGammaBlock *>(data);

                // the standardized point
                double zm = (x - log_gamma->nu) / log_gamma->lambda;

                // the derivative
                return (std::exp(zm) - log_gamma->alpha) / log_gamma->lambda;
            }

            // Combine f and df in one function
            static void significance_function_fdf(double x, void * data, double * f, double * df)
            {
                LogGammaBlock * log_gamma = static_cast<LogGammaBlock *>(data);

                // the standardized points
                double zp = (log_gamma->cache[log_gamma->id] - log_gamma->nu) / log_gamma->lambda;
                double zm = (x - log_gamma->nu) / log_gamma->lambda;

                // the function value
                *f = log_gamma->alpha * (zp - zm) - std::exp(zp) + std::exp(zm);

                // the derivative
                *df = (std::exp(zm) - log_gamma->alpha) / log_gamma->lambda;
            }

            virtual TestStatistic primary_test_statistic() const
            {
                return test_statistics::Empty();
            }

            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                ObservablePtr observable = this->cache.observable(id)->clone(cache.parameters());

                return LogLikelihoodBlockPtr(new LogGammaBlock(cache, cache.add(observable),
                    central - sigma_lower, central, central + sigma_upper, alpha, lambda, _number_of_observations));
            }
        };

        struct AmorosoBlock :
            public LogLikelihoodBlock
        {
            ObservableCache cache;

            ObservableCache::Id id;

            const double physical_limit;

            const double theta, alpha, beta;

            double norm;

            unsigned _number_of_observations;

            AmorosoBlock(const ObservableCache & cache, ObservableCache::Id id, const double & physical_limit,
                         const double & theta, const double & alpha, const double & beta,
                         const unsigned & number_of_observations) :
                cache(cache),
                id(id),
                physical_limit(physical_limit),
                theta(theta),
                alpha(alpha),
                beta(beta),
                _number_of_observations(number_of_observations)
            {

                if (theta <= 0)
                    throw InternalError("LogLikelihoodBlock::Amoroso: scale parameter theta (" + stringify(theta) +
                            ") must be positive for an upper limit");
                if (alpha <= 0)
                    throw InternalError("LogLikelihoodBlock::Amoroso: shape parameter alpha (" + stringify(alpha) +
                            ") must be positive");
                if (beta <= 0)
                    throw InternalError("LogLikelihoodBlock::Amoroso: shape parameter beta (" + stringify(beta) +
                                ") must be positive");

                // calculate normalization factors that are independent of x
                norm = -1.0 * gsl_sf_lngamma(alpha) + std::log(std::fabs(beta / theta));
            }

            virtual ~AmorosoBlock()
            {
            }

            virtual std::string as_string() const
            {
                std::string name = cache.observable(id)->name().str();
                std::string result = "Amoroso limit: ";
                result += "mode at " + name + " = " + stringify(mode(), 5);
                result += " (a = " + stringify(physical_limit, 5) + ", theta = " + stringify(theta, 5);
                result += ", alpha = " + stringify(alpha, 5) + ", beta = " + stringify(beta, 5)+ ")";

                if (0 == _number_of_observations)
                    result += "; no observation";

                return result;
            }

            double cdf(const double & x) const
            {
                // Weibull transform
                double w = std::pow((x - physical_limit) / theta, beta);

                if (beta / theta < 0)
                    return gsl_sf_gamma_inc_Q(alpha, w);

                else
                    return 1.0 - gsl_sf_gamma_inc_Q(alpha, w);
            }

            virtual double evaluate() const
            {
                // standardized transform
                const double z = (cache[id] - physical_limit) / theta;

                return norm + (alpha * beta - 1) * std::log(z) - std::pow(z, beta);
            }

            inline double mode() const
            {
                return physical_limit + theta * std::pow(alpha - 1 / beta, 1 / beta);
            }

            virtual unsigned number_of_observations() const
            {
                return _number_of_observations;
            }

            /*
             * Draw from standard gamma.
             * Usually one would have to perform an inverse Weibull transform,
             * but when plugging it into the pdf, we would have to do a Weibull again
             * for the term in the exponential, so they cancel.
             * For the power term, need to remove the effect of beta.
             * The norm, where \alpha, \beta and \theta appear, is just the right one.
             */
            virtual double sample(gsl_rng * rng) const
            {
                const double w = gsl_ran_gamma(rng, alpha, 1.0);
                const double z = std::pow(w, 1 / beta);

                // compare with experimental distribution, not the prediction!
                // seems strange, but we just want the distribution of the test statistic
                // and hopefully it is independent of the best-fit parameters chosen
                return norm + (alpha * beta - 1) * std::log(z) - w;
            }

            virtual double significance() const
            {
                const double value = cache[id];

                // if mode at the boundary, the significance is just the cumulative at the point
                if (std::abs(alpha * beta - 1.0) < 1e-13)
                {
                    // find the probability between limit and current point
                    const double p = cdf(value);

                    // transform to standard Gaussian sigma units
                    return gsl_cdf_ugaussian_Pinv( (p + 1) / 2.0);
                }

                double x_min, x_max;
                if (value > mode())
                {
                    x_min = physical_limit;
                    x_max = mode();
                }
                else
                {
                    x_min = mode();
                    // increase upper boundary until it contains the point
                    x_max = x_min + (mode() - value);
                    while (significance_function_f(x_max, (void *)this) < 0)
                    {
                        x_max *= 2;
                    }
                }
                double estimate = (x_min + x_max) / 2;

                const gsl_root_fsolver_type * solver_type = gsl_root_fsolver_brent;
                gsl_root_fsolver * solver = gsl_root_fsolver_alloc (solver_type);
                gsl_function f;
                f.function = &significance_function_f;
                f.params = (void *)this;
                gsl_root_fsolver_set(solver, &f, x_min, x_max);

                int status;
                int iter = 0, max_iter = 400;

                do
                {
                    iter++;
                    status = gsl_root_fsolver_iterate(solver);
                    estimate = gsl_root_fsolver_root(solver);
                    x_min = gsl_root_fsolver_x_lower(solver);
                    x_max = gsl_root_fsolver_x_upper(solver);

                    status = gsl_root_test_interval(x_min, x_max, 0.0, 1e-7);
                }
                while (status == GSL_CONTINUE && iter < max_iter);

                if (GSL_SUCCESS != status)
                {
                    throw InternalError("Could not find the mirror point, stopped after "
                        + stringify(iter) + " iterations with f(" + stringify(estimate) + ") = " + stringify(f.function(estimate, f.params)));
                }

                gsl_root_fsolver_free(solver);

                // find probability of smaller excess ( 1-ordinary p value)
                const double p = std::fabs(cdf(value) - cdf(estimate));

                // transform to Gaussian sigmas (>= zero because p >= 0)
                const double abs_significance =  gsl_cdf_ugaussian_Pinv((p + 1) / 2.0);

                // determine sign: + if measured value (i.e. *mode()*) exceeds predicted value (i.e. *value*)
                return (mode() > value ? +1.0 : -1.0) * abs_significance;
            }

            /*
             *
             * For standardized coordinates z = (x - a) / \theta,
             * x_{-} = current estimate of root
             * x_{+} = fixed at current value of observable
             * f(z_{-}) = log(f(z_{+})) - log(f(z_{-}))
             * f(z_{-}) = \alpha (z_{+} - z_{-}) - (e^{z_{+}} - e^{z_{-}}).
             * This is the log of the pdf, up to constants.
             * Need this as interface for GSL root solving
             */
            static double significance_function_f(double x, void * data)
            {
                AmorosoBlock* a = static_cast<AmorosoBlock *>(data);

                // the standardized points
                const double zp = (a->cache[a->id] - a->physical_limit) / a->theta;
                const double zm = (x - a->physical_limit) / a->theta;

                // avoid infinity when zm is at physical_limit
                if (zm == 0.0)
                {
                    return std::numeric_limits<double>::max();
                }

                // the function value
                return (a->alpha * a->beta - 1) * (std::log(zp) - std::log(zm))
                    + std::pow(zm, a->beta) - std::pow(zp, a->beta);
            }

            virtual TestStatistic primary_test_statistic() const
            {
                return test_statistics::Empty();
            }

            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                ObservablePtr observable = this->cache.observable(id)->clone(cache.parameters());
                return LogLikelihoodBlockPtr(new AmorosoBlock(cache, cache.add(observable), physical_limit, theta, alpha, beta, _number_of_observations));
            }
        };

        struct MixtureBlock :
            public LogLikelihoodBlock
        {
            std::vector<LogLikelihoodBlockPtr> components;
            std::vector<double> weights;
            std::vector<std::array<double, 2>> test_stat;
            mutable std::vector<double> temp;

            MixtureBlock(const std::vector<LogLikelihoodBlockPtr> & components,
                         const std::vector<double> & weights,
                         const std::vector<std::array<double, 2>> & test_stat) :
                components(components),
                weights(weights),
                test_stat(test_stat),
                temp(weights.size(), 0.0)
            {
            }

            ~MixtureBlock()
            {
            }

            std::string as_string() const
            {
                std::string ret_val = "Mixture: \n";
                for (const auto & component : components)
                    ret_val += (*component).as_string() + '\n';
                return ret_val;
            }

            LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                std::vector<LogLikelihoodBlockPtr> clones;
                for (const auto & component : components)
                    clones.push_back((*component).clone(cache));

                return LogLikelihoodBlockPtr(new MixtureBlock(clones, weights, test_stat));
            }

            double evaluate() const
            {
                // find biggest element
                auto v = temp.begin();

                for (auto c = components.cbegin() ; c != components.cend() ; ++c, ++v)
                    *v = (**c).evaluate();

                auto max_val = std::max_element(temp.cbegin(), temp.cend());
                double ret_val = 0;

                // computed weighted sum, renormalize exponents
                v = temp.begin();
                for (auto w = weights.cbegin(); w != weights.cend() ; ++w, ++v)
                {
                    ret_val += *w * std::exp(*v - *max_val);
                }

                ret_val = std::log(ret_val) + *max_val;

                return ret_val;
            }

            unsigned number_of_observations() const
            {
                return components.front()->number_of_observations();
            }

            double sample(gsl_rng * /*rng*/) const
            {
                throw InternalError("LogLikelihoodBlock::MixtureBlock::sample() not implemented yet");
            }

            double significance() const
            {
                double value = -2.0 * evaluate();
                for (const auto & pair : test_stat)
                {
                    if (value <= pair[1]) return pair[0];
                }

                return std::numeric_limits<double>::quiet_NaN();
            }

            virtual TestStatistic primary_test_statistic() const
            {
                return test_statistics::ChiSquare(
                    gsl_cdf_chisq_Pinv(gsl_cdf_chisq_P(power_of<2>(significance()), 1), number_of_observations()),
                    number_of_observations());
            }
        };

        struct MultivariateGaussianBlock :
            public LogLikelihoodBlock
        {
            ObservableCache _cache;

            std::vector<ObservableCache::Id> _ids;

            const unsigned _dim_pred;
            const unsigned _dim_meas;

            // inputs
            gsl_vector * const _mean;
            gsl_matrix * const _covariance;
            gsl_matrix * const _response;
            const unsigned _number_of_observations;

            // the normalization constant of the density
            const double _norm;

            // cholesky matrix of covariance, and inverse of covariance
            gsl_matrix * _chol;
            gsl_matrix * _covariance_inv;

            // temporary storage for evaluation
            gsl_vector * _observables;
            gsl_vector * _measurements;
            gsl_vector * _measurements_2;

            MultivariateGaussianBlock(const ObservableCache & cache, const std::vector<ObservableCache::Id> && ids,
                    gsl_vector * mean, gsl_matrix * covariance, gsl_matrix * response, const unsigned & number_of_observations) :
                _cache(cache),
                _ids(ids),
                _dim_pred(ids.size()),
                _dim_meas(mean->size),
                _mean(mean),
                _covariance(covariance),
                _response(response),
                _number_of_observations(number_of_observations),
                _norm(compute_norm()),
                _chol(gsl_matrix_alloc(covariance->size1, covariance->size2)),
                _covariance_inv(gsl_matrix_alloc(covariance->size1, covariance->size2)),
                _observables(gsl_vector_alloc(_dim_pred)),
                _measurements(gsl_vector_alloc(_dim_meas)),
                _measurements_2(gsl_vector_alloc(_dim_meas))
            {
                if (_covariance->size1 != _covariance->size2)
                    throw InternalError("MultivariateGaussianBlock: covariance matrix is not a square matrix");

                if (_dim_meas != _covariance->size1)
                    throw InternalError("MultivariateGaussianBlock: number of measurements and dimension of covariance matrix are not identical");

                if (_dim_meas != _response->size1)
                    throw InternalError("MultivariateGaussianBlock: number of measurements and number of rows in response matrix are not identical");

                if (_dim_pred != _response->size2)
                    throw InternalError("MultivariateGaussianBlock: number of predictions and number of columns in response matrix are not identical");

                // cholesky decomposition (informally: the sqrt of the covariance matrix)
                // the GSL matrix contains both the cholesky and its transpose, see GSL reference, ch. 14.5
                cholesky();
                invert_covariance();

                // keep only the lower and diagonal parts, set upper parts to zero
                for (unsigned i = 0; i < _dim_meas ; ++i)
                {
                    for (unsigned j = i + 1 ; j < _dim_meas ; ++j)
                    {
                        gsl_matrix_set(_chol, i, j, 0.0);
                    }
                }
            }

            virtual ~MultivariateGaussianBlock()
            {
                gsl_matrix_free(_covariance_inv);
                gsl_matrix_free(_chol);
                gsl_matrix_free(_covariance);
                gsl_matrix_free(_response);

                gsl_vector_free(_measurements_2);
                gsl_vector_free(_measurements);
                gsl_vector_free(_observables);
                gsl_vector_free(_mean);
            }

            virtual std::string as_string() const
            {
                const auto k = _mean->size;

                std::string result = "Multivariate Gaussian: ";
                result += "means = ( ";
                for (std::size_t i = 0 ; i < k ; ++i)
                {
                    result += stringify(gsl_vector_get(_mean, i)) + " ";
                }
                result += "), covariance matrix = (";
                for (std::size_t i = 0 ; i < k ; ++i)
                {
                    result += "( ";
                    for (std::size_t j = 0 ; j < k ; ++j)
                    {
                        result += stringify(gsl_matrix_get(_covariance, i, j)) + " ";
                    }
                    result += ")";
                }
                result += "), inverse covariance matrix = (";
                for (std::size_t i = 0 ; i < k ; ++i)
                {
                    result += "( ";
                    for (std::size_t j = 0 ; j < k ; ++j)
                    {
                        result += stringify(gsl_matrix_get(_covariance_inv, i, j)) + " ";
                    }
                    result += ")";
                }
                result += " )";

                if (0 == _number_of_observations)
                    result += "; no observation";

                return result;
            }

            // compute cholesky decomposition of covariance matrix
            void cholesky()
            {
                // copy covariance matrix
                gsl_matrix_memcpy(_chol, _covariance);
                if (GSL_SUCCESS != GSL_LINALG_CHOLESKY_DECOMP(_chol))
                {
                    throw InternalError("MultivariateGaussianBlock: Cholesky decomposition failed");
                }
            }

            // invert covariance matrix based on previously obtained Cholesky decomposition
            void invert_covariance()
            {
                // copy cholesky matrix
                gsl_matrix_memcpy(_covariance_inv, _chol);

                // compute inverse matrix from cholesky
                if (GSL_SUCCESS != gsl_linalg_cholesky_invert(_covariance_inv))
                {
                    throw InternalError("MultivariateGaussianBlock: Cholesky inversion failed");
                }
            }


            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                const auto dim_meas = _mean->size;
                const auto dim_pred = _ids.size();

                std::vector<ObservableCache::Id> ids;

                // add observables to cache
                for (auto i = 0u ; i < dim_pred ; ++i)
                {
                    ids.push_back(cache.add(this->_cache.observable(this->_ids[i])->clone(cache.parameters())));
                }

                gsl_vector * mean = gsl_vector_alloc(dim_meas);
                gsl_vector_memcpy(mean, _mean);

                gsl_matrix * covariance = gsl_matrix_alloc(dim_meas, dim_meas);
                gsl_matrix_memcpy(covariance, _covariance);

                gsl_matrix * response = gsl_matrix_alloc(dim_meas, dim_pred);
                gsl_matrix_memcpy(response, _response);

                return LogLikelihoodBlockPtr(new MultivariateGaussianBlock(cache, std::move(ids), mean, covariance, response, _number_of_observations));
            }

            // compute the normalization constant on log scale
            // -k/2 * log 2 Pi - 1/2 log(abs(det(V^{-1})))
            double compute_norm()
            {
                // copy covariance matrix
                gsl_matrix * M = gsl_matrix_alloc(_dim_meas, _dim_meas);
                gsl_matrix_memcpy(M, _covariance);

                // find LU decomposition
                int signum = 0;
                gsl_permutation * p = gsl_permutation_alloc(_dim_meas);
                gsl_linalg_LU_decomp(M, p, &signum);

                // calculate determinant
                const double log_det = gsl_linalg_LU_lndet(M);

                gsl_permutation_free(p);
                gsl_matrix_free(M);

                return -0.5 * _dim_meas * std::log(2 * M_PI) - 0.5 * log_det;
            }

            double chi_square() const
            {
                // read observable values from cache, and subtract mean
                for (auto i = 0u ; i < _dim_pred ; ++i)
                {
                    gsl_vector_set(_observables, i, _cache[_ids[i]]);
                }

                // prepare for centering
                //   measurements <- mean
                gsl_vector_memcpy(_measurements, _mean);

                // apply response matrix and center the gaussian:
                //   measurements <- R * observables - measurements
                gsl_blas_dgemv(CblasNoTrans, 1.0, _response, _observables, -1.0, _measurements);

                // observables <- inv(covariance) * measurements
                gsl_blas_dgemv(CblasNoTrans, 1.0, _covariance_inv, _measurements, 0.0, _measurements_2);

                double result;
                gsl_blas_ddot(_measurements, _measurements_2, &result);

                return result;
            }

            virtual double evaluate() const
            {
                return _norm - 0.5 * chi_square();
            }

            virtual unsigned number_of_observations() const
            {
                return _number_of_observations;
            }

            virtual double sample(gsl_rng * rng) const
            {
                // generate standard normals in observables
                for (auto i = 0u ; i < _dim_meas ; ++i)
                {
                    gsl_vector_set(_measurements, i, gsl_ran_ugaussian(rng));
                }

                // transform: observables2 <- _chol * observables
                gsl_blas_dgemv(CblasNoTrans, 1.0, _chol, _measurements, 0.0, _measurements_2);

                // To be consistent with the univariate Gaussian, we would center observables around theory,
                // then compare to theory. Hence we can forget about theory, and stay centered on zero.
                // transform: observables <- inv(covariance) * observables2
                gsl_blas_dgemv(CblasNoTrans, 1.0, _covariance_inv, _measurements_2, 0.0, _measurements);

                double result;
                gsl_blas_ddot(_measurements, _measurements_2, &result);
                result *= -0.5;
                result += _norm;

                return result;
            }

            virtual double significance() const
            {
                const auto chi_squared = this->chi_square();

                // find probability of this excess or less ( 1 - usual p-value)
                const double p = gsl_cdf_chisq_P(chi_squared, _mean->size);

                // transform to standard Gaussian sigma units
                // return  significance is >= 0, since p >= 0
                // and a negative significance is ruled out by definition
                return gsl_cdf_ugaussian_Pinv((p + 1) / 2.0);
            }

            virtual TestStatistic primary_test_statistic() const
            {
                return test_statistics::ChiSquare(this->chi_square(), this->_number_of_observations);
            }
        };

        struct UniformBoundBlock :
            public LogLikelihoodBlock
        {
            ObservableCache cache;

            std::vector<ObservableCache::Id> ids;

            const unsigned number_of_observables;

            const double bound;
            const double uncertainty;

            UniformBoundBlock(const ObservableCache & cache, const std::vector<ObservableCache::Id> && ids,
                              const double & bound, const double & uncertainty) :
                cache(cache),
                ids(ids),
                number_of_observables(ids.size()),
                bound(bound),
                uncertainty(uncertainty)
            {
            }

            virtual ~UniformBoundBlock()
            {
            }

            virtual std::string as_string() const
            {
                std::string result = "UniformBound: ";
                result += "bound = " + stringify(bound) + " +- " + stringify(uncertainty);

                return result;
            }

            virtual double evaluate() const
            {
                double saturation = 0.0;

                for (auto i : ids)
                {
                    saturation += cache[i];
                }

                if (saturation < 0.0)
                {
                    throw InternalError("Contribution to the uniform bound must be positive; found to be negative!");
                }
                else if ((0.0 <= saturation) && (saturation < bound))
                {
                    return 0.0;
                }
                else
                {
                    if (uncertainty == 0.0)
                    {
                        return - std::numeric_limits<double>::infinity();
                    }
                    else
                    {
                        // add a gaussian like penalty
                        return - 0.5 * power_of<2>((saturation - bound) / uncertainty);
                    }
                }
            }

            virtual unsigned number_of_observations() const
            {
                return 0.0;
            }

            virtual double sample(gsl_rng * /*rng*/) const
            {
                return 0.0;
            }

            virtual double significance() const
            {
                double saturation = 0.0;

                for (auto i : ids)
                {
                    saturation += cache[i];
                }

                if (saturation < 0.0)
                {
                    throw InternalError("Contribution to the uniform bound must be positive; found to be negative!");
                }
                else if ((0.0 <= saturation) && (saturation < bound))
                {
                    return 0.0;
                }
                else
                {
                    if (uncertainty == 0.0)
                    {
                        return - std::numeric_limits<double>::infinity();
                    }
                    else
                    {
                        // add a gaussian like penalty
                        return (saturation - bound) / uncertainty;
                    }
                }
            }

            virtual TestStatistic primary_test_statistic() const
            {
                double chi = significance();
                return test_statistics::ChiSquare(chi * chi, 1, chi);
            }

            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                std::vector<ObservableCache::Id> ids;

                // add observables to cache
                for (auto i = 0u ; i < number_of_observables ; ++i)
                {
                    ids.push_back(cache.add(this->cache.observable(this->ids[i])->clone(cache.parameters())));
                }

                return LogLikelihoodBlockPtr(new UniformBoundBlock(cache, std::move(ids), bound, uncertainty));
            }
        };
    }

    LogLikelihoodBlock::~LogLikelihoodBlock()
    {
    }

    LogLikelihoodBlockPtr
    LogLikelihoodBlock::Gaussian(ObservableCache cache, const ObservablePtr & observable,
            const double & min, const double & central, const double & max,
            const unsigned & number_of_observations)
    {
        // check input
        if (min >= central)
            throw InternalError("LogLikelihoodBlock::Gaussian: min value >= central value");

        if (max <= central)
            throw InternalError("LogLikelihoodBlock::Gaussian: max value <= central value");

        unsigned index = cache.add(observable);

        return LogLikelihoodBlockPtr(new implementation::GaussianBlock(cache, index, min, central, max, number_of_observations));
    }

    LogLikelihoodBlockPtr
    LogLikelihoodBlock::LogGamma(ObservableCache cache, const ObservablePtr & observable,
            const double & min, const double & central, const double & max,
            const double & alpha, const double & lambda,
            const unsigned & number_of_observations)
    {

        // check input
        if (min >= central)
            throw InternalError("LogLikelihoodBlock::LogGamma: min value >= central value");

        if (max <= central)
            throw InternalError("LogLikelihoodBlock::LogGamma: max value <= central value");

        if (alpha <= 0)
            throw InternalError("LogLikelihoodBlock::LogGamma: shape parameter alpha (" + stringify(alpha) +
                                ") must be positive");

        unsigned index = cache.add(observable);

        return LogLikelihoodBlockPtr(new implementation::LogGammaBlock(cache, index, min, central, max, alpha, lambda, number_of_observations));
    }

    LogLikelihoodBlockPtr
    LogLikelihoodBlock::Amoroso(ObservableCache cache, const ObservablePtr & observable,
                                const double & physical_limit, const double & upper_limit_10,
                                const double & upper_limit_50, const double & upper_limit_90,
                                const double & theta, const double & alpha, const double & beta,
                                const unsigned & number_of_observations)
    {
        // check input
        if (upper_limit_10 <= physical_limit)
            throw InternalError("LogLikelihoodBlock::AmorosoLimit: upper_limit_10 <= physical_limit");

        if (upper_limit_50 <= physical_limit)
            throw InternalError("LogLikelihoodBlock::AmorosoLimit: upper_limit_50 <= physical_limit");

        if (upper_limit_90 <= upper_limit_50)
            throw InternalError("LogLikelihoodBlock::AmorosoLimit: upper_limit_90 <= upper_limit_50");

        unsigned index = cache.add(observable);
        implementation::AmorosoBlock * a = new implementation::AmorosoBlock(cache, index, physical_limit, theta, alpha, beta, number_of_observations);

        // check consistency
        if (std::abs(a->cdf(upper_limit_10) - 0.10) > 1e-4)
        {
            throw InternalError("LogLikelihood::Amoroso.ctor: For the current parameter values, cdf(x_10) = "
                + stringify(a->cdf(upper_limit_10)) + " deviates from 10%.");
        }
        if (std::abs(a->cdf(upper_limit_50) - 0.50) > 1e-4)
        {
            throw InternalError("LogLikelihood::Amoroso.ctor: For the current parameter values, cdf(x_50) = "
                + stringify(a->cdf(upper_limit_50)) + " deviates from 50%.");
        }
        if (std::abs(a->cdf(upper_limit_90) - 0.90) > 1e-4)
        {
            throw InternalError("LogLikelihood::Amoroso.ctor: For the current parameter values, cdf(x_90) = "
                + stringify(a->cdf(upper_limit_90)) + " deviates from 90%.");
        }

        return LogLikelihoodBlockPtr(a);
    }

    LogLikelihoodBlockPtr
    LogLikelihoodBlock::Amoroso(ObservableCache cache, const ObservablePtr & observable,
        const double & physical_limit, const double & theta, const double & alpha, const double & beta,
        const unsigned & number_of_observations)
    {
        unsigned index = cache.add(observable);
        return LogLikelihoodBlockPtr(new implementation::AmorosoBlock(cache, index, physical_limit, theta, alpha, beta, number_of_observations));
    }

    LogLikelihoodBlockPtr
    LogLikelihoodBlock::Mixture(const std::vector<LogLikelihoodBlockPtr> & components,
                                const std::vector<double> & weights,
                                const std::vector<std::array<double, 2>> & test_stat)
    {
        if (components.size() != weights.size())
            throw InternalError("LogLikelihoodBlock::Mixture(): components and weights don't match");

        // normalize weights
        double sum = 0; // std::accumulate(weights.cbegin(), weights.cend(), 0);
        for (double weight : weights)
            sum += weight;

        std::vector<double> norm_weights(weights);

        std::transform(weights.cbegin(), weights.cend(), norm_weights.begin(),
                       std::bind(std::multiplies<double>(), 1.0 / sum, std::placeholders::_1));
        Log::instance()->message("MixtureBlock()", ll_debug)
            << "sum = " << sum <<  ", norm. weights " << stringify_container(norm_weights);

        return LogLikelihoodBlockPtr(new implementation::MixtureBlock(components, norm_weights, test_stat));
    }

    LogLikelihoodBlockPtr
    LogLikelihoodBlock::MultivariateGaussian(ObservableCache cache, const std::vector<ObservablePtr> & observables,
            gsl_vector * mean, gsl_matrix * covariance, gsl_matrix * response, const unsigned & number_of_observations)
    {
        std::vector<unsigned> indices;
        for (auto & o : observables)
        {
            indices.push_back(cache.add(o));
        }

        return LogLikelihoodBlockPtr(new implementation::MultivariateGaussianBlock(cache, std::move(indices), mean, covariance, response, number_of_observations));
    }

    LogLikelihoodBlockPtr
    LogLikelihoodBlock::UniformBound(ObservableCache cache, const std::vector<ObservablePtr> & observables,
            const double & bound, const double & uncertainty)
    {
        std::vector<unsigned> indices;
        for (auto & o : observables)
        {
            indices.push_back(cache.add(o));
        }

        return LogLikelihoodBlockPtr(new implementation::UniformBoundBlock(cache, std::move(indices), bound, uncertainty));
    }

    template <>
    struct WrappedForwardIteratorTraits<LogLikelihood::ConstraintIteratorTag>
    {
        using UnderlyingIterator = std::vector<Constraint>::iterator;
    };
    template class WrappedForwardIterator<LogLikelihood::ConstraintIteratorTag, Constraint>;

    template <>
    struct Implementation<LogLikelihood>
    {
        Parameters parameters;

        // Cache observable predictions
        ObservableCache cache;

        // Container for all named constraints
        std::vector<Constraint> constraints;

        // Container for all external likelihood blocks
        std::vector<LogLikelihoodBlockPtr> external_blocks;

        Implementation(const Parameters & parameters) :
            parameters(parameters),
            cache(parameters)
        {
        }

        std::pair<double, double>
        bootstrap_p_value(const unsigned & datasets)
        {
            // Algorithm:
            // 1. For fixed parameters, create data sets under the model.
            // 2. Use the likelihood as test statistic, T=L, calculate it for each data set.
            // 3. Compare with likelihood of "observed" data set to define p-value
            //      p = #llh < llh(obs) / #trials

            // observed value
            double t_obs = 0;

            // set up for sampling
            for (const auto & constraint : constraints)
            {
                for (auto b = constraint.begin_blocks(), b_end = constraint.end_blocks() ; b != b_end ; ++b)
                {
                    if (! (*b)->number_of_observations())
                        continue;
                    t_obs += (*b)->evaluate();
                }
            }

            Log::instance()->message("log_likelihood.bootstrap_pvalue", ll_informational)
                                     << "The value of the test statistic (total likelihood) "
                                     << "for the current parameters is = " << t_obs;

            // count data sets with smaller likelihood
            unsigned n_low = 0;

            // test value
            double t;

            gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, datasets);

            Log::instance()->message("log_likelihood.bootstrap_pvalue", ll_informational)
                                     << "Begin sampling " << datasets << " simulated "
                                     << "values of the likelihood";

            // collect samples
            for (unsigned i = 0 ; i < datasets ; ++i)
            {
                t = 0.0;
                for (const auto & constraint : constraints)
                {
                    for (auto b = constraint.begin_blocks(), b_end = constraint.end_blocks() ; b != b_end ; ++b)
                    {
                        t += (*b)->sample(rng);
                    }
                }

                if (t < t_obs)
                {
                    ++n_low;
                }
            }

            // mode of binomial posterior
            double p = n_low / double(datasets);

            // determine uncertainty of p-value
            // Just the variance of a binomial posterior
            double p_expected = double(n_low + 1) / double(datasets + 2);
            double uncertainty = std::sqrt(p_expected * (1 - p_expected) / double(datasets + 3));

            Log::instance()->message("log_likelihood.bootstrap_pvalue", ll_informational)
                                     << "The simulated p-value is " << p
                                     << " with uncertainty " << uncertainty;

            gsl_rng_free(rng);

            return std::make_pair(p, uncertainty);
        }

        double log_likelihood() const
        {
            double result = 0.0;

            // loop over all constraint-based likelihood blocks
            for (const auto & constraint : constraints)
            {
                for (auto b = constraint.begin_blocks(), b_end = constraint.end_blocks() ; b != b_end ; ++b)
                {
                    double llh = (*b)->evaluate();
                    if (! std::isfinite(llh))
                        return -std::numeric_limits<double>::infinity();

                    result += llh;
                }
            }

            // loop over all external likelihood blocks
            for (const auto & block : external_blocks)
            {
                double llh = block->evaluate();
                if (! std::isfinite(llh))
                    return -std::numeric_limits<double>::infinity();

                result += llh;
            }

            return result;
        }
    };

    LogLikelihood::LogLikelihood(const Parameters & parameters) :
        PrivateImplementationPattern<LogLikelihood>(new Implementation<LogLikelihood>(parameters))
    {
    }

    LogLikelihood::~LogLikelihood()
    {
    }

    void
    LogLikelihood::add(const ObservablePtr & observable, const double & min, const double & central, const double & max,
            const unsigned & number_of_observations)
    {
        LogLikelihoodBlockPtr b = LogLikelihoodBlock::Gaussian(_imp->cache, observable, min, central, max, number_of_observations);
        _imp->constraints.push_back(
            Constraint(observable->name(), std::vector<ObservablePtr>{ observable }, std::vector<LogLikelihoodBlockPtr>{ b }));
    }

    void
    LogLikelihood::add(const Constraint & constraint)
    {
        std::vector<ObservablePtr> observables;
        std::vector<LogLikelihoodBlockPtr> blocks;

        for (auto b = constraint.begin_blocks(), b_end = constraint.end_blocks() ; b != b_end ; ++b)
        {
            // Clone each LogLikelihoodBlock onto our ObservableCache
            blocks.push_back((*b)->clone(_imp->cache));
        }

        std::copy(constraint.begin_observables(), constraint.end_observables(), std::back_inserter(observables));

        // retain a proper copy of the constraint to iterate over
        _imp->constraints.push_back(Constraint(constraint.name(), observables, blocks));
    }

    void
    LogLikelihood::add(const LogLikelihoodBlockPtr & block)
    {
        _imp->external_blocks.push_back(block->clone(_imp->cache));
    }

    LogLikelihood::ConstraintIterator
    LogLikelihood::begin() const
    {
        return ConstraintIterator(_imp->constraints.begin());
    }

    LogLikelihood::ConstraintIterator
    LogLikelihood::end() const
    {
        return ConstraintIterator(_imp->constraints.end());
    }

    std::pair<double, double>
    LogLikelihood::bootstrap_p_value(const unsigned & datasets)
    {
        return _imp->bootstrap_p_value(datasets);
    }

    LogLikelihood
    LogLikelihood::clone() const
    {
        LogLikelihood result(_imp->parameters.clone());
        result._imp->cache = _imp->cache.clone(result._imp->parameters);

        for (const auto & constraint : _imp->constraints)
        {
            result.add(constraint);
        }

        for (const auto & block : _imp->external_blocks)
        {
            result.add(block);
        }

        return result;
    }

    unsigned
    LogLikelihood::number_of_observations() const
    {
        unsigned result = 0;
        for (const auto & constraint : _imp->constraints)
        {
            for (auto b = constraint.begin_blocks(), b_end = constraint.end_blocks() ; b != b_end ; ++b)
            {
                result += (**b).number_of_observations();
            }
        }

        return result;
    }

    Parameters
    LogLikelihood::parameters() const
    {
        return _imp->parameters;
    }

    ObservableCache
    LogLikelihood::observable_cache() const
    {
        return _imp->cache;
    }

    double
    LogLikelihood::operator() () const
    {
        _imp->cache.update();

        return _imp->log_likelihood();
    }
}
