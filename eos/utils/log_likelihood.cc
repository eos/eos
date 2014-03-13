/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2014 Danny van Dyk
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

#include <eos/utils/log_likelihood.hh>
#include <eos/utils/equation_solver.hh>
#include <eos/utils/log.hh>
#include <eos/utils/observable_cache.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/verify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>

namespace eos
{
    namespace implementation
    {
        struct GaussianBlock :
            public LogLikelihoodBlock
        {
            ObservableCache cache;

            ObservableCache::Id id;

            double central, sigma_lower, sigma_upper;

            // the probability covered to the left of the central value
            double prob_lower;

            // coefficients needed for sampling from asymmetric Gaussian on finite support
            // the cumulative is a piecewise function
            // CDF(x) = CDF_lower(x, sigma_lower) if x < central, else CDF_upper(x, sigma_upper)
            // To ensure that cumulative is
            // a) continuous at the central value
            // b) zero when x < x_min
            // c) one when  x > x_max
            // need to fix the coefficients in
            // @f$CDF_{lower}(x) = c_{lower} * \right( \Phi(x, \sigma_{lower}) - \Phi(x_{min},\sigma_{lower}) \left)$@f
            // @f$CDF_{upper}(x) = c_{upper} * \right( \Phi(x, \sigma_{upper}) + P_{lower}/c_{upper} - 1/2\left)$@f.
            double c_lower, c_upper;
            double phi_min, phi_max;

            unsigned _number_of_observations;

            GaussianBlock(const ObservableCache & cache, ObservableCache::Id id, const double & min, const double & central, const double & max,
                    const unsigned & number_of_observations) :
                cache(cache),
                id(id),
                central(central),
                sigma_lower(central - min),
                sigma_upper(max - central),
                _number_of_observations(number_of_observations)
            {
            }

            virtual ~GaussianBlock()
            {
            }

            virtual std::string as_string() const
            {
                std::string result = "Gaussian: ";
                result += stringify(central);
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
                if (value > central)
                    sigma = sigma_upper;
                else
                    sigma = sigma_lower;

                double chi = (value - central) / sigma;
                double norm = 1.0 / std::sqrt(2.0 * M_PI) / sigma;

                return std::log(norm) - chi * chi / 2.0;
            }

            virtual unsigned number_of_observations() const
            {
                return _number_of_observations;
            }

            virtual void prepare_sampling()
            {
                // fix the model, thus prediction of central value
                // since we can't explicitly change the data, we make the bootstrap
                // assumption: the measured is the true distribution. This creates
                // a bias towards higher p-values which cannot be easily estimated.
                // Would rather use cache[id], but we don't know what the experimental distribution
                // would look like, so use central. In this case, no difference.
                double mu = central;

                // allow a range of observables of three sigmas around model prediction
                // in order to avoid possible unphysical values
                double range_min = mu - 3.0 * sigma_lower;
                double range_max = mu + 3.0 * sigma_upper;

                // following constants needed for proper sampling
                phi_min = gsl_cdf_gaussian_P(range_min - mu, sigma_lower);
                phi_max = gsl_cdf_gaussian_P(range_max - mu, sigma_upper);
                prob_lower = (0.5 - phi_min) / (phi_max - phi_min);
                c_lower = 2.0 * prob_lower / (1.0 - phi_min);
                c_upper = (1.0 - prob_lower) / (phi_max - 0.5);
            }

            virtual double sample(gsl_rng * rng) const
            {
                // find out if sample in upper or lower part
                double u = gsl_rng_uniform(rng);

                // get a sample  observable using inverse transform method
                double obs, sigma;
                if (u < prob_lower)
                {
                    obs = gsl_cdf_gaussian_Pinv(u / c_lower + phi_min, sigma_lower) + central;
                    sigma = sigma_lower;
                }
                else
                {
                    obs = gsl_cdf_gaussian_Pinv((u - prob_lower) / c_upper + 0.5, sigma_upper) + central;
                    sigma = sigma_upper;
                }

                // calculate the properly normalized log likelihood
                // note that we generate from gaussian around cache[id],
                // but compare with central
                return -std::log(std::sqrt(2.0 * M_PI) * sigma) - power_of<2>((central - obs) / sigma) / 2.0;
            }

            virtual double significance() const
            {
                double value = cache[id];
                double sigma = 0.0;

                // allow for asymmetric Gaussian uncertainty
                if (value > central)
                    sigma = sigma_upper;
                else
                    sigma = sigma_lower;

                // return positive significance if measured value exceeds predictions
                return (central - value) / sigma;
            }

            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                ObservablePtr observable = this->cache.observable(id)->clone(cache.parameters());

                return LogLikelihoodBlockPtr(new GaussianBlock(cache, cache.add(observable), central - sigma_lower, central, central + sigma_upper, _number_of_observations));
            }
        };

        // For more details on the LogGamma distribution, see [C2004].
        struct LogGammaBlock :
            public LogLikelihoodBlock
        {
            ObservableCache cache;

            ObservableCache::Id id;

            double central, sigma_lower, sigma_upper;

            double nu, lambda, alpha;

            double norm;

            unsigned _number_of_observations;

            LogGammaBlock(const ObservableCache & cache, ObservableCache::Id id, const double & min, const double & central, const double & max,
                    const unsigned & number_of_observations) :
                cache(cache),
                id(id),
                central(central),
                sigma_lower(central - min),
                sigma_upper(max - central),
                _number_of_observations(number_of_observations)
            {
                // standardize scales, such that lower is one, thus fixing the sign of lambda
                const double sigma_plus( (sigma_upper > sigma_lower)? sigma_upper / sigma_lower : sigma_lower / sigma_upper);
                const double sigma_minus(1.0);
                if (sigma_plus < 1 + 6e-2)
                {
                    Log::instance()->message("LogLikelihoodBlock::LogGamma.ctor", ll_warning)
                        << "For nearly symmetric uncertainties (" << sigma_lower << " vs " << sigma_upper
                        << "), this procedure may fail to find the correct parameter values. "
                           "Please use a Gaussian block instead.";
                }

                // for positive skew, \lambda is negative
                // in the fit, \lambda always considered negative, so it only changes sign for negative skew
                const double lambda_scale_factor = sigma_upper > sigma_lower ? sigma_lower / sigma_minus : -1.0 * sigma_upper / sigma_minus;

                /* find the parameters using good starting values. Assume upper > lower=1, and fix sign at the end */

                // functions only good up to 10 % in region where uncertainties differ by 3-100 %. Not accurate for very [a]symmetric cases.
                // They are found by playing around in mathematica, using that
                // 1. \alpha depends only on \sigma_{+}
                // 2. \lambda is a scale parameter, so we solve the problem for the standard case, and then rescale
                double lambda_initial = -56 + 55.0 * gsl_cdf_gaussian_P(sigma_plus - 1.0, 0.05);
                double alpha_initial = std::pow(1.13 / (sigma_plus - 1), 1.3);

                EquationSolver solver(EquationSolver::Config::Default());

                //add free Parameter: initial value, error
                solver.add("lambda", lambda_initial, lambda_initial / 10.0);
                //add positive Parameter: initial value, error, min, max. \alpha for 5% asymmetry at 500, so 1000 is well above
                solver.add("alpha", alpha_initial, alpha_initial / 5.0, 0.0, 1000);

                // add constraints for standardized problem
                solver.add(std::bind(&LogGammaBlock::constraint, *this, std::placeholders::_1, sigma_plus, sigma_minus));

                // check errors manually, then restore default behavior later
                gsl_error_handler_t * default_gsl_error_handler = gsl_set_error_handler_off();

                // find the solution
                auto solution = solver.solve();

                // global minimum at zero value, often Minuit claims to not have found it while it actually did
                if (! solution.valid && solution.value > 1e-3)
                {
                    Log::instance()->message("LogLikelihood::LogGamma.ctor", ll_error)
                        << "Solution of constraints failed";
                }

                // now we have all values
                lambda = lambda_scale_factor * solution.parameters[0];
                alpha  = solution.parameters[1];
                nu = central - lambda * std::log(alpha);

                // calculate normalization factors that are independent of x
                norm = -1.0 * gsl_sf_lngamma(alpha) - std::log(std::fabs(lambda));

                // restore default error handler
                gsl_set_error_handler(default_gsl_error_handler);
            }

            LogGammaBlock(const ObservableCache & cache, ObservableCache::Id id,
                          const double & min, const double & central, const double & max,
                          const double & lambda, const double & alpha,
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
                static const double eps = 1e-4;
                if (std::abs(cdf(central + sigma_upper) - cdf(central - sigma_lower) - 0.68268949213708585 ) > eps)
                {
                    throw InternalError("LogLikelihoodBlock::LogGamma.ctor: For the current parameter values, "
                                        "the interval [lower, upper] doesn't contain approx. 68%");

                }
                const double z_plus = (central + sigma_upper - nu) / lambda;
                const double z_minus = (central - sigma_lower - nu) / lambda;
                if (std::fabs(alpha * z_plus - std::exp(z_plus) - alpha * z_minus + std::exp(z_minus)) > eps)
                {
                    throw InternalError("LogLikelihoodBlock::LogGamma.ctor: For the current parameter values, "
                                        "the probability density at lower is not equal to the probability density at upper");
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

            // Implements the constraint that the cumulative(x = \mu +- sigma_{+-} | \lambda, \alpha) = 0.84 [0.16]
            double constraint(const std::vector<double> & parameters, const double & sigma_plus, const double & sigma_minus)
            {
                if (parameters.size() != 2)
                {
                    throw InternalError("LogLikelihoodBlock::LogGamma.constraint: parameter dimension is "
                        + stringify(parameters.size()) + ", should be 2.");
                }

                const double & lambda = parameters[0];
                const double & alpha  = parameters[1];

                // standardized mode at 0
                const double nu = 0.0 - lambda * std::log(alpha);

                gsl_sf_result result;

                // standardized coordinates at plus/minus
                const double z_plus = (sigma_plus - nu) / lambda;
                const double z_minus = (-sigma_minus  - nu) / lambda;

                // first constraint: pdf's should be equal, neglect prefactors
                const double first = std::fabs(alpha * z_plus - std::exp(z_plus) - alpha * z_minus + std::exp(z_minus));

                // second constraint: 68% interval
                int ret_code = gsl_sf_gamma_inc_Q_e(alpha, std::exp(z_plus), &result);
                if ( ! GSL_SUCCESS == ret_code)
                    throw InternalError("LogLikelihoodBlock::LogGamma: cannot evaluate cumulative at lambda = " + stringify(lambda)
                                        + ", alpha = " + stringify(alpha) + ". GSL reports: " + gsl_strerror(ret_code)
                                        + ". Perhaps the input is too [a]symmetric?");
                const double cdf_plus = result.val;

                ret_code = gsl_sf_gamma_inc_Q_e(alpha, std::exp(z_minus), &result);
                if ( ! GSL_SUCCESS == ret_code)
                    throw InternalError("LogLikelihoodBlock::LogGamma: cannot evaluate cumulative at lambda = " + stringify(lambda)
                                        + ", alpha = " + stringify(alpha) + ". GSL reports: " + gsl_strerror(ret_code)
                                        + ". Perhaps the input is too [a]symmetric?");
                const double cdf_minus = result.val;

                const double second = std::fabs((cdf_plus - cdf_minus) - 0.68268949213708585);

                return first + second;
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

                if ( ! status == GSL_SUCCESS)
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

            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                ObservablePtr observable = this->cache.observable(id)->clone(cache.parameters());

                return LogLikelihoodBlockPtr(new LogGammaBlock(cache, cache.add(observable),
                    central - sigma_lower, central, central + sigma_upper, lambda, alpha, _number_of_observations));
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
                std::string name = cache.observable(id)->name();
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
                // and hopefully it is independent of the best fit parameters chosen
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

                if ( ! status == GSL_SUCCESS)
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
            mutable std::vector<double> temp;

            MixtureBlock(const std::vector<LogLikelihoodBlockPtr> & components, 
                         const std::vector<double> & weights) :
                components(components),
                weights(weights),
                temp(weights.size(), 0.0)
            {
            }

            ~MixtureBlock()
            {
            }

            std::string as_string() const
            {
                std::string ret_val = "Mixture: \n";
                for (auto c = components.cbegin() ; c != components.cend() ; ++c)
                    ret_val += (**c).as_string() + '\n';
                return ret_val;
            }

            LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                std::vector<LogLikelihoodBlockPtr> clones;
                for (auto c = components.cbegin() ; c != components.cend() ; ++c)
                    clones.push_back((**c).clone(cache));

                return LogLikelihoodBlockPtr(new MixtureBlock(clones, weights));
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
                for (auto w = weights.cbegin(); w != weights.end() ; ++w, ++v)
                {
                    ret_val += *w * std::exp(*v - *max_val);
                }

                ret_val = std::log(ret_val) + *max_val;

                return ret_val;
            }

            unsigned number_of_observations() const
            {
                unsigned ret_val = 0;
                for (auto c = components.cbegin() ; c != components.cend() ; ++c)
                    ret_val += (**c).number_of_observations();
                return ret_val;
            }

            double sample(gsl_rng * /*rng*/) const
            {
                throw InternalError("LogLikelihoodBlock::MixtureBlock::sample() not implemented yet");
            }

            double significance() const
            {
                throw InternalError("LogLikelihoodBlock::MixtureBlock::significance() not implemented yet");
            }
        };

        template <std::size_t n_>
        struct MultivariateGaussianBlock :
            public LogLikelihoodBlock
        {
            ObservableCache cache;

            std::array<ObservableCache::Id, n_> ids;

            std::array<double, n_> mean;
            std::array<std::array<double, n_>, n_> covariance;
            std::array<std::array<double, n_>, n_> covariance_inv;

            // the normalization constant of the density
            double norm;

            // cholesky matrix of covariance
            std::array<std::array<double, n_>, n_> chol;

            unsigned _number_of_observations;

            MultivariateGaussianBlock(const ObservableCache & cache, const std::array<ObservableCache::Id, n_> & ids,
                                     const std::array<double, n_> & mean, const std::array<std::array<double, n_>, n_> & covariance,
                                     const unsigned & number_of_observations) :
                cache(cache),
                ids(ids),
                mean(mean),
                covariance(covariance),
                norm(compute_norm()),
                _number_of_observations(number_of_observations)
            {
                unsigned k = mean.size();

                if ( ids.size() != k || covariance.size() != k|| covariance.front().size() != k)
                    throw InternalError("MultivariateGaussianBlock.ctor: dimensions of observables, mean and covariance not aligned");

                // cholesky decomposition (informally: the sqrt of the covariance matrix)
                // the GSL matrix contains both the cholesky and its transpose, see gsl ref, ch. 14.5
                gsl_matrix * chol = cholesky();
                invert_covariance(chol);

                // copy only lower and diagonal part
                // set upper to zero
                for (unsigned i = 0 ; i < k ; ++i)
                {
                    for (unsigned j = 0 ; j < k ; ++j)
                    {
                        this->chol[i][j] = (j>i) ? 0 : gsl_matrix_get(chol, i, j);
                    }
                }
                gsl_matrix_free(chol);
            }

            virtual ~MultivariateGaussianBlock()
            {
            }

            virtual std::string as_string() const
            {
                std::string result = "Multivariate Gaussian: ";
                result += "means = ( ";
                for (std::size_t i = 0 ; i < n_ ; ++i)
                {
                    result += stringify(mean[i]) + " ";
                }
                result += "), covariance matrix = (";
                for (std::size_t i = 0 ; i < n_ ; ++i)
                {
                    result += "( ";
                    for (std::size_t j = 0 ; j < n_ ; ++j)
                    {
                        result += stringify(covariance[i][j]) + " ";
                    }
                    result += ")";
                }
                result += "), inverse covariance matrix = (";
                for (std::size_t i = 0 ; i < n_ ; ++i)
                {
                    result += "( ";
                    for (std::size_t j = 0 ; j < n_ ; ++j)
                    {
                        result += stringify(covariance_inv[i][j]) + " ";
                    }
                    result += ")";
                }
                result += " )";

                if (0 == _number_of_observations)
                    result += "; no observation";

                return result;
            }

            // compute cholesky decomposition of covariance matrix
            gsl_matrix * cholesky()
            {
                // dimensionality of parameter space
                unsigned k = mean.size();

                // copy covariance matrix
                gsl_matrix * L = gsl_matrix_alloc (k, k);
                for (unsigned i = 0 ; i < k ; ++i)
                {
                    for (unsigned j = 0 ; j < k ; ++j)
                        gsl_matrix_set(L, i, j, covariance[i][j]);
                }

                gsl_linalg_cholesky_decomp(L);

                return L;
            }

            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                std::array<ObservableCache::Id, n_> indices;

                // add observables to cache
                for (auto i = 0u ; i < n_ ; ++i)
                {
                    indices[i] = cache.add(this->cache.observable(this->ids[i])->clone(cache.parameters()));
                }

                return LogLikelihoodBlockPtr(new implementation::MultivariateGaussianBlock<n_>(cache, indices, mean, covariance, _number_of_observations));
            }

            // compute the normalization constant on log scale
            // -k/2 * log 2 Pi - 1/2 log(abs(det(V^{-1})))
            double compute_norm()
            {
                // dimensionality of parameter space
                unsigned k = mean.size();

                // copy covariance matrix
                gsl_matrix * m = gsl_matrix_alloc (k, k);
                for (unsigned i = 0 ; i < k ; ++i)
                {
                    for (unsigned j = 0 ; j < k ; ++j)
                        gsl_matrix_set(m, i, j, covariance[i][j]);
                }

                // find LU decomposition
                gsl_permutation * p = gsl_permutation_alloc(k);
                int signum = 0;
                gsl_linalg_LU_decomp(m, p, &signum);

                // calculate determinant
                double log_det = gsl_linalg_LU_lndet(m);

                gsl_permutation_free(p);
                gsl_matrix_free(m);

                return -0.5 * k * std::log(2 * M_PI) - 0.5 * log_det;
            }

            virtual double evaluate() const
            {
                std::array<double, n_> observables;

                // read out observable values
                auto o = observables.begin();
                for (auto i = ids.begin(), i_end = ids.end() ; i != i_end ; ++i, ++o)
                {
                    *o = cache[*i];
                }

                // center the gaussian
                observables = observables - mean;

                return norm - dot(observables, covariance_inv * observables) / 2.0;
            }

            void invert_covariance(gsl_matrix * chol)
            {
                if ( ! chol)
                    throw InternalError("MultivariateGaussianBlock.invert_covariance: cholesky decomposition undefined.");

                // dimensionality of parameter space
                unsigned k = mean.size();

                // copy cholesky matrix
                gsl_matrix * inverse = gsl_matrix_alloc (k, k);
                gsl_matrix_memcpy(inverse, chol);

                // compute inverse matrix from cholesky
                gsl_linalg_cholesky_invert (inverse);

                // copy elements
                for (unsigned i = 0 ; i < k ; ++i)
                {
                    for (unsigned j = 0 ; j < k ; ++j)
                        covariance_inv[i][j] = gsl_matrix_get(inverse, i, j);
                }

                gsl_matrix_free(inverse);
            }

            virtual unsigned number_of_observations() const
            {
                return n_;
            }

            virtual double sample(gsl_rng * rng) const
            {
                // Holds the samples in n dimensions
                std::array<double, n_> random_samples;

                // generate standard normals
                std::generate(random_samples.begin(), random_samples.end(), std::bind(gsl_ran_ugaussian, rng));

                // transform
                random_samples = chol * random_samples;

                return norm - 0.5 * dot(random_samples, covariance_inv * random_samples);
            }

            virtual double significance() const
            {
                std::array<double, n_> observables;

                // read out observable values
                auto o = observables.begin();
                for (auto i = ids.begin(), i_end = ids.end() ; i != i_end ; ++i, ++o)
                {
                    *o = cache[*i];
                }

                // center the gaussian
                observables = observables - mean;

                double chi_squared = dot(observables, covariance_inv * observables);

                // find probability of this excess or less ( 1 - usual p-value)
                double p = gsl_cdf_chisq_P(chi_squared, n_);

                // transform to standard Gaussian sigma units
                // return  significance is >= 0, since p >= 0
                // and a negative significance is ruled out by definition
                return gsl_cdf_ugaussian_Pinv((p + 1) / 2.0);
            }
        };

        // explicit instantiation
        template struct MultivariateGaussianBlock<2>;
        template struct MultivariateGaussianBlock<3>;
        template struct MultivariateGaussianBlock<4>;
        template struct MultivariateGaussianBlock<6>;
    }

    LogLikelihoodBlock::~LogLikelihoodBlock()
    {
    }

    void
    LogLikelihoodBlock::prepare_sampling()
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
            const unsigned & number_of_observations)
    {
        // check input
        if (min >= central)
            throw InternalError("LogLikelihoodBlock::LogGamma: min value >= central value");

        if (max <= central)
            throw InternalError("LogLikelihoodBlock::LogGamma: max value <= central value");

        unsigned index = cache.add(observable);

        return LogLikelihoodBlockPtr(new implementation::LogGammaBlock(cache, index, min, central, max, number_of_observations));
    }

    LogLikelihoodBlockPtr
    LogLikelihoodBlock::LogGamma(ObservableCache cache, const ObservablePtr & observable,
            const double & min, const double & central, const double & max,
            const double & lambda, const double & alpha,
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

        return LogLikelihoodBlockPtr(new implementation::LogGammaBlock(cache, index, min, central, max, lambda, alpha, number_of_observations));
    }

    LogLikelihoodBlockPtr
    LogLikelihoodBlock::AmorosoLimit(ObservableCache cache, const ObservablePtr & observable,
                    const double & physical_limit, const double & upper_limit_90, const double & upper_limit_95,
                    const double & theta, const double & alpha,
                    const unsigned & number_of_observations)
    {
        // check input
        if (upper_limit_90 <= physical_limit)
            throw InternalError("LogLikelihoodBlock::AmorosoLimit: upper_limit_90 <= physical_limit");

        if (upper_limit_95 <= physical_limit)
            throw InternalError("LogLikelihoodBlock::AmorosoLimit: upper_limit_95 <= physical_limit");

        if (upper_limit_95 <= upper_limit_90)
            throw InternalError("LogLikelihoodBlock::AmorosoLimit: upper_limit_95 <= upper_limit_90");

        unsigned index = cache.add(observable);

        implementation::AmorosoBlock * a = new implementation::AmorosoBlock(cache, index, physical_limit, theta, alpha, 1 / alpha, number_of_observations);

        // check consistency
        if (std::abs(a->cdf(upper_limit_90) - 0.90) > 1e-4)
        {
            throw InternalError("LogLikelihood::AmorosoLimit.ctor: For the current parameter values, cdf(x_90) = "
                + stringify(a->cdf(upper_limit_90)) + " deviates from 90%.");

        }
        if (std::abs(a->cdf(upper_limit_95) - 0.95) > 1e-4)
        {
            throw InternalError("LogLikelihood::AmorosoLimit.ctor: For the current parameter values, cdf(x_95) = "
                + stringify(a->cdf(upper_limit_95)) + " deviates from 95%.");
        }

        return LogLikelihoodBlockPtr(a);
    }

    LogLikelihoodBlockPtr
    LogLikelihoodBlock::AmorosoMode(ObservableCache cache, const ObservablePtr & observable,
                                const double & physical_limit, const double & mode,
                                const double & upper_limit_90, const double & upper_limit_95,
                                const double & theta, const double & alpha, const double & beta,
                                const unsigned & number_of_observations)
    {
        // check input
        if (mode <= physical_limit)
            throw InternalError("LogLikelihoodBlock::AmorosoLimit: upper_limit_5 <= physical_limit");

        if (upper_limit_90 <= physical_limit)
            throw InternalError("LogLikelihoodBlock::AmorosoLimit: upper_limit_90 <= physical_limit");

        if (upper_limit_95 <= upper_limit_90)
            throw InternalError("LogLikelihoodBlock::AmorosoLimit: upper_limit_95 <= upper_limit_90");

        unsigned index = cache.add(observable);

        implementation::AmorosoBlock * a = new implementation::AmorosoBlock(cache, index, physical_limit, theta, alpha, beta, number_of_observations);

        // check consistency
        if (std::abs(a->mode() - mode) > 1e-4)
        {
            throw InternalError("LogLikelihood::Amoroso.ctor: For the current parameter values, Amoroso::mode() = "
                + stringify(a->mode()) + " deviates from mode supplied " + stringify(mode));
        }
        if (std::abs(a->cdf(upper_limit_90) - 0.90) > 1e-4)
        {
            throw InternalError("LogLikelihood::Amoroso.ctor: For the current parameter values, cdf(x_90) = "
                + stringify(a->cdf(upper_limit_90)) + " deviates from 90%.");
        }
        if (std::abs(a->cdf(upper_limit_95) - 0.95) > 1e-4)
        {
            throw InternalError("LogLikelihood::Amoroso.ctor: For the current parameter values, cdf(x_95) = "
                + stringify(a->cdf(upper_limit_95)) + " deviates from 95%.");
        }

        return LogLikelihoodBlockPtr(a);
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
                                const std::vector<double> & weights)
    {
        if (components.size() != weights.size())
            throw InternalError("LogLikelihoodBlock::Mixture(): components and weights don't match");

        // normalize weights
        double sum = 0; // std::accumulate(weights.cbegin(), weights.cend(), 0);
        for (auto w = weights.begin() ; w != weights.end() ; ++w)
            sum += *w;

        std::vector<double> norm_weights(weights);

        std::transform(weights.cbegin(), weights.cend(), norm_weights.begin(),
                       std::bind(std::multiplies<double>(), 1.0 / sum, std::placeholders::_1));
        Log::instance()->message("MixtureBlock()", ll_debug)
            << "sum = " << sum <<  ", norm. weights " << stringify_container(norm_weights);

        return LogLikelihoodBlockPtr(new implementation::MixtureBlock(components, norm_weights));
    }

    template <std::size_t n_>
    LogLikelihoodBlockPtr
    LogLikelihoodBlock::MultivariateGaussian(ObservableCache cache, const std::array<ObservablePtr, n_> & observables,
                                             const std::array<double, n_> & mean, const std::array<std::array<double, n_>, n_> & covariance,
                                             const unsigned & number_of_observations)
    {
        std::array<ObservableCache::Id, n_> indices;

        // add observables to cache
        for (auto i = 0u ; i < n_ ; ++i)
        {
            indices[i] = cache.add(observables[i]);
        }

        return LogLikelihoodBlockPtr(new implementation::MultivariateGaussianBlock<n_>(cache, indices, mean, covariance, number_of_observations));
    }

    // explicit instantiation
    template LogLikelihoodBlockPtr LogLikelihoodBlock::MultivariateGaussian<2>(ObservableCache cache, const std::array<ObservablePtr, 2> & observables,
                                             const std::array<double, 2> & mean, const std::array<std::array<double, 2>, 2> & covariance,
                                             const unsigned & number_of_observations = 2u);
    template LogLikelihoodBlockPtr LogLikelihoodBlock::MultivariateGaussian<3>(ObservableCache cache, const std::array<ObservablePtr, 3> & observables,
                                             const std::array<double, 3> & mean, const std::array<std::array<double, 3>, 3> & covariance,
                                             const unsigned & number_of_observations = 3u);
    template LogLikelihoodBlockPtr LogLikelihoodBlock::MultivariateGaussian<4>(ObservableCache cache, const std::array<ObservablePtr, 4> & observables,
                                             const std::array<double, 4> & mean, const std::array<std::array<double, 4>, 4> & covariance,
                                             const unsigned & number_of_observations = 4u);
    template LogLikelihoodBlockPtr LogLikelihoodBlock::MultivariateGaussian<6>(ObservableCache cache, const std::array<ObservablePtr, 6> & observables,
                                             const std::array<double, 6> & mean, const std::array<std::array<double, 6>, 6> & covariance,
                                             const unsigned & number_of_observations = 6u);

    template <std::size_t n_>
    LogLikelihoodBlockPtr
    LogLikelihoodBlock::MultivariateGaussian(ObservableCache cache, const std::array<ObservablePtr, n_> & observables,
                                                      const std::array<double, n_> & mean, const std::array<double, n_> & variances,
                                                      const std::array<std::array<double, n_>, n_> & correlation,
                                                      const unsigned & number_of_observations)
    {
        // convert correlation matrix to covariance matrix
        std::array<std::array<double, n_>, n_> covariance(correlation);
        for (auto i = 0u ; i < n_ ; ++i)
        {
            // diagonal elements are just the variances
            covariance[i][i] = variances[i];

            for (auto j = i + 1 ; j < n_ ; ++j)
            {
                covariance[i][j] = covariance[j][i] = correlation[i][j] * std::sqrt(variances[i] * variances[j]);
            }
        }

        return MultivariateGaussian(cache, observables, mean, covariance, number_of_observations);
    }

    // explicit instantiation
    template LogLikelihoodBlockPtr LogLikelihoodBlock::MultivariateGaussian<2>(ObservableCache cache, const std::array<ObservablePtr, 2> & observables,
                                             const std::array<double, 2> & mean, const std::array<double, 2> & variances,
                                             const std::array<std::array<double, 2>, 2> & correlation,
                                             const unsigned & number_of_observations = 2u);
    template LogLikelihoodBlockPtr LogLikelihoodBlock::MultivariateGaussian<3>(ObservableCache cache, const std::array<ObservablePtr, 3> & observables,
                                             const std::array<double, 3> & mean, const std::array<double, 3> & variances,
                                             const std::array<std::array<double, 3>, 3> & correlation,
                                             const unsigned & number_of_observations = 3u);
    template LogLikelihoodBlockPtr LogLikelihoodBlock::MultivariateGaussian<4>(ObservableCache cache, const std::array<ObservablePtr, 4> & observables,
                                             const std::array<double, 4> & mean, const std::array<double, 4> & variances,
                                             const std::array<std::array<double, 4>, 4> & correlation,
                                             const unsigned & number_of_observations = 4u);
    template LogLikelihoodBlockPtr LogLikelihoodBlock::MultivariateGaussian<6>(ObservableCache cache, const std::array<ObservablePtr, 6> & observables,
                                             const std::array<double, 6> & mean, const std::array<double, 6> & variances,
                                             const std::array<std::array<double, 6>, 6> & correlation,
                                             const unsigned & number_of_observations = 6u);

    template <>
    struct Implementation<LogLikelihood>
    {
        Parameters parameters;

        // Cache observable predictions
        ObservableCache cache;

        // Container for all named constraints
        std::vector<Constraint> constraints;

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

            // set up for sampling
            for (auto c = constraints.cbegin(), c_end = constraints.cend() ; c != c_end ; ++c)
            {
                for (auto b = c->begin_blocks(), b_end = c->end_blocks() ; b != b_end ; ++b)
                {
                    (*b)->prepare_sampling();
                }
            }

            // observed value
            double t_obs = this->log_likelihood();

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
                for (auto c = constraints.cbegin(), c_end = constraints.cend() ; c != c_end ; ++c)
                {
                    for (auto b = c->begin_blocks(), b_end = c->end_blocks() ; b != b_end ; ++b)
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

        double log_likelihood()
        {
            double result = 0.0;

            // loop over all likelihood blocks
            for (auto c = constraints.cbegin(), c_end = constraints.cend() ; c != c_end ; ++c)
            {
                for (auto b = c->begin_blocks(), b_end = c->end_blocks() ; b != b_end ; ++b)
                {
                    result += (*b)->evaluate();
                }
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

    template class WrappedForwardIterator<LogLikelihood::ConstraintIteratorTag, Constraint>;

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

        for (auto c = _imp->constraints.cbegin(), c_end = _imp->constraints.cend() ; c != c_end ; ++c)
        {
            result.add(*c);
        }

        return result;
    }

    unsigned
    LogLikelihood::number_of_observations() const
    {
        unsigned result = 0;
        for (auto c = _imp->constraints.cbegin(), c_end = _imp->constraints.cend() ; c != c_end ; ++c)
        {
            for (auto b = c->begin_blocks(), b_end = c->end_blocks() ; b != b_end ; ++b)
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
    LogLikelihood::operator() ()
    {
        _imp->cache.update();

        return _imp->log_likelihood();
    }

    double
    LogLikelihood::operator() (const Parameter::Id & id)
    {
        _imp->cache.update(id);

        return _imp->log_likelihood();
    }

    void
    LogLikelihood::reset()
    {
        _imp->cache.reset();
    }
}
