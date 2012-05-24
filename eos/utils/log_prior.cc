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

#include <eos/utils/log_prior.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/equation_solver.hh>
#include <eos/utils/log.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <limits>
#include <numeric>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_spline.h>

namespace eos
{
    template class WrappedForwardIterator<LogPrior::IteratorTag, ParameterDescription>;

    namespace priors
    {
        struct RangeError :
            public Exception
        {
            RangeError(const std::string & message) throw () :
                Exception("Range Error: " + message)
            {
            }
        };

        struct UnknownPriorError :
            public Exception
        {
            UnknownPriorError(const std::string & message) throw () :
                Exception("Unknown prior error: " + message)
            {
            }
        };

        /*!
         * Flat or uniform prior
         */
        class Flat :
            public LogPrior
        {
            private:
                std::string _name;

                ParameterRange _range;

                // The flat prior always returns this value
                double _value;

            public:
                Flat(const Parameters & parameters, const std::string & name, const ParameterRange & range) :
                    LogPrior(parameters),
                    _name(name),
                    _range(range),
                    _value(std::log(1.0 / (range.max - range.min)))
                {
                    if (range.min >= range.max)
                    {
                        throw RangeError("LogPrior::Flat(" + _name +"): minimum (" + stringify(range.min)
                                          + ") must be smaller than maximum (" + stringify(range.max) + ")");
                    }
                    _parameter_descriptions.push_back(ParameterDescription{ _parameters[name], range.min, range.max, false, false });
                }

                virtual ~Flat()
                {
                }

                virtual std::string as_string() const
                {
                    std::string result = "Parameter: " + _name + ", prior type: flat, range: [" + stringify(_range.min) + "," + stringify(_range.max) + "]";
                    return result;
                }

                virtual double operator()() const
                {
                    return _value;
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    return LogPriorPtr(new priors::Flat(parameters, _name, _range));
                }

                virtual double sample(gsl_rng * rng) const
                {
                    return gsl_rng_uniform(rng) * (_range.max - _range.min) + _range.min;
                }

                virtual double mean() const
                {
                    return (_range.max - _range.min) / 2.0;
                }

                virtual double variance() const
                {
                    return power_of<2>(_range.max - _range.min) / 12.0;
                }
        };

        /*!
         * [asymmetric] Gaussian or Normal prior distribution
         */
        class Gauss :
            public LogPrior
        {
            private:
                std::string _name;

                ParameterRange _range;

                double _lower, _central, _upper;

                double _sigma_lower, _sigma_upper;

                double _norm_lower, _norm_upper;

                // the probability covered to the left of the central value
                double _prob_lower;

                // coefficients needed for sampling from asymmetric Gaussian on finite support
                // the cumulative is a piecewise function
                // CDF(x) = CDF_lower(x, sigma_lower) if x < central, else CDF_upper(x, sigma_upper)
                // To ensure that cumulative is
                // a) continuous at the central value
                // b) zero when x < x_min
                // c) one when  x > x_max
                // d) relative prob of upper vs. lower part is given by ratio of standard Gaussian cumulatives from [x_min, x_c] and [x_c, x_max]
                // need to fix the coefficients in
                // @f$CDF_{lower}(x) = c_{lower} * ( \Phi(x - x_c / \sigma_{lower}) + b_lower $@f
                double _b_lower, _c_lower;

            public:
                Gauss(const Parameters & parameters, const std::string & name, const ParameterRange & range,
                        const double & lower, const double & central, const double & upper) :
                    LogPrior(parameters),
                    _name(name),
                    _range(range),
                    _lower(lower),
                    _central(central),
                    _upper(upper),
                    _sigma_lower(central - lower),
                    _sigma_upper(upper - central)
                {
                    if (range.min >= range.max)
                    {
                        throw RangeError("LogPrior::Gauss(" + _name +"): minimum (" + stringify(range.min)
                                          + ") must be smaller than maximum (" + stringify(range.max) + ")");
                    }
                    _parameter_descriptions.push_back(ParameterDescription{ _parameters[name], range.min, range.max, false, false });

                    // scale factor takes finite range into account. For large range, it is 1
                    _c_lower = 1.0 / (gsl_cdf_gaussian_P(_range.max - _central, _sigma_upper)
                                    - gsl_cdf_gaussian_P(_range.min - _central, _sigma_lower));
                    _b_lower = -1.0 * gsl_cdf_gaussian_P(_range.min - _central, _sigma_lower) * _c_lower;

                    _norm_lower =  std::log(_c_lower / (std::sqrt(2.0 * M_PI) * _sigma_lower));
                    _norm_upper =  std::log(_c_lower/ (std::sqrt(2.0 * M_PI) * _sigma_upper));

                    _prob_lower = _c_lower / 2.0 + _b_lower;

                    // sanity check: cdf must be continuous at the central point
                    if (std::fabs((_c_lower / 2.0 + _b_lower) - (_c_lower / 2.0 + _b_lower)) > 1e-12)
                        throw InternalError("LogPrior::Gauss: cdf not continuous");
                }

                virtual ~Gauss()
                {
                }

                virtual std::string as_string() const
                {
                    std::string result = "Parameter: " + _name + ", prior type: Gaussian, range: [" + stringify(_range.min) + "," + stringify(_range.max) + "]";

                    result += ", x = " + stringify(_central);
                    if (std::abs(_sigma_upper - _sigma_lower) < 1e-15)
                    {
                        result += " +- " + stringify(_sigma_upper);
                    }
                    else
                    {
                        result += " + " + stringify(_sigma_upper) + " - " + stringify(_sigma_lower);
                    }
                    return result;
                }

                virtual double operator()() const
                {
                    double sigma = 0.0, norm = 0.0;

                    // read parameter's current value
                    double x = _parameter_descriptions.front().parameter;

                    if (x < _central)
                    {
                        sigma = _sigma_lower;
                        norm = _norm_lower;
                    }
                    else
                    {
                        sigma = _sigma_upper;
                        norm = _norm_upper;
                    }

                    return norm - power_of<2>((x - _central) / sigma) / 2.0;
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    return LogPriorPtr(new priors::Gauss(parameters, _name, _range, _lower, _central, _upper));
                }

                virtual double sample(gsl_rng * rng) const
                {
                    // find out if sample in upper or lower part
                    double u = gsl_rng_uniform(rng);

                    // get a sample from lower or upper part using inverse transform method
                    // CDF = c \Phi(x - x_{central} / \sigma) + b
                    if (u < _prob_lower)
                        return gsl_cdf_ugaussian_Pinv((u - _b_lower) / _c_lower) * _sigma_lower + _central;
                    else
                        return gsl_cdf_ugaussian_Pinv((u - _b_lower) / _c_lower) * _sigma_upper + _central;
                }

                virtual double mean() const
                {
                    return _central;
                }

                ///Only true if parameter range is the whole real line
                virtual double variance() const
                {
                    return (power_of<2>(_sigma_upper) + power_of<2>(_sigma_lower)) / 2.0;
                }
        };

        /*!
         * [asymmetric] log-gamma prior distribution
         * Useful to input information from another paper stating a quantity is known to be
         * x = 1 + 0.20 - 0.15
         *
         * For symmetric uncertainties, one should always use the Gaussian distribution,
         * but in the asymmetric case, it may be desirable to have a smooth distribution everywhere,
         * whose cumulative, F(x), satisfies the following conditions:
         *
         * 1. F(1 + 0.20) == 0.84134
         * 2. F(1 - 0.15) == 0.15865
         *
         * This way, the interval contains the ominous 68.32 % probability,
         * as a Gaussian distribution does in [\mu - \sigma, \mu + \sigma].
         *
         * In addition, the mode is required to be at 1
         *
         * 3. f'(x=1) = 0
         *
         * Given the three conditions, the three parameters \nu, \lambda, \alpha are uniquely determined,
         * they are found by numerical optimization.
         *
         * More details on the distribution can be found in [C2004], Sec. 2.
         */
        class LogGamma :
            public LogPrior
        {
            private:
                std::string _name;

                ParameterRange _range;

                double _central, _sigma_lower, _sigma_upper;
                double _sigma_plus, _sigma_minus;

                double _nu, _lambda, _alpha;

                double _norm;

                LogGamma(const Parameters & parameters) :
                    LogPrior(parameters)
                {
                }

            public:
                LogGamma(const Parameters & parameters, const std::string & name, const ParameterRange & range,
                        const double & lower, const double & central, const double & upper) :
                    LogPrior(parameters),
                    _name(name),
                    _range(range),
                    _central(central),
                    _sigma_lower(central - lower),
                    _sigma_upper(upper - central),
                    _sigma_plus( (_sigma_upper > _sigma_lower)? _sigma_upper / _sigma_lower : _sigma_lower / _sigma_upper),
                    _sigma_minus(1.0)
                {
                    if (range.min >= range.max)
                    {
                        throw RangeError("LogPrior::LogGamma(" + _name +"): minimum (" + stringify(range.min)
                                          + ") must be smaller than maximum (" + stringify(range.max) + ")");
                    }
                    _parameter_descriptions.push_back(ParameterDescription{ _parameters[name], range.min, range.max, false, false });

                    // avoid extrapolation from polynomial
                    if (_sigma_plus < 1.03)
                    {
                        throw InternalError("LogPrior::LogGammaFor nearly symmetric uncertainties (" + stringify(_sigma_lower) + " vs " + stringify(_sigma_upper)
                            + "), this procedure fails to find the correct parameter values. Please use a Gaussian prior instead.");
                    }

                    // for positive skew, \lambda is negative
                    // in the fit, \lambda always considered negative, so it only changes sign for negative skew
                    const double lambda_scale_factor = _sigma_upper > _sigma_lower ? _sigma_lower / _sigma_minus : -1.0 * _sigma_upper / _sigma_minus;

                    /* find the parameters using good starting values. Assume upper > lower=1, and fix sign at the end */

                    // solved constraints for particular values of sigma_plus, interpolate linearly to find good starting position
                    // the points are called knots in the spline jargon
                    static const std::vector<double> knots_sigma  {  1.03,  1.04,  1.05,  1.06,  1.1 , 1.15,  1.2,  1.3 ,  1.6 ,  1.8 ,  1.9 ,  2.0 ,  2.5 ,  3.2 ,  4.0 ,  5.0  , 10.0};
                    static const std::vector<double> knots_lambda { -12.4, -8.70, -7.00, -5.90, -3.67, -2.6, -2.0, -1.44, -0.88, -0.73, -0.69, -0.65, -0.53, -0.45, -0.39, -0.35 , -0.27 };
                    static const std::vector<double> knots_alpha  {  127.,  72.4,  46.9,  32.9,  12.4,  5.9,  3.5,  1.78,  0.64,  0.44,  0.38,  0.33,  0.21,  0.15,  0.10,  0.073,  0.029 };

                    if (_sigma_plus > knots_sigma.back())
                    {
                        Log::instance()->message("LogPrior::LogGamma::ctor", ll_warning)
                            << "Asymmetry " << _sigma_plus << " very large. Extrapolation required. GSL behavior undefined";
                    }
                    gsl_interp * spline = gsl_interp_alloc(gsl_interp_linear, knots_sigma.size());
                    gsl_interp_accel * accelerator = NULL;
                    gsl_interp_init(spline, &knots_sigma[0], &knots_lambda[0], knots_sigma.size());
                    const double lambda_initial = gsl_interp_eval(spline, &knots_sigma[0], &knots_lambda[0], _sigma_plus, accelerator);

                    gsl_interp_init(spline, &knots_sigma[0], &knots_alpha[0], knots_sigma.size());
                    const double alpha_initial  = gsl_interp_eval(spline, &knots_sigma[0], &knots_alpha[0] , _sigma_plus, accelerator);
                    gsl_interp_free(spline);

                    EquationSolver solver(EquationSolver::Config::Default());

                    //add free Parameter: initial value, error
                    solver.add("lambda", lambda_initial, lambda_initial / 5.0, -30.0, 0.0);
                    //add positive Parameter: initial value, error, min, max. \alpha for 5% asymmetry at 500, so 1000 is well above
                    solver.add("alpha", alpha_initial, alpha_initial / 5.0, 0.0, 1000);

                    // add constraints for standardized problem
                    solver.add(std::bind(&LogGamma::constraint, *this, std::placeholders::_1));

                    // check errors manually, then restore default behavior later
                    gsl_error_handler_t * default_gsl_error_handler = gsl_set_error_handler_off();

                    // find the solution
                    auto solution = solver.solve();

                    // global minimum at zero value, often Minuit claims to not have found it while it actually did
                    if ( ((! solution.valid) && solution.value > 1e-4) || solution.parameters[1] > 500)
                    {
                        Log::instance()->message("LogPrior::LogGamma.ctor", ll_informational)
                            << "Standardized: nu = " <<  solution.parameters[0] * std::log(solution.parameters[1])
                            << ", lambda = " << solution.parameters[0]
                            << ", alpha = " << solution.parameters[1]
                            << ", solution = " <<  solution.value
                            << ", valid = " << solution.valid;

                        throw InternalError("Solution of constraints for '" + _name + "' failed");
                    }

                    // now we have all values
                    _lambda = lambda_scale_factor * solution.parameters[0];
                    _alpha  = solution.parameters[1];
                    _nu = central - _lambda * std::log(_alpha);

                    // account for finite range. Multiply with (CDF(max) - CDF(min))^{-1}
                    // CDF (\lambda > 0) = 1 - CDF(\lambda < 0)
                    _norm =   gsl_sf_gamma_inc_Q(_alpha, std::exp( (_range.max - _nu) / _lambda))
                            - gsl_sf_gamma_inc_Q(_alpha, std::exp( (_range.min - _nu) / _lambda));
                    if (lambda_scale_factor < 0)
                        _norm *= -1.0;
                    _norm = -1.0 * std::log(_norm);

                    // calculate normalization factors that are independent of x
                    _norm += -1.0 * gsl_sf_lngamma(_alpha) - std::log(std::fabs(_lambda));

                    // restore default error handler
                    gsl_set_error_handler(default_gsl_error_handler);
                }

                virtual ~LogGamma()
                {
                }

                virtual std::string as_string() const
                {
                    std::string result = "Parameter: " + _name + ", prior type: LogGamma, range: [" + stringify(_range.min) + "," + stringify(_range.max) + "]";

                    result += ", x = " + stringify(_central);
                    result += " + " + stringify(_sigma_upper) + " - " + stringify(_sigma_lower);
                    result += ", nu: " + stringify(_nu) + ", lambda: " + stringify(_lambda);
                    result += ", alpha: " + stringify(_alpha);

                    return result;
                }

                // optimize parameters such that the two constraints are met
                double constraint(const std::vector<double> & parameter_values) const
                {
                    gsl_sf_result result;
                    const double & lambda = parameter_values[0];
                    const double & alpha  = parameter_values[1];

                    // standardized mode at 0
                    const double nu = 0.0 - lambda * std::log(alpha);

                    // standardized coordinates at plus/minus
                    const double z_plus = (_sigma_plus - nu) / lambda;
                    const double z_minus = (-1.0 *_sigma_minus - nu) / lambda;

                    // first constraint: pdf's should be equal, neglect prefactors
                    const double first = std::fabs(alpha * z_plus - std::exp(z_plus) - alpha * z_minus + std::exp(z_minus));

                    // second constraint: 68% interval
                    int ret_code = gsl_sf_gamma_inc_Q_e(alpha, std::exp(z_plus), &result);
                    if ( ! GSL_SUCCESS == ret_code)
                        throw InternalError("LogPrior::LogGamma: cannot evaluate cumulative for '" + _name + "' at lambda = " + stringify(lambda)
                                            + ", alpha = " + stringify(alpha) + ". GSL reports: " + gsl_strerror(ret_code)
                                            + ". Perhaps the input is too [a]symmetric?");
                    const double cdf_plus = result.val;

                    ret_code = gsl_sf_gamma_inc_Q_e(alpha, std::exp(z_minus), &result);
                    if ( ! GSL_SUCCESS == ret_code)
                        throw InternalError("LogPrior::LogGamma: cannot evaluate cumulative for '" + _name + "' at lambda = " + stringify(lambda)
                                            + ", alpha = " + stringify(alpha) + ". GSL reports: " + gsl_strerror(ret_code)
                                            + ". Perhaps the input is too [a]symmetric?");
                    const double cdf_minus = result.val;

                    const double second = std::fabs((cdf_plus - cdf_minus) - 0.68268949213708585);

                    return first + second;
                }

                virtual double operator()() const
                {
                    const double z = (_parameter_descriptions.front().parameter - _nu) / _lambda;
                    return _norm + _alpha * z - std::exp(z);
                }

                // change private members by hand. Saves time on optimization
                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    priors::LogGamma * log_gamma = new priors::LogGamma(parameters);
                    log_gamma->_alpha = _alpha;
                    log_gamma->_central = _central;
                    log_gamma->_lambda = _lambda;
                    log_gamma->_name = _name;
                    log_gamma->_norm = _norm;
                    log_gamma->_nu = _nu;
                    log_gamma->_parameter_descriptions.push_back(ParameterDescription{ parameters[_name], _range.min, _range.max, false, false });
                    log_gamma->_range = _range;
                    log_gamma->_sigma_lower = _sigma_lower;
                    log_gamma->_sigma_minus = _sigma_minus;
                    log_gamma->_sigma_plus = _sigma_plus;
                    log_gamma->_sigma_upper = _sigma_upper;
                    return LogPriorPtr(log_gamma);
                }

                /*!
                 * Use the fact that if x' ~ StdLogGamma(\alpha), then
                 * x = \nu + \lambda * x' ~ LogGamma(\nu, \lambda, \alpha)
                 */
                virtual double sample(gsl_rng * rng) const
                {
                    double x = _range.min;

                    while(true)
                    {
                        x = _lambda * std::log(gsl_ran_gamma(rng, _alpha, 1.0)) + _nu;

                        if (_range.min < x && x < _range.max)
                            break;
                    }

                    return x;
                }

                virtual double mean() const
                {
                    gsl_sf_result result;
                    if ( ! GSL_SUCCESS == gsl_sf_psi_e(_alpha, &result))
                    {
                        Log::instance()->message("LogPrior::LogGamma.variance", ll_error)
                            << "Error in evaluating the trigamma function in GSL";
                    }
                    return _nu + _lambda * result.val;
                }

                ///Only true if parameter range is the whole real line
                virtual double variance() const
                {
                    gsl_sf_result result;
                    if ( ! GSL_SUCCESS == gsl_sf_psi_1_e(_alpha, &result))
                    {
                        Log::instance()->message("LogPrior::LogGamma.variance", ll_error)
                            << "Error in evaluating the trigamma function in GSL";
                    }
                    return power_of<2>(_lambda) * result.val;
                }
        };

        class Discrete :
            public LogPrior
        {
            private:
                std::string _name;

                std::vector<double> _values;

                double _prob;

            public:
                Discrete(const Parameters & parameters, const std::string & name, const std::set<double> & values) :
                    LogPrior(parameters),
                    _name(name),
                    _values(values.begin(), values.end()),
                    _prob(1.0 / _values.size())
                {
                    _parameter_descriptions.push_back(ParameterDescription{ _parameters[name], _values.front(), _values.back(), false, true });
                }

                virtual ~Discrete()
                {
                }

                virtual std::string as_string() const
                {
                    std::string result = "Parameter: " + _name + ", prior type: discrete, values = ";
                    result += stringify(_values.cbegin(), _values.cend());

                    return result;
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    return LogPriorPtr(new priors::Discrete(parameters, _name, std::set<double>(_values.begin(), _values.end())));
                }

                virtual double operator() () const
                {
                    return _prob;
                }

                virtual double sample(gsl_rng * rng) const
                {
                    return _values[gsl_rng_uniform_int(rng, _values.size())];
                }

                virtual double mean() const
                {
                    throw InternalError("Attempting to obtain mean of a discrete prior");
                    return 0.0;
                }

                virtual double variance() const
                {
                    throw InternalError("Attempting to obtain variance of a discrete prior");
                    return 0.0;
                }
        };
    }

    LogPrior::LogPrior(const Parameters & parameters) :
        _parameters(parameters)
    {
    }

    LogPrior::Iterator
    LogPrior::begin()
    {
        return LogPrior::Iterator(_parameter_descriptions.begin());
    }

    LogPrior::Iterator
    LogPrior::end()
    {
        return LogPrior::Iterator(_parameter_descriptions.end());
    }

    LogPriorPtr
    LogPrior::Discrete(const Parameters & parameters, const std::string & name, const std::set<double> & values)
    {
        LogPriorPtr prior = std::make_shared<eos::priors::Discrete>(parameters, name, values);

        return prior;
    }

    LogPriorPtr
    LogPrior::Flat(const Parameters & parameters, const std::string & name, const ParameterRange & range)
    {
        LogPriorPtr prior = std::make_shared<eos::priors::Flat>(parameters, name, range);

        return prior;
    }

    LogPriorPtr
    LogPrior::Gauss(const Parameters & parameters, const std::string & name, const ParameterRange & range,
            const double & lower, const double & central, const double & upper)
    {
        // check input
        if (lower >= central)
            throw InternalError("LogPrior::Gauss: lower value (" + stringify(lower) + ") >= central value (" + stringify(central) + ")");

        if (upper <= central)
            throw InternalError("LogPrior::Gauss: upper value (" + stringify(upper) + ") <= central value (" + stringify(central) + ")");

        LogPriorPtr prior = std::make_shared<eos::priors::Gauss>(parameters, name, range, lower, central, upper);

        return prior;
    }

    LogPriorPtr
    LogPrior::LogGamma(const Parameters & parameters, const std::string & name, const ParameterRange & range,
                       const double & lower, const double & central, const double & upper)
    {
        // check input
        if (lower >= central)
            throw InternalError("LogPrior::LogGamma: lower value (" + stringify(lower) + ") >= central value (" + stringify(central) + ")");

        if (upper <= central)
            throw InternalError("LogPrior::LogGamma: upper value (" + stringify(upper) + ") <= central value (" + stringify(central) + ")");

        LogPriorPtr prior = std::make_shared<eos::priors::LogGamma>(parameters, name, range, lower, central, upper);

        return prior;
    }

    LogPriorPtr
    LogPrior::Make(const Parameters & parameters, const std::string & s)
    {
        // extract parameter name
        auto loc1 = s.find_first_of(':');
        auto loc2 = s.find_first_of(',');

        std::string par_name = s.substr(loc1 + 2, loc2 - loc1 - 2);

        // extract type
        loc1 = s.find_first_of(':', loc2 + 1);
        loc2 = s.find_first_of(',', loc2 + 1);

        std::string prior_type = s.substr(loc1 + 2, loc2 - loc1 - 2);

        // extract range
        loc1 = s.find_first_of('[', loc2 + 1);
        loc2 = s.find_first_of(',', loc2 + 1);
        std::string range_min = s.substr(loc1 + 1, loc2 - loc1 - 1);

        loc1 = s.find_first_of(',', loc1 + 1);
        loc2 = s.find_first_of(']', loc1 + 1);
        std::string range_max = s.substr(loc1 + 1, loc2 - loc1 - 1);

        ParameterRange range { destringify<double>(range_min), destringify<double>(range_max) };

        if (prior_type == "flat")
        {
            return Flat(parameters, par_name, range);
        }
        if (prior_type == "Gaussian" || prior_type == "LogGamma")
        {
            // extract central value
            loc1 = s.find_first_of('=', loc2 + 1);
            loc2 = s.find_first_of('+', loc2 + 1);

            double central = destringify<double>(s.substr(loc1 + 2, loc2 - loc1 - 2));

            double sigma_upper, sigma_lower;

            // extract sigma_upper, lower
            if ( s[loc2 + 1] == '-')
            {
                sigma_upper = destringify<double>(s.substr(loc2 + 2));
                sigma_lower = sigma_upper;
            }
            else
            {
                loc1 = loc2;
                loc2 = s.find_first_of('-', loc2 + 1);

                sigma_upper = destringify<double>(s.substr(loc1 + 1, loc2 - loc1 - 1));

                loc1 = loc2;
                loc2 = s.find_first_of(',', loc2 + 1);

                // Gaussian
                if (loc2 == std::string::npos)
                {
                    sigma_lower = destringify<double>(s.substr(loc1 + 1));
                }
                // LogGamma
                // don't parse until end of string, but until next comma
                else
                {
                    sigma_lower = destringify<double>(s.substr(loc1 + 1, loc2 - loc1 -1));
                }
            }

            if (prior_type == "Gaussian")
                return Gauss(parameters, par_name, range, central - sigma_lower, central, central + sigma_upper);
            else
                return LogGamma(parameters, par_name, range, central - sigma_lower, central, central + sigma_upper);
        }

        throw priors::UnknownPriorError("Cannot construct prior from '" + s + "'");
    }
}
