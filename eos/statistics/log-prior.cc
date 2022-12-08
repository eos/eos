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

#include <eos/statistics/log-prior.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/log.hh>
#include <eos/maths/power-of.hh>
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
    template <>
    struct WrappedForwardIteratorTraits<LogPrior::IteratorTag>
    {
        using UnderlyingIterator = std::vector<ParameterDescription>::iterator;
    };

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
                Parameter _parameter;

                std::string _name;

                ParameterRange _range;

                // The flat prior always returns this value
                double _value;

            public:
                Flat(const Parameters & parameters, const std::string & name, const ParameterRange & range) :
                    LogPrior(parameters),
                    _parameter(parameters[name]),
                    _name(name),
                    _range(range),
                    _value(std::log(1.0 / (range.max - range.min)))
                {
                    if (range.min >= range.max)
                    {
                        throw RangeError("LogPrior::Flat(" + _name +"): minimum (" + stringify(range.min)
                                          + ") must be smaller than maximum (" + stringify(range.max) + ")");
                    }
                    _parameter_descriptions.push_back(ParameterDescription{ _parameters[name].clone(), range.min, range.max, false });
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

                virtual void sample()
                {
                    _parameter.set(_parameter.evaluate_generator() * (_range.max - _range.min) + _range.min);
                }

                virtual bool informative() const
                {
                    return false;
                }
        };

        /*!
         * [asymmetric] Gaussian or Normal prior distribution with finite support
         */
        class CurtailedGauss :
            public LogPrior
        {
            private:
                Parameter _parameter;

                const std::string _name;

                const ParameterRange _range;

                const double _lower, _central, _upper;

                const double _sigma_lower, _sigma_upper;

                // coefficients needed for sampling from asymmetric Gaussian on finite support
                // the pdf is a piecewise function of y given x^{+a}_{-b}
                // P(y|x,a,b) = \theta(y-x) c_a \mathcal{N}(y|x,a)  + \theta(x-y) c_b \mathcal{N}(y|x,b)
                // Fix c_a, c_b by requiring that the PDF
                // a) be continuous at x
                // b) integrate to one over the range
                const double _c_a, _c_b;

                // The probability covered to the left of the central value. It is just
                // c_b (1/2 - \Phi(y_-|x,b)
                // but for the sampling it is faster to precompute the call to Gaussian cumulative.
                const double _prob_lower;

                // PDF normalization factor precomputed for operator()
                const double _norm_lower, _norm_upper;

            public:
                CurtailedGauss(const Parameters & parameters, const std::string & name, const ParameterRange & range,
                        const double & lower, const double & central, const double & upper) :
                    LogPrior(parameters),
                    _parameter(parameters[name]),
                    _name(name),
                    _range(range),
                    _lower(lower),
                    _central(central),
                    _upper(upper),
                    _sigma_lower(central - lower),
                    _sigma_upper(upper - central),
                    _c_a(1.0 / ((_sigma_lower / _sigma_upper) * (0.5 - gsl_cdf_gaussian_P(range.min - central, _sigma_lower))
                                 + gsl_cdf_gaussian_P(range.max - central, _sigma_upper) - 0.5)),
                    _c_b(_sigma_lower/ _sigma_upper * _c_a),
                    _prob_lower(_c_b * (0.5 - gsl_cdf_gaussian_P(range.min - central, _sigma_lower))),
                    _norm_lower(std::log(_c_b / std::sqrt(2 * M_PI) / _sigma_lower)),
                    _norm_upper(std::log(_c_a / std::sqrt(2 * M_PI) / _sigma_upper))
                {
                    if (range.min >= range.max)
                    {
                        throw RangeError("LogPrior::Gauss(" + _name +"): minimum (" + stringify(range.min)
                                          + ") must be smaller than maximum (" + stringify(range.max) + ")");
                    }
                    _parameter_descriptions.push_back(ParameterDescription{ _parameters[name].clone(), range.min, range.max, false });
                }

                virtual ~CurtailedGauss()
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
                    double x = _parameter_descriptions.front().parameter->evaluate();

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
                    return norm - 0.5 * power_of<2>((x - _central) / sigma);
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    return LogPriorPtr(new priors::CurtailedGauss(parameters, _name, _range, _lower, _central, _upper));
                }

                virtual void sample()
                {
                    // CDF = c \Phi(x - x_{central} / \sigma) + b

                    // find out if sample in upper or lower part
                    const auto p = _parameter.evaluate_generator();

                    if (p < _prob_lower)
                       _parameter.set(gsl_cdf_gaussian_Pinv((p - _prob_lower) / _c_b + 0.5,  _sigma_lower) + _central);
                    else
                       _parameter.set(gsl_cdf_gaussian_Pinv((p - _prob_lower) / _c_a + 0.5,  _sigma_upper) + _central);
                }

                virtual bool informative() const
                {
                    return true;
                }
        };

        /*!
         * Prior distribution for renormalization scales
         */
        class Scale :
            public LogPrior
        {
            private:
                Parameter _parameter;

                const std::string _name;

                const ParameterRange _range;

                const double _mu_0, _lambda;

                const double _min, _max;

                const double _ln_lambda;

            public:
                Scale(const Parameters & parameters, const std::string & name, const ParameterRange & range, const double & mu_0, const double & lambda) :
                    LogPrior(parameters),
                    _parameter(parameters[name]),
                    _name(name),
                    _range(range),
                    _mu_0(mu_0),
                    _lambda(lambda),
                    _min(mu_0 / lambda),
                    _max(mu_0 * lambda),
                    _ln_lambda(std::log(lambda))
                {
                    _parameter_descriptions.push_back(ParameterDescription{ _parameters[name].clone(), mu_0 / lambda, mu_0 * lambda, false });
                }

                virtual ~Scale()
                {
                }

                virtual std::string as_string() const
                {
                    std::string result = "Parameter: " + _name + ", prior type: Scale, range: [" + stringify(_mu_0 / _lambda) + "," + stringify(_mu_0 * _lambda) + "]";
                    result += ", mu_0 = " + stringify(_mu_0) + ", lambda = " + stringify(_lambda);

                    return result;
                }

                virtual double operator()() const
                {
                    // read parameter's current value
                    double x = _parameter_descriptions.front().parameter->evaluate();

                    if ((x < _min) || (_max < x))
                        return -std::numeric_limits<double>::infinity();

                    return 1.0 / (2.0 * _ln_lambda * x);
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    return LogPriorPtr(new priors::Scale(parameters, _name, _range, _mu_0, _lambda));
                }

                virtual void sample()
                {
                    // CDF: p = [\ln x - \ln \mu_0 + \ln \lambda] / (2.0 \ln \lambda)
                    // inverse CDF: x = \mu_0 * \lambda^(2 p - 1)

                    _parameter.set(_mu_0 * std::pow(_lambda, 2.0 * _parameter.evaluate_generator() - 1.0));
                }

                virtual bool informative() const
                {
                    return true;
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
    LogPrior::Flat(const Parameters & parameters, const std::string & name, const ParameterRange & range)
    {
        LogPriorPtr prior = std::make_shared<eos::priors::Flat>(parameters, name, range);

        return prior;
    }

    LogPriorPtr
    LogPrior::CurtailedGauss(const Parameters & parameters, const std::string & name, const ParameterRange & range,
            const double & lower, const double & central, const double & upper)
    {
        // check input
        if (lower >= central)
            throw InternalError("LogPrior::Gauss: lower value (" + stringify(lower) + ") >= central value (" + stringify(central) + ")");

        if (upper <= central)
            throw InternalError("LogPrior::Gauss: upper value (" + stringify(upper) + ") <= central value (" + stringify(central) + ")");

        LogPriorPtr prior = std::make_shared<eos::priors::CurtailedGauss>(parameters, name, range, lower, central, upper);

        return prior;
    }

    LogPriorPtr
    LogPrior::Scale(const Parameters & parameters, const std::string & name, const ParameterRange & range,
            const double & mu_0, const double & lambda)
    {
        // check input
        if (mu_0 <= 0.0)
            throw InternalError("LogPrior::Scale: default value mu_0 must be strictly positive");

        if (lambda <= 1.0)
            throw InternalError("LogPrior::Scale: scale factor lambda must be strictly larger than 1");

        LogPriorPtr prior = std::make_shared<eos::priors::Scale>(parameters, name, range, mu_0, lambda);

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
        if (prior_type == "Gaussian")
        {
            // extract central value
            loc1 = s.find_first_of('=', loc2 + 1);
            loc2 = s.find_first_of('+', loc2 + 1);

            double central = destringify<double>(s.substr(loc1 + 2, loc2 - loc1 - 3));

            double sigma_upper, sigma_lower;

            // extract sigma_upper, lower
            if ( s[loc2 + 1] == '-')
            {
                sigma_upper = destringify<double>(s.substr(loc2 + 3));
                sigma_lower = sigma_upper;
            }
            else
            {
                loc1 = loc2;
                loc2 = s.find_first_of('-', loc2 + 1);

                sigma_upper = destringify<double>(s.substr(loc1 + 2, loc2 - loc1 - 3));

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
                    sigma_lower = destringify<double>(s.substr(loc1 + 2, loc2 - loc1 - 2));
                }
            }

            if (prior_type == "Gaussian")
                return CurtailedGauss(parameters, par_name, range, central - sigma_lower, central, central + sigma_upper);
        }

        throw priors::UnknownPriorError("Cannot construct prior from '" + s + "'");
    }

    template class WrappedForwardIterator<LogPrior::IteratorTag, ParameterDescription>;
}
