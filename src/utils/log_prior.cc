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

#include <src/utils/log_prior.hh>
#include <src/utils/power_of.hh>
#include <src/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <gsl/gsl_cdf.h>

namespace eos
{
    template class WrappedForwardIterator<LogPrior::IteratorTag, ParameterDescription>;

    namespace priors
    {
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
                    _parameter_descriptions.push_back(ParameterDescription{ _parameters[name], range.min, range.max, false, false });
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
                // need to fix the coefficients in
                // @f$CDF_{lower}(x) = c_{lower} * \right( \Phi(x, \sigma_{lower}) - \Phi(x_{min},\sigma_{lower}) \left)$@f
                // @f$CDF_{upper}(x) = c_{upper} * \right( \Phi(x, \sigma_{upper}) + P_{lower}/c_{upper} - 1/2\left)$@f.
                double _c_lower, _c_upper;
                double _phi_min, _phi_max;

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
                    _parameter_descriptions.push_back(ParameterDescription{ _parameters[name], range.min, range.max, false, false });

                    /*
                     * calculate anything that doesn't depend on parameter x only once.
                     * the normalization constant is found from
                     * $1=c\cdot\int_{x_{min}}^{x_{max}}p_{0}\left(x\right)=
                     *    c\int_{x_{min}}^{\mu}f\left(x|\mu,\sigma_{lower}\right)+c\int_{\mu}^{x_{max}}f\left(x|\mu,\sigma_{upper}\right)$
                     * where f denotes the Normal distribution
                     * $f\left(x|\mu,\sigma\right)=\frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{(x-\mu)^{2}}{2\sigma^{2}}}$
                     */

                    _phi_min = gsl_cdf_gaussian_P(range.min - _central, _sigma_lower);
                    _phi_max = gsl_cdf_gaussian_P(range.max - _central, _sigma_upper);

                    double overall_normalization = 1.0 / (_phi_max - _phi_min);

                    _prob_lower = (0.5 - _phi_min) * overall_normalization;

                    _norm_lower = std::log(overall_normalization / (std::sqrt(2.0 * M_PI) * _sigma_lower));
                    _norm_upper = std::log(overall_normalization / (std::sqrt(2.0 * M_PI) * _sigma_upper));

                    _c_lower = 2.0 * _prob_lower / (1.0 - _phi_min);
                    _c_upper = (1.0 - _prob_lower) / (_phi_max - 0.5);
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
                    if (u < _prob_lower)
                        return gsl_cdf_gaussian_Pinv(u / _c_lower + _phi_min, _sigma_lower) + _central;
                    else
                        return gsl_cdf_gaussian_Pinv((u - _prob_lower) / _c_upper + 0.5, _sigma_upper) + _central;
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
        LogPriorPtr prior = std::make_shared<eos::priors::Gauss>(parameters, name, range, lower, central, upper);

        return prior;
    }
}
