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

#include <src/utils/log.hh>
#include <src/utils/log_likelihood.hh>
#include <src/utils/log_prior.hh>
#include <src/utils/observable_cache.hh>
#include <src/utils/power_of.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include <gsl/gsl_cdf.h>

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

            GaussianBlock(const ObservableCache & cache, ObservableCache::Id id, const double & min, const double & central, const double & max) :
                cache(cache),
                id(id),
                central(central),
                sigma_lower(central - min),
                sigma_upper(max - central)
            {
            }

            ~GaussianBlock()
            {
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

            virtual LogLikelihoodBlockPtr clone(ObservableCache cache) const
            {
                ObservablePtr observable = this->cache.observable(id)->clone(cache.parameters());

                return LogLikelihoodBlockPtr(new GaussianBlock(cache, cache.add(observable), central - sigma_lower, central, central + sigma_upper));
            }
        };
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
            const double & min, const double & central, const double & max)
    {
        // check input
        if (min >= central)
            throw InternalError("LogLikelihoodBlock::Gaussian: min value >= central value");

        if (max <= central)
            throw InternalError("LogLikelihoodBlock::Gaussian: max value <= central value");

        unsigned index = cache.add(observable);

        return LogLikelihoodBlockPtr(new implementation::GaussianBlock(cache, index, min, central, max));
    }

    template <>
    struct Implementation<LogLikelihood>
    {
        Parameters parameters;

        // Cache observable predictions
        ObservableCache cache;

        // Container for all existing independent (i.e. uncorrelated) block
        std::vector<LogLikelihoodBlockPtr> blocks;

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
            for (auto b = blocks.cbegin(), b_end = blocks.cend() ; b != b_end ; ++b)
            {
                (*b)->prepare_sampling();
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
                for (auto b = blocks.cbegin(), b_end = blocks.cend() ; b != b_end ; ++b)
                {
                    t += (*b)->sample(rng);
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

            for (auto b = blocks.cbegin(), b_end = blocks.cend() ; b != b_end ; ++b)
            {
                result += (*b)->evaluate();
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
    LogLikelihood::add(const ObservablePtr & observable, const double & min, const double & central, const double & max)
    {
        _imp->blocks.push_back(LogLikelihoodBlock::Gaussian(_imp->cache, observable, min, central, max));
    }

    void
    LogLikelihood::add(const Constraint & constraint)
    {
        for (auto b = constraint.begin_blocks(), b_end = constraint.end_blocks() ; b != b_end ; ++b)
        {
            // Clone each LogLikelihoodBlock onto our ObservableCache
            _imp->blocks.push_back((*b)->clone(_imp->cache));
        }
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

        for (auto b = _imp->blocks.cbegin(), b_end = _imp->blocks.cend() ; b != b_end ; ++b)
        {
            result._imp->blocks.push_back((*b)->clone(result._imp->cache));
        }

        return result;
    }

    unsigned
    LogLikelihood::number_of_observations() const
    {
        return _imp->blocks.size();
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
