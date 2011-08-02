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

namespace eos
{
    template <>
    struct Implementation<LogLikelihood>
    {
        Parameters parameters;

        // Cache observable predictions
        ObservableCache cache;

        // <index into cache, min, central, max>
        std::vector<std::tuple<ObservableCache::Id, double, double, double>> observations;

        double chi_squared;

        Implementation(const Parameters & parameters) :
            parameters(parameters),
            cache(parameters),
            chi_squared(0)
        {
        }

        void add(const ObservablePtr & observable, const double & min, const double & central, const double & max)
        {
            auto result = cache.add(observable);
            observations.push_back(std::make_tuple(result, min, central, max));
        }

        std::pair<double, double>
        bootstrap_p_value(const unsigned & datasets)
        {
            // create a Gaussian prior for each observation
            // in order to sample possible observations
            // take exp. uncertainty fixed, but sample around predictions
            // for the (best fit) parameters
            // TODO can we make sure that no unphysical observations are generated. Is three sigma safe?
            // TODO asymmetry not correctly captured, as the discontinuity in the pdf shouldn't be at central,
            // but rather at x_obs
            Parameters pars = parameters.clone();
            std::vector<LogPriorPtr> generators;

            // store predictions to compare with, duplicate if necessary to avoid index lookup
            // in data generation loop
            std::vector<double> preds;
            for (auto i = observations.cbegin(), i_end = observations.cend() ; i != i_end ; ++i)
            {
                double central     = cache[std::get<0>(*i)]; //fixed parameters => fixed predictions
                double sigma_lower = std::get<2>(*i) - std::get<1>(*i);
                double sigma_upper = std::get<3>(*i) - std::get<2>(*i);
                generators.push_back( LogPrior::Gauss(pars, "mass::b(MSbar)",
                                ParameterRange{ central - 3 * sigma_lower, central + 3 * sigma_upper },
                                central - sigma_lower, central, central + sigma_upper));
                preds.push_back(cache[std::get<0>(*i)]);
            }

            // save the current chi^2 for the observed data points
            double chi_squared_obs = chi_squared;

            Log::instance()->message("log_likelihood.bootstrap_pvalue", ll_informational)
                                                 << "chi^2_obs = " << chi_squared_obs;

            gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, 1536);

            // keep track how often chi2 is bigger in pseudo data
            unsigned counter = 0;

            for (unsigned i=0; i < datasets; ++i)
            {
                double pseudo_chi2 = 0;
                double pseudo_obs = 0;
                double sigma = 0;

                // generate one pseudo data set
                unsigned i = 0;
                for (auto g = generators.begin(), g_end = generators.end() ; g != g_end ; ++g, ++i)
                {
                    pseudo_obs = (**g).sample(rng);
                    if (pseudo_obs > std::get<2>(observations[i]) )
                        sigma = std::get<3>(observations[i]) - std::get<2>(observations[i]);
                    else
                        sigma = std::get<2>(observations[i]) - std::get<1>(observations[i]);
                    pseudo_chi2 += power_of<2>( (pseudo_obs - preds[i]) / sigma);
                }

                // compare with predictions
                if (pseudo_chi2 > chi_squared_obs)
                    counter++;
            }

            gsl_rng_free(rng);

            // naive estimate: mode of binomial posterior
            double p = double(counter)/double(datasets);

            double p_expected = double(counter+1)/double(datasets+2);
            double uncertainty = sqrt(p_expected * ( 1 - p_expected) / double(datasets + 3));

            return std::make_pair(p, uncertainty);
        }

        double operator() ()
        {
           /*
            * proper normalization
            * Could calculate only once upon initialization, but if
            * uncertainty is parameter dependent, that's ruled out.
            */
           double norm = 1.0;

           chi_squared = 0.0;
           double sigma = 1.0;

           // compare predictions with observations
           for (auto i = observations.cbegin(), i_end = observations.cend() ; i != i_end ; ++i)
           {
               // get prediction
               double value = cache[std::get<0>(*i)];

               // get most likely observed value
               double central = std::get<2>(*i);

               // allow for asymmetric Gaussian uncertainty
               if (value > central)
                   sigma = std::get<3>(*i) - central;
               else
                   sigma = central - std::get<1>(*i);

               double chi = (value - central) / sigma;

               chi_squared += chi * chi;

               norm /= std::sqrt(2.0 * M_PI) * sigma;
           }

           return std::log(norm) - chi_squared / 2.0;
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
        _imp->add(observable, min, central, max);
    }

    std::pair<double, double>
    LogLikelihood::bootstrap_p_value(const unsigned & datasets)
    {
        return _imp->bootstrap_p_value(datasets);
    }

    double
    LogLikelihood::chi_squared() const
    {
        return _imp->chi_squared;
    }

    LogLikelihood
    LogLikelihood::clone() const
    {
        Parameters parameters(_imp->parameters.clone());
        LogLikelihood result(parameters);

        for (auto o = _imp->observations.cbegin(), o_end = _imp->observations.cend() ; o != o_end ; ++o)
        {
            // create clone with independent parameters object
            ObservablePtr ptr = _imp->cache.observable(std::get<0>(*o))->clone(parameters);
            result.add(ptr, std::get<1>(*o), std::get<2>(*o), std::get<3>(*o));
        }

        return result;
    }

    unsigned
    LogLikelihood::number_of_observations() const
    {
        return _imp->observations.size();
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

        return (*_imp)();
    }

    double
    LogLikelihood::operator() (const Parameter::Id & id)
    {
        _imp->cache.update(id);

        return (*_imp)();
    }

    void
    LogLikelihood::reset()
    {
        _imp->cache.reset();
    }
}
