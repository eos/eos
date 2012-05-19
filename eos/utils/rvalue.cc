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

#include <eos/utils/exception.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/rvalue.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/log.hh>

#include <cmath>
#include <limits>

namespace eos
{
    double
    RValue::gelman_rubin(const std::vector<double> & chain_means, const std::vector<double> & chain_variances,
            const unsigned & chain_length)
    {
        if (chain_means.size() != chain_variances.size())
            throw InternalError("RValue::gelman_rubin: chain means and chain variances are not aligned!");

        double n = chain_length;
        double m = chain_means.size();

        if (m <= 1)
            throw InternalError("RValue::gelman_rubin: Need at least two chains to compute R-value!");

        // init
        double variance_of_means = 0.0;
        double mean_of_means = 0.0;
        double mean_of_variances = 0.0;
        double variance_of_variances = 0.0;

        // use Welford's method here as well with temporary variables
        double previous_mean_of_means(0);
        double previous_mean_of_variances(0);
        for (unsigned i = 0 ; i < m ; ++i)
        {
            if (0 == i)
            {
                mean_of_means = chain_means.front();
                variance_of_means = 0;
                mean_of_variances = chain_variances.front();
                variance_of_variances = 0;

                continue;
            }

            // temporarily store previous mean of step (i-1)
            previous_mean_of_means = mean_of_means;
            previous_mean_of_variances = mean_of_variances;

            // update step
            mean_of_means += (chain_means[i] - previous_mean_of_means) / (i + 1.0);
            variance_of_means += (chain_means[i] - previous_mean_of_means) * (chain_means[i] - mean_of_means);
            mean_of_variances += (chain_variances[i] - previous_mean_of_variances) / (i + 1.0);
            variance_of_variances += (chain_variances[i] - previous_mean_of_variances) * (chain_variances[i] - mean_of_variances);
        }

        variance_of_means /= m - 1.0;
        variance_of_variances /= m - 1.0;

        //use Gelman/Rubin notation
        double B = variance_of_means * n;
        double W = mean_of_variances;
        double sigma_squared = (n - 1.0) / n * W + B / n;

        // avoid NaN due to divide by zero
        if (0.0 == W)
        {
            Log::instance()->message("Rvalue.gelman_rubin", ll_debug) << "W = 0. Avoiding R = NaN.";
            return std::numeric_limits<double>::max();
        }

        //estimated scale reduction
        double R = 0.0;

        // compute covariances using the means from above
        double covariance_22 = 0.0; // cov(s_i^2, \bar{x_i}^2
        double covariance_21 = 0.0; // cov(s_i^2, \bar{x_i}

        for (unsigned i = 0 ; i < m ; ++i)
        {
            covariance_21 += (chain_variances[i] - mean_of_variances) * (chain_means.at(i) - mean_of_means);
            covariance_22 += (chain_variances[i] - mean_of_variances) * (power_of<2>(chain_means[i]) - power_of<2>(mean_of_means));
        }

        covariance_21 /= m - 1.0;
        covariance_22 /= m - 1.0;

        // scale of t-distribution
        double V = sigma_squared + B / (m * n);

        // estimation of scale variance
        double a = (n - 1.0) * (n - 1.0) / (n * n * m) * variance_of_variances;
        double b = (m + 1) * (m + 1) / (m * n * m * n) * 2.0 / (m - 1) * B * B;
        double c = 2.0 * (m + 1.0) * (n - 1.0) / (m * n * n) * n / m * (covariance_22 - 2.0 * mean_of_means * covariance_21);
        double variance_of_V = a + b + c;

        // degrees of freedom of t-distribution
        double df = 2.0 * V * V / variance_of_V;

        if (df <= 2)
        {
            Log::instance()->message("Rvalue.gelman_rubin", ll_debug)
                << "DoF (" << df << ") below 2. Avoiding R = NaN.";
            return std::numeric_limits<double>::max();;
        }

        // sqrt of estimated scale reduction if sampling were continued
        R = std::sqrt(V / W * df / (df - 2.0));

        // R smaller, but close to 1 is OK.
        if ((R < 0.99) && (n > 100))
            throw InternalError("MarkovChainSampler::compute_rvalue: R-value " + stringify(R, 4) + " < 0.99. Check for a bug in the implementation!");

        return R;
    }

    double
    RValue::approximation(const std::vector<double> & chain_means, const std::vector<double> & chain_variances,
            const unsigned & chain_length)
    {
        if (chain_means.size() != chain_variances.size())
            throw InternalError("RValue::gelman_rubin: chain means and chain variances are not aligned!");

        double n = chain_length;
        double m = chain_means.size();

        if (m <= 1)
            throw InternalError("RValue::gelman_rubin: Need at least two chains to compute R-value!");

        // init
        double variance_of_means = 0.0;
        double mean_of_means = 0.0;
        double mean_of_variances = 0.0;
        double variance_of_variances = 0.0;

        // use Welford's method here as well with temporary variables
        double previous_mean_of_means(0);
        double previous_mean_of_variances(0);
        for (unsigned i = 0 ; i < m ; ++i)
        {
            if (0 == i)
            {
                mean_of_means = chain_means.front();
                variance_of_means = 0;
                mean_of_variances = chain_variances.front();
                variance_of_variances = 0;

                continue;
            }

            // temporarily store previous mean of step (i-1)
            previous_mean_of_means = mean_of_means;
            previous_mean_of_variances = mean_of_variances;

            // update step
            mean_of_means += (chain_means[i] - previous_mean_of_means) / (i + 1.0);
            variance_of_means += (chain_means[i] - previous_mean_of_means) * (chain_means[i] - mean_of_means);
            mean_of_variances += (chain_variances[i] - previous_mean_of_variances) / (i + 1.0);
            variance_of_variances += (chain_variances[i] - previous_mean_of_variances) * (chain_variances[i] - mean_of_variances);
        }

        variance_of_means /= m - 1.0;
        variance_of_variances /= m - 1.0;

        //use Gelman/Rubin notation
        double B = variance_of_means * n;
        double W = mean_of_variances;
        double sigma_squared = (n - 1.0) / n * W + B / n;

        // avoid NaN due to divide by zero
        if (0.0 == W)
        {
            return std::numeric_limits<double>::max();
        }

        //estimated scale reduction
        double R = std::sqrt(sigma_squared / W);

        return R;
    }
}
