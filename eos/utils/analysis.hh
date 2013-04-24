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

#ifndef EOS_GUARD_SRC_UTILS_ANALYSIS_HH
#define EOS_GUARD_SRC_UTILS_ANALYSIS_HH 1

#include <eos/utils/analysis-fwd.hh>
#include <eos/utils/hdf5-fwd.hh>
#include <eos/utils/log_likelihood.hh>
#include <eos/utils/log_prior.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/verify.hh>

#include <vector>
#include <set>

#include <gsl/gsl_multimin.h>

namespace ROOT
{
    namespace Minuit2
    {
        class FunctionMinimum;
    }
}

namespace eos
{
    class Analysis :
        public PrivateImplementationPattern<Analysis>
    {
        public:
            friend struct Implementation<Analysis>;

            struct OptimizationOptions;
            struct Output;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * Extracts parameters, observables from LogLikelihood.
             * The default prior (flat) is assumed for all parameters.
             *
             * @param log_likelihood  The LogLikelihood functor which shall be analysed.
             *
             * @note Analysis assumes ownership of log_likelihood
             */
            Analysis(const LogLikelihood & log_likelihood);

            /// Destructor.
            ~Analysis();

            /// Clone this Analysis
            AnalysisPtr clone() const;
            ///@}

            ///@name Accessors
            ///@{

            /// Retrieve a set of all parameters, including ranges
            const std::vector<ParameterDescription> & parameter_descriptions() const;

            /*!
             * Retrieve a parameter by index.
             *
             * @param index The index of the parameter.
             */
            Parameter operator[] (const unsigned & index);

            /// Retrieve our associates Parameters object
            Parameters parameters() const;

            /*!
             * Add one or more parameters and associated prior density
             *
             * @param prior    The logarithmic prior density.
             * @param nuisance False for a parameter of interest
             */
            bool add(const LogPriorPtr & prior, bool nuisance = false);

            /*!
             * Write parameter descriptions, constraints, observables
             * into the hdf5 file under the given group name.
             */
            void dump_descriptions(hdf5::File & file, std::string data_set_base = "/descriptions") const;

            /*!
             * Read in parameter descriptions from a previous dump.
             *
             * @param sample_file The HDF5 file with the information.
             * @param base The base directory where to look for the data set.
             * @return The descriptions, one per parameter
             */
            static std::vector<ParameterDescription> read_descriptions(const hdf5::File & file, std::string data_set_base = "/descriptions");

            /*!
             * Calculate the p-value based on the @f$\chi^2 @f$
             * test statistic for fixed parameter_values
             * @param parameter_values
             * @param simulate if true, simulate data sets to estimate distribution of test statistic,
             * else use the cumulative of the @f$\chi^2 @f$-distribution with (N-k) degrees-of-freedom,
             * where N is the number of observations and k is the number of fitted parameters
             * @return < @f$\chi^2 @f$, p>
             */
            std::pair<double, double>
            goodness_of_fit(const std::vector<double> & parameter_values, const unsigned & simulated_datasets, std::string output_file = "");

            /// Retrieve the overall Log(likelihood) for this analysis.
            LogLikelihood log_likelihood();

            /// Retrieve the overall Log(prior) for this analysis.
            double log_prior();

            /*!
             * Find the prior for a given parameter
             */
            LogPriorPtr log_prior(const std::string & name) const;

            /// Retrieve the overall Log(posterior) for this analysis.
            /// Incorporate normalization constant, the evidence here in getter if available.
            double log_posterior();

            /*!
             * Check if a given parameter is a nuisance parameter for this Analysis.
             *
             * @param name The name of the parameter we are interested in.
             */
            bool nuisance(const std::string & name) const;
            ///@}

            /*!
             * Optimize the posterior using the Nelder-Mead simplex algorithm.
             * @param initial_guess Starting point for simplex construction
             * @param options If no tuning desired, use Analysis::OptimizationOptions::Defaults()
             * @return <parameter values at mode, posterior value at mode>
             */
            std::pair<std::vector<double>, double>
            optimize(const std::vector<double> & initial_guess, const OptimizationOptions & options);

            const ROOT::Minuit2::FunctionMinimum &
            optimize_minuit(const std::vector<double> & initial_guess, const OptimizationOptions & options);

            /*!
             * Restrict to a subrange of a given parameter.
             * @param name
             * @param min
             * @param max
             */
            void restrict(const std::string & name, const double & min, const double & max);
    };

        struct Analysis::OptimizationOptions
        {
                /// Options are: "migrad", "minimize", "scan", "simplex" from minuit2
                std::string algorithm;

                /// Keep the value of nuisance parameters with a flat prior fixed at the current value during optimization,
                /// to avoid flat directions that cause Migrad to fail.
                bool fix_flat_nuisance;

                /// Fraction of parameter range, in [0,1].
                /// Useful only for simplex method
                VerifiedRange<double> initial_step_size;

                /// If algorithm doesn't converge before, quit
                /// after maximum_iterations.
                unsigned maximum_iterations;

                /*!
                 * If non-zero, perform MCMC iterations first,
                 * before Minuit2 is invoked from the last point of the chain.
                 *
                 * @note This is only useful when call from MarkovChainSampler,
                 *       further control of the chain is taken from the MarkovChainSampler::Config object.
                 */
                bool mcmc_pre_run;

                /*!
                 *  Once the algorithm has shrunk the probe
                 *  simplex below this size, convergence is declared.
                 *
                 *  For minuit, it is just their tolerance parameter
                 */
                VerifiedRange<double> tolerance;

                /*!
                 * When comparing two modes found by minuit to decide whether they correspond
                 * to the same mode, the splitting_tolerance decides how far
                 * in relative units their distance may be.
                 */
                VerifiedRange<double> splitting_tolerance;

                /// 0 - low, 1 - medium, 2 - high precision
                VerifiedRange<unsigned> strategy_level;

                static OptimizationOptions Defaults();

            private:
                OptimizationOptions();
        };

     struct Analysis::Output
     {
         typedef hdf5::Composite<hdf5::Scalar<const char *>, hdf5::Scalar<double>, hdf5::Scalar<double>,
                                 hdf5::Scalar<int>, hdf5::Scalar<const char *>> DescriptionType;
         static DescriptionType description_type();
         static std::tuple<const char *, double, double, int, const char *> description_record();
     };
}

#endif
