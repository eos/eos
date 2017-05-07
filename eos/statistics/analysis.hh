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

#ifndef EOS_GUARD_SRC_STATISTICS_ANALYSIS_HH
#define EOS_GUARD_SRC_STATISTICS_ANALYSIS_HH 1

#include <eos/statistics/analysis-fwd.hh>
#include <eos/statistics/log-likelihood.hh>
#include <eos/statistics/log-prior.hh>
#include <eos/utils/density.hh>
#include <eos/utils/hdf5-fwd.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/verify.hh>

#include <set>
#include <vector>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

namespace ROOT
{
    namespace Minuit2
    {
        class FunctionMinimum;
    }
}

namespace eos
{
    struct MinuitAdapter;

    class Analysis :
        public Density
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
            virtual ~Analysis();

            /// Clone this Analysis
            // todo remove
            AnalysisPtr old_clone() const;

            virtual DensityPtr clone() const;

            virtual double evaluate() const;

            virtual Iterator begin() const;
            virtual Iterator end() const;
            ///@}

            ///@name Accessors
            ///@{
            // todo remove
            /// Retrieve a set of all parameters, including ranges
            const std::vector<ParameterDescription> & parameter_descriptions() const;

            // todo remove as functionality available from begin(), end()
            /*!
             * Retrieve a parameter by index.
             *
             * @param index The index of the parameter.
             */
            MutablePtr operator[] (const unsigned & index) const;

            // todo remove
            /// Retrieve our associated Parameters object
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
            virtual void dump_descriptions(hdf5::File & file, const std::string & data_set_base) const;

            /*!
             * Read in parameter descriptions from a previous dump.
             *
             * @param sample_file The HDF5 file with the information.
             * @param base The base directory where to look for the data set.
             * @return The descriptions, one per parameter
             */
            static std::vector<ParameterDescription> read_descriptions(const hdf5::File & file, std::string data_set_base = "/descriptions");

            /*!
<<<<<<< HEAD:eos/statistics/analysis.hh
             * Read the description part of chain's prerun from hdf5 file.
             *
             * @param file
             * @param data_base_name The directory in the file under which the data is parsed.
             * @param descriptions All parameter ranges etc. Beware, the association to the underlying Parameters object is independent.
             * @param priors The string representation of a prior distribution.
             * @param constraints The string representation of an individual constraint.
             * @param hash The EOS version used to create the file.
             */
            static void read_descriptions(hdf5::File & file, const std::string & data_base_name,
                                          std::vector<ParameterDescription>& descriptions,
                                          std::vector<std::string> & priors,
                                          std::vector<std::string> & constraints,
                                          std::string & hash);

            /*!
             * Calculate the p-value based on the @f$\chi^2 @f$
             * test statistic for fixed parameter_values
=======
             * Calculate two p values based on the @f$\chi^2 @f$
             * test statistic for fixed parameter_values.
             *
             * The first is based on pulls, or significances,
             * defined in units of Gaussian standard deviations. Their squared sum
             * is a @f$\chi^2 @f$. The second uses the log-likelihood as a test statistic,
             * and empirically generates data from the likelihood blocks to simulate the statistic's
             * distribution and a p value. With N observations, this p value is transformed into a @f$\chi^2 @f$
             * using the inverse cumulative of the @f$\chi^2 @f$ distribution with N degrees of freedom.
             *
             * Both @f$\chi^2 @f$ values are transformed into a p value through
             * cumulative of the @f$\chi^2 @f$-distribution with (N-k) degrees-of-freedom,
             * where N is the number of observations and k is the number of fitted parameters.
             *
>>>>>>> 40893a8... [utils] Improve Analysis.goodness_of_fit documentation and error messages:eos/utils/analysis.hh
             * @param parameter_values
             * @param simulated_datasets The number of data sets used to estimate the distribution of log-likelihood test statistic.
             * @param output_file If given, store pulls and @f$\chi^2 @f$ in an HDF5 file.
             * @return < @f$\chi^2 @f$, p>
             *
             * @note Nuisance parameters are assumed to have an informative prior counted as one @e observation.
             *       They therefore cancel in computing the degrees of freedom.
             */
            std::pair<double, double>
            goodness_of_fit(const std::vector<double> & parameter_values, const unsigned & simulated_datasets, std::string output_file = "");

            /// Retrieve the overall Log(likelihood) for this analysis.
            LogLikelihood log_likelihood() const;

            /// Retrieve the overall Log(prior) for this analysis.
            double log_prior() const;

            /*!
             * Find the prior for a given parameter
             */
            LogPriorPtr log_prior(const std::string & name) const;

            /// Retrieve the overall Log(posterior) for this analysis.
            /// Incorporate normalization constant, the evidence here in getter if available.
            double log_posterior() const;

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

        private:
            /*!
             * Find index of definition of parameter
             * @param name
             * @return index if found, _parameter_descriptions.size() if not found
             */
            unsigned index(const std::string & name) const;

            /*!
             * routine needed for optimization. It returns the negative
             * log(posterior) at parameters
             * @param parameters the values of the parameters
             * @param data pointer to an Analysis object
             * @return
             */
            static double
            negative_log_posterior(const gsl_vector * pars, void * data);

            Analysis *
            private_clone() const;

            LogLikelihood _log_likelihood;

            Parameters _parameters;

            /// prior in N dimensions can decouple
            /// at most into N 1D priors
            std::vector<LogPriorPtr> _priors;

            /// Parameter, minimum, maximum, nuisance
            std::vector<ParameterDescription> _parameter_descriptions;

            /// names of all parameters. prevent using a parameter twice
            std::set<std::string> _parameter_names;

            /// Adapter to let minuit operate on posterior
            MinuitAdapter * _minuit;
    };

        // todo move optimization into separate class
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

     /*!
      * Compute an initial guess of the proposal covariance matrix.
      *
      * The variance of each parameter is taken from the prior distribution
      * and scaled if desired for higher efficiency. Zero correlation is assumed a priori.
      *
      * @param analysis The analysis supplying the prior information.
      * @param scale_reduction Value by which sqrt(variance) of parameters is divided.
      * @param scale_nuisance Decide whether only scan parameters or all parameters are scaled.
      * @return The covariance matrix in row major format.
      */
     std::vector<double> proposal_covariance(const Analysis & analysis,
                                             double scale_reduction=1,
                                             bool scale_nuisance=true);
}

#endif
