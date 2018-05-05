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

#include <eos/statistics/log-posterior_TEST.hh>
#include <eos/statistics/density-wrapper_TEST.hh>
#include <eos/statistics/population-monte-carlo-sampler.hh>
#include <eos/statistics/markov-chain-sampler.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <test/test.hh>

using namespace test;
using namespace eos;

class PopulationMonteCarloSamplerTest :
    public TestCase
{
    public:
        PopulationMonteCarloSamplerTest() :
            TestCase("population_monte_carlo_sampler_test")
        {
        }

        void wrapped_density() const
        {
            /* create prerun output */
            DensityWrapper density = make_multivariate_unit_normal(2);
            MarkovChainSampler::Config mcmc_config = MarkovChainSampler::Config::Default();
            mcmc_config.need_main_run = false;
            mcmc_config.number_of_chains = 2;
            mcmc_config.output_file = EOS_BUILDDIR "/eos/statistics/pmc_sampler_TEST-density-prerun.hdf5";
            mcmc_config.parallelize = false;
            mcmc_config.prerun_iterations_update = 300;
            mcmc_config.prerun_iterations_max = 5000;
            mcmc_config.prerun_iterations_min = 500;
            mcmc_config.seed = 1246122;
            MarkovChainSampler sampler(density.clone(), mcmc_config);
            sampler.run();

            /* run PMC */
            PopulationMonteCarloSampler::Config pmc_config = PopulationMonteCarloSampler::Config::Default();
            pmc_config.max_updates = 5;
            pmc_config.samples_per_component = 400;
            pmc_config.final_samples = 50000;
            pmc_config.output_file = EOS_BUILDDIR "/eos/statistics/pmc_sampler_TEST-density.hdf5";;
            pmc_config.parallelize = true;
            pmc_config.seed = 23;
            pmc_config.store = true;
            pmc_config.skip_initial = 0.2;
            pmc_config.patch_length = 100;
            pmc_config.target_ncomponents = 2;
            pmc_config.group_by_r_value = 1.2;
            PopulationMonteCarloSampler pmc_sampler(density.clone(), hdf5::File::Open(mcmc_config.output_file), pmc_config);
            pmc_sampler.run();
            TEST_CHECK(pmc_sampler.status().converged);
}

        virtual void run() const
        {
            wrapped_density();
            // initialize from a MCMC prerun
            {
                /* setup bimodal distribution */

                Parameters p = Parameters::Defaults();
                Kinematics k;
                std::array<ObservablePtr, 2> obs;
                obs[0] = ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k));
                obs[1] = ObservablePtr(new AbsoluteTestObservable(p, k, "mass::c"));

                ObservableCache cache(p);

                // multivariate gaussian with no correlation
                std::array<double, 2> mean{{ +5, +5 }};
                std::array<std::array<double, 2>, 2> covariance;
                covariance[0][0] = 0.1 * 0.1;
                covariance[1][1] = 0.05 * 0.05;
                covariance[0][1] = covariance[1][0] = 0.0;

                auto block = LogLikelihoodBlock::MultivariateGaussian<2>(cache, obs, mean, covariance);

                LogLikelihood llh(p);
                llh.add(Constraint("test::correlated-gaussian-m_b-and-m_c", std::vector<ObservablePtr>(obs.begin(), obs.end()), { block }));

                LogPosterior log_posterior(llh);

                // parameter that affects the observable
                log_posterior.add(LogPrior::Flat(p, "mass::b(MSbar)", ParameterRange{ -10, 10 }));

                // the 2nd parameter
                // NOTE: it is crucial for the accuracy of the covariance estimate that the range be large enough
                //       to cover the gaussian peak out to many sigmas
                log_posterior.add(LogPrior::Flat(p, "mass::c", ParameterRange{ -10, 10 }));

                /* setup the MCMC sampler for the prerun to create the proposal */

                static const std::string mcmc_file_name = EOS_BUILDDIR "/eos/statistics/pmc_sampler_TEST-mcmc-prerun.hdf5";

                std::remove(mcmc_file_name.c_str());

                {
                    MarkovChainSampler::Config config = MarkovChainSampler::Config::Default();
                    config.chunk_size = 1;
                    config.chunks = 1;
                    config.number_of_chains = 10;
                    config.parallelize = true;
                    config.prerun_iterations_update = 650;
                    config.prerun_iterations_max = 2000;
                    config.prerun_iterations_min = 5000;
                    config.proposal_initial_covariance = proposal_covariance(log_posterior, 10);
                    config.output_file = mcmc_file_name;
                    config.seed = 784213135;
                    config.skip_initial = 0.2;

                    std::vector<std::tuple<std::string, double, double> > part;

                    //TODO set starting point
                    //part.push_back(std::make_tuple(std::string("mass::c"), -5.5, -4.5));
                    //part.push_back(std::make_tuple(std::string("mass::c"), +4.5, +5.5));

                    Log::instance()->set_log_level(ll_silent);
                    MarkovChainSampler sampler(log_posterior.clone(), config);
                    sampler.run();
                    Log::instance()->set_log_level(ll_debug);
                }

                /* initialize PMC from MCMC */
                static const std::string pmc_output = EOS_BUILDDIR "/eos/statistics/pmc_sampler_TEST-output.hdf5";
                static const std::string pmc_output_components = EOS_BUILDDIR "/eos/statistics/pmc_sampler_TEST-output-components.hdf5";
                static const std::string pmc_output_hc = EOS_BUILDDIR "/eos/statistics/pmc_sampler_TEST-output-hc.hdf5";
                static const std::string pmc_output_resume = EOS_BUILDDIR "/eos/statistics/pmc_sampler_TEST-output-resume.hdf5";
                static const std::string pmc_output_split = EOS_BUILDDIR "/eos/statistics/pmc_sampler_TEST-output-split.hdf5";

                PopulationMonteCarloSampler::Config pmc_config = PopulationMonteCarloSampler::Config::Default();
                pmc_config.max_updates = 5;
                pmc_config.samples_per_component = 2000;
                pmc_config.final_samples = 5000;
                pmc_config.output_file = pmc_output;
                pmc_config.parallelize = true;
                pmc_config.seed = 23;
                pmc_config.store = true;

                // perform complete run
                {
                	PopulationMonteCarloSampler::Config temp_config(pmc_config);
                    temp_config.skip_initial = 0.2;
                    temp_config.patch_length = 400;
                    temp_config.target_ncomponents = 2;
                    PopulationMonteCarloSampler pmc_sampler(log_posterior.clone(), hdf5::File::Open(mcmc_file_name), temp_config);
                    pmc_sampler.run();
                    TEST_CHECK(pmc_sampler.status().converged);
                }

                // save initial status for later resumption
                {
                	PopulationMonteCarloSampler::Config temp_config(pmc_config);
                    temp_config.output_file = pmc_output_components;
                    temp_config.skip_initial = 0.2;
                    temp_config.patch_length = 400;
                    temp_config.target_ncomponents = 2;
                    PopulationMonteCarloSampler pmc_sampler(log_posterior.clone(), hdf5::File::Open(mcmc_file_name), temp_config);
                    pmc_sampler.draw_samples();
                }

                // resuming from previous step
                {
                    pmc_config.output_file = pmc_output_resume;
                    PopulationMonteCarloSampler pmc_sampler(log_posterior.clone(), hdf5::File::Open(pmc_output_components), pmc_config);
                    pmc_sampler.run();
                    TEST_CHECK(pmc_sampler.status().converged);
                }

                // splitting up the calculation
                {
                    pmc_config.samples_per_component = 3001;
                    pmc_config.output_file = pmc_output_split;
                    PopulationMonteCarloSampler pmc_sampler(log_posterior.clone(), hdf5::File::Open(pmc_output_components), pmc_config);
                    pmc_sampler.calculate_weights(pmc_output_components, 0, pmc_config.samples_per_component - 1);
                }

                // hierarchical clustering: integrate over a subdomain only, eliminate one peak with half of probability
                // but limit prior ranges
                {
                    pmc_config.samples_per_component = 1000;
                    pmc_config.group_by_r_value = 1.1;
                    pmc_config.output_file = pmc_output_hc;
                    pmc_config.patch_length = 400;
                    pmc_config.skip_initial = 0.2;
                    pmc_config.target_ncomponents = 2;

                    LogPosterior ana(llh);
                    ana.add(LogPrior::Flat(p, "mass::b(MSbar)", ParameterRange{ -10, 10 }));
                    ana.add(LogPrior::Flat(p, "mass::c", ParameterRange{ -10, 0 }));
                    PopulationMonteCarloSampler pmc_sampler(ana.clone(), hdf5::File::Open(mcmc_file_name), pmc_config);
                    pmc_sampler.run();
                    TEST_CHECK(pmc_sampler.status().converged);
                }

                // read in results and check
                {
                    hdf5::File file = hdf5::File::Open(pmc_output);
                    hdf5::File file_resume = hdf5::File::Open(pmc_output_resume);
                    hdf5::File file_split = hdf5::File::Open(pmc_output_split);
                    hdf5::File file_hc = hdf5::File::Open(pmc_output_hc);

                    auto data_set = file.open_data_set("/data/final/statistics", PopulationMonteCarloSampler::Output::statistics_type());
                    auto data_set_resume = file_resume.open_data_set("/data/final/statistics", PopulationMonteCarloSampler::Output::statistics_type());
                    auto data_set_resume0 = file_resume.open_data_set("/data/0/samples", PopulationMonteCarloSampler::Output::sample_type(2));
                    auto data_set_split = file_split.open_data_set("/data/weights", PopulationMonteCarloSampler::Output::weight_type());
                    auto data_set_hc = file_hc.open_data_set("/data/final/statistics", PopulationMonteCarloSampler::Output::statistics_type());
                    data_set.end();
                    data_set_resume.end();

                    auto record = PopulationMonteCarloSampler::Output::statistics_record();
                    auto record_resume = PopulationMonteCarloSampler::Output::statistics_record();
                    auto record_hc = PopulationMonteCarloSampler::Output::statistics_record();

                    data_set >> record;
                    data_set_resume >> record_resume;
                    data_set_hc >> record_hc;

                    // evidence = 2 (from absolute value for par mass::c) times 1 / (prior ranges)
                    TEST_CHECK_RELATIVE_ERROR(2.0 / 20.0 / 20.0, std::get<2>(record), 5e-3);
                    TEST_CHECK_RELATIVE_ERROR(2.0 / 20.0 / 20.0, std::get<2>(record_resume), 5e-3);
                    TEST_CHECK_RELATIVE_ERROR(1.0 / 20.0 / 10.0, std::get<2>(record_hc), 5e-3);

                    // check that initial components are identical
                    auto data_set_comp = file.open_data_set("/data/initial/components", PopulationMonteCarloSampler::Output::component_type(2));
                    auto data_set_comp_resume = file_resume.open_data_set("/data/initial/components", PopulationMonteCarloSampler::Output::component_type(2));

                    auto record_comp = PopulationMonteCarloSampler::Output::component_record(2);
                    auto record_comp_resume = PopulationMonteCarloSampler::Output::component_record(2);
                    data_set_comp >> record_comp;
                    data_set_comp_resume >> record_comp_resume;

                    TEST_CHECK_EQUAL(std::get<0>(record_comp), std::get<0>(record_comp_resume));
                    TEST_CHECK_EQUAL(std::get<1>(record_comp)[0], std::get<1>(record_comp_resume)[0]);
                    TEST_CHECK_EQUAL(std::get<2>(record_comp)[0], std::get<2>(record_comp_resume)[0]);

                    // check random record
                    auto record_resume0 = PopulationMonteCarloSampler::Output::sample_record(2);
                    auto record_split = PopulationMonteCarloSampler::Output::weight_record();

                    for (unsigned i = 0 ; i < data_set_split.records() ; ++i)
                    {
                        data_set_resume0 >> record_resume0;
                        data_set_split >> record_split;

                        TEST_CHECK_NEARLY_EQUAL(record_resume0[3], std::get<0>(record_split), 1e-17); // posterior
                        // todo For some reason, values don't agree bitwise. Why?
                        TEST_CHECK_RELATIVE_ERROR(record_resume0[4], std::get<1>(record_split), 2e-14); // weight
                    }
                }
            }
        }
} population_monte_carlo_sampler_test;
