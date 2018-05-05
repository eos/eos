/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Frederik Beaujean
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

#include <eos/statistics/markov-chain-sampler.hh>

#include <test/test.hh>
#include <eos/statistics/density-wrapper_TEST.hh>
#include <eos/statistics/histogram.hh>
#include <eos/statistics/log-posterior_TEST.hh>
#include <eos/statistics/proposal-functions.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/power_of.hh>

using namespace test;
using namespace eos;

template <typename T_>
void bin_data_set(hdf5::DataSet<T_ > & data_set, Histogram<1> & hist, const unsigned & dimension,
                  const double & mu, const double & sigma, double & chi_squared)
{
    chi_squared = 0;

    // reset pointer to start reading from beginning again
    data_set.set_index(0);

    std::vector<double> record;
    for (unsigned i = 0 ; i < data_set.records() ; ++i)
    {
        data_set >> record;
        hist.insert(record[dimension]);
        chi_squared += power_of<2>((record[dimension] - mu) / sigma);
    }
}

class MarkovChainSamplerTest :
    public TestCase
{
    public:
        MarkovChainSamplerTest() :
            TestCase("markov_chain_sampler_test")
        {
        }

        virtual void run() const
        {
            // check MarkovChainSampler::Config
            TEST_SECTION("config",
            {
                MarkovChainSampler::Config conf = MarkovChainSampler::Config::Default();
                TEST_CHECK_THROWS(VerifiedRangeUnderflow, conf.min_efficiency = -0.1);

                MarkovChainSampler::Config conf_quick = MarkovChainSampler::Config::Quick();
                TEST_CHECK_THROWS(VerifiedRangeOverflow, conf.max_efficiency = 23.1);
            });


            /* sample from a density */
            {
                DensityWrapper density = make_multivariate_unit_normal(2);
                MarkovChainSampler::Config config = MarkovChainSampler::Config::Default();
                config.chunk_size = 5000;
                config.chunks = 10;
                config.number_of_chains = 2;
                config.output_file = EOS_BUILDDIR "/eos/statistics/markov-chain-sampler_TEST_density.hdf5";
                config.parallelize = false;
                config.seed = 1246122;
                MarkovChainSampler sampler(density.clone(), config);
                sampler.run();
            }

            // check pre run, main run and HDF5 storage
            {
                static const std::string file_name(EOS_BUILDDIR "/eos/statistics/markov-chain-sampler_TEST.hdf5");
                std::remove(file_name.c_str());

                // store to HDF5
                TEST_SECTION("run-and-store",
                {
                    LogPosterior log_posterior(make_log_posterior(true));

                    MarkovChainSampler::Config config = MarkovChainSampler::Config::Quick();
                    config.chunk_size = 100;
                    config.chunks = 6;
                    config.max_efficiency = 0.75;
                    config.min_efficiency = 0.20;
                    config.need_prerun = true;
                    config.number_of_chains = 3;
                    config.output_file = file_name;
                    config.parallelize = true;
                    config.prerun_iterations_update = 500;
                    config.prerun_iterations_min = 1000;
                    config.proposal_initial_covariance = proposal_covariance(log_posterior, 2);
                    config.rvalue_criterion_param = 1.1;
                    config.scale_automatic = true;
                    config.seed = 1346;
                    config.store = true;
                    config.store_prerun = true;
                    config.use_posterior_rvalue = true;
                    config.use_strict_rvalue_definition = true;

                    MarkovChainSampler sampler(log_posterior.clone(), config);
                    sampler.run();

                    MarkovChainSampler::PreRunInfo pre_info(sampler.pre_run_info());

                    TEST_CHECK(pre_info.converged);
                    TEST_CHECK_EQUAL(pre_info.iterations_at_convergence, pre_info.iterations);
                    TEST_CHECK_EQUAL(pre_info.iterations_at_convergence, 1000);
                    TEST_CHECK_NEARLY_EQUAL(pre_info.rvalue_parameters[0], 1.0, 5e-3);
                });

                // check sizes of data sets
                {
                    auto f = hdf5::File::Open(file_name);
                    hdf5::Array<1, double> sample_type
                    {
                        "samples",
                        { 1 + 1 },
                    };

                    {
                        auto data_set = f.open_data_set("/prerun/chain #0/proposal/meta", proposal_functions::meta_type());
                        TEST_CHECK_EQUAL(data_set.records(), 1);

                        auto meta_record = proposal_functions::meta_record();
                        data_set >> meta_record;

                        TEST_CHECK_EQUAL(std::get<0>(meta_record), std::string("MultivariateGaussian"));
                        TEST_CHECK_EQUAL(std::get<1>(meta_record), 1u);
                    }
                    {
                        auto data_set_pre = f.open_data_set("/prerun/chain #1/samples", sample_type);
                        TEST_CHECK_EQUAL(data_set_pre.records(), 1000);

                        auto data_set_main = f.open_data_set("/main run/chain #1/samples", sample_type);
                        TEST_CHECK_EQUAL(data_set_main.records(), 600);
                    }
                    {
                        auto data_set = f.open_data_set("/prerun/chain #0/stats/mode", sample_type);
                        // mode found during 2x500 iterations
                        TEST_CHECK_EQUAL(data_set.records(), 2);

                        std::vector<double> record(2);
                        data_set.end();
                        data_set >> record;
                        TEST_CHECK_RELATIVE_ERROR(record[0], 4.2, 1e-4);
                        TEST_CHECK_RELATIVE_ERROR(record[1], 1.201325, 1e-4);
                    }
                    {
                        hdf5::Composite<hdf5::Scalar<const char *>, hdf5::Scalar<double>,
                        hdf5::Scalar<double>, hdf5::Scalar<int>, hdf5::Scalar<const char *>> parameter_descriptions_type
                        {
                            "parameter description",
                            hdf5::Scalar<const char *>("name"),
                            hdf5::Scalar<double>("min"),
                            hdf5::Scalar<double>("max"),
                            hdf5::Scalar<int>("nuisance"),
                            hdf5::Scalar<const char *>("prior"),
                        };
                        auto data_set_pre = f.open_data_set("/descriptions/prerun/chain #2/parameters", parameter_descriptions_type);
                        TEST_CHECK_EQUAL(data_set_pre.records(), 1);

                        auto record_pre = std::make_tuple("parameter_name", 1.0, 2.0, 3, "prior");
                        data_set_pre >> record_pre;
                        TEST_CHECK_EQUAL(std::get<0>(record_pre), std::string("mass::b(MSbar)"));
                        TEST_CHECK_EQUAL(std::get<1>(record_pre), 3.7);
                        TEST_CHECK_EQUAL(std::get<2>(record_pre), 4.9);
                        TEST_CHECK_EQUAL(std::get<3>(record_pre), false);
                        TEST_CHECK_EQUAL(std::get<4>(record_pre), std::string("Parameter: mass::b(MSbar), prior type: flat, range: [3.7,4.9]"));

                        auto data_set_main = f.open_data_set("/descriptions/main run/chain #2/parameters", parameter_descriptions_type);
                        TEST_CHECK_EQUAL(data_set_main.records(), 1);

                        auto record_main = std::make_tuple("parameter_name", 1.0, 2.0, 3, "prior");
                        data_set_main >> record_main;

                        TEST_CHECK_EQUAL(std::string(std::get<0>(record_pre)), std::string(std::get<0>(record_main)));
                        TEST_CHECK_EQUAL(std::get<1>(record_pre), std::get<1>(record_main));
                        TEST_CHECK_EQUAL(std::get<2>(record_pre), std::get<2>(record_main));
                        TEST_CHECK_EQUAL(std::get<3>(record_pre), std::get<3>(record_main));
                        TEST_CHECK_EQUAL(std::string(std::get<4>(record_pre)), std::string(std::get<4>(record_main)));
                    }
                    {
                        hdf5::Composite<hdf5::Scalar<const char *>> constraint_type
                        {
                            "constraints",
                            hdf5::Scalar<const char *>("name"),
                        };
                        auto data_set_pre = f.open_data_set("/descriptions/prerun/chain #1/constraints", constraint_type);
                        auto record_pre = std::make_tuple("parameter_name");
                        data_set_pre >> record_pre;
                        TEST_CHECK_EQUAL(std::get<0>(record_pre), std::string("mass::b(MSbar)"));

                        auto data_set_main = f.open_data_set("/descriptions/main run/chain #1/constraints", constraint_type);
                        auto record_main = std::make_tuple("parameter_name");
                        data_set_main >> record_main;
                        TEST_CHECK_EQUAL(std::string(std::get<0>(record_pre)), std::string(std::get<0>(record_main)));
                    }
                    {
                        hdf5::Array<1, double> covariance_type
                        {
                            "samples",
                            { 1 * 1 },
                        };
                        auto data_set_pre_0 = f.open_data_set("/prerun/chain #0/proposal/covariance", covariance_type);
                        std::vector<double> record_0(1);
                        data_set_pre_0 >> record_0;

                        auto data_set_pre_1 = f.open_data_set("/prerun/chain #1/proposal/covariance", covariance_type);
                        std::vector<double> record_1(1);
                        data_set_pre_1 >> record_1;

                        auto data_set_pre_2 = f.open_data_set("/prerun/chain #2/proposal/covariance", covariance_type);
                        std::vector<double> record_2(1);
                        data_set_pre_2 >> record_2;

                        /* covariances identical in first round.
                         * values from variance of uniform prior 1.2**2 / 12.0 = 0.12 and scale reduction 2**2
                         * and the scaling inside the proposal, if it is taken into account, 2.38**2.
                         */
                        TEST_CHECK_RELATIVE_ERROR(record_0[0], 0.03 * 2.38 * 2.38, 1e-15);
                        TEST_CHECK_EQUAL(record_0[0], record_1[0]);
                        TEST_CHECK_EQUAL(record_0[0], record_2[0]);
                        TEST_CHECK_EQUAL(record_1[0], record_2[0]);

                        data_set_pre_0 >> record_0;
                        data_set_pre_1 >> record_1;
                        data_set_pre_2 >> record_2;

                        // but different after first adapt
                        TEST_CHECK(record_0[0] != record_1[0]);
                        TEST_CHECK(record_0[0] != record_2[0]);
                        TEST_CHECK(record_1[0] != record_2[0]);
                    }
                }
            }
        }
} markov_chain_sampler_test;
