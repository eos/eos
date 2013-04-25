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

#include <test/test.hh>
#include <eos/utils/analysis_TEST.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/histogram.hh>
#include <eos/utils/markov_chain_sampler.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/proposal_functions.hh>

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


            // check pre run, main run and HDF5 storage
            {
                static const std::string file_name(EOS_BUILDDIR "/eos/utils/markov_chain_sampler_TEST.hdf5");
                std::remove(file_name.c_str());

                // store to HDF5
                TEST_SECTION("run-and-store",
                {
                    Analysis analysis(make_analysis(true));

                    MarkovChainSampler::Config config = MarkovChainSampler::Config::Quick();
                    config.chunk_size = 100;
                    config.chunks = 6;
                    config.max_efficiency = 0.75;
                    config.min_efficiency = 0.20;
                    config.need_prerun = true;
                    config.number_of_chains = 3;
                    config.output_file = file_name;
                    config.parallelize = true;
                    config.find_modes = true;
                    config.prerun_iterations_update = 500;
                    config.prerun_iterations_min = 1000;
                    config.rvalue_criterion_param = 1.1;
                    config.scale_automatic = true;
                    config.scale_reduction = 2.0;
                    config.seed = 1346;
                    config.store = true;
                    config.store_prerun = true;
                    config.use_posterior_rvalue = true;
                    config.use_strict_rvalue_definition = true;

                    MarkovChainSampler sampler(analysis, config);
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
                        // mode found during 2x500 iterations + the one from minuit
                        TEST_CHECK_EQUAL(data_set.records(), 3);

                        std::vector<double> record(2);
                        data_set.end();
                        data_set >> record;
                        TEST_CHECK_RELATIVE_ERROR(record[0], 4.2, 1e-7);
                        TEST_CHECK_RELATIVE_ERROR(record[1], 1.201325, 1e-5);
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
                        TEST_CHECK_EQUAL(std::get<0>(record_pre), std::string("test-observable[mass::b(MSbar)]"));

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

            // check local/global proposal function
            {
                /* setup bimodal distribution */

                Parameters p = Parameters::Defaults();
                Kinematics k;
                std::array<ObservablePtr, 2> obs;
                obs[0] = ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)"));
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
                llh.add(Constraint("Correlated Gaussian for m_b and m_c", std::vector<ObservablePtr>(obs.begin(), obs.end()), { block }));

                Analysis analysis(llh);

                // an interesting parameter that affects the observable
                analysis.add(LogPrior::Flat(p, "mass::b(MSbar)", ParameterRange{ -10, 10 }));

                // the 2nd parameter
                // NOTE: it is crucial for the accuracy of the covariance estimate that the range be large enough
                //       to cover the gaussian peak out to many sigmas
                analysis.add(LogPrior::Flat(p, "mass::c", ParameterRange{ -10, 10 }));


                /* setup the sampler */
                static const std::string file_name = EOS_BUILDDIR "/eos/utils/markov_chain_sampler_TEST-global-local.hdf5";
                static const std::string file_name_resume = EOS_BUILDDIR "/eos/utils/markov_chain_sampler_TEST-global-local-resume.hdf5";
                std::remove(file_name.c_str());

                MarkovChainSampler::Config config = MarkovChainSampler::Config::Default();
                config.number_of_chains = 2;
                config.chunk_size = 10000;
                config.chunks = 3;
                config.seed = 784213135;
                config.output_file = file_name;
                config.scale_reduction = 10;
                config.parallelize = true;
                config.find_modes = true;
                config.prerun_iterations_update = 650;
                config.prerun_iterations_max = 2000;
                config.prerun_iterations_min = 5000;
                config.scale_reduction = 10;
                config.skip_initial = 0.2;
                config.store_prerun = true;

                // more history points slow down the program significantly when likelihood is fast
                auto gl_config = proposal_functions::GlobalLocal::Config::Default();
                gl_config.join_chains_symmetrically = true;
                gl_config.local_jump_probability = 0.95;
                gl_config.perform_clustering = true;
                gl_config.skip_initial = config.skip_initial;
                config.global_local_config.reset(new proposal_functions::GlobalLocal::Config(gl_config));
                {
                    std::vector<std::tuple<std::string, double, double> > part;

                    part.push_back(std::make_tuple(std::string("mass::c"), -5.5, -4.5));
                    config.partitions.push_back(part);

                    part.clear();
                    part.push_back(std::make_tuple(std::string("mass::c"), +4.5, +5.5));
                    config.partitions.push_back(part);

                    MarkovChainSampler sampler(analysis, config);
                    sampler.run();
                }
                /* open HDF5 file and run checks on data */
                {
                    // read file back in
                    hdf5::File file = hdf5::File::Open(file_name, H5F_ACC_RDONLY);
                    hdf5::Array<1, double> sample_type
                    {
                        "samples",
                        { 2 + 1 },
                    };
                    auto data_set = file.open_data_set("/main run/chain #0/samples", sample_type);

                    unsigned n = data_set.records();
                    TEST_CHECK_EQUAL(n, config.chunks * config.chunk_size);

                    /* analyze histogram for mass::b, should be that of a Gaussian */

                    unsigned n_bins = 60;
                    Histogram<1> hist_b = Histogram<1>::WithEqualBinning(4, 6, n_bins);

                    double chi_squared = 0.0;

                    bin_data_set(data_set, hist_b, 0, mean[0], std::sqrt(covariance[0][0]), chi_squared);

                    // use sqrt(var) from Chi2 distribution as uncertainty
                    TEST_CHECK_RELATIVE_ERROR((double)n, (double)chi_squared, std::sqrt(2.0 * n));

                    std::vector<double> pdf_vec;

                    for (auto b = hist_b.begin(), b_end = hist_b.end() ; b != b_end ; ++b)
                    {
                        pdf_vec.push_back(b->value);
                    }

                    auto cumulative_hist = estimate_cumulative_distribution(hist_b);

                    std::vector<double> cumulative_vec;
                    for (auto b = cumulative_hist.begin(), b_end = cumulative_hist.end() ; b != b_end ; ++b)
                    {
                        cumulative_vec.push_back(b->value);
                    }

                    TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(0.5 * n_bins) - 1], 0.5,  2.5e-2);

                    // one sigma interval should be [4.9, 5.1], whole range covers n_sigmas
                    double n_sigmas = (6.0 - 4.0) / 0.1;
                    TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(11.0 / n_sigmas * n_bins) - 1] - cumulative_vec[unsigned(9.0 / n_sigmas * n_bins) - 1], 0.68,  2e-2);

                    // two sigma interval should be [4.8, 5.2], less precise
                    TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(12.0 / n_sigmas * n_bins) - 1] - cumulative_vec[unsigned(8.0 / n_sigmas * n_bins) - 1], 0.95,  3e-2);

                    /* check mass::c, should be two Gaussian of same height */

                    unsigned n_bins_c = 700;
                    Histogram<1> hist_c = Histogram<1>::WithEqualBinning(-6, +6, n_bins_c);

                    bin_data_set(data_set, hist_c, 1, mean[1], std::sqrt(covariance[1][1]), chi_squared);

                    pdf_vec.clear();

                    for (auto b = hist_c.begin(), b_end = hist_c.end() ; b != b_end ; ++b)
                    {
                        pdf_vec.push_back(b->value);
                    }

                    auto cumulative_hist_c = estimate_cumulative_distribution(hist_c);

                    cumulative_vec.clear();

                    for (auto b = cumulative_hist_c.begin(), b_end = cumulative_hist_c.end() ; b != b_end ; ++b)
                    {
                        cumulative_vec.push_back(b->value);
                    }

                    // same probability on either side of 0, both peaks equally strong
                    TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(0.5 * n_bins_c) - 1], 0.5,  3e-2);

                    // one sigma interval should be [4.95, 5.05] and [-5.05, -4.95], whole range covers 20 sigmas
                    n_sigmas = (6.0 - (-6.0)) / 0.05;
                    TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(21.0 / n_sigmas * n_bins_c) - 1] - cumulative_vec[unsigned(19.0 / n_sigmas * n_bins_c) - 1], 0.68 / 2.0,  5e-2);
                    TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned((1 - (19.0 / n_sigmas)) * n_bins_c) - 1] - cumulative_vec[unsigned((1 - (21.0 / n_sigmas)) * n_bins_c) - 1], 0.68 / 2.0,  5e-2);

                    /*  check pre run */

                    data_set = file.open_data_set("/prerun/chain #0/samples", sample_type);
                    TEST_CHECK(data_set.records() <=
                               std::max(config.prerun_iterations_min, config.prerun_iterations_max) + config.prerun_iterations_update);
                    TEST_CHECK(data_set.records() >= config.prerun_iterations_min);
                }

                // check proposal I/O in HDF5
                {
                    static const std::string file_name_build = EOS_BUILDDIR "/eos/utils/markov_chain_sampler_TEST-build-global-local.hdf5";

                    // read preruns, and store global local to disk
                    TEST_SECTION("build global local from disk",
                    {
                        proposal_functions::GlobalLocal::Config & gl_config = *config.global_local_config;
                        gl_config.join_chains_symmetrically = true;

                        std::vector<std::shared_ptr<hdf5::File>> input_files;
                        auto input_file = std::make_shared<hdf5::File>(hdf5::File::Open(file_name, H5F_ACC_RDONLY));
                        input_files.push_back(input_file);
                        TEST_CHECK(input_files.front()->group_exists("/prerun/chain #0"));
                        MarkovChainSampler::build_global_local(file_name_build, input_files, gl_config);
                    });

                    TEST_SECTION("read global local from disk and repeat main run",
                    {
                        hdf5::File file_build = hdf5::File::Open(file_name_build);

                        // check that it compiles and runs
                        proposal_functions::Factory::make(file_build, "/global local", "GlobalLocal", 2);

                        config.output_file = file_name_resume;
                        // avoid extending file name with every check just to avoid overwriting
                        std::remove(file_name_resume.c_str());
                        MarkovChainSampler sampler(analysis, config);
                        sampler.resume(file_build);
                    });
                }

                // do results agree, whether I resume or not?
                {
                    hdf5::File f = hdf5::File::Open(file_name);
                    hdf5::File g = hdf5::File::Open(file_name_resume);

                    //todo create factories with dimension argument to create all HDF5 output type => one definition!
                    hdf5::Array<1, double> sample_type
                    {
                        "samples",
                        { 2 + 1 },
                    };
                    auto data_set_f = f.open_data_set("/main run/chain #0/samples", sample_type);
                    auto data_set_g = g.open_data_set("/main run/chain #0/samples", sample_type);

                    std::vector<double> record_f(3);
                    std::vector<double> record_g(3);

                    data_set_f >> record_f;
                    data_set_g >> record_g;

                    TEST_CHECK_EQUAL(record_f, record_g);

                    data_set_f.end();
                    data_set_g.end();

                    data_set_f >> record_f;
                    data_set_g >> record_g;

                    TEST_CHECK_EQUAL(record_f, record_g);

                    data_set_f = f.open_data_set("/main run/chain #1/samples", sample_type);
                    data_set_g = g.open_data_set("/main run/chain #1/samples", sample_type);

                    data_set_f >> record_f;
                    data_set_g >> record_g;

                    TEST_CHECK_EQUAL(record_f, record_g);

                    data_set_f.end();
                    data_set_g.end();

                    data_set_f >> record_f;
                    data_set_g >> record_g;

                    TEST_CHECK_EQUAL(record_f, record_g);
                }
            }
        }
} markov_chain_sampler_test;
