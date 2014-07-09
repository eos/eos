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

#include <eos/utils/population_monte_carlo_sampler.hh>
#include <eos/utils/analysis_TEST.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/markov_chain_sampler.hh>
#include <test/test.hh>

using namespace test;
using namespace eos;

namespace
{
struct ThesisExample
{
    struct Config
    {
        // extra dimensions on top of two
        unsigned n_extra_dim;

        // flat dimensions
        unsigned n_flat_dim;

        // relative weight of modes
        VerifiedRange<double> weight;

        static Config Default()
        {
            return Config();
        }

    private:
        Config() :
            n_extra_dim(0),
            n_flat_dim(0),
            weight(1, std::numeric_limits<double>::max(), 1)
        {
        }
    };

    Kinematics k;
    Parameters p;
    ObservableCache cache;

    // defines the position of the mode
    const double nu; // = 10;
    const double theta;
    const double alpha;

    // extract proper '1sigma' boundaries by solving the usual constraints with nu, theta, alpha fixed and a,b variable
    const double a;
    const double b;

    ThesisExample::Config config;
    Analysis analysis;

    ThesisExample(const ThesisExample::Config & config) :
        p(Parameters::Defaults()),
        cache(p),
        nu(10),
        theta(1),
        alpha(1),
        a(-1.315856893331208),
        b(0.916552003733646),
        config(config),
        analysis(make_analysis())
    {
    }

    Analysis make_analysis()
    {
        // index offset into available observables / parameters
        unsigned offset = 0;

        std::vector<LogLikelihoodBlockPtr> loggamma_components
        {
            LogLikelihoodBlock::LogGamma(cache, observables(0), -nu + a, -nu, -nu + b, theta, alpha),
                LogLikelihoodBlock::LogGamma(cache, observables(0), nu + a, nu, nu + b, theta, alpha)
        };
        std::vector<double> loggamma_weights{ 1, 1 };

        std::vector<LogLikelihoodBlockPtr> blocks;

        blocks.push_back(LogLikelihoodBlock::Mixture(loggamma_components, loggamma_weights));

        std::vector<LogLikelihoodBlockPtr> gaussian_components
        {
            LogLikelihoodBlock::Gaussian(cache, observables(1), -nu - 1, -nu, -nu + 1),
                LogLikelihoodBlock::Gaussian(cache, observables(1),  nu - 1, nu, nu +1)
        };
        std::vector<double> weights { 1, config.weight };
        std::cout << config.weight << std::endl;

        blocks.push_back(LogLikelihoodBlock::Mixture(gaussian_components, weights));

        LogLikelihood llh(p);

        llh.add(Constraint("Bimodal LogGamma distribution",
                           std::vector<ObservablePtr> { observables(0) },
                           std::vector<LogLikelihoodBlockPtr> { blocks[0] }));
        ++offset;

        llh.add(Constraint("Bimodal Gaussian distribution",
                           std::vector<ObservablePtr> { observables(1) },
                           std::vector<LogLikelihoodBlockPtr> { blocks[1] }));
        ++offset;

        // add extra dimensions
        for (unsigned j = 0 ; j < config.n_extra_dim; ++j)
        {
            auto block = LogLikelihoodBlock::LogGamma(cache, observables(offset + j), nu + a, nu, nu + b, theta, alpha);
            llh.add(Constraint("Single LogGamma",
                           std::vector<ObservablePtr> { observables(offset + j) },
                           std::vector<LogLikelihoodBlockPtr> { block }));
        }
        offset += config.n_extra_dim;

        for (unsigned j = 0 ; j < config.n_extra_dim ; ++j)
        {
            auto block = LogLikelihoodBlock::Gaussian(cache, observables(offset + j),  nu - 1, nu, nu +1);
            llh.add(Constraint("Single Gaussian",
                           std::vector<ObservablePtr> { observables(offset + j) },
                           std::vector<LogLikelihoodBlockPtr> { block }));
        }
        offset += config.n_extra_dim;

        // add flat dimensions
        offset += config.n_flat_dim;

        Analysis analysis(llh);

        // add all necessary parameters, all nuisance except first two
        for (unsigned j = 0 ; j < offset ; ++j)
        {
            bool nuisance = (j >= 2) ? true : false;
            // cut on parameter name from "test-observable[B->K::b^t_1@KMPW2010]"
            std::string name = observables(j)->name();
            analysis.add(LogPrior::Flat(p, name.substr(16, name.size() - 17), ParameterRange{ -3 * nu, 3 * nu }), nuisance);
        }
        return analysis;
    }

    ObservablePtr observables(const unsigned & index)
    {
        static std::vector<ObservablePtr> observables
        {
            ObservablePtr(new TestObservable(p, k, "mass::e")),
                    ObservablePtr(new TestObservable(p, k, "mass::mu")),
                    ObservablePtr(new TestObservable(p, k, "mass::tau")),
                    ObservablePtr(new TestObservable(p, k, "mass::s(2GeV)")),
                    ObservablePtr(new TestObservable(p, k, "mass::c")),
                    ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")),
                    ObservablePtr(new TestObservable(p, k, "mass::t(pole)")),
                    ObservablePtr(new TestObservable(p, k, "mass::K_d")),
                    ObservablePtr(new TestObservable(p, k, "mass::K^*_d")),
                    ObservablePtr(new TestObservable(p, k, "mass::B_d")),
                    ObservablePtr(new TestObservable(p, k, "mass::B_u")),
                    ObservablePtr(new TestObservable(p, k, "mass::B_s")),
                    ObservablePtr(new TestObservable(p, k, "mass::W")),
                    ObservablePtr(new TestObservable(p, k, "mass::Z")),
                    ObservablePtr(new TestObservable(p, k, "decay-constant::B_d")),
                    ObservablePtr(new TestObservable(p, k, "decay-constant::B_u")),
                    ObservablePtr(new TestObservable(p, k, "decay-constant::B_s")),
                    ObservablePtr(new TestObservable(p, k, "life_time::B_d")),
                    ObservablePtr(new TestObservable(p, k, "life_time::B_u")),
                    ObservablePtr(new TestObservable(p, k, "life_time::B_s")),
                    ObservablePtr(new TestObservable(p, k, "B->K^*::F^V(0)@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K^*::F^A0(0)@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K^*::F^A1(0)@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K^*::F^A2(0)@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K^*::b^V_1@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K^*::b^A0_1@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K^*::b^A1_1@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K^*::b^A2_1@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::F^p(0)@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::F^0(0)@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::F^t(0)@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::b^p_1@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::b^0_1@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::b^t_1@KMPW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::alpha^V0_0@BFW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::alpha^V0_1@BFW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::alpha^Vt_0np@BFW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::alpha^Vt_1np@BFW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::alpha^T0_0@BFW2010")),
                    ObservablePtr(new TestObservable(p, k, "B->K::alpha^T0_1@BFW2010")),
                    ObservablePtr(new TestObservable(p, k, "CKM::A")),
                    ObservablePtr(new TestObservable(p, k, "CKM::lambda")),
                    ObservablePtr(new TestObservable(p, k, "CKM::rhobar")),
                    ObservablePtr(new TestObservable(p, k, "CKM::etabar")),
                    ObservablePtr(new TestObservable(p, k, "QED::alpha_e(m_b)")),
                    ObservablePtr(new TestObservable(p, k, "QCD::mu_t")),
                    ObservablePtr(new TestObservable(p, k, "QCD::alpha_s(MZ)")),
                    ObservablePtr(new TestObservable(p, k, "QCD::mu_b")),
                    ObservablePtr(new TestObservable(p, k, "QCD::mu_c")),
                    ObservablePtr(new TestObservable(p, k, "QCD::Lambda")),
        };

        if (index > observables.size())
            throw InternalError("ThesisExample::observable(): index" + stringify(index) + " out of bounds.");

        return observables[index];
    }

    void mcmc(const MarkovChainSampler::Config & mcmc_config)
    {
        std::remove(mcmc_config.output_file.c_str());
        MarkovChainSampler sampler(analysis, mcmc_config);
        sampler.run();

        // simple check if results stored
        std::vector<std::shared_ptr<hdf5::File>> input_files;
        auto input_file = std::make_shared<hdf5::File>(hdf5::File::Open(mcmc_config.output_file, H5F_ACC_RDONLY));
        input_files.push_back(input_file);
        TEST_CHECK(input_files.front()->group_exists("/prerun/chain #0"));
    }

    void pmc(const PopulationMonteCarloSampler::Config & pmc_config, const std::string & mcmc_file_name)
    {
        PopulationMonteCarloSampler pmc_sampler(analysis, hdf5::File::Open(mcmc_file_name), pmc_config);
        pmc_sampler.run();
    }
};
}

class PopulationMonteCarloSamplerTest :
    public TestCase
{
    public:
        PopulationMonteCarloSamplerTest() :
            TestCase("population_monte_carlo_sampler_test")
        {
        }

        void thesis_examples() const
        {
            /* set up default config options geared at 2D */

            static const std::string output_dir = "/.th/pcl128c/scratch/beaujean/diss/";

            MarkovChainSampler::Config mcmc_config = MarkovChainSampler::Config::Default();
            {
                mcmc_config.chunk_size = 1;
                mcmc_config.chunks = 1;
                mcmc_config.find_modes = false;
                mcmc_config.need_main_run = false;
                mcmc_config.parallelize = true;
                mcmc_config.prerun_chains_per_partition = 20;
                mcmc_config.prerun_iterations_update = 500;
                mcmc_config.prerun_iterations_max = 2000;
                mcmc_config.prerun_iterations_min = 40000;
                mcmc_config.output_file = "MISSING FILENAME";
                mcmc_config.scale_reduction = 10;
                mcmc_config.seed = 78421313;
                mcmc_config.skip_initial = 0.2;
                mcmc_config.store_prerun = true;
            }

            PopulationMonteCarloSampler::Config pmc_config = PopulationMonteCarloSampler::Config::Default();
            {
                pmc_config.chunks = 20;
                pmc_config.chunk_size = 500;
                pmc_config.degrees_of_freedom = 12;
                pmc_config.final_chunk_size = 5e5;
                pmc_config.group_by_r_value = 1.2;
                //               pmc_config.ignore_groups = std::vector<unsigned>{ 1,2,3 };
                pmc_config.maximum_relative_std_deviation = 0.02;
                pmc_config.minimum_overlap = 0.02;
                pmc_config.output_file = "MISSING FILENAME";
                pmc_config.parallelize = true;
                pmc_config.seed = 23;
                pmc_config.sliding_window = 100;
                pmc_config.skip_initial = 0.2;
                pmc_config.store = true;
                pmc_config.super_clusters = 3;
            }

            Log::instance()->set_log_level(ll_informational);

            // can only run one at a time. So never accidentally use setting by a previous example.

#if 1
            ::ThesisExample::Config config = ::ThesisExample::Config::Default();
            static const std::string id = output_dir + "simple-2D-full-";
            mcmc_config.output_file = id + "prerun.hdf5";
            pmc_config.output_file = id + "pmc.hdf5";
            ::ThesisExample ex(config);
            mcmc_config.prerun_iterations_min = 200000;
            ex.mcmc(mcmc_config);
            ex.pmc(pmc_config, mcmc_config.output_file);
#endif

#if 0
            ::ThesisExample::Config config = ::ThesisExample::Config::Default();
            config.weight = 1e5;
            static const std::string id = output_dir + "simple-2D-half-";

            mcmc_config.output_file = id + "prerun.hdf5";
            pmc_config.output_file = id + "pmc.hdf5";
            ::ThesisExample ex(config);
//            ex.mcmc(mcmc_config);
            ex.pmc(pmc_config, mcmc_config.output_file);
#endif

#if 0
            ::ThesisExample::Config config = ::ThesisExample::Config::Default();
            config.weight = 1e5;
            config.n_extra_dim = 10;
            ::ThesisExample ex(config);

            static const std::string id = output_dir + "simple-ND-half-";
            mcmc_config.output_file = id + "prerun.hdf5";
            pmc_config.output_file = id + "pmc.hdf5";

            // onle one mode
            ex.analysis.restrict("mass::e", 0, 3 * ex.nu);
            ex.analysis.restrict("mass::mu", 0, 3 * ex.nu);

            mcmc_config.prerun_chains_per_partition = 5;
            mcmc_config.prerun_iterations_min = 5e4;
            mcmc_config.prerun_iterations_update = 1000;
            ex.mcmc(mcmc_config);

            pmc_config.chunk_size = 1500;
//            pmc_config.degrees_of_freedom = -1;
            pmc_config.ignore_eff_sample_size = true;
            pmc_config.sliding_window = 150;
            pmc_config.super_clusters = ex.analysis.parameter_descriptions().size() + 10;
            ex.pmc(pmc_config, mcmc_config.output_file);
#endif

#if 0
            ::ThesisExample::Config config = ::ThesisExample::Config::Default();
            config.weight = 1e5;
            config.n_extra_dim = 4;
            config.n_flat_dim = 2;
            ::ThesisExample ex(config);

            static const std::string id = output_dir + "simple-ND-flat-";

            mcmc_config.output_file = id + "prerun.hdf5";
            pmc_config.output_file = id + "pmc.hdf5";

            // onle one mode
            ex.analysis.restrict("mass::e", 0, 3 * ex.nu);
            ex.analysis.restrict("mass::mu", 0, 3 * ex.nu);

            mcmc_config.prerun_chains_per_partition = 5;
            mcmc_config.prerun_iterations_min = 5e4;
            mcmc_config.prerun_iterations_update = 1000;
            ex.mcmc(mcmc_config);

            pmc_config.chunk_size = 1000;
//            pmc_config.degrees_of_freedom = 5;
            pmc_config.ignore_eff_sample_size = true;
            pmc_config.sliding_window = 200;
            pmc_config.super_clusters = ex.analysis.parameter_descriptions().size() + 5;

            ex.pmc(pmc_config, mcmc_config.output_file);
#endif

        }

        virtual void run() const
        {
            thesis_examples();
#if 0
            // convergence checking
            {
                PopulationMonteCarloSampler::Config config = PopulationMonteCarloSampler::Config::Default();
                config.output_file = EOS_BUILDDIR "/eos/utils/pmc_sampler_TEST-convergence.hdf5";
                config.convergence_eff_sample_size= 0.9;
                config.convergence_perplexity = 0.9;
                config.minimum_eff_sample_size = 0.5;
                config.minimum_perplexity = 0.5;
                config.minimum_steps = 3;
                config.maximum_relative_std_deviation = 1e-2;

                bool flat = true;
                Analysis analysis(make_analysis(flat));
                PopulationMonteCarloSampler pmc_sampler(analysis, config);

                /* direct convergence based on current status only */

                PopulationMonteCarloSampler::Status new_status;
                new_status.eff_sample_size = 0.18;
                new_status.perplexity = 0.98;

                static const bool check_for_convergence = true;
                bool converged = pmc_sampler.status(new_status, check_for_convergence);
                TEST_CHECK(! converged);

                new_status.eff_sample_size = 0.98;
                converged = pmc_sampler.status(new_status);
                TEST_CHECK(converged);

                /* indirect convergence based on previous, mock-up values */

                {
                    auto file = hdf5::File::Create(config.output_file);

                    std::vector<double> eff_sample_sizes { 0.5712, 0.5739, 0.5698 } ;
                    std::vector<double> perplexities { 0.7942, 0.7876, 0.7956 };

                    TEST_CHECK_EQUAL(eff_sample_sizes.size(), config.minimum_steps);
                    TEST_CHECK_EQUAL(perplexities.size(), config.minimum_steps);

                    for (unsigned i = 0 ; i < config.minimum_steps ; ++i)
                    {
                        auto statistics = file.create_data_set("/data/" + stringify(i) + "/statistics",
                            PopulationMonteCarloSampler::Output::statistics_type());
                        auto statistics_record = std::make_tuple(perplexities[i], eff_sample_sizes[i], 11.11);
                        statistics << statistics_record;
                    }
                }

                new_status.eff_sample_size = 0.1;
                converged = pmc_sampler.status(new_status);
                TEST_CHECK(converged);

                {
                    auto file = hdf5::File::Create(config.output_file);

                    std::vector<double> eff_sample_sizes { 0.5712, 0.61, 0.68 } ;
                    std::vector<double> perplexities { 0.7942, 0.82, 0.853 };

                    TEST_CHECK_EQUAL(eff_sample_sizes.size(), config.minimum_steps);
                    TEST_CHECK_EQUAL(perplexities.size(), config.minimum_steps);

                    for (unsigned i = 0 ; i < config.minimum_steps ; ++i)
                    {
                        auto statistics = file.create_data_set("/data/" + stringify(i) + "/statistics",
                            PopulationMonteCarloSampler::Output::statistics_type());
                        auto statistics_record = std::make_tuple(perplexities[i], eff_sample_sizes[i], 11.11);
                        statistics << statistics_record;
                    }
                }

                new_status.eff_sample_size = 0.1;
                converged = pmc_sampler.status(new_status, check_for_convergence);
                TEST_CHECK(! converged);
            }

            // initialization with minuit
            {
                bool flat = true;
                Analysis analysis(make_analysis(flat));

                PopulationMonteCarloSampler::Config config = PopulationMonteCarloSampler::Config::Default();
                config.chunks = 2;
                config.chunk_size = 500;
                config.component_weights = std::vector<double>(2, 1/2.0);
                config.final_chunk_size = 500;
                config.mode_distance = 1e-3;
                config.output_file = "/tmp/pmc_parallel.hdf5";
                config.parallelize = true;
                config.random_start = false;
                config.seed = 23;
                config.starting_points = 15;
                config.store = true;

                PopulationMonteCarloSampler pop(analysis, config);
                pop.run();
            }
            // initialize from a MCMC prerun
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

                // parameter that affects the observable
                analysis.add(LogPrior::Flat(p, "mass::b(MSbar)", ParameterRange{ -10, 10 }));

                // the 2nd parameter
                // NOTE: it is crucial for the accuracy of the covariance estimate that the range be large enough
                //       to cover the gaussian peak out to many sigmas
                analysis.add(LogPrior::Flat(p, "mass::c", ParameterRange{ -10, 10 }));

                /* setup the MCMC sampler for the prerun to create the global local proposal */

                static const std::string mcmc_file_name = EOS_BUILDDIR "/eos/utils/pmc_sampler_TEST-mcmc-prerun.hdf5";
                static const std::string global_local_input = EOS_BUILDDIR "/eos/utils/pmc_sampler_TEST-global-local.hdf5";

                std::remove(mcmc_file_name.c_str());

                {
                    MarkovChainSampler::Config config = MarkovChainSampler::Config::Default();
                    config.chunk_size = 1;
                    config.chunks = 1;
                    config.find_modes = true;
                    config.number_of_chains = 2;
                    config.parallelize = true;
                    config.prerun_iterations_update = 650;
                    config.prerun_iterations_max = 2000;
                    config.prerun_iterations_min = 5000;
                    config.output_file = mcmc_file_name;
                    config.scale_reduction = 10;
                    config.seed = 784213135;
                    config.scale_reduction = 10;
                    config.skip_initial = 0.2;
                    config.store_prerun = true;

                    // more history points slow down the program significantly when likelihood is fast
                    auto gl_config = proposal_functions::GlobalLocal::Config::Default();
                    gl_config.history_points = 50;
                    gl_config.history_points_local_covariance_size = 500;
                    gl_config.join_chains_symmetrically = true;
                    gl_config.perform_clustering = true;
                    gl_config.skip_initial = config.skip_initial;
                    config.global_local_config.reset(new proposal_functions::GlobalLocal::Config(gl_config));
                    std::vector<std::tuple<std::string, double, double> > part;

                    part.push_back(std::make_tuple(std::string("mass::c"), -5.5, -4.5));
                    config.partitions.push_back(part);

                    part.clear();
                    part.push_back(std::make_tuple(std::string("mass::c"), +4.5, +5.5));
                    config.partitions.push_back(part);

                    Log::instance()->set_log_level(ll_silent);
                    MarkovChainSampler sampler(analysis, config);
                    sampler.run();
                    Log::instance()->set_log_level(ll_debug);

                    std::vector<std::shared_ptr<hdf5::File>> input_files;
                    auto input_file = std::make_shared<hdf5::File>(hdf5::File::Open(mcmc_file_name, H5F_ACC_RDONLY));
                    input_files.push_back(input_file);
                    TEST_CHECK(input_files.front()->group_exists("/prerun/chain #0"));
                    MarkovChainSampler::build_global_local(global_local_input, input_files, gl_config);
                }

                /* initialize PMC from MCMC */
                static const std::string pmc_output = EOS_BUILDDIR "/eos/utils/pmc_sampler_TEST-output.hdf5";
                static const std::string pmc_output_components = EOS_BUILDDIR "/eos/utils/pmc_sampler_TEST-output-components.hdf5";
                static const std::string pmc_output_hc = EOS_BUILDDIR "/eos/utils/pmc_sampler_TEST-output-hc.hdf5";
                static const std::string pmc_output_resume = EOS_BUILDDIR "/eos/utils/pmc_sampler_TEST-output-resume.hdf5";
                static const std::string pmc_output_split = EOS_BUILDDIR "/eos/utils/pmc_sampler_TEST-output-split.hdf5";

                PopulationMonteCarloSampler::Config pmc_config = PopulationMonteCarloSampler::Config::Default();
                pmc_config.chunks = 5;
                pmc_config.chunk_size = 2000;
                pmc_config.components_per_cluster = 1;
                pmc_config.final_chunk_size = 5000;
                pmc_config.output_file = pmc_output;
                pmc_config.parallelize = true;
                pmc_config.random_start = false;
                pmc_config.seed = 23;
                pmc_config.store = true;

                // save initial status for later resumption
                {
                    std::string temp = pmc_config.output_file;
                    pmc_config.output_file = pmc_output_components;
                    PopulationMonteCarloSampler pmc_sampler(analysis, hdf5::File::Open(global_local_input), pmc_config);
                    pmc_sampler.draw_samples();
                    pmc_config.output_file = temp;
                }

                {
                    PopulationMonteCarloSampler pmc_sampler(analysis, hdf5::File::Open(global_local_input), pmc_config);
                    pmc_sampler.run();
                }

                // resuming from previous step
                {
                    pmc_config.output_file = pmc_output_resume;
                    PopulationMonteCarloSampler pmc_sampler_resume(analysis, hdf5::File::Open(pmc_output_components), pmc_config);
                    pmc_sampler_resume.run();
                }

                // splitting up the calculation
                {
                    pmc_config.chunk_size = 3001;
                    pmc_config.output_file = pmc_output_split;
                    PopulationMonteCarloSampler pmc_sampler(analysis, hdf5::File::Open(pmc_output_components), pmc_config);
                    pmc_sampler.calculate_weights(pmc_output_components, 0, pmc_config.chunk_size - 1);
                }

                // hierarchical clustering: integrate over a subdomain only, eliminate one peak with half of probability
                {
                    pmc_config.chunk_size = 1000;
                    pmc_config.group_by_r_value = 1.1;
                    pmc_config.minimum_overlap = 0.02;
                    pmc_config.output_file = pmc_output_hc;
                    pmc_config.patch_around_local_mode = true;
                    pmc_config.sliding_window = 400;
                    pmc_config.skip_initial = 0.2;
                    pmc_config.super_clusters = 2;

                    Analysis ana(llh);
                    ana.add(LogPrior::Flat(p, "mass::b(MSbar)", ParameterRange{ -10, 10 }));
                    ana.add(LogPrior::Flat(p, "mass::c", ParameterRange{ -10, 10 }));
                    ana.restrict("mass::c", -10, 0);
                    PopulationMonteCarloSampler pmc_sampler(ana, hdf5::File::Open(mcmc_file_name), pmc_config);
                    pmc_sampler.run();
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
                    TEST_CHECK_RELATIVE_ERROR(1.0 / 20.0 / 20.0, std::get<2>(record_hc), 5e-3);

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

                        TEST_CHECK_EQUAL(record_resume0[3], std::get<0>(record_split)); // posterior
                        // todo For some reason, values don't agree bitwise. Why?
                        TEST_CHECK_RELATIVE_ERROR(record_resume0[4], std::get<1>(record_split), 2e-14); // weight
                    }
                }
            }
#endif
        }
} population_monte_carlo_sampler_test;
