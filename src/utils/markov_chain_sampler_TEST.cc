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
#include <src/utils/analysis_TEST.hh>
#include <src/utils/markov_chain_sampler.hh>

using namespace test;
using namespace eos;

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
            static const double eps = 1e-14;

            // check MarkovChainSampler::Config
            {
                MarkovChainSampler::Config conf = MarkovChainSampler::Config::Default();
                TEST_CHECK_THROWS(VerifiedRangeUnderflow, conf.min_efficiency = -0.1);

                MarkovChainSampler::Config conf_quick = MarkovChainSampler::Config::Quick();
                TEST_CHECK_THROWS(VerifiedRangeOverflow, conf.max_efficiency = 23.1);
            }

            // empty Analysis
            {
                AnalysisPtr analysis = std::make_shared<Analysis>(std::make_shared<LogLikelihood>(Parameters::Defaults()));
                MarkovChain chain(analysis, 13);

                TEST_CHECK_THROWS(InternalError, chain.run(2));
            }

            // check pre run, main run and HDF5 storage
            {
                static const std::string file_name(EOS_BUILDDIR "/src/utils/markov_chain_sampler_TEST.hdf5");
                std::remove(file_name.c_str());

                // store to HDF5
                {
                    Analysis analysis(make_analysis(true));

                    MarkovChainSampler::Config config = MarkovChainSampler::Config::Quick();
                    config.need_prerun = true;
                    config.number_of_chains = 3;
                    config.min_efficiency = 0.20;
                    config.rvalue_criterion_param = 1.1;
                    config.use_strict_rvalue_definition = true;
                    config.use_posterior_rvalue = true;
                    config.parallelize = true;
                    config.prerun_iterations_update = 500;
                    config.chunk_size = 100;
                    config.chunks = 5;
                    config.seed = 1346;
                    config.output_file.reset(new ScanFile(ScanFile::Create(file_name, "markov_chain_sampler_TEST")));

                    MarkovChainSampler sampler(analysis, config);
                    sampler.run();

                    MarkovChainSampler::PreRunInfo pre_info(sampler.pre_run_info());

                    TEST_CHECK(pre_info.converged);
                    TEST_CHECK_EQUAL(pre_info.iterations_at_convergence, pre_info.iterations);
                    TEST_CHECK_EQUAL(pre_info.iterations_at_convergence, 1500);
                    TEST_CHECK_NEARLY_EQUAL(pre_info.rvalue_parameters.front(), 1.007285609117614777, eps);
                    TEST_CHECK_RELATIVE_ERROR(pre_info.rvalue_posterior, 1.042859920455756262, eps);
                }

                // open HDF5 file and run checks on data
                {
                    // read file from disk
                    ScanFile file = ScanFile::Open(file_name);

                    // main run of third chain
                    ScanFile::DataSet data_set = file["chain #2"];

                    /*
                     * check in ipython using
                     *   import h5py
                     *   f = h5py.File("test.hdf5"); f["data"]['chain #2'].shape
                     */
                    TEST_CHECK_EQUAL(data_set.records(), 500);

                    // did we store one param + log(posterior)?
                    TEST_CHECK_EQUAL(data_set.fields(), 2);

                    auto f = data_set.begin_fields();

                    // proper parameter?
                    TEST_CHECK_EQUAL(f->name(), "mass::b(MSbar)");
                    TEST_CHECK_EQUAL(f->get("min", 0.0), 3.7);
                    TEST_CHECK_EQUAL(f->get("max", 0.0), 4.9);
                    TEST_CHECK_EQUAL(f->get("nuisance", 17.0), false);

                    ++f;

                    // only log(posterior)?
                    TEST_CHECK_EQUAL(f->name(), "posterior");
                    // posterior has currently no further attributes
                }
            }
        }
} markov_chain_sampler_test;
