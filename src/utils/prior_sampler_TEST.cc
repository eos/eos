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

#include <src/utils/prior_sampler.hh>
#include <src/utils/log_prior.hh>
#include <src/utils/analysis_TEST.hh>
#include <test/test.hh>

using namespace test;
using namespace eos;

class PriorSamplerTest :
    public TestCase
{
    public:
        PriorSamplerTest() :
            TestCase("observable_evaluator_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-15;

            // create simple PriorSampler
            {
                Parameters p = Parameters::Defaults();

                // load observables
                ObservableSet o;
                o.add(ObservablePtr(new TestObservable(p, Kinematics(), "mass::b(MSbar)")));
                o.add(ObservablePtr(new TestObservable(p, Kinematics(), "mass::c")));
                o.add(ObservablePtr(new TestObservable(p, Kinematics(), "mass::s")));

                PriorSampler::Config config = PriorSampler::Config::Default();
                config.chunks = 2;
                config.chunk_size = 2;
                config.seed = 1;
                config.parallelize = true;
                config.store_parameters = true;
                config.output_file.reset(new ScanFile(ScanFile::Create("/tmp/prior_sampler_test.hdf5", "prior_sampler_test")));

                PriorSampler sampler(o, config);

                sampler.add(LogPrior::Gauss(p, "mass::b(MSbar)", ParameterRange{3.5, 4.5}, 4.1, 4.2, 4.3));
                sampler.add(LogPrior::Gauss(p, "mass::c", ParameterRange{1, 2}, 1.1, 1.2, 1.3));
                sampler.add(LogPrior::Flat(p, "mass::s", ParameterRange{0, 0.1}));

                sampler.run();
            }

            // open HDF5 file and run checks on data
            {
                // read file from disk
                ScanFile file = ScanFile::Open("/tmp/prior_sampler_test.hdf5");

                // main run of third chain
                ScanFile::DataSet data_set = file["data"];

                TEST_CHECK_EQUAL(data_set.records(), 4);

                // did we store 3 param + 3 observables
                TEST_CHECK_EQUAL(data_set.fields(), 6);

                auto f = data_set.begin_fields();

                // proper parameter?
                TEST_CHECK_EQUAL(f->name(), "mass::b(MSbar)");
                TEST_CHECK_EQUAL(f->get("min", 0.0), 3.5);
                TEST_CHECK_EQUAL(f->get("max", 0.0), 4.5);
                TEST_CHECK_EQUAL(f->get("nuisance", 1.0), false);

                auto record = data_set[0];

                // parameters = observables
                TEST_CHECK_NEARLY_EQUAL(record[0], 4.178903953642595  , eps);
                TEST_CHECK_NEARLY_EQUAL(record[1], 1.477605489382485  , eps);
                TEST_CHECK_NEARLY_EQUAL(record[2], 0.07203244892880321, eps);
                TEST_CHECK_NEARLY_EQUAL(record[3], 4.178903953642595  , eps);
                TEST_CHECK_NEARLY_EQUAL(record[4], 1.477605489382485  , eps);
                TEST_CHECK_NEARLY_EQUAL(record[5], 0.07203244892880321, eps);
            }
        }
} prior_sampler_test;
