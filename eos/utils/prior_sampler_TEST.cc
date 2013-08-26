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

#include <eos/utils/prior_sampler.hh>
#include <eos/utils/analysis_TEST.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log_prior.hh>
#include <test/test.hh>

using namespace test;
using namespace eos;

class PriorSamplerTest :
    public TestCase
{
    public:
        PriorSamplerTest() :
            TestCase("prior_sampler_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-14;

            static const std::string file_name(EOS_BUILDDIR "/eos/utils/prior_sampler_TEST.hdf5");

            PriorSampler::Config config = PriorSampler::Config::Default();
            config.n_samples = 4;
            config.seed = 1;
            config.parallelize = true;
            config.store_parameters = true;
            config.output_file.reset(new hdf5::File(hdf5::File::Create(file_name)));

            const unsigned n_samples = config.n_samples;// create simple PriorSampler
            {
                Parameters p = Parameters::Defaults();

                // load observables
                ObservableSet o;
                o.add(ObservablePtr(new TestObservable(p, Kinematics(), "mass::b(MSbar)")));
                o.add(ObservablePtr(new TestObservable(p, Kinematics(), "mass::c")));
                o.add(ObservablePtr(new TestObservable(p, Kinematics(), "mass::s(2GeV)")));

                PriorSampler sampler(o, config);

                sampler.add(LogPrior::Gauss(p, "mass::b(MSbar)", ParameterRange{3.5, 4.5}, 4.1, 4.2, 4.3));
                sampler.add(LogPrior::Gauss(p, "mass::c", ParameterRange{1, 2}, 1.1, 1.2, 1.3));
                sampler.add(LogPrior::Flat(p, "mass::s(2GeV)", ParameterRange{0, 0.1}));

                sampler.run();
            }

            // open HDF5 file and run checks on data
            {
                // read file from disk
                auto file = hdf5::File::Open(file_name);

                auto data_obs = file.open_data_set("/data/observables", PriorSampler::observables_type(3));

                TEST_CHECK_EQUAL(data_obs.records(), n_samples);

                std::vector<double> obs_record(3);
                data_obs >> obs_record;

                // parameters == observables
                TEST_CHECK_NEARLY_EQUAL(obs_record[0], 4.17890395364246,   eps);
                TEST_CHECK_NEARLY_EQUAL(obs_record[1], 1.47760548938249,   eps);
                TEST_CHECK_NEARLY_EQUAL(obs_record[2], 0.0720324489288032, eps);

                hdf5::Composite<hdf5::Scalar<double>, hdf5::Scalar<double>> par_type
                {
                    "parameter description",
                    hdf5::Scalar<double>("min"),
                    hdf5::Scalar<double>("max"),
                };
                auto data_par = file.open_data_set("/descriptions/parameters/0", par_type);
                auto par_record = std::make_tuple(0.0, 0.0);
                data_par >> par_record;

                // proper parameter?
                TEST_CHECK_EQUAL(std::get<0>(par_record), 3.5);
                TEST_CHECK_EQUAL(std::get<1>(par_record), 4.5);
            }
        }
} prior_sampler_test;
