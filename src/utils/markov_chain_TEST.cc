/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

#include <config.h>
#include <src/utils/markov_chain.hh>
#include <src/utils/analysis_TEST.hh>
#include <test/test.hh>

using namespace test;
using namespace eos;

class MarkovChainTest :
    public TestCase
{
    public:
        MarkovChainTest() :
            TestCase("markov_chain_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-14;
            Analysis analysis = make_analysis(false);

            // check that clones are independent
            {
                auto clone1 = analysis.clone();
                auto clone2 = analysis.clone();

                // create two chains and compare results
                MarkovChain chain1(clone1, 13);
                MarkovChain chain2(clone2, 13134);

                double mB1_before(chain1.parameter_descriptions().front().parameter);
                double mB2_before(chain2.parameter_descriptions().front().parameter);

                // running chain 1 shouldn't affect chain 2
                // chain1 should have some moves accepted
                chain1.run(300);
                if (chain1.parameter_descriptions().front().parameter() == mB1_before)
                    TEST_CHECK_FAILED("chain1 did not move");

                // running chain1 shouldn't affect chain2
                TEST_CHECK_EQUAL(mB2_before, chain2.parameter_descriptions().front().parameter());
            }

            // step by step analysis of proposed moves. Out of range, too improbable, accepted... all in here
            // assume scale = 1.0
            {
                auto clone = analysis.clone();
                MarkovChain chain(clone, 1313);

                //initial points
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().point.front(),4.8910029008984566, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().log_likelihood,-22.490603892714894, eps);

                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().log_prior,-10.6705455880927786, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().log_posterior,-33.1611494808076728, eps);

                // next point is not accepted, out-of-range
                chain.run(1);
                TEST_CHECK_EQUAL(chain.iterations_last_run(), 1);
                TEST_CHECK(chain.proposal_accepted() ==false);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().point.front(),4.8910029008984566, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_proposal().log_posterior,-33.1611494808076728, eps); //not evaluated if proposal out of range


                // again out of range
                chain.run(1);
                TEST_CHECK(chain.proposal_accepted() ==false);
                TEST_CHECK_EQUAL(chain.iterations_last_run(), 1);
                TEST_CHECK_EQUAL(chain.statistics().iterations_total,2);

                // now we accept it
                chain.run(1);
                TEST_CHECK(chain.proposal_accepted() ==true);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().point.front(),4.2610145541614974, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().log_likelihood,1.1975077688130611, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().log_prior,0.417799138695899186, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().log_posterior,1.61530690750896033, eps); //reevaluated

                // in range, but not accepted
                chain.run(1);
                TEST_CHECK(chain.proposal_accepted() ==false);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().point.front(),4.2610145541614974, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().log_likelihood,1.1975077688130611, eps);


                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().log_prior,0.417799138695899186, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().log_posterior,1.61530690750896033, eps);

                TEST_CHECK_RELATIVE_ERROR(chain.info_at_proposal().point.front(),4.5447255790251928, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_proposal().log_likelihood,-4.5581396819233833, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_proposal().log_prior,0.336372185233400678, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_proposal().log_posterior,-4.2217674966899823, eps);

                // jump to next accepted point.
                chain.run(53);
                TEST_CHECK(chain.proposal_accepted() ==true);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().point.front(),4.3827129748249529, eps);
            }

            // this verifies that chain.run(1) 57 times will end at same point as chain.run(57)
            // it is the second accepted move
            {
                auto clone = analysis.clone();
                MarkovChain chain1(clone, 1313);

                unsigned N = 57;
                for(unsigned i=0; i<N;++i)
                    chain1.run(1);

                auto clone2 = analysis.clone();
                MarkovChain chain2(clone2, 1313);

                chain2.run(N);

                TEST_CHECK_EQUAL(chain1.info_at_current().point.front(), chain2.info_at_current().point.front());
                TEST_CHECK_EQUAL(chain1.statistics().iterations_total, N);
                TEST_CHECK_EQUAL(chain2.statistics().iterations_total, N);
                TEST_CHECK_EQUAL(chain1.statistics().mean_of_parameters.front(), chain2.statistics().mean_of_parameters.front());
                TEST_CHECK(chain1.statistics().iterations_accepted.front() != chain2.statistics().iterations_accepted.front()); //efficiencies reset in each call to run
            }

            // same as before, but now check statistics calculations
            // used open-office Welford.ods to check the  values by hand
            {
                auto clone = analysis.clone();
                MarkovChain chain(clone, 1313);

                chain.run(1);
                TEST_CHECK_EQUAL(chain.statistics().mean_of_parameters.front(), chain.info_at_current().point.front()); //mean of one iteration is just the current point
                TEST_CHECK_EQUAL(chain.statistics().variance_of_parameters.front(),0);
                TEST_CHECK_EQUAL(chain.statistics().mean_of_posterior, chain.info_at_current().log_posterior);
                TEST_CHECK_EQUAL(chain.statistics().variance_of_posterior, 0);

                chain.run(1);
                TEST_CHECK_EQUAL(chain.statistics().mean_of_parameters.front(), chain.info_at_current().point.front()); //mean of one iteration is just the current point
                TEST_CHECK_EQUAL(chain.statistics().variance_of_parameters.front(),0);
                TEST_CHECK_EQUAL(chain.statistics().mean_of_posterior, chain.info_at_current().log_posterior);
                TEST_CHECK_EQUAL(chain.statistics().variance_of_posterior, 0);

                chain.run(1);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mean_of_parameters.front(), 4.6810067853194699, eps); //mean of one iteration is just the current point
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_parameters.front(), 0.13229510567478892, eps);
                TEST_CHECK_RELATIVE_ERROR (chain.statistics().mean_of_posterior, -21.568997351368793, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_posterior, 403.13397297616291, eps);

                chain.run(1);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mean_of_parameters.front(), 4.5760087275299766, eps); //mean of one iteration is just the current point
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_parameters.front(), 0.13229510567478892, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mean_of_posterior, -15.772921286649355, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_posterior, 403.13397297616081, eps);

                // jump to next accepted move.  But first had higher posterior
                chain.run(53);
                TEST_CHECK_RELATIVE_ERROR(chain.info_at_current().point.front(), 4.3827129748249529, eps);

                TEST_CHECK_RELATIVE_ERROR(chain.statistics().parameters_at_mode.front(), 4.2610145541614974, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mode_of_posterior, 1.61530690750896033, eps);
                TEST_CHECK(chain.info_at_current().log_posterior != chain.statistics().mode_of_posterior);

            }

            // test statistics: mean and variance. Should be around 4.2 and 0.1^2 for many iterations
            {
                auto analysis = std::make_shared<Analysis>(make_analysis(true));
                MarkovChain chain(analysis, 13);
                chain.set_scale(0, 0.2);
                chain.run(6);

                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mean_of_parameters.front(),     +4.5698409084012441,   eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_parameters.front(), +0.024118869364409479, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().parameters_at_mode.front(),     +4.2528309889691,      eps);

                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mean_of_posterior,              -6.6427427635412242,   eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_posterior,         +14.246281595803778,    eps);
            }

            // check that efficiency is correctly recorded
            // first accepted move after 3 iterations, second after 57 iterations
            {
                auto clone1 = analysis.clone();
                MarkovChain chain1(clone1, 1313);

                auto clone2 = analysis.clone();
                MarkovChain chain2(clone2, 1313);

                chain1.run(10);
                TEST_CHECK_EQUAL(chain1.statistics().iterations_total,10);
                TEST_CHECK_EQUAL(chain1.statistics().iterations_accepted.front(),1);
                TEST_CHECK_EQUAL(chain1.statistics().iterations_rejected.front(),9);

                chain2.run(57);
                TEST_CHECK_EQUAL(chain2.statistics().iterations_total, 57);
                TEST_CHECK_EQUAL(chain2.statistics().iterations_accepted.front(),2);
                TEST_CHECK_EQUAL(chain2.statistics().iterations_rejected.front(),55);
            }

            // check 3D as well
            {
                auto clone = analysis.clone();
                clone->add(LogPrior::Flat(clone->parameters(), "Abs{c9}", ParameterRange{0.5, 4}), false);
                clone->add(LogPrior::Flat(clone->parameters(), "mass::e", ParameterRange{5.1e-4, 5.11e-4}), true);
                MarkovChain chain(clone, 1313);

                // check initial
                TEST_CHECK_EQUAL(chain.info_at_current().point, chain.info_at_proposal().point);
                TEST_CHECK_EQUAL(chain.info_at_current().log_posterior, chain.info_at_proposal().log_posterior);

                chain.run(45);
                chain.run(37);
                TEST_CHECK_EQUAL(chain.statistics().iterations_total, 45+37);

                // if last move accepted, then current and proposal should agree
                TEST_CHECK(chain.proposal_accepted());
                TEST_CHECK_EQUAL(chain.info_at_current().point, chain.info_at_proposal().point);
                TEST_CHECK_EQUAL(chain.info_at_current().log_posterior, chain.info_at_proposal().log_posterior);
                TEST_CHECK_EQUAL(chain.statistics().iterations_accepted.back(), 10);

                // if last rejected, current and proposal should disagree, but only in last parameter
                // but posterior is arbitrary, since proposal could be out-of-range
                // and thus posterior not evaluated/updated
                chain.run(2);
                TEST_CHECK(! chain.proposal_accepted());
                TEST_CHECK_EQUAL(chain.statistics().iterations_accepted.back(), 1);
                TEST_CHECK(chain.info_at_current().point != chain.info_at_proposal().point);
                TEST_CHECK_EQUAL(chain.info_at_current().point[0], chain.info_at_proposal().point[0]);
                TEST_CHECK_EQUAL(chain.info_at_current().point[1], chain.info_at_proposal().point[1]);
            }

            // check likelihood short cut: define an irrelevant parameter first
            {
                Parameters parameters = Parameters::Defaults();

                LogLikelihoodPtr llh(new LogLikelihood(parameters));
                llh->add(ObservablePtr(new TestObservable(parameters, Kinematics(),
                            "mass::b(MSbar)")), 4.1, 4.2, 4.3);

                AnalysisPtr ana = std::make_shared<Analysis>(llh);
                // the irrelevant parameter
                ana->add(LogPrior::Flat(parameters, "CKM::etabar", ParameterRange{ 3.7, 4.9 }));
                // an interesting parameter that affects the observable
                ana->add(LogPrior::Flat(parameters, "mass::b(MSbar)", ParameterRange{ 3.7, 4.9 }));
                // another irrelevant parameters
                ana->add(LogPrior::Flat(parameters, "B->K^*::a1_uncertainty@BZ2004", ParameterRange{ 3.7, 4.9 }));

                MarkovChain chain(ana, 12345);

                // if this evaluates to NaN, initialization went wrong
                TEST_CHECK_RELATIVE_ERROR( (*llh)(), -14.758100406210971, eps);
            }

            // try sampling with discrete parameters only
            {
                Parameters parameters = Parameters::Defaults();
                LogLikelihoodPtr llh(new LogLikelihood(parameters));
                llh->add(ObservablePtr(new TestObservable(parameters, Kinematics(), "mass::b(MSbar)")), 4.1, 4.2, 4.3);

                AnalysisPtr ana = std::make_shared<Analysis>(llh);
                ana->add(LogPrior::Discrete(parameters, "mass::b(MSbar)", std::set<double>{ 4.15, 4.4 }));

                MarkovChain chain(ana, 12346);

                unsigned samples = 1000;
                chain.run(samples);

                static const std::string file_name(EOS_BUILDDIR "/src/utils/markov_chain_TEST_discrete.hdf5");
                std::remove(file_name.c_str());
                std::shared_ptr<ScanFile> file(new ScanFile(ScanFile::Create(file_name, "markov_chain_TEST")));
                std::shared_ptr<ScanFile::DataSet> data_set(new ScanFile::DataSet(file->add("chain #0", 1 + 1)));

                chain.dump_history(*data_set);

                // make sure there are only two values
                unsigned counter = 0;
                for (unsigned i = 0 ; i < samples ; ++i)
                {
                    double param = (*data_set)[i][0];

                    if (std::abs(param - 4.15) < eps)
                    {
                        ++counter;
                        continue;
                    }

                    if (std::abs(param - 4.4) < eps)
                        continue;

                    TEST_CHECK_FAILED("Parameter in step "+ stringify(i)+" of value "+stringify(param)+", which is neither 4.15 nor 4.4!");
                }

                // pdf at 4.15 is about 6.5208 higher than at 4.4
                // but accuracy is low for 1000 samples
                TEST_CHECK_NEARLY_EQUAL(1.0 * counter / (samples - counter), 6.5208, 2e-1);

                // remove the HDF5 file
                std::remove(file_name.c_str());
            }
        }
} markov_chain_test;
