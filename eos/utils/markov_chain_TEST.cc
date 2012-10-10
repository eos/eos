/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

#include <config.h>
#include <eos/utils/analysis_TEST.hh>
#include <eos/utils/markov_chain.hh>
#include <eos/utils/proposal_functions.hh>
#include <test/test.hh>

#include <algorithm>

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

            // empty Analysis
            TEST_SECTION("empty-analysis",
            {
                Analysis analysis(LogLikelihood(Parameters::Defaults()));
                std::shared_ptr<MarkovChain::ProposalFunction> ppf(new proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.01 }));
                TEST_CHECK_THROWS(InternalError, MarkovChain chain(analysis, 13, ppf));
            });

            // check if empty proposal function throws InternalError
            TEST_SECTION("empty-proposal-function",
            {
                TEST_CHECK_THROWS(InternalError, MarkovChain chain(analysis, 13, std::shared_ptr<MarkovChain::ProposalFunction>()));
            });

            // check that clones are independent
            TEST_SECTION("clones",
            {
                std::shared_ptr<MarkovChain::ProposalFunction> ppf1(new proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.01 }));
                std::shared_ptr<MarkovChain::ProposalFunction> ppf2(new proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.01 }));

                // create two chains and compare results
                MarkovChain chain1(analysis, 13, ppf1);
                MarkovChain chain2(analysis, 13134, ppf2);

                double mB1_before(chain1.parameter_descriptions().front().parameter);
                double mB2_before(chain2.parameter_descriptions().front().parameter);

                // running chain 1 shouldn't affect chain 2
                // chain1 should have some moves accepted
                chain1.run(300);
                if (chain1.parameter_descriptions().front().parameter() == mB1_before)
                    TEST_CHECK_FAILED("chain1 did not move");

                // running chain1 shouldn't affect chain2
                TEST_CHECK_EQUAL(mB2_before, chain2.parameter_descriptions().front().parameter());
            });
#if 0
            // step by step analysis of proposed moves. Out of range, too improbable, accepted... all in here
            TEST_SECTION("step-by-step",
            {
                std::shared_ptr<MarkovChain::ProposalFunction> ppf(new proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.03 }));
                MarkovChain chain(analysis, 1313, ppf);

                std::cout << "initial:" << std::endl;
                std::cout << chain.current_state().point[0] << std::endl;
                std::cout << chain.current_state().log_likelihood << std::endl;
                std::cout << chain.current_state().log_prior << std::endl;
                std::cout << chain.current_state().log_posterior << std::endl;

                //initial points
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().point.front(),    4.8910029008984566, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_likelihood, -22.4906038927148940, eps);

                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_prior,      -10.6705455880927786, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_posterior,  -33.1611494808076728, eps);

                // next point is not accepted, out-of-range
                chain.run(1);
                std::cout << chain.current_state().point[0] << std::endl;
                std::cout << chain.current_state().log_likelihood << std::endl;
                std::cout << chain.current_state().log_prior << std::endl;
                std::cout << chain.current_state().log_posterior << std::endl;
                TEST_CHECK_EQUAL(chain.iterations_last_run(), 1);
                TEST_CHECK_EQUAL(chain.proposal_accepted(), false);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().point.front(),    4.8910029008984566, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.proposed_state().log_posterior, -33.1611494808076728, eps); //not evaluated if proposal out of range

                // again out of range
                chain.run(1);
                std::cout << chain.current_state().point[0] << std::endl;
                std::cout << chain.current_state().log_likelihood << std::endl;
                std::cout << chain.current_state().log_prior << std::endl;
                std::cout << chain.current_state().log_posterior << std::endl;
                TEST_CHECK_EQUAL(chain.proposal_accepted(), false);
                TEST_CHECK_EQUAL(chain.iterations_last_run(), 1);
                TEST_CHECK_EQUAL(chain.statistics().iterations_total,2);

                // now we accept it
                chain.run(1);
                std::cout.precision(16);
                std::cout << chain.current_state().point[0] << std::endl;
                std::cout << chain.current_state().log_likelihood << std::endl;
                std::cout << chain.current_state().log_prior << std::endl;
                std::cout << chain.current_state().log_posterior << std::endl;
                TEST_CHECK(chain.proposal_accepted() ==true);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().point.front(),    4.872512050011867, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_likelihood, -21.229976310768970, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_prior,       -9.779735023878651, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_posterior,   -31.00971133464762, eps); // reevaluated

                // in range, but not accepted
                chain.run(1);
                TEST_CHECK_EQUAL(chain.proposal_accepted(), false);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().point.front(),    4.872512050011867, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_likelihood, -21.229976310768970, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_prior,       -9.779735023878651, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_posterior,   -31.00971133464762, eps); // not reevaluated

                chain.run(1);
                TEST_CHECK_EQUAL(chain.proposal_accepted(), true);
                TEST_CHECK_RELATIVE_ERROR(chain.proposed_state().point.front(),    4.812810423973447, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.proposed_state().log_likelihood, -17.393184226736560, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.proposed_state().log_prior,       -7.136975460614643, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.proposed_state().log_posterior,  -24.530159687351200, eps);

                chain.run(45);
                TEST_CHECK_EQUAL(chain.statistics().iterations_total, 50); // 5x1 + 1x45
                TEST_CHECK_EQUAL(chain.statistics().iterations_accepted, 21); // 21 acceptions in the last 45 iterations
                TEST_CHECK_EQUAL(chain.statistics().iterations_rejected, 24); // 24 acceptions in the last 45 iterations
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().point[0],         4.3920766161687550, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_likelihood,  -0.4610247641525864, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_prior,        1.3805078458753980, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_posterior,    0.9194830817228119, eps);
            });

            // this verifies that chain.run(1) 50 times will end at same point as chain.run(50)
            TEST_SECTION("iterations",
            {
                std::shared_ptr<MarkovChain::ProposalFunction> ppf1(new proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.01 }));
                std::shared_ptr<MarkovChain::ProposalFunction> ppf2(new proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.01 }));
                MarkovChain chain1(analysis, 1313, ppf1);
                MarkovChain chain2(analysis, 1313, ppf2);

                static const unsigned N = 57;
                for (unsigned i = 0 ; i < N; ++i)
                {
                    chain1.run(1);
                }
                chain2.run(N);

                TEST_CHECK_EQUAL(chain1.current_state().point.front(), chain2.current_state().point.front());
                TEST_CHECK_EQUAL(chain1.statistics().iterations_total, N);
                TEST_CHECK_EQUAL(chain2.statistics().iterations_total, N);
                TEST_CHECK_EQUAL(chain1.statistics().mean_of_parameters.front(), chain2.statistics().mean_of_parameters.front());
                TEST_CHECK(chain1.statistics().iterations_accepted != chain2.statistics().iterations_accepted); //efficiencies reset in each call to run
            });

            // same as before, but now check statistics calculations
            // used open-office spreadsheet Welford.ods to check the values by hand
            TEST_SECTION("statistics",
            {
                std::shared_ptr<MarkovChain::ProposalFunction> ppf(new proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.05 }));
                MarkovChain chain(analysis, 1313, ppf);

                chain.run(1);
                TEST_CHECK_EQUAL(chain.statistics().mean_of_parameters.front(), chain.current_state().point.front()); // mean of one iteration is just the current point
                TEST_CHECK_EQUAL(chain.statistics().variance_of_parameters.front(),0);
                TEST_CHECK_EQUAL(chain.statistics().mean_of_posterior, chain.current_state().log_posterior);
                TEST_CHECK_EQUAL(chain.statistics().variance_of_posterior, 0);

                chain.run(1);
                TEST_CHECK_EQUAL(chain.statistics().mean_of_parameters.front(), chain.current_state().point.front()); // 2nd iteration has no accepted points.
                TEST_CHECK_EQUAL(chain.statistics().variance_of_parameters.front(),0);
                TEST_CHECK_EQUAL(chain.statistics().mean_of_posterior, chain.current_state().log_posterior);
                TEST_CHECK_EQUAL(chain.statistics().variance_of_posterior, 0);

                chain.run(1);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mean_of_parameters[0],       4.8830457056161270000, eps); // mean of first two accepted points
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_parameters[0],   0.0001899508702833833, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mean_of_posterior,         -32.2395994688615600000, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_posterior,       2.5477632735536560000, eps);

                chain.run(7);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mean_of_parameters[0],     4.804728197447452000, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_parameters[0], 0.005006338824707367, eps);
            });

            // test statistics: mean and variance. Should be around 4.2 and 0.1^2 for many iterations
            TEST_SECTION("statistics2",
            {
                std::shared_ptr<MarkovChain::ProposalFunction> ppf(new proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.01 }));
                MarkovChain chain(make_analysis(true), 13, ppf);

                chain.run(1000);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mean_of_parameters.front(),     4.18292971650058200, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_parameters.front(), 0.01072836817339981, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().parameters_at_mode.front(),     4.19990405758096100, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().mean_of_posterior,              0.65087328379658350, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.statistics().variance_of_posterior,          0.57268622365364590, eps);
            });

            // check that efficiency is correctly recorded
            TEST_SECTION("efficiency",
            {
                std::shared_ptr<MarkovChain::ProposalFunction> ppf(new proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.03 }));

                MarkovChain chain(analysis, 1313, ppf);

                chain.run(50);
                TEST_CHECK_EQUAL(chain.statistics().iterations_total,    50);
                TEST_CHECK_EQUAL(chain.statistics().iterations_accepted, 23);
                TEST_CHECK_EQUAL(chain.statistics().iterations_rejected, 27);
            });

            // check 3D as well
            TEST_SECTION("3d",
            {
                auto clone = analysis.clone();
                clone->add(LogPrior::Flat(clone->parameters(), "Abs{c9}", ParameterRange{0.5, 4}), false);
                clone->add(LogPrior::Flat(clone->parameters(), "mass::e", ParameterRange{5.1e-4, 5.11e-4}), true);
                std::shared_ptr<MarkovChain::ProposalFunction> ppf(new proposal_functions::MultivariateGaussian(3,
                            std::vector<double>
                            {
                                0.01, 0.0,   0.0,
                                0.0,  0.025, 0.0,
                                0.0,  0.0,   2.5e-11
                            }));
                MarkovChain chain(*clone, 1313, ppf);

                // check initial
                TEST_CHECK_EQUAL(chain.current_state().point, chain.proposed_state().point);
                TEST_CHECK_EQUAL(chain.current_state().log_posterior, chain.proposed_state().log_posterior);

                chain.run(81);
                TEST_CHECK_EQUAL(chain.statistics().iterations_total, 45+36);

                // if last move is accepted, then current and proposal should agree
                TEST_CHECK(chain.proposal_accepted());
                TEST_CHECK_EQUAL(chain.current_state().point, chain.proposed_state().point);
                TEST_CHECK_EQUAL(chain.current_state().log_posterior, chain.proposed_state().log_posterior);
                TEST_CHECK_EQUAL(chain.statistics().iterations_accepted, 10);

                // if last move is rejected, current and proposal should disagree
                // posterior is arbitrary, since proposal could be out-of-range
                // and thus posterior not evaluated/updated
                chain.run(1);
                TEST_CHECK(! chain.proposal_accepted());
                TEST_CHECK_EQUAL(chain.statistics().iterations_accepted, 0);
                TEST_CHECK(chain.current_state().point != chain.proposed_state().point);
            });
#endif
#if 0
            // check likelihood short cut: define an irrelevant parameter first
            {
                Parameters parameters = Parameters::Defaults();

                LogLikelihood llh(parameters);
                llh.add(ObservablePtr(new TestObservable(parameters, Kinematics(),
                            "mass::b(MSbar)")), 4.1, 4.2, 4.3);

                Analysis analysis(llh);
                // the irrelevant parameter
                analysis.add(LogPrior::Flat(parameters, "CKM::etabar", ParameterRange{ 3.7, 4.9 }));
                // an interesting parameter that affects the observable
                analysis.add(LogPrior::Flat(parameters, "mass::b(MSbar)", ParameterRange{ 3.7, 4.9 }));
                // another irrelevant parameters
                analysis.add(LogPrior::Flat(parameters, "B->K^*::a1_uncertainty@BZ2004", ParameterRange{ 3.7, 4.9 }));

                std::shared_ptr<MarkovChain::ProposalFunction> ppf(new proposal_functions::MultivariateGaussian(3,
                            std::vector<double>
                            {
                                0.10, 0.0,   0.0,
                                0.0,  0.09,  0.0,
                                0.0,  0.0,   0.08
                            }));
                MarkovChain chain(analysis, 12345, ppf);

                // if this evaluates to NaN, initialization went wrong
                TEST_CHECK_RELATIVE_ERROR(llh(), -14.758100406210971, eps);
            }
#endif
#if 0
            // try sampling with discrete parameters only
            {
                Parameters parameters = Parameters::Defaults();
                LogLikelihood llh(parameters);
                llh.add(ObservablePtr(new TestObservable(parameters, Kinematics(), "mass::b(MSbar)")), 4.1, 4.2, 4.3);

                Analysis analysis(llh);
                analysis->add(LogPrior::Discrete(parameters, "mass::b(MSbar)", std::set<double>{ 4.15, 4.4 }));

                MarkovChain chain(ana, 12346);

                unsigned samples = 1000;
                chain.run(samples);

                static const std::string file_name(EOS_BUILDDIR "/eos/utils/markov_chain_TEST_discrete.hdf5");
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
#endif

            // changing the point of a chain by hand
            TEST_SECTION("set-point",
            {
                std::shared_ptr<MarkovChain::ProposalFunction> ppf(new proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.01 }));
                auto analysis = make_analysis(false);
                MarkovChain chain(analysis, 13, ppf);

                MarkovChain::HyperParameter h { 0 } ;
                chain.set_point(std::vector<double>{ 4.3}, h);
                TEST_CHECK_EQUAL(chain.current_state().point.front(), 4.3);

                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_likelihood, 0.88364655978937656, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_prior, 0.883646846442260436, eps);
                TEST_CHECK_RELATIVE_ERROR(chain.current_state().log_posterior,  0.88364655978937656 + 0.883646846442260436, eps);
            });
            TEST_SECTION("Multivariate::adapt",
            {
                Parameters parameters = Parameters::Defaults();

                LogLikelihood llh(parameters);
                llh.add(ObservablePtr(new TestObservable(parameters, Kinematics(),
                            "mass::b(MSbar)")), 4.1, 4.2, 4.3);

                Analysis analysis(llh);
                // the irrelevant parameter
                analysis.add(LogPrior::Flat(parameters, "mass::c", ParameterRange{ 1.2, 1.34 }));
                // an interesting parameter that affects the observable
                analysis.add(LogPrior::Flat(parameters, "mass::b(MSbar)", ParameterRange{ 3.7, 4.9 }));

                // _covariance with zero correlation
                std::vector<double> cov(4, 0.0);
                cov[0] = 0.0049;
                cov[3] = 0.01;

                auto mvg = new proposal_functions::MultivariateGaussian(2, cov);
                mvg->covariance_scale = 2.38 * 2.38 / 2.0;

                std::shared_ptr<MarkovChain::ProposalFunction> ppf(mvg);
                MarkovChain chain(analysis, 12345, ppf);

                chain.run(5e4);
                double efficiency = 0.24;
                std::cout << "Efficiency " << efficiency << std::endl;
                ppf->adapt(chain.history().states.cbegin(), chain.history().states.cend(), efficiency, 0.2, 0.35);

                //cloning
                auto prop = mvg->clone();
                auto mvg2 = static_cast<proposal_functions::MultivariateGaussian *>(prop.get());

                // copy correct
                TEST_CHECK_EQUAL(gsl_matrix_get(mvg->covariance(), 0, 0), gsl_matrix_get(mvg2->covariance(), 0, 0));
                TEST_CHECK_EQUAL(gsl_matrix_get(mvg->covariance(), 0, 1), gsl_matrix_get(mvg2->covariance(), 0, 1));
                TEST_CHECK_EQUAL(gsl_matrix_get(mvg->covariance(), 1, 1), gsl_matrix_get(mvg2->covariance(), 1, 1));

                chain.run(500);
                ppf->adapt(chain.history().states.cbegin(), chain.history().states.cend(), efficiency, 0.2, 0.35);
                TEST_CHECK(gsl_matrix_get(mvg->covariance(), 0, 0) != gsl_matrix_get(mvg2->covariance(), 0, 0));

                // check if internal sample covariance() copied correctly as well
                mvg2->adapt(chain.history().states.cbegin(), chain.history().states.cend(), efficiency, 0.2, 0.35);
                TEST_CHECK_EQUAL(gsl_matrix_get(mvg->covariance(), 0, 0), gsl_matrix_get(mvg2->covariance(), 0, 0));
                TEST_CHECK_EQUAL(gsl_matrix_get(mvg->covariance(), 0, 1), gsl_matrix_get(mvg2->covariance(), 0, 1));
                TEST_CHECK_EQUAL(gsl_matrix_get(mvg->covariance(), 1, 1), gsl_matrix_get(mvg2->covariance(), 1, 1));
            });

            // adaptation with full correlation
            {
                Parameters p = Parameters::Defaults();
                Kinematics k;
                std::array<ObservablePtr, 2> obs;
                obs[0] = ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)"));
                obs[1] = ObservablePtr(new TestObservable(p, k, "mass::c"));

                ObservableCache cache(p);

                // multivariate gaussian with full (\rho = 1) correlation
                std::array<double, 2> mean{{ 4.3, 1.1 }};
                std::array<std::array<double, 2>, 2> _covariance;
                _covariance[0][0] = 0.1 * 0.1;
                _covariance[1][1] = 0.05 * 0.05;
                _covariance[0][1] = _covariance[1][0] = 0.0048;

                auto block = LogLikelihoodBlock::MultivariateGaussian<2>(cache, obs, mean, _covariance);

                LogLikelihood llh(p);
                llh.add(Constraint("Correlated Gaussian", std::vector<ObservablePtr>(obs.begin(), obs.end()), { block }));

                Analysis analysis(llh);

                // an interesting parameter that affects the observable
                analysis.add(LogPrior::Flat(p, "mass::b(MSbar)", ParameterRange{ 3.7, 4.9 }));

                // the 2nd parameter
                // NOTE: it is crucial for the accuracy of the _covariance estimate that the range be large enough
                //       to cover the gaussian peak out to many sigmas
                analysis.add(LogPrior::Flat(p, "mass::c", ParameterRange{ 0.7, 1.34 }));

                //start with uncorrelated gaussian
                std::vector<double> cov_initial
                {
                    0.01,   0.0045,
                    0.0045, 0.0025
                };

                const double scale = 2.38 * 2.38 / 2.0;
                bool automatic_scaling = true;
                TEST_SECTION("adapt correlated gaussian",
                {
                    proposal_functions::MultivariateGaussian* mvg = new proposal_functions::MultivariateGaussian(2, cov_initial, automatic_scaling);
                    std::shared_ptr<MarkovChain::ProposalFunction> ppf(mvg);
                    MarkovChain chain(analysis, 12345, ppf);

                    chain.run(1e4);

                    double efficiency = 1.0 * chain.statistics().iterations_accepted / (chain.statistics().iterations_accepted + chain.statistics().iterations_rejected);
                    ppf->adapt(chain.history().states.cbegin(), chain.history().states.cend(), efficiency, 0.2, 0.35);

                    // we get input _covariance up to ca. 10%, but with the usual scaling factor
                    double high_eps = 7e-2;
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvg->covariance(), 0, 0) / scale, 0.0100, high_eps);
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvg->covariance(), 0, 1) / scale, 0.0048, high_eps);
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvg->covariance(), 1, 0) / scale, 0.0048, high_eps);
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvg->covariance(), 1, 1) / scale, 0.0025, high_eps);

                    chain.run(5e4);
                    efficiency = 1.0 * chain.statistics().iterations_accepted / (chain.statistics().iterations_accepted + chain.statistics().iterations_rejected) ;
                    ppf->adapt(chain.history().states.cbegin(), chain.history().states.cend(), efficiency, 0.2, 0.35);
                    // sample _covariance is estimated more accurately after second step
                });

                TEST_SECTION("adapt correlated student t",
                {
                    double dof = 8;
                    proposal_functions::MultivariateStudentT * mvt = new proposal_functions::MultivariateStudentT(2, cov_initial, dof);

                    std::shared_ptr<MarkovChain::ProposalFunction> ppf(mvt);
                    MarkovChain chain(analysis, 12345, ppf);

                    chain.run(12e4);

                    double efficiency = 1.0 * chain.statistics().iterations_accepted / (chain.statistics().iterations_accepted + chain.statistics().iterations_rejected);

                    ppf->adapt(chain.history().states.cbegin(), chain.history().states.cend(), efficiency, 0.2, 0.35);

                    chain.run(3e4);
                    ppf->adapt(chain.history().states.cbegin(), chain.history().states.cend(), efficiency, 0.2, 0.35);

                    // with student-t, get same efficiency and more accurate result in less than half the trials
                    double low_eps = 4e-2;
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvt->covariance(), 0, 0) / scale, 0.0100, low_eps);
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvt->covariance(), 0, 1) / scale, 0.0048, low_eps);
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvt->covariance(), 1, 0) / scale, 0.0048, low_eps);
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvt->covariance(), 1, 1) / scale, 0.0025, low_eps);
                });

                TEST_SECTION("adapt correlated student t in gaussian limit",
                {
                    double dof = 1000;
                    proposal_functions::MultivariateStudentT * mvt = new proposal_functions::MultivariateStudentT(2, cov_initial, dof);
                    proposal_functions::MultivariateGaussian * mvg = new proposal_functions::MultivariateGaussian(2, cov_initial);

                    MarkovChain::State current;
                    current.point.push_back(4.3);
                    current.point.push_back(1.1);

                    // some extreme values
                    MarkovChain::State proposal;
                    proposal.point.push_back(4.1);
                    proposal.point.push_back(1.26);

                    TEST_CHECK_RELATIVE_ERROR(mvt->evaluate(current, proposal), mvg->evaluate(current, proposal), 3e-2);

                    std::shared_ptr<MarkovChain::ProposalFunction> ppf(mvt);
                    MarkovChain chain(analysis, 12345, ppf);

                    chain.run(3e4);

                    // hope that with this efficiency, scale is not changed automatically
                    double efficiency = 0.25;
                    ppf->adapt(chain.history().states.cbegin(), chain.history().states.cend(), efficiency, 0.2, 0.35);

                    // with student-t, get same efficiency
                    double low_eps = 4e-2;
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvt->covariance(), 0, 0) / scale, 0.0100, low_eps);
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvt->covariance(), 0, 1) / scale, 0.0048, low_eps);
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvt->covariance(), 1, 0) / scale, 0.0048, low_eps);
                    TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(mvt->covariance(), 1, 1) / scale, 0.0025, low_eps);

                    // now compare with Gaussian distribution
                    // both have same _covariance matrix
                    mvg->adapt(chain.history().states.cbegin(), chain.history().states.cend(), efficiency, 0.2, 0.35);

                    TEST_CHECK_RELATIVE_ERROR(mvt->evaluate(chain.current_state(), chain.proposed_state()),
                                              mvg->evaluate(chain.current_state(), chain.proposed_state()), 2e-2);
                });

//                TEST_SECTION("integration",
                {
                    proposal_functions::MultivariateGaussian* mvg = new proposal_functions::MultivariateGaussian(2, cov_initial);
                    mvg->covariance_scale = 2.38 * 2.38 / 2.0;
                    std::shared_ptr<MarkovChain::ProposalFunction> ppf(mvg);

                    MarkovChain chain(analysis, 13, ppf);
                    MarkovChain::HyperParameter h;
                    h.component = 0;
                    std::vector<double> point_initial { 4.3, 1.15 };
                    chain.set_point(point_initial, h);

                    unsigned iterations = 1e5;

                    chain.run(iterations);

                    auto stats  = chain.statistics();
                    double efficiency = stats.iterations_accepted;
                    efficiency /= double(stats.iterations_accepted + stats.iterations_rejected);
//                    TEST_CHECK_NEARLY_EQUAL(0.25, efficiency, 1e-4);
                    TEST_CHECK( 0.2 < efficiency && efficiency < 0.3);

                    std::tuple<double, double> result;
                    chain.normalized_density(result, point_initial, iterations);
                    const double & enumerator = std::get<0>(result);
                    const double & denominator = std::get<1>(result);

                    // likelihood at the point
                    p["mass::b(MSbar)"] = point_initial[0];
                    p["mass::c"] = point_initial[1];
//
//                    // evidence = 1/range = prior
//                    // normalized posterior = likelihood
                    // accuracy is quite low, even though it it only a 2D problem
                    TEST_CHECK_RELATIVE_ERROR(analysis.log_likelihood()(), std::log(enumerator / denominator), 5.5e-2);

                } //);
#if 0
                //todo sample _covariance not positive definite. What to do? Tried initial_scale, but didn't help
                // global-local
                {
                    double initial_scale = 2.38 * 2.38 / 2.0;
                    // run two chains
                    proposal_functions::MultivariateGaussian* mvg1 = new proposal_functions::MultivariateGaussian(2, cov_initial);
                    mvg1->covariance_scale = initial_scale;
                    std::shared_ptr<MarkovChain::ProposalFunction> ppf1(mvg1);
                    proposal_functions::MultivariateGaussian* mvg2 = new proposal_functions::MultivariateGaussian(2, cov_initial);
                    mvg2->covariance_scale = initial_scale;
                    std::shared_ptr<MarkovChain::ProposalFunction> ppf2(mvg2);

                    std::vector<MarkovChain> chains;
                    chains.push_back(MarkovChain(analysis, 1, ppf1));
                    chains.push_back(MarkovChain(analysis, 2, ppf2));

                    unsigned iterations = 5e3;
                    for (auto c= chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
                    {
                        c->run(iterations);
                    }

                    proposal_functions::GlobalLocal gl(chains, 0.5, 2);

                    gsl_rng * rng;
                    rng = gsl_rng_alloc(gsl_rng_mt19937);
                    gsl_rng_set(rng, 4);

                    MarkovChain::State proposal(chains[0].proposed_state());

                    // pick 2nd component
                    gl.propose(proposal, chains[0].current_state(), rng);
                    TEST_CHECK_NEARLY_EQUAL(gl.evaluate(proposal, chains[0].current_state()), 0.717062241009793, eps);

                    // nonlocal, then pick 1st component
                    gl.propose(proposal, chains[0].current_state(), rng);

                    gl.propose(proposal, chains[0].current_state(), rng);

                    // local jump
                    gl.propose(proposal, chains[0].current_state(), rng);
                    TEST_CHECK_NEARLY_EQUAL(gl.evaluate(proposal, chains[0].current_state()), 1.5584608222464997, eps);

                    // local jump
                    gl.propose(proposal, chains[0].current_state(), rng);
                }
#endif
            }

        // Monte Carlo integration
            {
                // 1D Gaussian integration
                Analysis analysis = make_analysis(true);

                std::vector<double> cov_initial { 0.5 * 0.5 };
                proposal_functions::MultivariateGaussian* mvg = new proposal_functions::MultivariateGaussian(1, cov_initial);
                // automatic scaling is not perfect very low number of dimension, so revert its effect
                mvg->rescale(1.0 / mvg->covariance_scale);
                std::shared_ptr<MarkovChain::ProposalFunction> ppf(mvg);

                MarkovChain chain(analysis, 13, ppf);
                MarkovChain::HyperParameter h { 0 };
                std::vector<double> point_initial{ 4.2 };
                chain.set_point(point_initial, h);

                unsigned iterations = 1e5;

                chain.run(iterations);

                auto stats  = chain.statistics();
                double efficiency = stats.iterations_accepted;
                efficiency /= double(stats.iterations_accepted + stats.iterations_rejected);
                TEST_CHECK( 0.2 < efficiency && efficiency < 0.3);

                std::tuple<double, double> result;
                chain.normalized_density(result, point_initial, iterations);
                const double & enumerator = std::get<0>(result);
                const double & denominator = std::get<1>(result);

                // evidence = 1/range = prior
                // normalized posterior = likelihood
                // don't forget to exponentiate log likelihood output
                analysis.parameters()["mass::b(MSbar)"] = 4.2;
                TEST_CHECK_RELATIVE_ERROR(std::exp(analysis.log_likelihood()()), enumerator / denominator , 6e-3);

                // At a point away from the maximum, the precision isn't necessarily lower!
                point_initial[0] = 4.3;
                chain.run(iterations);
                chain.normalized_density(result, point_initial, iterations);

                analysis.parameters()["mass::b(MSbar)"] = 4.3;
                TEST_CHECK_RELATIVE_ERROR(std::exp(analysis.log_likelihood()()), enumerator / denominator , 2.5e-3);
            }

            // test History
            {
                MarkovChain::History history;
                std::vector<double> means, variances;

                MarkovChain::State s;

                s.point = std::vector<double> { 1.2, 3.3 };
                history.states.push_back(s);

                history.mean_and_variance(history.states.begin(), history.states.end(), means, variances);
                TEST_CHECK_EQUAL(means[0], 1.2);
                TEST_CHECK_EQUAL(means[1], 3.3);
                TEST_CHECK_EQUAL(variances[0], 0.0);
                TEST_CHECK_EQUAL(variances[1], 0.0);

                s.point = std::vector<double> { 2.3, 4.5 };
                history.states.push_back(s);
                history.mean_and_variance(history.states.begin(), history.states.end(), means, variances);
                TEST_CHECK_EQUAL(means[0], 1.75);
                TEST_CHECK_EQUAL(means[1], 3.9);
                TEST_CHECK_RELATIVE_ERROR(variances[0], 0.605, eps);
                TEST_CHECK_RELATIVE_ERROR(variances[1], 0.72, eps);

                s.point = std::vector<double> { 2.8, 4.1 };
                history.states.push_back(s);
                history.mean_and_variance(history.states.begin(), history.states.end(), means, variances);
                TEST_CHECK_RELATIVE_ERROR(means[0], 2.1, eps);
                TEST_CHECK_RELATIVE_ERROR(means[1], 11.9 / 3.0, eps);
                TEST_CHECK_RELATIVE_ERROR(variances[0], 0.67, eps);
                TEST_CHECK_RELATIVE_ERROR(variances[1], 0.37 + 1 / 300.0, eps);

                // skip first two elements
                MarkovChain::State::Iterator it = history.states.begin();
                it += 2;
                history.mean_and_variance(it, history.states.end(), means, variances);
                TEST_CHECK_EQUAL(means[0], 2.8);
                TEST_CHECK_EQUAL(means[1], 4.1);
                TEST_CHECK_EQUAL(variances[0], 0);
                TEST_CHECK_EQUAL(variances[1], 0);

                TEST_CHECK_THROWS(InternalError, history.mean_and_variance(it, it, means, variances));
            }
            // random index
          {
                gsl_rng * rng;
                rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, 46);

                std::vector<double> cumulative { 0.22, 0.4, 0.6, 0.8, 1 };
                TEST_CHECK_EQUAL(proposal_functions::random_index(cumulative, rng), 3);
                TEST_CHECK_EQUAL(proposal_functions::random_index(cumulative, rng), 1);
                TEST_CHECK_EQUAL(proposal_functions::random_index(cumulative, rng), 3);
                TEST_CHECK_EQUAL(proposal_functions::random_index(cumulative, rng), 3);
                TEST_CHECK_EQUAL(proposal_functions::random_index(cumulative, rng), 1);
                TEST_CHECK_EQUAL(proposal_functions::random_index(cumulative, rng), 4);
                TEST_CHECK_EQUAL(proposal_functions::random_index(cumulative, rng), 3);
                TEST_CHECK_EQUAL(proposal_functions::random_index(cumulative, rng), 0);
            }
        }
} markov_chain_test;
