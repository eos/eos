/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Frederik Beaujean
 * Copyright (c) 2011, 2012, 2013, 2015, 2016 Danny van Dyk
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

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>

using namespace test;
using namespace eos;

class LogPosteriorTest :
    public TestCase
{
    public:
        LogPosteriorTest() :
            TestCase("log_posterior_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-13;

            // check cloning and values
            {
                LogPosterior log_posterior = make_log_posterior(false);

                auto clone1 = log_posterior.old_clone();
                auto clone2 = log_posterior.old_clone();

                // make sure observable's value is not equal to central value
                MutablePtr p = (*clone1)[0];
                p->set(4.3); //posterior mode
                p = (*clone2)[0];
                p->set(4.4); //log_prior mode

                // for comparison used ipython's log(scipy.stats.norm.pdf(4.3, loc=4.4, scale=0.1))
                // value at center of both Gaussian distributions. so pdf the same
                TEST_CHECK_RELATIVE_ERROR(clone1->log_likelihood()(), +0.88364655978936768, eps);
                TEST_CHECK_RELATIVE_ERROR(clone1->log_prior(), 0.883646846442260436, eps);

                // almost, but not quite identical
                TEST_CHECK_RELATIVE_ERROR(clone1->log_likelihood()(), clone1->log_prior(), 1e-6);

                TEST_CHECK_RELATIVE_ERROR(clone2->log_likelihood()(), -0.61635344021063077, eps);

                TEST_CHECK_RELATIVE_ERROR(clone2->log_prior(), 1.38364684644226932, eps);


                // now change a parameter which is not scanned
                TEST_CHECK(log_posterior.parameters()["b->s::Re{c7}"] != 2.599);
                log_posterior.parameters()["b->s::Re{c7}"] = 2.599;
                LogPosteriorPtr clone3 = log_posterior.old_clone();

                TEST_CHECK_EQUAL(double(log_posterior.parameters()["b->s::Re{c7}"] ),
                                 double( clone3->parameters()["b->s::Re{c7}"] ));

            }

            // smart parameter adding
            {
                Parameters parameters = Parameters::Defaults();

                LogLikelihood llh(parameters);
                llh.add(ObservablePtr(new ObservableStub(parameters, "mass::b(MSbar)")), 4.1, 4.2, 4.3);
                LogPosterior log_posterior(llh);

                // store a clone with no parameters
                auto clone_bare = log_posterior.old_clone();


                // 4.4 +- 0.1
                log_posterior.add(LogPrior::Gauss(parameters, "mass::b(MSbar)",
                        ParameterRange
                            { 3.7, 4.9 }, 4.3, 4.4, 4.5));



                MutablePtr p = log_posterior[0];
                p->set(4.3); //posterior mode


                TEST_CHECK_NEARLY_EQUAL(log_posterior.log_likelihood()(), +0.88364655978936768, eps);
                TEST_CHECK_NEARLY_EQUAL(log_posterior.log_prior(), 0.883646846442260436, eps);
                // slightly different due to normalization of prior
                TEST_CHECK(log_posterior.log_likelihood()() != log_posterior.log_prior());

                // now check cloning

                auto clone = log_posterior.old_clone();
                MutablePtr p2 = (*clone)[0];

                TEST_CHECK_EQUAL(p->evaluate(), p2->evaluate());


                // change clone only
                p2->set(4.112);
                TEST_CHECK(log_posterior.log_likelihood()() != clone->log_likelihood()());
                TEST_CHECK(log_posterior.log_prior() != clone->log_prior());

                // same value for clone and original
                p2->set(4.3);

                TEST_CHECK_EQUAL(log_posterior.log_likelihood()(), clone->log_likelihood()());
                TEST_CHECK_EQUAL(log_posterior.log_prior(), clone->log_prior());
            }

            // nuisance properties.nuisance())
            {
                LogPosterior log_posterior = make_log_posterior(false);

                TEST_CHECK(!log_posterior.nuisance("mass::b(MSbar)"));

                log_posterior.add(LogPrior::Flat(log_posterior.parameters(), "mass::c", ParameterRange{ 1.4, 2.2}), true);

                TEST_CHECK(log_posterior.nuisance("mass::c"));
            }

            // stop if prior undefined
            {
                Parameters parameters = Parameters::Defaults();

                LogLikelihood llh(parameters);
                llh.add(ObservablePtr(new ObservableStub(parameters, "mass::b(MSbar)")), 4.1, 4.2, 4.3);
                LogPosterior log_posterior(llh);

                TEST_CHECK_THROWS(InternalError, log_posterior.log_prior());
            }

            // 1D optimization
            {
                LogPosterior log_posterior = make_log_posterior(false);

                std::vector<double> initial_guess(1, 4.161345);
                LogPosterior::OptimizationOptions options = LogPosterior::OptimizationOptions::Defaults();
                options.tolerance = 1e-5;
                options.initial_step_size = 0.1;
                auto pair = log_posterior.optimize(initial_guess, options);

                auto best_fit_parameter = pair.first;

                TEST_CHECK_NEARLY_EQUAL(best_fit_parameter.front(), 4.3, 1e-5);
                TEST_CHECK_NEARLY_EQUAL(pair.second, 1.7672934062316281 , 1e-8);
            }

            // 5D optimization
            {
                Parameters parameters = Parameters::Defaults();

                LogLikelihood llh(parameters);
                llh.add(ObservablePtr(new ObservableStub(parameters, "mass::b(MSbar)")),  4.1,      4.2,    4.3);
                llh.add(ObservablePtr(new ObservableStub(parameters, "mass::c")),         1.15,     1.2,    1.25);
                llh.add(ObservablePtr(new ObservableStub(parameters, "mass::s(2GeV)")),   5e-3,    10e-3,  15e-3);
                llh.add(ObservablePtr(new ObservableStub(parameters, "mass::t(pole)")), 171,      172,    173);
                llh.add(ObservablePtr(new ObservableStub(parameters, "mass::e")),       510.5e-6, 511e-6, 511.5e-6);

                LogPosterior log_posterior(llh);

                log_posterior.add(LogPrior::Flat(log_posterior.parameters(), "mass::b(MSbar)", ParameterRange{ 4     , 4.5   }));
                log_posterior.add(LogPrior::Flat(log_posterior.parameters(), "mass::c",        ParameterRange{ 1     , 2     }));
                log_posterior.add(LogPrior::Flat(log_posterior.parameters(), "mass::s(2GeV)",  ParameterRange{ 1e-3  , 25e-3 }));
                log_posterior.add(LogPrior::Flat(log_posterior.parameters(), "mass::t(pole)",  ParameterRange{ 168   , 177   }));
                log_posterior.add(LogPrior::Flat(log_posterior.parameters(), "mass::e",        ParameterRange{ 500e-6, 520e-6}));

                std::vector<double> initial_guess{ 4.1001, 1.90014, 3.00045e-3, 174.6345, 515.51e-6 };

                LogPosterior::OptimizationOptions options = LogPosterior::OptimizationOptions::Defaults();
                options.tolerance = 1e-5;
                options.initial_step_size = 0.1;
                std::vector<double> optimum = log_posterior.optimize(initial_guess, options).first;

                TEST_CHECK_NEARLY_EQUAL(optimum[0], 4.2   , 1e-5);
                TEST_CHECK_NEARLY_EQUAL(optimum[1], 1.2   , 1e-5);
                TEST_CHECK_NEARLY_EQUAL(optimum[2], 1e-2  , 1e-5);
                TEST_CHECK_NEARLY_EQUAL(optimum[3], 172   , 2e-5);
                TEST_CHECK_NEARLY_EQUAL(optimum[4], 511e-6, 1e-5);

                /* try again with Minuit */

                //somehow minuit doesn't coverge with 4.1001
                initial_guess[0] = 3.2;

                // use lowest accuracy
                auto config = LogPosterior::OptimizationOptions::Defaults();
                config.strategy_level = 0;
                const ROOT::Minuit2::FunctionMinimum & data_at_min = log_posterior.optimize_minuit(initial_guess, config);

                // check parameters at mode
                auto u_par = data_at_min.UserParameters();
                TEST_CHECK_NEARLY_EQUAL(u_par.Value(0), 4.2   , 1e-4);
                TEST_CHECK_NEARLY_EQUAL(u_par.Value(1), 1.2   , 1e-4);
                TEST_CHECK_NEARLY_EQUAL(u_par.Value(2), 1e-2  , 1e-4);
                TEST_CHECK_NEARLY_EQUAL(u_par.Value(3), 172   , 1e-4);
                TEST_CHECK_NEARLY_EQUAL(u_par.Value(4), 511e-6, 1e-4);

                // should find input uncertainties
                auto u_cov = data_at_min.UserCovariance();
                TEST_CHECK_NEARLY_EQUAL(sqrt(u_cov(0,0)) , 0.10   , 5e-3);
                TEST_CHECK_NEARLY_EQUAL(sqrt(u_cov(1,1)) , 0.05   , 5e-3);
                TEST_CHECK_NEARLY_EQUAL(sqrt(u_cov(2,2)) , 5e-3   , 5e-5);
                TEST_CHECK_NEARLY_EQUAL(sqrt(u_cov(3,3)) , 1      , 5e-2);
                TEST_CHECK_NEARLY_EQUAL(sqrt(u_cov(4,4)) , 5e-7   , 5e-9);

                // no correlation present
                TEST_CHECK_NEARLY_EQUAL(u_cov(0,1) / sqrt(fabs(u_cov(0,0) * u_cov(1,1))), 0   , 5e-3);
                TEST_CHECK_NEARLY_EQUAL(u_cov(1,3) / sqrt(fabs(u_cov(1,1) * u_cov(3,3))), 0   , 2e-2);
            }

            // goodness_of_fit
            {
                Parameters parameters = Parameters::Defaults();

                LogLikelihood llh(parameters);
                llh.add(ObservablePtr(new ObservableStub(parameters, "mass::c")), 1.182, 1.192, 1.202);
                llh.add(ObservablePtr(new ObservableStub(parameters, "mass::c")), 1.19, 1.2, 1.21);

                LogPosterior log_posterior(llh);

                log_posterior.add(LogPrior::Flat(log_posterior.parameters(), "mass::c", ParameterRange{ 1 , 2 }));

                // in middle of both observations
                std::vector<double> best_fit_parameter(1, 1.196);

                // each observation 0.4 sigma away from mode
//                auto ret = log_posterior.goodness_of_fit(best_fit_parameter);
//                TEST_CHECK_NEARLY_EQUAL(ret.first, 0.32, eps);
//                TEST_CHECK_NEARLY_EQUAL(ret.first, 0.57160764495333116, eps);

                // now use simulation.
                // p-value _not_ corrected for DoF, thus biased towards p=1
                auto ret = log_posterior.goodness_of_fit(best_fit_parameter, 5e4);
                TEST_CHECK_NEARLY_EQUAL(ret.first, 0.852143788, 5e-3);
                TEST_CHECK_NEARLY_EQUAL(ret.second, 0.57160, 5e-3);
            }
        }
} log_posterior_test;
