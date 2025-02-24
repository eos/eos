/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Frederik Beaujean
 * Copyright (c) 2011-2023 Danny van Dyk
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

#include <eos/statistics/log-posterior_TEST.hh>

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

                auto clone1 = log_posterior.clone();
                auto clone2 = log_posterior.clone();

                // make sure observable's value is not equal to central value
                Parameter p1 = (*clone1)[0];
                p1.set(4.3); //posterior mode
                Parameter p2 = (*clone2)[0];
                p2.set(4.4); //log_prior mode

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
                LogPosteriorPtr clone3 = log_posterior.clone();

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
                auto clone_bare = log_posterior.clone();


                // 4.4 +- 0.1
                log_posterior.add(LogPrior::CurtailedGauss(parameters, "mass::b(MSbar)", 3.7, 4.9, 4.3, 4.4, 4.5));



                Parameter p = log_posterior[0];
                p.set(4.3); //posterior mode


                TEST_CHECK_NEARLY_EQUAL(log_posterior.log_likelihood()(), +0.88364655978936768, eps);
                TEST_CHECK_NEARLY_EQUAL(log_posterior.log_prior(), 0.883646846442260436, eps);
                // slightly different due to normalization of prior
                TEST_CHECK(log_posterior.log_likelihood()() != log_posterior.log_prior());

                // now check cloning

                auto clone = log_posterior.clone();
                Parameter p2 = (*clone)[0];

                TEST_CHECK_EQUAL(p.evaluate(), p2.evaluate());


                // change clone only
                p2.set(4.112);
                TEST_CHECK(log_posterior.log_likelihood()() != clone->log_likelihood()());
                TEST_CHECK(log_posterior.log_prior() != clone->log_prior());

                // same value for clone and original
                p2.set(4.3);

                TEST_CHECK_EQUAL(log_posterior.log_likelihood()(), clone->log_likelihood()());
                TEST_CHECK_EQUAL(log_posterior.log_prior(), clone->log_prior());
            }

            // stop if prior undefined
            {
                Parameters parameters = Parameters::Defaults();

                LogLikelihood llh(parameters);
                llh.add(ObservablePtr(new ObservableStub(parameters, "mass::b(MSbar)")), 4.1, 4.2, 4.3);
                LogPosterior log_posterior(llh);

                TEST_CHECK_THROWS(InternalError, log_posterior.log_prior());
            }
        }
} log_posterior_test;
