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

#include <test/test.hh>
#include <eos/utils/analysis_TEST.hh>
#include <eos/utils/log_likelihood.hh>
#include <eos/utils/power_of.hh>
#include <algorithm>

using namespace test;
using namespace eos;

namespace eos
{
    class LogLikelihoodTest :
        public TestCase
    {
        public:
            LogLikelihoodTest() :
                TestCase("log_likelihood_test")
            {
            }

            virtual void run() const
            {
                Parameters p = Parameters::Defaults();

                Kinematics k;
                k.declare("s", 15.0);

                static double eps = 1e-14;

                // symmetric gaussian test
                {
                    LogLikelihood llh(p);
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +4.2, +4.3, +4.4);

                    p["mass::b(MSbar)"] = 4.2;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +0.88364655978937656 , eps);

                    p["mass::b(MSbar)"] = 4.4;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +0.88364655978937656 , eps);
                }

                // asymmetric gaussian
                // values differ at one sigma from mode
                {
                    LogLikelihood llh(p);
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +4.24, +4.25, +4.3);


                    p["mass::b(MSbar)"] = 4.2;
                    TEST_CHECK_NEARLY_EQUAL(llh(), -8.813768347217003, eps);

                    p["mass::b(MSbar)"] = 4.24;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +3.18623165278344, eps);


                    p["mass::b(MSbar)"] = 4.3;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +1.5767937403493217, eps);
                }

                // Multiple test
                // Just the sum of individual log(likelihood) terms
                {
                    LogLikelihood llh(p);
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +4.24, +4.25, +4.3);
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::c")),        +1.33, +1.82, +1.90);
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::tau")),      +1.85, +2.00, +2.18);

                    p["mass::b(MSbar)"] = 4.2;
                    p["mass::c"] = 1.5;
                    p["mass::tau"] = 2.28;

                    TEST_CHECK_NEARLY_EQUAL(llh(), -9.646618122332889, eps);
                }

                // clone test
                {
                    LogLikelihood llh1(p);
                    llh1.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +4.2, +4.3, +4.4);

                    p["mass::b(MSbar)"] = 4.2;
                    TEST_CHECK_NEARLY_EQUAL(llh1(), 0.88364655978937656 , eps);

                    LogLikelihood llh2 = llh1.clone();
                    TEST_CHECK_EQUAL(llh1(), llh2());

                    //change parameters of ll1, but not of llh2
                    p["mass::b(MSbar)"] = 4.30;
                    TEST_CHECK_NEARLY_EQUAL(llh1(), 1.383646559789377, eps);
                    TEST_CHECK_NEARLY_EQUAL(llh2(), 0.88364655978937656, eps);

                    llh2.parameters()["mass::b(MSbar)"] = 4.60;
                    TEST_CHECK_NEARLY_EQUAL(llh1(), 1.383646559789377, eps);
                    TEST_CHECK_NEARLY_EQUAL(llh2(), -3.116353440210579, eps);
                }

                // multiple instances of same observable, to mimic results from different experiments
                {
                    LogLikelihood llh(p);
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +4.1, +4.2, +4.3);
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +4.3, +4.4, +4.5);

                    p["mass::b(MSbar)"] = 4.30;
                    TEST_CHECK_NEARLY_EQUAL(llh(), 2 * 0.88364655978937656, eps);

                    // only one prediction
                    TEST_CHECK_EQUAL(llh.observable_cache().size(), 1);
                }

                // observables vary only by kinematic, but identical in name => they should be different _predictions
                {
                    Kinematics kin;
                    kin.declare("s", 1.0);

                    LogLikelihood llh(p);
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")),   +4.1, +4.2, +4.3);
                    llh.add(ObservablePtr(new TestObservable(p, kin, "mass::b(MSbar)")), +4.3, +4.4, +4.5);

                    p["mass::b(MSbar)"] = 4.30;
                    TEST_CHECK_NEARLY_EQUAL(llh(), 2 * 0.88364655978937656, eps);

                    // two different _predictions
                    TEST_CHECK_EQUAL(llh.observable_cache().size(), 2);
                }

                // observables vary only by option, but identical in name => they should be different predicitions
                {
                    LogLikelihood llh(p);
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +4.1, +4.2, +4.3);

                    TestObservable * obs = new TestObservable(p, k, "mass::b(MSbar)");
                    obs->set_option("opt", "har");
                    llh.add(ObservablePtr(obs), +4.3, +4.4, +4.5);

                    p["mass::b(MSbar)"] = 4.30;
                    TEST_CHECK_NEARLY_EQUAL(llh(), 2 * 0.88364655978937656, eps);

                    // two different _predictions
                    TEST_CHECK_EQUAL(llh.observable_cache().size(), 2);
                }

                //verify that evaluation of observables only when needed works
                {
                    LogLikelihood llh(p);
                    ObservablePtr obs(new TestObservable(p, k, "mass::b(MSbar)"));
                    llh.add(obs, +4.2, +4.3, +4.4);

                    //add another observable
                    ObservablePtr obs2(new TestObservable(p, k, "Abs{c10'}"));
                    llh.add(obs2, 5, +5.3, +5.4);
                    p["Abs{c10'}"] = 5.3;



                    Parameter par = p["mass::b(MSbar)"];
                    par = 4.20;


                    TEST_CHECK_NEARLY_EQUAL(llh(), 1.16868083091064 , eps);

                    Parameter irrelevant_par = p["mass::e"];
                    irrelevant_par = 21;

                    //is it really irrelevant?
                    auto result = std::find(obs->begin(), obs->end(), irrelevant_par.id() );
                    TEST_CHECK(result == obs->end() );

                    //another irrelevant parameter
                    result = std::find(obs->begin(), obs->end(), p["Abs{c10}"].id() );
                    TEST_CHECK(result == obs->end() );

                    //only one is relevant though
                    result = std::find(obs->begin(), obs->end(), par.id() );
                    TEST_CHECK(result != obs->end() );

                    unsigned length = 0;
                    for ( auto dep = obs->begin(), dep_end = obs->end(); dep != dep_end; ++dep_end )
                        ++length;

                    TEST_CHECK_EQUAL(length, 1);

                    //no reevaluation
                    TEST_CHECK_NEARLY_EQUAL(llh(irrelevant_par.id()), 1.16868083091064 , eps);

                    //need to reevaluate
                    p["mass::b(MSbar)"] = 4.30;

                    TEST_CHECK_NEARLY_EQUAL(llh(par.id()), 1.66868083091064 , eps);

                    //evaluate
                    llh();

                    TEST_CHECK_EQUAL(llh.observable_cache()[1], 5.3);

                    //change value through the backdoor
                    p["Abs{c10'}"] = 6;

                    //same as before
                    llh(p["mass::b(MSbar)"].id());
                    TEST_CHECK_EQUAL(llh.observable_cache()[1], 5.3);

                    //now observable should be reevaluated
                    llh(p["Abs{c10'}"].id());
                    TEST_CHECK_EQUAL(llh.observable_cache()[1], 6);
                }

                //prevent user from using stale cache
                {
                    LogLikelihood llh(p);
                    ObservablePtr obs(new TestObservable(p, k, "mass::b(MSbar)"));
                    llh.add(obs, +4.2, +4.3, +4.4);

                    llh(p["mass::b(MSbar)"].id());

                    //add another observable
                    ObservablePtr obs2(new TestObservable(p, k, "Abs{c10'}"));
                    TEST_CHECK_THROWS(InternalError, llh.add(obs2, 5, +5.3, +5.4) );
                }

                // check single Gaussian block likelihood
                {
                    ObservablePtr obs(new TestObservable(p, k, "mass::b(MSbar)"));
                    ObservableCache cache(p);

                    auto block = LogLikelihoodBlock::Gaussian(cache, obs, +4.2, +4.3, +4.4);

                    // the model prediction
                    p["mass::b(MSbar)"] = 4.35;
                    cache.update();
                    block->prepare_sampling();

                    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                    gsl_rng_set(rng, 1243);

                    double sample = block -> sample(rng);
                    // likelihood from chi2 = (4.3 - 4.52)^2 / 0.1^2
                    TEST_CHECK_NEARLY_EQUAL(sample, -1.214746202801726, eps );

                    // now generate multiple samples for perfect fit
                    // thus a chi^2 distribution with 1 DoF
                    p["mass::b(MSbar)"] = 4.3;
                    cache.update();
                    block->prepare_sampling();

                    sample = block->sample(rng);

                    // likelihood from chi2 = (4.3 - 4.358)^2 / 0.1^2
                    TEST_CHECK_NEARLY_EQUAL(sample, 1.213893687750542, eps );

                    gsl_rng_set(rng, 15458);

                    unsigned n = 1e5;

                    double mean = 0;
                    unsigned n_in = 0;

                    for (unsigned i = 0; i < n; ++i)
                    {
                        sample = block->sample(rng);

                        // transform from llh to chi^2
                        sample = (sample - 1.383646559789373) * (-2);

                        mean += (sample - mean) / double(i + 1);

                        // check how many within one sigma
                        if (sample <= 1)
                            ++n_in;
                    }

                    // slight excess of samples: 68.4 % > 68.27
                    // where does it come from?
                    // checked different RNG algorithms, didn't make a difference
                    // We only allowed variations of the data within 3 sigma,
                    // in order to avoid possible unphysical values. This
                    // slightly pushed the bulk up towards the central value.
                    TEST_CHECK_NEARLY_EQUAL(n_in / double(n), 0.684 , 3e-3);

                    // as a consequence, mean always too low (<1), seems to converge around 0.972
                    TEST_CHECK_NEARLY_EQUAL(mean, 0.97, 1e-2);

                    // do not allow wrong input
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::Gaussian(cache, obs, +4.2, +4.3, +1.2));
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::Gaussian(cache, obs, +10., +4.3, +4.4));

                    gsl_rng_free(rng);
                }

                // multivariate gaussian
                {
                    std::vector<ObservablePtr> obs;
                    obs.push_back(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")));
                    obs.push_back(ObservablePtr(new TestObservable(p, k, "mass::c")));

                    ObservableCache cache(p);

                    // start with two uncorrelated Gaussians
                    std::array<double, 2> mean {{4.3, 1.1}};
                    std::array<std::array<double, 2>, 2> covariance;
                    covariance[0][0] = 0.1 * 0.1;
                    covariance[1][1] = 0.05 * 0.05;
                    covariance[0][1] = covariance[1][0] = 0;

                    auto block = LogLikelihoodBlock::MultivariateGaussian<2>(cache, obs, mean, covariance);

                    // create two one dim. gaussians to compare with
                    auto block1 = LogLikelihoodBlock::Gaussian(cache, obs[0], +4.20, +4.30, +4.40);
                    auto block2 = LogLikelihoodBlock::Gaussian(cache, obs[1], +1.05, +1.10, +1.15);

                    // update the common cache so observable now have values different from nan
                    p["mass::b(MSbar)"] = 4.35;
                    p["mass::c"] = 1.2;
                    cache.update();

                    // log of product of single pdfs is just the combined log
                    TEST_CHECK_NEARLY_EQUAL(block1->evaluate() + block2->evaluate(), block->evaluate(), 1e-13);

                    // with correlation, results are slightly inaccurate due to matrix inversion and determinant
                    covariance[0][1] = covariance[1][0] = 0.003;
                    block = LogLikelihoodBlock::MultivariateGaussian<2>(cache, obs, mean, covariance);

                    TEST_CHECK_NEARLY_EQUAL(block->evaluate(), 1.30077135, 1e-8);

                    /* test sampling */
                    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                    gsl_rng_set(rng, 1243);

                    block->prepare_sampling();
                    double sample = block -> sample(rng);

                    unsigned n2_in = 0;
                    unsigned n3_in = 0;
                    unsigned n = 1e5;

                    // prefactor of multivariate pdf
                    double normalization = -log(2 * M_PI) - 0.5 * log(1.6e-5);

                    for (unsigned i = 0; i < n; ++i)
                    {
                        sample = block->sample(rng);

                        // transform from llh to chi^2
                        sample = (sample - normalization) * (-2);

                        // check how many within a given quantile
                        if (sample <= 2)
                            ++n2_in;
                        // check how many within a given quantile
                        if (sample <= 3)
                            ++n3_in;
                    }

                    // compare with \chi^2 distribution with 2 DoF
                    TEST_CHECK_NEARLY_EQUAL(n2_in / double(n), 0.63212055882855767, 3e-3);
                    TEST_CHECK_NEARLY_EQUAL(n3_in / double(n), 0.77686983985157021, 3e-3);

                    /* cloning */
                    Parameters new_pars = p.clone();
                    ObservableCache new_cache(new_pars);
                    auto block_clone = block->clone(new_cache);
                    new_cache.update();

                    double old_value = block->evaluate();
                    TEST_CHECK_EQUAL(old_value, block_clone->evaluate());

                    // with updated parameters, results should differ
                    new_pars["mass::c"] = 1.232;
                    new_cache.update();
                    TEST_CHECK(block_clone->evaluate() != block->evaluate());
                }

                // bootstrap p-value calculation
                {
                    Parameters parameters  = Parameters::Defaults();
                    LogLikelihood llh(parameters);
                    llh.add(ObservablePtr(new TestObservable(parameters, Kinematics(),
                        "mass::c")), 1.182, 1.192, 1.202);
                    llh.add(ObservablePtr(new TestObservable(parameters, Kinematics(),
                        "mass::c")), 1.19, 1.2, 1.21);

                    parameters["mass::c"] = 1.196;
                    llh();

                    double p_value = llh.bootstrap_p_value(5e4).first;
                    // p-value from chi^2=0.32 and two degrees-of-freedom
                    // since data restricted to three sigma around central value,
                    // p-value should be slightly biased upwards
                    TEST_CHECK_NEARLY_EQUAL(p_value, 0.852143788, 5e-3);
                }
            }
    } log_likelihood_test;
}
