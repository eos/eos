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
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +4.24, +4.25, +4.30);
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

                // iteration
                {
                    LogLikelihood llh(p);

                    // add blocks manually
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")), +4.1, +4.2, +4.3);
                    llh.add(ObservablePtr(new TestObservable(p, k, "mass::c")), +1.15, +1.2, +1.25);

                    // now add a  constraint
                    auto obs = ObservablePtr(new TestObservable(p, k, "mass::e"));
                    llh.add(Constraint("harr", std::vector<ObservablePtr>{ obs },
                        std::vector<LogLikelihoodBlockPtr>{ LogLikelihoodBlock::LogGamma(llh.observable_cache(), obs, 0.1, 0.11, 0.13) } ));

                    // remember to evaluate likelihood to fill the cache
                    p["mass::b(MSbar)"] = 4.25;
                    p["mass::c"] = 1.3;
                    p["mass::e"] = 0.115;
                    llh();

                    // use cached values
                    TEST_CHECK_EQUAL(4.250, llh.observable_cache()[0]);
                    TEST_CHECK_EQUAL(1.300, llh.observable_cache()[1]);
                    TEST_CHECK_EQUAL(0.115, llh.observable_cache()[2]);

                    // check significances
                    auto c = llh.begin();
                    TEST_CHECK_RELATIVE_ERROR(-0.5, (**c->begin_blocks()).significance(), eps);
                    ++c;
                    TEST_CHECK_RELATIVE_ERROR(-2, (**c->begin_blocks()).significance(), eps);

                    std::string observable_values;
                    auto cache = llh.observable_cache();
                    for (unsigned i = 0 ; i < cache.size() ; ++i)
                    {
                        observable_values += cache.observable(i)->name() + " = " + stringify(double(cache[i])) + "; ";
                    }
                    TEST_CHECK_EQUAL("test-observable[mass::b(MSbar)] = 4.25; test-observable[mass::c] = 1.3; test-observable[mass::e] = 0.115; ", observable_values);

                    // check that looping gives proper results as well
                    std::string constraints_significances;
                    for (auto c = llh.begin(), c_end = llh.end() ; c != c_end ; ++c)
                    {
                        for (auto b = c->begin_blocks(), b_end = c->end_blocks() ; b != b_end ; ++b)
                        {
                                constraints_significances += c->name() + ": " + stringify((**b).significance()) + "; ";
                        }
                    }

                    // check observations
                    TEST_CHECK_EQUAL(3, llh.number_of_observations());
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

                // observables vary only by option, but identical in name => they should be different predictions
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

                    /* significance */
                    block = LogLikelihoodBlock::Gaussian(cache, obs, +4.1, +4.3, +4.4);

                    p["mass::b(MSbar)"] = 4.35;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(-0.5, block->significance(), eps);

                    p["mass::b(MSbar)"] = 4.25;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(+0.25, block->significance(), eps);
                }

                // LogGamma
                {
                    static const double low_eps = 5e-4;

                    ObservablePtr obs(new TestObservable(p, k, "mass::b(MSbar)"));
                    ObservableCache cache(p);

                    double min = 0.34;
                    double central = 0.53;
                    double max = 0.63;
                    auto log_gamma = LogLikelihoodBlock::LogGamma(cache, obs, min, central, max);

                    // the model prediction
                    p["mass::b(MSbar)"] = 0.57;
                    cache.update();

                    TEST_CHECK_RELATIVE_ERROR(log_gamma->evaluate(), +1.005543554, low_eps);

                    // pdf value at one sigma border
                    p["mass::b(MSbar)"] = central + 0.2;
                    cache.update();
                    double pdf_max = log_gamma->evaluate();

                    p["mass::b(MSbar)"] = central - 0.2;
                    cache.update();
                    double pdf_min = log_gamma->evaluate();

                    // away from the mode, pdf falls more rapidly where uncertainty is smaller
                    TEST_CHECK(pdf_max < pdf_min);

                    // construct with known parameters (expect no exception)
                    auto log_gamma_manual = LogLikelihoodBlock::LogGamma(cache, obs, min, central, max, 6.8790736808e-02, 3.8305604649e-01);
                    TEST_CHECK_RELATIVE_ERROR(log_gamma->evaluate(), log_gamma_manual->evaluate(), low_eps);

                    // cloning
                    ObservableCache cache_clone(Parameters::Defaults());
                    auto log_gamma_clone = log_gamma->clone(cache_clone);

                    // is clone independent?
                    cache_clone.parameters()["mass::b(MSbar)"] = 15;
                    cache_clone.update();
                    TEST_CHECK(log_gamma->evaluate() != log_gamma_clone->evaluate());

                    // does clone have same state?
                    cache_clone.parameters()["mass::b(MSbar)"] = double(p["mass::b(MSbar)"]);
                    cache_clone.update();
                    TEST_CHECK_RELATIVE_ERROR(log_gamma_clone->evaluate(), log_gamma->evaluate(), 1e-14);

                    /* significance */
                    // last two examples show that when the point is at the one sigma interval boundary
                    // we get one sigma deviation as desired
                    cache_clone.parameters()["mass::b(MSbar)"] = min;
                    cache_clone.update();
                    TEST_CHECK_RELATIVE_ERROR(+1.0, log_gamma_clone->significance(), 2.1e-6);
                    cache_clone.parameters()["mass::b(MSbar)"] = max;
                    cache_clone.update();
                    TEST_CHECK_RELATIVE_ERROR(-1.0, log_gamma_clone->significance(), 2.1e-6);
                }

                // compare gaussian and loggamma
                {
                    ObservablePtr obs(new TestObservable(p, k, "mass::b(MSbar)"));
                    ObservableCache cache(p);

                    double min = 0.42;
                    double central = 0.53;
                    double max = 0.63;
                    auto log_gamma = LogLikelihoodBlock::LogGamma(cache, obs, min, central, max);
                    auto gauss = LogLikelihoodBlock::Gaussian(cache, obs, min, central, max);

                    /* agreement not very precise due to slight asymmetry of uncertainties */

                    // the model prediction
                    p["mass::b(MSbar)"] = 0.53;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(log_gamma->evaluate(), gauss->evaluate(), 3.8e-2);

                    p["mass::b(MSbar)"] = 0.63;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(log_gamma->evaluate(), gauss->evaluate(), 7e-2);

                    p["mass::b(MSbar)"] = 0.34;
                    cache.update();
                    TEST_CHECK_NEARLY_EQUAL(log_gamma->evaluate(), gauss->evaluate(), 0.13);

                    min = 0.425;
                    central = 0.53;
                    max = 0.63;
                    log_gamma = LogLikelihoodBlock::LogGamma(cache, obs, min, central, max);
                    gauss = LogLikelihoodBlock::Gaussian(cache, obs, min, central, max);

                    /* agreement not very precise due to slight asymmetry of uncertainties */

                    // the model prediction
                    p["mass::b(MSbar)"] = 0.53;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(log_gamma->evaluate(), gauss->evaluate(), 3.7e-2);

                    p["mass::b(MSbar)"] = 0.63;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(log_gamma->evaluate(), gauss->evaluate(), 3.5e-2);

                    p["mass::b(MSbar)"] = 0.34;
                    cache.update();
                    TEST_CHECK_NEARLY_EQUAL(log_gamma->evaluate(), gauss->evaluate(), 1e-1);
                }

                // loggamma sampling
                {
                    // sampling almost same as for gaussian, but with slight asymmetry
                    ObservablePtr obs(new TestObservable(p, k, "mass::b(MSbar)"));
                    ObservableCache cache(p);

                    double min = 4.195;
                    double central = 4.3;
                    double max = 4.4;
                    auto log_gamma = LogLikelihoodBlock::LogGamma(cache, obs, min, central, max);

                    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                    gsl_rng_set(rng, 2022);

                    log_gamma->prepare_sampling();

                    unsigned n = 1e4;
                    unsigned n_in = 0;

                    // pdf value at one sigma border
                    p["mass::b(MSbar)"] = max;
                    cache.update();
                    double pdf_max = log_gamma->evaluate();

                    p["mass::b(MSbar)"] = min;
                    cache.update();
                    double pdf_min = log_gamma->evaluate();

                    // both value should be close, so average seems justified for a rough approximation
                    double pdf_avg = 0.5 * (pdf_min + pdf_max);

                    for (unsigned i = 0; i < n; ++i)
                    {
                        double sample = log_gamma->sample(rng);

                        // check how many within one sigma
                        if (sample > pdf_avg)
                            ++n_in;
                    }

                    // should have 1 sigma probability in interval
                    TEST_CHECK_RELATIVE_ERROR(n_in / double(n), 0.684, 4e-3);

                    gsl_rng_free(rng);
                }

                // Amoroso limit
                {
                    ObservablePtr obs(new TestObservable(p, k, "mass::b(MSbar)"));
                    ObservableCache cache(p);

                    double mode = 0.0;
                    double x_90 = 0.6;
                    double x_95 = 0.8;

                    // construction checks
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::AmorosoLimit(cache, obs, mode, x_90, x_95, 0.2, 2));

                    // construction with correct parameters
                    auto amoroso = LogLikelihoodBlock::AmorosoLimit(cache, obs, mode, x_90, x_95, 0.19963939, 1.16193254);

                    // evaluation
                    p["mass::b(MSbar)"] = 0.53;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(amoroso->evaluate(), -0.782449964008941, 1e-8);

                    p["mass::b(MSbar)"] = 0.9101;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(amoroso->evaluate(), -2.155377494863876, 1e-8);

                    // cloning
                    ObservableCache cache_clone(Parameters::Defaults());
                    auto amoroso_clone = amoroso->clone(cache_clone);

                    // is clone independent?
                    cache_clone.parameters()["mass::b(MSbar)"] = 15;
                    cache_clone.update();
                    TEST_CHECK(amoroso->evaluate() != amoroso_clone->evaluate());

                    // does clone have same state?
                    cache_clone.parameters()["mass::b(MSbar)"] = 0.9101;
                    cache_clone.update();
                    TEST_CHECK_EQUAL(amoroso->evaluate(), amoroso_clone->evaluate());

                    /* sampling */
                    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                    gsl_rng_set(rng, 1357);

                    amoroso->prepare_sampling();

                    unsigned n = 1e4;
                    unsigned n_in_90 = 0;
                    unsigned n_in_95 = 0;

                    // pdf value at one sigma border
                    p["mass::b(MSbar)"] = x_90;
                    cache.update();
                    double pdf_90 = amoroso->evaluate();

                    p["mass::b(MSbar)"] = x_95;
                    cache.update();
                    double pdf_95 = amoroso->evaluate();

                    for (unsigned i = 0; i < n; ++i)
                    {
                        double sample = amoroso->sample(rng);

                        // check how many within one sigma
                        if (sample > pdf_90)
                            ++n_in_90;
                        if (sample > pdf_95)
                            ++n_in_95;
                    }

                    // should have 1 sigma probability in interval
                    TEST_CHECK_RELATIVE_ERROR(n_in_90 / double(n), 0.9, 1e-3);
                    TEST_CHECK_RELATIVE_ERROR(n_in_95 / double(n), 0.95, 1e-3);

                    gsl_rng_free(rng);

                    /* significance */

                    // limits at 0.9 and 1.0
                    amoroso = LogLikelihoodBlock::AmorosoLimit(cache, obs, 0.0, 0.9, 1.0, 0.99739274853, 0.19554265174);

                    p["mass::b(MSbar)"] = 0.1;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(0.137109, amoroso->significance(), 3e-6);

                    // the one sigma distance
                    p["mass::b(MSbar)"] = 0.63616;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(1.0, amoroso->significance(), 2e-7);
                }

                // Amoroso
                {
                    ObservablePtr obs(new TestObservable(p, k, "mass::b(MSbar)"));
                    ObservableCache cache(p);

                    // use 2011 LHCb/CMS limit on B_s -> mu mu
                    {
                        static const double physical_limit = 0.0;
                        static const double x_10 = 0.132749474699 * 10;
                        static const double x_50 = 0.446663009589 * 10;
                        static const double x_90 = 0.932149816388 * 10;
                        static const double theta = 6.4184393253;
                        static const double alpha = 8.1583565997e-01;
                        static const double beta = 1.8230347158;

                        // construction checks
                        TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::Amoroso(cache, obs, physical_limit, x_10, x_50, x_90, 0.2, 2, 3));

                        // construction with correct parameters
                        auto amoroso = LogLikelihoodBlock::Amoroso(cache, obs, physical_limit, x_10, x_50, x_90, theta, alpha, beta);

                        // evaluation at mode
                        p["mass::b(MSbar)"] = 3.112559;
                        cache.update();
                        TEST_CHECK_RELATIVE_ERROR(amoroso->evaluate(), std::log(1.332261877086652e-01), 1e-8);
                    }

                    // use 2012 LHCblimit on B_s -> mu mu
                    {
                        static const double physical_limit = 0.0;
                        static const double x_10 = 0.558367940293;
                        static const double x_50 = 2.03115589965;
                        static const double x_90 = 4.4528950788;
                        static const double theta = 2.9708273062;
                        static const double alpha = 8.2392613044e-01;
                        static const double beta = 1.6993290032;

                        // construction checks
                        TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::Amoroso(cache, obs, physical_limit, x_10, x_50, x_90, 0.2, 2, 3));

                        // construction with correct parameters
                        auto amoroso = LogLikelihoodBlock::Amoroso(cache, obs, physical_limit, x_10, x_50, x_90, theta, alpha, beta);
                        std::cout << amoroso->as_string() << std::endl;
                        // evaluation at mode
                        p["mass::b(MSbar)"] = 1.268439;
                        cache.update();
                        TEST_CHECK_RELATIVE_ERROR(amoroso->evaluate(), std::log(2.824624787700217e-01), 1e-8);

                        // significance
                        p["mass::b(MSbar)"] = 0.516344136;
                        cache.update();
                        TEST_CHECK_RELATIVE_ERROR(amoroso->significance(), +0.639662, 1e-5);
                        std::cout << "significance: " << amoroso->significance() << std::endl;

                        p["mass::b(MSbar)"] = 4.016344136;
                        cache.update();
                        TEST_CHECK_RELATIVE_ERROR(amoroso->significance(), -1.45552, 1e-5);
                    }
                }

                // multivariate gaussian
                {
                    std::array<ObservablePtr, 2> obs
                    {{
                        ObservablePtr(new TestObservable(p, k, "mass::b(MSbar)")),
                        ObservablePtr(new TestObservable(p, k, "mass::c"))
                    }};

                    ObservableCache cache(p);

                    // start with two uncorrelated Gaussians
                    std::array<double, 2> mean{{ 4.3, 1.1 }};
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

                    /* significance */
                    TEST_CHECK_RELATIVE_ERROR(1.557158038223962, block->significance(), eps);
                    TEST_CHECK_EQUAL(2, block->number_of_observations());

                    // with correlation, results are slightly inaccurate due to matrix inversion and determinant
                    covariance[0][1] = covariance[1][0] = 0.003;
                    block = LogLikelihoodBlock::MultivariateGaussian<2>(cache, obs, mean, covariance);

                    TEST_CHECK_NEARLY_EQUAL(block->evaluate(), 1.30077135, 1e-8);

                    // chi^2 now bigger with correlation!!
                    TEST_CHECK_RELATIVE_ERROR(1.6834363845158217821, block->significance(), eps);

                    p["mass::b(MSbar)"] = 4.6;
                    p["mass::c"] = 1.3;
                    cache.update();
                    TEST_CHECK_NEARLY_EQUAL(block->evaluate(),-4.597666149, 1e-8);

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
                    TEST_CHECK_RELATIVE_ERROR(old_value, block_clone->evaluate(), eps);

                    // with updated parameters, results should differ
                    new_pars["mass::c"] = 1.232;
                    new_cache.update();
                    TEST_CHECK(block_clone->evaluate() != block->evaluate());

                    /* interface with correlation */
                    mean = std::array<double, 2> {{ -0.32, 0.2 }};
                    std::array<double, 2> variances{{ 0.1321, 0.0601 }};
                    std::array<std::array<double, 2>, 2> correlation;
                    correlation[0][0] = correlation[1][1] = 1.0;
                    correlation[0][1] = correlation[1][0] = 0.08;

                    // calculate covariance by hand
                    // covariance matrix = (( 0.1321 0.007128179571 )( 0.007128179571 0.0601 ))
                    covariance[0][0] = 0.1321;
                    covariance[1][1] = 0.0601;
                    covariance[0][1] = covariance[1][0] = correlation[0][1] * std::sqrt(variances[0] * variances[1]);

                    auto mvg_correlation = LogLikelihoodBlock::MultivariateGaussian(cache, obs, mean, variances, correlation);
                    auto mvg_covariance = LogLikelihoodBlock::MultivariateGaussian(cache, obs, mean, covariance);
                    TEST_CHECK_RELATIVE_ERROR(mvg_covariance->evaluate(), mvg_correlation->evaluate(), eps);
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

                // mixture density
                {
                    ObservableCache cache(p);

                    std::vector<LogLikelihoodBlockPtr> components
                    {
                        LogLikelihoodBlock::Gaussian(cache, ObservablePtr(new TestObservable(p, Kinematics(), "mass::c")),
                                                     -5, -4, -3),
                        LogLikelihoodBlock::Gaussian(cache, ObservablePtr(new TestObservable(p, Kinematics(), "mass::c")),
                                                     3, 4, 5)
                    };
                    std::vector<double> weights { 0.9, 0.1 };

                    auto m = LogLikelihoodBlock::Mixture(components, weights);

                    p["mass::c"] = 4;
                    cache.update();
                    const double pdf_suppressed = m->evaluate();

                    // modes far from each other, so only single pdf matters. Lower precision.
                    TEST_CHECK_RELATIVE_ERROR(pdf_suppressed, std::log(weights[1]) + components[1]->evaluate() , 5e-14);

                    p["mass::c"] = -4;
                    cache.update();
                    const double pdf_favored = m->evaluate();

                    // modes far from each other, so only single pdf matters. High precision.
                    TEST_CHECK_RELATIVE_ERROR(pdf_favored, std::log(weights[0]) + components[0]->evaluate() , 2e-15);

                    // ratio of pdfs at mode given by weight ratio
                    TEST_CHECK_RELATIVE_ERROR(pdf_favored, pdf_suppressed + std::log(weights[0] / weights[1]), 1e-12);
                }
            }
    } log_likelihood_test;
}
