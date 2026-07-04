/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2026 Danny van Dyk
 * Copyright (c) 2011      Frederik Beaujean
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
#include <eos/statistics/log-likelihood.hh>
#include <eos/statistics/log-posterior_TEST.hh>
#include <eos/statistics/test-statistic-impl.hh>
#include <eos/maths/power-of.hh>
#include <algorithm>
#include <cmath>

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
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)), +4.2, +4.3, +4.4);

                    p["mass::b(MSbar)"] = 4.2;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +0.88364655978937656 , eps);

                    p["mass::b(MSbar)"] = 4.4;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +0.88364655978937656 , eps);
                }

                // asymmetric gaussian
                // values differ at one sigma from mode
                {
                    LogLikelihood llh(p);
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)), +4.24, +4.25, +4.3);

                    p["mass::b(MSbar)"] = 4.2;
                    TEST_CHECK_NEARLY_EQUAL(llh(), -9.912380635885128, eps);

                    p["mass::b(MSbar)"] = 4.24;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +2.087619364115315, eps);

                    p["mass::b(MSbar)"] = 4.3;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +2.0876193641153127, eps);
                }

                // Multiple test
                // Just the sum of individual log(likelihood) terms
                {
                    LogLikelihood llh(p);
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)), +4.24, +4.25, +4.30);
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::c",        k)), +1.33, +1.82, +1.90);
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::tau",      k)), +1.85, +2.00, +2.18);

                    p["mass::b(MSbar)"] = 4.2;
                    p["mass::c"] = 1.5;
                    p["mass::tau"] = 2.28;

                    TEST_CHECK_NEARLY_EQUAL(llh(), -10.11630282317536, eps);
                }

                // clone test
                {
                    LogLikelihood llh1(p);
                    llh1.add(ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)), +4.2, +4.3, +4.4);

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
                    std::cout << "FOO" << std::endl;
                    LogLikelihood llh(p);

                    // add blocks manually
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)), +4.1,  +4.2, +4.3);
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::c",        k)), +1.15, +1.2, +1.25);

                    // now add a  constraint
                    auto obs = ObservablePtr(new ObservableStub(p, "mass::e", k));
                    llh.add(Constraint("test::electron-mass", std::vector<ObservablePtr>{ obs },
                        std::vector<LogLikelihoodBlockPtr>{ LogLikelihoodBlock::LogGamma(llh.observable_cache(), obs, 0.1, 0.11, 0.13, 0.338082, -0.00649023) } ));

                    // remember to evaluate likelihood to fill the cache
                    p["mass::b(MSbar)"] = 4.25;
                    p["mass::c"] = 1.3;
                    p["mass::e"] = 0.115;
                    llh();

                    // use cached values
                    TEST_CHECK_EQUAL(4.250, llh.observable_cache()[ObservableCache::ObservableId(0)]);
                    TEST_CHECK_EQUAL(1.300, llh.observable_cache()[ObservableCache::ObservableId(1)]);
                    TEST_CHECK_EQUAL(0.115, llh.observable_cache()[ObservableCache::ObservableId(2)]);

                    // check significances
                    auto c = llh.begin();
                    TEST_CHECK_RELATIVE_ERROR((**c->begin_blocks()).significance(), -0.5, eps);
                    ++c;
                    TEST_CHECK_RELATIVE_ERROR((**c->begin_blocks()).significance(), -2, eps);

                    std::string observable_values;
                    auto cache = llh.observable_cache();
                    for (unsigned i = 0 ; i < cache.size() ; ++i)
                    {
                        auto id = ObservableCache::ObservableId(i);
                        observable_values += cache.observable(id)->name().str() + " = " + stringify(double(cache[id])) + "; ";
                    }
                    TEST_CHECK_EQUAL_STR("mass::b(MSbar) = 4.25; mass::c = 1.3; mass::e = 0.115; ", observable_values);

                    // check that looping gives proper results as well
                    std::string constraints_significances;
                    for (auto c = llh.begin(), c_end = llh.end() ; c != c_end ; ++c)
                    {
                        for (auto b = c->begin_blocks(), b_end = c->end_blocks() ; b != b_end ; ++b)
                        {
                                constraints_significances += c->name().str() + ": " + stringify((**b).significance()) + "; ";
                        }
                    }

                    // check observations
                    TEST_CHECK_EQUAL(3, llh.number_of_observations());
                }
                // multiple instances of same observable, to mimic results from different experiments
                {
                    LogLikelihood llh(p);
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)), +4.1, +4.2, +4.3);
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)), +4.3, +4.4, +4.5);

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
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)),   +4.1, +4.2, +4.3);
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", kin)), +4.3, +4.4, +4.5);

                    p["mass::b(MSbar)"] = 4.30;
                    TEST_CHECK_NEARLY_EQUAL(llh(), 2 * 0.88364655978937656, eps);

                    // two different _predictions
                    TEST_CHECK_EQUAL(llh.observable_cache().size(), 2);
                }

                // observables vary only by option, but identical in name => they should be different predictions
                {
                    LogLikelihood llh(p);
                    llh.add(ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)), +4.1, +4.2, +4.3);

                    ObservableStub * obs = new ObservableStub(p, "mass::b(MSbar);opt=har", k);
#if 0
                    obs->set_option("opt", "har");
#endif
                    llh.add(ObservablePtr(obs), +4.3, +4.4, +4.5);

                    p["mass::b(MSbar)"] = 4.30;
                    TEST_CHECK_NEARLY_EQUAL(llh(), 2 * 0.88364655978937656, eps);

                    // two different _predictions
                    TEST_CHECK_EQUAL(llh.observable_cache().size(), 2);
                }

                // check single Gaussian block likelihood
                {
                    ObservablePtr obs(new ObservableStub(p, "mass::b(MSbar)", k));
                    ObservableCache cache(p);

                    auto block = LogLikelihoodBlock::Gaussian(cache, obs, +4.2, +4.3, +4.4);

                    // the model prediction
                    p["mass::b(MSbar)"] = 4.35;
                    cache.update();

                    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                    gsl_rng_set(rng, 1243);

                    double sample = block -> sample(rng);
                    // likelihood from chi2 = (4.35 - 4.2975)^2 / 0.1^2
                    double target = 1.383332873466108;
                    TEST_CHECK_NEARLY_EQUAL(sample, target, eps );

                    // theory is completely irrelevant, it is added and subtracted internally
                    p["mass::b(MSbar)"] = 11234.35;
                    cache.update();
                    gsl_rng_set(rng, 1243);
                    TEST_CHECK_NEARLY_EQUAL(block->sample(rng), target, eps );

                    // likelihood from chi2 = (4.35 - 4.27296)^2 / 0.1^2
                    p["mass::b(MSbar)"] = 4.35;
                    cache.update();
                    TEST_CHECK_NEARLY_EQUAL(block->sample(rng), 1.086906027470852, eps );

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
                    // 1 sigma interval
                    TEST_CHECK_NEARLY_EQUAL(n_in / double(n), 0.68268949213708585, 1.0 / std::sqrt(n));

                    // chi^2 distribution with 1 dof has mean 1
                    TEST_CHECK_NEARLY_EQUAL(mean, 1, 1.0 / std::sqrt(n));

                    // do not allow wrong input
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::Gaussian(cache, obs, +4.2, +4.3, +1.2));
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::Gaussian(cache, obs, +10., +4.3, +4.4));
                    gsl_rng_free(rng);

                    /* significance */
                    block = LogLikelihoodBlock::Gaussian(cache, obs, +4.1, +4.3, +4.4);

                    p["mass::b(MSbar)"] = 4.35;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(block->significance(), -0.5,  eps);

                    p["mass::b(MSbar)"] = 4.25;
                    cache.update();
                    TEST_CHECK_RELATIVE_ERROR(block->significance(), +0.25, eps);
                }

                // LogGamma
                {
                    static const double low_eps = 5e-4;

                    ObservablePtr obs(new ObservableStub(p, "mass::b(MSbar)", k));
                    ObservableCache cache(p);

                    double min = 0.34;
                    double central = 0.53;
                    double max = 0.63;
                    auto log_gamma = LogLikelihoodBlock::LogGamma(cache, obs, min, central, max, 0.383056, 0.0687907);

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
                    auto log_gamma_manual = LogLikelihoodBlock::LogGamma(cache, obs, min, central, max, 3.8305604649e-01, 6.8790736808e-02);
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
                    TEST_CHECK_RELATIVE_ERROR(log_gamma_clone->significance(), +1.0, 1e-5);
                    cache_clone.parameters()["mass::b(MSbar)"] = max;
                    cache_clone.update();
                    TEST_CHECK_RELATIVE_ERROR(log_gamma_clone->significance(), -1.0, 1e-5);
                }

                // compare gaussian and loggamma
                {
                    ObservablePtr obs(new ObservableStub(p, "mass::b(MSbar)", k));
                    ObservableCache cache(p);

                    double min = 0.42;
                    double central = 0.53;
                    double max = 0.63;
                    auto log_gamma = LogLikelihoodBlock::LogGamma(cache, obs, min, central, max, 11.867, 0.358334);
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
                    TEST_CHECK_NEARLY_EQUAL(log_gamma->evaluate(), gauss->evaluate(), 0.132);

                    min = 0.425;
                    central = 0.53;
                    max = 0.63;
                    log_gamma = LogLikelihoodBlock::LogGamma(cache, obs, min, central, max, 80.2465, 0.916982);
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
                    ObservablePtr obs(new ObservableStub(p, "mass::b(MSbar)", k));
                    ObservableCache cache(p);

                    double min = 4.195;
                    double central = 4.3;
                    double max = 4.4;
                    auto log_gamma = LogLikelihoodBlock::LogGamma(cache, obs, min, central, max, 46.8496, 0.699917);

                    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                    gsl_rng_set(rng, 2022);

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

                // Amoroso
                {
                    ObservablePtr obs(new ObservableStub(p, "mass::b(MSbar)", k));
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
                        ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)),
                        ObservablePtr(new ObservableStub(p, "mass::c",        k))
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
                    TEST_CHECK_RELATIVE_ERROR(block->significance(), 1.557158038223962, eps);
                    TEST_CHECK_EQUAL(2, block->number_of_observations());

                    // with correlation, results are slightly inaccurate due to matrix inversion and determinant
                    covariance[0][1] = covariance[1][0] = 0.003;
                    block = LogLikelihoodBlock::MultivariateGaussian<2>(cache, obs, mean, covariance);

                    TEST_CHECK_NEARLY_EQUAL(block->evaluate(), 1.30077135, 1e-8);

                    // chi^2 now bigger with correlation!!
                    TEST_CHECK_RELATIVE_ERROR(block->significance(), 1.6834363845158217821, eps);

                    p["mass::b(MSbar)"] = 4.6;
                    p["mass::c"] = 1.3;
                    cache.update();
                    TEST_CHECK_NEARLY_EQUAL(block->evaluate(),-4.597666149, 1e-8);

                    /* test sampling */
                    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                    gsl_rng_set(rng, 1243);

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
                    llh.add(ObservablePtr(new ObservableStub(parameters, "mass::c")), 1.182, 1.192, 1.202);
                    llh.add(ObservablePtr(new ObservableStub(parameters, "mass::c")), 1.19, 1.2, 1.21);

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
                        LogLikelihoodBlock::Gaussian(cache, ObservablePtr(new ObservableStub(p, "mass::c")), -5, -4, -3),
                        LogLikelihoodBlock::Gaussian(cache, ObservablePtr(new ObservableStub(p, "mass::c")),  3,  4,  5)
                    };
                    std::vector<double> weights { 0.9, 0.1 };
                    std::vector<std::array<double, 2>> test_stat {};

                    auto m = LogLikelihoodBlock::Mixture(components, weights, test_stat);

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

                // A: UniformBound block (the factory is not exercised elsewhere)
                {
                    ObservablePtr   obs(new ObservableStub(p, "mass::c", k));
                    ObservableCache cache(p);
                    auto            block = LogLikelihoodBlock::UniformBound(cache, std::vector<ObservablePtr>{ obs }, 1.0, 0.1);

                    TEST_CHECK(block->as_string().find("UniformBound") != std::string::npos);
                    TEST_CHECK_EQUAL(block->number_of_observations(), 0u);

                    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
                    TEST_CHECK_EQUAL(block->sample(rng), 0.0);
                    gsl_rng_free(rng);

                    // saturation below the bound -> no penalty
                    p["mass::c"] = 0.5;
                    cache.update();
                    TEST_CHECK_EQUAL(block->evaluate(), 0.0);
                    TEST_CHECK_EQUAL(block->significance(), 0.0);

                    // saturation above the bound -> Gaussian penalty and non-zero significance
                    p["mass::c"] = 1.2;
                    cache.update();
                    TEST_CHECK_NEARLY_EQUAL(block->evaluate(), -0.5 * power_of<2>((1.2 - 1.0) / 0.1), eps);
                    TEST_CHECK_NEARLY_EQUAL(block->significance(), (1.2 - 1.0) / 0.1, eps);
                    TEST_CHECK_NO_THROW(block->primary_test_statistic());

                    // clone reproduces the parent (both saturated at 1.2)
                    Parameters      np = p.clone();
                    ObservableCache ncache(np);
                    auto            clone = block->clone(ncache);
                    np["mass::c"] = 1.2;
                    ncache.update();
                    TEST_CHECK_NEARLY_EQUAL(clone->evaluate(), block->evaluate(), eps);

                    // negative saturation is an error
                    p["mass::c"] = -0.5;
                    cache.update();
                    TEST_CHECK_THROWS(InternalError, block->evaluate());
                    TEST_CHECK_THROWS(InternalError, block->significance());

                    // with a zero uncertainty, saturation above the bound yields -infinity
                    auto block0 = LogLikelihoodBlock::UniformBound(cache, std::vector<ObservablePtr>{ ObservablePtr(new ObservableStub(p, "mass::c", k)) }, 1.0, 0.0);
                    p["mass::c"] = 1.2;
                    cache.update();
                    TEST_CHECK(! std::isfinite(block0->evaluate()));
                    TEST_CHECK(! std::isfinite(block0->significance()));
                }

                // B: as_string(), primary_test_statistic(), and clone() for the remaining block types
                {
                    // LogGamma
                    ObservablePtr   obs(new ObservableStub(p, "mass::b(MSbar)", k));
                    ObservableCache cache(p);
                    auto            lg = LogLikelihoodBlock::LogGamma(cache, obs, 0.34, 0.53, 0.63, 0.383056, 0.0687907);
                    p["mass::b(MSbar)"] = 0.53;
                    cache.update();
                    TEST_CHECK(! lg->as_string().empty());
                    TEST_CHECK_NO_THROW(lg->primary_test_statistic());
                    TEST_CHECK_NO_THROW(lg->clone(cache));
                }
                {
                    // Mixture
                    ObservableCache                    cache(p);
                    std::vector<LogLikelihoodBlockPtr> comps{
                        LogLikelihoodBlock::Gaussian(cache, ObservablePtr(new ObservableStub(p, "mass::c")), -5, -4, -3),
                        LogLikelihoodBlock::Gaussian(cache, ObservablePtr(new ObservableStub(p, "mass::c")),  3,  4,  5)
                    };
                    std::vector<double>                weights{ 0.9, 0.1 };
                    std::vector<std::array<double, 2>> test_stat{};
                    auto                               mix = LogLikelihoodBlock::Mixture(comps, weights, test_stat);
                    p["mass::c"] = -4;
                    cache.update();
                    TEST_CHECK(! mix->as_string().empty());
                    TEST_CHECK_NO_THROW(mix->primary_test_statistic());
                    TEST_CHECK_NO_THROW(mix->clone(cache));
                }
                {
                    // MultivariateGaussian
                    std::array<ObservablePtr, 2> obs{
                        { ObservablePtr(new ObservableStub(p, "mass::b(MSbar)", k)), ObservablePtr(new ObservableStub(p, "mass::c", k)) }
                    };
                    ObservableCache cache(p);
                    std::array<double, 2>                mean{ { 4.3, 1.1 } };
                    std::array<std::array<double, 2>, 2> covariance;
                    covariance[0][0] = 0.01;
                    covariance[1][1] = 0.0025;
                    covariance[0][1] = covariance[1][0] = 0.0;
                    auto mvg = LogLikelihoodBlock::MultivariateGaussian<2>(cache, obs, mean, covariance);
                    p["mass::b(MSbar)"] = 4.35;
                    p["mass::c"]        = 1.1;
                    cache.update();
                    TEST_CHECK(! mvg->as_string().empty());
                    TEST_CHECK_NO_THROW(mvg->primary_test_statistic());
                }

                // C: MixtureBlock::sample() is not implemented and must throw
                {
                    ObservableCache                    cache(p);
                    std::vector<LogLikelihoodBlockPtr> comps{ LogLikelihoodBlock::Gaussian(cache, ObservablePtr(new ObservableStub(p, "mass::c")), -1, 0, 1) };
                    std::vector<double>                weights{ 1.0 };
                    std::vector<std::array<double, 2>> test_stat{};
                    auto                               mix = LogLikelihoodBlock::Mixture(comps, weights, test_stat);

                    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
                    TEST_CHECK_THROWS(InternalError, mix->sample(rng));
                    gsl_rng_free(rng);
                }

                // D: factory error paths
                {
                    ObservablePtr   obs(new ObservableStub(p, "mass::c", k));
                    ObservableCache cache(p);

                    // LogGamma
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::LogGamma(cache, obs, 0.6, 0.5, 0.7, 0.3, 0.06));  // min >= central
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::LogGamma(cache, obs, 0.3, 0.5, 0.4, 0.3, 0.06));  // max <= central
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::LogGamma(cache, obs, 0.3, 0.5, 0.7, -1.0, 0.06)); // alpha <= 0

                    // Amoroso: mis-ordered limits (checked before the cdf-consistency checks)
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::Amoroso(cache, obs, 1.0, 0.5, 2.0, 3.0, 1.0, 1.0, 1.0)); // x_10 <= physical_limit
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::Amoroso(cache, obs, 1.0, 2.0, 0.5, 3.0, 1.0, 1.0, 1.0)); // x_50 <= physical_limit
                    TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::Amoroso(cache, obs, 0.0, 1.0, 2.0, 1.5, 1.0, 1.0, 1.0)); // x_90 <= x_50

                    // Mixture: components and weights sizes disagree
                    {
                        std::vector<LogLikelihoodBlockPtr> comps{ LogLikelihoodBlock::Gaussian(cache, obs, -1, 0, 1) };
                        std::vector<double>                weights{ 0.5, 0.5 };
                        std::vector<std::array<double, 2>> test_stat{};
                        TEST_CHECK_THROWS(InternalError, LogLikelihoodBlock::Mixture(comps, weights, test_stat));
                    }

                    // Unbinned1D
                    {
                        std::vector<Kinematics> kinematics{ Kinematics{} };
                        std::vector<double>     resolution{ 0.1 };
                        std::vector<Kinematics> observations{ Kinematics{} };
                        TEST_CHECK_THROWS(InternalError,
                                          LogLikelihoodBlock::Unbinned1D(cache, "TestLegendre1D::P(z)", std::vector<Kinematics>{}, Options{}, resolution, observations)); // empty kinematics
                        TEST_CHECK_THROWS(InternalError,
                                          LogLikelihoodBlock::Unbinned1D(cache, "TestLegendre1D::P(z)", kinematics, Options{}, std::vector<double>{ 0.1, 0.2 }, observations)); // size mismatch
                        TEST_CHECK_THROWS(InternalError,
                                          LogLikelihoodBlock::Unbinned1D(cache, "TestLegendre1D::P(z)", kinematics, Options{}, resolution, std::vector<Kinematics>{})); // empty observations
                    }
                }

                // E: externally-added likelihood blocks (LogLikelihood::add(block) and the external-block loops)
                {
                    Parameters    pe = Parameters::Defaults();
                    LogLikelihood llh(pe);

                    // add a standalone block directly (not via a constraint)
                    llh.add(LogLikelihoodBlock::Gaussian(llh.observable_cache(), ObservablePtr(new ObservableStub(pe, "mass::c")), 1.1, 1.2, 1.3));

                    pe["mass::c"] = 1.2;
                    TEST_CHECK(std::isfinite(llh()));

                    // cloning must carry the external block along
                    LogLikelihood llh_clone = llh.clone();
                    TEST_CHECK_NEARLY_EQUAL(llh_clone(), llh(), eps);

                    // an external block evaluating to -infinity short-circuits the total likelihood
                    llh.add(LogLikelihoodBlock::UniformBound(llh.observable_cache(), std::vector<ObservablePtr>{ ObservablePtr(new ObservableStub(pe, "mass::c")) }, 0.5, 0.0));
                    pe["mass::c"] = 1.2;
                    TEST_CHECK(! std::isfinite(llh()));
                }
            }
    } log_likelihood_test;

    class UnbinnedLogLikelihoodTest :
        public TestCase
    {
        public:
            UnbinnedLogLikelihoodTest() :
                TestCase("unbinned_log_likelihood_test")
            {
            }

            virtual void run() const
            {
                // The convolution of the test PDF 'TestLegendre1D::P(z)' with a Gaussian resolution
                // function (mean 0, standard deviation 0.1), evaluated on a uniform grid of N = 2^8
                // points spanning the smeared domain [z_smeared_min, z_smeared_max] = [-1.5, 5.5].
                //
                // The reference values below are the resolution-smeared density on that grid, computed
                // independently with Mathematica. They were generated with a grid spacing of 7/(N-1)
                // (inclusive endpoints), whereas the periodic grid underlying the discrete Fourier
                // transform has spacing 7/N; the factor (N-1)/N = 255/256 below corrects for this.
                static const std::size_t N = 256;
                static const double      z_smeared_min = -1.5;
                static const double      z_smeared_max = +5.5;
                static const double      dz            = (z_smeared_max - z_smeared_min) / N; // = 7/256
                static const double      sigma         = 0.1;

                static const std::array<double, N> reference{{
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 9.8096750507335319e-11, 4.8994547536473863e-10,
                2.2808818915303313e-09, 9.9003955209872559e-09, 4.0081628049196628e-08, 1.5140779927891626e-07,
                5.3388787147320642e-07, 1.7581873495209658e-06, 5.4105051871749915e-06, 1.5568520639437843e-05,
                4.1919296052406931e-05, 0.00010570750226186107, 0.00024988876108756496, 0.00055440287373774097,
                0.0011558617225092187, 0.0022679740834220794, 0.0041953727246893521, 0.0073309847622660407,
                0.012128152463587722, 0.019044987204231922, 0.028468940703678169, 0.04064010963427709,
                0.055596606477380954, 0.073160581334740057, 0.092970184593649696, 0.11454685534998099,
                0.1373763355568024, 0.1609803899497354, 0.18496366568521261, 0.2090315355161908,
                0.23298444582654129, 0.25669891391220345, 0.28010488476152967, 0.30316582942180786,
                0.32586420951757183, 0.34819226096975547, 0.37014681131077348, 0.3917266630519648,
                0.41293139836227799, 0.43376088240169014, 0.45421507490612867, 0.47429396474657715,
                0.49399754907484072, 0.51332582721581699, 0.53227879902127051, 0.55085646446104319,
                0.56905882352944925, 0.58688587622549537, 0.60433762254902035, 0.62141406250000009,
                0.63811519607843126, 0.65444102328431364, 0.67039154411764701, 0.68596675857843137,
                0.7011666666666666, 0.71599126838235294, 0.73044056372549016, 0.74451455269607825,
                0.75821323529411755, 0.77153661151960795, 0.78448468137254912, 0.79705744485294128,
                0.80925490196078442, 0.82107705269607856, 0.83252389705882368, 0.84359543504901957,
                0.85429166666666667, 0.86461259191176476, 0.87455821078431406, 0.88412852328431379,
                0.89332352941176485, 0.90214322916666678, 0.9105876225490197, 0.91865670955882361,
                0.92635049019607851, 0.9336689644607844, 0.9406121323529415, 0.94717999387254914,
                0.95337254901960788, 0.95918979779411773, 0.96463174019607867, 0.96969837622549004,
                0.97438970588235285, 0.97870572916666654, 0.98264644607843143, 0.98621185661764721,
                0.98940196078431364, 0.99221675857843128, 0.99465625000000024, 0.9967204350490193,
                0.99840931372549024, 0.99972288602941162, 1.0006611519607844, 1.0012241115196079,
                1.0014117647058822, 1.0012241115196079, 1.0006611519607844, 0.99972288602941184,
                0.99840931372549013, 0.99672043504901953, 0.99465625000000024, 0.9922167585784315,
                0.98940196078431386, 0.98621185661764721, 0.98264644607843143, 0.97870572916666676,
                0.97438970588235307, 0.96969837622549027, 0.96463174019607867, 0.95918979779411773,
                0.95337254901960788, 0.94717999387254914, 0.9406121323529415, 0.93366896446078418,
                0.92635049019607829, 0.91865670955882339, 0.9105876225490197, 0.90214322916666678,
                0.89332352941176463, 0.88412852328431357, 0.87455821078431406, 0.86461259191176443,
                0.85429166666666678, 0.84359543504901946, 0.83252389705882368, 0.82107705269607856,
                0.80925490196078442, 0.79705744485294128, 0.78448468137254912, 0.77153661151960784,
                0.75821323529411766, 0.74451455269607847, 0.73044056372549027, 0.71599126838235305,
                0.70116666666666672, 0.68596675857843148, 0.67039154411764712, 0.65444102328431386,
                0.63811519607843137, 0.62141406250000009, 0.60433762254902035, 0.58688587622549537,
                0.56905882352944925, 0.55085646446104319, 0.53227879902127051, 0.51332582721581699,
                0.49399754907484072, 0.47429396474657715, 0.45421507490612856, 0.43376088240169008,
                0.41293139836227799, 0.3917266630519648, 0.37014681131077348, 0.34819226096975542,
                0.32586420951757178, 0.30316582942180792, 0.28010488476152967, 0.2566989139122034,
                0.23298444582654129, 0.2090315355161908, 0.18496366568521261, 0.16098038994973554,
                0.13737633555680248, 0.11454685534998099, 0.092970184593649793, 0.073160581334740071,
                0.055596606477381003, 0.040640109634277083, 0.028468940703678156, 0.019044987204231999,
                0.012128152463587784, 0.0073309847622661031, 0.0041953727246893035, 0.0022679740834221449,
                0.0011558617225092449, 0.00055440287373768361, 0.00024988876108765246, 0.00010570750226200931,
                4.1919296052469313e-05, 1.5568520639573728e-05, 5.4105051872746119e-06, 1.7581873495215873e-06,
                5.3388787161756879e-07, 1.5140779946748343e-07, 4.0081628093235466e-08, 9.9003956910643878e-09,
                2.2808819007805453e-09, 4.8994546284965436e-10, 9.8096741187619583e-11, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0
                }};

                Parameters    p = Parameters::Defaults();
                ObservableCache cache(p);

                // Build the grid of kinematics at which the PDF is evaluated, and the resolution function
                // sampled on the same grid in wrap-around ("FFT-native") order: index 0 is the zero offset,
                // ascending indices are positive offsets, the high indices are the negative offsets.
                // The resolution evaluations span [-3.5, +3.5] (= [-N/2, N/2] grid steps from zero).
                //
                // The factor 1/4 = 1/max(z (4 - z)) peak-normalizes the *unnormalized* test PDF that the
                // block evaluates via SignalPDF::evaluate_linear() (whose maximum is 4 at z = 2), so that
                // the convolved grid matches the unit-peak reference values.
                std::vector<Kinematics> kinematics;
                std::vector<double>     resolution(N);
                for (std::size_t i = 0 ; i < N ; ++i)
                {
                    Kinematics k;
                    k.declare("z",     z_smeared_min + i * dz);
                    k.declare("z_min", z_smeared_min);
                    k.declare("z_max", z_smeared_max);
                    kinematics.push_back(k);

                    const double offset = (i <= N / 2 ? static_cast<double>(i) : static_cast<double>(static_cast<long>(i) - static_cast<long>(N))) * dz;
                    resolution[i] = std::exp(-offset * offset / (2.0 * sigma * sigma)) / (sigma * std::sqrt(2.0 * M_PI)) * dz / 4.0;
                }

                // Place one observation at each grid point at which the smeared density is safely positive,
                // so that interpolation returns the grid value and std::log() stays well-conditioned.
                std::vector<Kinematics> observations;
                double                  expected = 0.0;
                for (std::size_t i = 0 ; i < N ; ++i)
                {
                    const double value = reference[i] * 255.0 / 256.0;
                    if (value < 1.0e-3)
                        continue;

                    Kinematics k;
                    k.declare("z",     z_smeared_min + i * dz);
                    k.declare("z_min", z_smeared_min);
                    k.declare("z_max", z_smeared_max);
                    observations.push_back(k);

                    expected += std::log(value);
                }

                auto block = LogLikelihoodBlock::Unbinned1D(cache, "TestLegendre1D::P(z)", kinematics, Options{}, resolution, observations);

                // The block reads the unnormalized PDF values off the grid from the cache, so the cache
                // must be updated before evaluation (as LogLikelihood::operator()() does in production).
                cache.update();

                // The block convolves the PDF with the resolution and sums the logarithm of the smeared
                // density evaluated at each observation; this must reproduce the reference values.
                TEST_CHECK_NEARLY_EQUAL(block->evaluate(), expected, 1.0e-6);

                // exercise the remaining Unbinned1DLikelihoodBlock methods
                TEST_CHECK(block->as_string().find("Unbinned") != std::string::npos);
                TEST_CHECK_NO_THROW(block->primary_test_statistic());
                TEST_CHECK_NO_THROW(block->clone(cache));

                // sample() and significance() are not implemented and must throw
                gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
                TEST_CHECK_THROWS(InternalError, block->sample(rng));
                gsl_rng_free(rng);
                TEST_CHECK_THROWS(InternalError, block->significance());
            }
    } unbinned_log_likelihood_test;
}
