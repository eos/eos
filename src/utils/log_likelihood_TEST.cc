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
#include <src/utils/analysis_TEST.hh>
#include <src/utils/log_likelihood.hh>
#include <src/utils/power_of.hh>
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
                    TEST_CHECK_EQUAL(llh.predictions().size(), 1);
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
                    TEST_CHECK_EQUAL(llh.predictions().size(), 2);
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
                    TEST_CHECK_EQUAL(llh.predictions().size(), 2);
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

                    TEST_CHECK_EQUAL(llh.predictions()[1], 5.3);

                    //change value through the backdoor
                    p["Abs{c10'}"] = 6;

                    //same as before
                    llh(p["mass::b(MSbar)"].id());
                    TEST_CHECK_EQUAL(llh.predictions()[1], 5.3);

                    //now observable should be reevaluated
                    llh(p["Abs{c10'}"].id());
                    TEST_CHECK_EQUAL(llh.predictions()[1], 6);
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

                // bootstrap p-value calculation
                {
                    Parameters parameters  = Parameters::Defaults();
                    LogLikelihoodPtr llh(new LogLikelihood(parameters));
                    llh->add(ObservablePtr(new TestObservable(parameters, Kinematics(),
                        "mass::c")), 1.182, 1.192, 1.202);
                    llh->add(ObservablePtr(new TestObservable(parameters, Kinematics(),
                        "mass::c")), 1.19, 1.2, 1.21);

                    parameters["mass::c"] = 1.196;
                    (*llh)();

                    double p_value = llh->bootstrap_p_value(5e4).first;
                    //p-value from chi^2=0.32 and two degrees-of-freedom
                    TEST_CHECK_NEARLY_EQUAL(p_value, 0.852143788, 1e-2);
                }

            }
    } log_likelihood_test;
}
