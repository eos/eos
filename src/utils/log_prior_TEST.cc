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

#include <test/test.hh>
#include <src/utils/log_prior.hh>

using namespace test;
using namespace eos;

class LogPriorTest :
    public TestCase
{
    public:
        LogPriorTest() :
            TestCase("log_prior_test")
        {
        }

        virtual void run() const
        {
            Parameters parameters = Parameters::Defaults();

            static const double eps = 1e-12;
            static const double central = 4.3;
            static const double sig_lower = 0.1;
            static const double sig_upper = 0.2;

            // parameter range
            {
                ParameterRange range{ -1.0, 1.0 };

                TEST_CHECK_EQUAL(range.min, -1.0);
                TEST_CHECK_EQUAL(range.max, +1.0);
            }

            // flat prior
            {
                // use factory
                LogPriorPtr flat_prior = LogPrior::Flat(parameters, "mass::b(MSbar)", ParameterRange{ 4.2, 4.5 });
                TEST_CHECK_NEARLY_EQUAL((*flat_prior)(), 1.2039728043259361, eps);

                //a continuous parameter of interest
                TEST_CHECK( ! flat_prior->begin()->nuisance);
                TEST_CHECK( ! flat_prior->begin()->discrete);
                TEST_CHECK_EQUAL(flat_prior->begin()->parameter.name(), "mass::b(MSbar)");
            }


            /*
             * gaussian prior
             * compare with
             *   import scipy.stats
             *   a = scipy.stats.norm.cdf(4.57, loc=4.3, scale=0.2) - 0.5
             *   b = -scipy.stats.norm.cdf(4.15, loc=4.3, scale=0.1) + 0.5
             *   c = 1/(a+b) #overall normalization to have \int p_0(x) =1 over allowed range
             *   log(c*scipy.stats.norm.pdf(4.389, loc=4.3, scale=0.2))
             */
            {
                // use factory
                LogPriorPtr gauss_prior = LogPrior::Gauss(parameters, "mass::b(MSbar)", ParameterRange{ 4.15, 4.57 },
                        central - sig_lower, central, central+sig_upper);

                // set m_b = 4.2
                parameters["mass::b(MSbar)"] = 4.2;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 1.0524382901193807, eps);

                //set m_b = 4.25
                parameters["mass::b(MSbar)"] = 4.25;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 1.4274382901193783, eps);

                //set m_b = 4.389
                parameters["mass::b(MSbar)"] = 4.389;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 0.76027860955943105, eps);
            }

            // cloning
            {
                Parameters independent = Parameters::Defaults();
                LogPriorPtr gauss_prior1 = LogPrior::Gauss(parameters, "mass::b(MSbar)", ParameterRange{ 4.15, 4.57 },
                        central - sig_lower, central, central+sig_upper);
                LogPriorPtr gauss_prior2 = gauss_prior1->clone(independent);

                parameters["mass::b(MSbar)"] = 4.389;
                independent["mass::b(MSbar)"] = 4.25;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior1)(), 0.76027860955943105, eps);
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior2)(), 1.42743829011937830, eps);
            }

            // vary one sigma interval
            {
                LogPriorPtr gauss_prior = LogPrior::Gauss(parameters, "mass::b(MSbar)", ParameterRange{ 3.7, 4.9 },
                        4.3, 4.4, 4.5);

                parameters["mass::b(MSbar)"] = 4.2;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), -0.6163531535577218, eps);

                parameters["mass::b(MSbar)"] = 4.3;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), +0.883646846442260436, eps);
            }

            //test sampling
            {
                LogPriorPtr gauss_prior = LogPrior::Gauss(parameters, "mass::b(MSbar)", ParameterRange{ 4.15, 4.57 },
                                        4.2, 4.3, 4.5);

                gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, 1243);

                //calculate manually what will happen for the first three numbers.

                TEST_CHECK_NEARLY_EQUAL(gauss_prior -> sample(rng), 4.559828404943599, eps );
                TEST_CHECK_NEARLY_EQUAL(gauss_prior -> sample(rng), 4.390994101050554, eps );
                TEST_CHECK_NEARLY_EQUAL(gauss_prior -> sample(rng), 4.286086713649191, eps );

                gsl_rng_free(rng);
            }

            //test discrete prior
            {
                LogPriorPtr discrete = LogPrior::Discrete(Parameters::Defaults(), "Arg{c9}", std::set<double>{0, M_PI} );

                // spits out 1 1 0 1 1 on range {0,1}
                gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, 1243);

                TEST_CHECK_EQUAL(M_PI, discrete->sample(rng));
                TEST_CHECK_EQUAL(M_PI, discrete->sample(rng));
                TEST_CHECK_EQUAL(0   , discrete->sample(rng));
                TEST_CHECK_EQUAL(M_PI, discrete->sample(rng));
                TEST_CHECK_EQUAL(M_PI, discrete->sample(rng));

                TEST_CHECK_EQUAL(discrete->begin()->parameter.name(), "Arg{c9}");
                TEST_CHECK(!discrete->begin()->nuisance);
                TEST_CHECK(discrete->begin()->discrete);

                LogPriorPtr discrete2 = LogPrior::Discrete(Parameters::Defaults(), "Arg{c10'}", std::set<double>{0, 0.75, 1, 0.25, 0.5} );

                // spits out indices 4 0 1 4 1 on range {0...4}
                gsl_rng* rng2 = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, 984);

                // note that samples are ordered according to size
                // so element[4] is the largest value
                TEST_CHECK_EQUAL(1   , discrete2->sample(rng2));
                TEST_CHECK_EQUAL(0   , discrete2->sample(rng2));
                TEST_CHECK_EQUAL(0.25, discrete2->sample(rng2));
                TEST_CHECK_EQUAL(1   , discrete2->sample(rng2));
                TEST_CHECK_EQUAL(0.25, discrete2->sample(rng2));

                gsl_rng_free(rng);
                gsl_rng_free(rng2);
            }
        }
} log_prior_test;
