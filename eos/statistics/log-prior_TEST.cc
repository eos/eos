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
#include <eos/statistics/log-prior.hh>
#include <eos/maths/power-of.hh>

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
            static const double mu_0 = 4.18;
            static const double lambda = 2.0;

            // parameter range
            {
                ParameterRange range{ -1.0, 1.0 };

                TEST_CHECK_EQUAL(range.min, -1.0);
                TEST_CHECK_EQUAL(range.max, +1.0);
            }

            const auto inverse_cdf = [](const LogPriorPtr & prior, Parameter & param, const double & p)
            {
                param.set_generator(p);
                prior->sample();
                return param.evaluate();
            };

            // flat prior
            {
                // use factory
                LogPriorPtr flat_prior = LogPrior::Flat(parameters, "mass::b(MSbar)", ParameterRange{ 4.2, 4.5 });
                Parameter param        = parameters["mass::b(MSbar)"];
                TEST_CHECK_NEARLY_EQUAL((*flat_prior)(),                     1.2039728043259361, eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(flat_prior, param, 0.0), 4.2,                eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(flat_prior, param, 0.5), 4.35,               eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(flat_prior, param, 1.0), 4.5,                eps);

                // a continuous parameter of interest
                TEST_CHECK(! flat_prior->begin()->nuisance);
                TEST_CHECK_EQUAL(flat_prior->begin()->parameter->name(), "mass::b(MSbar)");
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
                LogPriorPtr gauss_prior = LogPrior::CurtailedGauss(parameters, "mass::b(MSbar)", ParameterRange{ 4.15, 4.57 },
                        central - sig_lower, central, central + sig_upper);

                parameters["mass::b(MSbar)"] = 4.2;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 0.6555737246958322 /*1.0524382901193807*/, eps);

                parameters["mass::b(MSbar)"] = 4.25;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 1.0305737246958295, eps);

                parameters["mass::b(MSbar)"] = 4.389;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 1.0565612246958278, eps);

                // continuity at zero
                parameters["mass::b(MSbar)"] = 4.3 - 1e-7;
                const double lower_limit = (*gauss_prior)();

                parameters["mass::b(MSbar)"] = 4.3 + 1e-7;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), lower_limit, eps);

                // no CDF checks here: mean and median do not coincide for asymmetric uncertainties
            }

            // cloning
            {
                Parameters independent = Parameters::Defaults();
                LogPriorPtr gauss_prior1 = LogPrior::CurtailedGauss(parameters, "mass::b(MSbar)", ParameterRange{ 4.15, 4.57 },
                        central - sig_lower, central, central+sig_upper);
                LogPriorPtr gauss_prior2 = gauss_prior1->clone(independent);

                parameters["mass::b(MSbar)"] = 4.389;
                independent["mass::b(MSbar)"] = 4.25;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior1)(), 1.0565612246958278, eps);
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior2)(), 1.0305737246958295, eps);
            }

            // vary one sigma interval
            {
                LogPriorPtr gauss_prior = LogPrior::CurtailedGauss(parameters, "mass::b(MSbar)", ParameterRange{ 3.7, 4.9 },
                        4.3, 4.4, 4.5);

                parameters["mass::b(MSbar)"] = 4.2;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), -0.616353153557734281, eps);

                parameters["mass::b(MSbar)"] = 4.3;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 0.883646846442265719, eps);
            }

            // asymmetric
            {
                LogPriorPtr gauss_prior = LogPrior::CurtailedGauss(parameters, "mass::b(MSbar)", ParameterRange{ 0.2, 0.55 },
                            0.319, 0.369, 0.485);

                parameters["mass::b(MSbar)"] = 0.32;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 1.176587791815339, eps);

                parameters["mass::b(MSbar)"] = 0.44;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 1.4694735825406655, eps);
            }

            // Scale prior
            {
                // use factory
                LogPriorPtr scale_prior = LogPrior::Scale(parameters, "mass::b(MSbar)", ParameterRange{ 2.0, 10.0 },
                    mu_0, lambda);

                Parameter param = parameters["mass::b(MSbar)"];

                param = 3.0;
                TEST_CHECK_NEARLY_EQUAL((*scale_prior)(), 0.2404491734814939, eps);

                param = 4.0;
                TEST_CHECK_NEARLY_EQUAL((*scale_prior)(), 0.1803368801111204, eps);

                param = 7.0;
                TEST_CHECK_NEARLY_EQUAL((*scale_prior)(), 0.1030496457777831, eps);

                // central value at p = 0.5
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(scale_prior, param, 0.5), mu_0,          eps);
                // lower boundary at p = 0.0
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(scale_prior, param, 0.0), mu_0 / lambda, eps);
                // upper boundary at p = 1.0
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(scale_prior, param, 1.0), mu_0 * lambda, eps);
            }

            //Make
            {
                Parameters p = Parameters::Defaults();

                std::string s_flat("Parameter: b->smumu::Re{c10}, prior type: flat, range: [-15,15]");
                LogPriorPtr prior_flat = LogPrior::Make(p, s_flat);
                TEST_CHECK_EQUAL(s_flat, prior_flat->as_string());

                std::string s_gauss("Parameter: CKM::A, prior type: Gaussian, range: [0.774,0.834], x = 0.804 +- 0.01");
                LogPriorPtr prior_gauss = LogPrior::Make(p, s_gauss);
                TEST_CHECK_EQUAL(s_gauss, prior_gauss->as_string());
            }
        }
} log_prior_test;
