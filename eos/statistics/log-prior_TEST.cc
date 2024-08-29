/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
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

#include <test/test.hh>
#include <eos/constraint.hh>
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

            const auto cdf = [](const LogPriorPtr & prior, Parameter & param, const double & x)
            {
                param.set(x);
                prior->compute_cdf();
                return param.evaluate_generator();
            };

            const auto inverse_cdf = [](const LogPriorPtr & prior, Parameter & param, const double & p)
            {
                param.set_generator(p);
                prior->sample();
                return param.evaluate();
            };

            // flat prior
            {
                // use factory
                LogPriorPtr flat_prior = LogPrior::Flat(parameters, "mass::b(MSbar)", 4.2, 4.5 );
                Parameter param        = parameters["mass::b(MSbar)"];
                TEST_CHECK_NEARLY_EQUAL((*flat_prior)(),                     1.2039728043259361, eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(flat_prior, param, 0.0), 4.2,                eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(flat_prior, param, 0.5), 4.35,               eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(flat_prior, param, 1.0), 4.5,                eps);

                // a continuous parameter of interest
                TEST_CHECK_EQUAL(flat_prior->begin()->name(), "mass::b(MSbar)");
            }


            /*
             * curtailed gaussian prior
             * compare with
             *   import scipy.stats
             *   a = scipy.stats.norm.cdf(4.57, loc=4.3, scale=0.2) - 0.5
             *   b = -scipy.stats.norm.cdf(4.15, loc=4.3, scale=0.1) + 0.5
             *   c = 1/(a+b) #overall normalization to have \int p_0(x) =1 over allowed range
             *   log(c*scipy.stats.norm.pdf(4.389, loc=4.3, scale=0.2))
             */
            {
                // use factory
                LogPriorPtr gauss_prior = LogPrior::CurtailedGauss(parameters, "mass::b(MSbar)", 4.15, 4.57,
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
                LogPriorPtr gauss_prior1 = LogPrior::CurtailedGauss(parameters, "mass::b(MSbar)", 4.15, 4.57,
                        central - sig_lower, central, central+sig_upper);
                LogPriorPtr gauss_prior2 = gauss_prior1->clone(independent);

                parameters["mass::b(MSbar)"] = 4.389;
                independent["mass::b(MSbar)"] = 4.25;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior1)(), 1.0565612246958278, eps);
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior2)(), 1.0305737246958295, eps);
            }

            // vary one sigma interval
            {
                LogPriorPtr gauss_prior = LogPrior::CurtailedGauss(parameters, "mass::b(MSbar)", 3.7, 4.9,
                        4.3, 4.4, 4.5);

                parameters["mass::b(MSbar)"] = 4.2;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), -0.616353153557734281, eps);

                parameters["mass::b(MSbar)"] = 4.3;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 0.883646846442265719, eps);
            }

            // asymmetric
            {
                LogPriorPtr gauss_prior = LogPrior::CurtailedGauss(parameters, "mass::b(MSbar)", 0.2, 0.55,
                            0.319, 0.369, 0.485);

                parameters["mass::b(MSbar)"] = 0.32;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 1.176587791815339, eps);

                parameters["mass::b(MSbar)"] = 0.44;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 1.4694735825406655, eps);
            }

            // Scale prior
            {
                // use factory
                LogPriorPtr scale_prior = LogPrior::Scale(parameters, "mass::b(MSbar)", 2.0, 10.0,
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

            // Gaussian prior (w/ infinite support)
            {
                // use factory
                static const double mu    = 0.5;
                static const double sigma = 0.25;
                LogPriorPtr gaussian_prior = LogPrior::Gaussian(parameters, "mass::b(MSbar)", mu, sigma);

                Parameter param = parameters["mass::b(MSbar)"];

                param = -0.5;
                TEST_CHECK_NEARLY_EQUAL((*gaussian_prior)(), -7.532644172085, eps);

                param =  0.5;
                TEST_CHECK_NEARLY_EQUAL((*gaussian_prior)(), +0.467355827915, eps);

                param = +1.5;
                TEST_CHECK_NEARLY_EQUAL((*gaussian_prior)(), -7.532644172085, eps);

                // CDF
                TEST_CHECK_NEARLY_EQUAL(cdf(gaussian_prior, param, -0.5), 0.000031671242, eps);
                TEST_CHECK_NEARLY_EQUAL(cdf(gaussian_prior, param, +0.5), 0.500000000000, eps);
                TEST_CHECK_NEARLY_EQUAL(cdf(gaussian_prior, param, +1.5), 0.999968328758, eps);

                // inverse CDF
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(gaussian_prior, param, 0.005), -0.1439573258872, eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(gaussian_prior, param, 0.250),  0.3313775624510, eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(gaussian_prior, param, 0.900),  0.8203878913862, eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(gaussian_prior, param, 0.950),  0.9112134067379, eps);
            }

            // Multivariate Gaussian
            {
                // look up constraint entry
                auto entry = Constraints()["B->K::FormFactors[parametric,LCSR]@GKvD:2018A"];
                // create prior
                auto prior = entry->make_prior(parameters, Options());

                // set generator values to 0.5 and sample
                parameters["B->K::alpha^f+_0@BSZ2015"].set_generator(0.5);
                parameters["B->K::alpha^f+_1@BSZ2015"].set_generator(0.5);
                parameters["B->K::alpha^f+_2@BSZ2015"].set_generator(0.5);
                parameters["B->K::alpha^f0_1@BSZ2015"].set_generator(0.5);
                parameters["B->K::alpha^f0_2@BSZ2015"].set_generator(0.5);
                parameters["B->K::alpha^fT_0@BSZ2015"].set_generator(0.5);
                parameters["B->K::alpha^fT_1@BSZ2015"].set_generator(0.5);
                parameters["B->K::alpha^fT_2@BSZ2015"].set_generator(0.5);
                prior->sample();

                // check that the parameters are at the median
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^f+_0@BSZ2015"], +0.2655528728950212, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^f+_1@BSZ2015"], -0.6466140804171657, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^f+_2@BSZ2015"], -0.1337677825631754, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^f0_1@BSZ2015"], +0.3841222758978433, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^f0_2@BSZ2015"], -0.6628825163091753, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^fT_0@BSZ2015"], +0.2510192324158927, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^fT_1@BSZ2015"], -0.6508680050905388, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^fT_2@BSZ2015"], +0.0999901466869552, eps);

                // set generator values to 0.5 and sample
                parameters["B->K::alpha^f+_0@BSZ2015"].set_generator(0.15865525393145702);
                parameters["B->K::alpha^f+_1@BSZ2015"].set_generator(0.15865525393145702);
                parameters["B->K::alpha^f+_2@BSZ2015"].set_generator(0.15865525393145702);
                parameters["B->K::alpha^f0_1@BSZ2015"].set_generator(0.15865525393145702);
                parameters["B->K::alpha^f0_2@BSZ2015"].set_generator(0.15865525393145702);
                parameters["B->K::alpha^fT_0@BSZ2015"].set_generator(0.15865525393145702);
                parameters["B->K::alpha^fT_1@BSZ2015"].set_generator(0.15865525393145702);
                parameters["B->K::alpha^fT_2@BSZ2015"].set_generator(0.15865525393145702);
                prior->sample();

                // check that the parameters match the reference implementation
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^f+_0@BSZ2015"], +0.185453816796, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^f+_1@BSZ2015"], -0.802422750729, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^f+_2@BSZ2015"], +2.133896762012, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^f0_1@BSZ2015"], +0.013676201228, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^f0_2@BSZ2015"], -0.697392274325, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^fT_0@BSZ2015"], +0.162651391807, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^fT_1@BSZ2015"], -0.790991231378, eps);
                TEST_CHECK_NEARLY_EQUAL(parameters["B->K::alpha^fT_2@BSZ2015"], +4.223260161888, eps);
            }

            // Poisson prior: k = 10.0
            {
                // use factory
                static const double k = 10.0;

                LogPriorPtr poisson_prior = LogPrior::Poisson(parameters, "mass::b(MSbar)", k);

                Parameter param = parameters["mass::b(MSbar)"];

                param = 3.0 / k;
                TEST_CHECK_NEARLY_EQUAL((*poisson_prior)(), -4.815704593400,  eps);

                param = 7.0 / k;
                TEST_CHECK_NEARLY_EQUAL((*poisson_prior)(), -0.3427259895280, eps);

                param = 10.0 / k;
                TEST_CHECK_NEARLY_EQUAL((*poisson_prior)(),  0.2240234498590, eps);

                param = 20.0 / k;
                TEST_CHECK_NEARLY_EQUAL((*poisson_prior)(), -2.844504744542,  eps);

                // CDF
                TEST_CHECK_NEARLY_EQUAL(cdf(poisson_prior, param,  5 / k), 0.013695268598, eps);
                TEST_CHECK_NEARLY_EQUAL(cdf(poisson_prior, param,  9 / k), 0.294011679659, eps);
                TEST_CHECK_NEARLY_EQUAL(cdf(poisson_prior, param, 11 / k), 0.540111297306, eps);
                TEST_CHECK_NEARLY_EQUAL(cdf(poisson_prior, param, 20 / k), 0.989188281173, eps);

                // inverse CDF
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(poisson_prior, param, 0.05) * k,  6.169007289395323, eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(poisson_prior, param, 0.25) * k,  8.619809702379529, eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(poisson_prior, param, 0.90) * k, 15.40664117197652,  eps);
                TEST_CHECK_NEARLY_EQUAL(inverse_cdf(poisson_prior, param, 0.95) * k, 16.9622192357219,   eps);
            }

            // Transform prior
            {
                // use factory
                std::vector<double> shift     = {0.0,0.0};
                std::vector<std::vector<double>> transform = {{ 0.707106, 0.707106 }, {-0.707106, 0.707106}};
                std::vector<double> min       = {-2.0,-2.0};
                std::vector<double> max       = {2.0,2.0};
                LogPriorPtr transform_prior = LogPrior::Transform(parameters, {"scnuee::Re{cVL}","scnuee::Re{cVR}"}, shift, transform, min, max);

                parameters["scnuee::Re{cVL}"] = 0.0;
                parameters["scnuee::Re{cVR}"] = 0.0;
                TEST_CHECK_NEARLY_EQUAL((*transform_prior)(), -2.77259, 1.0e-5);

                parameters["scnuee::Re{cVL}"] = -3.0;
                parameters["scnuee::Re{cVR}"] = 0.0;
                TEST_CHECK(! std::isfinite((*transform_prior)()));
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
