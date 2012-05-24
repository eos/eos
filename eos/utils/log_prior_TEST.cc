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
#include <eos/utils/log_prior.hh>
#include <eos/utils/histogram.hh>
#include <eos/utils/power_of.hh>

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
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), -0.616353153557734281, eps);

                parameters["mass::b(MSbar)"] = 4.3;
                TEST_CHECK_NEARLY_EQUAL((*gauss_prior)(), 0.883646846442265719, eps);
            }

            //test sampling
            {
                ParameterRange range { 4.15, 4.57 };
                LogPriorPtr gauss_prior = LogPrior::Gauss(parameters, "mass::b(MSbar)", range, 4.2, 4.3, 4.5);

                gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, 1243);

                /* symmetric Gaussian with huge range */
                range  = ParameterRange{ 3, 6 };
                gauss_prior = LogPrior::Gauss(parameters, "mass::b(MSbar)", range, 4.2, 4.3, 4.4);

                gsl_rng_set(rng, 1243);

                unsigned n_bins = 80;
                Histogram<1> hist = Histogram<1>::WithEqualBinning(3.8, 4.8, n_bins);

                double chi_sq = 0.0;
                unsigned N = 1e5;
                for (unsigned i = 0 ; i < N ; ++i)
                {
                    double x = gauss_prior->sample(rng);
                    hist.insert(x);
                    chi_sq += power_of<2>(x - 4.3) / 0.01;
                }

                std::vector<double> pdf_vec;

                for (auto b = hist.begin(), b_end = hist.end() ; b != b_end ; ++b)
                {
                    pdf_vec.push_back(b->value);
                }

                auto cumulative_hist = estimate_cumulative_distribution(hist);

                std::vector<double> cumulative_vec;
                for (auto b = cumulative_hist.begin(), b_end = cumulative_hist.end() ; b != b_end ; ++b)
                {
                    cumulative_vec.push_back(b->value);
                }

                TEST_CHECK_RELATIVE_ERROR(chi_sq / N, 1.0, 1.1e-2);

                TEST_CHECK_RELATIVE_ERROR(pdf_vec[unsigned(0.5 * n_bins) - 1], pdf_vec[unsigned(0.5 * n_bins)], 1.4e-2);

                // median in the middle
                TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(0.5 * n_bins) - 1], 0.5,  1e-2);

                // one sigma interval should be [4.2, 4.4]
                TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(0.6 * n_bins) - 1] - cumulative_vec[unsigned(0.4 * n_bins) - 1], 0.68,  5e-2);

                // two sigma interval should be [4.1, 4.5]
                TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(0.7 * n_bins) - 1] - cumulative_vec[unsigned(0.3 * n_bins) - 1], 0.95,  1.4e-2);

                gsl_rng_free(rng);
            }

                /* symmetric Gaussian with small range */
            {
                ParameterRange range  = ParameterRange{ 4.1, 4.5 };
                LogPriorPtr gauss_prior = LogPrior::Gauss(parameters, "mass::b(MSbar)", range, 4.2, 4.3, 4.4);

                gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, 1243);

                unsigned n_bins = 80;
                Histogram<1> hist = Histogram<1>::WithEqualBinning(3.8, 4.8, n_bins);

                double chi_sq = 0.0;
                unsigned N = 1e5;
                for (unsigned i = 0 ; i < N ; ++i)
                {
                    double x = gauss_prior->sample(rng);
                    hist.insert(x);
                    chi_sq += power_of<2>(x - 4.3) / 0.01;
                }

                std::vector<double> pdf_vec;

                for (auto b = hist.begin(), b_end = hist.end() ; b != b_end ; ++b)
                {
                    pdf_vec.push_back(b->value);
                }

                auto cumulative_hist = estimate_cumulative_distribution(hist);

                std::vector<double> cumulative_vec;
                for (auto b = cumulative_hist.begin(), b_end = cumulative_hist.end() ; b != b_end ; ++b)
                {
                    cumulative_vec.push_back(b->value);
                }

                // account for finite range: only a fraction of usual Gaussian is contained
                static const double scale_factor = 1/ 0.9545;

                TEST_CHECK_RELATIVE_ERROR(pdf_vec[unsigned(0.5 * n_bins) - 1], pdf_vec[unsigned(0.5 * n_bins)], 1.4e-2);

                // median in the middle
                TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(0.5 * n_bins) - 1], 0.5,  1e-2);

                // one sigma interval should be [4.2, 4.4]
                TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(0.6 * n_bins) - 1] - cumulative_vec[unsigned(0.4 * n_bins) - 1], 0.68 * scale_factor,  5e-2);

                // two sigma interval should be [4.1, 4.5]
                TEST_CHECK_RELATIVE_ERROR(cumulative_vec[unsigned(0.7 * n_bins) - 1] - cumulative_vec[unsigned(0.3 * n_bins) - 1], 0.95 * scale_factor,  1.4e-2);

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

            // LogGamma prior
            {
                static const double low_eps = 5e-5;

                LogPriorPtr log_gamma = LogPrior::LogGamma(parameters, "mass::b(MSbar)", ParameterRange{ -10, 10 }, 0.34, 0.53, 0.63);

                parameters["mass::b(MSbar)"] = 0.57;
                TEST_CHECK_RELATIVE_ERROR((*log_gamma)(), +1.005543554, low_eps);
                parameters["mass::b(MSbar)"] = 0.3;
                TEST_CHECK_RELATIVE_ERROR((*log_gamma)(), +0.1737006298, low_eps);

                // lead to trouble when alpha was allowed to be negative
                // so make sure ctor doesn't throw
                log_gamma = LogPrior::LogGamma(parameters, "mass::c", ParameterRange{ 1.09, 1.41 }, +1.18, +1.27, +1.34);

                // bad starting point?
                log_gamma = LogPrior::LogGamma(parameters, "B->K::F^p(0)@KMPW2010", ParameterRange{ 0.3, 0.44 }, 0.32, 0.34, 0.39);

                // barely soluble
                log_gamma = LogPrior::LogGamma(parameters, "mass::c", ParameterRange{ -10, 10 }, -1, 0, +1.03);

                // completely symmetric
                TEST_CHECK_THROWS(InternalError, LogPrior::LogGamma(parameters, "mass::c", ParameterRange{ 1.09, 1.41 }, +1.18, +1.27, +1.36));

                // almost symmetric
                TEST_CHECK_THROWS(InternalError, LogPrior::LogGamma(parameters, "mass::c", ParameterRange{ 1.09, 1.41 }, +0.99, +1.00, +1.01002));

                // shouldn't throw
                LogPrior::LogGamma(parameters, "mass::c", ParameterRange{ 1.09, 1.41 }, +0.99, +1.00, +1.0108);

                // back to useful example
                log_gamma = LogPrior::LogGamma(parameters, "mass::b(MSbar)", ParameterRange{ -10, 10 }, 0.66, 1.20, 1.8);

                parameters["mass::b(MSbar)"] = 0.57;
                TEST_CHECK_NEARLY_EQUAL((*log_gamma)(), -1.05802105, 4e-4);
                parameters["mass::b(MSbar)"] = 0.92;
                TEST_CHECK_NEARLY_EQUAL((*log_gamma)(), -0.4842376414, 4e-4);

                // cloning
                auto p = Parameters::Defaults();
                p["mass::b(MSbar)"] = 1.7;
                auto clone = log_gamma->clone(p);
                TEST_CHECK((*log_gamma)() != (*clone)());
                p["mass::b(MSbar)"] = 0.92;
                TEST_CHECK_EQUAL((*log_gamma)(), (*clone)());

                // sampling
                gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, 2022);

                double mean = 0.0;
                unsigned N = 1e4;
                for (unsigned i = 0 ; i < N ; ++i)
                {
                    mean += log_gamma->sample(rng);
                }

                mean /= double(N);

                // mean = \nu + \lambda * \psi(\alpha)
                TEST_CHECK_RELATIVE_ERROR(mean, 1.28979, 3e-3);

                gsl_rng_free(rng);

                // finite range increases pdf
                log_gamma = LogPrior::LogGamma(parameters, "mass::b(MSbar)", ParameterRange{ 0.2, 0.7 }, 0.34, 0.53, 0.63);
                parameters["mass::b(MSbar)"] = 0.57;
                TEST_CHECK_RELATIVE_ERROR((*log_gamma)(), 1.139778733, low_eps);
            }

            //Make
            {
                Parameters p = Parameters::Defaults();

                std::string s_flat("Parameter: Re{c10}, prior type: flat, range: [-15,15]");
                LogPriorPtr prior_flat = LogPrior::Make(p, s_flat);
                TEST_CHECK_EQUAL(s_flat, prior_flat->as_string());

                std::string s_gauss("Parameter: CKM::A, prior type: Gaussian, range: [0.774,0.834], x = 0.804 +- 0.01");
                LogPriorPtr prior_gauss = LogPrior::Make(p, s_gauss);
                TEST_CHECK_EQUAL(s_gauss, prior_gauss->as_string());

                std::string s_gamma("Parameter: mass::c, prior type: LogGamma, range: [1,1.48], x = 1.27 + 0.07 - 0.09, nu: 1.201633227, lambda: 0.1046632085, alpha: 1.921694426");
                LogPriorPtr prior_gamma1 = LogPrior::Make(p, s_gamma);
                LogPriorPtr prior_gamma2 = LogPrior::LogGamma(p, "mass::c", ParameterRange{ 1.0, 1.48 }, 1.18, 1.27, 1.34);

                // string comparison should work
                TEST_CHECK_EQUAL(s_gamma, prior_gamma1->as_string());
                TEST_CHECK_EQUAL(s_gamma, prior_gamma2->as_string());
            }
        }
} log_prior_test;
