/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
 * Copyright (c) 2011 Danny van Dyk
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
#include <eos/utils/hdf5.hh>
#include <eos/utils/proposal_functions.hh>
#include <eos/utils/power_of.hh>
#include <algorithm>

using namespace test;
using namespace eos;
using namespace eos::proposal_functions;

class ProposalFunctionsTest :
    public TestCase
{
    public:
        ProposalFunctionsTest() :
            TestCase("proposal_functions_test")
        {
        }

        virtual void run() const
        {
            // proposal_functions::sliding_window()
            {
                unsigned j, j_min, j_max;

                static const unsigned K = 2500;
                static const unsigned size = 1000;

                j = 1300;
                proposal_functions::sliding_window(K, size, j, j_min, j_max);
                TEST_CHECK_EQUAL(j_min,  801);
                TEST_CHECK_EQUAL(j_max, 1801);

                j = 2200;
                proposal_functions::sliding_window(K, size, j, j_min, j_max);
                TEST_CHECK_EQUAL(j_min, 1500);
                TEST_CHECK_EQUAL(j_max, 2500);

                j = 200;
                proposal_functions::sliding_window(K, size, j, j_min, j_max);
                TEST_CHECK_EQUAL(j_min, 0);
                TEST_CHECK_EQUAL(j_max, 1000);

                TEST_CHECK_THROWS(InternalError, proposal_functions::sliding_window(size / 2, size, j, j_min, j_max));
            }

            // check single Gaussian proposal function
            {
                MultivariateGaussian ppf(1, std::vector<double>{ 0.01 }, false);

                MarkovChain::State current;
                current.point.push_back(4.35);

                MarkovChain::State proposal;
                proposal.point.push_back(4.0);

                //TODO check new values
                //TEST_CHECK_RELATIVE_ERROR(-4.74135344021063, ppf.evaluate(current, proposal), eps);
                //TEST_CHECK_RELATIVE_ERROR(-4.74135344021063, ppf.evaluate(proposal, current), eps);
            }
            // check multivariate Gaussian proposal function
            {
                auto ppf = proposal_functions::MultivariateGaussian(2,
                        std::vector<double>
                        {
                            0.01,  0.003,
                            0.003, 0.0025
                        },
                        false);

                MarkovChain::State current;
                current.point.push_back(4.3);
                current.point.push_back(1.1);

                MarkovChain::State proposal;
                proposal.point.push_back(4.35);
                proposal.point.push_back(1.2);

                TEST_CHECK_RELATIVE_ERROR(1.30077135, ppf.evaluate(current, proposal), 1e-8);
                TEST_CHECK_RELATIVE_ERROR(1.30077135, ppf.evaluate(proposal, current), 1e-8);
            }

            // check multivariate Gauss with zero efficiency
            {
                // covariance is not invertible
                TEST_CHECK_THROWS(InternalError,
                    proposal_functions::MultivariateGaussian(3,
                    std::vector<double>
                    {
                        0.0,  0.0, 0.0,
                        0.0, 0.0025, 0.0,
                        0.0, 0.0,  0.6
                    },
                    false)
                );

                // negative off-diagonal elements are OK
                    proposal_functions::MultivariateGaussian(3,
                    std::vector<double>
                    {
                        0.1,  0.0, -0.01,
                        0.0, 0.0025, 0.0,
                        -0.01, 0.0,  0.6
                    },
                    false);
            }

            // check multivariate Gauss and StudentT
            {
                // _covariance with zero correlation
                std::vector<double> cov(4, 0.0);
                cov[0] = 0.0049;
                cov[3] = 0.01;

                bool automatic_scaling = false;
                auto mvg = new MultivariateGaussian(2, cov, automatic_scaling);

                double dof = 5;
                auto mvt = new MultivariateStudentT(2, cov, dof, automatic_scaling);

                // evaluate
                MarkovChain::State current;

                current.point = std::vector<double> (1, 1.25);
                current.point.push_back(4.3);

                MarkovChain::State proposal = current;
                proposal.point[0] = 1.3;
                proposal.point[1] = 4.4;

                TEST_CHECK_RELATIVE_ERROR(mvg->evaluate(proposal, current), 2.368866023, 1e-9);
                TEST_CHECK_RELATIVE_ERROR(mvt->evaluate(proposal, current), 2.200202941, 1e-9);

                proposal.point[0] = 1.26;
                proposal.point[1] = 4.424;
                TEST_CHECK_RELATIVE_ERROR(mvg->evaluate(proposal, current), 2.344963982, 1e-9);
                TEST_CHECK_RELATIVE_ERROR(mvt->evaluate(proposal, current), 2.174596526, 1e-9);
            }

            // cloning
            {
                // _covariance with zero correlation
                std::vector<double> cov(4, 0.0);
                cov[0] = 0.0049;
                cov[3] = 0.01;

                auto mvg = new MultivariateGaussian(2, cov);
                mvg->cooling_power = 3.0;
                mvg->covariance()->data[0] = 0.0050;

                double dof = 5;
                auto mvt = new MultivariateStudentT(2, cov, dof);
                mvt->adaptations = 8;
                mvt->covariance()->data[1] = mvt->covariance()->data[2] = 0.0001;

                // don't use get() on temporary!
                ProposalFunctionPtr p_mvg_clone = mvg->clone();
                ProposalFunctionPtr p_mvt_clone = mvt->clone();

                auto mvg_clone = dynamic_cast<proposal_functions::MultivariateGaussian *>(p_mvg_clone.get());
                auto mvt_clone = dynamic_cast<proposal_functions::MultivariateStudentT *>(p_mvt_clone.get());

                TEST_CHECK_EQUAL(mvg->cooling_power, mvg_clone->cooling_power);
                TEST_CHECK_EQUAL(mvg->covariance()->data[0], mvg_clone->covariance()->data[0]);
                TEST_CHECK_EQUAL(mvg->covariance()->data[3], mvg_clone->covariance()->data[3]);
                TEST_CHECK_EQUAL(mvg->adaptations, mvg_clone->adaptations);

                TEST_CHECK_EQUAL(mvt->cooling_power, mvt_clone->cooling_power);
                TEST_CHECK_EQUAL(mvt->covariance()->data[0], mvt_clone->covariance()->data[0]);
                TEST_CHECK_EQUAL(mvt->covariance()->data[1], mvt_clone->covariance()->data[1]);
                TEST_CHECK_EQUAL(mvt->covariance()->data[3], mvt_clone->covariance()->data[3]);
                TEST_CHECK_EQUAL(mvt->adaptations, mvt_clone->adaptations);
            }

            // check sampling of single gaussian
            {
                auto ppf = proposal_functions::MultivariateGaussian(1, std::vector<double>{ 0.01 }, false);

                gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, 1243);

                double mean = 0.0, chi_squared = 0.0;

                MarkovChain::State current;
                current.point.push_back(4.3);

                MarkovChain::State proposal;
                proposal.point.push_back(0.0);

                static const unsigned N = 50000;
                for (unsigned i = 0 ; i < N ; ++i)
                {
                    ppf.propose(proposal, current, rng);
                    mean += proposal.point[0];
                    chi_squared += power_of<2>(proposal.point[0] - 4.3) / 0.01;
                }
                mean /= N;
                chi_squared /= N;

                TEST_CHECK_RELATIVE_ERROR(mean,        4.3, 0.001);
                TEST_CHECK_RELATIVE_ERROR(chi_squared, 1.0, 0.005);

                gsl_rng_free(rng);
            }

            // check sampling of multivariate gaussian
            {
                auto ppf = proposal_functions::MultivariateGaussian(2,
                        std::vector<double>
                        {
                            0.01,  0.003,
                            0.003, 0.0025
                        },
                        false);

                gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, 1243);

                std::vector<double> mean{ 0.0, 0.0 };
                double chi_squared = 0.0;

                MarkovChain::State current;
                current.point.push_back(4.3);
                current.point.push_back(1.1);

                MarkovChain::State proposal;
                proposal.point.push_back(0.0);
                proposal.point.push_back(0.0);

                double normalization = -log(2 * M_PI) - 0.5 * log(1.6e-5);

                static const unsigned N = 50000;
                for (unsigned i = 0 ; i < N ; ++i)
                {
                    ppf.propose(proposal, current, rng);
                    mean[0] += proposal.point[0];
                    mean[1] += proposal.point[1];
                    chi_squared += -2.0 * (ppf.evaluate(proposal, current) - normalization);
                }
                mean[0] /= N;
                mean[1] /= N;
                chi_squared /= N;

                TEST_CHECK_RELATIVE_ERROR(mean[0],     4.3, 0.001);
                TEST_CHECK_RELATIVE_ERROR(mean[1],     1.1, 0.001);
                TEST_CHECK_RELATIVE_ERROR(chi_squared, 2.0, 0.005);

                gsl_rng_free(rng);
            }

            // read and write to HDF5
            {
                static const std::string file_name = EOS_BUILDDIR "/eos/utils/proposal_functions_TEST-rdwr.hdf5";
                {
                    hdf5::File file = hdf5::File::Create(file_name);
                }

                /* Gaussian */
                {
                    // covariance() with zero correlation
                    std::vector<double> cov(4, 0.0);
                    cov[0] = 0.0049;
                    cov[3] = 0.01;

                    auto mvg = new MultivariateGaussian(2, cov);

                    // write
                    hdf5::File file = hdf5::File::Open(file_name, H5F_ACC_RDWR);
                    mvg->dump_state(file, "/Gaussian/");

                    // and read back in
                    ProposalFunctionPtr p_mvg_copy = Factory::make(file, "/Gaussian", "MultivariateGaussian", 2);
                    auto mvg_copy = dynamic_cast<proposal_functions::MultivariateGaussian *>(p_mvg_copy.get());

                    TEST_CHECK_EQUAL(mvg->cooling_power, mvg_copy->cooling_power);
                    TEST_CHECK_EQUAL(mvg->adaptations, mvg_copy->adaptations);
                    TEST_CHECK_EQUAL(mvg->covariance()->data[0], mvg_copy->covariance()->data[0]);
                    TEST_CHECK_EQUAL(mvg->covariance()->data[1], mvg_copy->covariance()->data[2]);
                    TEST_CHECK_EQUAL(mvg->covariance()->data[3], mvg_copy->covariance()->data[3]);
                }

                /* Student T */
                {
                    std::vector<double> cov(4, 0.0);
                    cov[0] = 0.0049;
                    cov[3] = 0.01;

                    double dof = 2.423;

                    auto mvt = new MultivariateStudentT(2, cov, dof);

                    // write
                    hdf5::File file = hdf5::File::Open(file_name, H5F_ACC_RDWR);
                    mvt->dump_state(file, "/StudentT");
                    // and read back in
                    ProposalFunctionPtr p_mvt_copy = Factory::make(file, "/StudentT", "MultivariateStudentT", 2);
                    auto mvt_copy = dynamic_cast<proposal_functions::MultivariateStudentT *>(p_mvt_copy.get());

                    TEST_CHECK_EQUAL(mvt->cooling_power, mvt_copy->cooling_power);
                    TEST_CHECK_EQUAL(mvt->adaptations, mvt_copy->adaptations);
                    TEST_CHECK_EQUAL(mvt->dof, mvt_copy->dof);
                    TEST_CHECK_EQUAL(mvt->covariance()->data[0], mvt_copy->covariance()->data[0]);
                    TEST_CHECK_EQUAL(mvt->covariance()->data[1], mvt_copy->covariance()->data[2]);
                    TEST_CHECK_EQUAL(mvt->covariance()->data[3], mvt_copy->covariance()->data[3]);
                }
            }

            // BlockDecomposition
            {
                std::vector<double> cov(4, 0.0);
                cov[0] = 0.0049;
                cov[3] = 0.01;
                MultivariateProposalPtr mv(new MultivariateGaussian(2, cov));

                Parameters p = Parameters::Defaults();
                LogPriorPtr flat = LogPrior::Flat(p, "mass::c", ParameterRange{ 1.0, 3.0 });

                static const std::string file_name(EOS_BUILDDIR "/eos/utils/proposal_functions_TEST-block-decomposition.hdf5");

                // one Multivariate
                {
                    BlockDecomposition bd;

                    bd.add(mv);

                    MarkovChain::State current;
                    current.point = std::vector<double>{4.3, 1.1};

                    MarkovChain::State proposal;
                    proposal.point = std::vector<double>{4.2, 1.13};

                    TEST_CHECK_RELATIVE_ERROR(bd.evaluate(proposal, current), mv->evaluate(proposal, current), 1e-15);

                    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
                    gsl_rng_set(rng, 1346);

                    auto proposal_bd = proposal;

                    bd.propose(proposal_bd, current, rng);

                    gsl_rng_set(rng, 1346);
                    mv->propose(proposal, current, rng);

                    TEST_CHECK_EQUAL(proposal_bd.point[0], proposal.point[0]);
                    TEST_CHECK_EQUAL(proposal_bd.point[1], proposal.point[1]);
                }
                //one prior
                {
                    BlockDecomposition bd;

                    bd.add(flat);

                    MarkovChain::State current;
                    current.point = std::vector<double>{ 1.1 };

                    MarkovChain::State proposal;
                    proposal.point = std::vector<double>{ 1.13 };

                    TEST_CHECK_NEARLY_EQUAL(bd.evaluate(proposal, current), std::log(0.5), 1e-15);
                }
                // multivariate and prior
                {
                    BlockDecomposition bd;

                    bd.add(mv);
                    bd.add(flat);

                    MarkovChain::State current;
                    current.point = std::vector<double>{4.3, 1.1, 1.8};

                    MarkovChain::State proposal;
                    proposal.point = std::vector<double>{4.2, 1.13, 1.46};

                    MarkovChain::State current_2d;
                    current_2d.point = std::vector<double>{4.3, 1.1};

                    MarkovChain::State proposal_2d;
                    proposal_2d.point = std::vector<double>{4.2, 1.13};

                    TEST_CHECK_NEARLY_EQUAL(bd.evaluate(proposal, current),
                                            mv->evaluate(proposal_2d, current_2d) + std::log(0.5), 1e-15);

                    auto file = hdf5::File::Create(file_name);
                    bd.dump_state(file, "/block decomposition");
                }
                // read/write
                {
                    auto file = hdf5::File::Open(file_name);

                    ProposalFunctionPtr bd = Factory::make(file, "/block decomposition", "BlockDecomposition", 3);

                    MarkovChain::State current;
                    current.point = std::vector<double>{4.3, 1.1, 1.8};

                    MarkovChain::State proposal;
                    proposal.point = std::vector<double>{4.2, 1.13, 1.46};

                    MarkovChain::State current_2d;
                    current_2d.point = std::vector<double>{4.3, 1.1};

                    MarkovChain::State proposal_2d;
                    proposal_2d.point = std::vector<double>{4.2, 1.13};

                    TEST_CHECK_NEARLY_EQUAL(bd->evaluate(proposal, current),
                                            mv->evaluate(proposal_2d, current_2d) + std::log(0.5), 1e-15);
                }
            }
        }
} proposal_functions_test;
