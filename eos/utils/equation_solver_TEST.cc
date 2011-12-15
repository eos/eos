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

#include <eos/utils/equation_solver.hh>
#include <eos/utils/log.hh>
#include <test/test.hh>

#include <gsl/gsl_sf_gamma.h>

using namespace test;
using namespace eos;

class EquationSolverTest :
    public TestCase
{
    public:
        EquationSolverTest() :
            TestCase("equation_solver_test")
        {
        }

        double simple_constraint(const std::vector<double> & parameters)
        {
            return parameters[0];
        }

        double difficult_constraint(const std::vector<double> & parameters, double central, double x, double cdf)
        {
            const double & lambda = parameters[0];
            const double & alpha  = parameters[1];

            double result = gsl_sf_gamma_inc_Q(alpha, std::exp((x - (central - lambda * std::log(alpha))) / lambda));

            return result - cdf;
        }

        virtual void run() const
        {
            static const double eps = 1e-13;

            // construction
            {
                EquationSolver solver(EquationSolver::Config::Default());
                solver.add("z", 15, 1);
                solver.add(std::bind(&EquationSolverTest::simple_constraint, *this, std::placeholders::_1));

                auto solution = solver.solve();
                TEST_CHECK(solution.valid);
                TEST_CHECK_NEARLY_EQUAL(0.0, solution.value, eps);
                TEST_CHECK_NEARLY_EQUAL(0.0, solution.parameters[0], 1e-6);
            }

            // a more complicated problem
            {
                EquationSolver solver(EquationSolver::Config::Default());
                solver.add("lambda", -7, 2);
                solver.add("alpha", 1.5, 0.2);

                double central = 0.53;
                double sigma_upper = 1.9;
                double sigma_lower = 1.0;

                solver.add(std::bind(&EquationSolverTest::difficult_constraint, *this, std::placeholders::_1, central, central + sigma_upper, 0.84134));
                solver.add(std::bind(&EquationSolverTest::difficult_constraint, *this, std::placeholders::_1, central, central - sigma_lower, 0.158655));

                auto solution = solver.solve();

                TEST_CHECK_NEARLY_EQUAL(0.0, solution.value, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(-2.1677, solution.parameters[0], 1e-4);
                TEST_CHECK(solution.valid);
            }
        }
} equation_solver_test;
