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

#ifndef EOS_GUARD_EOS_UTILS_EQUATION_SOLVER_HH
#define EOS_GUARD_EOS_UTILS_EQUATION_SOLVER_HH 1

#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/verify.hh>

#include <functional>
#include <vector>

namespace eos
{
    class EquationSolver :
        public PrivateImplementationPattern<EquationSolver>
    {
        public:
            typedef std::function<double (const std::vector<double> &)> Equation;

            struct Config;
            struct Solution;

            EquationSolver(const EquationSolver::Config & config);

            ~EquationSolver();

            void add(const Equation & equation);

            void add(const std::string & name, const double & initial_value, const double & error);

            void add(const std::string & name, const double & initial_value, const double & error, const double & min, const double & max);

            Solution solve();
    };

    struct EquationSolver::Config
    {
        private:
            Config();

        public:
            /// Maximum number of evaluations
            unsigned maximum_steps;

            /// [0 = low, 1 = medium, 2 = high] precision
            VerifiedRange<unsigned> strategy;

            /// Required precision.
            VerifiedRange<double> tolerance;

            static Config Default();
    };

    struct EquationSolver::Solution
    {
        /// The parameter values giving the best solution to the problem.
        std::vector<double> parameters;

        /// Value of the function to be minimized at the minimum found.
        double value;

        /// True if the solution finding succeeded.
        bool valid;
    };
}

#endif
