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
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>


#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>

namespace eos
{
    namespace minuit
    {
        struct SolverAdapter :
            public ROOT::Minuit2::FCNBase
        {
            // each equation of form f(x1, x2, ... ) = 0
            std::vector<EquationSolver::Equation> equations;

            ROOT::Minuit2::MnUserParameters user_parameters;

            std::shared_ptr<ROOT::Minuit2::FunctionMinimum> minimum;

            virtual double
            Up() const
            {
                return 0.5;
            }

            // minimum at f = 0. Add up squared terms, each one with
            // the minimum at zero
            virtual double
            operator()(const std::vector<double> & parameter_values) const
            {
                double result = 0.0;
                double squared_sum = 0.0;

                for (auto e = equations.begin(), e_end = equations.end() ; e != e_end ; ++e)
                {
                    result = (*e)(parameter_values);
                    squared_sum += power_of<2>(result);
                }

                return squared_sum;
            }

            virtual ~SolverAdapter()
            {
            }
        };
    }

    template <>
    struct Implementation<EquationSolver>
    {
        EquationSolver::Config config;

        minuit::SolverAdapter mn;

        Implementation(const EquationSolver::Config & config) :
            config(config)
        {
        }

        void
        add(const EquationSolver::Equation & equation)
        {
            mn.equations.push_back(equation);
        }

        void
        add(const std::string & name, const double & initial_value, const double & error)
        {
            mn.user_parameters.Add(name, initial_value, error);
        }

        void
        add(const std::string & name, const double & initial_value, const double & error, const double & min, const double & max)
        {
            mn.user_parameters.Add(name, initial_value, error, min, max);
        }

        EquationSolver::Solution
        solve()
        {
            // create MIGRAD minimizer with strategy
            ROOT::Minuit2::MnMigrad migrad(mn, mn.user_parameters, config.strategy);

            // minimize and save results
            ROOT::Minuit2::FunctionMinimum minimum(migrad(config.maximum_steps, config.tolerance));

            // copy values
            EquationSolver::Solution solution { minimum.UserParameters().Params(), minimum.Fval(), minimum.IsValid() };

            return solution;
        }
    };

    EquationSolver::EquationSolver(const EquationSolver::Config & config) :
        PrivateImplementationPattern<EquationSolver>(new Implementation<EquationSolver>(config))
    {
    }

    EquationSolver::~EquationSolver()
    {
    }

    void
    EquationSolver::add(const Equation & equation)
    {
        _imp->add(equation);
    }

    void
    EquationSolver::add(const std::string & name, const double & initial_value, const double & error)
    {
        _imp->add(name, initial_value, error);
    }

    void
    EquationSolver::add(const std::string & name, const double & initial_value, const double & error, const double & min, const double & max)
    {
        _imp->add(name, initial_value, error, min, max);
    }

    EquationSolver::Solution
    EquationSolver::solve()
    {
        return _imp->solve();
    }

    EquationSolver::Config::Config() :
        maximum_steps(2000),
        strategy(0, 2, 2),
        tolerance(0.0, 1.0, 1e-12)
    {
    }

    EquationSolver::Config
    EquationSolver::Config::Default()
    {
        return Config();
    }
}
