/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
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

#include <eos/observable.hh>
#include <eos/utils/destringify.hh>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <list>
#include <tuple>

using namespace eos;

struct DoUsage
{
    std::string what;

    DoUsage(const std::string & what) :
        what(what)
    {
    }
};

int
main(int argc, char * argv[])
{
    try
    {
        Parameters parameters(Parameters::Defaults());
        Kinematics kinematics;
        kinematics.declare("s_min");
        kinematics.declare("s_max");

        double s_min(0.0), s_max(0.0);
        std::list<std::tuple<std::string, std::list<Parameter>>> budgets;
        std::list<ObservablePtr> observables;

        for (char ** a(argv + 1), ** a_end(argv + argc) ; a != a_end ; ++a)
        {
            std::string argument(*a);
            if ("--parameter" == argument)
            {
                std::string name(*(++a));
                double value(destringify<double>(*(++a)));
                parameters.set(name, value);
                continue;
            }

            if ("--range" == argument)
            {
                s_min = destringify<double>(*(++a));
                s_max = destringify<double>(*(++a));
                std::cerr << "Range: " << s_min << " .. " << s_max << std::endl;
                continue;
            }

            if ("--budget" == argument)
            {
                std::string name(*(++a));
                budgets.push_back(std::make_tuple(name, std::list<Parameter>()));
                continue;
            }

            if ("--vary" == argument)
            {
                if (budgets.empty())
                    throw DoUsage("Specify a budget before any variations");

                std::string name(*(++a));
                std::get<1>(budgets.back()).push_back(parameters[name]);
                continue;
            }

            if ("--observable" != argument)
                throw DoUsage("Unknown option: '" + argument + "'");

            std::string name(*(++a));
            ObservablePtr ptr(Observable::make(name, parameters, kinematics, Options()));
            if (! ptr)
                throw DoUsage("Unknown observable: '" + name + "'");

            observables.push_back(ptr);
        }

        if (observables.empty())
            throw DoUsage("Need at least one observable");

        kinematics.set("s_min", s_min);
        kinematics.set("s_max", s_max);

        for (auto o(observables.begin()), o_end(observables.end()) ; o != o_end ; ++o)
        {
            double central = (*o)->evaluate();
            double delta_min = 0.0, delta_max = 0.0;

            std::list<std::tuple<std::string, double, double>> uncertainties;

            for (auto b(budgets.begin()), b_end(budgets.end()) ; b != b_end ; ++b)
            {
                std::string name = std::get<0>(*b);
                std::list<Parameter> variations = std::get<1>(*b);

                double budget_max = 0.0, budget_min = 0.0;

                for (auto p(variations.begin()), p_end(variations.end()) ; p != p_end ; ++p)
                {
                    double old_p = *p;
                    double max = 0.0, min = 0.0, value;

                    *p = p->min();
                    value = (*o)->evaluate();
                    if (value > central)
                        max = value - central;

                    if (value < central)
                        min = central - value;

                    *p = p->max();
                    value = (*o)->evaluate();
                    if (value > central)
                        max = std::max(max, value - central);

                    if (value < central)
                        min = std::max(min, central - value);

                    *p = old_p;

                    delta_min += min * min;
                    delta_max += max * max;

                    budget_min += min * min;
                    budget_max += max * max;
                }

                uncertainties.push_back(std::make_tuple(name, std::sqrt(budget_min), std::sqrt(budget_max)));
            }

            delta_max = std::sqrt(delta_max);
            delta_min = std::sqrt(delta_min);

            std::cout << std::setprecision(7) << std::scientific;
            std::cout << (*o)->name() << '[' << s_min << ".." << s_max << "]:" << std::endl;
            std::cout << '\t'
                << central << ' '
                << '-' << delta_min << ' '
                << '+' << delta_max << '\t'
                << '-' << std::abs(delta_min / central) * 100 << " %  "
                << '+' << std::abs(delta_max / central) * 100 << " %  "
                << std::endl;

            if (uncertainties.empty())
                continue;

            std::cout << "Budgets:" << std::endl;
            for (auto b(uncertainties.begin()) ; b != uncertainties.end() ; ++b)
            {
                std::cout
                    << '\t' << '-' << std::get<1>(*b) << " +" << std::get<2>(*b)
                    << " [" << std::get<0>(*b) << "] " << std::endl;
            }
            std::cout << std::endl;
        }
    }
    catch (DoUsage & e)
    {
        std::cerr << e.what << std::endl;
        std::cerr << "Usage: observables --range SMIN SMAX [--parameter NAME VALUE]* [--vary NAME]* [--observable NAME]+" << std::endl;
        return EXIT_FAILURE;
    }
    catch (Exception & e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (std::exception & e)
    {
        std::cerr << "STL Exception; " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (std::string & e)
    {
        std::cerr << "String Exception: " << e << std::endl;
        return EXIT_FAILURE;
    }
    catch (...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
