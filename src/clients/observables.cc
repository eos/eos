/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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
#include <iostream>
#include <list>

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
        kinematics.declare("s");
        double s_low(0.0), s_high(19.21);
        unsigned points = 50;
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
                std::cerr << "Parameter: " << name << " -> " << value << std::endl;
                continue;
            }

            if ("--range" == argument)
            {
                s_low = destringify<double>(*(++a));
                s_high = destringify<double>(*(++a));
                std::cerr << "Range: " << s_low << " .. " << s_high << std::endl;
                continue;
            }

            if ("--points" == argument)
            {
                points = destringify<unsigned>(*(++a));
                std::cerr << "Points: " << points << std::endl;
                continue;
            }

            if ("--budget" == argument)
            {
                std::string name(*(++a));
                budgets.push_back(std::make_tuple(name, std::list<Parameter>()));
                std::cerr << "Budget: " << name << std::endl;
                continue;
            }

            if ("--vary" == argument)
            {
                if (budgets.empty())
                    throw DoUsage("Specify a budget before any variations");

                std::string name(*(++a));
                std::get<1>(budgets.back()).push_back(parameters[name]);
                std::cerr << "Vary: " << name << std::endl;
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

        std::cout << "## Observables ##" << std::endl;
        for (auto o(observables.begin()), o_end(observables.end()) ; o != o_end ; ++o)
        {
            std::cout << "# " << (*o)->name() << std::endl;
        }

        std::cout << "## Data ##" << std::endl;
        for (auto o(observables.begin()), o_end(observables.end()) ; o != o_end ; ++o)
        {
            std::cout << "# " << (*o)->name() << std::endl;
            for (unsigned j = 0 ; j <= points ; ++j)
            {
                double s = s_low + j * (s_high - s_low) / points;

                std::cout << s << std::flush;
                kinematics.set("s", s);

                double central = (*o)->evaluate();
                double delta_min = 0.0, delta_max = 0.0;

                std::cout << '\t' << central << std::flush;

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

                    std::cout << '\t' << std::sqrt(budget_min) << '\t' << std::sqrt(budget_max) << std::flush;
                }

                delta_max = std::sqrt(delta_max);
                delta_min = std::sqrt(delta_min);

                std::cout << '\t' << delta_min << '\t' << delta_max << std::endl;
            }

            std::cout << std::endl << std::endl;
        }
    }
    catch (DoUsage & e)
    {
        std::cerr << e.what << std::endl;
        std::cerr << "Usage: observables --range SMIN SMAX [--points N] [--parameter NAME VALUE]* [--vary NAME]* [--observable NAME]+" << std::endl;
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
    catch (...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
