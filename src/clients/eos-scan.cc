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
#include <eos/utils/cartesian-product.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/lock.hh>
#include <eos/utils/log.hh>
#include <eos/utils/mutex.hh>
#include <eos/utils/thread_pool.hh>

#include <cmath>
#include <config.h>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <utility>
#include <vector>

using namespace eos;

struct Input
{
        const double min;

        const double max;

        const double o_min;

        const double o;

        const double o_max;

        const std::string o_name;
};

struct ScanData
{
        std::string name;
        unsigned    points;
        double      min;
        double      max;
};

class WilsonScan
{
    public:
        Mutex * const mutex;

        std::list<std::pair<Input, ObservablePtr>> bins;

        std::list<ScanData> scan_data;

        std::list<Input> inputs;

        std::list<std::string> variation_names;

        double theory_uncertainty;

        std::list<std::pair<std::vector<double>, double>> results;

        WilsonScan(const std::list<ScanData> & scan_data, const std::list<Input> & inputs, const std::list<std::pair<std::string, double>> & param_changes,
                   const std::list<std::string> & variation_names, const double & theory_uncertainty) :
            mutex(new Mutex),
            scan_data(scan_data),
            inputs(inputs),
            variation_names(variation_names),
            theory_uncertainty(theory_uncertainty)
        {
            Parameters parameters = Parameters::Defaults();
            Kinematics kinematics;
            kinematics.declare("s_min");
            kinematics.declare("s_max");

            for (const auto & param_change : param_changes)
            {
                parameters[param_change.first] = param_change.second;
            }

            for (const auto & input : inputs)
            {
                bins.push_back(std::make_pair(input, Observable::make(input.o_name, parameters.clone(), kinematics.clone(), Options())));
            }
        }

        void
        calc_chi_square(const Input & input, const ObservablePtr & observable, const CartesianProduct<std::vector<double>>::Iterator & wc_iterator)
        {
            Kinematics k = observable->kinematics();
            k.set("s_min", input.min);
            k.set("s_max", input.max);

            ObservablePtr o      = observable->clone();
            Parameters    params = o->parameters();

            auto                sd        = scan_data.cbegin();
            std::vector<double> wc_values = *wc_iterator;
            for (auto w = wc_values.cbegin(); wc_values.cend() != w; ++w, ++sd)
            {
                params[sd->name] = *w;
            }

            // Calculate chi^2
            double central   = o->evaluate();
            double delta_min = 0.0, delta_max = 0.0;
            for (auto & variation_name : variation_names)
            {
                Parameter p     = params[variation_name];
                double    old_p = p();
                double    max = 0.0, min = 0.0, value;

                p     = p.min();
                value = o->evaluate();
                if (value > central)
                {
                    max = value - central;
                }

                if (value < central)
                {
                    min = central - value;
                }

                p     = p.max();
                value = o->evaluate();
                if (value > central)
                {
                    max = std::max(max, value - central);
                }

                if (value < central)
                {
                    min = std::max(min, central - value);
                }

                p = old_p;

                delta_min += min * min;
                delta_max += max * max;
            }

            delta_min += pow(central * theory_uncertainty, 2);
            delta_max += pow(central * theory_uncertainty, 2);

            delta_max = std::sqrt(delta_max);
            delta_min = std::sqrt(delta_min);

            double chi = 0.0;
            if (input.o - central > delta_max)
            {
                chi = input.o - central - delta_max;
            }
            else if (central - input.o > delta_min)
            {
                chi = central - input.o - delta_min;
            }

            chi                /= (input.o_max - input.o_min);
            double chi_squared  = chi * chi;

            {
                Lock l(*mutex);
                results.push_back(std::make_pair(wc_values, chi_squared));
            }
        }

        void
        scan()
        {
            std::cout << "# Generated by eos-scan (" EOS_GITHEAD ")" << std::endl;
            std::cout << "# Scan data" << std::endl;

            CartesianProduct<std::vector<double>> cp;
            for (auto sd = scan_data.cbegin(); scan_data.cend() != sd; ++sd)
            {
                std::vector<double> set;
                double              delta = (sd->max - sd->min) / sd->points;

                for (unsigned i = 0; i <= sd->points; ++i)
                {
                    set.push_back(sd->min + delta * i);
                }
                cp.over(set);

                std::cout << "#   " << sd->name << ": [" << sd->min << ", " << sd->max << "], increment = " << delta << std::endl;
            }

            std::cout << "# Inputs" << std::endl;
            for (const auto & input : inputs)
            {
                std::cout << "#   " << input.o_name << "[" << input.min << ", " << input.max << "]" << " = (" << input.o_min << ", " << input.o << ", " << input.o_max << ")"
                          << std::endl;
            }

            TicketList    tickets;
            unsigned long jobs = 0;
            for (auto bin = bins.begin(); bins.end() != bin; ++bin)
            {
                for (auto w = cp.begin(); cp.end() != w; ++w)
                {
                    ThreadPool::instance()->wait_for_free_capacity();
                    tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&WilsonScan::calc_chi_square, this, bin->first, bin->second, w)));
                    ++jobs;
                    if (jobs % 100 == 0)
                    {
                        std::cerr << '[' << jobs << '/' << cp.size() << ']' << std::endl;
                    }
                }
            }

            tickets.wait();

            std::cout << std::scientific << std::setprecision(7);
            for (const auto & result : results)
            {
                for (auto w = result.first.cbegin(), w_end = result.first.cend(); w != w_end; ++w)
                {
                    std::cout << *w << '\t';
                }

                std::cout << result.second << std::endl;
            }
        }
};

class DoUsage
{
    private:
        std::string _what;

    public:
        DoUsage(const std::string & what) :
            _what(what)
        {
        }

        const std::string &
        what() const
        {
            return _what;
        }
};

int
main(int argc, char * argv[])
{
    try
    {
        std::list<ScanData>                       scan_data;
        std::list<Input>                          input;
        std::list<std::string>                    variation_names;
        std::list<std::pair<std::string, double>> param_changes;
        double                                    theory_uncertainty = 0.0;

        Log::instance()->set_program_name("eos-scan");

        for (char **a(argv + 1), **a_end(argv + argc); a != a_end; ++a)
        {
            std::string argument(*a);
            if ("--scan" == argument)
            {
                std::string name   = std::string(*(++a));
                unsigned    points = destringify<unsigned>(*(++a));
                double      min    = destringify<double>(*(++a));
                double      max    = destringify<double>(*(++a));
                scan_data.push_back(ScanData{ name, points, min, max });
                continue;
            }

            if ("--parameter" == argument)
            {
                std::string name  = std::string(*(++a));
                double      value = destringify<double>(*(++a));
                param_changes.push_back(std::make_pair(name, value));
                continue;
            }

            if ("--input" == argument)
            {
                std::string observable(*(++a));
                double      k1      = destringify<double>(*(++a));
                double      k2      = destringify<double>(*(++a));
                double      min     = destringify<double>(*(++a));
                double      central = destringify<double>(*(++a));
                double      max     = destringify<double>(*(++a));

                input.push_back(Input{ k1, k2, min, central, max, observable });

                continue;
            }

            if ("--vary" == argument)
            {
                std::string variation_name(*(++a));
                variation_names.push_back(variation_name);

                continue;
            }

            if ("--theory-uncertainty" == argument)
            {
                theory_uncertainty = destringify<double>(*(++a));

                continue;
            }

            throw DoUsage("Unknown command line argument: " + argument);
        }

        if (scan_data.empty())
        {
            throw DoUsage("Need at least one scan parameter");
        }

        if (input.empty())
        {
            throw DoUsage("Need at least one input");
        }

        WilsonScan scanner(scan_data, input, param_changes, variation_names, theory_uncertainty);
        scanner.scan();
    }
    catch (DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-scan" << std::endl;
        std::cout << "  [--vary PARAMETER]*" << std::endl;
        std::cout << "  [--parameter NEWVALUE]*" << std::endl;
        std::cout << "  [--input NAME SMIN SMAX MIN CENTRAL MAX]+" << std::endl;
        std::cout << "  [--scan PARAMETER POINTS MIN MAX]+" << std::endl;
        std::cout << "  [--theory-uncertainty PERCENT]" << std::endl;
    }
    catch (Exception & e)
    {
        std::cerr << "Caught exception: '" << e.what() << "'" << std::endl;
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cerr << "Aborting after unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
