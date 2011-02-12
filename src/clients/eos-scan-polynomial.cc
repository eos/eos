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

#include <config.h>
#include <src/observable.hh>
#include <src/utils/cartesian-product.hh>
#include <src/utils/chi-squared.hh>
#include <src/utils/destringify.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/one-of.hh>
#include <src/utils/power_of.hh>
#include <src/utils/stringify.hh>
#include <src/utils/wilson-polynomial.hh>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>

using namespace eos;

class DoUsage
{
    private:
        std::string _what;

    public:
        DoUsage(const std::string & what) :
            _what(what)
        {
        }

        const std::string & what() const
        {
            return _what;
        }
};

struct ObservableInput
{
    ObservablePtr observable;

    double min, central, max;
};

struct ObservableRatioInput
{
    ObservablePtr numerator, denominator;

    double min, central, max;
};

typedef OneOf<ObservableInput, ObservableRatioInput> Input;

struct ScanData
{
    std::string name;

    unsigned points;

    double min;

    double max;
};

struct ObservableChiSquared
{
    WilsonPolynomial observable_central;
    std::vector<WilsonPolynomial> observable_max, observable_min;

    double experiment_min, experiment_central, experiment_max;

    double operator() ()
    {
        static WilsonPolynomialEvaluator evaluator;

        double theory_central = observable_central.accept_returning<double>(evaluator);
        double delta_min = 0.0, delta_max = 0.0;

        // Handle parameters with lowered values
        for (auto i = observable_min.cbegin(), i_end = observable_min.cend() ; i != i_end ; ++i)
        {
            double value = i->accept_returning<double>(evaluator);

            if (value < theory_central)
            {
                delta_min += power_of<2>(theory_central - value);
            }
            else if (value > theory_central)
            {
                delta_max += power_of<2>(theory_central - value);
            }
        }

        // Handle parameters with raised values
        for (auto i = observable_max.cbegin(), i_end = observable_max.cend() ; i != i_end ; ++i)
        {
            double value = i->accept_returning<double>(evaluator);

            if (value < theory_central)
            {
                delta_min += power_of<2>(theory_central - value);
            }
            else if (value > theory_central)
            {
                delta_max += power_of<2>(theory_central - value);
            }
        }

        return ChiSquared::with_theory_offset(theory_central - std::sqrt(delta_min), theory_central, theory_central + std::sqrt(delta_max), experiment_min, experiment_central, experiment_max);
    }
};

struct ObservableRatioChiSquared
{
    WilsonPolynomial numerator_central;
    std::vector<WilsonPolynomial> numerator_max, numerator_min;

    WilsonPolynomial denominator_central;
    std::vector<WilsonPolynomial> denominator_max, denominator_min;

    double experiment_min, experiment_central, experiment_max;

    double operator() ()
    {
        static WilsonPolynomialEvaluator evaluator;

        double theory_central = numerator_central.accept_returning<double>(evaluator)
            / denominator_central.accept_returning<double>(evaluator);
        double delta_min = 0.0, delta_max = 0.0;

        // Handle parameters with lowered values
        for (auto n = numerator_min.cbegin(), d = denominator_min.cbegin(), n_end = numerator_min.cend() ; n != n_end ; ++n, ++d)
        {
            double value = n->accept_returning<double>(evaluator) / d->accept_returning<double>(evaluator);

            if (value < theory_central)
            {
                delta_min += power_of<2>(theory_central - value);
            }
            else if (value > theory_central)
            {
                delta_max += power_of<2>(theory_central - value);
            }
        }

        // Handle parameters with raised values
        for (auto n = numerator_max.cbegin(), d = denominator_max.cbegin(), n_end = numerator_max.cend() ; n != n_end ; ++n, ++d)
        {
            double value = n->accept_returning<double>(evaluator) / d->accept_returning<double>(evaluator);

            if (value < theory_central)
            {
                delta_min += power_of<2>(theory_central - value);
            }
            else if (value > theory_central)
            {
                delta_max += power_of<2>(theory_central - value);
            }
        }

        return ChiSquared::with_theory_offset(theory_central - std::sqrt(delta_min), theory_central, theory_central + std::sqrt(delta_max), experiment_min, experiment_central, experiment_max);
    }
};

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        Parameters parameters;

        std::list<ScanData> scans;

        std::list<std::string> coefficients;

        std::list<std::string> variations;

        std::list<Input> inputs;

        CommandLine() :
            parameters(Parameters::Defaults())
        {
        }

        void parse(int argc, char ** argv)
        {
            std::shared_ptr<Kinematics> kinematics(new Kinematics);

            for (char ** a(argv + 1), ** a_end(argv + argc) ; a != a_end ; ++a)
            {
                std::string argument(*a);
                if ("--scan-abs" == argument)
                {
                    std::string coefficient = std::string(*(++a));
                    std::string name = "Abs{" + coefficient + "}";
                    unsigned points = destringify<unsigned>(*(++a));
                    double min = destringify<double>(*(++a));
                    double max = destringify<double>(*(++a));
                    scans.push_back(ScanData{name, points, min, max});
                    if (coefficients.cend() == std::find(coefficients.cbegin(), coefficients.cend(), coefficient))
                        coefficients.push_back(coefficient);

                    continue;
                }
                if ("--scan-arg" == argument)
                {
                    std::string coefficient = std::string(*(++a));
                    std::string name = "Arg{" + coefficient + "}";
                    unsigned points = destringify<unsigned>(*(++a));
                    double min = destringify<double>(*(++a));
                    double max = destringify<double>(*(++a));
                    scans.push_back(ScanData{name, points, min, max});
                    if (coefficients.cend() == std::find(coefficients.cbegin(), coefficients.cend(), coefficient))
                        coefficients.push_back(coefficient);

                    continue;
                }

                if ("--kinematics" == argument)
                {
                    std::string name = std::string(*(++a));
                    double value = destringify<double>(*(++a));
                    kinematics->declare(name);
                    kinematics->set(name, value);

                    continue;
                }

                if ("--observable" == argument)
                {
                    std::string observable_name(*(++a));

                    ObservableInput input;
                    input.observable = Observable::make(observable_name, parameters, *kinematics, Options());
                    if (! input.observable)
                        throw DoUsage("Unknown observable '" + observable_name + "'");

                    input.min = destringify<double>(*(++a));
                    input.central = destringify<double>(*(++a));
                    input.max = destringify<double>(*(++a));

                    inputs.push_back(input);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--ratio" == argument)
                {
                    std::string numerator_name(*(++a)), denominator_name(*(++a));

                    ObservableRatioInput input;
                    input.numerator = Observable::make(numerator_name, parameters, *kinematics, Options());
                    if (! input.numerator)
                        throw DoUsage("Unknown observable '" + numerator_name + "'");

                    input.denominator = Observable::make(denominator_name, parameters, *kinematics, Options());
                    if (! input.numerator)
                        throw DoUsage("Unknown observable '" + denominator_name + "'");

                    input.min = destringify<double>(*(++a));
                    input.central = destringify<double>(*(++a));
                    input.max = destringify<double>(*(++a));

                    inputs.push_back(input);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--vary" == argument)
                {
                    std::string variation_name(*(++a));
                    Parameter p = parameters[variation_name];
                    variations.push_back(variation_name);

                    continue;
                }

                throw DoUsage("Unknown command line argument: " + argument);
            }
        }
};

class WilsonScannerPolynomial
{
    private:
        CartesianProduct<std::vector<double>> _points;

        std::list<Parameter> _scan_parameters;

        std::list<std::function<double ()>> _chi_squared_functions;

    public:
        WilsonScannerPolynomial()
        {
            std::cout << std::scientific;
            std::cout << "# Scan generated by eos-scan-polynomial (" EOS_GITHEAD ")" << std::endl;
            std::cout << "# Coefficients:" << std::endl;
            for (auto c = CommandLine::instance()->coefficients.cbegin(), c_end = CommandLine::instance()->coefficients.cend() ; c != c_end ; ++c)
            {
                std::cout << "#   " << *c << std::endl;
            }

            std::cout << "# Scans:" << std::endl;
            CartesianProduct<std::vector<double>> cp;
            for (auto s = CommandLine::instance()->scans.cbegin(), s_end = CommandLine::instance()->scans.cend() ; s != s_end ; ++s)
            {
                double delta = (s->max - s->min) / s->points;

                std::vector<double> scan_set;
                for (unsigned i = 0 ; i <= s->points ; ++i)
                {
                    scan_set.push_back(s->min + delta * i);
                }

                _points.over(scan_set);
                _scan_parameters.push_back(CommandLine::instance()->parameters[s->name]);

                std::cout << "#   " << s->name << ": [" << s->min << ", " << s->max << "], increment = " << delta << std::endl;
            }

            std::cout << "# Inputs:" << std::endl;
            for (auto i = CommandLine::instance()->inputs.cbegin(), i_end = CommandLine::instance()->inputs.cend() ; i != i_end ; ++i)
            {
                i->accept(*this);
            }

        }

        void visit(const ObservableInput & i)
        {
            Parameters parameters = CommandLine::instance()->parameters;

            std::cout << "#   " << i.observable->name() << '[' << i.observable->kinematics().as_string() << "] = (" << i.min << ", " << i.central << ", " << i.max << ")" << std::endl;

            ObservableChiSquared chi_squared;
            chi_squared.observable_central = make_polynomial(i.observable, CommandLine::instance()->coefficients);
            for (auto v = CommandLine::instance()->variations.cbegin(), v_end = CommandLine::instance()->variations.cend() ; v != v_end ; ++v)
            {
                Parameter p = parameters[*v];
                double old_p = p();

                p = p.min();
                chi_squared.observable_min.push_back(make_polynomial(i.observable, CommandLine::instance()->coefficients));

                p = p.max();
                chi_squared.observable_max.push_back(make_polynomial(i.observable, CommandLine::instance()->coefficients));

                p = old_p;
            }
            chi_squared.experiment_min = i.min;
            chi_squared.experiment_central = i.central;
            chi_squared.experiment_max = i.max;

            _chi_squared_functions.push_back(std::function<double ()>(chi_squared));
        }

        void visit(const ObservableRatioInput & i)
        {
            Parameters parameters = i.numerator->parameters();

            std::cout << "#   " << i.numerator->name() << "[Kinematics]" << " / " << i.denominator->name() << "[Kinematics]" << " = (" << i.min << ", " << i.central << ", " << i.max << ")" << std::endl;

            ObservableRatioChiSquared chi_squared;
            chi_squared.numerator_central = make_polynomial(i.numerator, CommandLine::instance()->coefficients);
            chi_squared.denominator_central = make_polynomial(i.denominator, CommandLine::instance()->coefficients);
            for (auto v = CommandLine::instance()->variations.cbegin(), v_end = CommandLine::instance()->variations.cend() ; v != v_end ; ++v)
            {
                Parameter p = parameters[*v];
                double old_p = p();

                p = p.min();
                chi_squared.numerator_min.push_back(make_polynomial(i.numerator, CommandLine::instance()->coefficients));
                chi_squared.denominator_min.push_back(make_polynomial(i.denominator, CommandLine::instance()->coefficients));

                p = p.max();
                chi_squared.numerator_max.push_back(make_polynomial(i.numerator, CommandLine::instance()->coefficients));
                chi_squared.denominator_max.push_back(make_polynomial(i.denominator, CommandLine::instance()->coefficients));

                p = old_p;
            }
            chi_squared.experiment_min = i.min;
            chi_squared.experiment_central = i.central;
            chi_squared.experiment_max = i.max;

            _chi_squared_functions.push_back(std::function<double ()>(chi_squared));
        }

        void scan()
        {
            for (auto i = _points.begin(), i_end = _points.end() ; i != i_end ; ++i)
            {
                std::vector<double> v = *i;
                std::vector<double> result;

                auto p = _scan_parameters.begin();
                for (auto j = v.begin(), j_end = v.end() ; j != j_end ; ++j, ++p)
                {
                    (*p) = *j;
                    result.push_back(*j);
                }

                double chi_squared = 0.0;
                for (auto c = _chi_squared_functions.cbegin(), c_end = _chi_squared_functions.cend() ; c != c_end ; ++c)
                {
                    chi_squared += (*c)();
                }
                result.push_back(chi_squared);


                auto r = result.cbegin(), r_end = result.cend();
                if (r != r_end)
                    std::cout << *r;

                ++r;
                for ( ; r != r_end ; ++r)
                {
                    std::cout << '\t' << *r;
                }

                std::cout << '\n';
            }
        }
};

int
main(int argc, char * argv[])
{
    try
    {
        CommandLine::instance()->parse(argc, argv);

        if (CommandLine::instance()->inputs.empty())
            throw DoUsage("No inputs specified");

        WilsonScannerPolynomial scanner;
        scanner.scan();
    }
    catch(DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-scan-polynomial" << std::endl;
        std::cout << "  [--vary PARAMETER]*" << std::endl;
        std::cout << "  [[--kinematics NAME VALUE]* --observable NAME MIN CENTRAL MAX]+" << std::endl;
        std::cout << "  [[--scan-abs COEFFICIENT POINTS MIN MAX] | [--scan-arg COEFFICIENT POINTS MIN MAX]]+" << std::endl;
    }
    catch(Exception & e)
    {
        std::cerr << "Caught exception: '" << e.what() << "'" << std::endl;
        return EXIT_FAILURE;
    }
    catch(...)
    {
        std::cerr << "Aborting after unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
