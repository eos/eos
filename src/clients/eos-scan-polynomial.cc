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
#include <eos/observable.hh>
#include <eos/utils/cartesian-product.hh>
#include <eos/utils/chi-squared.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/one-of.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/scan_file.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/thread_pool.hh>
#include <eos/utils/wilson-polynomial.hh>

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

struct ObservableHTLikeRatioInput
{
    ObservablePtr numerator, denominator1, denominator2;

    double min, central, max;
};

typedef OneOf<ObservableInput, ObservableRatioInput, ObservableHTLikeRatioInput> Input;

struct ScanData
{
    std::string name;

    unsigned points;

    double min;

    double max;
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

        std::string output;

        std::string creator;

        double theory_uncertainty;

        CommandLine() :
            parameters(Parameters::Defaults()),
            theory_uncertainty(0.0)
        {
        }

        void parse(int argc, char ** argv)
        {
            Log::instance()->set_program_name("eos-scan-polynomial");

            std::shared_ptr<Kinematics> kinematics(new Kinematics);

            creator = std::string(argv[0]);
            for (int i = 1 ; i < argc ; ++i)
            {
                creator += ' ' + std::string(argv[i]);
            }

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
                    std::string name(*(++a));

                    ObservableInput input;
                    input.min = destringify<double>(*(++a));
                    input.central = destringify<double>(*(++a));
                    input.max = destringify<double>(*(++a));
                    input.observable = Observable::make(name, parameters, *kinematics, Options());
                    if (! input.observable)
                        throw DoUsage("Unknown observable '" + name + "'");

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

                if ("--ht-like-ratio" == argument)
                {
                    std::string numerator_name(*(++a)), denominator1_name(*(++a)), denominator2_name(*(++a));

                    ObservableHTLikeRatioInput input;
                    input.numerator = Observable::make(numerator_name, parameters, *kinematics, Options());
                    if (! input.numerator)
                        throw DoUsage("Unknown observable '" + numerator_name + "'");

                    input.denominator1 = Observable::make(denominator1_name, parameters, *kinematics, Options());
                    if (! input.numerator)
                        throw DoUsage("Unknown observable '" + denominator1_name + "'");

                    input.denominator2 = Observable::make(denominator2_name, parameters, *kinematics, Options());
                    if (! input.numerator)
                        throw DoUsage("Unknown observable '" + denominator2_name + "'");

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

                if ("--output" == argument)
                {
                    output = std::string(*(++a));
                    continue;
                }

                if ("--parameter" == argument)
                {
                    std::string name(*(++a));
                    double value = destringify<double>(*(++a));

                    parameters[name] = value;

                    continue;
                }

                if ("--theory-uncertainty" == argument)
                {
                    theory_uncertainty = destringify<double>(*(++a));

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

        std::vector<Parameter> _scan_parameters;

        std::vector<Parameter> _variations;

        std::vector<std::tuple<ObservablePtr, double, double, double, std::vector<std::tuple<ObservablePtr, ObservablePtr>>>> _observables;

        std::vector<Ticket> _tickets;

        ScanFile _output;

        std::vector<ScanFile::DataSet> _data_sets;

    public:
        WilsonScannerPolynomial() :
            _output(ScanFile::Create(CommandLine::instance()->output, "eos-scan-polynomial"))
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

            std::cout << "# Variations:" << std::endl;
            for (auto v = CommandLine::instance()->variations.cbegin(), v_end = CommandLine::instance()->variations.cend() ; v != v_end ; ++v)
            {
                Parameter p = CommandLine::instance()->parameters[*v];
                std::cout << "#   " << p.name() << ": " << p.min() << " < " << p << " < " << p.max() << std::endl;
                _variations.push_back(p);
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

            ObservablePtr observable = make_polynomial_observable(make_polynomial(i.observable, CommandLine::instance()->coefficients), parameters);
            std::vector<std::tuple<ObservablePtr, ObservablePtr>> varied_observables;
            for (auto v = _variations.begin(), v_end = _variations.end() ; v != v_end ; ++v)
            {
                double old_v = *v;

                *v = v->max();
                ObservablePtr raised = make_polynomial_observable(make_polynomial(i.observable, CommandLine::instance()->coefficients), parameters);

                *v = v->min();
                ObservablePtr lowered = make_polynomial_observable(make_polynomial(i.observable, CommandLine::instance()->coefficients), parameters);

                *v = old_v;

                varied_observables.push_back(std::make_tuple(raised, lowered));
            }
            _observables.push_back(std::make_tuple(observable, i.min, i.central, i.max, varied_observables));
        }

        void visit(const ObservableRatioInput & i)
        {
            Parameters parameters = CommandLine::instance()->parameters;

            std::cout << "#   " << i.numerator->name() << "[Kinematics]" << " / " << i.denominator->name() << "[Kinematics]" << " = (" << i.min << ", " << i.central << ", " << i.max << ")" << std::endl;

            ObservablePtr observable = make_polynomial_ratio(make_polynomial(i.numerator, CommandLine::instance()->coefficients),
                    make_polynomial(i.denominator, CommandLine::instance()->coefficients),
                    parameters);

            std::vector<std::tuple<ObservablePtr, ObservablePtr>> varied_observables;
            for (auto v = _variations.begin(), v_end = _variations.end() ; v != v_end ; ++v)
            {
                double old_v = *v;

                *v = v->max();
                ObservablePtr raised = make_polynomial_ratio(make_polynomial(i.numerator, CommandLine::instance()->coefficients),
                        make_polynomial(i.denominator, CommandLine::instance()->coefficients),
                        parameters);

                *v = v->min();
                ObservablePtr lowered = make_polynomial_ratio(make_polynomial(i.numerator, CommandLine::instance()->coefficients),
                        make_polynomial(i.denominator, CommandLine::instance()->coefficients),
                        parameters);

                *v = old_v;

                varied_observables.push_back(std::make_tuple(raised, lowered));
            }
            _observables.push_back(std::make_tuple(observable, i.min, i.central, i.max, varied_observables));
        }

        void visit(const ObservableHTLikeRatioInput & i)
        {
            Parameters parameters = CommandLine::instance()->parameters;

            std::cout << "#   " << i.numerator->name() << "[Kinematics]"
                << " / Sqrt(" << i.denominator1->name() << " * " << i.denominator2->name() << ")"
                << "[Kinematics]" << " = (" << i.min << ", " << i.central << ", " << i.max << ")"
                << std::endl;

            ObservablePtr observable = make_polynomial_ht_like_ratio(make_polynomial(i.numerator, CommandLine::instance()->coefficients),
                    make_polynomial(i.denominator1, CommandLine::instance()->coefficients),
                    make_polynomial(i.denominator2, CommandLine::instance()->coefficients),
                    parameters);

            std::vector<std::tuple<ObservablePtr, ObservablePtr>> varied_observables;
            for (auto v = _variations.begin(), v_end = _variations.end() ; v != v_end ; ++v)
            {
                double old_v = *v;

                *v = v->max();
                ObservablePtr raised = make_polynomial_ht_like_ratio(make_polynomial(i.numerator, CommandLine::instance()->coefficients),
                        make_polynomial(i.denominator1, CommandLine::instance()->coefficients),
                        make_polynomial(i.denominator2, CommandLine::instance()->coefficients),
                        parameters);

                *v = v->min();
                ObservablePtr lowered = make_polynomial_ht_like_ratio(make_polynomial(i.numerator, CommandLine::instance()->coefficients),
                        make_polynomial(i.denominator1, CommandLine::instance()->coefficients),
                        make_polynomial(i.denominator2, CommandLine::instance()->coefficients),
                        parameters);

                *v = old_v;

                varied_observables.push_back(std::make_tuple(raised, lowered));
            }
            _observables.push_back(std::make_tuple(observable, i.min, i.central, i.max, varied_observables));
        }

        void scan_range(const CartesianProduct<std::vector<double>>::Iterator & begin, const CartesianProduct<std::vector<double>>::Iterator & end,
                unsigned index)
        {
            // Clone
            Parameters parameters = CommandLine::instance()->parameters.clone();

            std::vector<Parameter> scan_parameters;
            for (auto p = _scan_parameters.cbegin(), p_end = _scan_parameters.cend() ; p != p_end ; ++p)
            {
                scan_parameters.push_back(parameters[p->name()]);
            }

            std::vector<std::tuple<ObservablePtr, double, double, double, std::vector<std::tuple<ObservablePtr, ObservablePtr>>>> observables;
            for (auto o = _observables.cbegin(), o_end = _observables.cend() ; o != o_end ; ++o)
            {
                std::vector<std::tuple<ObservablePtr, ObservablePtr>> varied_observables;
                for (auto v = std::get<4>(*o).begin(), v_end = std::get<4>(*o).end() ; v != v_end ; ++v)
                {
                    varied_observables.push_back(std::make_tuple(std::get<0>(*v)->clone(parameters), std::get<1>(*v)->clone(parameters)));
                }

                observables.push_back(std::make_tuple(std::get<0>(*o)->clone(parameters), std::get<1>(*o), std::get<2>(*o), std::get<3>(*o), varied_observables));
            }

            // Allocate a write buffer
            ScanFile::WriteBuffer buffer(_scan_parameters.size() + 1);

            // Scan our range
            for (auto i = begin, i_end = end ; i != i_end ; ++i)
            {
                std::vector<double> v = *i;
                std::vector<double> result;

                auto p = scan_parameters.begin();
                for (auto j = v.begin(), j_end = v.end() ; j != j_end ; ++j, ++p)
                {
                    (*p) = *j;
                    result.push_back(*j);
                }

                // Calculate chi^2
                double chi_squared = 0.0;
                for (auto o = observables.begin(), o_end = observables.end() ; o != o_end ; ++o)
                {
                    ObservablePtr observable = std::get<0>(*o);
                    double central = observable->evaluate();
                    double delta_min = 0.0, delta_max = 0.0;
                    std::vector<std::tuple<ObservablePtr, ObservablePtr>> varied_observables = std::get<4>(*o);

                    for (auto v = varied_observables.begin(), v_end = varied_observables.end() ; v != v_end ; ++v)
                    {
                        double max = 0.0, min = 0.0, value;

                        // Handle parameters with lowered values
                        value = std::get<1>(*v)->evaluate();

                        if (value > central)
                            max = value - central;

                        if (value < central)
                            min = central - value;

                        // Handle parameters with raised values
                        value = std::get<0>(*v)->evaluate();

                        if (value > central)
                            max = std::max(max, value - central);

                        if (value < central)
                            min = std::max(min, central -value);

                        delta_min += min * min;
                        delta_max += max * max;
                    }

                    delta_min += power_of<2>(central * CommandLine::instance()->theory_uncertainty);
                    delta_max += power_of<2>(central * CommandLine::instance()->theory_uncertainty);

                    chi_squared += ChiSquared::with_theory_offset(central - std::sqrt(delta_min), central, central + std::sqrt(delta_max),
                            std::get<1>(*o), std::get<2>(*o), std::get<3>(*o));
                }

                result.push_back(chi_squared);

                buffer << result;

                if (buffer.capacity() == buffer.size())
                {
                    _data_sets[index] << buffer;
                    buffer.clear();
                }
            }

            _data_sets[index] << buffer;
            buffer.clear();
        }

        void scan()
        {
            // Find chunk size for N chunks
            unsigned chunk_size = _points.size() / ThreadPool::instance()->number_of_threads();

            // Enqueue N - 1 jobs
            auto c = _points.begin();
            unsigned i = 0;
            for ( ; i < ThreadPool::instance()->number_of_threads() - 1 ; ++i, c += chunk_size)
            {
                auto begin = c, end = c;
                end += chunk_size;

                _data_sets.push_back(_output.add("chunk #" + stringify(i), _scan_parameters.size() + 1));

                auto f = _data_sets.back().begin_fields();
                for (auto p = _scan_parameters.cbegin(), p_end = _scan_parameters.cend() ; p != p_end ; ++p, ++f)
                {
                    f->name(p->name());
                }
                f->name("posterior");

                _tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&WilsonScannerPolynomial::scan_range, this, begin, end, i)));
            }

            // Enqueue an extra job for the remains
            _data_sets.push_back(_output.add("chunk #" + stringify(i), _scan_parameters.size() + 1));
            auto f = _data_sets.back().begin_fields();
            for (auto p = _scan_parameters.cbegin(), p_end = _scan_parameters.cend() ; p != p_end ; ++p, ++f)
            {
                f->name(p->name());
            }
            f->name("posterior");
            _tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&WilsonScannerPolynomial::scan_range, this, c, _points.end(), i)));

            // Wait for job completion
            for (auto t = _tickets.begin(), t_end = _tickets.end() ; t != t_end ; ++t)
            {
                t->wait();
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
            throw DoUsage("Need to specify at least one input!");

        if (CommandLine::instance()->output.empty())
            throw DoUsage("Need to specify output!");

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
