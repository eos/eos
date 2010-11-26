/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/factory.hh>
#include <src/utils/destringify.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/one-of.hh>
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

    Kinematics kinematics;
};

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        Parameters parameters;

        std::list<std::string> coefficients;

        std::list<ObservableInput> inputs;

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
                if ("--coefficient" == argument)
                {
                    std::string coefficient = std::string(*(++a));
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
                    input.kinematics = *kinematics;
                    kinematics.reset(new Kinematics);
                    input.observable = RareBFactory::make(observable_name, parameters, ObservableOptions());
                    if (! input.observable)
                        throw DoUsage("Unknown observable '" + observable_name + "'");

                    inputs.push_back(input);

                    continue;
                }

                throw DoUsage("Unknown command line argument: " + argument);
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
            throw DoUsage("No input specified");

        WilsonPolynomialEvaluator evaluator;
        WilsonPolynomialPrinter printer;
        for (auto i = CommandLine::instance()->inputs.cbegin(), i_end = CommandLine::instance()->inputs.cend() ; i != i_end ; ++i)
        {
            WilsonPolynomial polynomial = make_polynomial(i->observable, i->kinematics, CommandLine::instance()->coefficients);

            std::cout << i->observable->name() + "[";

            auto c = CommandLine::instance()->coefficients.cbegin(), c_end = CommandLine::instance()->coefficients.cend();
            if (c != c_end)
            {
                std::cout << *c;
                ++c;
            }

            for ( ; c != c_end ; ++c)
            {
                std::cout << ", " << *c;
            }
            std::cout << "] = " << polynomial.accept_returning<std::string>(printer) << std::endl;
            std::cout << "polynomial = " << polynomial.accept_returning<double>(evaluator) << std::endl;
            std::cout << "direct     = " << i->observable->evaluate(i->kinematics) << std::endl;
        }
    }
    catch(DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-print-polynomial" << std::endl;
        std::cout << "  [--coefficient WILSONCOEFFICIENT]*" << std::endl;
        std::cout << "  [[--kinematics NAME VALUE]* --observable NAME]+" << std::endl;
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
