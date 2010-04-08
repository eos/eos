/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/factory.hh>
#include <src/utils/destringify.hh>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>

using namespace wf;

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
        bool range(false);
        std::list<Parameter> variations;
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
                range = true;
                std::cerr << "Range: " << s_min << " .. " << s_max << std::endl;
                continue;
            }

            if ("--vary" == argument)
            {
                std::string name(*(++a));
                variations.push_back(parameters[name]);
                continue;
            }

            if ("--observable" != argument)
                throw DoUsage("Unknown option: '" + argument + "'");

            ObservableOptions options;
            argument = *(++a);

            std::string::size_type pos;
            while (std::string::npos != (pos = argument.rfind(',')))
            {
                std::string::size_type sep(argument.find('=', pos + 1));
                if (std::string::npos == sep)
                    throw DoUsage("Invalid observable option: '" + argument.substr(pos + 1) + "'");

                std::string key(argument.substr(pos + 1, sep - pos - 1));
                std::string value(argument.substr(sep + 1));

                options.set(key, value);
                argument.erase(pos);
            }

            ObservablePtr ptr(RareBFactory::make(argument, parameters, options));
            if (! ptr)
                throw DoUsage("Unknown observable: '" + argument + "'");

            observables.push_back(ptr);
        }

        if (observables.empty())
            throw DoUsage("Need at least one observable");

        const unsigned points = 50;

        kinematics.set("s_min", s_min);
        kinematics.set("s_max", s_max);

        for (auto o(observables.begin()), o_end(observables.end()) ; o != o_end ; ++o)
        {
            double central = (*o)->evaluate(kinematics);
            double delta_min = 0.0, delta_max = 0.0;

            for (auto p(variations.begin()), p_end(variations.end()) ; p != p_end ; ++p)
            {
                double old_p = *p;
                double max = 0.0, min = 0.0, value;

                *p = p->min();
                value = (*o)->evaluate(kinematics);
                if (value > central)
                    max = value - central;

                if (value < central)
                    min = central - value;

                *p = p->max();
                value = (*o)->evaluate(kinematics);
                if (value > central)
                    max = std::max(max, value - central);

                if (value < central)
                    min = std::max(min, central - value);

                *p = old_p;

                delta_min += min * min;
                delta_max += max * max;
            }

            delta_max = std::sqrt(delta_max);
            delta_min = std::sqrt(delta_min);

            std::cout << (*o)->name() << " = " << central << " -" << delta_min << " +"  << delta_max << std::endl;
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
    catch (...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
