/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>

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
        if (2 > argc)
            throw DoUsage("Need at least one observable");

        Parameters parameters(Parameters::Defaults());
        Kinematics kinematics;
        kinematics.declare("s");

        std::list<ObservablePtr> observables;

        for (char ** o(argv + 1), ** o_end(argv + argc) ; o != o_end ; ++o)
        {
            ObservableOptions options;
            std::string observable(*o);

            std::string::size_type pos;
            while (std::string::npos != (pos = observable.rfind(',')))
            {
                std::string::size_type sep(observable.find('=', pos + 1));
                if (std::string::npos == sep)
                    throw DoUsage("Invalid option: '" + observable.substr(pos + 1) + "'");

                std::string key(observable.substr(pos + 1, sep - pos - 1));
                std::string value(observable.substr(sep + 1));

                options.set(key, value);
                observable.erase(pos);
            }

            ObservablePtr ptr(BToKstarDileptonFactory::make(observable, parameters, options));
            if (! ptr)
                throw DoUsage("Unknown observable: '" + observable + "'");

            observables.push_back(ptr);
        }

        const unsigned points = 300;

        std::cout << "#\ts";
        for (auto o(observables.begin()), o_end(observables.end()) ; o != o_end ; ++o)
        {
            std::cout << '\t' << (*o)->name();
        }
        std::cout << std::endl;

        for (unsigned j = 0 ; j <= points ; ++j)
        {
            const double s_low = 0.0;
            const double s_high = 19.211;
            double s = s_low + j * (s_high - s_low) / points;

            std::cout << s << std::flush;
            kinematics.set("s", s);

            for (auto o(observables.begin()), o_end(observables.end()) ; o != o_end ; ++o)
            {
                std::cout << '\t' << (*o)->evaluate(kinematics);
            }

            std::cout << std::endl;
        }
    }
    catch (DoUsage & e)
    {
        std::cout << e.what << std::endl;
        std::cout << "Usage: observables OBSERVABLE [OBSERVABLE [...]]" << std::endl;
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
