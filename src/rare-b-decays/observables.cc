/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>

#include <cstdlib>
#include <iostream>
#include <vector>

using namespace wf;

int
main(int argc, char * argv[])
{
    try
    {
        Parameters parameters(Parameters::Defaults());
        Kinematics kinematics;
        kinematics.declare("s");
        ObservableOptions options;
        std::vector<ObservablePtr> observables
        {
            BToKstarDileptonFactory::make("dBR/ds", parameters, options),
            BToKstarDileptonFactory::make("A_FB(s)", parameters, options),
            BToKstarDileptonFactory::make("F_L(s)", parameters, options)
        };

        const unsigned points = 100;

        std::cout << "#s(GeV^2) dGamma A_FB F_L" << std::endl;
        for (unsigned j = 0 ; j <= points ; ++j)
        {
            const double s_low = 13.932;
            const double s_high = 19.211;
            double s = s_low + j * (s_high - s_low) / points;

            std::cout << s << std::flush;
            kinematics.set("s", s);

            for (auto k(observables.begin()), k_end(observables.end()) ; k != k_end ; ++k)
            {
                if (*k)
                    std::cout << '\t' << (*k)->evaluate(kinematics);
            }
            std::cout << std::endl;
        }
    }
    catch (Exception & e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    catch (...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
