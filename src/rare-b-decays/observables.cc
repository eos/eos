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
        std::vector<ObservablePtr> observables
        {
            BToKstarDileptonFactory::make("dBR/ds", parameters),
            BToKstarDileptonFactory::make("dA_FB/ds", parameters),
            BToKstarDileptonFactory::make("dF_L/ds", parameters)
        };

        const unsigned points = 100;
        std::vector<double> scales{ 2.0, 4.0, 8.0 };

        for (auto i = scales.begin() ; i != scales.end() ; ++i)
        {
            std::cout << "# mu = " << *i << std::endl;
            parameters.set("mu", *i);
            BToKstarDilepton<LowRecoil> decay(parameters);

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
            std::cout << std::endl << std::endl;
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
