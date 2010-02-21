/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>

#include <cstdlib>
#include <iostream>

using namespace wf;

int
main(int argc, char * argv[])
{
    try
    {
        BToKstarDilepton<LowRecoil> decay(Parameters::Defaults());

        std::cout << "#s(GeV^2) dGamma A_FB F_L" << std::endl;

        const unsigned N = 100;
        const static double Gamma = 6.58211899e-22 * 1e-3 / 1.53e-12; // cf. [PDG2008], hbar / tau_B
        for (unsigned i = 0 ; i <= N ; ++i)
        {
            const double s_low = 13.932;
            const double s_high = 19.211;
            double s = s_low + i * (s_high - s_low) / N;

            std::cout
                << s << "\t"
                << decay.a_long(left_handed, s).absolute_squared() << "\t"
                << decay.a_long(right_handed,s).absolute_squared() << "\t"
                << decay.a_perp(left_handed, s).absolute_squared() << "\t"
                << decay.a_perp(right_handed,s).absolute_squared() << "\t"
                << decay.a_par(left_handed,  s).absolute_squared() << "\t"
                << decay.a_par(right_handed, s).absolute_squared() << "\t"
                << std::endl;
        }
    }
    catch(...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
