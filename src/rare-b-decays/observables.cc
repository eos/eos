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
        Decay<BToKstarDilepton> decay(4.0);

        std::cout << "#s(GeV^2) dGamma A_FB F_L" << std::endl;

        const unsigned N = 100;
        const static double Gamma = 6.58211899e-22 * 1e-3 / 1.52e-12; // cf. [PDG2006], hbar / tau_B
        for (unsigned i = 0 ; i <= N ; ++i)
        {
            const double s_low = 13.932;
            const double s_high = 19.211;
            double s = s_low + i * (s_high - s_low) / N;

            double dGamma = decay.a_long(left_handed, s).absolute_squared()
                + decay.a_long(right_handed, s).absolute_squared()
                + decay.a_perp(left_handed, s).absolute_squared()
                + decay.a_perp(right_handed, s).absolute_squared()
                + decay.a_par(left_handed, s).absolute_squared()
                + decay.a_par(right_handed, s).absolute_squared();

            double J_6 = 1.5 * (
                    (decay.a_par(left_handed, s) * decay.a_perp(left_handed, s).conjugate()).real()
                    -(decay.a_par(right_handed, s) * decay.a_perp(right_handed, s).conjugate()).real());

            double num_F_L = decay.a_long(left_handed, s).absolute_squared()
                + decay.a_long(right_handed, s).absolute_squared();

            std::cout
                << s << "\t"
                << dGamma / Gamma << "\t"
                << J_6 / dGamma << "\t"
                << num_F_L / dGamma
                << std::endl;
        }
    }
    catch(...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
