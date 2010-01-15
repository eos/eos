/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/low-recoil/low-recoil.hh>

#include <cstdlib>
#include <iostream>

using namespace wf;

int
main(int argc, char * argv[])
{
    try
    {
        Decay<BToKstarDilepton> decay(4.0);

        std::cout << "#s(GeV^2) dGamma" << std::endl;

        const unsigned N = 100;
        for (unsigned i = 0 ; i < N ; ++i)
        {
            const double s_low = 17;
            const double s_high = 19.25;
            double s = s_low + i * (s_high - s_low) / N;

            double dGamma = decay.a_long(left_handed, s).absolute_squared()
                + decay.a_long(right_handed, s).absolute_squared()
                + decay.a_perp(left_handed, s).absolute_squared()
                + decay.a_perp(right_handed, s).absolute_squared()
                + decay.a_par(left_handed, s).absolute_squared()
                + decay.a_par(right_handed, s).absolute_squared();

            std::cout
                << s << "\t"
                << dGamma
                << std::endl;
        }
    }
    catch(...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
