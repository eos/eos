/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>

#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace wf;

int
main(int argc, char * argv[])
{
    try
    {
        ObservableOptions options;
        options.set("form-factors", "BZ2004");
        BToKstarDilepton<LowRecoil> decay(Parameters::Defaults(), options);

        std::cout << "#s(GeV^2) A_0^L A_0^R A_pp^L A_pp^R A_pa^L A_pa^R" << std::endl;

        const unsigned N = 20;
        for (unsigned i = 0 ; i <= N ; ++i)
        {
            const double s_low = 14.0;
            const double s_high = 19.21;
            double s = s_low + i * (s_high - s_low) / N;

            std::cout
                << std::setprecision(5)
                << std::scientific
                << s << "\t"
                << real(decay.a_long(left_handed, s)) << "\t"
                << imag(decay.a_long(left_handed, s)) << "\t"
                << real(decay.a_long(right_handed,s)) << "\t"
                << imag(decay.a_long(right_handed,s)) << "\t"
                << real(decay.a_perp(left_handed, s)) << "\t"
                << imag(decay.a_perp(left_handed, s)) << "\t"
                << real(decay.a_perp(right_handed,s)) << "\t"
                << imag(decay.a_perp(right_handed,s)) << "\t"
                << real(decay.a_par(left_handed,  s)) << "\t"
                << imag(decay.a_par(left_handed,  s)) << "\t"
                << real(decay.a_par(right_handed, s)) << "\t"
                << imag(decay.a_par(right_handed, s)) << "\t"
                << std::endl;
        }
    }
    catch(...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
