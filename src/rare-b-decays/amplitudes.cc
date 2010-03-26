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
        ObservableOptions options;
        options.set("form-factors", "ABHH1999::2");
        BToKstarDilepton<LargeRecoil> decay(Parameters::Defaults(), options);

        std::cout << "#s(GeV^2) A_0^L A_0^R A_pp^L A_pp^R A_pa^L A_pa^R" << std::endl;

        const unsigned N = 100;
        const static double Gamma = 6.58211899e-22 * 1e-3 / 1.53e-12; // cf. [PDG2008], hbar / tau_B
        for (unsigned i = 0 ; i <= N ; ++i)
        {
            const double s_low = 1.0;
            const double s_high = 6.0;
            double s = s_low + i * (s_high - s_low) / N;

            std::cout
                << s << "\t"
                << decay.a_long(left_handed, s).real() << "\t"
                << decay.a_long(left_handed, s).imaginary() << "\t"
                << decay.a_long(right_handed,s).real() << "\t"
                << decay.a_long(right_handed,s).imaginary() << "\t"
                << decay.a_perp(left_handed, s).real() << "\t"
                << decay.a_perp(left_handed, s).imaginary() << "\t"
                << decay.a_perp(right_handed,s).real() << "\t"
                << decay.a_perp(right_handed,s).imaginary() << "\t"
                << decay.a_par(left_handed,  s).real() << "\t"
                << decay.a_par(left_handed,  s).imaginary() << "\t"
                << decay.a_par(right_handed, s).real() << "\t"
                << decay.a_par(right_handed, s).imaginary() << "\t"
                << std::endl;
        }
    }
    catch(...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
