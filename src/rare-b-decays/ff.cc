/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <src/rare-b-decays/form-factors.hh>

#include <cstdlib>
#include <iostream>

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
        if (2 != argc)
            throw DoUsage("Need at exactly one set of form factors");

        std::shared_ptr<FormFactors<BToKstar>> set(FormFactorFactory<BToKstar>::create(argv[1], Parameters::Defaults()));
        if (! set)
            throw DoUsage("Unknown set of form factors: '" + std::string(argv[1]) + "'");

        const unsigned points = 300;

        std::cout << "#s\tV\tA0\tA1\tA2\txi_perp\txi_par" << std::endl;

        for (unsigned j = 0 ; j <= points ; ++j)
        {
            const double s_low = 0;
            const double s_high = 19.211;
            double s = s_low + j * (s_high - s_low) / points;
            double s_hat = s / 5.28 / 5.28;

            std::cout
                << s << '\t'
                << set->v(s_hat) << '\t'
                << set->a_0(s_hat) << '\t'
                << set->a_1(s_hat) << '\t'
                << set->a_2(s_hat) << '\t'
                << set->v(s_hat) * 5.28 / (5.28 + 0.896) << '\t'
                << (5.28 + 0.896) * 5.28 / (5.28*5.28 + 0.896*0.896 - s) * set->a_1(s_hat) - (1.0 - 0.896/5.28) * set->a_2(s_hat) << '\t'
                << std::endl;
        }
    }
    catch (DoUsage & e)
    {
        std::cout << e.what << std::endl;
        std::cout << "Usage: ff FORMFACTORSET" << std::endl;
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
