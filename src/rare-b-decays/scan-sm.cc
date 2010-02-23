/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>
#include <map>

using namespace wf;

class Bin
{
    public:
        const double min;

        const double max;

        const double o_min;

        const double o;

        const double o_max;

        const std::string o_name;

        const std::string o_options;
};

struct Result
{
    double c9;

    double c10;

    double likelihood;

    bool end;
};

int
main(int argc, char * argv[])
{
    try
    {
        // max(s) = (m_B - m_Kstar)^2 = 19.211
        std::list<Bin> data = {
            // [BaBar2006] data
            Bin{10.24, 19.21, -0.38, -0.72, -1.08, "A_FB", ""},
            Bin{10.24, 19.21, 0.21e-7, 0.37e-6, 0.55e-6, "BR", ""},

            // [Belle2009] data
            Bin{14.18, 16.00, -0.96, -0.70, -0.38, "A_FB", ""},
            Bin{16.00, 19.21, -0.81, -0.66, -0.46, "A_FB", ""},
            Bin{14.18, 16.00, 0.71e-7, 1.05e-7, 1.42e-7, "BR", ""},
            Bin{16.00, 19.21, 1.64e-7, 2.04e-7, 2.47e-7, "BR", ""},

            // [CDF2010] data
            Bin{14.18, 16.00, -0.67, -0.42, -0.17, "A_FB", ""},
            Bin{16.00, 19.21, -0.96, -0.70, -0.35, "A_FB", ""},
            Bin{14.18, 16.00, 1.02e-7, 1.51e-7, 2.00e-7, "BR", ""},
            Bin{16.00, 19.21, 0.86e-7, 1.35e-7, 1.84e-7, "BR", ""},
        };

        Parameters parameters(Parameters::Defaults());
        Parameter c9 = parameters["c9"];
        Parameter c10 = parameters["c10"];

        Kinematics kinematics;
        kinematics.declare("s_min");
        kinematics.declare("s_max");

        std::list<std::pair<Bin, ObservablePtr>> bins;
        for (auto d(data.begin()) ; d != data.end() ; ++d)
        {
            // TODO: Read options from o_options
            ObservableOptions options;
            bins.push_back(std::make_pair(*d, BToKstarDileptonFactory::make(d->o_name, parameters, options)));
        }

        double max_likelihood(-1e7);
        std::list<Result> results;
        for (int i(-50) ; i <= 50 ; ++i)
        {
            c9 = 1.0 * i / 3.0;

            for (int j(-50) ; j <= 50 ; ++j)
            {
                c10 = 1 * j / 3.0;

                double chi_squared(0.0);

                for (auto bin = bins.begin() ; bin != bins.end() ; ++bin)
                {
                    kinematics.set("s_min", bin->first.min);
                    kinematics.set("s_max", bin->first.max);

                    double value = bin->second->evaluate(kinematics) / (bin->first.max - bin->first.min);
                    double chi = (value - bin->first.o) / (bin->first.o_max - bin->first.o_min);
                    chi_squared += chi * chi;
                }
                double likelihood = std::exp(-0.5 * chi_squared);
                max_likelihood = std::max(max_likelihood, likelihood);
                results.push_back(Result{c9(), c10(), likelihood, false});
            }

            results.back().end = true;
        }

        std::cout << "# max_likelihood = " << max_likelihood << std::endl;
        std::map<double, unsigned> distribution;
        for (auto r(results.begin()), r_end(results.end()) ; r != r_end ; ++r)
        {
            double likelihood = 0.1 * std::floor(10 * r->likelihood / max_likelihood);
            std::cout << r->c9 << '\t' << r->c10 << '\t' << likelihood << std::endl;
            distribution[likelihood] += 1;

            if (r->end)
                std::cout << std::endl;
        }

        std::cout << std::endl;
        std::cout << "# Distribution" << std::endl;
        for (auto d(distribution.begin()) ; d != distribution.end() ; ++d)
        {
            std::cout << "# " << d->first << " : " << d->second << std::endl;
        }
    }
    catch(std::string & e)
    {
        std::cout << "Caught exception: '" << e << "'" << std::endl;
        return EXIT_FAILURE;
    }
    catch(...)
    {
        std::cerr << "Aborting after unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
