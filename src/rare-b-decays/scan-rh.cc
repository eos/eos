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
};

struct Result
{
    double c9prime;

    double c10prime;

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
            // [BaBar2008] data
            Bin{10.24, 19.21, -0.44, -0.76, -1.28, "A_FB"},

            // [Belle2009] data
            Bin{14.18, 16.00, -0.96, -0.70, -0.38, "A_FB"},
            Bin{16.00, 19.21, -0.81, -0.66, -0.46, "A_FB"},
            Bin{14.18, 16.00, 0.71e-7, 1.05e-7, 1.42e-7, "BR"},
            Bin{16.00, 19.21, 1.64e-7, 2.04e-7, 2.47e-7, "BR"},

            // [CDF2010] data
            Bin{14.18, 16.00, -0.67, -0.42, -0.17, "A_FB"},
            Bin{16.00, 19.21, -0.96, -0.70, -0.35, "A_FB"},
            Bin{14.18, 16.00, 1.02e-7, 1.51e-7, 2.00e-7, "BR"},
            Bin{16.00, 19.21, 0.86e-7, 1.35e-7, 1.84e-7, "BR"},
        };

        Parameters parameters(Parameters::Defaults());
        Parameter c7prime = parameters.declare("c7prime");
        Parameter c9prime = parameters.declare("c9prime");
        Parameter c10prime = parameters.declare("c10prime");

        Kinematics kinematics;
        kinematics.declare("s_min");
        kinematics.declare("s_max");

        std::list<std::pair<Bin, ObservablePtr>> bins;
        for (auto d(data.begin()) ; d != data.end() ; ++d)
        {
            bins.push_back(std::make_pair(*d, BToKstarDileptonFactory::make(d->o_name, parameters)));
        }

        double max_likelihood(-1e7);
        std::list<Result> results;
        for (int i(-50) ; i <= 50 ; ++i)
        {
            c9prime = 1.0 * i / 3.0;

            for (int j(-50) ; j <= 50 ; ++j)
            {
                c10prime = 1 * j / 3.0;

                double likelihood(0.0);

                for (auto bin = bins.begin() ; bin != bins.end() ; ++bin)
                {
                    kinematics.set("s_min", bin->first.min);
                    kinematics.set("s_max", bin->first.max);

                    double value = bin->second->evaluate(kinematics) / (bin->first.max - bin->first.min);
                    double chi = (value - bin->first.o) / (bin->first.o_max - bin->first.o_min);
                    likelihood += -2.0 * std::log(std::abs(chi));
                }

                max_likelihood = std::max(max_likelihood, likelihood);
                results.push_back(Result{c9prime(), c10prime(), likelihood, false});
            }

            results.back().end = true;
        }

        std::cout << "# max_likelihood = " << max_likelihood << std::endl;
        std::map<signed, unsigned> distribution;
        for (auto r(results.begin()), r_end(results.end()) ; r != r_end ; ++r)
        {
            double likelihood = std::ceil((r->likelihood - max_likelihood) / bins.size());
            std::cout << r->c9prime << '\t' << r->c10prime << '\t' << likelihood << std::endl;
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
    catch(...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
