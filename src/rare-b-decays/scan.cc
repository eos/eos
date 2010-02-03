/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>
#include <tr1/functional>

using namespace wf;

struct BinData { double min, max, o_min, o, o_max; std::string name; };

struct Bin { BinData data; ObservablePtr o; };

int
main(int argc, char * argv[])
{
    try
    {
        // TODO: Upper bin q^2 is 19.30 instead of 19.20. Change after more recent
        // PDG data has been included.
        std::list<BinData> bin_data = {
            // [BaBar2008] data
            //BinData{10.24, 19.20, -0.44, -0.76, -1.28, "A_FB"},

            // [Belle2009] data
            BinData{14.18, 16.00, -0.96, -0.70, -0.38, "A_FB"},
            BinData{16.00, 19.20, -0.81, -0.66, -0.46, "A_FB"},
            BinData{14.18, 16.00, 0.71e-7, 1.05e-7, 1.42e-7, "BR"},
            BinData{16.00, 19.20, 1.64e-7, 2.04e-7, 2.47e-7, "BR"},

            // [CDF2010] data
            BinData{14.18, 16.00, -0.67, -0.42, -0.17, "A_FB"},
            BinData{16.00, 19.20, -0.96, -0.70, -0.35, "A_FB"},
            BinData{14.18, 16.00, 1.02e-7, 1.51e-7, 2.00e-7, "BR"},
            BinData{16.00, 19.20, 0.86e-7, 1.35e-7, 1.84e-7, "BR"},
        };

        Parameters parameters(Parameters::StandardModell());
        Parameter c9 = parameters["c9"];
        Parameter c10 = parameters["c10"];

        Kinematics kinematics;
        kinematics.declare("s_min");
        kinematics.declare("s_max");

        std::list<Bin> bins;
        for (auto bd(bin_data.begin()) ; bd != bin_data.end() ; ++bd)
        {
            bins.push_back(Bin{*bd, BToKstarDileptonFactory::make(bd->name, parameters)});
        }

        for (int i(-50) ; i <= 50 ; ++i)
        {
            c9 = 1.0 * i / 3.0;

            for (int j(-50) ; j <= 50 ; ++j)
            {
                c10 = 1 * j / 3.0;

                double likelyhood(0.0);

                for (auto bin = bins.begin() ; bin != bins.end() ; ++bin)
                {
                    kinematics.set("s_min", bin->data.min);
                    kinematics.set("s_max", bin->data.max);

                    double value = bin->o->evaluate(kinematics) / (bin->data.max - bin->data.min);
                    double chi = (value - bin->data.o) / (bin->data.o_max - bin->data.o_min);
                    likelyhood += -2.0 * std::log(std::abs(chi));
                }

                std::cout << c9 << '\t' << c10 << '\t' << likelyhood << std::endl;
            }

            std::cout << std::endl;
        }
    }
    catch(...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
