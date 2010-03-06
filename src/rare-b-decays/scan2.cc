/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <src/utils/lock.hh>
#include <src/utils/mutex.hh>
#include <src/utils/thread_pool.hh>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <tr1/functional>
#include <utility>
#include <vector>

using namespace wf;

struct Input
{
    const double min;

    const double max;

    const double o_min;

    const double o;

    const double o_max;

    const std::string o_name;

    const std::string o_options;
};

class Scan2
{
    public:
        Mutex mutex;

        std::list<std::pair<Input, ObservablePtr>> bins;

        double max_likelihood;

        std::map<std::pair<double, double>, double> results;

        std::string x_name;

        std::string y_name;

        Scan2(const std::string & x_name, const std::string & y_name, const std::initializer_list<Input> & inputs) :
            max_likelihood(0.0),
            x_name(x_name),
            y_name(y_name)
        {
            for (auto i(inputs.begin()), i_end(inputs.end()) ; i != i_end ; ++i)
            {
                //TODO: Create options from i->o_opions!
                ObservableOptions options;
                bins.push_back(std::make_pair(*i, BToKstarDileptonFactory::make(i->o_name, Parameters::Defaults(), options)));
            }
        }

        void calc_likelihood(const double & x, const double & y)
        {
            double chi_squared(0.0);
            Kinematics k;
            k.declare("s_min");
            k.declare("s_max");

            for (auto bin = bins.begin() ; bin != bins.end() ; ++bin)
            {
                ObservablePtr o(bin->second->clone());

                k.set("s_min", bin->first.min);
                k.set("s_max", bin->first.max);
                o->parameters().set(x_name, x);
                o->parameters().set(y_name, y);

                double value = o->evaluate(k) / (bin->first.max - bin->first.min);
                double chi = (value - bin->first.o) / (bin->first.o_max - bin->first.o_min);
                chi_squared += chi * chi;
            }

            double likelihood = std::exp(-0.5 * chi_squared);

            {
                Lock l(mutex);
                results[std::make_pair(x, y)] = likelihood;
                max_likelihood = std::max(max_likelihood, likelihood);
            }
        }

        void scan()
        {
            TicketList tickets;

            for (int i(-50) ; i <= 50 ; ++i)
            {
                double x = 1.0 * i / 5.0;

                for (int j(-50) ; j <= 50 ; ++j)
                {
                    double y = 1.0 * j / 5.0;

                    tickets.push_back(ThreadPool::instance()->enqueue(std::tr1::bind(std::tr1::mem_fn(&Scan2::calc_likelihood), this, x, y)));
                }
            }

            tickets.wait();

            std::vector<double> likelihoods;
            likelihoods.reserve(results.size());
            double integral = 0.0;

            double index(results.begin()->first.second);
            for (auto r(results.begin()), r_end(results.end()) ; r != r_end ; ++r)
            {
                if (r->first.second < index)
                    std::cout << std::endl;

                index = r->first.second;
                std::cout << r->first.first << '\t' << r->first.second << '\t' << r->second / max_likelihood << std::endl;
                likelihoods.push_back(r->second / max_likelihood);
                integral += r->second / max_likelihood;
            }

            std::vector<double> upper(4, 1.0), lower(4, 0.0), value(4, 0.5);
            std::vector<double> ratio = { 0.683, 0.900, 0.950, 0.954 };
            for (unsigned i(0) ; i < 10 ; ++i)
            {
                std::vector<double> partial(4, 0.0);

                for (auto l(likelihoods.begin()), l_end(likelihoods.end()) ; l != l_end ; ++l)
                {
                    for (unsigned j(0) ; j < partial.size() ; ++j)
                        if (*l > value[j])
                            partial[j] += *l;
                }

                for (unsigned j(0) ; j < partial.size() ; ++j)
                {
                    if ((partial[j] / integral) > ratio[j])
                        lower[j] = value[j];
                    else
                        upper[j] = value[j];

                    value[j] = (upper[j] + lower[j]) / 2.0;
                }
            }

            std::cout << "# Confidence Levels" << std::endl;
            for (unsigned j(0) ; j < ratio.size() ; ++j)
            {
                std::cout << "# " << ratio[j] << " -> " << value[j] << std::endl;
            }

            std::cout << "# max(likelihood) = " << max_likelihood << std::endl;
        }
};

class DoUsage
{
    private:
        std::string _what;

    public:
        DoUsage(const std::string & what) :
            _what(what)
        {
        }

        const std::string & what() const
        {
            return _what;
        }
};

int
main(int argc, char * argv[])
{
    try
    {
        if (argc < 3)
            throw DoUsage("Need two parameter names");

        std::string x(argv[1]), y(argv[2]);

        // max(s) = (m_B - m_Kstar)^2 = 19.211
        std::initializer_list<Input> input = {
            // [BaBar2006] data
            //Input{00.10, 08.41,  0.12e-6,  0.27e-6,  0.44e-6, "BR@LargeRecoil", ""},
            //Input{10.24, 19.21,  0.21e-6,  0.37e-6,  0.55e-6, "BR@LowRecoil", ""},

            // [BaBar2008] data
            //Input{00.10, 06.25, +0.01,    -0.24,    -0.42,    "A_FB@LargeRecoil", ""},
            //Input{10.24, 19.21, -0.44,    -0.76,    -1.28,    "A_FB@LowRecoil", ""},

            // [Belle2009] data
            Input{00.10, 02.00, -0.12,    -0.47,    -0.76,    "A_FB@LargeRecoil", ""},
            Input{02.00, 04.30, +0.30,    -0.11,    -0.47,    "A_FB@LargeRecoil", ""},
            Input{04.30, 08.68, -0.11,    -0.45,    -0.73,    "A_FB@LargeRecoil", ""},
            Input{14.18, 16.00, -0.96,    -0.70,    -0.38,    "A_FB@LowRecoil", ""},
            Input{16.00, 19.21, -0.81,    -0.66,    -0.46,    "A_FB@LowRecoil", ""},
            Input{00.10, 02.00,  0.99e-7,  1.46e-7,  1.98e-7, "BR@LargeRecoil", ""},
            Input{02.00, 04.30,  0.52e-7,  0.86e-7,  1.24e-7, "BR@LargeRecoil", ""},
            Input{04.30, 08.68,  0.83e-7,  1.37e-7,  1.96e-7, "BR@LargeRecoil", ""},
            Input{14.18, 16.00,  0.71e-7,  1.05e-7,  1.42e-7, "BR@LowRecoil", ""},
            Input{16.00, 19.21,  1.64e-7,  2.04e-7,  2.47e-7, "BR@LowRecoil", ""},

            // [CDF2010] data
            Input{00.10, 02.00, +0.87,    -0.13,    -2.03,    "A_FB@LargeRecoil", ""},
            Input{02.00, 04.30, +0.36,    -0.19,    -0.73,    "A_FB@LargeRecoil", ""},
            Input{04.30, 08.68, +0.39,    +0.06,    -0.29,    "A_FB@LargeRecoil", ""},
            Input{14.18, 16.00, -0.67,    -0.42,    -0.17,    "A_FB@LowRecoil", ""},
            Input{16.00, 19.21, -0.96,    -0.70,    -0.35,    "A_FB@LowRecoil", ""},
            Input{00.10, 02.00,  0.49e-7,  0.98e-7,  1.47e-7, "BR@LargeRecoil", ""},
            Input{02.00, 04.30,  0.53e-7,  1.00e-7,  1.47e-7, "BR@LargeRecoil", ""},
            Input{04.30, 08.68,  0.97e-7,  1.69e-7,  1.41e-7, "BR@LargeRecoil", ""},
            Input{14.18, 16.00,  1.02e-7,  1.51e-7,  2.00e-7, "BR@LowRecoil", ""},
            Input{16.00, 19.21,  0.86e-7,  1.35e-7,  1.84e-7, "BR@LowRecoil", ""},
        };

        Scan2 scanner(x, y, input);
        scanner.scan();
    }
    catch(DoUsage & e)
    {
        std::cout << e.what() << std::endl;
    }
    catch(Exception & e)
    {
        std::cout << "Caught exception: '" << e.what() << "'" << std::endl;
        return EXIT_FAILURE;
    }
    catch(...)
    {
        std::cerr << "Aborting after unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
