/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <src/utils/destringify.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/histogram.hh>
#include <src/utils/scan_file.hh>
#include <src/utils/stringify.hh>

#include <cmath>
#include <list>
#include <iostream>
#include <string>

using namespace eos;

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

/* Marginalisation Classes */
double marginalise_by_sum(const double & previous, const double & next)
{
    return previous + next;
}

double marginalise_by_max_likelihood(const double & previous, const double & next)
{
    return std::max(previous, std::exp(-0.5 * next));
}

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        std::list<std::string> files;

        unsigned x_index, y_index;

        std::array<double, 2> start;

        std::array<double, 2> end;

        std::array<unsigned, 2> count;

        std::function<double (const double &, const double &)> marginalise;

        CommandLine() :
            x_index(0),
            y_index(0),
            start(std::array<double, 2>{{ 0.0, 0.0 }}),
            end(std::array<double, 2>{{ 15.0, 15.0 }}),
            count(std::array<unsigned, 2>{{ 60, 60 }}),
            marginalise(marginalise_by_sum)
        {
        }

        void parse(int argc, char ** argv)
        {
            for (char ** a(argv + 1), ** a_end(argv + argc) ; a != a_end ; ++a)
            {
                std::string argument(*a);

                if ("--max-exp" == argument)
                {
                    marginalise = &marginalise_by_max_likelihood;

                    continue;
                }

                if ("--file" == argument)
                {
                    files.push_back(std::string(*(++a)));

                    continue;
                }

                if ("--x" == argument)
                {
                    x_index = destringify<unsigned>(*(++a));
                    start[0] = destringify<double>(*(++a));
                    end[0] = destringify<double>(*(++a));
                    count[0] = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--y" == argument)
                {
                    y_index = destringify<unsigned>(*(++a));
                    start[1] = destringify<double>(*(++a));
                    end[1] = destringify<double>(*(++a));
                    count[1] = destringify<unsigned>(*(++a));

                    continue;
                }

                throw DoUsage("Unknown command line argument: " + argument);
            }
        }
};

int
main(int argc, char * argv[])
{
    try
    {
        CommandLine::instance()->parse(argc, argv);

        if ((0 == CommandLine::instance()->x_index) || (0 == CommandLine::instance()->y_index))
            throw DoUsage("Need to specify --x-index and --y-index with non-zero values!");

        unsigned x_index = CommandLine::instance()->x_index - 1;
        unsigned y_index = CommandLine::instance()->y_index - 1;

        Histogram<2> histogram = Histogram<2>::WithEqualBinning(CommandLine::instance()->start,
                CommandLine::instance()->end, CommandLine::instance()->count);

        for (auto f = CommandLine::instance()->files.cbegin(), f_end = CommandLine::instance()->files.cend() ; f != f_end ; ++f)
        {
            try
            {
                std::cout << "# " << *f << std::endl;
                ScanFile file = ScanFile::Open(*f);

                std::cout << "#   Creator:     " << file.creator() << std::endl;
                std::cout << "#   EOS Version: " << file.eos_version() << std::endl;

                for (auto d = file.begin(), d_end = file.end() ; d != d_end ; ++d)
                {
                    unsigned posterior_index = d->tuple_size() - 1;

                    if (posterior_index <= x_index)
                        throw DoUsage("X index '" + stringify(CommandLine::instance()->x_index) + "' exceeds number of parameters!");

                    if (posterior_index <= y_index)
                        throw DoUsage("Y index '" + stringify(CommandLine::instance()->y_index) + "' exceeds number of parameters!");

                    ScanFile::Tuple tuple = (*d)[0];

                    std::cout << "#   Data set '" << d->name() << "' with " << d->tuples() << " of " << d->tuple_size() << " elements each" << std::endl;
                    for (unsigned i = 0 ; i < d->tuples() ; ++i, ++tuple)
                    {
                        std::array<double, 2> coords{{ tuple[x_index], tuple[y_index] }};
                        auto b = histogram.find(coords);
                        if (histogram.end() == b)
                        {
                            std::cerr << "Did not find bin suitable for '(" << coords[0] << ", " << coords[1] << ")'. You might need to adjust the histogram configuration!" << std::endl;
                            continue;
                        }

                        b->value = CommandLine::instance()->marginalise(b->value, tuple[posterior_index]);
                    }
                }
            }
            catch (ScanFileError & e)
            {
                std::cout << "#   Error reading " << *f << std::endl;
            }

            auto b = histogram.cbegin(), b_end = histogram.cend();
            double last_y = b->lower[1];
            for ( ; b != b_end ; ++b)
            {
                if (b->lower[1] < last_y)
                    std::cout << std::endl;

                last_y = b->lower[1];

                std::cout << b->lower[0] << '\t' << b->lower[1] << '\t' << b->upper[0] << '\t' << b->upper[1] << '\t' << b->value << std::endl;
            }
        }
    }
    catch(DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-marginalise" << std::endl;
        std::cout << "  --x INDEX MIN MAX COUNT" << std::endl;
        std::cout << "  --y INDEX MIN MAX COUNT" << std::endl;
        std::cout << "  [--file NAME]+" << std::endl;
        std::cout << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  eos-marginalise --x 1 0.0 15.0 60 --y 2 0.0 6.28 32 --file input1.hdf5 --file input2.hdf5" << std::endl;
    }
    catch(Exception & e)
    {
        std::cerr << "Caught exception: '" << e.what() << "'" << std::endl;
        return EXIT_FAILURE;
    }
    catch(...)
    {
        std::cerr << "Aborting after unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
