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
#include <limits>
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

        std::string x_name, y_name;

        std::array<double, 2> start;

        std::array<double, 2> end;

        std::array<unsigned, 2> count;

        std::function<double (const double &, const double &)> marginalise;

        CommandLine() :
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
                    x_name = *(++a);
                    start[0] = destringify<double>(*(++a));
                    end[0] = destringify<double>(*(++a));
                    count[0] = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--y" == argument)
                {
                    y_name = *(++a);
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

        if (CommandLine::instance()->x_name.empty() || CommandLine::instance()->y_name.empty())
            throw DoUsage("Need to specify no-empty field names for --x and --y!");

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
                    unsigned x_index = d->find_field_index(CommandLine::instance()->x_name),
                             y_index = d->find_field_index(CommandLine::instance()->y_name),
                             posterior_index = d->find_field_index("posterior");

                    ScanFile::Record record = (*d)[0];

                    std::cout << "#   Data set '" << d->name() << "' with " << d->records() << " records of " << d->fields() << " elements each" << std::endl;
                    for (unsigned i = 0 ; i < d->records() ; ++i, ++record)
                    {
                        std::array<double, 2> coords{{ record[x_index], record[y_index] }};
                        auto b = histogram.find(coords);
                        if (histogram.end() == b)
                        {
                            std::cerr << "Did not find bin suitable for '(" << coords[0] << ", " << coords[1] << ")'. You might need to adjust the histogram configuration!" << std::endl;
                            continue;
                        }

                        b->value = CommandLine::instance()->marginalise(b->value, record[posterior_index]);
                    }
                }
            }
            catch (ScanFileError & e)
            {
                std::cout << "#   Error reading " << *f << std::endl;
            }

            auto b = histogram.cbegin(), b_end = histogram.cend();
            double last_y = b->lower[1];
            double min_value = std::numeric_limits<double>::max(), max_value = -std::numeric_limits<double>::max();
            double integral = 0.0;
            for ( ; b != b_end ; ++b)
            {
                if (b->lower[1] < last_y)
                    std::cout << std::endl;

                last_y = b->lower[1];

                std::cout << b->lower[0] << '\t' << b->lower[1] << '\t' << b->upper[0] << '\t' << b->upper[1] << '\t' << b->value << std::endl;
                min_value = std::min(min_value, b->value);
                max_value = std::max(max_value, b->value);
                integral += b->value;
            }

            std::array<double, 3> ratios{{ 0.683, 0.954, 0.997 }};
            std::array<double, 3> partials;
            std::array<double, 3> thresholds, upper_bounds, lower_bounds;
            thresholds.fill((min_value + max_value) / 2.0);
            upper_bounds.fill(max_value);
            lower_bounds.fill(min_value);

            for (auto i = 0 ; i < 10 ; ++i)
            {
                partials.fill(0.0);

                for (auto b = histogram.cbegin(), b_end = histogram.cend() ; b != b_end ; ++b)
                {
                    for (unsigned j = 0 ; j < thresholds.size() ; ++j)
                    {
                        if (b->value >= thresholds[j])
                            partials[j] += b->value;
                    }
                }

                for (unsigned j = 0 ; j < thresholds.size() ; ++j)
                {
                    double ratio = partials[j] / integral;

                    if (ratio > ratios[j])
                        lower_bounds[j] = thresholds[j];

                    if (ratio < ratios[j])
                        upper_bounds[j] = thresholds[j];

                    thresholds[j] = (upper_bounds[j] + lower_bounds[j]) / 2.0;
                }
            }

            for (unsigned j = 0 ; j < thresholds.size() ; ++j)
            {
                std::cout << "# " << ratios[j] << " -> " << thresholds[j] << " @ " << partials[j] / integral << std::endl;
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
