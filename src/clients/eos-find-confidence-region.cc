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
#include <src/utils/log.hh>
#include <src/utils/scan_file.hh>
#include <src/utils/stringify.hh>

#include <cmath>
#include <limits>
#include <iostream>
#include <string>
#include <vector>

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

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        std::vector<std::string> files;

        unsigned max_iterations;

        CommandLine() :
            max_iterations(10)
        {
        }

        void parse(int argc, char ** argv)
        {
            Log::instance()->set_program_name("eos-find-confidence-region");

            for (char ** a(argv + 1), ** a_end(argv + argc) ; a != a_end ; ++a)
            {
                std::string argument(*a);

                if ("--file" == argument)
                {
                    files.push_back(std::string(*(++a)));

                    continue;
                }

                if ("--max-iterations" == argument)
                {
                    max_iterations = destringify<unsigned>(*(++a));

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

        if (CommandLine::instance()->files.empty())
            throw DoUsage("Need to specify at least one input file!");

        if (10 > CommandLine::instance()->max_iterations)
            throw DoUsage("Need at least 10 iterations for any amount of accuracy!");

        double integral = 0.0;
        double max_value = -std::numeric_limits<double>::max(), min_value = std::numeric_limits<double>::max();
        std::vector<ScanFile> files;

        for (auto f = CommandLine::instance()->files.cbegin(), f_end = CommandLine::instance()->files.cend() ; f != f_end ; ++f)
        {
            try
            {
                std::cout << "# " << *f << std::endl;
                ScanFile file = ScanFile::Open(*f);
                files.push_back(file);

                std::cout << "#   Creator:     " << file.creator() << std::endl;
                std::cout << "#   EOS Version: " << file.eos_version() << std::endl;

                for (auto d = file.begin(), d_end = file.end() ; d != d_end ; ++d)
                {
                    std::cout << "#   Data set '" << d->name() << "' with " << d->records() << " records of " << d->fields() << " elements each" << std::endl;

                    unsigned posterior_index = d->find_field_index("posterior");
                    ScanFile::Record record = (*d)[0];

                    for (unsigned i = 0 ; i < d->records() ; ++i, ++record)
                    {
                        double value = std::exp(-1.0 * record[posterior_index]);
                        max_value = std::max(max_value, value);
                        min_value = std::min(min_value, value);
                        integral += value;
                    }
                }
            }
            catch (ScanFileError & e)
            {
                std::cout << "#   Error reading " << *f << std::endl;
            }
        }

        // Find 68.3% CR, 95.4% CR, 99.7% CR
        std::array<double, 3> ratios{{ 0.683, 0.954, 0.997 }};
        std::array<double, 3> partials;
        std::array<double, 3> thresholds, upper_bounds, lower_bounds;
        thresholds.fill((min_value + max_value) / 2.0);
        upper_bounds.fill(max_value);
        lower_bounds.fill(min_value);

        for (unsigned i = 0 ; i < CommandLine::instance()->max_iterations ; ++i)
        {
            std::cout << "# Iteration #" << (i + 1) << std::endl;

            partials.fill(0.0);

            for (auto f = files.begin(), f_end = files.end() ; f != f_end ; ++f)
            {
                try
                {
                    std::cout << "#   " << f->file_name() << std::endl;

                    for (auto d = f->begin(), d_end = f->end() ; d != d_end ; ++d)
                    {
                        std::cout << "#     " << d->name() << std::endl;

                        unsigned posterior_index = d->find_field_index("posterior");
                        ScanFile::Record record = (*d)[0];

                        for (unsigned j = 0 ; j < d->records() ; ++j, ++record)
                        {
                            double value = std::exp(-1.0 * record[posterior_index]);

                            for (unsigned k = 0 ; k < ratios.size() ; ++k)
                            {
                                if (value < thresholds[k])
                                    continue;

                                partials[k] += value;
                            }
                        }
                    }
                }
                catch (ScanFileError & e)
                {
                    std::cout << "#   Error reading " << f->file_name() << std::endl;
                }
            }

            bool converged = true;
            for (unsigned j = 0 ; j < ratios.size() ; ++j)
            {
                double ratio = partials[j] / integral;
                std::cout << "# Target: " << ratios[j] << ":: Threshold " << thresholds[j] << " yields " << ratio << std::endl;

                if (std::abs(ratio - ratios[j]) < 0.001)
                {
                    std::cout << "#   converged!" << std::endl;

                    continue;
                }
                else
                {
                    converged = false;
                }

                if (ratio > ratios[j])
                    lower_bounds[j] = thresholds[j];

                if (ratio < ratios[j])
                    upper_bounds[j] = thresholds[j];

                thresholds[j] = (upper_bounds[j] + lower_bounds[j]) / 2.0;
                std::cout << "#   New threshold: " << thresholds[j] << std::endl;
            }

            if (converged)
                break;
        }

        std::cout << "# Results:" << std::endl;
        for (unsigned i = 0 ; i < thresholds.size() ; ++i)
        {
            std::cout << "#  Target " << ratios[i] << ": Threshold " << thresholds[i] << " yields " << partials[i] / integral << std::endl;
        }
    }
    catch(DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-marginalise" << std::endl;
        std::cout << "  --max-iterations NUMBER" << std::endl;
        std::cout << "  [--file NAME]+" << std::endl;
        std::cout << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  eos-marginalise --max-iterations 10 --file input1.hdf5 --file input2.hdf5" << std::endl;
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
