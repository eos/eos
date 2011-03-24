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
#include <src/utils/scan_file.hh>

#include <list>
#include <iomanip>
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

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        std::list<std::string> filenames;

        CommandLine()
        {
        }

        void parse(int argc, char ** argv)
        {
            for (char ** a(argv + 1), ** a_end(argv + argc) ; a != a_end ; ++a)
            {
                std::string argument(*a);

                if ("--file" == argument)
                {
                    filenames.push_back(std::string(*(++a)));

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

        if (CommandLine::instance()->filenames.empty())
            throw DoUsage("Need to specify at least one filename to inspect!");

        for (auto f = CommandLine::instance()->filenames.cbegin(), f_end = CommandLine::instance()->filenames.cend() ; f != f_end ; ++f)
        {
            std::cout << "# File: " << *f << std::endl;
            try
            {
                ScanFile file = ScanFile::Open(*f);

                std::cout << "#   Creator:     " << file.creator() << std::endl;
                std::cout << "#   EOS Version: " << file.eos_version() << std::endl;

                for (auto d = file.begin(), d_end = file.end() ; d != d_end ; ++d)
                {
                    std::cout << "#   Dataset '" + d->name() + "': (" << d->records() << ", " << d->fields() << ')' << std::endl;
                    for (auto f = d->begin_fields(), f_end = d->end_fields() ; f != f_end ; ++f)
                    {
                        std::cout << "#     Field '" << *f << "'" << std::endl;
                    }

                    std::cout << std::scientific << std::setprecision(9);
                    ScanFile::Record record = (*d)[0];
                    for (unsigned i = 0 ; i < d->records() ; ++i, ++record)
                    {
                        for (unsigned j = 0 ; j < d->fields() ; ++j)
                        {
                            std::cout << record[j] << '\t';
                        }

                        std::cout << std::endl;
                    }
                }
            }
            catch (ScanFileError & e)
            {
                std::cout << "#   Error reading " << *f << std::endl;
            }
        }
    }
    catch(DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-inspect" << std::endl;
        std::cout << "  [--file NAME]+" << std::endl;
        std::cout << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  eos-inspect --file input1.hdf5 --file input2.hdf5" << std::endl;
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
