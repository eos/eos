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

        for (auto f = CommandLine::instance()->filenames.cbegin(), f_end = CommandLine::instance()->filenames.cend() ; f != f_end ; ++f)
        {
            std::cout << "# " << *f << std::endl;
            try
            {
                ScanFile file = ScanFile::Open(*f);

                std::cout << "#   Creator:          " << file.creator() << std::endl;
                std::cout << "#   EOS Version:      " << file.eos_version() << std::endl;
                std::cout << "#   Number of tuples: " << file.scan_size() << std::endl;
                std::cout << "#   Tuple size:       " << file.tuple_size() << std::endl;

                ScanFile::Tuple tuple = file[0];
                for (auto i = 0 ; i < file.scan_size() ; ++i, ++tuple)
                {
                    tuple.read();

                    for (auto j = 0 ; j < file.tuple_size() ; ++j)
                    {
                        std::cout << tuple[j] << '\t';
                    }

                    std::cout << std::endl;
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
        std::cout << "Usage: eos-list-parameters" << std::endl;
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
