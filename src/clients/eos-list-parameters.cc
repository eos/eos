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

#include <src/factory.hh>
#include <src/utils/destringify.hh>
#include <src/utils/histogram.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/one-of.hh>
#include <src/utils/random_number_engine.hh>
#include <src/utils/stringify.hh>
#include <src/utils/wilson-polynomial.hh>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <list>
#include <random>
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
        Parameters parameters;

        CommandLine() :
            parameters(Parameters::Defaults())
        {
        }

        void parse(int argc, char ** argv)
        {
            std::shared_ptr<Kinematics> kinematics(new Kinematics);

            for (char ** a(argv + 1), ** a_end(argv + argc) ; a != a_end ; ++a)
            {
                std::string argument(*a);

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

        std::string::size_type max_name_length = 20;
        for (auto p = CommandLine::instance()->parameters.begin(), p_end = CommandLine::instance()->parameters.end() ; p != p_end ; ++p)
        {
            max_name_length = std::max(max_name_length, p->name().length());
        }

        for (auto p = CommandLine::instance()->parameters.begin(), p_end = CommandLine::instance()->parameters.end() ; p != p_end ; ++p)
        {
            std::cout
                << std::setw(max_name_length) << std::setiosflags(std::ios::right)
                << p->name() << '\t'
                << std::setw(7) << std::scientific << std::setprecision(4) << std::setiosflags(std::ios::left | std::ios::showpos)
                << p->min() << '\t' << (*p)() << '\t' << p->max()
                << std::endl;
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
