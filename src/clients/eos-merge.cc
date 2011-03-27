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
#include <set>
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
        std::list<std::string> files;

        std::string output;

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
                    files.push_back(std::string(*(++a)));

                    continue;
                }

                if ("--output" == argument)
                {
                    output = std::string(*(++a));

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
            throw DoUsage("Need to specify at least one input file");

        std::set<std::string> set_names;
        std::vector<ScanFile> files;

        try
        {
            for (auto f = CommandLine::instance()->files.cbegin(), f_end = CommandLine::instance()->files.cend() ; f != f_end ; ++f)
            {
                std::cout << "# Merging these files" << std::endl;
                std::cout << "# " << *f << std::endl;
                ScanFile file = ScanFile::Open(*f);

                std::cout << "#   Creator:     " << file.creator() << std::endl;
                std::cout << "#   EOS Version: " << file.eos_version() << std::endl;

                files.push_back(file);
                for (auto d = file.begin(), d_end = file.end() ; d != d_end ; ++d)
                {
                    set_names.insert(d->name());
                }
            }

            std::vector<std::vector<ScanFile::DataSet>> all_data_sets(set_names.size());
            for (auto f = files.begin(), f_end = files.end() ; f != f_end ; ++f)
            {
                auto a = all_data_sets.begin();
                for (auto d = f->begin(), d_end = f->end() ; d != d_end ; ++d, ++a)
                {
                    a->push_back(*d);
                }
            }

            ScanFile output = ScanFile::Create(CommandLine::instance()->output, "eos-merge");

            for (auto a = all_data_sets.begin(), a_end = all_data_sets.end() ; a != a_end ; ++a)
            {
                std::vector<ScanFile::DataSet> data_sets(*a);

                std::string name = data_sets[0].name();
                unsigned fields = data_sets[0].fields();
                unsigned records = data_sets[0].records();
                std::vector<ScanFile::FieldInfo> field_infos(data_sets[0].begin_fields(), data_sets[0].end_fields());

                ScanFile::DataSet output_set = output.add(name, fields);
                std::copy(data_sets[0].begin_fields(), data_sets[0].end_fields(), output_set.begin_fields());
                ScanFile::WriteBuffer output_buffer(fields);

                std::vector<ScanFile::Record> all_records;
                for (auto d = data_sets.begin(), d_end = data_sets.end() ; d != d_end ; ++d)
                {
                    all_records.push_back((*d)[0]);
                }

                std::vector<double> values;
                for (unsigned i = 0 ; i < records ; ++i)
                {
                    values = all_records[0].data();

                    for (auto r = all_records.begin(), r_end = all_records.end() ; r != r_end ; ++r)
                    {
                        values[fields - 1] += (*r)[fields - 1];
                        ++(*r);
                    }

                    output_buffer << values;
                    if (output_buffer.capacity() == output_buffer.size())
                    {
                        output_set << output_buffer;
                        output_buffer.clear();
                    }
                }

                output_set << output_buffer;
            }
        }
        catch (ScanFileError & e)
        {
            std::cout << "# Encountered ScanFileError: '" << e.what() << "'" << std::endl;
        }
    }
    catch(DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-merge" << std::endl;
        std::cout << "  [--file NAME]+" << std::endl;
        std::cout << "  --output NAME" << std::endl;
        std::cout << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  eos-merge --file input1.hdf5 --file input2.hdf5 --output output.hdf5" << std::endl;
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

