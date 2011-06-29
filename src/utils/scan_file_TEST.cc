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

#include <config.h>
#include <test/test.hh>
#include <src/utils/scan_file.hh>

#include <cmath>
#include <cstdio>
#include <string>
#include <set>
#include <tuple>

#include <iostream>

using namespace test;
using namespace eos;

class ScanFileTest :
    public TestCase
{
    public:
        ScanFileTest() :
            TestCase("scan_file_test")
        {
        }

        virtual void run() const
        {
            static const std::string filename(EOS_BUILDDIR "/src/utils/scan_file_TEST.hdf5");

            std::remove(filename.c_str());

            // Create file
            {
                ScanFile test_file = ScanFile::Create(filename, "scan_file_TEST");

                // Set 'result #1'
                {
                    ScanFile::DataSet test_set = test_file.add("result #1", 3);
                    unsigned idx = 1;
                    for (auto f = test_set.begin_fields(), f_end = test_set.end_fields() ; f != f_end ; ++f, ++idx)
                    {
                        f->name("result #1, field #" + stringify(idx));
                        f->set("min", 0.0);
                        f->set("max", idx * 1.0);
                        f->set("nuisance", false);
                    }

                    test_set << std::vector<double>{ 3.0, 2.0, 1.0 };

                    ScanFile::WriteBuffer test_buffer(3);
                    for (auto i = 0 ; i < 1002 ; ++i)
                    {
                        test_buffer << std::vector<double>{ 4.0, 5.0, 6.0 };
                    }
                    test_set << test_buffer;

                    test_set << std::vector<double>{ 7.0, 8.0, 9.0 };
                }

                // Set 'result #2'
                {
                    ScanFile::DataSet test_set = test_file.add("result #2", 7);
                    test_set << std::vector<double>{ 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 };
                }

                // Set 'result #3'
                {
                    ScanFile::DataSet test_set = test_file.add("result #3", 5);
                    ScanFile::WriteBuffer test_buffer(5);

                    for (auto i = 0 ; i < 1023 ; ++i) // 1024 currently max capacity
                    {
                        test_buffer << std::vector<double>{ 4.0, 5.0, 6.0, 7.0, 8.0 };
                    }
                    test_buffer << std::vector<double>{ 5.0, 6.0, 7.0, 8.0, 9.0 };

                    test_set << test_buffer;
                }

                // Set 'result #4'
                {
                    ScanFile::DataSet test_set = test_file.add("result #4", 3);
                    ScanFile::WriteBuffer test_buffer(3);

                    ScanFile::Record r = test_file["result #1"][0];
                    for (auto i = 0 ; i < 1004 ; ++i, ++r)
                    {
                        test_buffer << r;
                    }

                    test_set << test_buffer;
                }
            }

            // Open file
            {
                ScanFile test_file = ScanFile::Open(filename);
                TEST_CHECK_EQUAL(test_file.creator(),     "scan_file_TEST");
                TEST_CHECK_EQUAL(test_file.eos_version(), EOS_GITHEAD);

                static const std::vector<std::tuple<std::string, unsigned, unsigned>> references
                {
                    std::tuple<std::string, unsigned, unsigned>(std::string("result #1"), 3, 1004),
                    std::tuple<std::string, unsigned, unsigned>(std::string("result #2"), 7, 1),
                    std::tuple<std::string, unsigned, unsigned>(std::string("result #3"), 5, 1024),
                    std::tuple<std::string, unsigned, unsigned>(std::string("result #4"), 3, 1004),
                };

                auto r = references.cbegin();
                for (ScanFile::Iterator d = test_file.begin(), d_end = test_file.end() ; d != d_end ; ++d, ++r)
                {
                    TEST_CHECK_EQUAL(std::get<0>(*r), d->name());
                    TEST_CHECK_EQUAL(std::get<1>(*r), d->fields());
                    TEST_CHECK_EQUAL(std::get<2>(*r), d->records());
                }

                // "result #1"
                {
                    ScanFile::DataSet test_set = test_file["result #1"];
                    ScanFile::Record test_record = test_set[0];

                    TEST_CHECK_EQUAL(3.0, test_record[0]);
                    TEST_CHECK_EQUAL(2.0, test_record[1]);
                    TEST_CHECK_EQUAL(1.0, test_record[2]);

                    unsigned idx = 1;
                    for (auto f = test_set.begin_fields(), f_end = test_set.end_fields() ; f != f_end ; ++f, ++idx)
                    {
                        std::cout << "f->name = " << f->name() << std::endl;

                        // test the name
                        TEST_CHECK_EQUAL("result #1, field #" + stringify(idx), f->name());

                        // test direct attribute access
                        TEST_CHECK_EQUAL(0.0,                                   f->get("min", 1.0));
                        TEST_CHECK_EQUAL(1.0 * idx,                             f->get("max", 0.5));
                        TEST_CHECK_EQUAL(false,                                 f->get("nuisance", 17.0));

                        // test iteration over attributes
                        const std::set<std::pair<const std::string, double>> references
                        {
                            std::make_pair(std::string("min"),       0.0),
                            std::make_pair(std::string("max"),       1.0 * idx),
                            std::make_pair(std::string("nuisance"),  false),
                        };

                        for (auto i = f->begin_attributes(), i_end = f->end_attributes() ; i != i_end ; ++i)
                        {
                            TEST_CHECK(references.end() != references.find(*i));
                        }
                        TEST_CHECK_EQUAL(std::distance(f->begin_attributes(), f->end_attributes()), references.size());
                    }
                }

                // "result #2"
                {
                    ScanFile::DataSet test_set = test_file["result #2"];
                    ScanFile::Record test_record = test_set[0];

                    TEST_CHECK_EQUAL(7.0, test_record[0]);
                    TEST_CHECK_EQUAL(6.0, test_record[1]);
                    TEST_CHECK_EQUAL(5.0, test_record[2]);
                    TEST_CHECK_EQUAL(4.0, test_record[3]);
                    TEST_CHECK_EQUAL(3.0, test_record[4]);
                    TEST_CHECK_EQUAL(2.0, test_record[5]);
                    TEST_CHECK_EQUAL(1.0, test_record[6]);
                }

                // "result #3"
                {
                    ScanFile::DataSet test_set = test_file["result #3"];
                    ScanFile::Record test_record = test_set[0];

                    TEST_CHECK_EQUAL(4.0, test_record[0]);
                    TEST_CHECK_EQUAL(5.0, test_record[1]);
                    TEST_CHECK_EQUAL(6.0, test_record[2]);
                    TEST_CHECK_EQUAL(7.0, test_record[3]);
                    TEST_CHECK_EQUAL(8.0, test_record[4]);

                    test_record = test_set[33];

                    TEST_CHECK_EQUAL(4.0, test_record[0]);
                    TEST_CHECK_EQUAL(5.0, test_record[1]);
                    TEST_CHECK_EQUAL(6.0, test_record[2]);
                    TEST_CHECK_EQUAL(7.0, test_record[3]);
                    TEST_CHECK_EQUAL(8.0, test_record[4]);

                    test_record = test_set[1023];

                    TEST_CHECK_EQUAL(5.0, test_record[0]);
                    TEST_CHECK_EQUAL(6.0, test_record[1]);
                    TEST_CHECK_EQUAL(7.0, test_record[2]);
                    TEST_CHECK_EQUAL(8.0, test_record[3]);
                    TEST_CHECK_EQUAL(9.0, test_record[4]);
                }

                // "result #4"
                {
                    ScanFile::DataSet test_set = test_file["result #1"];

                    ScanFile::Record r = test_set[0];
                    TEST_CHECK_EQUAL(3.0, r[0]);
                    TEST_CHECK_EQUAL(2.0, r[1]);
                    TEST_CHECK_EQUAL(1.0, r[2]);
                    ++r;

                    for (auto i = 0 ; i < 1002 ; ++i, ++r)
                    {
                        TEST_CHECK_EQUAL(4.0, r[0]);
                        TEST_CHECK_EQUAL(5.0, r[1]);
                        TEST_CHECK_EQUAL(6.0, r[2]);
                    }

                    TEST_CHECK_EQUAL(7.0, r[0]);
                    TEST_CHECK_EQUAL(8.0, r[1]);
                    TEST_CHECK_EQUAL(9.0, r[2]);
                }
            }
        }
} scan_file_test;

class ScanFileWriteBufferTest :
    public TestCase
{
    public:
        ScanFileWriteBufferTest() :
            TestCase("scan_file_write_buffer_test")
        {
        }

        virtual void run() const
        {
            // empty buffer
            {
                ScanFile::WriteBuffer buffer(4, 0);

                TEST_CHECK_EQUAL(0, buffer.capacity());
                TEST_CHECK_EQUAL(0, buffer.size());
            }

            // random capacity
            {
                ScanFile::WriteBuffer buffer(3, 478594);
                TEST_CHECK_EQUAL(478594, buffer.capacity());
                TEST_CHECK_EQUAL(0,      buffer.size());
            }

            // writing to empty buffer
            {
                ScanFile::WriteBuffer buffer(2, 0);
                for (unsigned i = 0 ; i < 1300 ; ++i)
                {
                    buffer << std::vector<double>{ 9.9, 7.7 };
                }

                TEST_CHECK_EQUAL(2048, buffer.capacity());
                TEST_CHECK_EQUAL(1300, buffer.size());
            }
        }
} scan_file_write_buffer_test;
