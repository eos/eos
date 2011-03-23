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
#include <tuple>

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
            static const std::string filename("/tmp/scan_file_TEST.scan");

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
                        *f = "result #1, field #" + stringify(idx);
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
                };

                auto r = references.cbegin();
                for (ScanFile::Iterator d = test_file.begin(), d_end = test_file.end() ; d != d_end ; ++d, ++r)
                {
                    TEST_CHECK_EQUAL(std::get<0>(*r), d->name());
                    TEST_CHECK_EQUAL(std::get<1>(*r), d->tuple_size());
                    TEST_CHECK_EQUAL(std::get<2>(*r), d->tuples());
                }

                // "result #1"
                {
                    ScanFile::DataSet test_set = test_file["result #1"];
                    ScanFile::Tuple test_tuple = test_set[0];

                    TEST_CHECK_EQUAL(3.0, test_tuple[0]);
                    TEST_CHECK_EQUAL(2.0, test_tuple[1]);
                    TEST_CHECK_EQUAL(1.0, test_tuple[2]);

                    unsigned idx = 1;
                    for (auto f = test_set.begin_fields(), f_end = test_set.end_fields() ; f != f_end ; ++f, ++idx)
                    {
                        TEST_CHECK_EQUAL("result #1, field #" + stringify(idx), *f);
                    }
                }

                // "result #2"
                {
                    ScanFile::DataSet test_set = test_file["result #2"];
                    ScanFile::Tuple test_tuple = test_set[0];

                    TEST_CHECK_EQUAL(7.0, test_tuple[0]);
                    TEST_CHECK_EQUAL(6.0, test_tuple[1]);
                    TEST_CHECK_EQUAL(5.0, test_tuple[2]);
                    TEST_CHECK_EQUAL(4.0, test_tuple[3]);
                    TEST_CHECK_EQUAL(3.0, test_tuple[4]);
                    TEST_CHECK_EQUAL(2.0, test_tuple[5]);
                    TEST_CHECK_EQUAL(1.0, test_tuple[6]);
                }

                // "result #3"
                {
                    ScanFile::DataSet test_set = test_file["result #3"];
                    ScanFile::Tuple test_tuple = test_set[0];

                    TEST_CHECK_EQUAL(4.0, test_tuple[0]);
                    TEST_CHECK_EQUAL(5.0, test_tuple[1]);
                    TEST_CHECK_EQUAL(6.0, test_tuple[2]);
                    TEST_CHECK_EQUAL(7.0, test_tuple[3]);
                    TEST_CHECK_EQUAL(8.0, test_tuple[4]);

                    test_tuple = test_set[33];

                    TEST_CHECK_EQUAL(4.0, test_tuple[0]);
                    TEST_CHECK_EQUAL(5.0, test_tuple[1]);
                    TEST_CHECK_EQUAL(6.0, test_tuple[2]);
                    TEST_CHECK_EQUAL(7.0, test_tuple[3]);
                    TEST_CHECK_EQUAL(8.0, test_tuple[4]);

                    test_tuple = test_set[1023];

                    TEST_CHECK_EQUAL(5.0, test_tuple[0]);
                    TEST_CHECK_EQUAL(6.0, test_tuple[1]);
                    TEST_CHECK_EQUAL(7.0, test_tuple[2]);
                    TEST_CHECK_EQUAL(8.0, test_tuple[3]);
                    TEST_CHECK_EQUAL(9.0, test_tuple[4]);
                }
            }
        }
} scan_file_test;
