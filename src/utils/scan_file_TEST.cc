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
                ScanFile test_file = ScanFile::Create(filename, "scan_file_TEST", 3, 5);

                ScanFile::Tuple tuple = test_file[0];
                tuple[0] = 3.0;
                tuple[1] = 2.0;
                tuple[2] = 1.0;
                tuple.write();
            }

            // Open file
            {
                ScanFile test_file = ScanFile::Open(filename);
                TEST_CHECK_EQUAL(test_file.creator(),     "scan_file_TEST");
                TEST_CHECK_EQUAL(test_file.eos_version(), EOS_GITHEAD);
                TEST_CHECK_EQUAL(test_file.tuple_size(),  3);
                TEST_CHECK_EQUAL(test_file.scan_size(),   5);
            }
        }
} scan_file_test;
