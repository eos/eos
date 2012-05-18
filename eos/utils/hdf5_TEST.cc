/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
 * Copyright (c) 2011 Frederik Beaujean
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
#include <eos/utils/hdf5.hh>

using namespace test;
using namespace eos;

class HDF5FileTest :
    public TestCase
{
    public:
        HDF5FileTest() :
            TestCase("hdf5_file_test")
        {
        }

        virtual void run() const
        {
            static const std::string filename(EOS_BUILDDIR "/eos/utils/hdf5_TEST-file.hdf5");

            std::remove(filename.c_str());

            hdf5::Composite<hdf5::Scalar<double>, hdf5::Array<1, double>> type
            {
                "component",
                hdf5::Scalar<double>("weight"),
                hdf5::Array<1, double>("means", { 4 }),
            };

            // Create file with a data set of composite values
            {
                hdf5::File file = hdf5::File::Create(filename);

                auto data_set = file.create_data_set("/data/1/components", type);

                data_set << std::make_tuple(17.0, std::vector<double>{ 1.0, 2.0, 3.0, 4.0 });

                TEST_CHECK(file.group_exists("/data"));
                TEST_CHECK(file.group_exists("/data/1"));
                TEST_CHECK_EQUAL(file.number_of_objects("/data"), 1);
                TEST_CHECK_EQUAL(file.number_of_objects("/data/1"), 1);
            }

            // Open file
            {
                hdf5::File file = hdf5::File::Open(filename, H5F_ACC_RDWR);
                auto data_set = file.open_data_set("/data/1/components", type);

                TEST_CHECK_EQUAL(data_set.records(), 1);

                std::tuple<double, std::vector<double>> record;
                data_set >> record;

                TEST_CHECK_EQUAL(17.0, std::get<0>(record));
                TEST_CHECK_EQUAL(1.0,  std::get<1>(record)[0]);
                TEST_CHECK_EQUAL(2.0,  std::get<1>(record)[1]);
                TEST_CHECK_EQUAL(3.0,  std::get<1>(record)[2]);
                TEST_CHECK_EQUAL(4.0,  std::get<1>(record)[3]);

                // append to the end
                data_set << std::make_tuple(-17.0, std::vector<double>{ -1.0, -2.0, -3.0, -4.0 });
            }

            // Open again
            {
                hdf5::File file = hdf5::File::Open(filename, H5F_ACC_RDONLY);
                auto data_set = file.open_data_set("/data/1/components", type);

                TEST_CHECK_EQUAL(data_set.records(), 2);

                std::tuple<double, std::vector<double>> record;

                data_set >> record;
                TEST_CHECK_EQUAL(17.0, std::get<0>(record));

                data_set >> record;
                TEST_CHECK_EQUAL(-17.0, std::get<0>(record));
            }
            // Non-existing file
            {
                H5E_BEGIN_TRY
                {
                    TEST_CHECK_THROWS(HDF5Error, hdf5::File::Open(filename + ".har", H5F_ACC_RDONLY));
                }
                H5E_END_TRY;
            }

            // create a copy of a group in another file
            {
                hdf5::File file = hdf5::File::Open(filename, H5F_ACC_RDONLY);

                static const std::string filename_copy(EOS_BUILDDIR "/eos/utils/hdf5_TEST-copy.hdf5");

                hdf5::File file_copy = hdf5::File::Create(filename_copy);

                static const std::string directory = "/data/1";

                file.copy(directory, file_copy);
                TEST_CHECK(file_copy.group_exists(directory));

                auto data_set = file_copy.open_data_set("/data/1/components", type);

                std::tuple<double, std::vector<double>> record;

                data_set >> record;
                TEST_CHECK_EQUAL(17.0, std::get<0>(record));

                data_set >> record;
                TEST_CHECK_EQUAL(-17.0, std::get<0>(record));
            }
        }
} hdf5_file_test;

class HDF5AttributeTest:
    public TestCase
{
    public:
        HDF5AttributeTest() :
            TestCase("hdf5_attribute_test")
        {
        }

        virtual void run() const
        {
            static const std::string filename(EOS_BUILDDIR "/eos/utils/hdf5_TEST-attribute.hdf5");
            std::remove(filename.c_str());

            hdf5::Composite<hdf5::Scalar<double>, hdf5::Array<1, double>> record_type
            {
                "component",
                hdf5::Scalar<double>("weight"),
                hdf5::Array<1, double>("means", { 4 }),
            };

            hdf5::Scalar<int> attribute_type("converged");

            // Create data set with a data set of composite values
            {
                hdf5::File file = hdf5::File::Create(filename);
                {
                    auto data_set = file.create_data_set("/components", record_type);
                    auto attr = data_set.create_or_open_attribute("converged", attribute_type);
                    attr = 1;
                }
                {
                    auto data_set = file.open_data_set("/components", record_type);
                    auto attr = data_set.open_attribute("converged", attribute_type);
                    TEST_CHECK_EQUAL(1, attr.value());
                }
            }
        }
} hdf5_attribute_test;

