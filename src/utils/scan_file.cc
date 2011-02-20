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
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/scan_file.hh>
#include <src/utils/stringify.hh>

#include <hdf5.h>

#include <iostream>

namespace eos
{
    ScanFileError::ScanFileError(const std::string & message) :
        Exception(message)
    {
    }

    ScanFileHDF5Error::ScanFileHDF5Error(const std::string & function, int code) :
        ScanFileError("Error when calling '" + function + "': error code is " + stringify(code))
    {
    }

    template <>
    struct Implementation<ScanFile>
    {
        const hid_t file_id;

        hid_t group_id_metadata;

        hid_t space_id_scan;

        hid_t set_id_scan;

        std::string creator;

        std::string eos_version;

        hsize_t tuple_size;

        hsize_t scan_size;

        // Create a new HDF5 file
        Implementation(const std::string & filename, const std::string & creator,
                const unsigned & tuple_size, const unsigned & scan_size) :
            file_id(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)),
            group_id_metadata(H5I_INVALID_HID),
            space_id_scan(H5I_INVALID_HID),
            set_id_scan(H5I_INVALID_HID),
            creator(creator),
            eos_version(EOS_GITHEAD),
            tuple_size(tuple_size),
            scan_size(scan_size)
        {
            if (H5I_INVALID_HID == file_id)
                throw ScanFileHDF5Error("H5Fcreate", file_id);

            group_id_metadata = H5Gcreate2(file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (H5I_INVALID_HID == group_id_metadata)
                throw ScanFileHDF5Error("H5Gcreate2", group_id_metadata);

            add_metadata_field("/metadata/creator", creator);
            add_metadata_field("/metadata/eos_version", EOS_GITHEAD);

            hsize_t dimensions[2] = { scan_size, tuple_size };

            space_id_scan = H5Screate_simple(2, dimensions, 0);
            if (H5I_INVALID_HID == space_id_scan)
                throw ScanFileHDF5Error("H5Screate_simple", space_id_scan);

            set_id_scan = H5Dcreate2(file_id, "/scan", H5T_IEEE_F64LE, space_id_scan, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (H5I_INVALID_HID == set_id_scan)
                throw ScanFileHDF5Error("H5Dcreate2", set_id_scan);
        }

        Implementation(const std::string & filename) :
            file_id(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT))
        {
            if (H5I_INVALID_HID == file_id)
                throw ScanFileHDF5Error("H5Fopen", file_id);

            group_id_metadata = H5Gopen2(file_id, "/metadata", H5P_DEFAULT);
            if (H5I_INVALID_HID == group_id_metadata)
                throw ScanFileHDF5Error("H5Gopen2", group_id_metadata);

            creator = read_metadata_field("/metadata/creator");
            eos_version = read_metadata_field("/metadata/eos_version");

            set_id_scan = H5Dopen2(file_id, "/scan", H5P_DEFAULT);
            if (H5I_INVALID_HID == set_id_scan)
                throw ScanFileHDF5Error("H5Dopen2", set_id_scan);

            space_id_scan = H5Dget_space(set_id_scan);
            if (H5I_INVALID_HID == space_id_scan)
                throw ScanFileHDF5Error("H5Dget_space", space_id_scan);

            if (0 >= H5Sis_simple(space_id_scan))
                throw ScanFileError("File at hand is not an EOS scan file: '/scan' is not associated to a simple data space: " + stringify(H5Sis_simple(space_id_scan)));

            if (2 != H5Sget_simple_extent_ndims(space_id_scan))
                throw ScanFileError("File at hand is not an EOS scan file: '/scan' is not of rank 2");

            hsize_t dimensions[2] = { -1, -1 };
            if (2 != H5Sget_simple_extent_dims(space_id_scan, dimensions, 0))
                throw ScanFileError("File at hand is fishy");

            scan_size = dimensions[0];
            tuple_size = dimensions[1];
        }

        ~Implementation()
        {
            herr_t ret;

            ret = H5Dclose(set_id_scan);
            if (0 > ret)
                throw ScanFileHDF5Error("H5Dclose", ret);

            ret = H5Sclose(space_id_scan);
            if (0 > ret)
                throw ScanFileHDF5Error("H5Sclose", ret);

            ret = H5Gclose(group_id_metadata);
            if (0 > ret)
                throw ScanFileHDF5Error("H5Gclose", ret);

            ret = H5Fclose(file_id);
            if (0 > ret)
                throw ScanFileHDF5Error("H5Fclose", ret);
        }

        void add_metadata_field(const std::string & field_name, const std::string & value)
        {
            hid_t type_id = H5Tcopy(H5T_C_S1);
            if (H5I_INVALID_HID == type_id)
                throw ScanFileHDF5Error("H5Tcopy", type_id);

            hid_t space_id = H5Screate(H5S_SCALAR);
            if (H5I_INVALID_HID == space_id)
                throw ScanFileHDF5Error("H5Screate", space_id);

            herr_t ret = H5Tset_size(type_id, value.size());
            if (ret > 0)
                throw ScanFileHDF5Error("H5Tset_size", ret);

            hid_t set_id = H5Dcreate(group_id_metadata, field_name.c_str(), type_id, space_id, H5P_DEFAULT);
            if (H5I_INVALID_HID == set_id)
                throw ScanFileHDF5Error("H5Dcreate", set_id);

            ret = H5Dwrite(set_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value[0]);

            ret = H5Dclose(set_id);
            if (ret > 0)
                throw ScanFileHDF5Error("H5Dclose", ret);
        }

        std::string read_metadata_field(const std::string & field_name)
        {
            hid_t set_id = H5Dopen(group_id_metadata, field_name.c_str());
            if (H5I_INVALID_HID == set_id)
                throw ScanFileHDF5Error("H5Dopen", set_id);

            hsize_t set_size = H5Dget_storage_size(set_id);
            std::string result = std::string(set_size, '\0');

            hid_t type_id = H5Tcopy(H5T_C_S1);
            if (H5I_INVALID_HID == type_id)
                throw ScanFileHDF5Error("H5Tcopy", type_id);

            herr_t ret = H5Tset_size(type_id, set_size);
            if (ret > 0)
                throw ScanFileHDF5Error("H5Tset_size", ret);

            ret = H5Dread(set_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &result[0]);
            if (ret > 0)
                throw ScanFileHDF5Error("H5Dread", ret);

            ret = H5Dclose(set_id);
            if (ret > 0)
                throw ScanFileHDF5Error("H5Dclose", ret);

            return result;
        }
    };

    ScanFile::ScanFile(Implementation<ScanFile> * imp) :
        PrivateImplementationPattern<ScanFile>(imp)
    {
    }

    ScanFile::~ScanFile()
    {
    }

    ScanFile
    ScanFile::Create(const std::string & filename, const std::string & creator,
            const unsigned & tuple_size, const unsigned & scan_size)
    {
        return ScanFile(new Implementation<ScanFile>(filename, creator, tuple_size, scan_size));
    }

    ScanFile
    ScanFile::Open(const std::string & filename)
    {
        return ScanFile(new Implementation<ScanFile>(filename));
    }

    const std::string &
    ScanFile::creator() const
    {
        return _imp->creator;
    }

    const std::string &
    ScanFile::eos_version() const
    {
        return _imp->eos_version;
    }

    int
    ScanFile::tuple_size() const
    {
        return _imp->tuple_size;
    }

    int
    ScanFile::scan_size() const
    {
        return _imp->scan_size;
    }

    ScanFile::Tuple
    ScanFile::operator[] (const unsigned & index)
    {
        return ScanFile::Tuple(_imp, index);
    }

    template <>
    struct Implementation<ScanFile::Tuple>
    {
        std::shared_ptr<Implementation<ScanFile>> file_imp;

        const unsigned tuple_size;

        unsigned index;

        double * buffer;

        hid_t set_id_scan;

        hid_t space_id_memory;

        hid_t space_id_file;

        Implementation(const std::shared_ptr<Implementation<ScanFile>> & file_imp, const unsigned & index) :
            file_imp(file_imp),
            tuple_size(file_imp->tuple_size),
            index(index),
            set_id_scan(file_imp->set_id_scan)
        {
            buffer = new double[tuple_size];
            hsize_t dimensions[2] = { 1, tuple_size };
            space_id_memory = H5Screate_simple(2, dimensions, 0);
            if (H5I_INVALID_HID == space_id_memory)
                throw ScanFileHDF5Error("H5Screate_simple", space_id_memory);

            space_id_file = H5Dget_space(set_id_scan);
            if (H5I_INVALID_HID == space_id_file)
                throw ScanFileHDF5Error("H5Dget_space", space_id_file);

            select();
        }

        ~Implementation()
        {
            delete[] buffer;

            herr_t ret = H5Sclose(space_id_memory);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sclose", ret);

            ret = H5Sclose(space_id_file);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sclose", ret);
        }

        void select()
        {
            hsize_t offset[2] = { index, 0 };
            hsize_t count[2] = { 1, tuple_size };

            herr_t ret = H5Sselect_hyperslab(space_id_file, H5S_SELECT_SET, offset, 0, count, 0);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sselect_hyperslab", ret);
        }

        void read()
        {
            herr_t ret = H5Dread(set_id_scan, H5T_IEEE_F64LE, space_id_memory, space_id_file, H5P_DEFAULT, buffer);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Dread", ret);
        }

        void write()
        {
            herr_t ret = H5Dwrite(set_id_scan, H5T_IEEE_F64LE, space_id_memory, space_id_file, H5P_DEFAULT, buffer);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Dwrite", ret);
        }
    };

    ScanFile::Tuple::Tuple(const std::shared_ptr<Implementation<ScanFile>> & file_imp, const unsigned & index) :
        PrivateImplementationPattern<ScanFile::Tuple>(new Implementation<ScanFile::Tuple>(file_imp, index))
    {

    }

    ScanFile::Tuple::~Tuple()
    {
    }

    ScanFile::Tuple &
    ScanFile::Tuple::operator++ ()
    {
        _imp->index += 1;
        _imp->select();

        return *this;
    }

    double &
    ScanFile::Tuple::operator[] (const unsigned & index)
    {
        return _imp->buffer[index];
    }

    void
    ScanFile::Tuple::read()
    {
        _imp->read();
    }

    void
    ScanFile::Tuple::write()
    {
        _imp->write();
    }
}
