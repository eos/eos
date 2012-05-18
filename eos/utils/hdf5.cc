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
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/stringify.hh>

#include <hdf5.h>

namespace eos
{
    HDF5Error::HDF5Error(const std::string & message) :
        Exception(message)
    {
    }

    template <> struct Implementation<hdf5::FileHandle>
    {
        const hid_t file_id;

        bool read_only;

        Implementation(const hid_t & file_id, bool read_only) :
            file_id(file_id),
            read_only(read_only)
        {
        }

        ~Implementation()
        {
            herr_t ret;

            if (! read_only)
            {
                ret = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
                if (0 > ret)
                    throw HDF5Error("H5Fflush failed and returned " + stringify(ret));
            }

            ret = H5Fclose(file_id);
            if (0 > ret)
                throw HDF5Error("H5Fclose failed and returned " + stringify(ret));
        }
    };

    template <> struct Implementation<hdf5::DataSetHandle>
    {
        hdf5::FileHandle file_handle;

        hid_t data_set_id;

        hid_t space_id_file;

        hid_t space_id_memory_element;

        hid_t type_id;

        hsize_t size;

        hsize_t capacity;

        Implementation(const hdf5::FileHandle & file_handle, const hid_t & data_set_id, const hid_t & space_id_file, hsize_t size) :
            file_handle(file_handle),
            data_set_id(data_set_id),
            space_id_file(space_id_file),
            space_id_memory_element(H5Screate(H5S_SCALAR)),
            type_id(H5Dget_type(data_set_id)),
            size(size),
            capacity(H5Sget_simple_extent_npoints(space_id_file))
        {
        }

        ~Implementation()
        {
            herr_t ret;

            // truncate the data set to its actual size
            if (! file_handle.read_only())
            {
                hsize_t dimension = size;
                ret = H5Dset_extent(data_set_id, &dimension);
                if (0 > ret)
                    throw HDF5Error("H5Dset_extent failed and returned " + stringify(ret));
            }

            ret = H5Sclose(space_id_memory_element);
            if (0 > ret)
                throw HDF5Error("H5Sclose failed and returned " + stringify(ret));

            ret = H5Sclose(space_id_file);
            if (0 > ret)
                throw HDF5Error("H5Sclose failed and returned " + stringify(ret));

            ret = H5Tclose(type_id);
            if (0 > ret)
                throw HDF5Error("H5Tclose failed and returned " + stringify(ret));

            ret = H5Dclose(data_set_id);
            if (0 > ret)
                throw HDF5Error("H5Dclose failed and returned " + stringify(ret));
        }
    };

    template <> struct Implementation<hdf5::AttributeHandle>
    {
        hdf5::DataSetHandle data_set_handle;

        hid_t attribute_id;

        hid_t type_id;

        hid_t space_id_file;

        hid_t space_id_memory;

        Implementation(const hdf5::DataSetHandle & data_set_handle, const hid_t & attribute_id) :
            data_set_handle(data_set_handle),
            attribute_id(attribute_id),
            type_id(H5Aget_type(attribute_id)),
            space_id_file(H5Aget_space(attribute_id)),
            space_id_memory(H5Screate(H5S_SCALAR))
        {
        }

        ~Implementation()
        {
            herr_t ret;

            ret = H5Sclose(space_id_memory);
            if (0 > ret)
                throw HDF5Error("H5Sclose failed and returned " + stringify(ret));

            ret = H5Sclose(space_id_file);
            if (0 > ret)
                throw HDF5Error("H5Sclose failed and returned " + stringify(ret));

            ret = H5Tclose(type_id);
            if (0 > ret)
                throw HDF5Error("H5Tclose failed and returned " + stringify(ret));

            ret = H5Aclose(attribute_id);
            if (0 > ret)
                throw HDF5Error("H5Aclose failed and returned " + stringify(ret));
        }
    };

    namespace hdf5
    {
        /* Type System */
        Type::~Type()
        {
        }

        /* Representation Classes */

        FileHandle::FileHandle(const hid_t & file_id, bool read_only) :
            PrivateImplementationPattern<hdf5::FileHandle>(new Implementation<hdf5::FileHandle>(file_id, read_only))
        {
        }

        FileHandle::~FileHandle()
        {
        }

        hid_t
        FileHandle::id() const
        {
            return _imp->file_id;
        }

        bool
        FileHandle::read_only() const
        {
            return _imp->read_only;
        }

        DataSetHandle::DataSetHandle(const FileHandle & file_handle, const hid_t & data_set_id, const hid_t & space_id_file, hsize_t size) :
            PrivateImplementationPattern<hdf5::DataSetHandle>(new Implementation<DataSetHandle>(file_handle, data_set_id, space_id_file, size))
        {
        }

        DataSetHandle::~DataSetHandle()
        {
        }

        hid_t
        DataSetHandle::set_id() const
        {
            return _imp->data_set_id;
        }

        hid_t
        DataSetHandle::space_id() const
        {
            return _imp->space_id_file;
        }

        hid_t
        DataSetHandle::type_id() const
        {
            return _imp->type_id;
        }

        hsize_t
        DataSetHandle::size() const
        {
            return _imp->size;
        }

        void
        DataSetHandle::select(hsize_t start, hsize_t count)
        {
            H5Sselect_hyperslab(_imp->space_id_file, H5S_SELECT_SET, &start, 0, &count, 0);
        }

        void
        DataSetHandle::write_one(const void * buffer)
        {
            if (_imp->size + 5 >= _imp->capacity)
            {
                hsize_t new_capacity = _imp->capacity + 1000;
                hsize_t max_capacity = H5S_UNLIMITED;

                herr_t ret = H5Sset_extent_simple(_imp->space_id_file, 1, &new_capacity, &max_capacity);
                if (0 > ret)
                    throw HDF5Error("H5Sset_extent_simple failed and returned " + stringify(ret));

                select(_imp->size, 1);

                ret = H5Dset_extent(_imp->data_set_id, &new_capacity);
                if (0 > ret)
                    throw HDF5Error("H5Dset_extent failed and returned " + stringify(ret));

                _imp->capacity = new_capacity;
            }

            H5Dwrite(_imp->data_set_id, _imp->type_id, _imp->space_id_memory_element, _imp->space_id_file, H5P_DEFAULT, buffer);
            ++_imp->size;
        }

        void
        DataSetHandle::read_one(void * buffer)
        {
            H5Dread(_imp->data_set_id, _imp->type_id, _imp->space_id_memory_element, _imp->space_id_file, H5P_DEFAULT, buffer);
        }

        AttributeHandle
        DataSetHandle::create_attribute(const std::string & name, const hid_t & type_id)
        {
            hid_t space_id = H5Screate(H5S_SCALAR);

            hid_t attr_id = H5Acreate2(_imp->data_set_id, name.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
            if (H5I_INVALID_HID == attr_id)
                throw HDF5Error("H5Acreate2 failed to create '" + name + "'");

            return AttributeHandle(*this, attr_id);
        }

        AttributeHandle
        DataSetHandle::open_attribute(const std::string & name, const hid_t & type_id)
        {
            hid_t attr_id = H5Aopen(_imp->data_set_id, name.c_str(), H5P_DEFAULT);
            if (H5I_INVALID_HID == attr_id)
                throw HDF5Error("H5Acreate2 failed to open '" + name + "'");

            herr_t ret = H5Tequal(H5Aget_type(attr_id), type_id);
            if (0 > ret)
                throw InternalError("Trying to open attribute '" + name + "' of incompatible type");

            return AttributeHandle(*this, attr_id);
        }

        AttributeHandle::AttributeHandle(const DataSetHandle & data_set_handle, const hid_t & attribute_id) :
            PrivateImplementationPattern<hdf5::AttributeHandle>(new Implementation<hdf5::AttributeHandle>(data_set_handle, attribute_id))
        {
        }

        AttributeHandle::~AttributeHandle()
        {
        }

        hid_t
        AttributeHandle::attribute_id() const
        {
            return _imp->attribute_id;
        }

        hid_t
        AttributeHandle::type_id() const
        {
            return _imp->type_id;
        }

        void
        AttributeHandle::write(const void * buffer)
        {
            herr_t ret = H5Awrite(_imp->attribute_id, _imp->type_id, buffer);
            if (0 > ret)
                throw HDF5Error("H5Awrite failed and returned " + stringify(ret));
        }

        void
        AttributeHandle::read(void * buffer) const
        {
            herr_t ret = H5Aread(_imp->attribute_id, _imp->type_id, buffer);
            if (0 > ret)
                throw HDF5Error("H5Awrite failed and returned " + stringify(ret));
        }

        File::File(const FileHandle & handle, const std::string & file_name) :
            _handle(handle),
            _name(file_name)
        {
        }

        File::~File()
        {
        }

        File
        File::Create(const std::string & file_name)
        {
            H5dont_atexit();
            hid_t file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            if (H5I_INVALID_HID == file_id)
                throw HDF5Error("H5Fcreate failed and returned " + stringify(file_id));

            return File(FileHandle(file_id, false), file_name);
        }

        // Open an existing HDF5 file
        File
        File::Open(const std::string & file_name, unsigned mode)
        {
            H5dont_atexit();
            hid_t file_id = H5Fopen(file_name.c_str(), mode, H5P_DEFAULT);
            if (H5I_INVALID_HID == file_id)
                throw HDF5Error("H5Fopen failed to open '" + file_name + "' and returned " + stringify(file_id));

            return File(FileHandle(file_id, mode == H5F_ACC_RDONLY), file_name);
        }

        bool
        File::Exists(const std::string & file_name)
        {
            H5E_BEGIN_TRY
            {
                hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
                if (H5I_INVALID_HID == file_id)
                    return false;
            }
            H5E_END_TRY;

            return true;
        }

        DataSetHandle
        File::_create_data_set(const std::string & name, const TypePtr & type)
        {
            static const unsigned capacity = 10;
            hid_t space_id_file, dcpl_id, lcpl_id, set_id;

            // create space id for in-file representation
            {
                hsize_t dimensions = capacity;
                hsize_t max_dimensions = H5S_UNLIMITED;
                space_id_file = H5Screate_simple(1, &dimensions, &max_dimensions);
                if (H5I_INVALID_HID == space_id_file)
                    throw HDF5Error("H5Screate_simple failed and returned " + stringify(space_id_file));
            }

            // create a property list for in-file data
            {
                dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
                if (H5I_INVALID_HID == dcpl_id)
                    throw HDF5Error("H5Pcreate failed and returned " + stringify(dcpl_id));

                hsize_t chunk_size = capacity;
                herr_t ret = H5Pset_chunk(dcpl_id, 1, &chunk_size);
                if (0 > ret)
                    throw HDF5Error("H5Pset_chunk failed and returned " + stringify(dcpl_id));

                lcpl_id = H5Pcreate(H5P_LINK_CREATE);
                if (H5I_INVALID_HID == lcpl_id)
                    throw HDF5Error("H5Pcreate failed and returned " + stringify(lcpl_id));

                ret = H5Pset_create_intermediate_group(lcpl_id, 1);
                if (0 > ret)
                    throw HDF5Error("H5Pset_create_intermediate_group failed and returned " + stringify(lcpl_id));
            }

            // create set id for in-file data
            {
                set_id = H5Dcreate2(_handle.id(), name.c_str(), type->type_id(), space_id_file, lcpl_id, dcpl_id, H5P_DEFAULT);
                if (H5I_INVALID_HID == set_id)
                    throw HDF5Error("H5Dcreate2 failed to create '" + name + "' and returned " + stringify(set_id));
            }

            return DataSetHandle(_handle, set_id, space_id_file, 0);
        }

        DataSetHandle
        File::_open_data_set(const std::string & name, const TypePtr & type) const
        {
            hid_t space_id_file, set_id;
            hsize_t size = 0;

            // open set by name and retrieve its space id
            {
                set_id = H5Dopen2(_handle.id(), name.c_str(), H5P_DEFAULT);
                if (H5I_INVALID_HID == set_id)
                throw HDF5Error("H5Dopen2 failed to open '" + name + "' and returned " + stringify(set_id));

                space_id_file = H5Dget_space(set_id);
                if (H5I_INVALID_HID == space_id_file)
                    throw HDF5Error("H5Dget_space failed and returned " + stringify(space_id_file));

                hid_t set_type_id = H5Dget_type(set_id);
                if (0 >= H5Tequal(set_type_id, type->type_id()))
                    throw HDF5Error("DataSet<>'s type id and HDF5's type_id for data set '" + name + "' do not match");

                if (1 != H5Sget_simple_extent_ndims(space_id_file))
                    throw HDF5Error("H5Dget_simple_extent_ndims reports rank != 1 for DataSet<>");

                if (0 >= H5Sget_simple_extent_dims(space_id_file, &size, 0))
                    throw HDF5Error("H5Sget_simple_extent_dims failed");
            }

            return DataSetHandle(_handle, set_id, space_id_file, size);
        }

        void
        File::copy(const std::string& source_directory, File& destination_file, const std::string& destination_directory) const
        {
            std::string dest = (destination_directory.empty()) ? source_directory : destination_directory;

            // create intermediate groups if needed
            hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);
            H5Pset_create_intermediate_group(lcpl_id, 1);

            herr_t ret = H5Ocopy(_handle.id(), source_directory.c_str(), destination_file._handle.id(),
                                 dest.c_str(), H5P_OBJECT_COPY_DEFAULT, lcpl_id);

            if (0 > ret)
                throw HDF5Error("H5Ocopy failed to copy " + source_directory + " to "
                                + destination_directory + " in file " + destination_file.name()
                                + " and returned " + stringify(ret));
        }

        bool
        File::group_exists(const std::string & name)
        {
            H5E_BEGIN_TRY
            {
            hid_t group_id = H5Gopen2(_handle.id(), name.c_str(), H5P_DEFAULT);
            if (H5I_INVALID_HID == group_id)
                return false;
            }
            H5E_END_TRY;

            return true;
        }

        unsigned
        File::number_of_objects(const std::string & name)
        {
            H5G_info_t info;

            herr_t ret = H5Gget_info_by_name(_handle.id(), name.c_str(), &info, H5P_DEFAULT);
            if (0 > ret)
                throw HDF5Error("H5Gget_info_by_name(" + name + ") failed and returned " + stringify(ret));

            return info.nlinks;
        }
    }
}
