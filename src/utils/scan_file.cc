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
#include <src/utils/wrapped_forward_iterator-impl.hh>

#include <algorithm>

#include <hdf5.h>
#include <hdf5_hl.h>

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

    template class WrappedForwardIterator<ScanFile::IteratorTag, ScanFile::DataSet>;

    template <>
    struct Implementation<ScanFile>
    {
        const hid_t file_id;

        hid_t group_id_metadata;

        hid_t group_id_data;

        std::string creator;

        std::string eos_version;

        std::vector<ScanFile::DataSet> data_sets;

        bool read_only;

        // Create a new HDF5 file
        Implementation(const std::string & filename, const std::string & creator) :
            file_id(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)),
            group_id_metadata(H5I_INVALID_HID),
            group_id_data(H5I_INVALID_HID),
            creator(creator),
            eos_version(EOS_GITHEAD),
            read_only(false)
        {
            if (H5I_INVALID_HID == file_id)
                throw ScanFileHDF5Error("H5Fcreate", file_id);

            group_id_metadata = H5Gcreate2(file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (H5I_INVALID_HID == group_id_metadata)
                throw ScanFileHDF5Error("H5Gcreate2", group_id_metadata);

            add_metadata_field("/metadata/creator", creator);
            add_metadata_field("/metadata/eos_version", EOS_GITHEAD);

            group_id_data = H5Gcreate2(file_id, "/data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (H5I_INVALID_HID == group_id_data)
                throw ScanFileHDF5Error("H5Gcreate2", group_id_data);
        }

        // Callback function for iteration over data sets
        static herr_t add_data_set(hid_t /* unused */, const char * name, const H5L_info_t * /* unused */, void * _imp)
        {
            Implementation * imp = static_cast<Implementation<ScanFile> *>(_imp);

            imp->data_sets.push_back(ScanFile::DataSet(imp, std::string(name)));

            return 0;
        }

        // Open an existing HDF5 file
        Implementation(const std::string & filename) :
            file_id(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)),
            group_id_metadata(H5I_INVALID_HID),
            group_id_data(H5I_INVALID_HID),
            read_only(true)
        {
            if (H5I_INVALID_HID == file_id)
                throw ScanFileHDF5Error("H5Fopen", file_id);

            group_id_metadata = H5Gopen2(file_id, "/metadata", H5P_DEFAULT);
            if (H5I_INVALID_HID == group_id_metadata)
                throw ScanFileHDF5Error("H5Gopen2", group_id_metadata);

            creator = read_metadata_field("/metadata/creator");
            eos_version = read_metadata_field("/metadata/eos_version");

            group_id_data = H5Gopen2(file_id, "/data", H5P_DEFAULT);
            if (H5I_INVALID_HID == group_id_data)
                throw ScanFileHDF5Error("H5Gopen2", group_id_data);

            hsize_t idx = 0;
            herr_t ret = H5Literate(group_id_data, H5_INDEX_NAME, H5_ITER_INC, &idx, &Implementation<ScanFile>::add_data_set, this);
            if (0 > ret)
                throw ScanFileHDF5Error("H5Literate", ret);
        }

        ~Implementation()
        {
            herr_t ret;

            // Close all data sets
            data_sets.clear();

            ret = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
            if (0 > ret)
                throw ScanFileHDF5Error("H5Fflush", ret);

            ret = H5Gclose(group_id_data);
            if (0 > ret)
                throw ScanFileHDF5Error("H5Gclose", ret);

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

            hid_t set_id = H5Dcreate2(group_id_metadata, field_name.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (H5I_INVALID_HID == set_id)
                throw ScanFileHDF5Error("H5Dcreate2", set_id);

            ret = H5Dwrite(set_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value[0]);

            ret = H5Dclose(set_id);
            if (ret > 0)
                throw ScanFileHDF5Error("H5Dclose", ret);
        }

        std::string read_metadata_field(const std::string & field_name)
        {
            hid_t set_id = H5Dopen2(group_id_metadata, field_name.c_str(), H5P_DEFAULT);
            if (H5I_INVALID_HID == set_id)
                throw ScanFileHDF5Error("H5Dopen2", set_id);

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

        static bool same_name_as(const std::string & name, const ScanFile::DataSet & data_set)
        {
            return name == data_set.name();
        }

        std::vector<ScanFile::DataSet>::iterator find_data_set(const std::string & name)
        {
            return std::find_if(data_sets.begin(), data_sets.end(), std::bind(Implementation<ScanFile>::same_name_as, name, std::placeholders::_1));
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
    ScanFile::Create(const std::string & filename, const std::string & creator)
    {
        H5dont_atexit();
        return ScanFile(new Implementation<ScanFile>(filename, creator));
    }

    ScanFile
    ScanFile::Open(const std::string & filename)
    {
        H5dont_atexit();
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

    ScanFile::DataSet
    ScanFile::operator[] (const std::string & name)
    {
        auto d = _imp->find_data_set(name);

        if (_imp->data_sets.end() == d)
            throw ScanFileError("No such data set '" + name + "'");

        return *d;
    }

    ScanFile::DataSet
    ScanFile::add(const std::string & name, unsigned tuple_size)
    {
        ScanFile::DataSet result = ScanFile::DataSet(_imp.get(), name, tuple_size);
        _imp->data_sets.push_back(result);

        return result;
    }

    ScanFile::Iterator
    ScanFile::begin()
    {
        return ScanFile::Iterator(_imp->data_sets.begin());
    }

    ScanFile::Iterator
    ScanFile::end()
    {
        return ScanFile::Iterator(_imp->data_sets.end());
    }

    template <>
    struct Implementation<ScanFile::WriteBuffer>
    {
        static const unsigned chunk_size = 1024;

        unsigned tuples;

        unsigned tuple_size;

        std::vector<double> buffer;

        Implementation(const unsigned & tuple_size) :
            tuples(0),
            tuple_size(tuple_size)
        {
            buffer.reserve(tuple_size * chunk_size);
        }

        void clear()
        {
            tuples = 0;
            buffer.clear(); // clear only affects size, not capacity
        }
    };

    ScanFile::WriteBuffer::WriteBuffer(const unsigned & tuple_size) :
        PrivateImplementationPattern<ScanFile::WriteBuffer>(new Implementation<ScanFile::WriteBuffer>(tuple_size))
    {
    }

    ScanFile::WriteBuffer::~WriteBuffer()
    {
    }

    void
    ScanFile::WriteBuffer::clear()
    {
        _imp->clear();
    }

    unsigned
    ScanFile::WriteBuffer::capacity() const
    {
        return Implementation<ScanFile::WriteBuffer>::chunk_size;
    }

    unsigned
    ScanFile::WriteBuffer::size() const
    {
        return _imp->tuples;
    }

    ScanFile::WriteBuffer &
    operator<< (ScanFile::WriteBuffer & lhs, const std::vector<double> & rhs)
    {
        if (Implementation<ScanFile::WriteBuffer>::chunk_size == lhs._imp->tuples)
            throw InternalError("Extending WriteBuffer capacity is not yet implemented");

        for (auto i = rhs.cbegin(), i_end = rhs.cend() ; i != i_end ; ++i)
            lhs._imp->buffer.push_back(*i);

        ++lhs._imp->tuples;

        return lhs;
    }

    template class WrappedForwardIterator<ScanFile::DataSet::FieldIteratorTag, std::string>;

    template <>
    struct Implementation<ScanFile::DataSet>
    {
        Implementation<ScanFile> * file_imp;

        std::string name;

        unsigned tuple_size;

        unsigned current_index;

        unsigned tuples;

        unsigned capacity;

        static const unsigned chunk_size = 1024;

        hid_t group_id_data;

        hid_t set_id;

        hid_t dcpl_id;

        hid_t space_id_memory;

        hid_t space_id_file;

        hid_t space_id_file_writing;

        std::vector<std::string> field_names;

        bool truncate_on_destroy;

        // Constructor to create a new data set
        Implementation(Implementation<ScanFile> * file_imp, const std::string & name, unsigned tuple_size) :
            file_imp(file_imp),
            name(name),
            tuple_size(tuple_size),
            current_index(0),
            tuples(0),
            capacity(chunk_size),
            group_id_data(file_imp->group_id_data),
            set_id(H5I_INVALID_HID),
            space_id_memory(H5I_INVALID_HID),
            space_id_file(H5I_INVALID_HID),
            space_id_file_writing(H5I_INVALID_HID),
            field_names(tuple_size, ""),
            truncate_on_destroy(true)
        {
            // create space id for in-memory representation
            {
                hsize_t dimensions[2] = { 1, tuple_size };
                hsize_t max_dimensions[2] = { capacity, tuple_size };
                space_id_memory = H5Screate_simple(2, dimensions, max_dimensions);
                if (H5I_INVALID_HID == space_id_memory)
                    throw ScanFileHDF5Error("H5Screate_simple", space_id_memory);
            }

            // create space id for in-file representation
            {
                hsize_t dimensions[2] = { capacity, tuple_size };
                hsize_t max_dimensions[2] = { H5S_UNLIMITED, tuple_size };
                space_id_file = H5Screate_simple(2, dimensions, max_dimensions);
                if (H5I_INVALID_HID == space_id_file)
                    throw ScanFileHDF5Error("H5Screate_simple", space_id_file);

                space_id_file_writing = H5Scopy(space_id_file);
                if (H5I_INVALID_HID == space_id_file_writing)
                    throw ScanFileHDF5Error("H5Scopy", space_id_file_writing);
            }

            // create a property list for in-file data
            {
                dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
                if (H5I_INVALID_HID == dcpl_id)
                    throw ScanFileHDF5Error("H5Pcreate", dcpl_id);

                hsize_t chunk_size[2] = { capacity, tuple_size };
                herr_t ret = H5Pset_chunk(dcpl_id, 2, chunk_size);
                if (0 > ret)
                    throw ScanFileHDF5Error("H5Pset_chunk", ret);
            }

            // create set id for in-file data
            {
                set_id = H5Dcreate2(group_id_data, name.c_str(), H5T_IEEE_F64LE, space_id_file, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
                if (H5I_INVALID_HID == set_id)
                    throw ScanFileHDF5Error("H5Dcreate2", set_id);
            }
        }

        // Constructor to open an existing data set in append mode
        Implementation(Implementation<ScanFile> * file_imp, const std::string & name) :
            file_imp(file_imp),
            name(name),
            tuple_size(0),
            current_index(0),
            group_id_data(file_imp->group_id_data),
            set_id(H5I_INVALID_HID),
            space_id_memory(H5I_INVALID_HID),
            space_id_file(H5I_INVALID_HID),
            truncate_on_destroy(false)
        {
            // open set by name and retrieve its space id
            {
                set_id = H5Dopen2(group_id_data, name.c_str(), H5P_DEFAULT);
                if (H5I_INVALID_HID == set_id)
                    throw ScanFileHDF5Error("H5Dopen2", set_id);

                space_id_file = H5Dget_space(set_id);
                if (H5I_INVALID_HID == space_id_file)
                    throw ScanFileHDF5Error("H5Dget_space", space_id_file);

                space_id_file_writing = H5Scopy(space_id_file);
                if (H5I_INVALID_HID == space_id_file_writing)
                    throw ScanFileHDF5Error("H5Scopy", space_id_file_writing);

                if (0 >= H5Sis_simple(space_id_file))
                    throw ScanFileError("File at hand is not an EOS scan file: '/data/" + name + "' is not associated to a simple data space: " + stringify(H5Sis_simple(space_id_file)));

                if (2 != H5Sget_simple_extent_ndims(space_id_file))
                    throw ScanFileError("File at hand is not an EOS scan file: '/data/" + name + "' is not a data space of rank 2");

                hsize_t dimensions[2] = { -1, -1 };
                if (2 != H5Sget_simple_extent_dims(space_id_file, dimensions, 0))
                    throw ScanFileError("File at hand is fishy");

                current_index = capacity = tuples = dimensions[0];
                tuple_size = dimensions[1];
            }

            // create space id for in-memory representation
            {
                hsize_t dimensions[2] = { 1, tuple_size };
                hsize_t max_dimensions[2] = { chunk_size, tuple_size };
                space_id_memory = H5Screate_simple(2, dimensions, max_dimensions);
                if (H5I_INVALID_HID == space_id_memory)
                    throw ScanFileHDF5Error("H5Screate_simple", space_id_memory);
            }

            // read field names from attributes
            for (unsigned i = 0 ; i < tuple_size ; ++i)
            {
                std::string attr_name = "FIELD_" + stringify(i) + "_NAME";
                hsize_t attr_rank;
                H5T_class_t attr_class;
                size_t attr_size;

                // how large is the attribute's value?
                herr_t ret = H5LTget_attribute_info(group_id_data, name.c_str(), attr_name.c_str(), &attr_rank, &attr_class, &attr_size);
                if (ret < 0)
                    throw ScanFileHDF5Error("H5LTget_attribute_info", ret);

                // allocate a suitable buffer
                // fuck the HDF5 Lite API!
                std::vector<char> attr_data(attr_size + 1, 0);

                // get the attribute
                ret = H5LTget_attribute_string(group_id_data, name.c_str(), attr_name.c_str(), &attr_data[0]);
                if (ret < 0)
                    field_names.push_back("(UNLABELED)");
                else
                    field_names.push_back(std::string(&attr_data[0]));
            }
        }

        ~Implementation()
        {
            // truncate the data set to its actual size
            if (truncate_on_destroy && (! file_imp->read_only))
            {
                hsize_t dimensions[2] = { tuples, tuple_size };

                herr_t ret = H5Dset_extent(set_id, dimensions);
                if (0 > ret)
                    throw ScanFileHDF5Error("H5Sset_extent_simple", ret);
            }

            // create attributes for field info
            if (! file_imp->read_only)
            {
                unsigned i = 0;
                for (auto f = field_names.cbegin(), f_end = field_names.cend() ; f != f_end ; ++f, ++i)
                {
                    std::string attr_name = "FIELD_" + stringify(i) + "_NAME";
                    herr_t ret = H5LTset_attribute_string(group_id_data, name.c_str(), attr_name.c_str(), f->c_str());
                    if (ret < 0)
                        throw ScanFileHDF5Error("H5LTset_attribute_string", ret);
                }
            }

            herr_t ret = H5Fflush(set_id, H5F_SCOPE_LOCAL);
            if (0 > ret)
                throw ScanFileHDF5Error("H5Fflush", ret);

            ret = H5Dclose(set_id);
            if (0 > ret)
                throw ScanFileHDF5Error("H5Dclose", ret);

            ret = H5Sclose(space_id_memory);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sclose", ret);

            ret = H5Sclose(space_id_file_writing);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sclose", ret);

            ret = H5Sclose(space_id_file);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sclose", ret);
        }

        void append(const std::vector<double> & data)
        {
            unsigned buffer_length = data.size() / tuple_size;
            if (buffer_length > chunk_size)
                throw InternalError("Flushing a WriteBuffer with capacity > " + stringify(chunk_size) + " is not implemented yet");

            // can we store the data? if not, extend the data set
            if (capacity <= current_index + buffer_length)
            {
                capacity += chunk_size;

                hsize_t dimensions[2] = { capacity, tuple_size };
                hsize_t max_dimensions[2] = { H5S_UNLIMITED, tuple_size };

                herr_t ret = H5Sset_extent_simple(space_id_file_writing, 2, dimensions, max_dimensions);
                if (0 > ret)
                    throw ScanFileHDF5Error("H5Sset_extent_simple", ret);

                ret = H5Sset_extent_simple(space_id_file, 2, dimensions, max_dimensions);
                if (0 > ret)
                    throw ScanFileHDF5Error("H5Sset_extent_simple", ret);

                ret = H5Dset_extent(set_id, dimensions);
                if (0 > ret)
                    throw ScanFileHDF5Error("H5Dset_extent", ret);
            }

            // actually store the data
            hsize_t offset[2] = { current_index, 0 };
            hsize_t count[2] = { buffer_length, tuple_size };

            herr_t ret = H5Sselect_hyperslab(space_id_file_writing, H5S_SELECT_SET, offset, 0, count, 0);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sselect_hyperslab", ret);

            hsize_t max_dimensions[2] = { chunk_size, tuple_size };
            ret = H5Sset_extent_simple(space_id_memory, 2, count, max_dimensions);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sset_extent_simple", ret);

            ret = H5Dwrite(set_id, H5T_IEEE_F64LE, space_id_memory, space_id_file_writing, H5P_DEFAULT, &data[0]);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Dwrite", ret);

            current_index += buffer_length;
            tuples += buffer_length;

        }
    };

    ScanFile::DataSet::DataSet(Implementation<ScanFile> * imp, const std::string & name) :
        PrivateImplementationPattern<ScanFile::DataSet>(new Implementation<ScanFile::DataSet>(imp, name))
    {
    }

    ScanFile::DataSet::DataSet(Implementation<ScanFile> * imp, const std::string & name, unsigned tuple_size) :
        PrivateImplementationPattern<ScanFile::DataSet>(new Implementation<ScanFile::DataSet>(imp, name, tuple_size))
    {
    }

    ScanFile::DataSet::~DataSet()
    {
    }

    std::string
    ScanFile::DataSet::name() const
    {
        return _imp->name;
    }

    unsigned
    ScanFile::DataSet::tuple_size() const
    {
        return _imp->tuple_size;
    }

    unsigned
    ScanFile::DataSet::tuples() const
    {
        return _imp->tuples;
    }

    ScanFile::DataSet::FieldIterator
    ScanFile::DataSet::begin_fields()
    {
        return FieldIterator(_imp->field_names.begin());
    }

    ScanFile::DataSet::FieldIterator
    ScanFile::DataSet::end_fields()
    {
        return FieldIterator(_imp->field_names.end());
    }

    unsigned
    ScanFile::DataSet::find_field_index(const std::string & name) const
    {
        auto f = std::find(_imp->field_names.cbegin(), _imp->field_names.cend(), name);

        if (_imp->field_names.cend() == f)
            throw ScanFileError("Dataset '" + _imp->name + "' does not have a field named '" + name + "'");

        return std::distance(_imp->field_names.cbegin(), f);
    }

    ScanFile::Tuple
    ScanFile::DataSet::operator[] (const unsigned & index)
    {
        return ScanFile::Tuple(_imp, index);
    }

    ScanFile::DataSet &
    operator<< (ScanFile::DataSet & lhs, const std::vector<double> & rhs)
    {
        if (lhs._imp->tuple_size != rhs.size())
            throw InternalError("Trying to write a tuple of size '" + stringify(rhs.size()) + "' to a ScanFile::DataSet of width '" + stringify(lhs._imp->tuple_size) + "'");

        lhs._imp->append(rhs);

        return lhs;
    }

    ScanFile::DataSet &
    operator<< (ScanFile::DataSet & lhs, const ScanFile::WriteBuffer & rhs)
    {
        if (lhs._imp->tuple_size != rhs._imp->tuple_size)
            throw InternalError("Trying to write a buffer of size '" + stringify(rhs._imp->tuple_size) + "' to a ScanFile::DataSet of width '" + stringify(lhs._imp->tuple_size) + "'");

        lhs._imp->append(rhs._imp->buffer);

        return lhs;
    }

    template <>
    struct Implementation<ScanFile::Tuple>
    {
        std::shared_ptr<Implementation<ScanFile::DataSet>> data_set_imp;

        const unsigned tuple_size;

        unsigned tuples;

        unsigned index;

        double * buffer;

        hid_t set_id;

        hid_t space_id_memory;

        hid_t space_id_file_reading;

        Implementation(const std::shared_ptr<Implementation<ScanFile::DataSet>> & data_set_imp, const unsigned & index) :
            data_set_imp(data_set_imp),
            tuple_size(data_set_imp->tuple_size),
            tuples(data_set_imp->tuples),
            index(index),
            set_id(data_set_imp->set_id),
            space_id_memory(H5I_INVALID_HID),
            space_id_file_reading(H5Scopy(data_set_imp->space_id_file))
        {
            if (H5I_INVALID_HID == space_id_file_reading)
                throw ScanFileHDF5Error("H5Scopy", space_id_file_reading);

            hsize_t dimensions[2] = { 1, tuple_size };
            space_id_memory = H5Screate_simple(2, dimensions, 0);
            if (H5I_INVALID_HID == space_id_memory)
                throw ScanFileHDF5Error("H5Screate_simple", space_id_memory);

            buffer = new double[tuple_size];

            read();
        }

        ~Implementation()
        {
            delete[] buffer;

            herr_t ret = H5Sclose(space_id_memory);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sclose", ret);

            ret = H5Sclose(space_id_file_reading);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sclose", ret);
        }

        void read()
        {
            if (index >= tuples)
                return;

            hsize_t offset[2] = { index, 0 };
            hsize_t count[2] = { 1, tuple_size };

            herr_t ret = H5Sselect_hyperslab(space_id_file_reading, H5S_SELECT_SET, offset, 0, count, 0);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sselect_hyperslab", ret);

            ret = H5Dread(set_id, H5T_IEEE_F64LE, space_id_memory, space_id_file_reading, H5P_DEFAULT, buffer);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Dread", ret);
        }
    };

    ScanFile::Tuple::Tuple(const std::shared_ptr<Implementation<ScanFile::DataSet>> & data_set_imp, const unsigned & index) :
        PrivateImplementationPattern<ScanFile::Tuple>(new Implementation<ScanFile::Tuple>(data_set_imp, index))
    {
    }

    ScanFile::Tuple::~Tuple()
    {
    }

    ScanFile::Tuple &
    ScanFile::Tuple::operator++ ()
    {
        _imp->index += 1;
        _imp->read();

        return *this;
    }

    double
    ScanFile::Tuple::operator[] (const unsigned & index) const
    {
        return _imp->buffer[index];
    }
}
