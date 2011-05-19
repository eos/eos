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
#include <limits>

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

    struct HDF5File
    {
        std::string file_name;

        const hid_t file_id;

        hid_t group_id_metadata;

        hid_t group_id_data;

        std::string creator;

        std::string eos_version;

        bool read_only;

        // Create a new HDF5 file
        HDF5File(const std::string & file_name, const std::string & creator) :
            file_name(file_name),
            file_id(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)),
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

        // Open an existing HDF5 file
        HDF5File(const std::string & file_name) :
            file_name(file_name),
            file_id(H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)),
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
        }

        ~HDF5File()
        {
            herr_t ret;

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
    };

    template class WrappedForwardIterator<ScanFile::IteratorTag, ScanFile::DataSet>;

    template <>
    struct Implementation<ScanFile>
    {
        std::shared_ptr<HDF5File> hdf5_file;

        std::vector<ScanFile::DataSet> data_sets;

        // Create a new HDF5 file
        Implementation(const std::string & file_name, const std::string & creator) :
            hdf5_file(new HDF5File(file_name, creator))
        {
        }

        // Callback function for iteration over data sets
        static herr_t add_data_set(hid_t /* unused */, const char * name, const H5L_info_t * /* unused */, void * _imp)
        {
            Implementation<ScanFile> * imp = static_cast<Implementation<ScanFile> *>(_imp);

            imp->data_sets.push_back(ScanFile::DataSet(imp->hdf5_file, std::string(name)));

            return 0;
        }

        // Open an existing HDF5 file
        Implementation(const std::string & file_name) :
            hdf5_file(new HDF5File(file_name))
        {
            hsize_t idx = 0;
            herr_t ret = H5Literate(hdf5_file->group_id_data, H5_INDEX_NAME, H5_ITER_INC, &idx, &Implementation<ScanFile>::add_data_set, this);
            if (0 > ret)
                throw ScanFileHDF5Error("H5Literate", ret);
        }

        ~Implementation()
        {
            // first close our handles to the data_sets
            data_sets.clear();

            // then close our handle to the file
            hdf5_file.reset();
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
    ScanFile::Create(const std::string & file_name, const std::string & creator)
    {
        H5dont_atexit();
        return ScanFile(new Implementation<ScanFile>(file_name, creator));
    }

    ScanFile
    ScanFile::Open(const std::string & file_name)
    {
        H5dont_atexit();
        return ScanFile(new Implementation<ScanFile>(file_name));
    }

    const std::string &
    ScanFile::creator() const
    {
        return _imp->hdf5_file->creator;
    }

    const std::string &
    ScanFile::eos_version() const
    {
        return _imp->hdf5_file->eos_version;
    }

    const std::string &
    ScanFile::file_name() const
    {
        return _imp->hdf5_file->file_name;
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
    ScanFile::add(const std::string & name, unsigned fields)
    {
        ScanFile::DataSet result = ScanFile::DataSet(_imp->hdf5_file, name, fields);
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
        static const unsigned chunk_size;

        unsigned capacity;

        unsigned records;

        unsigned fields;

        std::vector<double> buffer;

        Implementation(const unsigned & fields, const unsigned & capacity) :
            capacity(capacity),
            records(0),
            fields(fields)
        {
            buffer.reserve(fields * capacity);
        }

        void clear()
        {
            records = 0;
            buffer.clear(); // clear only affects size, not capacity
        }

        void reserve(const unsigned & capacity)
        {
            buffer.reserve(fields * capacity);
            this->capacity = capacity;
        }
    };

    const unsigned Implementation<ScanFile::WriteBuffer>::chunk_size = 1024;

    ScanFile::WriteBuffer::WriteBuffer(const unsigned & fields, const unsigned & capacity) :
        PrivateImplementationPattern<ScanFile::WriteBuffer>(new Implementation<ScanFile::WriteBuffer>(fields, capacity))
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
        return _imp->capacity;
    }

    void
    ScanFile::WriteBuffer::reserve(const unsigned & capacity)
    {
        _imp->reserve(capacity);
    }

    unsigned
    ScanFile::WriteBuffer::size() const
    {
        return _imp->records;
    }

    ScanFile::WriteBuffer &
    operator<< (ScanFile::WriteBuffer & lhs, const std::vector<double> & rhs)
    {
        if (lhs._imp->records == lhs._imp->capacity)
            lhs._imp->reserve(lhs._imp->capacity + Implementation<ScanFile::WriteBuffer>::chunk_size);

        for (auto i = rhs.cbegin(), i_end = rhs.cend() ; i != i_end ; ++i)
        {
            lhs._imp->buffer.push_back(*i);
        }

        ++lhs._imp->records;

        return lhs;
    }

    /* ScanFile::FieldInfo */
    ScanFile::FieldInfo::FieldInfo(const std::string & name, const double & min, const double & max, bool nuisance, bool discrete) :
        name(name),
        min(min),
        max(max),
        nuisance(nuisance),
        discrete(discrete)
    {
    }

    /* ScanFile::DataSet */

    template class WrappedForwardIterator<ScanFile::DataSet::FieldIteratorTag, ScanFile::FieldInfo>;

    template <>
    struct Implementation<ScanFile::DataSet>
    {
        static const unsigned chunk_size;

        std::shared_ptr<HDF5File> file;

        std::string name;

        unsigned fields;

        unsigned current_index;

        unsigned records;

        unsigned capacity;

        hid_t group_id_data;

        hid_t set_id;

        hid_t dcpl_id;

        hid_t space_id_memory;

        hid_t space_id_file;

        hid_t space_id_file_writing;

        std::vector<ScanFile::FieldInfo> field_infos;

        bool truncate_on_destroy;

        // Constructor to create a new data set
        Implementation(const std::shared_ptr<HDF5File> & file, const std::string & name, unsigned fields) :
            file(file),
            name(name),
            fields(fields),
            current_index(0),
            records(0),
            capacity(chunk_size),
            group_id_data(file->group_id_data),
            set_id(H5I_INVALID_HID),
            space_id_memory(H5I_INVALID_HID),
            space_id_file(H5I_INVALID_HID),
            space_id_file_writing(H5I_INVALID_HID),
            field_infos(fields, ScanFile::FieldInfo("", 0.0, 0.0, false)),
            truncate_on_destroy(true)
        {
            // create space id for in-memory representation
            {
                hsize_t dimensions[2] = { 1, fields };
                hsize_t max_dimensions[2] = { capacity, fields };
                space_id_memory = H5Screate_simple(2, dimensions, max_dimensions);
                if (H5I_INVALID_HID == space_id_memory)
                    throw ScanFileHDF5Error("H5Screate_simple", space_id_memory);
            }

            // create space id for in-file representation
            {
                hsize_t dimensions[2] = { capacity, fields };
                hsize_t max_dimensions[2] = { H5S_UNLIMITED, fields };
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

                hsize_t chunk_size[2] = { capacity, fields };
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
        Implementation(const std::shared_ptr<HDF5File> & file, const std::string & name) :
            file(file),
            name(name),
            fields(0),
            current_index(0),
            group_id_data(file->group_id_data),
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

                current_index = capacity = records = dimensions[0];
                fields = dimensions[1];
            }

            // create space id for in-memory representation
            {
                hsize_t dimensions[2] = { 1, fields };
                hsize_t max_dimensions[2] = { H5S_UNLIMITED, fields };
                space_id_memory = H5Screate_simple(2, dimensions, max_dimensions);
                if (H5I_INVALID_HID == space_id_memory)
                    throw ScanFileHDF5Error("H5Screate_simple", space_id_memory);
            }

            // read field infos from attributes
            for (unsigned i = 0 ; i < fields ; ++i)
            {
                // name
                std::string name_attr_name = "FIELD_" + stringify(i) + "_NAME";
                hsize_t name_attr_rank;
                H5T_class_t name_attr_class;
                size_t name_attr_size;

                // how large is the attribute's value?
                herr_t ret = H5LTget_attribute_info(group_id_data, name.c_str(), name_attr_name.c_str(), &name_attr_rank, &name_attr_class, &name_attr_size);
                if (ret < 0)
                    throw ScanFileHDF5Error("H5LTget_attribute_info", ret);

                // allocate a suitable buffer
                // fuck the HDF5 Lite API!
                std::vector<char> name_attr_data(name_attr_size + 1, '\0');
                std::string name_attr;

                // get the name
                ret = H5LTget_attribute_string(group_id_data, name.c_str(), name_attr_name.c_str(), &name_attr_data[0]);
                if (ret < 0)
                    name_attr = "(UNLABELED)";
                else
                    name_attr = std::string(&name_attr_data[0]);

                // min
                std::string min_attr_name = "FIELD_" + stringify(i) + "_MIN";
                double min_attr;
                ret = H5LTget_attribute_double(group_id_data, name.c_str(), min_attr_name.c_str(), &min_attr);
                if (ret < 0)
                    min_attr = std::numeric_limits<double>::min();

                // max
                std::string max_attr_name = "FIELD_" + stringify(i) + "_MAX";
                double max_attr;
                ret = H5LTget_attribute_double(group_id_data, name.c_str(), max_attr_name.c_str(), &max_attr);
                if (ret < 0)
                    max_attr = std::numeric_limits<double>::max();

                // nuisance
                std::string nuisance_attr_name = "FIELD_" + stringify(i) + "_NUISANCE";
                unsigned short nuisance_attr;
                ret = H5LTget_attribute_ushort(group_id_data, name.c_str(), nuisance_attr_name.c_str(), &nuisance_attr);
                if (ret < 0)
                    nuisance_attr = false;

                // discrete
                std::string discrete_attr_name = "FIELD_" + stringify(i) + "_DISCRETE";
                unsigned short discrete_attr;
                ret = H5LTget_attribute_ushort(group_id_data, name.c_str(), discrete_attr_name.c_str(), &discrete_attr);
                if (ret < 0)
                    discrete_attr = false;

                field_infos.push_back(ScanFile::FieldInfo(name_attr, min_attr, max_attr, nuisance_attr, discrete_attr));
            }
        }

        ~Implementation()
        {
            // truncate the data set to its actual size
            if (truncate_on_destroy && (! file->read_only))
            {
                hsize_t dimensions[2] = { records, fields };

                herr_t ret = H5Dset_extent(set_id, dimensions);
                if (0 > ret)
                    throw ScanFileHDF5Error("H5Dset_extent", ret);
            }

            // create attributes for field info
            if (! file->read_only)
            {
                unsigned i = 0;
                for (auto f = field_infos.cbegin(), f_end = field_infos.cend() ; f != f_end ; ++f, ++i)
                {
                    std::string name_attr_name = "FIELD_" + stringify(i) + "_NAME",
                        max_attr_name = "FIELD_" + stringify(i) + "_MAX",
                        min_attr_name = "FIELD_" + stringify(i) + "_MIN",
                        nuisance_attr_name = "FIELD_" + stringify(i) + "_NUISANCE",
                        discrete_attr_name = "FIELD_" + stringify(i) + "_DISCRETE";

                    herr_t ret = H5LTset_attribute_string(group_id_data, name.c_str(), name_attr_name.c_str(), f->name.c_str());
                    if (ret < 0)
                        throw ScanFileHDF5Error("H5LTset_attribute_string", ret);

                    ret = H5LTset_attribute_double(group_id_data, name.c_str(), max_attr_name.c_str(), &f->max, 1);
                    if (ret < 0)
                        throw ScanFileHDF5Error("H5LTset_attribute_string", ret);

                    ret = H5LTset_attribute_double(group_id_data, name.c_str(), min_attr_name.c_str(), &f->min, 1);
                    if (ret < 0)
                        throw ScanFileHDF5Error("H5LTset_attribute_string", ret);

                    unsigned short nuisance = f->nuisance;
                    ret = H5LTset_attribute_ushort(group_id_data, name.c_str(), nuisance_attr_name.c_str(), &nuisance, 1);
                    if (ret < 0)
                        throw ScanFileHDF5Error("H5LTset_attribute_string", ret);

                    unsigned short discrete = f->discrete;
                    ret = H5LTset_attribute_ushort(group_id_data, name.c_str(), discrete_attr_name.c_str(), &discrete, 1);
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
                throw ScanFileHDF5Error("Implementation<ScanFile::DataSet>: H5Sclose(space_id_memory)", ret);

            ret = H5Sclose(space_id_file_writing);
            if (ret < 0)
                throw ScanFileHDF5Error("Implementation<ScanFile::DataSet>: H5Sclose(space_id_file_writing)", ret);

            ret = H5Sclose(space_id_file);
            if (ret < 0)
                throw ScanFileHDF5Error("Implementation<ScanFile::DataSet>: H5Sclose(space_id_file)", ret);
        }

        void append(const std::vector<double> & data)
        {
            unsigned buffer_length = data.size() / fields;
            if (0 == buffer_length)
                return;

            // can we store the data? if not, extend the data set
            if (capacity <= current_index + buffer_length)
            {
                unsigned need = buffer_length + current_index - capacity;
                unsigned chunks = need / chunk_size + 1;

                capacity += chunks * chunk_size;

                hsize_t dimensions[2] = { capacity, fields };
                hsize_t max_dimensions[2] = { H5S_UNLIMITED, fields };

                herr_t ret = H5Sset_extent_simple(space_id_file_writing, 2, dimensions, max_dimensions);
                if (0 > ret)
                    throw ScanFileHDF5Error("H5Sset_extent_simple(space_id_file_writing, 2, ["
                            + stringify(dimensions[0]) + "," + stringify(dimensions[1]) + "], ["
                            + stringify(max_dimensions[0]) + "," + stringify(max_dimensions[1]) + "])",
                            ret);

                ret = H5Sset_extent_simple(space_id_file, 2, dimensions, max_dimensions);
                if (0 > ret)
                    throw ScanFileHDF5Error("H5Sset_extent_simple(space_id_file, ...)", ret);

                ret = H5Dset_extent(set_id, dimensions);
                if (0 > ret)
                    throw ScanFileHDF5Error("H5Dset_extent", ret);
            }

            // actually store the data
            hsize_t offset[2] = { current_index, 0 };
            hsize_t count[2] = { buffer_length, fields };

            herr_t ret = H5Sselect_hyperslab(space_id_file_writing, H5S_SELECT_SET, offset, 0, count, 0);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sselect_hyperslab", ret);

            hsize_t max_dimensions[2] = { H5S_UNLIMITED, fields };
            ret = H5Sset_extent_simple(space_id_memory, 2, count, max_dimensions);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sset_extent_simple(space_id_memory, 2, [" + stringify(buffer_length) + "," + stringify(fields) + "], ["
                        + "UNLIMITED" + "," + stringify(fields) + ")", ret);

            ret = H5Dwrite(set_id, H5T_IEEE_F64LE, space_id_memory, space_id_file_writing, H5P_DEFAULT, &data[0]);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Dwrite", ret);

            current_index += buffer_length;
            records += buffer_length;

        }
    };

    const unsigned Implementation<ScanFile::DataSet>::chunk_size = 1024;

    ScanFile::DataSet::DataSet(const std::shared_ptr<HDF5File> & file, const std::string & name) :
        PrivateImplementationPattern<ScanFile::DataSet>(new Implementation<ScanFile::DataSet>(file, name))
    {
    }

    ScanFile::DataSet::DataSet(const std::shared_ptr<HDF5File> & file, const std::string & name, unsigned fields) :
        PrivateImplementationPattern<ScanFile::DataSet>(new Implementation<ScanFile::DataSet>(file, name, fields))
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
    ScanFile::DataSet::fields() const
    {
        return _imp->fields;
    }

    unsigned
    ScanFile::DataSet::records() const
    {
        return _imp->records;
    }

    ScanFile::DataSet::FieldIterator
    ScanFile::DataSet::begin_fields()
    {
        return FieldIterator(_imp->field_infos.begin());
    }

    ScanFile::DataSet::FieldIterator
    ScanFile::DataSet::end_fields()
    {
        return FieldIterator(_imp->field_infos.end());
    }

// TODO: Replace by lambda function!
namespace
{
    struct CompareFieldInfoByName
    {
        std::string _name;

        CompareFieldInfoByName(const std::string & name) :
            _name(name)
        {
        }

        bool operator() (const eos::ScanFile::FieldInfo & x)
        {
            return x.name == _name;
        }
    };
}

    unsigned
    ScanFile::DataSet::find_field_index(const std::string & name) const
    {
        auto f = std::find_if(_imp->field_infos.cbegin(), _imp->field_infos.cend(), CompareFieldInfoByName(name));

        if (_imp->field_infos.cend() == f)
            throw ScanFileError("Dataset '" + _imp->name + "' does not have a field named '" + name + "'");

        return std::distance(_imp->field_infos.cbegin(), f);
    }

    ScanFile::Record
    ScanFile::DataSet::operator[] (const unsigned & index)
    {
        return ScanFile::Record(_imp, index);
    }

    ScanFile::DataSet &
    operator<< (ScanFile::DataSet & lhs, const std::vector<double> & rhs)
    {
        if (lhs._imp->fields != rhs.size())
            throw InternalError("Trying to write a tuple of size '" + stringify(rhs.size()) + "' to a ScanFile::DataSet of width '" + stringify(lhs._imp->fields) + "'");

        lhs._imp->append(rhs);

        return lhs;
    }

    ScanFile::DataSet &
    operator<< (ScanFile::DataSet & lhs, const ScanFile::WriteBuffer & rhs)
    {
        if (lhs._imp->fields != rhs._imp->fields)
            throw InternalError("Trying to write a buffer of size '" + stringify(rhs._imp->fields) + "' to a ScanFile::DataSet of width '" + stringify(lhs._imp->fields) + "'");

        lhs._imp->append(rhs._imp->buffer);

        return lhs;
    }

    template <>
    struct Implementation<ScanFile::Record>
    {
        std::shared_ptr<Implementation<ScanFile::DataSet>> data_set_imp;

        const unsigned fields;

        unsigned records;

        unsigned index;

        std::vector<double> buffer;

        hid_t set_id;

        hid_t space_id_memory;

        hid_t space_id_file_reading;

        Implementation(const std::shared_ptr<Implementation<ScanFile::DataSet>> & data_set_imp, const unsigned & index) :
            data_set_imp(data_set_imp),
            fields(data_set_imp->fields),
            records(data_set_imp->records),
            index(index),
            buffer(fields, 0.0),
            set_id(data_set_imp->set_id),
            space_id_memory(H5I_INVALID_HID),
            space_id_file_reading(H5Scopy(data_set_imp->space_id_file))
        {
            if (H5I_INVALID_HID == space_id_file_reading)
                throw ScanFileHDF5Error("H5Scopy", space_id_file_reading);

            hsize_t dimensions[2] = { 1, fields };
            space_id_memory = H5Screate_simple(2, dimensions, 0);
            if (H5I_INVALID_HID == space_id_memory)
                throw ScanFileHDF5Error("H5Screate_simple", space_id_memory);

            read();
        }

        ~Implementation()
        {
            herr_t ret = H5Sclose(space_id_memory);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sclose", ret);

            ret = H5Sclose(space_id_file_reading);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sclose", ret);
        }

        void read()
        {
            if (index >= records)
                return;

            hsize_t offset[2] = { index, 0 };
            hsize_t count[2] = { 1, fields };

            herr_t ret = H5Sselect_hyperslab(space_id_file_reading, H5S_SELECT_SET, offset, 0, count, 0);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Sselect_hyperslab", ret);

            ret = H5Dread(set_id, H5T_IEEE_F64LE, space_id_memory, space_id_file_reading, H5P_DEFAULT, &buffer[0]);
            if (ret < 0)
                throw ScanFileHDF5Error("H5Dread", ret);
        }
    };

    ScanFile::Record::Record(const std::shared_ptr<Implementation<ScanFile::DataSet>> & data_set_imp, const unsigned & index) :
        PrivateImplementationPattern<ScanFile::Record>(new Implementation<ScanFile::Record>(data_set_imp, index))
    {
    }

    ScanFile::Record::~Record()
    {
    }

    const std::vector<double> &
    ScanFile::Record::data() const
    {
        return _imp->buffer;
    }

    ScanFile::Record &
    ScanFile::Record::operator++ ()
    {
        _imp->index += 1;
        _imp->read();

        return *this;
    }

    double
    ScanFile::Record::operator[] (const unsigned & index) const
    {
        return _imp->buffer[index];
    }

    ScanFile::WriteBuffer &
    operator<< (ScanFile::WriteBuffer & lhs, const ScanFile::Record & rhs)
    {
        if (Implementation<ScanFile::WriteBuffer>::chunk_size == lhs._imp->records)
            throw InternalError("Extending WriteBuffer capacity is not yet implemented");

        if (rhs._imp->fields != lhs._imp->fields)
            throw InternalError("Trying to write a ScanFile::Record of width '" + stringify(rhs._imp->fields) + "' to a ScanFile::WriteBuffer of width '" + stringify(lhs._imp->fields) + "'");

        lhs._imp->buffer.insert(lhs._imp->buffer.end(), rhs._imp->buffer.cbegin(), rhs._imp->buffer.cend());
        ++lhs._imp->records;

        return lhs;
    }
}
