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

#ifndef EOS_GUARD_SRC_UTILS_HDF5_HH
#define EOS_GUARD_SRC_UTILS_HDF5_HH 1

#include <eos/utils/exception.hh>
#include <eos/utils/hdf5-fwd.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/wrapped_forward_iterator.hh>
#include <hdf5.h>

#include <cstring>
#include <vector>

namespace eos
{
    /*!
     * HDF5Error is parent to all exceptions thrown by entities in namespace eos::hdf5, e.g when creating, opening or accessing hdf5 file data.
     */
    class HDF5Error :
        public Exception
    {
        public:
            /*!
             * Constructor.
             *
             * @param message The error message.
             */
            HDF5Error(const std::string & message);
    };

    namespace hdf5
    {
        /* Type System */
        typedef int hid_t;

        template <typename T_> struct DataType
        {
            static hid_t type_id();
        };

        template <> struct DataType<double>
        {
            static hid_t type_id() { return H5T_IEEE_F64LE; }
        };

        template <> struct DataType<unsigned>
        {
            static hid_t type_id() { return H5T_STD_U32LE; }
        };

        template <> struct DataType<int>
        {
            static hid_t type_id() { return H5T_STD_I32LE; }
        };

        template <> struct DataType<short>
        {
            static hid_t type_id() { return H5T_STD_I8LE; }
        };

        template <> struct DataType<const char *>
        {
            static hid_t type_id()
            {
                hid_t data_type = H5Tcopy(H5T_C_S1);
                H5Tset_size(data_type, H5T_VARIABLE);

                return data_type;
            }
        };

        /*
         * Type presents the C++ interface to HDF5 types handlers.
         */
        struct Type
        {
            virtual ~Type() = 0;

            virtual hid_t type_id() const = 0;

            virtual const std::string & name() const = 0;

            virtual hsize_t size() const = 0;

            virtual void copy_to_hdf5(const void * src, void * dest) const = 0;

            virtual void copy_from_hdf5(const void * src, void * dest) const = 0;
        };

        typedef std::shared_ptr<Type> TypePtr;

        /*!
         * Scalar represents a scalar built-in HDF5 data type.
         */
        template <typename T_> class Scalar :
            public Type
        {
            private:
                hid_t _type_id;

                std::string _name;

            public:
                typedef T_ Type;

                Scalar(const std::string & name) :
                    _type_id(hdf5::DataType<T_>::type_id()),
                    _name(name)
                {
                }

                ~Scalar()
                {
                }

                virtual hid_t type_id() const
                {
                    return _type_id;
                }

                virtual const std::string & name() const
                {
                    return _name;
                }

                virtual hsize_t size() const
                {
                    return H5Tget_size(_type_id);
                }

                virtual void copy_to_hdf5(const void * src, void * dest) const
                {
                    ::memcpy(dest, src, sizeof(T_));
                }

                virtual void copy_from_hdf5(const void * src, void * dest) const
                {
                    ::memcpy(dest, src, sizeof(T_));
                }
        };

        template <unsigned rank_, typename T_> class Array :
            public Type
        {
            private:
                hid_t _type_id;

                std::string _name;

                unsigned _elements;

            public:
                typedef std::vector<T_> Type;

                Array(const std::string & name, const std::initializer_list<hsize_t> & dimensions) :
                    _name(name)
                {
                    static_assert(rank_ > 0, "Array<0> is not supported");

                    if (dimensions.size() != rank_)
                        throw std::string("Array<> dimensions have wrong rank");

                    _type_id = H5Tarray_create2(hdf5::DataType<T_>::type_id(), rank_, &(*dimensions.begin()));

                    auto d = dimensions.begin();
                    _elements = *d; ++d;

                    for (auto d_end = dimensions.end() ; d != d_end ; ++d)
                    {
                        _elements *= *d;
                    }
                }

                virtual const std::string & name() const
                {
                    return _name;
                }

                virtual hid_t type_id() const
                {
                    return _type_id;
                }

                virtual hsize_t size() const
                {
                    return H5Tget_size(_type_id);
                }

                virtual void copy_to_hdf5(const void * src, void * dest) const
                {
                    const std::vector<T_> * _src = reinterpret_cast<const std::vector<T_> *>(src);
                    ::memcpy(dest, &(*_src)[0], _elements * sizeof(T_));
                }

                virtual void copy_from_hdf5(const void * src, void * dest) const
                {
                    std::vector<T_> * _dest = reinterpret_cast<std::vector<T_> *>(dest);
                    _dest->reserve(_elements);

                    ::memcpy(&(*_dest)[0], src, _elements * sizeof(T_));
                }
        };

        template <size_t index_, typename TypeTuple_, typename Tuple_> struct CompositeCopier
        {
            static std::ptrdiff_t copy_to_hdf5(const TypeTuple_ & types, const Tuple_ & src, char * dest)
            {
                std::ptrdiff_t offset = CompositeCopier<index_ - 1, TypeTuple_, Tuple_>::copy_to_hdf5(types, src, dest);

                std::get<index_>(types).copy_to_hdf5(&std::get<index_>(src), dest + offset);

                return offset + std::get<index_>(types).size();
            }

            static std::ptrdiff_t copy_from_hdf5(const TypeTuple_ & types, const char * src, Tuple_ & dest)
            {
                std::ptrdiff_t offset = CompositeCopier<index_ - 1, TypeTuple_, Tuple_>::copy_from_hdf5(types, src, dest);

                std::get<index_>(types).copy_from_hdf5(src + offset, &std::get<index_>(dest));

                return offset + std::get<index_>(types).size();
            }
        };

        template <typename TypeTuple_, typename Tuple_> struct CompositeCopier<0, TypeTuple_, Tuple_>
        {
            static std::ptrdiff_t copy_to_hdf5(const TypeTuple_ & types, const Tuple_ & src, char * dest)
            {
                std::get<0>(types).copy_to_hdf5(&std::get<0>(src), dest);

                return std::get<0>(types).size();
            }

            static std::ptrdiff_t copy_from_hdf5(const TypeTuple_ & types, const char * src, Tuple_ & dest)
            {
                std::get<0>(types).copy_from_hdf5(src, &std::get<0>(dest));

                return std::get<0>(types).size();
            }
        };

        template <typename ... T_> class Composite :
            public Type
        {
            private:
                hid_t _type_id;

                size_t _size;

                std::string _name;

                std::tuple<T_ ...> _types;

                template <typename R_> hid_t _insert(const hid_t & compound_id, size_t offset, const R_ & r)
                {
                    H5Tinsert(compound_id, r.name().c_str(), offset, r.type_id());

                    return compound_id;
                }

                template <typename R_, typename ... S_> hid_t _insert(const hid_t & compound_id, size_t offset, const R_ & r, const S_ & ... s)
                {
                    H5Tinsert(compound_id, r.name().c_str(), offset, r.type_id());

                    return _insert(compound_id, offset + r.size(), s ...);
                }

                template <typename R_> static size_t _compute_size(const R_ & r)
                {
                    return r.size();
                }

                template <typename R_, typename ... S_> static size_t _compute_size(const R_ & r, const S_ & ... s)
                {
                    return r.size() + _compute_size(s ...);
                }

            public:
                typedef std::tuple<typename T_::Type ...> Type;

                Composite(const std::string & name, const T_ & ... t) :
                    _type_id(_insert(H5Tcreate(H5T_COMPOUND, _compute_size(t ...)), 0, t ...)),
                    _size(H5Tget_size(_type_id)),
                    _name(name),
                    _types(t ...)
                {
                    if (_size != _compute_size(t ...))
                        throw HDF5Error("Composite::_compute_size returned different size than H5Tget_size");
                }

                ~Composite()
                {
                }

                virtual hid_t type_id() const
                {
                    return _type_id;
                }

                virtual const std::string & name() const
                {
                    return _name;
                }

                virtual hsize_t size() const
                {
                    return _size;
                }

                virtual void copy_to_hdf5(const void * src, void * dest) const
                {
                    const Type * _src = reinterpret_cast<const Type *>(src);
                    char * _dest = reinterpret_cast<char *>(dest);

                    CompositeCopier<sizeof...(T_) - 1, std::tuple<T_ ...>, Type>::copy_to_hdf5(_types, *_src, _dest);
                }

                virtual void copy_from_hdf5(const void * src, void * dest) const
                {
                    const char * _src = reinterpret_cast<const char *>(src);
                    Type * _dest = reinterpret_cast<Type *>(dest);

                    CompositeCopier<sizeof...(T_) - 1, std::tuple<T_ ...>, Type>::copy_from_hdf5(_types, _src, *_dest);
                }
        };

        /* Handle Classes */

        class FileHandle :
            public PrivateImplementationPattern<hdf5::FileHandle>
        {
            private:
                FileHandle(const hid_t & file_id, bool read_only);

            public:
                friend class File;

                ~FileHandle();

                hid_t id() const;

                bool read_only() const;
        };

        class DataSetHandle :
            public PrivateImplementationPattern<hdf5::DataSetHandle>
        {
            private:
                DataSetHandle(const FileHandle & file_handle, const hid_t & data_set_id, const hid_t & space_id_file, hsize_t size);

            public:
                friend class File;

                ~DataSetHandle();

                hid_t set_id() const;

                hid_t space_id() const;

                hid_t type_id() const;

                hsize_t size() const;

                void select(hsize_t start, hsize_t count);

                void write_one(const void * buffer);

                void read_one(void * buffer);

                AttributeHandle create_attribute(const std::string & name, const hid_t & type_id);

                AttributeHandle open_attribute(const std::string & name, const hid_t & type_id);
        };

        class AttributeHandle :
            PrivateImplementationPattern<AttributeHandle>
        {
            private:
                AttributeHandle(const DataSetHandle & data_set_handle, const hid_t & attribute_id);

            public:
                friend class DataSetHandle;

                ~AttributeHandle();

                hid_t attribute_id() const;

                hid_t type_id() const;

                void write(const void * buffer);

                void read(void * buffer) const;
        };

        class File
        {
            private:
                FileHandle _handle;

                std::string _name;

                /// Constructor.
                File(const FileHandle & handle, const std::string & name);

                DataSetHandle _open_data_set(const std::string & name, const TypePtr & type) const;

                DataSetHandle _create_data_set(const std::string & name, const TypePtr & type);

            public:
                ///@name Basic Functions
                ///@{
                /// Destructor.
                ~File();

                /*!
                 * Named constructor to create a new file.
                 *
                 * @param filename   File which shall be created.
                 * @param creator    Name of the creating program.
                 */
                static File Create(const std::string & filename);

                /*!
                 * Named constructor to open an existing file.
                 *
                 * @param file_name File which shall be opened.
                 */
                static File Open(const std::string & file_name, unsigned mode = H5F_ACC_RDONLY);

                static bool Exists(const std::string & filename);

                ///Access to file name
                const std::string & name() const
                {
                    return _name;
                }
                ///@}

                ///@name Data Access
                ///@{
                /*!
                 * Create a new data set by name.
                 *
                 * @param name   Absolute name of the new data set.
                 * @param t Instance of any of Scalar, Array or Composite that represents this data set's underlying data type.
                 */
                template <typename T_> DataSet<T_> create_data_set(const std::string & name, const T_ & t)
                {
                    TypePtr type = TypePtr(new T_(t));
                    auto result = DataSet<T_>(_create_data_set(name, type), type);

                    return result;
                }

                /*!
                 * Open an existing data set by name.
                 *
                 * @param name Absolute name of the data set.
                 * @param t Instance of any of Scalar, Array or Composite that represents this data set's underlying data type.
                 */
                template <typename T_> DataSet<T_> open_data_set(const std::string & name, const T_ & t)
                {
                    TypePtr type = TypePtr(new T_(t));
                    auto result = DataSet<T_>(_open_data_set(name, type), type);

                    return result;
                }

                /*!
                 * Open a data set by name. Create the data set if it doesn't exist.
                 *
                 * @note The error handling of the HDF5 library is suppressed for the duration of the execution.
                 * Do not use this method when parallel access to a file is expected.
                 *
                 * @param name Absolute name of the data set.
                 * @param t Instance of any of Scalar, Array or Composite that represents this data set's underlying data type.
                 */
                template <typename T_> DataSet<T_> create_or_open_data_set(const std::string & name, const T_ & t)
                {
                    TypePtr type = TypePtr(new T_(t));
                    H5E_BEGIN_TRY
                    {
                        // check if data set exists already, we suppress HDF5 error output
                        try
                        {
                            auto result = DataSet<T_>(_open_data_set(name, type), type);
                            return result;
                        }
                        catch (HDF5Error &)
                        {
                        }
                    }
                    H5E_END_TRY;
                    auto result = DataSet<T_>(_create_data_set(name, type), type);
                    return result;
                }

                ///@}

                ///@name Group Operations
                ///@{
                /*!
                 * Copy all objects under source directory in this file to the destination file.
                 * If the destination directory is omitted, the source directory is recreated
                 * in the destination file.
                 */
                void copy(const std::string& source_directory, File& destination_file, const std::string& destination_directory="") const;

                /*!
                 * Create a group within this HDF5 file.
                 *
                 * @param name Name of the group that shall be created.
                 */
                void create_group(const std::string & name);

                bool group_exists(const std::string & name);

                /// List how many objects, i.e. groups or data sets, are in a subdirectory.
                unsigned number_of_objects(const std::string & name);
                ///@}
        };

        /*!
         * DataSet<> represents one of the data sets within an hdf5 file.
         */
        template <typename T_> class DataSet
        {
            public:
                typedef typename T_::Type RecordType;

            private:
                DataSetHandle _handle;

                TypePtr _type;

                hsize_t _index;

                std::vector<char> _buffer;

                DataSet(const DataSetHandle & handle, const TypePtr & type) :
                    _handle(handle),
                    _type(type),
                    _index(0),
                    _buffer(type->size(), '\0')
                {
                }

                /// Append only to the end of the data set.
                void _append(const RecordType & record)
                {
                    _type->copy_to_hdf5(&record, &_buffer[0]);

                    _handle.select(_handle.size(), 1);
                    _handle.write_one(&_buffer[0]);
                }

                void _extract(RecordType & record)
                {
                    _handle.select(_index, 1);
                    _handle.read_one(&_buffer[0]);
                    ++_index;

                    _type->copy_from_hdf5(&_buffer[0], &record);
                }

            public:
                friend class hdf5::File;
                friend DataSet<T_> & operator<< <> (DataSet<T_> &, const typename T_::Type &);
                friend DataSet<T_> & operator>> <> (DataSet<T_> &, typename T_::Type &);

                ///@name Basic Functions
                ///@{
                /// Destructor.
                ~DataSet()
                {
                }
                ///@}

                ///@name Metadata Access
                ///@{
                /// Name of the data set
                std::string name() const;

                /// Number of rows in data set
                unsigned records() const
                {
                    return _handle.size();
                }

                ///@name Data Access
                ///@{

                /// Set index to last record.
                void end()
                {
                    _index = _handle.size() - 1;
                }

                /*!
                 * Retrieve a record.
                 *
                 * @param index Index of the record that shall be retrieved.
                 */
                RecordType operator[] (const unsigned & index);

                void set_index(const unsigned & index)
                {
                    _index = index;
                }

                ///@}

                ///@name Attribute Access
                ///@{
                template <typename U_> Attribute<U_> create_attribute(const std::string & name, const U_ & u)
                {
                    TypePtr type(new U_(u));
                    AttributeHandle handle = _handle.create_attribute(name, type->type_id());

                    return Attribute<U_>(handle, type);
                }

                template <typename U_> Attribute<U_> create_or_open_attribute(const std::string & name, const U_ & u)
                {
                    TypePtr type(new U_(u));

                    H5E_BEGIN_TRY
                    {
                        // check if data set exists already, we suppress HDF5 error output
                        try
                        {
                            AttributeHandle handle = _handle.open_attribute(name, type->type_id());
                            return Attribute<U_>(handle, type);
                        }
                        catch (HDF5Error &)
                        {
                        }
                    }
                    H5E_END_TRY;

                    AttributeHandle handle = _handle.create_attribute(name, type->type_id());
                    return Attribute<U_>(handle, type);
                }

                template <typename U_> Attribute<U_> open_attribute(const std::string & name, const U_ & u)
                {
                    TypePtr type(new U_(u));
                    AttributeHandle handle = _handle.open_attribute(name, type->type_id());

                    return Attribute<U_>(handle, type);
                }
                ///@}
        };

        /*!
         * Append a record to the end of an HDF5 DataSet<> .
         *
         * @param lhs The DataSet<> to which shall be written.
         * @param rhs The record which shall be appended.
         */
        template <typename T_> DataSet<T_> & operator<< (DataSet<T_> & lhs, const typename T_::Type & rhs)
        {
            lhs._append(rhs);

            return lhs;
        }

        /*!
         * Extract a record from an HDF5 DataSet.
         *
         * @param lhs The DataSet<> from which shall be read.
         * @param rhs The record to which the data shall be copied.
         */
        template <typename T_> DataSet<T_> & operator>> (DataSet<T_> & lhs, typename T_::Type & rhs)
        {
            lhs._extract(rhs);

            return lhs;
        }

        template <typename T_> class Attribute
        {
            public:
                typedef typename T_::Type RecordType;

            private:
                AttributeHandle _handle;

                TypePtr _type;

                mutable std::vector<char> _buffer;

                Attribute(const AttributeHandle & handle, TypePtr type) :
                    _handle(handle),
                    _type(type),
                    _buffer(type->size(), '\0')
                {
                }

            public:
                template <typename U_> friend class DataSet;

                RecordType value() const
                {
                    RecordType result;

                    _handle.read(&_buffer[0]);
                    _type->copy_from_hdf5(&_buffer[0], &result);

                    return result;
                }

                Attribute & operator= (const RecordType & rhs)
                {
                    _type->copy_to_hdf5(&rhs, &_buffer[0]);
                    _handle.write(&_buffer[0]);

                    return *this;
                }
        };
    }
}

#endif
