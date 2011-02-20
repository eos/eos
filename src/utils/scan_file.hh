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

#ifndef EOS_GUARD_SRC_UTILS_SCAN_FILE_HH
#define EOS_GUARD_SRC_UTILS_SCAN_FILE_HH 1

#include <src/utils/exception.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace eos
{
    /*!
     * ScanFileError is parent to all exceptions thrown by ScanFile when creating, opening or accessing scan file data.
     */
    class ScanFileError :
        public Exception
    {
        public:
            /*!
             * Constructor.
             *
             * @param message The error message.
             */
            ScanFileError(const std::string & message);
    };

    /*!
     * ScanFileHDF5Error is thrown when an error occurs in interfacing libhdf5 during any
     * of ScanFile's operations.
     */
    class ScanFileHDF5Error :
        public ScanFileError
    {
        public:
            /*!
             * Constructor.
             *
             * @param function Name of the HDF5 function that failed.
             * @param code     Error code returned by the failing function.
             */
            ScanFileHDF5Error(const std::string & function, int code);
    };

    /*!
     * ScanFile represents an HDF5 formatted file that contains EOS scan data.
     *
     * An HDF5 DDL representation of a scan file follows:
     *
     * @verbatim
     *   HDF5 "<FILE>" {
     *   GROUP "/" {
     *     GROUP "metadata" {
     *       DATASET "creator" {
     *         DATATYPE H5T_STRING {
     *           STRSIZE <CREATORLEN>
     *           STRPAD H5T_STR_NULLTERM
     *           CSET H5T_CSET_ASCII
     *           CTYPE H5T_C_S1
     *         }
     *         DATASPACE SCALAR
     *         DATA {
     *             "<CREATOR>"
     *         }
     *       }
     *       DATASET "eos_version" {
     *         DATATYPE H5T_STRING {
     *           STRSIZE <EOSVERSIONLEN>
     *           STRPAD H5T_STR_NULLTERM
     *           CSET H5T_CSET_ASCII
     *           CTYPE H5T_C_S1
     *         }
     *         DATASPACE SCALAR
     *         DATA {
     *             "<EOSVERSION>"
     *         }
     *       }
     *     }
     *     DATASET "scan" {
     *       DATATYPE H5T_IEEE_F64LE
     *       DATASPACE SIMPLE { ( <TUPLES>, <ELEMENTS> ) / ( <TUPLES>, <ELEMENTS> ) }
     *       DATA {
     *         ...
     *       }
     *     }
     *   }
     *   }
     * @endverbatim
     *
     * Here, <FILE> is the filename, <TUPLES> is the number of scan tuples of <ELEMENTS> elements
     * each. <CREATOR> is an identifier of the creating program, and <EOSVERSION> is a string
     * representation of the EOS version/revision that was used to create the file.
     */
    class ScanFile :
        public PrivateImplementationPattern<ScanFile>
    {
        private:
            /// Constructor.
            ScanFile(Implementation<ScanFile> * imp);

        public:
            class Tuple;

            ///@name Basic Functions
            ///@{
            /// Destructor.
            ~ScanFile();

            /*!
             * Named constructor to create a new file.
             *
             * @param filename   File which shall be created.
             * @param creator    Name of the creating program.
             * @param tuple_size Number of elements of each scan tuple.
             * @param scan_size  Number of scan tuples in file.
             */
            static ScanFile Create(const std::string & filename, const std::string & creator,
                    const unsigned & tuple_size, const unsigned & scan_size);

            /*!
             * Named constructor to open an existing file.
             *
             * @param filename File which shall be opened.
             */
            static ScanFile Open(const std::string & filename);
            ///@}

            ///@name Metadata Access
            ///@{
            /// Retrieve the file creator's name.
            const std::string & creator() const;

            /// Retrieve the file creator's version information.
            const std::string & eos_version() const;

            /// Retrieve the number of elements of a scan tuple.
            int tuple_size() const;

            /// Retrieve the number of scan tuples in the file.
            int scan_size() const;
            ///@}

            ///@name Data Access
            ///@{
            /*!
             * Retrieve a scan tuple by its index.
             *
             * @param index Index of the scan tuple that shall be retrieved.
             */
            Tuple operator[] (const unsigned & index);
            ///@}
    };

    /*!
     * ScanFile::Tuple represents one of the scan tuples within a ScanFile.
     */
    class ScanFile::Tuple :
        public PrivateImplementationPattern<ScanFile::Tuple>
    {
        private:
            Tuple(const std::shared_ptr<Implementation<ScanFile>> & imp, const unsigned & offset);

        public:
            friend class ScanFile;

            ///@name Basic Functions
            ///@{
            /// Destructor.
            ~Tuple();
            ///@}

            ///@name Data Access
            ///@{
            /// Advance the scan tuple by one within the parent ScanFile.
            Tuple & operator++ ();

            /*!
             * Retrieve a tuple element.
             *
             * @param index Index of the tuple element that shall be retrieved.
             */
            double & operator[] (const unsigned & index);

            /// Read the tuple data from the parent ScanFile
            void read();

            /// Write the tuple data to the parent ScanFile
            void write();
            ///@}
    };
}

#endif
