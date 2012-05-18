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


#ifndef EOS_GUARD_SRC_UTILS_HDF5_FWD_HH
#define EOS_GUARD_SRC_UTILS_HDF5_FWD_HH 1

namespace eos
{
    namespace hdf5
    {
        // Forward declarations.
        class AttributeHandle;

        template <typename T_> class Attribute;

        template <typename T_> class DataSet;
        template <typename T_> DataSet<T_> & operator<< (DataSet<T_> &, const typename T_::Type &);
        template <typename T_> DataSet<T_> & operator>> (DataSet<T_> &, typename T_::Type &);

        class File;

        struct Type;
        template <unsigned rank_, typename T_> class Array;
        template <typename ... T_> class Composite;
        template <typename T_> class Scalar;
    }
}

#endif
