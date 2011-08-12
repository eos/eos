/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_OPTIONS_HH
#define EOS_GUARD_SRC_UTILS_OPTIONS_HH 1

#include <eos/utils/exception.hh>
#include <eos/utils/private_implementation_pattern.hh>

#include <string>

namespace eos
{
    /*!
     * UnknownOptionError is thrown when no option of a given name could be
     * found.
     */
    struct UnknownOptionError :
        public Exception
    {
        UnknownOptionError(const std::string & key) throw ();
    };

    /*!
     * Options keeps the set of all string options for any Observable.
     */
    class Options :
        public PrivateImplementationPattern<Options>
    {
        public:
            friend Options operator+ (const Options &, const Options &);

            ///@name Basic Functions
            ///@{

            /// Constructor.
            Options();

            /*!
             * Constructor.
             *
             * Create an instance of Kinematics with a given set of initial
             * options.
             *
             * @param options The set of initial options from which this object shall be constructed.
             */
            Options(const std::initializer_list<std::pair<std::string, std::string>> & options);

            /// Destructor.
            ~Options();

            /// Equality comparison operator.
            bool operator== (const Options & rhs) const;

            /// Inequality comparison operator.
            bool operator!= (const Options & rhs) const;
            ///@}

            ///@name Access
            ///@{
            const std::string & operator[] (const std::string & key) const;

            bool has(const std::string & key) const;

            void set(const std::string & key, const std::string & value = "");

            std::string get(const std::string & key, const std::string & default_value = "") const;

            std::string as_string() const;
            ///@}
    };

    /// Merge operator.
    Options operator+ (const Options & lhs, const Options & rhs);
}

#endif
