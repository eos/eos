/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_REFERENCES_HH
#define EOS_GUARD_EOS_REFERENCES_HH 1

#include <eos/utils/reference-name.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

// forward declaration
namespace YAML
{
    class Emitter;
    class Node;
}

namespace eos
{
    // forward declaration
    class References;

    /*!
     * Reference is used to keep track of the known references.
     */
    class Reference :
        public PrivateImplementationPattern<Reference>
    {
        private:
            Reference(Implementation<Reference> * imp);

        public:
            friend class Implementation<References>;

            /// Copy constructor
            Reference(const Reference &) = default;
            Reference(Reference &&) = default;

            /// Destructor
            ~Reference();

            /// Return the reference's name
            const ReferenceName & name() const;

            /// Return the reference's authors
            const std::string & authors() const;

            /// Return the reference's title
            const std::string & title() const;

            /// Return the reference's eprint archive
            const std::string & eprint_archive() const;

            /// Return the reference's eprint id
            const std::string & eprint_id() const;

            /// Return the reference's inspire id
            const std::string & inspire_id() const;
    };

    using ReferencePtr = std::shared_ptr<const Reference>;

    /*!
     * Container around the known references
     */
    class References :
        public PrivateImplementationPattern<References>
    {
        public:
            /// Constructor.
            References();

            /// Destructor.
            ~References();

            ///@name Iteration over known constraints
            ///@{
            struct ReferenceIteratorTag;
            typedef WrappedForwardIterator<ReferenceIteratorTag, const std::pair<const ReferenceName, ReferencePtr>> ReferenceIterator;

            ReferenceIterator begin() const;
            ReferenceIterator end() const;
            ///@}

            /*!
             * Retrieve a Reference object by name.
             *
             * @param name  The name of the Reference that shall be retrieved.
             */
            ReferencePtr operator[] (const ReferenceName & name) const;
    };

    extern template class WrappedForwardIterator<References::ReferenceIteratorTag, const std::pair<const ReferenceName, ReferencePtr>>;

    /*!
     * UnknownReferenceError is thrown when References encounters an unknown constraint name.
     */
    struct UnknownReferenceError :
        public Exception
    {
        ///@name Basic Functions
        ///@{
        /*!
         * Constructor.
         *
         * @param name The offending reference name.
         */
        UnknownReferenceError(const ReferenceName & name);
        ///@}
    };

    /*!
     * ReferencesInputFileParseError is thrown when a malformed references file
     * cannot be parsed by libyaml-cpp.
     */
    struct ReferencesInputFileParseError :
        public Exception
    {
        ReferencesInputFileParseError(const std::string & file, const std::string & msg) throw ();
    };

    /*!
     * ReferencesInputFileNodeError is thrown when a malformed node is encountered
     * within the references file.
     */
    struct ReferencesInputFileNodeError :
        public Exception
    {
        ReferencesInputFileNodeError(const std::string & file, const std::string & node, const std::string & msg) throw ();
    };

    /*!
     * ReferencesInputDuplicateError is thrown when a duplicate parameter entry is encountered when parsing
     * the references file.
     */
    struct ReferencesInputDuplicateError :
        public Exception
    {
        ReferencesInputDuplicateError(const std::string & file, const std::string & msg) throw ();
    };

}

#endif
