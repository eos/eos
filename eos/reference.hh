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

#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <set>

// forward declaration
namespace YAML
{
    class Emitter;
    class Node;
} // namespace YAML

namespace eos
{
    // forward declaration
    class References;

    /*!
     * Reference is used to keep track of the known references.
     */
    class Reference : public PrivateImplementationPattern<Reference>
    {
        private:
            Reference(Implementation<Reference> * imp);

        public:
            friend class Implementation<References>;

            /// Copy constructor
            Reference(const Reference &) = default;
            Reference(Reference &&)      = default;

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
    class References : public PrivateImplementationPattern<References>
    {
        public:
            /// Constructor.
            References();

            /// Destructor.
            ~References();

            ///@name Iteration over known constraints
            ///@{
            struct ReferenceIteratorTag;
            using ReferenceIterator = WrappedForwardIterator<ReferenceIteratorTag, const std::pair<const ReferenceName, ReferencePtr>>;

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
     * Base class for all users of Reference objects.
     */
    class ReferenceUser
    {
        protected:
            std::set<ReferenceName> _references;

        public:
            ReferenceUser()  = default;
            ~ReferenceUser() = default;

            ///@name Iteration over references
            ///@{
            struct ReferenceIteratorTag;
            using ReferenceIterator = WrappedForwardIterator<ReferenceIteratorTag, const ReferenceName>;

            ReferenceIterator begin_references() const;
            ReferenceIterator end_references() const;

            ///@}

            ///@name Access
            ///@{
            /*!
             * Add a given reference to our list of used references.
             *
             * @param name   The reference name that we use.
             */
            inline void
            uses(const ReferenceName & name)
            {
                this->_references.insert(name);
            }

            /*!
             * Convenience access to add an entire set of used references.
             */
            inline void
            uses(const std::set<ReferenceName> & names)
            {
                for (const auto & name : names)
                {
                    uses(name);
                }
            }

            /*!
             * Copy reference names of another ReferenceUser to our list of used names.
             *
             * @param user The other ReferenceUser whose reference we are going to copy.
             */
            inline void
            uses(const ReferenceUser & user)
            {
                this->uses(user._references);
            }

            ///@}
    };

    extern template class WrappedForwardIterator<ReferenceUser::ReferenceIteratorTag, const ReferenceName>;

    /*!
     * UnknownReferenceError is thrown when References encounters an unknown constraint name.
     */
    struct UnknownReferenceError : public Exception
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
    struct ReferencesInputFileParseError : public Exception
    {
            ReferencesInputFileParseError(const std::string & file, const std::string & msg) throw();
    };

    /*!
     * ReferencesInputFileNodeError is thrown when a malformed node is encountered
     * within the references file.
     */
    struct ReferencesInputFileNodeError : public Exception
    {
            ReferencesInputFileNodeError(const std::string & file, const std::string & node, const std::string & msg) throw();
    };

    /*!
     * ReferencesInputDuplicateError is thrown when a duplicate parameter entry is encountered when parsing
     * the references file.
     */
    struct ReferencesInputDuplicateError : public Exception
    {
            ReferencesInputDuplicateError(const std::string & file, const std::string & msg) throw();
    };

} // namespace eos

#endif
