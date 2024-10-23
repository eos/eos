/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2005-2011 Ciaran McCreesh
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

#ifndef PALUDIS_GUARD_ARGS_ARGS_GROUP_HH
#define PALUDIS_GUARD_ARGS_ARGS_GROUP_HH 1

#include <eos/utils/private_implementation_pattern.hh>

#include <src/clients/cli_option.hh>
#include <string>

/** \file
 * Declarations for Group.
 */
namespace eos
{
    namespace cli
    {
        class Section;

        /**
         * Contains a related group of command line arguments.
         */
        class Group : public PrivateImplementationPattern<cli::Group>
        {
            private:
                const std::string _name;
                const std::string _description;

                Section * _section;

            public:
                /**
                 * Remove this group from our section.
                 */
                void remove();

                /**
                 * Fetch our section.
                 */
                Section *
                section() const
                {
                    return _section;
                }

                /**
                 * Add an Option instance (called by the Option
                 * constructor).
                 */
                void add(Option * const value);

                /**
                 * Remove an Option instance (called by
                 * Option::remove).  Calls Group::remove() if
                 * that would leave us with no Options.
                 */
                void remove(Option * const value);

                ///\name Iterate over our Options.
                ///\{

                struct ConstIteratorTag;
                using ConstIterator = WrappedForwardIterator<ConstIteratorTag, Option * const>;

                ConstIterator begin() const;
                ConstIterator end() const;

                ///\}

                ///\name Basic operations
                ///\{

                Group(Section * s, const std::string & name, const std::string & description);

                ~Group();

                Group(const Group &)              = delete;
                Group & operator= (const Group &) = delete;

                ///\}

                /**
                 * Fetch our name.
                 */
                const std::string &
                name() const
                {
                    return _name;
                }

                /**
                 * Fetch our description.
                 */
                const std::string &
                description() const
                {
                    return _description;
                }
        };
    } // namespace cli

    extern template class WrappedForwardIterator<cli::Group::ConstIteratorTag, cli::Option * const>;
} // namespace eos

#endif
