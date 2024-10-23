/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2009-2010 Ciaran McCreesh
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

#ifndef EOS_GUARD_SRC_CLIENTS_CLI_SECTION_HH
#define EOS_GUARD_SRC_CLIENTS_CLI_SECTION_HH 1

#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <string>

namespace eos
{
    namespace cli
    {
        class Group;
        class Handler;

        /**
         * Holds a number of Group instances.
         */
        class Section : public PrivateImplementationPattern<cli::Section>
        {
            public:
                Section(Handler * const, const std::string &);
                ~Section();

                struct GroupsConstIteratorTag;
                using GroupsConstIterator = WrappedForwardIterator<GroupsConstIteratorTag, const cli::Group>;
                GroupsConstIterator begin() const;
                GroupsConstIterator end() const;

                Handler *         handler() const;
                const std::string name() const;

                void add(Group * const);
                void remove(Group * const);
        };
    } // namespace cli

    extern template class WrappedForwardIterator<cli::Section::GroupsConstIteratorTag, const cli::Group>;
} // namespace eos

#endif
