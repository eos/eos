/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2009-2011 Ciaran McCreesh
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

#include <eos/utils/indirect-iterator-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <list>
#include <src/clients/cli_handler.hh>
#include <src/clients/cli_section.hh>

namespace eos
{
    template <> struct Implementation<cli::Section>
    {
            cli::Handler * const    handler;
            const std::string       name;
            std::list<cli::Group *> groups;

            Implementation(cli::Handler * const h, const std::string & s) :
                handler(h),
                name(s)
            {
            }
    };

    template <> struct WrappedForwardIteratorTraits<cli::Section::GroupsConstIteratorTag>
    {
            using UnderlyingIterator = IndirectIterator<std::list<cli::Group *>::const_iterator>;
    };

    template class WrappedForwardIterator<cli::Section::GroupsConstIteratorTag, const cli::Group>;

    namespace cli
    {
        Section::Section(cli::Handler * const h, const std::string & s) :
            PrivateImplementationPattern<cli::Section>(new Implementation<cli::Section>(h, s))
        {
            h->add(this);
        }

        Section::~Section() = default;

        Section::GroupsConstIterator
        Section::begin() const
        {
            return GroupsConstIterator(_imp->groups.cbegin());
        }

        Section::GroupsConstIterator
        Section::end() const
        {
            return GroupsConstIterator(_imp->groups.cend());
        }

        void
        Section::add(cli::Group * const g)
        {
            _imp->groups.push_back(g);
        }

        void
        Section::remove(cli::Group * const g)
        {
            _imp->groups.remove(g);
        }

        Handler *
        Section::handler() const
        {
            return _imp->handler;
        }

        const std::string
        Section::name() const
        {
            return _imp->name;
        }
    } // namespace cli
} // namespace eos
