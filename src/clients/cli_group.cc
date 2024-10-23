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

#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <list>
#include <src/clients/cli_group.hh>
#include <src/clients/cli_section.hh>

namespace eos
{
    /**
     * Imp data for Group.
     *
     * \ingroup grplibpaludisargs
     */
    template <> struct Implementation<cli::Group>
    {
            std::list<cli::Option *> args_options;
    };

    template <> struct WrappedForwardIteratorTraits<cli::Group::ConstIteratorTag>
    {
            using UnderlyingIterator = std::list<cli::Option *>::const_iterator;
    };

    namespace cli
    {
        Group::Group(Section * s, const std::string & our_name, const std::string & our_description) :
            PrivateImplementationPattern<cli::Group>(new Implementation<cli::Group>()),
            _name(our_name),
            _description(our_description),
            _section(s)
        {
            s->add(this);
        }

        void
        Group::remove()
        {
            _section->remove(this);
        }

        void
        Group::add(Option * const value)
        {
            /// \bug Should check for uniqueness of short and long names.
            _imp->args_options.push_back(value);
        }

        void
        Group::remove(Option * const value)
        {
            _imp->args_options.remove(value);
            if (_imp->args_options.empty())
            {
                remove();
            }
        }

        Group::~Group() = default;

        Group::ConstIterator
        Group::begin() const
        {
            return ConstIterator(_imp->args_options.begin());
        }

        Group::ConstIterator
        Group::end() const
        {
            return ConstIterator(_imp->args_options.end());
        }
    } // namespace cli

    template class WrappedForwardIterator<cli::Group::ConstIteratorTag, cli::Option * const>;
} // namespace eos
