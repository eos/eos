/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Danny van Dyk
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

#include <eos/utils/diagnostics.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <vector>

namespace eos
{
    template <> struct WrappedForwardIteratorTraits<Diagnostics::IteratorTag>
    {
            using UnderlyingIterator = std::vector<Diagnostics::Entry>::const_iterator;
    };
    template class WrappedForwardIterator<Diagnostics::IteratorTag, const Diagnostics::Entry>;

    template <> struct Implementation<Diagnostics>
    {
            std::vector<Diagnostics::Entry> entries;
    };

    Diagnostics::Diagnostics() :
        PrivateImplementationPattern<Diagnostics>(new Implementation<Diagnostics>())
    {
    }

    Diagnostics::~Diagnostics() {}

    void
    Diagnostics::add(const Entry & entry)
    {
        _imp->entries.push_back(entry);
    }

    Diagnostics::Iterator
    Diagnostics::begin() const
    {
        return Diagnostics::Iterator(_imp->entries.begin());
    }

    Diagnostics::Iterator
    Diagnostics::end() const
    {
        return Diagnostics::Iterator(_imp->entries.end());
    }

    unsigned
    Diagnostics::size() const
    {
        return _imp->entries.size();
    }
} // namespace eos
