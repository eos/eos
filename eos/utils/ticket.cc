/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2015 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the EOS program. EOS is free software;
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

#include <eos/utils/condition_variable.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/lock.hh>
#include <eos/utils/mutex.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/ticket.hh>

#include <list>
#include <memory>

namespace eos
{
    template <> struct Implementation<Ticket>
    {
            Mutex mutex;

            ConditionVariable completion;

            bool completed;

            Implementation() :
                completed(false)
            {
            }
    };

    Ticket::Ticket() :
        PrivateImplementationPattern<Ticket>(new Implementation<Ticket>)
    {
    }

    Ticket::~Ticket() {}

    void
    Ticket::mark()
    {
        Lock l(_imp->mutex);

        _imp->completed = true;
        _imp->completion.signal();
    }

    void
    Ticket::wait() const
    {
        Lock l(_imp->mutex);

        while (! _imp->completed)
        {
            _imp->completion.wait(_imp->mutex);
        }
    }

    template <> struct Implementation<TicketList>
    {
            std::list<std::shared_ptr<Implementation<Ticket>>> tickets;
    };

    TicketList::TicketList() :
        PrivateImplementationPattern<TicketList>(new Implementation<TicketList>)
    {
    }

    TicketList::~TicketList() {}

    void
    TicketList::push_back(const Ticket & ticket)
    {
        _imp->tickets.push_back(ticket._imp);
    }

    void
    TicketList::wait() const
    {
        while (! _imp->tickets.empty())
        {
            std::shared_ptr<Implementation<Ticket>> ticket(_imp->tickets.front());

            {
                Lock l(ticket->mutex);

                while (! ticket->completed)
                {
                    ticket->completion.wait(ticket->mutex);
                }
            }

            _imp->tickets.pop_front();
        }
    }
} // namespace eos
