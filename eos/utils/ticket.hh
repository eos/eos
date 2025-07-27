/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef EOS_GUARD_EOS_UTILS_TICKET_HH
#define EOS_GUARD_EOS_UTILS_TICKET_HH 1

#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    // Forward declaration for Ticket.
    class TicketList;

    /**
     * Ticket is used by all asynchronous function calls to relay/query
     * information on a function's completion status.
     */
    class Ticket : public PrivateImplementationPattern<Ticket>
    {
        public:
            /// \name Friends of Ticket
            /// \{

            friend class TicketList;

            /// \}

            /// \name Basic Operations
            /// \{

            /// Constructor.
            Ticket();

            /// Destructor.
            ~Ticket();

            /// \}

            /// Mark ticket as completed.
            void mark();

            /// Wait for ticket completion.
            void wait() const;
    };

    /**
     * TicketList is used when multiple tickets are needed.
     */
    class TicketList : public InstantiationPolicy<TicketList, NonCopyable>, public PrivateImplementationPattern<TicketList>
    {
        public:
            /// \name Basic Operations
            /// \{

            /// Constructor.
            TicketList();

            /// Destructor.
            ~TicketList();

            /// \}

            /// Push a ticket to the back of the list.
            void push_back(const Ticket & ticket);

            /// Wait for ticket completion.
            void wait() const;
    };
} // namespace eos

#endif
