/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2015 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_THREAD_POOL_HH
#define EOS_GUARD_EOS_UTILS_THREAD_POOL_HH 1

#include <eos/utils/instantiation_policy.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/ticket.hh>

#include <functional>

namespace eos
{
    class ThreadPool : public InstantiationPolicy<ThreadPool, Singleton>, public PrivateImplementationPattern<ThreadPool>
    {
        public:
            ThreadPool();

            ~ThreadPool();

            Ticket enqueue(const std::function<void(void)> & work);

            static ThreadPool * instance();

            void wait_for_free_capacity();

            unsigned number_of_threads() const;
    };
} // namespace eos

#endif
