/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2026 Danny van Dyk
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
#include <memory>

namespace eos
{
    class ThreadPool : public InstantiationPolicy<ThreadPool, Singleton>, public PrivateImplementationPattern<ThreadPool>
    {
        public:
            ThreadPool();

            ~ThreadPool();

            Ticket enqueue(const std::function<void(void)> & work);

            static ThreadPool * instance();

            /*!
             * A guard that a thread holds while it waits for the work it enqueued.
             *
             * An embedding runtime that serialises access to its own state -- the Python
             * interpreter, say -- must relinquish it for as long as a thread waits, since the
             * pool's threads need it in turn. Such a runtime derives from this class and
             * relinquishes that access for as long as one of its guards lives.
             */
            class WaitGuard
            {
                public:
                    virtual ~WaitGuard();
            };

            /// Sets the factory that supplies the wait guards.
            void set_wait_guard_factory(const std::function<std::unique_ptr<WaitGuard>()> & factory);

            /// Creates a wait guard, or an empty pointer if no factory has been set.
            std::unique_ptr<WaitGuard> wait_guard() const;

            void wait_for_free_capacity();

            unsigned number_of_threads() const;
    };
} // namespace eos

#endif
