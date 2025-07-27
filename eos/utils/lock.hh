/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Based upon 'mutex.hh' from Paludis, which is:
 *     Copyright (c) 2007 Ciaran McCreesh
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

#ifndef EOS_GUARD_EOS_UTILS_LOCK_HH
#define EOS_GUARD_EOS_UTILS_LOCK_HH 1

#include <eos/utils/instantiation_policy.hh>
#include <eos/utils/mutex.hh>

namespace eos
{
    class Lock : public InstantiationPolicy<Lock, NonCopyable>
    {
        private:
            /// Our mutex.
            Mutex * const _mutex;

        public:
            /// (Explicit) constructor.
            explicit Lock(Mutex &);

            /// Destructor.
            ~Lock();
    };

    class TryLock : public InstantiationPolicy<TryLock, NonCopyable>
    {
        private:
            /// Our mutex.
            Mutex * _mutex;

        public:
            /// (Explicit) constructor.
            explicit TryLock(Mutex & e);

            /// Destructor.
            ~TryLock();

            /// Return true if the lock worked.
            bool
            operator() () const
            {
                return 0 != _mutex;
            }
    };
} // namespace eos

#endif
