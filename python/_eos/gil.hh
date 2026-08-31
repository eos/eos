/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#include "eos/utils/thread_pool.hh"

#include <boost/python.hpp>

#ifndef EOS_PYTHON__EOS_GIL_HH
#  define EOS_PYTHON__EOS_GIL_HH 1

namespace impl
{
    /*!
     * Attaches the calling thread to the interpreter for the duration of the scope.
     *
     * Every use of the Python API from a thread that EOS created itself -- a thread of the
     * ``ThreadPool``, say -- must be wrapped in this guard. Nesting is permitted, both within
     * itself and within ``ScopedGILRelease``.
     *
     * ``PyGILState_Ensure`` and ``PyGILState_Release`` express thread-state attachment rather
     * than lock acquisition, which is what a free-threaded interpreter requires as well: there,
     * attachment carries no serialisation, so the same code parallelises instead of blocking.
     */
    class ScopedGILAcquire
    {
        private:
            PyGILState_STATE _state;

        public:
            ScopedGILAcquire() :
                _state(PyGILState_Ensure())
            {
            }

            ~ScopedGILAcquire() { PyGILState_Release(_state); }

            ScopedGILAcquire(const ScopedGILAcquire &)              = delete;
            ScopedGILAcquire & operator= (const ScopedGILAcquire &) = delete;
    };

    /*!
     * Detaches the calling thread from the interpreter for the duration of the scope.
     *
     * A thread that dispatches work to threads of EOS's own making must hold this guard while
     * it waits, so that those threads can attach themselves via ``ScopedGILAcquire``. The
     * calling thread must not touch the Python API while detached. Constructing the guard on a
     * thread that is already detached does nothing.
     *
     * This is the pool's ``WaitGuard`` for the interpreter, and doubles as a plain scope guard.
     */
    class ScopedGILRelease : public eos::ThreadPool::WaitGuard
    {
        private:
            PyThreadState * _state;

        public:
            ScopedGILRelease() :
                _state(PyGILState_Check() ? PyEval_SaveThread() : nullptr)
            {
            }

            ~ScopedGILRelease()
            {
                if (nullptr != _state)
                {
                    PyEval_RestoreThread(_state);
                }
            }

            ScopedGILRelease(const ScopedGILRelease &)              = delete;
            ScopedGILRelease & operator= (const ScopedGILRelease &) = delete;
    };
} // namespace impl

#endif // EOS_PYTHON__EOS_GIL_HH
