/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Based upon 'condition_variable.cc' from Paludis, which is:
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

#include <eos/utils/condition_variable.hh>
#include <eos/utils/instantiation_policy-impl.hh>

namespace eos
{
    ConditionVariable::ConditionVariable() :
        _cond(new pthread_cond_t)
    {
        pthread_cond_init(_cond, 0);
    }

    ConditionVariable::~ConditionVariable()
    {
        pthread_cond_destroy(_cond);
        delete _cond;
    }

    void
    ConditionVariable::broadcast()
    {
        pthread_cond_broadcast(_cond);
    }

    void
    ConditionVariable::signal()
    {
        pthread_cond_signal(_cond);
    }

    void
    ConditionVariable::wait(Mutex & m)
    {
        pthread_cond_wait(_cond, m.mutex());
    }
} // namespace eos
