/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Based upon 'mutex.cc' from Paludis, which is:
 *     Copyright (c) 2007 Ciaran McCreesh
 *
 * * This file is part of the EOS program. EOS is free software;
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

#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/mutex.hh>

using namespace eos;

Mutex::Mutex() :
    _attr(new pthread_mutexattr_t),
    _mutex(new pthread_mutex_t)
{
    pthread_mutexattr_init(_attr);
    pthread_mutexattr_settype(_attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init(_mutex, _attr);
}

Mutex::~Mutex()
{
    pthread_mutex_destroy(_mutex);
    pthread_mutexattr_destroy(_attr);

    delete _mutex;
    delete _attr;
}
