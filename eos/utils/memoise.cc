/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2011 Danny van Dyk
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

#include <eos/utils/memoise.hh>

namespace eos
{
    MemoisationControl::MemoisationControl() :
        _mutex(new Mutex)
    {
    }

    MemoisationControl::~MemoisationControl()
    {
        delete _mutex;
    }

    void
    MemoisationControl::register_clear_function(const std::function<void()> & clear_function)
    {
        Lock l(*_mutex);

        _clear_functions.push_back(clear_function);
    }

    void
    MemoisationControl::clear()
    {
        Lock l(*_mutex);

        for (auto & _clear_function : _clear_functions)
        {
            _clear_function();
        }
    }
} // namespace eos
