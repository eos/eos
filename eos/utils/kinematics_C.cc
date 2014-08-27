/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Frederik Beaujean, Christoph Bobeth
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

#include <eos/utils/kinematics_C.hh>

extern "C" {
    using namespace eos;

    Kinematics *
    EOS_Kinematics_new()
    {
        return new Kinematics();
    }

    void
    EOS_Kinematics_delete(Kinematics * kinematics)
    {
        delete kinematics;
        kinematics = nullptr;
    }

    void
    EOS_Kinematics_set(Kinematics * kinematics, const char * key, const double value)
    {
        kinematics->declare(key, value);
    }
}
