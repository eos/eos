/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Stephan Jahn
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

#include <eos/utils/options_C.hh>

extern "C" {
    using namespace eos;

    Options *
    EOS_Options_new()
    {
        return new Options();
    }

    void
    EOS_Options_delete(Options * options)
    {
        delete options;
        options = nullptr;
    }

    void
    EOS_Options_set(Options * options, const char * key, const char * value)
    {
        options->set(key, value);
    }
}
