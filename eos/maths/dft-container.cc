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

#include <eos/maths/dft-container-impl.hh>

 namespace eos
 {
    namespace dft
    {
        template class Container<double, 1>;
        template class Container<double, 2>;
        template class Container<double, 3>;
        template class Container<double, 4>;

        template class Container<std::complex<double>, 1>;
        template class Container<std::complex<double>, 2>;
        template class Container<std::complex<double>, 3>;
        template class Container<std::complex<double>, 4>;
    }
 }
