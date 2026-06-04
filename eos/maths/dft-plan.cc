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

#include <eos/maths/dft-plan-impl.hh>

namespace eos
{
    namespace dft
    {
        template class Plan<1, Direction::Forward>;
        template class Plan<1, Direction::Backward>;

        template class Plan<2, Direction::Forward>;
        template class Plan<2, Direction::Backward>;

        template class Plan<3, Direction::Forward>;
        template class Plan<3, Direction::Backward>;

        template class Plan<4, Direction::Forward>;
        template class Plan<4, Direction::Backward>;
    } // namespace eos::dft
} // namespace eos
