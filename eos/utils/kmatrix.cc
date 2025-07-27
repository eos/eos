/*
 * Copyright (c) 2019 Stephan Kürten
 * Copyright (c) 2019 Danny van Dyk
 * Copyright (c) 2021 Méril Reboud
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

#include <eos/maths/power-of.hh>
#include <eos/utils/kmatrix.hh>

namespace eos
{
    namespace kmatrix_utils
    {
        complex<double>
        blatt_weisskopf_factor(const unsigned & l, const complex<double> & z)
        {
            const complex<double> z2 = z * z;

            switch (l)
            {
                // We use the definitions of the PDG's resonance review
                // they agree with [CBHKSS:1995A], keeping only the denominator (numerators are accounted for in the K matrix definition)
                case 0:  return 1.0;
                case 1:  return std::sqrt(1.0 / (z2 + 1.0));
                case 2:  return std::sqrt(1.0 / (9.0 + z2 * (3.0 + z2)));
                case 3:  return std::sqrt(1.0 / (z2 * power_of<2>(z2 - 15.0) + 9.0 * power_of<2>(2.0 * z2 - 5.0))); // This one is wrong in [CBHKSS:1995A]
                case 4:  return std::sqrt(1.0 / (power_of<2>(power_of<2>(z2) - 45.0 * z2 + 105.0) + 25.0 * z2 * power_of<2>(2.0 * z2 - 21.0)));
                default: throw InternalError("Blatt-Weisskopf factors are not implemented for l > 4.");
            }
        }
    } // namespace kmatrix_utils
} // namespace eos
