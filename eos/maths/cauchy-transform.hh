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

#ifndef EOS_GUARD_EOS_MATHS_CAUCHY_TRANSFORM_HH
#define EOS_GUARD_EOS_MATHS_CAUCHY_TRANSFORM_HH 1

#include <eos/maths/complex.hh>

#include <memory>
#include <span>
#include <string>

namespace eos
{
    /*!
     * The Cauchy transform C_R(zeta) = int du R(u) / (zeta - u) of a resolution kernel R, for Im zeta > 0.
     *
     * Convolving the pole term Re[c / (x - z)], Im z < 0, with R yields Re[c C_R(x - z)].
     */
    class CauchyTransform
    {
        public:
            virtual ~CauchyTransform();

            // Set the kernel parameters, in the kernel's argument order after the offset.
            virtual void prepare(std::span<const double> parameters) = 0;

            virtual complex<double> evaluate(const complex<double> & zeta) const = 0;

            // The Cauchy transform of the kernel with the given name, such as 'Kernel::Gaussian', for
            // |Re zeta - mu| up to 'extent'. Throws InternalError for a kernel without a Cauchy transform.
            static std::unique_ptr<CauchyTransform> make(const std::string & kernel, const double & extent);
    };
} // namespace eos

#endif
