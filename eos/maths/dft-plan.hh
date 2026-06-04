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

#ifndef EOS_GUARD_EOS_MATHS_DFT_PLAN_HH
#define EOS_GUARD_EOS_MATHS_DFT_PLAN_HH 1

#include <eos/maths/dft-container.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    namespace dft
    {
        enum class Direction
        {
            Forward,
            Backward
        };

        template <std::size_t rank_, Direction direction_>
        class Plan : public PrivateImplementationPattern<dft::Plan<rank_, direction_>>
        {
            public:
                Plan(const std::array<std::size_t, rank_> & dimensions);

                // A Plan owns an FFTW plan and its scratch buffers. Copying would either alias them (via the
                // shared_ptr PIMPL) or require rebuilding the FFTW plan, so Plan is move-only.
                Plan(const Plan & other) = delete;
                Plan(Plan && other);
                Plan & operator=(const Plan & other) = delete;
                Plan & operator=(Plan && other);
                ~Plan();

                const std::array<std::size_t, rank_> & dimensions() const;

                // Forward: transforms from the time domain to the frequency domain.
                // Backward: transforms from the frequency domain to the time domain.
                void transform();

                Container<double, rank_> & time_domain_container();
                Container<std::complex<double>, rank_> & frequency_domain_container();
        };

        extern template class Plan<1, Direction::Forward>;
        extern template class Plan<1, Direction::Backward>;

        extern template class Plan<2, Direction::Forward>;
        extern template class Plan<2, Direction::Backward>;

        extern template class Plan<3, Direction::Forward>;
        extern template class Plan<3, Direction::Backward>;

        extern template class Plan<4, Direction::Forward>;
        extern template class Plan<4, Direction::Backward>;
    } // namespace eos::dft
} // namespace eos

#endif
