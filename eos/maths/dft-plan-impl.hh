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

#ifndef EOS_GUARD_EOS_MATHS_DFT_PLAN_IMPL_HH
#define EOS_GUARD_EOS_MATHS_DFT_PLAN_IMPL_HH 1

#include <eos/maths/dft-plan.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <format>

#include <fftw3.h>

namespace eos
{
    namespace dft
    {
        namespace impl
        {
            template <std::size_t rank_>
            inline std::array<std::size_t, rank_> frequency_dimensions(const std::array<std::size_t, rank_> & dimensions)
            {
                std::array<std::size_t, rank_> frequency_dimensions;

                // For a real-valued input of size N, the DFT is symmetric and only the first N/2 + 1 complex coefficients are unique.
                // Therefore, the size of the frequency domain container is N/2 + 1 in the last dimension, and the same as the input dimensions in the other dimensions.
                frequency_dimensions[rank_ - 1] = dimensions[rank_ - 1] / 2 + 1;

                for (std::size_t i = 0; i + 1 < rank_; ++i)
                {
                    frequency_dimensions[i] = dimensions[i];
                }

                return frequency_dimensions;
            }

            template <std::size_t rank_>
            inline void check_dimensions(const std::array<std::size_t, rank_> & dimensions)
            {
                // The real-to-complex (and complex-to-real) transforms implemented here rely on each dimension being even,
                // both for the size N/2 + 1 of the frequency-domain container and for the symmetry assumptions of the transform.
                for (std::size_t i = 0; i < rank_; ++i)
                {
                    if (dimensions[i] % 2 != 0)
                    {
                        throw InternalError(std::format("dft::Plan<{}>: dimension {} (= {}) must be even", rank_, i, dimensions[i]));
                    }
                }
            }

            template <std::size_t rank_, Direction direction_>
            struct PlanTraits;

            template <>
            struct PlanTraits<1, Direction::Forward>
            {
                static fftw_plan alloc(const std::array<std::size_t, 1> & dimensions, const dft::Container<double, 1> & time_domain_container, const dft::Container<std::complex<double>, 1> & frequency_domain_container)
                {
                    static const unsigned flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;

                    double * time_domain_data = const_cast<double *>(time_domain_container.data());
                    fftw_complex * frequency_domain_data = const_cast<fftw_complex *>(reinterpret_cast<const fftw_complex *>(frequency_domain_container.data()));

                    return fftw_plan_dft_r2c_1d(dimensions[0], time_domain_data, frequency_domain_data, flags);
                }
            };

            template <>
            struct PlanTraits<1, Direction::Backward>
            {
                static fftw_plan alloc(const std::array<std::size_t, 1> & dimensions, const dft::Container<double, 1> & time_domain_container, const dft::Container<std::complex<double>, 1> & frequency_domain_container)
                {
                    static const unsigned flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;

                    double * time_domain_data = const_cast<double *>(time_domain_container.data());
                    fftw_complex * frequency_domain_data = const_cast<fftw_complex *>(reinterpret_cast<const fftw_complex *>(frequency_domain_container.data()));

                    return fftw_plan_dft_c2r_1d(dimensions[0], frequency_domain_data, time_domain_data, flags);
                }
            };

            template <>
            struct PlanTraits<2, Direction::Forward>
            {
                static fftw_plan alloc(const std::array<std::size_t, 2> & dimensions, const dft::Container<double, 2> & time_domain_container, const dft::Container<std::complex<double>, 2> & frequency_domain_container)
                {
                    static const unsigned flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;

                    double * time_domain_data = const_cast<double *>(time_domain_container.data());
                    fftw_complex * frequency_domain_data = const_cast<fftw_complex *>(reinterpret_cast<const fftw_complex *>(frequency_domain_container.data()));

                    return fftw_plan_dft_r2c_2d(dimensions[0], dimensions[1], time_domain_data, frequency_domain_data, flags);
                }
            };

            template <>
            struct PlanTraits<2, Direction::Backward>
            {
                static fftw_plan alloc(const std::array<std::size_t, 2> & dimensions, const dft::Container<double, 2> & time_domain_container, const dft::Container<std::complex<double>, 2> & frequency_domain_container)
                {
                    static const unsigned flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;

                    double * time_domain_data = const_cast<double *>(time_domain_container.data());
                    fftw_complex * frequency_domain_data = const_cast<fftw_complex *>(reinterpret_cast<const fftw_complex *>(frequency_domain_container.data()));

                    return fftw_plan_dft_c2r_2d(dimensions[0], dimensions[1], frequency_domain_data, time_domain_data, flags);
                }
            };

            template <>
            struct PlanTraits<3, Direction::Forward>
            {
                static fftw_plan alloc(const std::array<std::size_t, 3> & dimensions, const dft::Container<double, 3> & time_domain_container, const dft::Container<std::complex<double>, 3> & frequency_domain_container)
                {
                    static const unsigned flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;

                    double * time_domain_data = const_cast<double *>(time_domain_container.data());
                    fftw_complex * frequency_domain_data = const_cast<fftw_complex *>(reinterpret_cast<const fftw_complex *>(frequency_domain_container.data()));

                    return fftw_plan_dft_r2c_3d(dimensions[0], dimensions[1], dimensions[2], time_domain_data, frequency_domain_data, flags);
                }
            };

            template <>
            struct PlanTraits<3, Direction::Backward>
            {
                static fftw_plan alloc(const std::array<std::size_t, 3> & dimensions, const dft::Container<double, 3> & time_domain_container, const dft::Container<std::complex<double>, 3> & frequency_domain_container)
                {
                    static const unsigned flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;

                    double * time_domain_data = const_cast<double *>(time_domain_container.data());
                    fftw_complex * frequency_domain_data = const_cast<fftw_complex *>(reinterpret_cast<const fftw_complex *>(frequency_domain_container.data()));

                    return fftw_plan_dft_c2r_3d(dimensions[0], dimensions[1], dimensions[2], frequency_domain_data, time_domain_data, flags);
                }
            };

            template <>
            struct PlanTraits<4, Direction::Forward>
            {
                static fftw_plan alloc(const std::array<std::size_t, 4> & dimensions, const dft::Container<double, 4> & time_domain_container, const dft::Container<std::complex<double>, 4> & frequency_domain_container)
                {
                    static const unsigned flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;

                    std::array<int, 4> n;
                    std::transform(dimensions.begin(), dimensions.end(), n.begin(),
                        [](std::size_t d) { return static_cast<int>(d); });

                    double * time_domain_data = const_cast<double *>(time_domain_container.data());
                    fftw_complex * frequency_domain_data = const_cast<fftw_complex *>(reinterpret_cast<const fftw_complex *>(frequency_domain_container.data()));

                    return fftw_plan_dft_r2c(4, n.data(), time_domain_data, frequency_domain_data, flags);
                }
            };

            template <>
            struct PlanTraits<4, Direction::Backward>
            {
                static fftw_plan alloc(const std::array<std::size_t, 4> & dimensions, const dft::Container<double, 4> & time_domain_container, const dft::Container<std::complex<double>, 4> & frequency_domain_container)
                {
                    static const unsigned flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;

                    std::array<int, 4> n;
                    std::transform(dimensions.begin(), dimensions.end(), n.begin(),
                        [](std::size_t d) { return static_cast<int>(d); });

                    double * time_domain_data = const_cast<double *>(time_domain_container.data());
                    fftw_complex * frequency_domain_data = const_cast<fftw_complex *>(reinterpret_cast<const fftw_complex *>(frequency_domain_container.data()));

                    return fftw_plan_dft_c2r(4, n.data(), frequency_domain_data, time_domain_data, flags);
                }
            };
        }
    }

    template <std::size_t rank_, dft::Direction direction_>
    struct Implementation<dft::Plan<rank_, direction_>>
    {
        dft::Container<double, rank_> time_domain_container;
        dft::Container<std::complex<double>, rank_> frequency_domain_container;

        fftw_plan plan;

        std::array<std::size_t, rank_> dimensions;

        Implementation(const std::array<std::size_t, rank_> & dimensions) :
            // validate the dimensions before allocating any containers or the FFTW plan
            time_domain_container((dft::impl::check_dimensions(dimensions), dimensions)),
            frequency_domain_container(dft::impl::frequency_dimensions(dimensions)),
            plan(dft::impl::PlanTraits<rank_, direction_>::alloc(dimensions, time_domain_container, frequency_domain_container)),
            dimensions(dimensions)
        {
        }

        ~Implementation()
        {
            fftw_destroy_plan(this->plan);
        }
    };

    namespace dft
    {
        template <std::size_t rank_, Direction direction_>
        Plan<rank_, direction_>::Plan(const std::array<std::size_t, rank_> & dimensions) :
            PrivateImplementationPattern<dft::Plan<rank_, direction_>>(new Implementation<dft::Plan<rank_, direction_>>(dimensions))
        {
        }

        template <std::size_t rank_, Direction direction_>
        Plan<rank_, direction_>::Plan(Plan && other) = default;

        template <std::size_t rank_, Direction direction_>
        Plan<rank_, direction_> & Plan<rank_, direction_>::operator=(Plan && other) = default;

        template <std::size_t rank_, Direction direction_>
        Plan<rank_, direction_>::~Plan() = default;

        template <std::size_t rank_, Direction direction_>
        const std::array<std::size_t, rank_> & Plan<rank_, direction_>::dimensions() const
        {
            return this->_imp->dimensions;
        }

        template <std::size_t rank_, Direction direction_>
        void Plan<rank_, direction_>::transform()
        {
            fftw_execute(this->_imp->plan);
        }

        template <std::size_t rank_, Direction direction_>
        Container<double, rank_> & Plan<rank_, direction_>::time_domain_container()
        {
            return this->_imp->time_domain_container;
        }

        template <std::size_t rank_, Direction direction_>
        Container<std::complex<double>, rank_> & Plan<rank_, direction_>::frequency_domain_container()
        {
            return this->_imp->frequency_domain_container;
        }
    } // namespace eos::dft
} // namespace eos

#endif
