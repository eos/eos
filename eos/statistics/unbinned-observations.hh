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

#ifndef EOS_GUARD_EOS_STATISTICS_UNBINNED_OBSERVATIONS_HH
#define EOS_GUARD_EOS_STATISTICS_UNBINNED_OBSERVATIONS_HH 1

#include <eos/utils/private_implementation_pattern.hh>

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace eos
{
    /*!
     * The observed events of an unbinned likelihood, flat and row-major with one coordinate per axis.
     *
     * Copies share the events: replacing them through one copy replaces them for all, e.g. for an
     * unbinned likelihood block and its clones, which rebuild what depends on the events at their
     * next evaluation.
     */
    class UnbinnedObservations : public PrivateImplementationPattern<UnbinnedObservations>
    {
        public:
            /*!
             * Generation identifies the events held by a set of observations.
             *
             * Two instances compare equal if and only if they have been obtained from the same
             * set of observations without an intervening set().
             */
            class Generation
            {
                private:
                    uint64_t _instance;

                    uint64_t _counter;

                    Generation(const uint64_t & instance, const uint64_t & counter) :
                        _instance(instance),
                        _counter(counter)
                    {
                    }

                    friend class UnbinnedObservations;

                public:
                    /// Create a generation that compares unequal to any generation of any set of observations.
                    Generation() :
                        _instance(0u),
                        _counter(0u)
                    {
                    }

                    bool operator== (const Generation &) const = default;
            };

            /*!
             * Constructor. Throws InternalError if the values are empty or not a multiple of the rank.
             *
             * @param values The observed events, flat and row-major.
             * @param rank   The number of coordinates per event.
             */
            UnbinnedObservations(const std::vector<double> & values, const std::size_t & rank);

            ~UnbinnedObservations();

            /// Replace the events, with the same checks as the constructor, and advance the generation.
            void set(const std::vector<double> & values);

            /// The events, flat and row-major; valid until the next set().
            std::span<const double> values() const;

            /// The number of coordinates per event.
            std::size_t rank() const;

            /// The number of events.
            std::size_t size() const;

            /// The current generation; it changes whenever set() succeeds.
            Generation generation() const;
    };
} // namespace eos

#endif
