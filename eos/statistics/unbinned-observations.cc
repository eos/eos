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

#include <eos/statistics/unbinned-observations.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/lock.hh>
#include <eos/utils/mutex.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <atomic>
#include <format>

namespace eos
{
    template <> struct Implementation<UnbinnedObservations>
    {
            std::size_t rank;

            // Identifies this set of observations; the counter advances once per successful set().
            uint64_t instance, counter;

            // Guards the exchange of `current` and the counter; the snapshots themselves are immutable.
            mutable Mutex mutex;

            std::shared_ptr<const UnbinnedObservations::Snapshot> current;

            Implementation(const std::size_t & rank) :
                rank(rank),
                instance(_next_instance()),
                counter(0u)
            {
            }

            static const std::vector<double> &
            check(const std::vector<double> & values, const std::size_t & rank)
            {
                if (0 == rank)
                {
                    throw InternalError("UnbinnedObservations: the rank must be positive");
                }

                if (values.empty())
                {
                    throw InternalError("UnbinnedObservations: the observations must not be empty");
                }

                if (values.size() % rank != 0)
                {
                    throw InternalError(std::format("UnbinnedObservations: {} coordinates are not a multiple of the rank {}", values.size(), rank));
                }

                return values;
            }

            static uint64_t
            _next_instance()
            {
                static std::atomic<uint64_t> next(1u);

                return next.fetch_add(1u, std::memory_order_relaxed);
            }
    };

    UnbinnedObservations::UnbinnedObservations(const std::vector<double> & values, const std::size_t & rank) :
        PrivateImplementationPattern<UnbinnedObservations>(new Implementation<UnbinnedObservations>(rank))
    {
        _imp->current = std::make_shared<const Snapshot>(Snapshot{ Implementation<UnbinnedObservations>::check(values, rank), Generation(_imp->instance, 0u) });
    }

    UnbinnedObservations::~UnbinnedObservations() {}

    void
    UnbinnedObservations::set(const std::vector<double> & values)
    {
        // copy the events before taking the lock, which then covers only the exchange
        auto next = std::make_shared<Snapshot>(Snapshot{ Implementation<UnbinnedObservations>::check(values, _imp->rank), Generation() });

        Lock lock(_imp->mutex);
        next->generation = Generation(_imp->instance, ++_imp->counter);
        _imp->current    = std::move(next);
    }

    std::shared_ptr<const UnbinnedObservations::Snapshot>
    UnbinnedObservations::snapshot() const
    {
        Lock lock(_imp->mutex);
        return _imp->current;
    }

    std::vector<double>
    UnbinnedObservations::values() const
    {
        return this->snapshot()->values;
    }

    std::size_t
    UnbinnedObservations::rank() const
    {
        return _imp->rank;
    }

    std::size_t
    UnbinnedObservations::size() const
    {
        return this->snapshot()->values.size() / _imp->rank;
    }

    UnbinnedObservations::Generation
    UnbinnedObservations::generation() const
    {
        return this->snapshot()->generation;
    }
} // namespace eos
