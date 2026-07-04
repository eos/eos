/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2026 Danny van Dyk
 * Copyright (c) 2011      Frederik Beaujean
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

#ifndef EOS_GUARD_EOS_UTILS_OBSERVABLE_CACHE_HH
#define EOS_GUARD_EOS_UTILS_OBSERVABLE_CACHE_HH 1

#include <eos/observable.hh>
#include <eos/utils/observable_set.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/strong-typedef.hh>

#include <span>
#include <vector>

namespace eos
{
    class ObservableCache : public PrivateImplementationPattern<ObservableCache>
    {
        public:
            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param parameters   The Parameters object common to all cached Observables.
             */
            ObservableCache(const Parameters & parameters);

            /// Destructor.
            ~ObservableCache();
            ///@}

            ///@name Access
            ///@{
            struct ObservableIdTag;

            class ObservableId : public StrongTypedef<unsigned, ObservableCache::ObservableIdTag>
            {
                public:
                    using StrongTypedef::StrongTypedef;
            };

            /*!
             * Add a given observable to the cache and return its unique ObservableId.
             *
             * @param observable The observable which shall be added to the cache.
             */
            ObservableId add(const ObservablePtr & observable);

            struct BatchIdTag;

            /*!
             * A BatchId identifies a batch of computations added to the cache via add_batch().
             *
             * In contrast to an ObservableId, which addresses a single prediction, a BatchId
             * addresses a contiguous block of predictions, one per observable in the batch.
             */
            class BatchId : public StrongTypedef<unsigned, ObservableCache::BatchIdTag>
            {
                public:
                    using StrongTypedef::StrongTypedef;
            };

            /*!
             * Add a batch of computations to the cache and return its unique BatchId.
             *
             * The observables of a batch are evaluated into a contiguous, pre-allocated block of
             * memory upon update(). Different batches may be evaluated in parallel; the computations
             * within a single batch are partitioned across the available threads. The results are
             * retrieved as a contiguous span via operator[](const BatchId &).
             *
             * @param computations The observables which shall be added to the cache as a single batch.
             */
            BatchId add_batch(std::vector<ObservablePtr> && computations);

            /// Update the predictions for all observables.
            void update();

            /// Retrieve the cache's common Parameters object.
            Parameters parameters() const;

            /*!
             * Retrieve a unique observable by its ObservableCache::ObservableId.
             *
             * @param id The ObservableCache::ObservableId whose associated ObservablePtr shall be retrieved.
             */
            ObservablePtr observable(const ObservableCache::ObservableId & id) const;

            /*!
             * Retrieve the prediction for a given observable from the cache.
             *
             * @param id The unique ObservableCache::ObservableId whose associated observable's prediction shall be retrieved.
             */
            double operator[] (const ObservableCache::ObservableId & id) const;

            /*!
             * Retrieve the predictions for a given batch of computations from the cache.
             *
             * @param id The unique ObservableCache::BatchId whose associated batch's predictions shall be retrieved.
             */
            std::span<const double> operator[] (const ObservableCache::BatchId & id) const;

            /// Retrieve the number of independent predictions from the cache.
            unsigned size() const;

            struct IteratorTag;
            using Iterator = WrappedForwardIterator<IteratorTag, ObservablePtr>;
            Iterator begin() const;
            Iterator end() const;
            ///@}

            /// Clone this cache whilst keeping the observables in the given order, i.e. all ids remain valid.
            ObservableCache clone(const Parameters & parameters) const;
    };

    extern template class WrappedForwardIterator<ObservableCache::IteratorTag, ObservablePtr>;
} // namespace eos

#endif
