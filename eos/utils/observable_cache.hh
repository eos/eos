/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
 * Copyright (c) 2011 Frederik Beaujean
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

#ifndef EOS_GUARD_SRC_UTILS_OBSERVABLE_CACHE_HH
#define EOS_GUARD_SRC_UTILS_OBSERVABLE_CACHE_HH 1

#include <eos/observable.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    class ObservableCache :
        public PrivateImplementationPattern<ObservableCache>
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
            typedef unsigned Id;

            /*!
             * Add a given observable to the cache and return its unique Id.
             *
             * @param observable The observable which shall be added to the cache.
             */
            Id add(const ObservablePtr & observable);

            /// Update the predictions for all observables.
            void update();

            /*!
             * Update the predicition for all observables that depend on a given
             * Parameter.
             *
             * @param id The Parameter::Id for which all depending observables shall be updated.
             */
            void update(const Parameter::Id & id);

            /// Retrieve the cache's common Parameters object.
            Parameters parameters() const;

            /// Reset the cache to the previous predictions.
            void reset();

            /*!
             * Retrieve a unique observable by its ObservableCache::Id.
             *
             * @param id The ObservableCache::Id whose associated ObservablePtr shall be retrieved.
             */
            ObservablePtr observable(const ObservableCache::Id & id) const;

            /*!
             * Retrieve the prediction for a given observable from the cache.
             *
             * @param id The unique ObservableCache::Id whose associated observable's prediction shall be retrieved.
             */
            double operator[] (const ObservableCache::Id & id) const;

            /// Retrieve the number of independent predictions from the cache.
            unsigned size() const;

            /// Clone this cache whilst keeping the observables in the given order, i.e. all ids remain valid.
            ObservableCache clone(const Parameters & parameters) const;
    };
}

#endif
