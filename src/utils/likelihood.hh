/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_LIKELIHOOD_HH
#define EOS_GUARD_SRC_UTILS_LIKELIHOOD_HH 1

#include <src/utils/observable.hh>
#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace eos
{
    /*!
     * Likelihood handles a set of ObservablePtr with associated measurement data.
     *
     * Access to any Likelihood is coherent, i.e., changes to one object will propagate
     * to every other object copy. To create an independent instance, use clone().
     */
    class Likelihood :
        public PrivateImplementationPattern<Likelihood>
    {
        public:
            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param parameters  The Parameters object to which all further ObservablePtr objects must be bound.
             */
            Likelihood(const Parameters & parameters);

            /*!
             * Destructor.
             */
            ~Likelihood();
            ///@}

            ///@name Access and evaluation
            ///@{
            /*!
             * Add an observable and its associated measurement.
             *
             * @param observable  The Observable that shall be added to the calculation of the likelihood.
             * @param min         The lower bound on the measurement.
             * @param central     The central value of the measurement.
             * @param max         The upper bound on the measurement.
             */
            void add(const ObservablePtr & observable, const double & min, const double & central, const double & max);

            /*!
             * Create an independent instance of this Likelihood that uses the same set of observables and measurements.
             */
            Likelihood clone() const;

            /*!
             * Retrieve the underlying Parameters object.
             */
            Parameters parameters() const;

            /*!
             * Evaluate the likelihood, i.e., return @f[\mathcal{L} = \exp\left[-\frac{1}{2}\sum_i \chi^2_i\right].@f]
             */
            double operator() () const;
            ///@}
    };
}

#endif
