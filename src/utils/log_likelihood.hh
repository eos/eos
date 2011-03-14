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

#ifndef EOS_GUARD_SRC_UTILS_LIKELIHOOD_HH
#define EOS_GUARD_SRC_UTILS_LIKELIHOOD_HH 1

#include <src/observable.hh>
#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern.hh>

#include <vector>

namespace eos
{
    // Forward declaration.
    class LogLikelihoodTest;


    /*!
     * LogLikelihood handles a set of ObservablePtr with associated measurement data.
     *
     * Access to any LogLikelihood is coherent, i.e., changes to one object will propagate
     * to every other object copy. To create an independent instance, use clone().
     */
    class LogLikelihood :
        public PrivateImplementationPattern<LogLikelihood>
    {
        private:
            /*!
             * Retrieve the observables and their values at last evaluation
             * @return
             */
            const std::vector<std::tuple<ObservablePtr, double>>&
            predictions() const;

        public:
            friend class LogLikelihoodTest;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param parameters  The Parameters object to which all further ObservablePtr objects must be bound.
             */
            LogLikelihood(const Parameters & parameters);

            /*!
             * Destructor.
             */
            ~LogLikelihood();
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
            void add(const ObservablePtr & observable, const double & min,
                    const double & central, const double & max);

            /*!
             * Create an independent instance of this LogLikelihood that uses the same set of observables and measurements.
             */
            LogLikelihood clone() const;

            /*!
             * Retrieve the underlying Parameters object.
             */
            Parameters parameters() const;

            /*!
             * Evaluate the log likelihood, i.e., return @f[2 \log \mathcal{L} = -\sum_i \chi^2_i + C\right]@f].
             * @note: all observables are recalculated
             */
            double operator()();

            /*!
             * Evaluate the likelihood, but recalculate only those observables which depend on par
             * with given id
             *
             * @warning Do not use this call before __all__ observables have been added to the likelihood!
             * Since the dependence is checked only the first time a parameter id is encountered, adding
             * an observable later will not establish the link between the id and the observable, leading
             * to incorrect values of the likelihood
             */
            double operator()(const Parameter::Id& id);

            /*!
             * Reload previous observable values.
             *
             * @note this is only useful in conjunction with operator()(const Parameter& par)
             * and the Markov Chain sampler
             */
            void reset();
            ///@}
    };

    typedef std::shared_ptr<LogLikelihood> LogLikelihoodPtr;
}

#endif
