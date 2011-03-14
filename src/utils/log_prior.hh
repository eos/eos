/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#ifndef EOS_GUARD_SRC_UTILS_PRIOR_HH
#define EOS_GUARD_SRC_UTILS_PRIOR_HH 1

#include <src/utils/log_prior-fwd.hh>
#include <src/utils/parameters.hh>
#include <src/utils/wrapped_forward_iterator.hh>

#include <vector>
#include <memory>

namespace eos
{
    /*!
     * Base class for log(prior) distributions.
     *
     * Has a container of subclasses describing independent 1..k dimensional prior distributions.
     * Taken together, they specify the full N dimensional prior. Any actual calculation is done
     * by the subclasses.
     */
    class LogPrior
    {
        protected:
            /// Our associated Parameters object.
            Parameters _parameters;

            /// All parameters for which this prior provides information.
            std::vector<ParameterDescription> _parameter_descriptions;

        public:
            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param parameters The Parameters objects from which we evaluate.
             */
            LogPrior(const Parameters & parameters);

            /*!
             * Create a clone (independent object) of this LogPrior object,
             *
             * @param parameters The Parameters object from which the clone shall evaluate.
             * @return The clone
             */
            virtual LogPriorPtr clone(const Parameters & parameters) const = 0;
            ///@}


            ///@name Iteration over descriptions
            ///@{
            struct IteratorTag;
            typedef WrappedForwardIterator<IteratorTag, ParameterDescription> Iterator;

            Iterator begin();
            Iterator end();
            ///@}

            ///@name Accessors
            ///@{
            /*!
             * Evaluate the natural logarithm of the prior.
             */
            virtual double operator() () const = 0;
            ///@}

            ///@name Named constructors for 1D prior distributions
            ///@{
            static LogPriorPtr Gauss(const Parameters & parameters, const std::string & name, const ParameterRange & range,
                    const double & lower, const double & central, const double & upper);
            static LogPriorPtr Flat(const Parameters & parameters, const std::string & name, const ParameterRange & range);
            ///@}
    };
}

#endif
