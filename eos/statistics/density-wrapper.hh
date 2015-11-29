/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Frederik Beaujean
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

#ifndef EOS_GUARD_EOS_STATISTICS_DENSITY_WRAPPER_HH
#define EOS_GUARD_EOS_STATISTICS_DENSITY_WRAPPER_HH

#include <eos/statistics/simple-parameters.hh>
#include <eos/utils/density.hh>

#include <vector>

namespace eos
{
   /*!
     * A wrapper around a multivariate scalar function.
     *
     * Note that initially no parameters are defined but they have to be added.
     * The WrappedDensity is called with a vector containing as many elements as defined parameters.
     *
     * Update parameter values either via the iterator interface, or via parameters() and then access parameters
     * either by name or by index.
     */
    class DensityWrapper :
        public Density
    {
        public:
            typedef double (* RawDensity)(const std::vector<double> &);
            typedef std::function< double (const std::vector<double> &)> WrappedDensity;

            /// Initialize with a WrappedDensity, that could point for example to a member method
            DensityWrapper(const WrappedDensity &);

            /// Initialize with a RawDensity typical of a free-standing function
            DensityWrapper(RawDensity);

            virtual ~DensityWrapper();

            /*!
             * Add a parameter to the density
             *
             * @param name Parameter name.
             * @param min Minimum allowed value. Note that this limit is purely informative and not enforced.
             * @param max Maximum allowed value. Note that this limit is purely informative and not enforced.
             * @param nuisance A nuisance parameter or not.
             */
            void add_parameter(const std::string & name, const double & min,
                               const double & max, bool nuisance=false);

            /// Access to container managing the underlying parameters.
            SimpleParameters & parameters();

            /// Evaluate the density function at the current parameter point on the log scale.
            virtual double evaluate() const;

            /// Create an independent copy of this density function.
            virtual DensityPtr clone() const;

            /// Iterate over the parameters relevant to this density function.
            ///@{
            virtual Iterator begin() const;
            virtual Iterator end() const;
            ///@}
        private:
            WrappedDensity _density;
            SimpleParameters _parameters;
    };
}
#endif
