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

#ifndef EOS_GUARD_EOS_STATISTICS_DENSITY_TEST_HH
#define EOS_GUARD_EOS_STATISTICS_DENSITY_TEST_HH

#include <eos/statistics/density.hh>

#include <vector>

namespace eos
{
    /*!
     * Wrapper class of a simple named parameter
     */
    class TestParameter :
        public Mutable
    {
        public:
            ///@name Basic Operatios
            ///@{
            /// Constructor
            TestParameter(const std::string & name, double value = 0);

            /// Destructor.
            virtual ~TestParameter();

            /// Make a copy of this Mutable.
            virtual MutablePtr clone() const;

            ///@name Access & Modification of the numeric Value
            ///@{
            /// Cast a Mutable to a double.
            virtual operator double () const;

            /// Retrieve a Mutable's numeric value.
            virtual double operator() () const;

            /// Retrieve a Mutable's numeric value.
            virtual double evaluate() const;

            /// Set a Mutable's numeric value.
            virtual const Mutable & operator= (const double &);

            /// Set a Mutable's numeric value.
            virtual void set(const double &);
            ///@}

            ///@name Access to Meta Data
            ///@{
            /// Retrieve the Parameter's name.
            virtual const std::string & name() const;
            ///@}
        private:
            std::string _name;
            double _value;
    };

    /*!
     * A wrapper around a multivariate scalar function
     */
    class TestDensity :
        public Density
    {
        public:
            typedef std::function<double (const std::vector<double> &)> WrappedDensity;

            TestDensity(const WrappedDensity &);

            virtual ~TestDensity();

            void add(const ParameterDescription & def);

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
            std::vector<ParameterDescription> _defs;
            mutable std::vector<double> _parameter_values;
    };

    TestDensity make_multivariate_unit_normal(const unsigned & ndim);
}

#endif // EOS_GUARD_EOS_STATISTICS_DENSITY_TEST_HH
