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
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

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

    class SimpleParameters;

    class SimpleParameter :
        public Mutable
    {
        public:
            friend class Implementation<SimpleParameters>;

            /*!
             * A unique number index of this parameter
             */
            typedef size_t Index;

            ///@name Basic Operations
            ///@{
            /// Destructor.
            virtual ~SimpleParameter();

            /*! Make a copy of this Mutable.
             * @note It is tied to an instance of SimpleParameters, hence not an independent copy.
             */
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
            /// Constructor
            SimpleParameter(const std::string & name, const Index & index,
                            const std::shared_ptr<std::vector<double>> & parameters);

            std::string _name;
            Index _index;
            std::shared_ptr<std::vector<double> > _parameters;
    };

    class SimpleParameters :
            public PrivateImplementationPattern<SimpleParameters>
    {
        public:
            ///@name Basic methods
            ///@{
            SimpleParameters();

            /*!
             * Destructor.
             */
            ~SimpleParameters();

            /// Yields and independent copy.
            SimpleParameters clone() const;

            ///@}

            ///@name Iteration
            ///@{
            struct IteratorTag;
            typedef WrappedForwardIterator<IteratorTag, ParameterDescription> Iterator;

            Iterator begin() const;
            Iterator end() const;
            ///@}

            ///@name Parameter access
            ///@{
            /*!
             * Declare a new parameter.
             *
             * @param name  Name of the new parameter to be declared.
             * @param value (Optional) value for the new parameter.
             */
            SimpleParameter & declare(const std::string & name, const double & min,
                                    const double & max, bool nuisance=false);

            /*!
             * Raw access to the values of all parameters.
             */
            const std::vector<double> & values() const;

            /*!
             * Retrieve a parameter's Parameter object by name.
             *
             * @param name  The name of the Parameter that shall be retrieved.
             */
            SimpleParameter & operator[] (const std::string & name) const;

            /*!
             * Retrieve a parameter's Parameter object by id.
             *
             * @param id    The id of the Parameter that shall be retrieved.
             */
            SimpleParameter & operator[] (const SimpleParameter::Index & id) const;
            ///@}

            /*!
             * Compare two instances of Parameters on inequality of their
             * underlying implementations.
             *
             * @param rhs   The right hand side of the binary != operator.
             */
            bool operator!= (const SimpleParameters & rhs) const;
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
