/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_STATISTICS_SIMPLE_PARAMETERS_HH
#define EOS_GUARD_EOS_STATISTICS_SIMPLE_PARAMETERS_HH

#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <vector>

namespace eos
{
    class SimpleParameters;

    class SimpleParameter :
        public Mutable
    {
        public:
            friend struct Implementation<SimpleParameters>;

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
            typedef WrappedForwardIterator<IteratorTag, const ParameterDescription> Iterator;

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

    extern template class WrappedForwardIterator<SimpleParameters::IteratorTag, const ParameterDescription>;
}
#endif
