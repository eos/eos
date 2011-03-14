/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_PARAMETERS_HH
#define EOS_GUARD_SRC_UTILS_PARAMETERS_HH 1

#include <src/utils/exception.hh>
#include <src/utils/private_implementation_pattern.hh>
#include <src/utils/random_number_engine.hh>
#include <src/utils/wrapped_forward_iterator.hh>

namespace eos
{
    /*!
     * UnknownParameterError is thrown when no parameter of a given
     * name could be found.
     */
    struct UnknownParameterError :
        public Exception
    {
        UnknownParameterError(const std::string & variable) throw ();
    };

    // Forward declaration.
    class Parameter;

    /*!
     * Parameters keeps the set of all numeric parameters for any Observable.
     *
     * Access to any Parameter or their values is coherent, i.e., changes to
     * a Parameter object will propagate to every other object with the same
     * parent Parameters and which handle the same parameter by name.
     */
    class Parameters :
        public PrivateImplementationPattern<Parameters>
    {
        private:
            ///@name Internal Data
            ///@{
            struct Data;
            ///@}

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor (private).
             *
             * Creates an instance of Parameters from a given implementation.
             * To create an instance of Parameters, see one of the named constructors.
             *
             * @param impl Implementation from which the Parameters object shall be constructed from.
             */
            Parameters(Implementation<Parameters> *);
            ///@}

        public:
            friend class Parameter;
            friend class Implementation<Parameter>;
            friend class Implementation<Parameters>;

            ///@name Basic Functions
            ///@{
            /*!
             * Named constructor.
             *
             * Creates an instance of Parameters with default values filled in.
             */
            static Parameters Defaults();

            Parameters clone() const;
            /*!
             * Destructor.
             */
            ~Parameters();
            ///@}

            ///@name Iteration
            ///@{
            struct IteratorTag;
            typedef WrappedForwardIterator<IteratorTag, Parameter> Iterator;

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
            Parameter declare(const std::string & name, const double & value = 0.0);

            /*!
             * Set a parameter's numeric value.
             *
             * @param name  The name of the parameter whose numeric value shall be changed.
             * @param value The parameter's new numeric value.
             */
            void set(const std::string & name, const double & value);

            /*!
             * Retrieve a parameter's Parameter object by name.
             *
             * @param name  The name of the Parameter that shall be retrieved.
             */
            Parameter operator[] (const std::string & name) const;
            ///@}

            /*!
             * Compare two instances of Parameters on inequality of their
             * underlying implementations.
             *
             * @param rhs   The right hand side of the binary != operator.
             */
            bool operator!= (const Parameters & rhs) const;
    };

    /*!
     * Parameter is the class that holds all information of one of Parameters' parameters.
     */
    class Parameter
    {
        private:
            struct Data;
            struct Template;

            ///@name Internal Data
            ///@{
            std::shared_ptr<Parameters::Data> _parameters_data;

            unsigned _index;
            ///@}

            ///@name Basic Functions
            ///@{
            Parameter(const std::shared_ptr<Parameters::Data> & imp, unsigned index);
            ///@}

        public:
            friend class Parameters;
            friend class Implementation<Parameters>;

            ///@name Basic Functions
            ///@{
            ~Parameter();
            ///@}

            ///@name Access & Modification of the Numeric Value
            ///@{
            /// Cast a Parameter's numeric value to a double.
            operator double () const;

            /// Retrieve a Parameter's numeric value.
            double operator() () const;

            /// Set a Parameter's numeric value.
            const Parameter & operator= (const double &);

            /*!
             * Retrieve a sample of the Parameter's numeric value based on its intrinsic probability distribution.
             *
             * @param engine The random number engine that shall be used to obtain the sample
             */
            double sample(RandomNumberEngine & engine) const;
            ///@}

            ///@name Access to Meta Data
            ///@{
            /// Retrieve the Parameter's name.
            const std::string & name() const;

            /// Retrieve the Parameter's (default) central value.
            const double & central() const;

            /// Retrieve the Parameter's (default) maximal value.
            const double & max() const;

            /// Retrieve the Parameter's (default) minimal value.
            const double & min() const;
            ///@}
    };

    struct ParameterDescription
    {
        Parameter parameter;

        double min;

        double max;

        bool nuisance;
    };

    struct ParameterRange
    {
            double min, max;
    };
}

#endif
