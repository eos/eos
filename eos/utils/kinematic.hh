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

#ifndef EOS_GUARD_SRC_UTILS_KINEMATIC_HH
#define EOS_GUARD_SRC_UTILS_KINEMATIC_HH 1

#include <eos/utils/exception.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*!
     * UnknownKinematicVariableError is thrown when no parameter of a given
     * name could be found.
     */
    struct UnknownKinematicVariableError :
        public Exception
    {
        UnknownKinematicVariableError(const std::string & variable) throw ();
    };

    // Forward declaration.
    class KinematicVariable;

    /*!
     * Kinematics keeps the set of all kinematic variables for any Observable.
     *
     * Access to any KinematicVariable or their values is coherent, i.e., changes to
     * a KinematicVariable object will propagate to every other object with the same
     * parent Kinematics and which handle the same variable by name.
     */
    class Kinematics :
        public PrivateImplementationPattern<Kinematics>
    {
        private:
            ///@name Internal Data
            ///@{
            struct Data;
            ///@}

        public:
            ///@name Basic Functions
            ///@{
            /// Constructor.
            Kinematics();

            /*!
             * Constructor.
             *
             * Create an instance of Kinematics with a given set of initial kinematic
             * variables.
             *
             * @param variables The set of initial kinematics variables from which this object shall be constructed.
             */
            Kinematics(const std::initializer_list<std::pair<std::string, double>> & variables);

            /// Destructor.
            ~Kinematics();

            /*!
             * Create an independent copy of this Kinematics object.
             */
            Kinematics clone() const;

            /// Equality comparison operator.
            bool operator== (const Kinematics & rhs) const;

            /// Inequality comparison operator.
            bool operator!= (const Kinematics & rhs) const;
            ///@}

            ///@name Variable access
            ///@{
            /*!
             * Declare a new kinematic variable.
             *
             * @param name  Name of the new variable to be declared.
             * @param value (Optional) value for the new variable.
             */
            KinematicVariable declare(const std::string & name, const double & value = 0.0);

            /*!
             * Set a kinematic variable's numeric value.
             *
             * @param name  The name of the variable whose numeric value shall be changed.
             * @param value The variable's new numeric value.
             */
            void set(const std::string & variable, const double & value);

            /*!
             * Retrieve a variable's KinematicVariable object by name.
             *
             * @param name  The name of the KinematicVariable that shall be retrieved.
             */
            KinematicVariable operator[] (const std::string & variable) const;
            ///@}

            ///@name Output
            ///@{
            /*!
             * Retrieve a string representation of the set of kinematic
             * variables.
             */
            std::string as_string() const;
            ///@}
    };


    /*!
     * KinematicVariable is the class that holds all information of one of KinematicVariables' parameters.
     */
    class KinematicVariable
    {
        private:
            ///@name Internal Data
            ///@{
            std::shared_ptr<Implementation<Kinematics>> _imp;

            unsigned _index;
            ///@}

            ///@name Basic Functions
            ///@{
            KinematicVariable(const std::shared_ptr<Implementation<Kinematics>> & imp, unsigned index);
            ///@}

        public:
            friend class Kinematics;
            friend struct Implementation<Kinematics>;

            ///@name Basic Functions
            ///@{
            ~KinematicVariable();
            ///@}

            ///@name Access & Modification of the Numeric Value
            ///@{
            /// Cast a KinematicVariable's numeric value to a double.
            operator double () const;

            /// Retrieve a KinematicVariable's numeric value.
            double operator() () const;

            /// Set a KinematicVariable's numeric value.
            const KinematicVariable & operator= (const double &);
            ///@}
    };

    inline double lambda(const double & a, const double & b, const double & c)
    {
        return a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c);
    }
}

#endif
