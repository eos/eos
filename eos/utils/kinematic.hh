/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2015, 2018 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_KINEMATIC_HH
#define EOS_GUARD_EOS_UTILS_KINEMATIC_HH 1

#include <eos/utils/exception.hh>
#include <eos/utils/mutable.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <set>

namespace eos
{
    /*!
     * UnknownKinematicVariableError is thrown when no parameter of a given
     * name could be found.
     */
    struct UnknownKinematicVariableError : public Exception
    {
            UnknownKinematicVariableError(const std::string & variable) throw();
    };

    /*!
     * DuplicateKinematicAliasError is thrown when an alias is defined twice.
     */
    struct DuplicateKinematicAliasError : public Exception
    {
            DuplicateKinematicAliasError(const std::string & alias, const std::string & variable) throw();
    };

    /*!
     * UnknownKinematicAliasError is thrown when no alias of a given
     * name could be found.
     */
    struct UnknownKinematicAliasError : public Exception
    {
            UnknownKinematicAliasError(const std::string & alias) throw();
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
    class Kinematics : public PrivateImplementationPattern<Kinematics>
    {
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

            Kinematics operator+ (const Kinematics & rhs) const;

            /// Equality comparison operator.
            bool operator== (const Kinematics & rhs) const;

            /// Inequality comparison operator.
            bool operator!= (const Kinematics & rhs) const;
            ///@}

            ///@name Variable access
            ///@{

            /*!
             * Create an alias of an existing kinematic variable, under a new name.
             *
             * @param alias Alternative name for the existing variable; must not exist!
             * @param name  Name of the existing variable; must exist!
             */
            void alias(const std::string & alias, const std::string & name);

            /*!
             * Remove an existing alias of a kinematic variable.
             *
             * @param alias Name of an existing alias that shall be removed.
             */
            void remove_alias(const std::string & alias);

            /*!
             * Reset all defined aliases.
             */
            void clear_aliases();

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

            ///@name Iteration over our kinematic variables
            ///@{
            struct KinematicVariableIteratorTag;
            using KinematicVariableIterator = WrappedForwardIterator<KinematicVariableIteratorTag, const KinematicVariable>;

            KinematicVariableIterator begin() const;
            KinematicVariableIterator end() const;
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

    extern template class WrappedForwardIterator<Kinematics::KinematicVariableIteratorTag, const KinematicVariable>;

    /*!
     * KinematicVariable is the class that holds all information of one of KinematicVariables' parameters.
     */
    class KinematicVariable : public Mutable
    {
        private:
            ///@name Internal Data
            ///@{
            std::shared_ptr<Implementation<Kinematics>> _imp;

            unsigned _index;

            bool _is_alias;
            ///@}

            ///@name Basic Functions
            ///@{
            KinematicVariable(const std::shared_ptr<Implementation<Kinematics>> & imp, unsigned index, bool is_alias);
            ///@}

        public:
            friend class Kinematics;
            friend struct Implementation<Kinematics>;

            /*!
             * A unique number that identifies this parameter at run time.
             */
            using Id = unsigned;

            ///@name Basic Functions
            ///@{
            ~KinematicVariable();

            /// Make a copy of this KinematicVariable as a MutablePtr.
            MutablePtr clone() const;
            ///@}

            ///@name Access & Modification of the Numeric Value
            ///@{
            /// Cast a KinematicVariable's numeric value to a double.
            operator double () const;

            /// Retrieve a KinematicVariable's numeric value.
            double operator() () const;

            /// Retrieve a KinematicVariable's numeric value.
            virtual double evaluate() const;

            /// Set a KinematicVariable's numeric value.
            const KinematicVariable & operator= (const double &);

            /// Retrieve the Parameter's id.
            Id id() const;

            /// Set a KinematicVariable's numeric value.
            virtual void set(const double &);
            ///@}

            ///@name Access to Meta Data
            ///@{
            /// Retrieve the Parameter's name.
            virtual const std::string & name() const;
            ///@}
    };

    /*!
     * Base class for all users of Kinematics objects.
     */
    class KinematicUser
    {
        protected:
            std::set<KinematicVariable::Id> _ids;

        public:
            ///@name Iteration over ids
            ///@{
            struct ConstIteratorTag;
            using ConstIterator = WrappedForwardIterator<ConstIteratorTag, const KinematicVariable::Id>;

            ConstIterator begin_kinematics() const;
            ConstIterator end_kinematics() const;
            ///@}

            ///@name Access
            ///@{
            /*!
             * Remove a kinematics variable from the set of used ids.
             */
            void drop(const KinematicVariable::Id & id);

            /*!
             * Add a given kinematics variable id to our list of used ids.
             *
             * @param id   The kinematics variable id that we use.
             */
            void uses_kinematic(const KinematicVariable::Id & id);

            /*!
             * Copy parameter ids of another KinematicUser to our list of used ids.
             *
             * @param user The other ParameterUser whose ids we are going to copy.
             */
            void uses_kinematic(const KinematicUser & user);
            ///@}
    };

    extern template class WrappedForwardIterator<KinematicUser::ConstIteratorTag, const KinematicVariable::Id>;

    /*!
     * Wrapper class to automate usage tracking of Kinematics objects.
     */
    class UsedKinematicVariable : public KinematicVariable
    {
        public:
            /*!
             * Constructor.
             *
             * Constructs a KinematicVariable object and registers its usage with a KinematicUser.
             *
             * @param variable The kinematics variable which is used.
             * @param user      The user of above kinematics variable.
             */
            UsedKinematicVariable(const KinematicVariable & variable, KinematicUser & user);

            using KinematicVariable::operator=;
    };

    template <typename T>
    inline T
    lambda(const T & a, const T & b, const T & c)
    {
        return a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c);
    }
} // namespace eos

#endif
