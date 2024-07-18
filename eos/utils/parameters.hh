/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2024 Danny van Dyk
 * Copyright (c) 2021 Philip Lüghausen
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

#ifndef EOS_GUARD_EOS_UTILS_PARAMETERS_HH
#define EOS_GUARD_EOS_UTILS_PARAMETERS_HH 1

#include <eos/utils/exception.hh>
#include <eos/utils/mutable.hh>
#include <eos/utils/parameters-fwd.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/units.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <set>
#include <limits>

namespace eos
{
    // Forward declarations
    class ParameterDefaults;

    /*!
     * UnknownParameterError is thrown when no parameter of a given
     * name could be found.
     */
    struct UnknownParameterError :
        public Exception
    {
        UnknownParameterError(const QualifiedName & variable) throw ();
    };

    /*!
     * ParameterInputFileParseError is thrown when a malformed parameter input
     * file cannot be parsed by libyaml-cpp.
     */
    struct ParameterInputFileParseError :
        public Exception
    {
        ParameterInputFileParseError(const std::string & file, const std::string & msg) throw ();
    };

    /*!
     * ParameterInputFileNodeError is thrown when a malformed node is encountered
     * within a parameter input file.
     */
    struct ParameterInputFileNodeError :
        public Exception
    {
        ParameterInputFileNodeError(const std::string & file, const std::string & node, const std::string & msg) throw ();
    };

    /*!
     * ParameterInputDuplicateError is thrown when a duplicate parameter entry is encountered when parsing
     * several input files.
     */
    struct ParameterInputDuplicateError :
        public Exception
    {
        ParameterInputDuplicateError(const std::string & file, const std::string & msg) throw ();
    };

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
            friend class ParameterDefaults;
            friend struct Implementation<Parameter>;
            friend struct Implementation<Parameters>;

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

            ///@name Iteration over parameters
            ///@{
            struct IteratorTag;
            using Iterator = WrappedForwardIterator<IteratorTag, Parameter>;

            Iterator begin() const;
            Iterator end() const;
            ///@}

            ///@name Iteration over parameter sections
            ///@{
            struct SectionIteratorTag;
            using SectionIterator = WrappedForwardIterator<SectionIteratorTag, const ParameterSection &>;

            SectionIterator begin_sections() const;
            SectionIterator end_sections() const;
            ///@}

            ///@name Access to default parameters
            ///@{
            /*!
             * Declare a new default parameter.
             *
             * @param name  Name of the new parameter to be declared.
             * @param latex LaTeX representation of the new parameter.
             * @param unit  Unit of the new parameter.
             * @param value (Optional) value for the new parameter.
             * @param min   (Optional) minimal value for the new parameter.
             * @param max   (Optional) maximal value for the new parameter.
             */
            static unsigned declare(const QualifiedName & name, const std::string & latex, Unit unit,
                const double & value = 0.0,
                const double & min = -std::numeric_limits<double>::max(),
                const double & max = +std::numeric_limits<double>::max());

            /*!
             * Redirect a parameter name to a different parameter id in the default set of parameters.
             *
             * The internal mapping of the parameter name will be redirected to the new id.
             * If the the parameter's previous id is not already aliased, it will become inaccessible.
             * This is useful for example to alias a parameter name to a different parameter object.
             *
             * @param name  Name of the parameter to be redirected.
             * @param id    The id of the parameter to which the name shall be redirected.
             */
            static void redirect(const QualifiedName & name, const unsigned & id);
            ///@}

            ///@name Parameter access
            ///@{
            /*!
             * Declare a previously undeclared parameter in the default set of parameters and insert it
             * into this parameter set.
             *
             * @param name  Name of the new parameter to be declared.
             * @param latex LaTeX representation of the new parameter.
             * @param unit  Unit of the new parameter.
             * @param value (Optional) value for the new parameter.
             * @param min   (Optional) minimal value for the new parameter.
             * @param max   (Optional) maximal value for the new parameter.
             */
            Parameter declare_and_insert(const QualifiedName & name, const std::string & latex, Unit unit,
                const double & value = 0.0,
                const double & min = -std::numeric_limits<double>::max(),
                const double & max = +std::numeric_limits<double>::max());

            /*!
             * Redirect a parameter name to a different parameter id in the default set of parameters
             * and apply the redirection to this parameter set.
             *
             * The internal mapping of the parameter name will be redirected to the new id.
             * If the the parameter's previous id is not already aliased, it will become inaccessible.
             * This is useful for example to alias a parameter name to a different parameter object.
             *
             * @param name  Name of the parameter to be redirected.
             * @param id    The id of the parameter to which the name shall be redirected.
             */
            void redirect_and_apply(const QualifiedName & name, const unsigned & id);

            /*!
             * Set a parameter's numeric value.
             *
             * @param name  The name of the parameter whose numeric value shall be changed.
             * @param value The parameter's new numeric value.
             */
            void set(const QualifiedName & name, const double & value);

            /*!
             * Verify if a parameter with a given name exists.
             *
             * @param name  The name to be checked against the known parameters.
             */
            bool has(const QualifiedName & name);

            /*!
             * Retrieve a parameter's Parameter object by name.
             *
             * @param name  The name of the Parameter that shall be retrieved.
             */
            Parameter operator[] (const QualifiedName & name) const;

            /*!
             * Retrieve a parameter's Parameter object by id.
             *
             * @param id    The id of the Parameter that shall be retrived.
             */
            Parameter operator[] (const unsigned & id) const;

            /*!
             * Override the parameter values from an external YAML file.
             *
             * @param file  The name of the YAML fie.
             */
            void override_from_file(const std::string & file);
            ///@}

            /*!
             * Compare two instances of Parameters on inequality of their
             * underlying implementations.
             *
             * @param rhs   The right hand side of the binary != operator.
             */
            bool operator!= (const Parameters & rhs) const;
    };

    extern template class WrappedForwardIterator<Parameters::IteratorTag, Parameter>;

    /*!
     * Parameter is the class that holds all information of one of Parameters' parameters.
     */
    class Parameter :
        public Mutable
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
            friend class ParameterDefaults;
            friend struct Implementation<Parameters>;

            /*!
             * A unique number that identifies this parameter at run time.
             */
            using Id = unsigned;

            ///@name Basic Functions
            ///@{
            /// Copy-constructor.
            Parameter(const Parameter & other);

            /// Destructor.
            ~Parameter();

            /// Make a copy of this Parameter as a MutablePtr.
            MutablePtr clone() const;
            ///@}

            ///@name Access & Modification of the Numeric Value
            ///@{
            /// Cast a Parameter's numeric value to a double.
            virtual operator double () const;

            /// Retrieve a Parameter's numeric value.
            virtual double operator() () const;

            /// Retrieve a Parameter's numeric value.
            virtual double evaluate() const;

            /// Retrieve a Parameter's generator value, used for prior sampling.
            virtual double evaluate_generator() const;

            /// Set a Parameter's numeric value.
            virtual const Parameter & operator= (const double &);

            /// Set a Parameter's numeric value.
            virtual void set(const double &);

            /// Set a Parameter's generator value, used for prior sampling.
            virtual void set_generator(const double &);
            ///@}

            ///@name Access to Meta Data
            ///@{
            /// Retrieve the Parameter's name.
            virtual const std::string & name() const;

            /// Retrieve the Parameter's (default) central value.
            const double & central() const;

            /// Retrieve the Parameter's (default) maximal value.
            const double & max() const;

            /// Set the Parameter's maximal value.
            void set_max(const double &);

            /// Retrieve the Parameter's (default) minimal value.
            const double & min() const;

            /// Set the Parameter's maximal value.
            void set_min(const double &);

            /// Retrieve the Parameter's id.
            Id id() const;

            /// Retrieve the Parameter's name as a LaTeX representation
            const std::string & latex() const;

            /// Retrieve the Parameter's unit
            Unit unit() const;
            ///@}
    };

    /**
     * ParameterSection is used to keep track of one or more ParameterGroup objects, and groups
     * them together under a common name. Examples of observable sections include SM & EFT parameters,
     * and form factor parameters.
     */
    class ParameterSection :
        public PrivateImplementationPattern<ParameterSection>
    {
        public:
            ParameterSection(Implementation<ParameterSection> *);

            ~ParameterSection();

            ///@name Iteration over groups
            ///@{
            struct GroupIteratorTag;
            using GroupIterator = WrappedForwardIterator<GroupIteratorTag, const ParameterGroup &>;

            GroupIterator begin() const;
            GroupIterator end() const;
            ///@}

            ///@name Meta data
            ///@{
            const std::string & name() const;
            const std::string & description() const;
            ///@}
    };
    extern template class WrappedForwardIterator<ParameterSection::GroupIteratorTag, const ParameterGroup &>;

    /*!
     * ParameterGroup is used to keep track of one or more Parameter objects, and groups
     * them together under a common name and description. Examples of Parameter Groups
     * include fermion mass parameters and B->D form factors parameters.
     */
    class ParameterGroup :
        public PrivateImplementationPattern<ParameterGroup>
    {
        public:
            ParameterGroup(Implementation<ParameterGroup> *);

            ~ParameterGroup();

            ///@name Iteration over parameters
            ///@{
            struct ParameterIteratorTag;
            using ParameterIterator = WrappedForwardIterator<ParameterIteratorTag, const Parameter>;

            ParameterIterator begin() const;
            ParameterIterator end() const;
            ///@}

            ///@name Meta data
            ///@{
            const std::string & name() const;
            const std::string & description() const;
            ///@}
    };
    extern template class WrappedForwardIterator<ParameterGroup::ParameterIteratorTag, const Parameter>;

    /*!
     * Base class for all users of Parameter objects.
     */
    class ParameterUser
    {
        protected:
            std::set<Parameter::Id> _ids;

        public:
            ///@name Iteration over ids
            ///@{
            struct ConstIteratorTag;
            using ConstIterator = WrappedForwardIterator<ConstIteratorTag, const Parameter::Id>;

            ConstIterator begin() const;
            ConstIterator end() const;
            ///@}

            ///@name Access
            ///@{
            /*!
             * Remove a parameter from the set of used ids.
             */
            void drop(const Parameter::Id & id);

            /*!
             * Add a given parameter id to our list of used ids.
             *
             * @param id   The parameter id that we use.
             */
            void uses(const Parameter::Id & id);

            /*!
             * Copy parameter ids of another ParameterUser to our list of used ids.
             *
             * @param user The other ParameterUser whose ids we are going to copy.
             */
            void uses(const ParameterUser & user);
            ///@}
    };

    extern template class WrappedForwardIterator<ParameterUser::ConstIteratorTag, const Parameter::Id>;

    /*!
     * Wrapper class to automate usage tracking of Parameter objects.
     */
    class UsedParameter :
        public Parameter
    {
        public:
            /*!
             * Constructor.
             *
             * Constructs a Parameter object and registers its usage with a ParameterUser.
             *
             * @param parameter The parameter which is used.
             * @param user      The user of above parameter.
             */
            UsedParameter(const Parameter & parameter, ParameterUser & user);

            using Parameter::operator=;
    };

    struct ParameterDescription
    {
        MutablePtr parameter;

        double min;

        double max;

        bool nuisance;
    };

    bool operator== (const ParameterDescription & lhs, const ParameterDescription & rhs);
}

#endif
