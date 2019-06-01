/* vim: set sw=4 sts=4 et tw=150 foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2016-2019 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_OBSERVABLE_HH
#define EOS_GUARD_SRC_UTILS_OBSERVABLE_HH 1

#include <eos/observable-fwd.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>

#include <string>

namespace eos
{
    /**
     * Observable is internally used to handle the creation,
     * evaluation and cloning of any (pseudo)observable quantities.
     */
    class Observable :
        public ParameterUser
    {
        public:
            virtual const QualifiedName & name() const = 0;

            virtual double evaluate() const = 0;

            virtual Kinematics kinematics() = 0;

            virtual Parameters parameters() = 0;

            virtual Options options() = 0;

            virtual ObservablePtr clone() const = 0;

            virtual ObservablePtr clone(const Parameters & parameters) const = 0;

            static ObservablePtr make(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options);
    };

    /**
     * ObservableSection is used to keep track of one or more ObservableGroup objects, and groups
     * them together under a common name. Examples of observable sections include semileptonic B decays and form factors.
     */
    class ObservableSection :
        public PrivateImplementationPattern<ObservableSection>
    {
        public:
            ObservableSection(Implementation<ObservableSection> *);

            ~ObservableSection();

            ///@name Iteration over groups
            ///@{
            struct GroupIteratorTag;
            typedef WrappedForwardIterator<GroupIteratorTag, const ObservableGroup &> GroupIterator;

            GroupIterator begin() const;
            GroupIterator end() const;
            ///@}

            ///@name Meta data
            ///@{
            const std::string & name() const;
            const std::string & description() const;
            ///@}
    };

    extern template class WrappedForwardIterator<ObservableSection::GroupIteratorTag, const ObservableGroup &>;

    /**
     * ObservableGroup is used to keep track of one or more ObservableEntry objects, and groups
     * them together under a common name and description. Examples of Observables Groups include B->pilnu observables and B->D form factors.
     */
    class ObservableGroup :
        public PrivateImplementationPattern<ObservableGroup>
    {
        public:
            ObservableGroup(Implementation<ObservableGroup> *);

            ~ObservableGroup();

            ///@name Iteration over observables
            ///@{
            struct ObservableIteratorTag;
            typedef WrappedForwardIterator<ObservableIteratorTag, const std::pair<const QualifiedName, ObservableEntryPtr>> ObservableIterator;

            ObservableIterator begin() const;
            ObservableIterator end() const;
            ///@}

            ///@name Meta data
            ///@{
            const std::string & name() const;
            const std::string & description() const;
            ///@}
    };
    extern template class WrappedForwardIterator<ObservableGroup::ObservableIteratorTag, const std::pair<const QualifiedName, ObservableEntryPtr>>;

    /**
     * ObservableEntry is internally used to keep track of the description and factory method
     * for any given Observable. This includes handling its construction (via the make() method), and
     * describing it (via the ostream & insert() method).
     */
    class ObservableEntry
    {
        public:
            friend std::ostream & operator<< (std::ostream &, const ObservableEntry &);

            ObservableEntry();

            virtual ~ObservableEntry();

            virtual ObservablePtr make(const Parameters &, const Kinematics &, const Options &) const = 0;

            virtual const QualifiedName & name() const = 0;

            virtual const std::string & latex() const = 0;

            ///@name Iteration over kinematic variables
            ///@{
            struct KinematicVariableIteratorTag;
            using KinematicVariableIterator = WrappedForwardIterator<KinematicVariableIteratorTag, const std::string &>;

            virtual KinematicVariableIterator begin_kinematic_variables() const = 0;
            virtual KinematicVariableIterator end_kinematic_variables() const = 0;
            ///@}

        protected:
            virtual std::ostream & insert(std::ostream & os) const = 0;
    };

    /*!
     * Output stream operator for ObservableEntry.
     */
    inline std::ostream & operator<< (std::ostream & os, const ObservableEntry & entry)
    {
        return entry.insert(os);
    }

    /*!
     * Container around the known and implemented signal PDFs
     */
    class Observables :
        public PrivateImplementationPattern<Observables>
    {
        public:
            /// Constructor.
            Observables();

            /// Destructor.
            ~Observables();

            ///@name Access of individual ObserableEntry instances
            ///@{
            ObservableEntryPtr operator[] (const QualifiedName &) const;
            ///@}

            ///@name Iteration over observables
            ///@{
            struct ObservableIteratorTag;
            typedef WrappedForwardIterator<ObservableIteratorTag, const std::pair<const QualifiedName, ObservableEntryPtr>> ObservableIterator;

            ObservableIterator begin() const;
            ObservableIterator end() const;
            ///@}

            ///@name Iteration over groups of observables
            ///@{
            struct SectionIteratorTag;
            typedef WrappedForwardIterator<SectionIteratorTag, const ObservableSection &> SectionIterator;

            SectionIterator begin_sections() const;
            SectionIterator end_sections() const;
            ///@}
    };

    extern template class WrappedForwardIterator<Observables::ObservableIteratorTag, const std::pair<const QualifiedName, ObservableEntryPtr>>;
    extern template class WrappedForwardIterator<Observables::SectionIteratorTag, const ObservableSection &>;
}

#endif
