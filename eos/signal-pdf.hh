/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2023 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_SIGNAL_PDF_HH
#define EOS_GUARD_EOS_UTILS_SIGNAL_PDF_HH 1

#include <eos/signal-pdf-fwd.hh>
#include <eos/utils/density.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/iterator-range.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>

#include <map>
#include <memory>
#include <string>

namespace eos
{
    class KinematicRange
    {
        public:
            const char * name;

            const double min, max;

            const std::string description;

            operator const char * () const { return name; }
    };

    class SignalPDF : public Density
    {
        public:
            virtual const QualifiedName & name() const = 0;

            virtual double evaluate() const = 0;

            virtual double normalization() const = 0;

            virtual Kinematics kinematics() = 0;

            virtual Parameters parameters() = 0;

            virtual Options options() = 0;

            virtual DensityPtr clone() const = 0;

            virtual DensityPtr clone(const Parameters & parameters) const = 0;

            static SignalPDFPtr make(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options);
    };

    /**
     * SignalPDFSection is used to keep track of one or more SignalPDFGroup objects, and groups
     * them together under a common name. Examples of observable sections include semileptonic B decays
     * or ee->hadrons.
     */
    class SignalPDFSection : public PrivateImplementationPattern<SignalPDFSection>
    {
        public:
            SignalPDFSection(Implementation<SignalPDFSection> *);

            ~SignalPDFSection();

            ///@name Iteration over groups
            ///@{
            struct GroupIteratorTag;
            using GroupIterator = WrappedForwardIterator<GroupIteratorTag, const SignalPDFGroup &>;

            GroupIterator begin() const;
            GroupIterator end() const;
            ///@}

            ///@name Meta data
            ///@{
            const std::string & name() const;
            const std::string & description() const;
            ///@}
    };

    extern template class WrappedForwardIterator<SignalPDFSection::GroupIteratorTag, const SignalPDFGroup &>;

    /**
     * SignalPDFGroup is used to keep track of one or more ObservableEntry objects, and groups
     * them together under a common name and description. Examples of Observables Groups include B->pilnu observables and B->D form factors.
     */
    class SignalPDFGroup : public PrivateImplementationPattern<SignalPDFGroup>
    {
        public:
            SignalPDFGroup(Implementation<SignalPDFGroup> *);

            ~SignalPDFGroup();

            ///@name Iteration over observables
            ///@{
            struct SignalPDFIteratorTag;
            using SignalPDFIterator = WrappedForwardIterator<SignalPDFIteratorTag, const std::pair<const QualifiedName, SignalPDFEntryPtr>>;

            SignalPDFIterator begin() const;
            SignalPDFIterator end() const;
            ///@}

            ///@name Meta data
            ///@{
            const std::string & name() const;
            const std::string & description() const;
            ///@}
    };
    extern template class WrappedForwardIterator<SignalPDFGroup::SignalPDFIteratorTag, const std::pair<const QualifiedName, SignalPDFEntryPtr>>;

    /*!
     * SignalPDFEntry is internally used to keep track of the SignalPDFDescription and the SignalPDFFactory
     * for any given SignalPDF. This includes handling its construction (via the make() method), and
     * describing it (via the ostream & insert() method).
     */
    class SignalPDFEntry
    {
        public:
            friend std::ostream & operator<< (std::ostream &, const SignalPDFEntry &);

            SignalPDFEntry();

            virtual ~SignalPDFEntry();

            virtual SignalPDFPtr make(const Parameters &, const Kinematics &, const Options &) const = 0;

            /// Return the SignalPDF name
            virtual const QualifiedName & name() const = 0;

            /// Return the SignalPDF description
            virtual const std::string & description() const = 0;

            ///@name Iteration over our kinematic ranges
            ///@{
            struct KinematicRangeIteratorTag;
            using KinematicRangeIterator = WrappedForwardIterator<KinematicRangeIteratorTag, const KinematicRange>;

            virtual KinematicRangeIterator begin_kinematic_ranges() const = 0;
            virtual KinematicRangeIterator end_kinematic_ranges() const   = 0;

            inline IteratorRange<KinematicRangeIterator>
            kinematic_ranges() const
            {
                return IteratorRange<KinematicRangeIterator>(begin_kinematic_ranges(), end_kinematic_ranges());
            }

            ///@}

        protected:
            virtual std::ostream & insert(std::ostream & os) const = 0;
    };

    extern template class WrappedForwardIterator<SignalPDFEntry::KinematicRangeIteratorTag, const KinematicRange>;

    /*!
     * Output stream operator for SignalPDFEntry.
     */
    inline std::ostream &
    operator<< (std::ostream & os, const SignalPDFEntry & entry)
    {
        return entry.insert(os);
    }

    /*!
     * Container around the known and implemented signal PDFs
     */
    class SignalPDFs : public PrivateImplementationPattern<SignalPDFs>
    {
        public:
            /// Constructor.
            SignalPDFs();

            /// Destructor.
            ~SignalPDFs();

            ///@name Access of individual ObserableEntry instances
            ///@{
            SignalPDFEntryPtr operator[] (const QualifiedName &) const;
            ///@}

            ///@name Iteration over observables
            ///@{
            struct SignalPDFIteratorTag;
            using SignalPDFIterator = WrappedForwardIterator<SignalPDFIteratorTag, const std::pair<const QualifiedName, SignalPDFEntryPtr>>;

            SignalPDFIterator begin() const;
            SignalPDFIterator end() const;
            ///@}

            ///@name Iteration over groups of observables
            ///@{
            struct SectionIteratorTag;
            using SectionIterator = WrappedForwardIterator<SectionIteratorTag, const SignalPDFSection &>;

            SectionIterator begin_sections() const;
            SectionIterator end_sections() const;
            ///@}
    };

    extern template class WrappedForwardIterator<SignalPDFs::SignalPDFIteratorTag, const std::pair<const QualifiedName, SignalPDFEntryPtr>>;
    extern template class WrappedForwardIterator<SignalPDFs::SectionIteratorTag, const SignalPDFSection &>;

    class SignalPDFEntries : public InstantiationPolicy<SignalPDFEntries, Singleton>
    {
        private:
            std::map<QualifiedName, std::shared_ptr<const SignalPDFEntry>> * _entries;

            SignalPDFEntries();

            ~SignalPDFEntries();

        public:
            friend class InstantiationPolicy<SignalPDFEntries, Singleton>;

            inline const std::map<QualifiedName, std::shared_ptr<const SignalPDFEntry>> &
            entries() const
            {
                return *_entries;
            }
    };

    /*!
     * SignalPDFNameError is thrown when SignalPDF::make encounters a malformed observable name.
     */
    struct SignalPDFNameError : public Exception
    {
            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param name The offending malformed observable name.
             */
            SignalPDFNameError(const std::string & name);
            ///@}
    };
} // namespace eos

#endif
