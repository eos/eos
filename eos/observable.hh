/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2016, 2017 Danny van Dyk
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

#include <eos/utils/exception.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>

#include <string>
#include <memory>

namespace eos
{
    class Observable;

    class Options;

    typedef std::shared_ptr<Observable> ObservablePtr;

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

            ///@name Iteration over known constraints
            ///@{
            struct ObservableIteratorTag;
            typedef WrappedForwardIterator<ObservableIteratorTag, const std::pair<const QualifiedName, const ObservableEntry *>> ObservableIterator;

            ObservableIterator begin() const;
            ObservableIterator end() const;
            ///@}
    };

    extern template class WrappedForwardIterator<Observables::ObservableIteratorTag, const std::pair<const QualifiedName, const ObservableEntry *>>;
}

#endif
