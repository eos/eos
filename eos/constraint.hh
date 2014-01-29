/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2012, 2013, 2014 Danny van Dyk
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


#ifndef EOS_GUARD_SRC_CONSTRAINT_HH
#define EOS_GUARD_SRC_CONSTRAINT_HH 1

#include <eos/observable.hh>
#include <eos/utils/log_likelihood-fwd.hh>
#include <eos/utils/observable_cache.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <vector>

namespace eos
{
    // Forward declaration.
    class Constraint;

    /*!
     * Constraint models experimental constraints via one or more LogLikelihoodBlock objects that
     * depend on one more Observable objects.
     */
    class Constraint :
        public PrivateImplementationPattern<Constraint>
    {
        public:
            /// Constructor.
            Constraint(const std::string & name, const std::vector<ObservablePtr> & observable,
                    const std::vector<LogLikelihoodBlockPtr> & blocks);

            /// Destructor.
            ~Constraint();

            ///@name Access
            ///@{
            std::string name() const;
            ///@}

            ///@name Iteration over Observables
            ///@{
            struct ObservableIteratorTag;
            typedef WrappedForwardIterator<ObservableIteratorTag, ObservablePtr> ObservableIterator;

            ObservableIterator begin_observables() const;
            ObservableIterator end_observables() const;
            ///@}

            ///@name Iteration over Blocks
            ///@{
            struct BlockIteratorTag;
            typedef WrappedForwardIterator<BlockIteratorTag, LogLikelihoodBlockPtr> BlockIterator;

            BlockIterator begin_blocks() const;
            BlockIterator end_blocks() const;
            ///@}

            /*!
             * Create one of the builtin constraints by name.
             *
             * Constraint names are structured as PROCESS::NAME@EXPERIMENT-YEAR
             *
             * @note You need to clone each LogLikelihoodBlock of a Constraint
             * to actually be able to use it.
             */
            static Constraint make(const std::string & name, const Options & options);
    };

    /*!
     * UnknownConstraintError is thrown when Constraint::make encounters an unknown constraint name.
     */
    struct UnknownConstraintError :
        public Exception
    {
        ///@name Basic Functions
        ///@{
        /*!
         * Constructor.
         *
         * @param name The offending constraint name.
         */
        UnknownConstraintError(const std::string & name);
        ///@}
    };
}

#endif
