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

#ifndef EOS_GUARD_SRC_UTILS_HISTROGRAM_HH
#define EOS_GUARD_SRC_UTILS_HISTROGRAM_HH 1

#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <array>

namespace eos
{
    template <std::size_t dimensions_> class Histogram;

    template <> class Histogram<1> :
        public PrivateImplementationPattern<Histogram<1>>
    {
        public:
            struct Bin;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * Creates an empty Histogram, i.e. a Histogram without any bins.
             */
            Histogram();

            /// Destructor.
            ~Histogram();

            /*!
             * Named constructor.
             *
             * Creates a Histogram with pre-existing bins of identical width.
             * @param start    Left-most value that shall be covered by the Histogram.
             * @param end      Right-most value that shall be covered by the Histogram.
             * @param count    Number of Bins in the Histogram.
             */
            static Histogram WithEqualBinning(const double & start, const double & end, const unsigned & count);
            ///@}

            ///@name Insertion and Retrieval
            ///@{
            /// Inserts a bin into the histogram.
            void insert(const Bin & bin);
            /// Inserts a value into the histogram's matching bin.
            void insert(const double & value);
            /// Returns the number of entries in the histogram.
            unsigned entries() const;
            ///@}

            ///@name Iteration
            ///@{
            struct ConstIteratorTag;
            typedef WrappedForwardIterator<ConstIteratorTag, const Bin> ConstIterator;

            ConstIterator begin() const;
            ConstIterator end() const;
            ///@}
    };

    /*!
     * Holds information on any one of Histogram<1>'s bins.
     */
    struct Histogram<1>::Bin
    {
        /// Lower limit of the Bin's interval.
        double lower;

        /// Width of the Bin's interval.
        double upper;

        /// Content of the Bin.
        double value;

        ///@name Basic Functions
        ///@{
        /*
         * Constructor.
         *
         * @param lower    Left-most (inclusive) border of the Bin's interval.
         * @param upper    Right-most (exclusive) border of the Bin's interval.
         * @param value    Initial value of the Bin.
         */
        Bin(const double & lower, const double & upper, const double & value = 0.0) :
            lower(lower),
            upper(upper),
            value(value)
        {
        }

        /// Destructor.
        ~Bin() { }
        ///@}

        /*!
         * Compare if this bin is left of another bin.
         *
         * @param other    The other Bin.
         */
        bool operator< (const Bin & other) const
        {
            return this->lower < other.lower;
        }
    };

    /*!
     * Compute the estimated cumulative distribution function for a given 1D distribution.
     *
     * @param distribution   The histogram of a distribution function.
     */
    Histogram<1> estimate_cumulative_distribution(const Histogram<1> & distribution);

    template <> class Histogram<2> :
        public PrivateImplementationPattern<Histogram<2>>
    {
        public:
            struct Bin;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * Creates an empty Histogram, i.e. a Histogram without any bins.
             */
            Histogram();

            /// Destructor.
            ~Histogram();

            /*!
             * Named constructor.
             *
             * Creates a Histogram with pre-existing bins of identical width.
             * @param start    Left-most value for each dimension that shall be covered by the Histogram.
             * @param end      Right-most value for each dimension that shall be covered by the Histogram.
             * @param count    Number of Bins for each dimension.
             */
            static Histogram WithEqualBinning(const std::array<double, 2> & start, const std::array<double, 2> & end, const std::array<unsigned, 2> & count);
            ///@}

            ///@name Insertion and Retrieval
            ///@{
            /// Inserts a bin into the histogram.
            void insert(const Bin & bin);

            /// Returns the number of entries in the histogram.
            unsigned entries() const;
            ///@}

            ///@name Iteration
            ///@{
            struct IteratorTag;
            typedef WrappedForwardIterator<IteratorTag, Bin> Iterator;

            Iterator begin();
            Iterator find(const std::array<double, 2> & coordinates);
            Iterator end();
            ///@}

            ///@name Constant Iteration
            ///@{
            struct ConstIteratorTag;
            typedef WrappedForwardIterator<ConstIteratorTag, const Bin> ConstIterator;

            ConstIterator cbegin() const;
            ConstIterator cend() const;
            ///@}
    };

    /*!
     * Holds information on any one of Histogram<2>'s bins.
     */
    struct Histogram<2>::Bin
    {
        /// Lower limit of the Bin's interval for each dimension.
        std::array<double, 2> lower;

        /// Upper limit of the Bin's interval for each dimension.
        std::array<double, 2> upper;

        /// Content of the Bin.
        double value;

        ///@name Basic Functions
        ///@{
        /*
         * Constructor.
         *
         * @param lower    Left-most (inclusive) border of the Bin's interval.
         * @param upper    Right-most (exclusive) border of the Bin's interval.
         * @param value    Initial value of the Bin.
         */
        Bin(const std::array<double, 2> & lower, const std::array<double, 2> & upper, const double & value = 0.0) :
            lower(lower),
            upper(upper),
            value(value)
        {
        }

        /// Destructor.
        ~Bin() { }
        ///@}

        /*!
         * Compare if this bin is of lower order than another bin.
         *
         * @param other    The other Bin.
         */
        bool operator< (const Bin & other) const
        {
            if (this->lower[0] < other.lower[0])
                return true;

            if (this->lower[0] > other.lower[0])
                return false;

            if (this->lower[1] < other.lower[1])
                return true;

            return false;
        }
    };
}

#endif
