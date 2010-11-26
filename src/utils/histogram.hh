/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_HISTROGRAM_HH
#define EOS_GUARD_SRC_UTILS_HISTROGRAM_HH 1

#include <src/utils/private_implementation_pattern.hh>
#include <src/utils/wrapped_forward_iterator.hh>

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
}

#endif
