#ifndef EOS_GUARD_EOS_STATISTICS_DENSITY_WRAPPER_HH
#define EOS_GUARD_EOS_STATISTICS_DENSITY_WRAPPER_HH

#include <eos/statistics/density.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <vector>

namespace eos
{
   /*!
     * A wrapper around a multivariate scalar function
     */
    class DensityWrapper :
        public Density
    {
        public:
            typedef double (* RawDensity)(const std::vector<double> &);
            typedef std::function< double (const std::vector<double> &)> WrappedDensity;

            /// todo document
            DensityWrapper(const WrappedDensity &);
            DensityWrapper(RawDensity);

            virtual ~DensityWrapper();

            void add(const ParameterDescription & def);

            /// Evaluate the density function at the current parameter point on the log scale.
            virtual double evaluate() const;

            /// Create an independent copy of this density function.
            virtual DensityPtr clone() const;

            /// Iterate over the parameters relevant to this density function.
            ///@{
            virtual Iterator begin() const;
            virtual Iterator end() const;
            ///@}
        private:
            WrappedDensity _density;
            std::vector<ParameterDescription> _defs;
            mutable std::vector<double> _parameter_values;
    };
}
#endif
