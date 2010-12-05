/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WILSON_FITTER_GUARD_SRC_UTILS_RANDOM_NUMBER_ENGINE_HH
#define WILSON_FITTER_GUARD_SRC_UTILS_RANDOM_NUMBER_ENGINE_HH 1

#include <src/utils/private_implementation_pattern.hh>

#include <cstdint>

namespace eos
{
    /*!
     * Generate pseudo-random floating point numbers
     * in the range 0 ... 1.
     *
     * RandomNumberEngine keeps its state/seed on copying.
     */
    class RandomNumberEngine :
        public PrivateImplementationPattern<RandomNumberEngine>
    {
        public:
            ///@{
            /// Basic functions.
            RandomNumberEngine();
            ~RandomNumberEngine();
            ///@}

            /// Obtain a pseudo-random number in the range 0.0 .. 1.0.
            double operator() ();

            /// Return the maximal value that can be drawn.
            double max() const { return 1.0; }

            /// Return the minimal value that can be drawn.
            double min() const { return 0.0; }
    };
}

#endif
