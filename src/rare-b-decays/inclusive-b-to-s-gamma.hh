/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_INCLUSIVE_B_TO_S_GAMMA_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_INCLUSIVE_B_TO_S_GAMMA_HH 1

#include <src/rare-b-decays/decays.hh>
#include <src/utils/complex.hh>
#include <src/utils/observable.hh>
#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace eos
{
    /*
     * Decay: b -> X_s gamma
     */

    // LO-NP Enhancement of [2006zs]
    struct Minimal
    {
    };

    template <>
    class BToXsGamma<Minimal> :
        public PrivateImplementationPattern<BToXsGamma<Minimal>>
    {
        public:
            BToXsGamma(const Parameters & parameters, const ObservableOptions & options);
            ~BToXsGamma();

            // Integrated Observables
            double integrated_branching_ratio() const;
    };
}
#endif
