/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_GAMMA_HH
#define WFITTER_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_GAMMA_HH 1

#include <src/rare-b-decays/decays.hh>
#include <src/utils/complex.hh>
#include <src/utils/observable.hh>
#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace wf
{
    /*
     * Decay: B -> K^* gamma
     */


    class BToKstarGamma :
        public PrivateImplementationPattern<BToKstarGamma>
    {
        public:
            BToKstarGamma(const Parameters & parameters, const ObservableOptions & options);
            ~BToKstarGamma();

            // Time dependent CP asymmetry S_K^*gamma
            double s_kstar_gamma() const;
    };
}

#endif
