/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_RARE_B_DECAYS_INCLUSIVE_B_TO_S_DILEPTON_HH
#define WFITTER_GUARD_SRC_RARE_B_DECAYS_INCLUSIVE_B_TO_S_DILEPTON_HH 1

#include <src/rare-b-decays/decays.hh>
#include <src/utils/complex.hh>
#include <src/utils/observable.hh>
#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace wf
{
    /*
     * Decay: B -> X_s l lbar
     */


    // As given in [HLMW2005]
    struct HLMW2005
    {
    };

    template <>
    class BToXsDilepton<HLMW2005> :
        public PrivateImplementationPattern<BToXsDilepton<HLMW2005>>
    {
        public:
            BToXsDilepton(const Parameters & parameters, const ObservableOptions & options);
            ~BToXsDilepton();

            // Differential Observables
            double differential_branching_ratio(const double & s) const;

            // Integrated Observables
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
    };
}

#endif
