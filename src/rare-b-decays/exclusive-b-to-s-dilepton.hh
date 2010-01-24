/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_HH
#define WFITTER_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_HH 1

#include <src/rare-b-decays/decays.hh>
#include <src/utils/complex.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace wf
{
    /*
     * Decay: B -> K l lbar
     */

    /*
     * Reference: [BHvD2010]
     */
    struct LowRecoil
    {
    };

    template <>
    class BToKstarDilepton<LowRecoil> :
        public PrivateImplementationPattern<BToKstarDilepton<LowRecoil>>
    {
        public:
            BToKstarDilepton(const double & mu);
            ~BToKstarDilepton();

            // [BHvD2010] Eqs. (??-??)
            Complex<double> a_long(const Helicity & h, const double & s);
            Complex<double> a_perp(const Helicity & h, const double & s);
            Complex<double> a_par(const Helicity & h, const double & s);
    };
}

#endif
