/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_LOW_RECOIL_LOW_RECOIL_HH
#define WFITTER_GUARD_SRC_LOW_RECOIL_LOW_RECOIL_HH 1

#include <src/utils/decay.hh>
#include <src/utils/complex.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace wf
{
    /*
     * Reference: arxiv:10XX.YYYY
     */

    struct BToKstarDilepton
    {
    };

    /*
     * @Decay: B -> K l lbar
     */
    template <>
    class Decay<BToKstarDilepton> :
        public PrivateImplementationPattern<Decay<BToKstarDilepton>>
    {
        public:
            Decay(const double & mu);
            ~Decay();

            // [BHvD2010] Eqs. (??-??)
            Complex<double> a_long(const Helicity & h, const double & s);
            Complex<double> a_perp(const Helicity & h, const double & s);
            Complex<double> a_par(const Helicity & h, const double & s);
    };
}

#endif
