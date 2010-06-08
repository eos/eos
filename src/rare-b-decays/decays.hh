/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_RARE_B_DECAYS_DECAYS_HH
#define WFITTER_GUARD_SRC_RARE_B_DECAYS_DECAYS_HH 1

namespace wf
{
    template <typename T_> class BToKstarDilepton;

    template <typename T_> class BToXsDilepton;

    template <typename T_> class BToXsGamma;

    enum Helicity
    {
        left_handed = -1,
        right_handed = +1
    };
}

#endif
