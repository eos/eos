/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_RARE_B_DECAYS_HARD_SCATTERING_HH
#define WFITTER_GUARD_SRC_RARE_B_DECAYS_HARD_SCATTERING_HH 1

#include <complex>

namespace wf
{
    using std::complex;

    struct HardScattering
    {
        // cf. [BFS2001], Eqs. (30)-(32), p. 8
        // @s    : Dilepton invariant mass
        // @u    : Relative contribution of the quark (versus ubar = 1-u for the
        //         antiquark) to the light meson's energy.
        // @m_q  : Mass of the internal loop quark
        // @m_B  : Mass of the B meson
        static complex<double> I1(const double & s, const double & u, const double & m_q, const double & m_B);
    };
}

#endif
