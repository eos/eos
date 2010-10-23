/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_CHARM_LOOPS_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_CHARM_LOOPS_HH 1

#include <src/utils/complex.hh>

namespace eos
{
    struct CharmLoops
    {
        /* One-loop functions */
        // cf. [BFS2001], Eq. (11), p. 4 with m_q -> 0
        static complex<double> h(const double & mu, const double & s);
        // cf. [BFS2001], Eq. (11), p. 4 with m_q -> 0
        static complex<double> h(const double & mu, const double & s, const double & m_q);

        /* Two-loop functions */
        // cf. [S2004], Eq. (29), p. 8
        static complex<double> A(const double & mu, const double & s, const double & m_b);
        // cf. [S2004], Eq. (30), pp. 8-9
        static complex<double> B(const double & mu, const double & s, const double & m_b);
        // cf. [S2004], Eq. (31), p. 9
        static complex<double> C(const double & mu, const double & s);

        /* Non-factorizing two loop contributions */
        // massless case, cf. [S2004], Eq. (16), p. 6
        static complex<double> F17_massless(const double & mu, const double & s, const double & m_b);
        static complex<double> F19_massless(const double & mu, const double & s, const double & m_b);
        static complex<double> F27_massless(const double & mu, const double & s, const double & m_b);
        static complex<double> F29_massless(const double & mu, const double & s, const double & m_b);
        // massless case, cf. [BFS2001], Eqs. (82)-(83), p. 30
        static complex<double> F87_massless(const double & mu, const double & s, const double & m_b);
        static complex<double> F89_massless(                   const double & s, const double & m_b);

        // massive case, cf. [ABGW2003], Eq. (7), p. 8
        static complex<double> F17_massive(const double & mu, const double & s, const double & m_b, const double & m_c);
        static complex<double> F19_massive(const double & mu, const double & s, const double & m_b, const double & m_c);
        static complex<double> F27_massive(const double & mu, const double & s, const double & m_b, const double & m_c);
        static complex<double> F29_massive(const double & mu, const double & s, const double & m_b, const double & m_c);

        // helper functions for F8j, cf. [BFS2001], Eqs. (29) and (84), pp. 8 and 30
        static complex<double> B0(const double & s, const double & m_q);
        static complex<double> C0(const double & s, const double & m_q);
    };
}

#endif
