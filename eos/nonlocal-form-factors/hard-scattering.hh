/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2016 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_SRC_NONLOCAL_FORM_FACTORS_HARD_SCATTERING_HH
#define EOS_GUARD_SRC_NONLOCAL_FORM_FACTORS_HARD_SCATTERING_HH 1

#include <complex>

namespace eos
{
    using std::complex;

    struct HardScattering
    {
        /*!
         * cf. [BFS:2001A], Eqs. (30)-(32), p. 8
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_q Mass of the internal loop quark
         * @param m_B Mass of the B meson
         */
        static complex<double> I1(const double & s, const double & u, const double & m_q, const double & m_B);

        /*!
         * Two-particle light-cone distribution amplitude up to 2nd order
         *         in Gegenbauer expansion of light mesons
         * cf. [BFS:2001A], Eq. (48)
         * @param u Relative contribution of the quark (versus ubar = 1-u for the
                 antiquark) to the light meson's energy.
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static double LCDA_2pt(const double & u, const double & a_1, const double & a_2);

        /*!
         * t_perp for q^2 = 0; cf. [BFS:2001A], Eq. (33)
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_q Mass of the internal loop quark
         * @param m_B Mass of the B meson
         */
        static complex<double> t_perp_s0(const double & u, const double & m_q, const double & m_B);

        /*!
         * t_perp; cf. [BFS:2001A], Eq. (27)
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_q Mass of the internal loop quark
         * @param m_B Mass of the B meson
         * @param m_M Mass of the K(K*) meson
         */
        static complex<double> t_perp(const double & s, const double & u, const double & m_q,
                                      const double & m_B, const double & m_M);

        /*!
         * t_par; cf. [BFS:2001A], Eq. (28)
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_q Mass of the internal loop quark
         * @param m_B Mass of the B meson
         * @param m_M Mass of the K(K*) meson
         */
        static complex<double> t_par(const double & s, const double & u, const double & m_q,
                                     const double & m_B, const double & m_M);

        /*!
         * [BFS:2001A], Eq. (26), term \propto C_8_eff
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_B Mass of the B meson
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static double j0(const double & s, const double & u, const double & m_B, const double & a_1,
                const double & a_2);

        /*!
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_q Mass of the internal loop quark
         * @param m_B Mass of the B meson
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static complex<double> j1(const double & s, const double & u, const double & m_q, const double & m_B,
                const double & a_1, const double & a_2);

        /*!
         * [BFS:2001A], Eq. (25)
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_q Mass of the internal loop quark
         * @param m_B Mass of the B meson
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static complex<double> j2(const double & s, const double & u, const double & m_q, const double & m_B,
                const double & a_1, const double & a_2);

        /*!
         * [BFS:2001A], Eq. (25)
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_B Mass of the B meson
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static complex<double> j2_massless(const double & s, const double & u, const double & m_B,
                const double & a_1, const double & a_2);

        /*!
         * [BFS:2001A], Eq. (25)
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_q Mass of the internal loop quark
         * @param m_B Mass of the B meson
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static complex<double> j3(const double & s, const double & u, const double & m_q, const double & m_B,
                const double & a_1, const double & a_2);

        /*!
         * [BFS:2001A], Eq. (25)
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_B Mass of the B meson
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static complex<double> j3_massless(const double & s, const double & u, const double & m_B,
                const double & a_1, const double & a_2);

        /*!
         * [BFS:2001A], Eq. (26), term \propto C1 - C6
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_q Mass of the internal loop quark
         * @param m_B Mass of the B meson
         * @param mu Renormalization scale
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static complex<double> j4(const double & s, const double & u, const double & m_q, const double & m_B,
                const double & mu, const double & a_1, const double & a_2);

        /*!
         * [BFS:2001A], Eq. (23), term \propto C_8_eff
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_q Mass of the internal loop quark
         * @param m_B Mass of the B meson
         * @param mu Renormalization scale
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static complex<double> j5(const double & s, const double & u, const double & m_q, const double & m_B,
                const double & mu, const double & a_1, const double & a_2);

        /*!
         * [BFS:2001A], Eq. (23), term \propto C_8_eff
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_q Mass of the internal loop quark
         * @param m_B Mass of the B meson
         * @param mu Renormalization scale
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static complex<double> j6(const double & s, const double & u, const double & m_q, const double & m_B,
                const double & mu, const double & a_1, const double & a_2);

        /*!
         * [BFS:2004A], Eq. (52), term \propto C_8_eff
         * @param s Dilepton invariant mass
         * @param u Relative contribution of the quark (versus ubar = 1-u for the antiquark) to the light meson's energy.
         * @param m_B Mass of the B meson
         * @param a_1 First Gegenbauer coefficient
         * @param a_2 Second Gegenbauer coefficient
         */
        static double j7(const double & s, const double & u, const double & m_B, const double & a_1, const double & a_2);
    };
}

#endif
