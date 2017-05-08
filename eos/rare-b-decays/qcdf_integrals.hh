/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_QCDF_INTEGRALS_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_QCDF_INTEGRALS_HH 1

#include <eos/utils/complex.hh>

namespace eos
{
    class QCDFIntegrals
    {
        public:
            struct Results;

            /*!
             * Return all QCDF Integrals for a b quark-antiquark loop with s = 0, according to [vD2011].
             *
             * @param m_b          Pole mass of the b quark.
             * @param m_B          Mass of the parent B quark.
             * @param mu           Renormalization scale mu.
             * @param a_1_perp     First Gegenbauer moment for the perpendicular amplitude.
             * @param a_2_perp     Second Gegenbauer moment for the perpendicular amplitude.
             * @param a_1_parallel First Gegenbauer moment for the paralle amplitude.
             * @param a_2_parallel Second Gegenbauer moment for the parallel amplitude.
             */
            static Results photon_bottom_case(const double & m_b, const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel);

            /*!
             * Return all QCDF Integrals for a c quark-antiquark loop with s = 0, according to [vD2011].
             *
             * @param m_c          Pole mass of the c quark.
             * @param m_B          Mass of the parent B quark.
             * @param mu           Renormalization scale mu.
             * @param a_1_perp     First Gegenbauer moment for the perpendicular amplitude.
             * @param a_2_perp     Second Gegenbauer moment for the perpendicular amplitude.
             * @param a_1_parallel First Gegenbauer moment for the paralle amplitude.
             * @param a_2_parallel Second Gegenbauer moment for the parallel amplitude.
             */
            static Results photon_charm_case(const double & m_c, const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel);

            /*!
             * Return all QCDF Integrals for a u,d,s (i.e. massless) quark-antiquark loops with s = 0, according to [vD2011].
             *
             * @param m_B          Mass of the parent B quark.
             * @param mu           Renormalization scale mu.
             * @param a_1_perp     First Gegenbauer moment for the perpendicular amplitude.
             * @param a_2_perp     Second Gegenbauer moment for the perpendicular amplitude.
             * @param a_1_parallel First Gegenbauer moment for the paralle amplitude.
             * @param a_2_parallel Second Gegenbauer moment for the parallel amplitude.
             */
            static Results photon_massless_case(const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel);

            /*!
             * Return all QCDF Integrals for a b quark-antiquark loop, according to [vD2011].
             *
             * @param s            Invariant quark-antiquark mass square.
             * @param m_b          Pole mass of the b quark.
             * @param m_B          Mass of the parent B quark.
             * @param mu           Renormalization scale mu.
             * @param a_1_perp     First Gegenbauer moment for the perpendicular amplitude.
             * @param a_2_perp     Second Gegenbauer moment for the perpendicular amplitude.
             * @param a_1_parallel First Gegenbauer moment for the paralle amplitude.
             * @param a_2_parallel Second Gegenbauer moment for the parallel amplitude.
             */
            static Results dilepton_bottom_case(const double & s, const double & m_b, const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel);

            /*!
             * Return all QCDF Integrals for a c quark-antiquark loop, according to [vD2011].
             *
             * @param s            Invariant quark-antiquark mass square.
             * @param m_c          Pole mass of the c quark.
             * @param m_B          Mass of the parent B quark.
             * @param mu           Renormalization scale mu.
             * @param a_1_perp     First Gegenbauer moment for the perpendicular amplitude.
             * @param a_2_perp     Second Gegenbauer moment for the perpendicular amplitude.
             * @param a_1_parallel First Gegenbauer moment for the paralle amplitude.
             * @param a_2_parallel Second Gegenbauer moment for the parallel amplitude.
             */
            static Results dilepton_charm_case(const double & s, const double & m_c, const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel);

            /*!
             * Return all QCDF Integrals for a u,d,s (i.e. massless) quark-antiquark loop for massless quarks, according to [vD2011].
             *
             * @param s            Invariant quark-antiquark mass square.
             * @param m_B          Mass of the parent B quark.
             * @param mu           Renormalization scale.
             * @param a_1_perp     First Gegenbauer moment for the perpendicular amplitude.
             * @param a_2_perp     Second Gegenbauer moment for the perpendicular amplitude.
             * @param a_1_parallel First Gegenbauer moment for the paralle amplitude.
             * @param a_2_parallel Second Gegenbauer moment for the parallel amplitude.
             */
            static Results dilepton_massless_case(const double & s, const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel);
    };

    /*!
     * Ensemble of the individual integral results.
     */
    struct QCDFIntegrals::Results
    {
        // cf. [vD2011], Eqs. (x), (y) for the perpendicular amplitudes
        complex<double> j0_perp;
        complex<double> j0bar_perp;
        complex<double> j1_perp;
        complex<double> j2_perp;
        complex<double> j4_perp;
        complex<double> j5_perp;
        complex<double> j6_perp; // This integral arises in perpendicular amplitudes, but depends on parallel Gegenbauer moments!
        double j7_perp;

        // cf. [vD2011], Eqs. (x), (y) for the parallel amplitudes
        complex<double> j0_parallel;
        complex<double> j1_parallel;
        complex<double> j3_parallel;
        complex<double> j4_parallel;

        /* For convenience, the integrals jtilde_{1,2} from j_{1,2,3} */
        complex<double> jtilde1_perp;
        complex<double> jtilde2_parallel;
    };
}

#endif
