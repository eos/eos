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

#ifndef EOS_GUARD_SRC_UTILS_LCDA_HH
#define EOS_GUARD_SRC_UTILS_LCDA_HH 1

#include <eos/utils/qcd.hh>

namespace eos
{
    /*!
     * Groups all functions related to LightCone Distribution Amplitudes
     * (LCDAs).
     */
    struct LCDA
    {
        /*!
         * Evolve the nth Gegenbauer moment a_n_0 from a scale mu_0 to a scale mu according to the LL RGE.
         *
         * Calculation according to [BBL:2006A].
         *
         * @param a_n_0 The Gegenbauer moment at the scale mu_0.
         * @param n     The index of the Gegenbauer moment.
         * @param eta   The ratio alpha_s(mu)/alpha_s(mu_0).
         * @param beta  The coefficients of the QCD beta function.
         */
        static double evolve_gegenbauer_moment(const double & a_n_0, const unsigned & n, const double & eta, const QCD::BetaFunction & beta);
    };
}

#endif
