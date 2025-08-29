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

#include <eos/form-factors/lcda.hh>

#include <cmath>

#include <gsl/gsl_sf_psi.h>

namespace eos
{
    double
    LCDA::evolve_gegenbauer_moment(const double & a_n_0, const unsigned & n, const double & eta, const QCD::BetaFunction & beta)
    {
        // Euler-Mascheroni constant.
        static const double gamma_e = 0.57721566490153286061;

        // cf. [BBL:2006], Eq. (2.14), p. 5
        double gamma_0 = 8 * QCD::casimir_f * (gsl_sf_psi_int(n + 2) + gamma_e - 0.75 - 0.5 / (n + 1) / (n + 2));

        // cf. [BBL:2006], Eq. (2.13), p. 5
        return std::pow(eta, gamma_0 / 2.0 / beta[0]) * a_n_0;
    }
}
