/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Philip Lüghausen
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_B_LCDAS_PARAM_HH
#define EOS_GUARD_EOS_FORM_FACTORS_B_LCDAS_PARAM_HH 1

#include <eos/form-factors/b-lcdas.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    namespace b_lcdas
    {
        /*!
        * Parametrization of the B-meson LCDAs
        */
        class Param:
            public BMesonLCDAs,
            public PrivateImplementationPattern<Param>
        {
            public:
                Param(const Parameters &, const Options &);
                ~Param();

                /*!
                *  B to gamma l nu
                */
                double L0() const;
                double L0inc(const double & Omega) const;
                double Binc(const double & Omega, const double & sigma) const;

                /*!
                * Leading twist two-particle LCDAs
                * omega: plus-component of the spectator momentum
                */
                double phi_plus(const double & omega) const;
                double phi_minus(const double & omega) const;
                double phi_bar(const double & omega) const;
                double phi_bar_d1(const double & omega) const;

                /*!
                * Next-to-leading twist two-particle LCDAs
                * omega: plus-component of the spectator momentum
                */
                double g_plus(const double & omega) const;
                double g_plus_d1(const double & omega) const;
                double g_plus_d2(const double & omega) const;

                double g_minusWW(const double & omega) const;
                double g_minusWW_d1(const double & omega) const;
                double g_minusWW_d2(const double & omega) const;

                double g_bar(const double & omega) const;
                double g_bar_d1(const double & omega) const;
                double g_bar_d2(const double & omega) const;
                double g_bar_d3(const double & omega) const;

                /*!
                * Leading power three-particle LCDAs
                * omega_1: plus-component of the spectator momentum
                * omega_2: plus-component of the gluon momentum
                * */
                double phi_3(const double & omega_1, const double & omega_2) const;
                double phi_4(const double & omega_1, const double & omega_2) const;

                double phi_bar_3(const double & omega_1, const double & omega_2) const;
                double phi_bar_4(const double & omega_1, const double & omega_2) const;

                double phi_bar2_3(const double & omega_1, const double & omega_2) const;
                double phi_bar2_4(const double & omega_1, const double & omega_2) const;

                double phi_bar_bar_3(const double & omega_1, const double & omega_2) const;
                double phi_bar_bar_4(const double & omega_1, const double & omega_2) const;

                double psi_bar_4(const double & omega_1, const double & omega_2) const;
                double chi_bar_4(const double & omega_1, const double & omega_2) const;

                double psi_bar_bar_4(const double & omega_1, const double & omega_2) const;
                double chi_bar_bar_4(const double & omega_1, const double & omega_2) const;
                /*!
                * Pseudo observables for the two-particle LCDAs
                */
                double inverse_lambda_plus() const;

                /*!
                * Leading power three-particle LCDAs
                * omega: plus-component of the spectator momentum
                * xi:    plus-component of the gluon momentum
                * */
                double psi_A(const double & omega, const double & xi) const;
                double psi_V(const double & omega, const double & xi) const;
                double X_A(const double & omega, const double & xi) const;
                double Y_A(const double & omega, const double & xi) const;

                /*!
                * Auxiliary functions for the three-particle LCDAs
                * See [KMO2006], below eq. (72), p. 28 for their definition.
                */
                double Xbar_A(const double & omega, const double & xi) const;
                double Ybar_A(const double & omega, const double & xi) const;

                /* Internal diagnostics */
                Diagnostics diagnostics() const;
        };
    }
}

#endif
