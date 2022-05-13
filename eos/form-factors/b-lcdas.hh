/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
 * Copyright (c) 2018 Nico Gubernari
 * Copyright (c) 2018 Ahmet Kokulu
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_B_LCDAS_HH
#define EOS_GUARD_EOS_FORM_FACTORS_B_LCDAS_HH 1

#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    /*!
     * Decomposition of B-meson to vacuum matrix elements of light-cone
     * dominated operators.
     *
     * The abtract base class defines the interface used in sum rules.
     * TODO: Add interface for B to gamma l nu
     * TODO: Default implementation of the methods that throws NotImplementedError?
     *
     */
    class BMesonLCDAs:
        public ParameterUser
    {
        public:
            // Can this work...? The logic should be implemented in the observables.
            // static std::shared_ptr<BMesonLCDAs> make(const std::string & model, const Parameters & parameters, const Options & options) {
            //     if (model == "exp") {
            //         return std::make_shared<b_lcdas::Exponential>(b_lcdas::Exponential(parameters, options));
            //     }
            //     else {
            //         throw std::logic_error("No other model than 'exp' implemented (WIP)");
            //     }
            // };

            /*!
             * Leading twist two-particle LCDAs
             *
             * omega: plus-component of the spectator momentum
             */
            virtual double phi_plus(const double & omega) const = 0;
            virtual double phi_minus(const double & omega) const = 0;
            virtual double phi_bar(const double & omega) const = 0;
            virtual double phi_bar_d1(const double & omega) const = 0;

            /*!
             * Next-to-leading twist two-particle LCDAs
             *
             * omega: plus-component of the spectator momentum
             */
            virtual double g_plus(const double & omega) const = 0;
            virtual double g_plus_d1(const double & omega) const = 0;
            virtual double g_plus_d2(const double & omega) const = 0;

            virtual double g_minusWW(const double & omega) const = 0;
            virtual double g_minusWW_d1(const double & omega) const = 0;
            virtual double g_minusWW_d2(const double & omega) const = 0;

            virtual double g_bar(const double & omega) const = 0;
            virtual double g_bar_d1(const double & omega) const = 0;
            virtual double g_bar_d2(const double & omega) const = 0;
            virtual double g_bar_d3(const double & omega) const = 0;

            /*!
             * Leading power three-particle LCDAs
             *
             * omega_1: plus-component of the spectator momentum
             * omega_2: plus-component of the gluon momentum
             * */
            virtual double phi_3(const double & omega_1, const double & omega_2) const = 0;
            virtual double phi_4(const double & omega_1, const double & omega_2) const = 0;

            virtual double phi_bar_3(const double & omega_1, const double & omega_2) const = 0;
            virtual double phi_bar_4(const double & omega_1, const double & omega_2) const = 0;

            virtual double phi_bar2_3(const double & omega_1, const double & omega_2) const = 0;
            virtual double phi_bar2_4(const double & omega_1, const double & omega_2) const = 0;

            virtual double phi_bar_bar_3(const double & omega_1, const double & omega_2) const = 0;
            virtual double phi_bar_bar_4(const double & omega_1, const double & omega_2) const = 0;

            virtual double psi_bar_4(const double & omega_1, const double & omega_2) const = 0;
            virtual double chi_bar_4(const double & omega_1, const double & omega_2) const = 0;

            virtual double psi_bar_bar_4(const double & omega_1, const double & omega_2) const = 0;
            virtual double chi_bar_bar_4(const double & omega_1, const double & omega_2) const = 0;
            /*!
             * Pseudo observables for the two-particle LCDAs
             */
            virtual double inverse_lambda_plus() const = 0;

            /*!
             * Leading power three-particle LCDAs
             *
             * omega: plus-component of the spectator momentum
             * xi:    plus-component of the gluon momentum
             * */
            virtual double psi_A(const double & omega, const double & xi) const = 0;
            virtual double psi_V(const double & omega, const double & xi) const = 0;
            virtual double X_A(const double & omega, const double & xi) const = 0;
            virtual double Y_A(const double & omega, const double & xi) const = 0;

            /*!
             * Auxiliary functions for the three-particle LCDAs
             *
             * See [KMO2006], below eq. (72), p. 28 for their definition.
             */
            virtual double Xbar_A(const double & omega, const double & xi) const = 0;
            virtual double Ybar_A(const double & omega, const double & xi) const = 0;

            // TODO: Does this method belong here? Is this for testing implementations only?
            /* Internal diagnostics */
            // virtual Diagnostics diagnostics() const = 0;
    };
}

#endif
