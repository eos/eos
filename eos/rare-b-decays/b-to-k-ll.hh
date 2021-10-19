/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2015, 2016 Danny van Dyk
 * Copyright (c) 2021 Méril Reboud
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_B_TO_K_LL_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_B_TO_K_LL_HH 1

#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{

    /*
     * Decay: B -> K l lbar.
     */
    class BToKDilepton :
        public ParameterUser,
        public PrivateImplementationPattern<BToKDilepton>
    {
        public:
            BToKDilepton(const Parameters & parameters, const Options & options);
            ~BToKDilepton();

            struct AngularCoefficients;
            struct Amplitudes;
            class AmplitudeGenerator;
            struct DipoleFormFactors;

            // Differential Observables
            double differential_branching_ratio(const double & s) const;
            double differential_flat_term(const double & s) const;
            double differential_forward_backward_asymmetry(const double & s) const;
            double differential_ratio_muons_electrons(const double & s) const;

            // Two Differential decay with in LHCb angular convention
            double two_differential_decay_width(const double & s, const double & c_theta_l) const;

            // Integrated Observables
            double integrated_decay_width(const double & s_min, const double & s_max) const;
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_branching_ratio_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_cp_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_flat_term(const double & s_min, const double & s_max) const;
            double integrated_flat_term_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_ratio_muons_electrons(const double & s_min, const double & s_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_s;
            static const std::string kinematics_description_c_theta_l;

            /*!
             * Auxilliary methods for unit tests and diagnostic purposes.
             */
            Amplitudes amplitudes(const double & q2) const;
            std::array<double, 3> angular_coefficients(const double & q2) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;
    };

    /*!
     * Amplitudes for the decay B -> K l lbar.
     */
    struct BToKDilepton::Amplitudes
    {
        complex<double> F_A;
        complex<double> F_V;
        complex<double> F_S;
        complex<double> F_P;
        complex<double> F_T;
        complex<double> F_T5;
    };
}

#endif
