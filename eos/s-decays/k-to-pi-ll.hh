/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_S_DECAYS_K_TO_PI_LL_HH
#define EOS_GUARD_SRC_S_DECAYS_K_TO_PI_LL_HH 1

#include <eos/maths/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{

    /*
     * Decay: K -> pi l lbar.
     */
    class KToPiDilepton :
        public ParameterUser,
        public PrivateImplementationPattern<KToPiDilepton>
    {
        public:
            KToPiDilepton(const Parameters & parameters, const Options & options);
            ~KToPiDilepton();

            struct AngularCoefficients;
            struct Amplitudes;
            class AmplitudeGenerator;
            struct DipoleFormFactors;

            // Differential Observables
            double differential_branching_ratio(const double & q2) const;
            double differential_flat_term(const double & q2) const;
            double differential_forward_backward_asymmetry(const double & q2) const;
            double differential_ratio_muons_electrons(const double & q2) const;

            // Two Differential decay with in LHCb angular convention
            double two_differential_decay_width(const double & q2, const double & c_theta_l) const;

            // Integrated Observables
            double integrated_decay_width(const double & q2_min, const double & q2_max) const;
            double integrated_branching_ratio(const double & q2_min, const double & q2_max) const;
            double integrated_flat_term(const double & q2_min, const double & q2_max) const;
            double integrated_forward_backward_asymmetry(const double & q2_min, const double & q2_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_q2;
            static const std::string kinematics_description_c_theta_l;

            /*!
             * Auxiliary methods for unit tests and diagnostic purposes.
             */
            Amplitudes amplitudes(const double & q2) const;
            std::array<double, 3> angular_coefficients(const double & q2) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };

    /*!
     * Amplitudes for the decay K -> pi l lbar.
     */
    struct KToPiDilepton::Amplitudes
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
