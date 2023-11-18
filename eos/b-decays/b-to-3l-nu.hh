/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Stephan Kürten
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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_3L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_3L_NU_HH 1

#include <cmath>
#include <eos/form-factors/mesonic.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>
#include <string>

namespace eos
{
    class BToThreeLeptonsNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToThreeLeptonsNeutrino>
    {
        private:
            double _asymmetry_numerator(const double & q2,const double & k2) const;

        public:
            BToThreeLeptonsNeutrino(const Parameters & parameters, const Options & options);
            ~BToThreeLeptonsNeutrino();

            /*!
             * observables of up to 5 kinematic variables
             * q2 is the invariant mass of the off-shell photon in the range 4 m_l'^2 <= q2 <= (m_B - m_l )^2
             * k2 is the invariant mass of the W-meson in the range m_l^2 <= k2 <= (m_B - sqrt(q2) )^2
             * z_gamma is the angle between the negatively charged lepton l' and the negative z-axis
             * z_w is the angle between the charged lepton l and the positive z-axis
             * phi is the angle between the q2 plane and the k2 plane
             */
            double double_differential_branching_ratio(const double & q2, const double & k2) const;
            double quintuple_differential_branching_ratio(const double & q2, const double & k2, const double & z_gamma, const double & z_w, const double & phi) const;
            double double_differential_decay_width(const double & q2, const double & k2) const;
            double quintuple_differential_decay_width(const double & q2, const double & k2, const double & z_gamma, const double & z_w, const double & phi) const;
            double integrated_branching_ratio(const double & q2_min, const double & q2_max, const double & k2_min, const double & k2_max) const;
            double double_differential_forward_backward_asymmetry(const double & q2, const double & k2) const;
            double integrated_forward_backward_asymmetry(const double & q2_min,const double & q2_max, const double & k2_min, const double & k2_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_q2;
            static const std::string kinematics_description_k2;
            static const std::string kinematics_description_z_gamma;
            static const std::string kinematics_description_z_w;
            static const std::string kinematics_description_phi;

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
}

#endif
