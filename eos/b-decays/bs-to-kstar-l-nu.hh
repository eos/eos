/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013, 2015, 2016 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_B_DECAYS_BS_TO_KSTAR_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_BS_TO_KSTAR_L_NU_HH 1

#include <eos/decays.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    class BsToKstarLeptonNeutrinoRatios;

    /*
     * Decay: B_s -> K^* l^- nubar, cf. [FMvD2014]
     */
    class BsToKstarLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BsToKstarLeptonNeutrino>
    {
        public:
            friend BsToKstarLeptonNeutrinoRatios;

            BsToKstarLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BsToKstarLeptonNeutrino();

            // Four Differential Observables, cf. [FMvD2015] Eq. (7)
            double four_differential_decay_width(const double & s, const double & c_theta_l, const double & c_theta_k, const double & phi) const;

            // Helicity Form Factors
            double Ftime(const double & s) const;
            double Flong(const double & s) const;
            double Fpara(const double & s) const;
            double Fperp(const double & s) const;

            // Single Differential Observables, cf. [FMvD2015]
            double differential_branching_ratio(const double & s) const;
            double differential_decay_width(const double & s) const;
            double differential_forward_backward_asymmetry(const double & s) const;
            double differential_longitudinal_polarisation(const double & s) const;
            double differential_transversal_polarisation(const double & s) const;
            double differential_transverse_asymmetry_2(const double & s) const;
            double differential_transverse_asymmetry_3(const double & s) const;
            double differential_transverse_asymmetry_4(const double & s) const;
            double differential_transverse_asymmetry_5(const double & s) const;
            double differential_transverse_asymmetry_re(const double & s) const;
            double differential_transverse_asymmetry_im(const double & s) const;
            double differential_h_1(const double & s) const;
            double differential_h_2(const double & s) const;
            double differential_h_3(const double & s) const;
            double differential_h_4(const double & s) const;
            double differential_h_5(const double & s) const;

            // Integrated Observables, cf. [FMvD2015]
            double integrated_decay_width(const double & s_min, const double & s_max) const;
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const;
            double integrated_transversal_polarisation(const double & s_min, const double & s_max) const;
            double integrated_h_1(const double & s_min, const double & s_max) const;
            double integrated_h_2(const double & s_min, const double & s_max) const;
            double integrated_h_3(const double & s_min, const double & s_max) const;
            double integrated_h_4(const double & s_min, const double & s_max) const;
            double integrated_h_5(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_5(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_re(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_im(const double & s_min, const double & s_max) const;
            double integrated_s_1s(const double & s_min, const double & s_max) const;
            double integrated_s_1c(const double & s_min, const double & s_max) const;
            double integrated_s_2s(const double & s_min, const double & s_max) const;
            double integrated_s_2c(const double & s_min, const double & s_max) const;
            double integrated_s_3(const double & s_min, const double & s_max) const;
            double integrated_s_4(const double & s_min, const double & s_max) const;
            double integrated_s_5(const double & s_min, const double & s_max) const;
            double integrated_s_6s(const double & s_min, const double & s_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_s;
            static const std::string kinematics_description_c_theta_l;
            static const std::string kinematics_description_c_theta_k;
            static const std::string kinematics_description_phi;
    };

    class BsToKstarLeptonNeutrinoRatios :
        public ParameterUser,
        public PrivateImplementationPattern<BsToKstarLeptonNeutrinoRatios>
    {
        public:
            BsToKstarLeptonNeutrinoRatios(const Parameters & parameters, const Options & options);
            ~BsToKstarLeptonNeutrinoRatios();

            double ratio_long() const;
            double ratio_para() const;
            double ratio_perp() const;
    };
}

#endif
