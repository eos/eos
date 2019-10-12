/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Danny van Dyk
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

#ifndef MASTER_GUARD_EOS_B_DECAYS_B_TO_D_PI_L_NU_HH
#define MASTER_GUARD_EOS_B_DECAYS_B_TO_D_PI_L_NU_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    class BToDPiLeptonNeutrino :
        public ParameterUser,
        public PrivateImplementationPattern<BToDPiLeptonNeutrino>
    {
        public:
            BToDPiLeptonNeutrino(const Parameters & parameters, const Options & options);
            ~BToDPiLeptonNeutrino();

            /*!
             * 1-dim. PDFs as functions of cos(theta_D), cost(theta_L), chi, q2, and w
             */
            double differential_pdf_chi(const double & chi) const;
            double differential_pdf_d(const double & c_d) const;
            double differential_pdf_l(const double & c_l) const;
            double differential_pdf_q2(const double & q2) const;
            double differential_pdf_w(const double & w) const;

            /*!
             * Partially-integrated 1-dim. PDFs for cos(theta_D), cost(theta_L), and chi
             */
            double integrated_pdf_d(const double & c_d_min, const double & c_d_max) const;
            double integrated_pdf_l(const double & c_l_min, const double & c_l_max) const;
            double integrated_pdf_chi(const double & chi_min, const double & chi_max) const;
            double integrated_pdf_w(const double & w_min, const double & w_max) const;

            double integrated_lepton_polarization(const double & q2_min, const double & q2_max) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_c_d;
            static const std::string kinematics_description_c_l;
            static const std::string kinematics_description_chi;
            static const std::string kinematics_description_q2;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;
    };
}

#endif
