/*
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_BS_TO_PHI_CHARMONIUM_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_BS_TO_PHI_CHARMONIUM_HH 1

#include <eos/maths/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /**
     * Decay: B_s -> phi psi
     * with psi a narrow charmonium: psi = { J/psi, psi(2S) }
     */
    class BsToPhiCharmonium :
        public ParameterUser,
        public PrivateImplementationPattern<BsToPhiCharmonium>
    {
        public:
            ///@name Basic operations.
            ///@{

            BsToPhiCharmonium(const Parameters & parameters, const Options & options);
            ~BsToPhiCharmonium();

            ///@}

            ///@name Observables
            ///@{

            /// Branching ratio
            double branching_ratio() const;

            /// Polarzation fractions and relative phases of the amplitudes
            double perp_polarization() const;
            double para_polarization() const;
            double long_polarization() const;

            double long_phase() const;
            double delta_perp_long() const;
            double delta_para_long() const;

            /// Angular observables as detected in the decay B_s -> phi psi (-> l^+ l^-)
            double S_1s_LHCb() const;
            double S_1c_LHCb() const;
            double S_3_LHCb() const;
            double S_4_LHCb() const;
            double S_8_LHCb() const;
            double S_9_LHCb() const;

            ///@}

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
