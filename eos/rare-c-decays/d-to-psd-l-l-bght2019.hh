/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026    Carolina Bolognani
 * Copyright (c) 2026    Dominik Suelmann
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

#ifndef MASTER_GUARD_EOS_RARE_C_DECAYS_D_TO_PSDEUDOSCALAR_LEPTON_LEPTON_BGHT2019_HH
#define MASTER_GUARD_EOS_RARE_C_DECAYS_D_TO_PSDEUDOSCALAR_LEPTON_LEPTON_BGHT2019_HH 1

#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/rare-c-decays/d-to-psd-l-l-base.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    template <> class DToPseudoscalarLeptonLeptonAmplitudes<tag::BGHT2019> : public DToPseudoscalarLeptonLepton::AmplitudeGenerator
    {
        public:
            UsedParameter m_rho;
            UsedParameter m_omega;
            UsedParameter m_phi;
            UsedParameter m_eta;
            UsedParameter m_eta_p;
            UsedParameter tau_rho;
            UsedParameter tau_omega;
            UsedParameter tau_phi;
            UsedParameter tau_eta;
            UsedParameter tau_eta_p;

            // resonance parameters
            UsedParameter a_rho;
            UsedParameter a_omega;
            UsedParameter a_phi;
            UsedParameter delta_rho;
            UsedParameter delta_omega_m_rho;
            UsedParameter delta_phi_m_rho;

            UsedParameter a_eta;
            UsedParameter a_eta_p;
            UsedParameter delta_eta;
            UsedParameter delta_eta_p_m_eta;

            DToPseudoscalarLeptonLeptonAmplitudes(const Parameters & p, const Options & o);
            ~DToPseudoscalarLeptonLeptonAmplitudes();

            virtual DToPseudoscalarLeptonLepton::Amplitudes amplitudes(const double & q2) const;
    };
} // namespace eos

#endif
