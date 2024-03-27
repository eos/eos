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

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_GVDV2020_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_GVDV2020_HH 1

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>
#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    template <>
    class BToKstarDileptonAmplitudes<tag::GvDV2020> :
        public BToKstarDilepton::AmplitudeGenerator
    {
        public:
            UsedParameter m_b_MSbar;
            UsedParameter m_s_MSbar;

            UsedParameter f_B;
            UsedParameter f_Kstar_par;
            UsedParameter lambda_B_p_inv;

            QuarkFlavorOption q;

            SwitchOption opt_nonlocal_formfactor;
            NonlocalFormFactorPtr<PToV> nonlocal_formfactor;

            static const std::vector<OptionSpecification> options;

            BToKstarDileptonAmplitudes(const Parameters & p, const Options & o);
            ~BToKstarDileptonAmplitudes() = default;

            BToKstarDilepton::FormFactorCorrections sb_contributions(const double & q2, const WilsonCoefficients<BToS> & wc) const;
            double m_b_PS() const;
            double mu_f() const;

            virtual double real_C9_perp(const double & s) const;
            virtual double real_C9_para(const double & s) const;
            virtual double imag_C9_perp(const double & s) const;
            virtual double imag_C9_para(const double & s) const;
            virtual double H_perp_corrections(const double & s) const;
            virtual double H_para_corrections(const double & s) const;
            virtual double H_long_corrections(const double & s) const;

            virtual BToKstarDilepton::Amplitudes amplitudes(const double & q2) const;
    };
}

#endif
