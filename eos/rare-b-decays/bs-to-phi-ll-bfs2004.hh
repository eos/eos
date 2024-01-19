/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_BS_TO_PHI_LL_BFS2004_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_BS_TO_PHI_LL_BFS2004_HH 1

#include <eos/rare-b-decays/bs-to-phi-ll-base.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>

namespace eos
{
    template <>
    class BsToPhiDileptonAmplitudes<tag::BFS2004> :
        public BsToPhiDilepton::AmplitudeGenerator
    {
        public:
            UsedParameter m_b_MSbar;
            UsedParameter m_c;
            UsedParameter m_s_MSbar;

            UsedParameter f_Bs;
            UsedParameter f_phi_par;
            UsedParameter f_phi_perp;
            UsedParameter lambda_B_p_inv;
            UsedParameter a_1_par;
            UsedParameter a_2_par;
            UsedParameter a_1_perp;
            UsedParameter a_2_perp;

            UsedParameter uncertainty_para;
            UsedParameter uncertainty_perp;
            UsedParameter uncertainty_long;

            UsedParameter uncertainty_xi_perp;
            UsedParameter uncertainty_xi_par;

            BooleanOption opt_ccbar_resonance;
            BooleanOption opt_use_nlo;

            bool ccbar_resonance;
            bool use_nlo;

            static const std::vector<OptionSpecification> options;

            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &)> qcdf_dilepton_massless_case;
            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &, const double &)> qcdf_dilepton_charm_case;
            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &, const double &)> qcdf_dilepton_bottom_case;

            std::string ff_relation;

            BsToPhiDileptonAmplitudes(const Parameters & p, const Options & o);
            ~BsToPhiDileptonAmplitudes();

            virtual BsToPhiDilepton::Amplitudes amplitudes(const double & q2) const;

            double m_b_PS() const;
            double mu_f() const;
            BsToPhiDilepton::DipoleFormFactors dipole_form_factors(const double & q2, const WilsonCoefficients<BToS> & wc) const;
            double norm(const double & q2) const;
            double xi_perp(const double & q2) const;
            double xi_par(const double & q2) const;

            virtual double real_C9_perp(const double & s) const;
            virtual double real_C9_para(const double & s) const;
            virtual double imag_C9_perp(const double & s) const;
            virtual double imag_C9_para(const double & s) const;
            virtual double H_perp_corrections(const double &) const
            {
                return 0.0;
            };
            virtual double H_para_corrections(const double &) const
            {
                return 0.0;
            };
            virtual double H_long_corrections(const double &) const
            {
                return 0.0;
            };
    };
}

#endif
