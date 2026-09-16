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

#include <eos/maths/power-of.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/rare-c-decays/d-to-psd-l-l-bght2019.hh>
#include <eos/utils/memoise.hh>

#include <gsl/gsl_sf.h>

namespace eos
{
    using namespace std::placeholders;

    DToPseudoscalarLeptonLeptonAmplitudes<tag::BGHT2019>::DToPseudoscalarLeptonLeptonAmplitudes(const Parameters & p, const Options & o) :
        AmplitudeGenerator(p, o),
        m_rho(p["mass::rho^0"], u),
        m_omega(p["mass::omega"], u),
        m_phi(p["mass::phi"], u),
        tau_rho(p["life_time::rho^0"], u),
        tau_omega(p["life_time::omega"], u),
        tau_phi(p["life_time::phi"], u),
        a_rho(p["::res_a_rho@GHM2021"], u),
        a_omega(p["::res_a_omega@GHM2021"], u),
        a_phi(p["::res_a_phi@GHM2021"], u),
        delta_rho(p["::res_delta_rho@GHM2021"], u),
        delta_omega_m_rho(p["::res_delta_omega_m_rho@GHM2021"], u),
        delta_phi_m_rho(p["::res_delta_phi_m_rho@GHM2021"], u),
        opt_nonlocal_formfactor(o, options, "nonlocal-formfactor"_ok),
        nonlocal_formfactor(NonlocalFormFactor<PToP>::make(_process() + "::" + opt_nonlocal_formfactor.value(), p, o))
    {
        Context ctx("When constructing B->Kll GvDV2020 amplitudes");
    }

    DToPseudoscalarLeptonLeptonAmplitudes<tag::BGHT2019>::~DToPseudoscalarLeptonLeptonAmplitudes() {}

    const std::vector<OptionSpecification> DToPseudoscalarLeptonLeptonAmplitudes<tag::BGHT2019>::options{
        {                   "q"_ok,                                            { "d"_ov, "u"_ov },        "d"_ov },
        { "nonlocal-formfactor"_ok, { "GvDV2020"_ov, "GRvDV2022order5"_ov, "GRvDV2022order6"_ov }, "GvDV2020"_ov }
    };

    /* Amplitudes */
    DToPseudoscalarLeptonLepton::Amplitudes
    DToPseudoscalarLeptonLeptonAmplitudes<tag::BGHT2019>::amplitudes(const double & q2) const
    {
        DToPseudoscalarLeptonLepton::Amplitudes result;

        const auto wc = model->wilson_coefficients_uc(opt_l.value(), opt_cp_conjugate.value());
        const double mcatmu = model->m_c_msbar(mu);
        //const double msatmu = model->m_s_msbar(mu);

        const double m_D2 = m_D * m_D, m_P2 = m_P * m_P;

        // cf. [GvDV:2020A] Eq. (A.5)
        const double calF_plus = form_factors->f_p(q2), calF_time = form_factors->f_0(q2), calF_T_plus = q2 / m_D / (m_D + m_P) * form_factors->f_t(q2);

        const complex<double> calH_plus = nonlocal_formfactor->H_plus(q2);

        double F_Tkin = calF_T_plus / calF_plus * 2.0 * std::sqrt(lambda(q2)) * beta_l(q2) * m_D / q2;
        double F_Skin = calF_time / calF_plus * 0.5 * (m_D2 - m_P2) / (mcatmu - m_s_MSbar);

        // Wilson coefficients
        const complex<double> c9_p = wc.c9() + wc.c9prime(), c10_p = wc.c10() + wc.c10prime(), c7_p = wc.c7() + wc.c7prime();

        // cf. [BHP:2007A], Eq. (3.2), p. 3 and 4 or [BKMS:2012A] (1205.5811)
        result.F_A  = c10_p;
        result.F_T  = F_Tkin * wc.cT();
        result.F_T5 = F_Tkin * wc.cT5();
        result.F_S  = F_Skin * (wc.cS() + wc.cSprime());
        result.F_P  = F_Skin * (wc.cP() + wc.cPprime()) + m_l() * c10_p * ((m_D2 - m_P2) / q2 * (calF_time / calF_plus - 1.0) - 1.0);
        result.F_V  = c9_p + 2.0 * m_b_MSbar() * m_D / q2 * c7_p * calF_T_plus / calF_plus
                     //+ 2.0 * m_b_PS() / m_D / xi_pseudo(q2) * (dff.calT - 16.0 * power_of<2>(M_PI) * power_of<3>(m_D()) / m_b_PS() / q2 * calH_plus)
                     + 8.0 * m_l * m_D / q2 * calF_T_plus / calF_plus * wc.cT();

        return result;
    }
} // namespace eos
