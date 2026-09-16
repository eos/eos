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
#include <eos/rare-c-decays/d-to-psd-l-l-bght2019.hh>
#include <eos/utils/memoise.hh>

#include <gsl/gsl_sf.h>

#include <iomanip>
#include <iostream>

namespace eos
{
    using namespace std::placeholders;

    DToPseudoscalarLeptonLeptonAmplitudes<tag::BGHT2019>::DToPseudoscalarLeptonLeptonAmplitudes(const Parameters & p, const Options & o) :
        AmplitudeGenerator(p, o),
        m_rho(p["mass::rho^0"], *this),
        m_omega(p["mass::omega"], *this),
        m_phi(p["mass::phi"], *this),
        m_eta(p["mass::eta"], *this),
        m_eta_p(p["mass::eta_prime"], *this),
        tau_rho(p["life_time::rho^0"], *this),
        tau_omega(p["life_time::omega"], *this),
        tau_phi(p["life_time::phi"], *this),
        tau_eta(p["life_time::eta"], *this),
        tau_eta_p(p["life_time::eta_prime"], *this),
        a_rho(p[_process_q() + "::res_a_rho@BGHT2019"], *this),
        a_omega(p[_process_q() + "::res_a_omega@BGHT2019"], *this),
        a_phi(p[_process_q() + "::res_a_phi@BGHT2019"], *this),
        delta_rho(p[_process_q() + "::res_delta_rho@BGHT2019"], *this),
        delta_omega_m_rho(p[_process_q() + "::res_delta_omega_m_rho@BGHT2019"], *this),
        delta_phi_m_rho(p[_process_q() + "::res_delta_phi_m_rho@BGHT2019"], *this),
        a_eta(p[_process_q() + "::res_a_eta@BGHT2019"], *this),
        a_eta_p(p[_process_q() + "::res_a_eta_prime@BGHT2019"], *this),
        delta_eta(p[_process_q() + "::res_delta_eta@BGHT2019"], *this),
        delta_eta_p_m_eta(p[_process_q() + "::res_delta_eta_prime_m_eta@BGHT2019"], *this)
    {
        Context ctx("When constructing D(s)->Pll BGHT2019 amplitudes");
    }

    DToPseudoscalarLeptonLeptonAmplitudes<tag::BGHT2019>::~DToPseudoscalarLeptonLeptonAmplitudes() {}

    /* Amplitudes */
    DToPseudoscalarLeptonLepton::Amplitudes
    DToPseudoscalarLeptonLeptonAmplitudes<tag::BGHT2019>::amplitudes(const double & q2) const
    {
        DToPseudoscalarLeptonLepton::Amplitudes result;

        const auto   wc     = model->wilson_coefficients_uc(opt_l.value(), opt_cp_conjugate.value());
        const double mcatmu = model->m_c_msbar(mu);

        const double m_D2 = m_D * m_D, m_P2 = m_P * m_P;

        // cf. [GvDV:2020A] Eq. (A.5)
        const double calF_plus = form_factors->f_p(q2), calF_time = form_factors->f_0(q2), calF_T_plus = q2 / m_D / (m_D + m_P) * form_factors->f_t(q2);

        // const complex<double> calH_plus = nonlocal_formfactor->H_plus(q2);

        double F_Tkin = calF_T_plus / calF_plus * 2.0 * std::sqrt(lambda(q2)) * beta_l(q2) * m_D / q2;
        double F_Skin = calF_time / calF_plus * 0.5 * (m_D2 - m_P2) / mcatmu;

        // Wilson coefficients
        const complex<double> c9_p = (wc.c9() + wc.c9prime()), c10_p = (wc.c10() + wc.c10prime()), c7_p = (wc.c7() + wc.c7prime());

        // [BGHT:2019A] Eq. (6) without isospin assumption
        const complex<double> c9R =
                (a_rho() / (complex<double>(q2 - m_rho * m_rho, m_rho * hbar / tau_rho))
                 + a_omega() * complex<double>(std::cos(delta_omega_m_rho()), std::sin(delta_omega_m_rho())) / (complex<double>(q2 - m_omega * m_omega, m_omega * hbar / tau_omega))
                 + a_phi() * complex<double>(std::cos(delta_phi_m_rho()), std::sin(delta_phi_m_rho())) / (complex<double>(q2 - m_phi * m_phi, m_phi * hbar / tau_phi)))
                * complex<double>(std::cos(delta_rho()), std::sin(delta_rho()));
        const complex<double> cPR = (a_eta() / (complex<double>(q2 - m_eta * m_eta, m_eta * hbar / tau_eta))
                                     + a_eta_p() * complex<double>(std::cos(delta_eta_p_m_eta()), std::sin(delta_eta_p_m_eta()))
                                               / (complex<double>(q2 - m_eta_p * m_eta_p, m_eta_p * hbar / tau_eta_p)))
                                    * complex<double>(std::cos(delta_eta()), std::sin(delta_eta()));

        std::cout << std::fixed;
        std::cout << std::setprecision(20);
        std::cout << "f_+: " << form_factors->f_p(q2) << std::endl;
        std::cout << "f_0: " << form_factors->f_0(q2) << std::endl;
        std::cout << "f_T: " << form_factors->f_t(q2) << std::endl;
        std::cout << "c9R: " << c9R << std::endl;
        std::cout << "c9+c9': " << c9_p << std::endl;
        std::cout << "c10+c10': " << c10_p << std::endl;
        std::cout << "c7+c7': " << c7_p << std::endl;
        std::cout << "cS+cS': " << wc.cS() + wc.cSprime() << std::endl;
        std::cout << "cP+cP': " << wc.cP() + wc.cPprime() << std::endl;
        std::cout << "cT': " << wc.cT() << std::endl;
        std::cout << "cT5': " << wc.cT5() << std::endl;
        std::cout << "cPR: " << cPR << std::endl;
        std::cout << "mc(mu): " << mcatmu << std::endl;

        // cf. [BHP:2007A], Eq. (3.2), p. 3 and 4 or [BKMS:2012A] (1205.5811)
        result.F_A  = c10_p;
        result.F_T  = F_Tkin * wc.cT();
        result.F_T5 = F_Tkin * wc.cT5();
        result.F_S  = F_Skin * (wc.cS() + wc.cSprime());
        result.F_P  = F_Skin * ((wc.cP() + wc.cPprime()) + cPR / ckm_factor_ds) + m_l() * c10_p * ((m_D2 - m_P2) / q2 * (calF_time / calF_plus - 1.0) - 1.0);
        result.F_V  = c9_p + 2.0 * mcatmu * m_D / q2 * c7_p * calF_T_plus / calF_plus
                     + c9R / ckm_factor_ds
                     //+ 2.0 * m_b_PS() / m_D / xi_pseudo(q2) * (dff.calT - 16.0 * power_of<2>(M_PI) * power_of<3>(m_D()) / m_b_PS() / q2 * calH_plus)
                     + 8.0 * m_l * m_D / q2 * calF_T_plus / calF_plus * wc.cT();

        return result;
    }
} // namespace eos
