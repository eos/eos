/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#include <eos/form-factors/form-factors.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/rare-b-decays/b-to-kstar-gamma-naive.hh>

#include <cmath>

namespace eos
{
    BToKstarGammaAmplitudes<tag::Naive>::BToKstarGammaAmplitudes(const Parameters & p, const Options & o) :
        AmplitudeGenerator(p, o)
    {
        Context ctx("When constructing B->K^*gamma Naive amplitudes");
    }

    BToKstarGamma::Amplitudes
    BToKstarGammaAmplitudes<tag::Naive>::amplitudes() const
    {
        const WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), LeptonFlavor::muon /*fake lepton flavor*/, cp_conjugate);

        // cf. [BFS:2001A], below Eq. (9), p. 4
        const complex<double> c7eff = wc.c7() - 1.0 / 3.0 * wc.c3() - 4.0 / 9.0 * wc.c4() - 20.0 / 3.0 * wc.c5() - 80.0 / 9.0 * wc.c6();

        const double ff_T1 = form_factors->t_1(0.0);

        const double calN = std::sqrt(alpha_e * power_of<3>(m_B) * power_of<3>(1.0 - m_Kstar * m_Kstar / (m_B * m_B)) / (32.0 * power_of<4>(M_PI))) * g_fermi
                            * model->m_b_msbar(mu()) * std::abs(model->ckm_tb() * conj(model->ckm_ts()));

        // cf. [BFS:2001A], Eq. (15), p. 5, at order alpha_s^0 and with the full tensor form factor in place of xi_perp
        const complex<double> a_left  = complex<double>(0.0, +1.0) * calN * ff_T1 * c7eff;
        const complex<double> a_right = complex<double>(0.0, -1.0) * calN * ff_T1 * wc.c7prime();

        return BToKstarGamma::Amplitudes{
            (a_left - a_right) / std::sqrt(2.0),
            (a_left + a_right) / std::sqrt(2.0),
        };
    }
} // namespace eos
