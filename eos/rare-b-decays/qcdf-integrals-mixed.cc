/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2012, 2016 Danny van Dyk
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

#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/nonlocal-form-factors/hard-scattering.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/rare-b-decays/qcdf-integrals-impl.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <functional>
#include <limits>

#include <iostream>

namespace eos
{
    using namespace std::placeholders;

    /* photon final state */

    // massless case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::photon_massless_case(const double & m_B,
            const double & m_V, const double & mu,
            const double & a_1_perp, const double & a_2_perp,
            const double & a_1_para, const double & a_2_para)
    {
        return QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::photon_massless_case(m_B, m_V, mu, a_1_perp, a_2_perp, a_1_para, a_2_para);
    }

    // charm case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::photon_charm_case(const double & m_c,
            const double & m_B, const double & m_V, const double & mu,
            const double & a_1_perp, const double & a_2_perp,
            const double & a_1_para, const double & a_2_para)
    {
        return QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::photon_charm_case(m_c, m_B, m_V, mu, a_1_perp, a_2_perp, a_1_para, a_2_para);
    }

    // bottom case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::photon_bottom_case(const double & m_b,
            const double & m_B, const double & m_V, const double & mu,
            const double & a_1_perp, const double & a_2_perp,
            const double & a_1_para, const double & a_2_para)
    {
        return QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::photon_bottom_case(m_b, m_B, m_V, mu, a_1_perp, a_2_perp, a_1_para, a_2_para);
    }


    /* dilepton final states */

    // massless case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_massless_case(const double & s,
            const double & m_B, const double & m_V, const double & mu,
            const double & a_1_perp, const double & a_2_perp,
            const double & a_1_para, const double & a_2_para)
    {
        return QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_massless_case(s, m_B, m_V, mu, a_1_perp, a_2_perp, a_1_para, a_2_para);
    }

    // charm case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_charm_case(const double & s,
            const double & m_c, const double & m_B, const double & m_V, const double & mu,
            const double & a_1_perp, const double & a_2_perp,
            const double & a_1_para, const double & a_2_para)
    {
        QCDFIntegrals<BToKstarDilepton> results;
        results = QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_charm_case(s, m_c, m_B, m_V, mu, a_1_perp, a_2_perp, a_1_para, a_2_para);

        // avoid NaN at u=1
        static const double u_min = 0.0 + 1e-5;
        static const double u_max = 1.0 - 1e-5;

        // we only need to correct integral J_1
        std::function<complex<double> (const double &)> j_1_perp = std::bind(&HardScattering::j1, s, _1, m_c, m_B, a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_1_para = std::bind(&HardScattering::j1, s, _1, m_c, m_B, a_1_para, a_2_para);
        cubature::Config cub_conf = cubature::Config().epsrel(1e-3);
        results.j1_perp     = integrate(j_1_perp, u_min, u_max, cub_conf);
        results.j1_parallel = integrate(j_1_para, u_min, u_max, cub_conf);

        // composite results
        const double sh = s / m_B / m_B;
        const double eh = (1.0 + power_of<2>(m_V / m_B) - sh) / 2.0;
        results.jtilde1_perp = 2.0 / eh * results.j1_perp + sh * results.j2_perp / (eh * eh);
        results.jtilde2_parallel = 2.0 / eh * results.j1_parallel + results.j3_parallel / (eh * eh);

        return results;
    }

    // bottom case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_bottom_case(const double & s,
            const double & m_b, const double & m_B, const double & m_V, const double & mu,
            const double & a_1_perp, const double & a_2_perp,
            const double & a_1_para, const double & a_2_para)
    {
        return QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_bottom_case(s, m_b, m_B, m_V, mu, a_1_perp, a_2_perp, a_1_para, a_2_para);
    }
}
