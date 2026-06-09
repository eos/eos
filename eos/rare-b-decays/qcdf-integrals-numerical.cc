/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2012, 2016 Danny van Dyk
 * Copyright (c) 2025 Florian Herren
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

#include <eos/maths/integrate-impl.hh>
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

    // massless case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::photon_massless_case(const double & /*m_B*/,
            const double & /*m_V*/, const double & /*mu*/,
            const double & /*a_1_perp*/, const double & /*a_2_perp*/,
            const double & /*a_1_para*/, const double & /*a_2_para*/)
    {
        throw InternalError("QCDFIntegralCalculator::photon_massless_case: Numerical integration of photon cases not supported");
    }

    // charm case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::photon_charm_case(const double & /*m_c*/,
            const double & /*m_B*/, const double & /*m_V*/, const double & /*mu*/,
            const double & /*a_1_perp*/, const double & /*a_2_perp*/,
            const double & /*a_1_para*/, const double & /*a_2_para*/)
    {
        throw InternalError("QCDFIntegralCalculator::photon_charm_case: Numerical integration of photon cases not supported");
    }

    // bottom case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::photon_bottom_case(const double & /*m_b*/,
            const double & /*m_B*/, const double & /*m_V*/, const double & /*mu*/,
            const double & /*a_1_perp*/, const double & /*a_2_perp*/,
            const double & /*a_1_para*/, const double & /*a_2_para*/)
    {
        throw InternalError("QCDFIntegralCalculator::photon_bottom_case: Numerical integration of photon cases not supported");
    }

    // massless case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_massless_case(const double & s,
            const double & m_B, const double & m_V, const double & mu,
            const double & a_1_perp, const double & a_2_perp,
            const double & a_1_para, const double & a_2_para)
    {
        QCDFIntegrals<BToKstarDilepton> results;

        // avoid NaN at u=1
        static const double u_min = 0.0 + 1e-5;
        static const double u_max = 1.0 - 1e-5;
        // We use the same regularising cut-off x ~= Lambda / m_B as in j7_zero as to ensure
        // a smooth transition B->K^*ll -> B->K^*gamma for s -> 0. (Lambda = 0.5 / GeV).
        // The relative error for j7  in the QCDF region 1 <= q^2 <= 6 is less than 25%.
        // Since j7 enters only via subleading terms, it amounts to a relative error of A_FB
        // in the SM of < 0.3%.
        static const double u_max_7 = 1.0 - 0.5 / m_B;

        cubature::Config cub_conf = cubature::Config().epsrel(1e-3);

        // perpendicular amplitude
        std::function<double (const double &)>          j_0_perp    = std::bind(&HardScattering::j0, s, _1, m_B, a_1_perp, a_2_perp);
        std::function<double (const double &)>          j_0bar_perp = std::bind(&HardScattering::j0, s, _1, m_B, -a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_1_perp    = std::bind(&HardScattering::j1, s, _1, 0.0, m_B, a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_2_perp    = std::bind(&HardScattering::j2_massless, s, _1, m_B, a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_4_perp    = std::bind(&HardScattering::j4, s, _1, 0.0, m_B, mu, a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_5_perp    = std::bind(&HardScattering::j5, s, _1, 0.0, m_B, mu, a_1_perp, a_2_perp);
        // This integral arises in perpendicular amplitudes, but depends on parallel Gegenbauer moments!
        std::function<complex<double> (const double &)> j_6_perp    = std::bind(&HardScattering::j6, s, _1, 0.0, m_B, mu, a_1_para, a_2_para);
        std::function<double (const double &)>          j_7_perp    = std::bind(&HardScattering::j7, s, _1, m_B, a_1_perp, a_2_perp);
        results.j0_perp    = integrate<1>(j_0_perp,    u_min, u_max, cub_conf);
        results.j0bar_perp = integrate<1>(j_0bar_perp, u_min, u_max, cub_conf);
        results.j1_perp    = integrate<1, 1, complex<double>>(j_1_perp,    u_min, u_max, cub_conf);
        results.j2_perp    = integrate<1, 1, complex<double>>(j_2_perp,    u_min, u_max, cub_conf);
        results.j4_perp    = integrate<1, 1, complex<double>>(j_4_perp,    u_min, u_max, cub_conf);
        results.j5_perp    = integrate<1, 1, complex<double>>(j_5_perp,    u_min, u_max, cub_conf);
        results.j6_perp    = integrate<1, 1, complex<double>>(j_6_perp,    u_min, u_max, cub_conf);
        results.j7_perp    = integrate<1>(j_7_perp,    u_min, u_max_7, cub_conf);

        // parallel amplitude
        std::function<double (const double &)>          j_0_para = std::bind(&HardScattering::j0, s, _1, m_B, a_1_para, a_2_para);
        std::function<complex<double> (const double &)> j_1_para = std::bind(&HardScattering::j1, s, _1, 0.0, m_B, a_1_para, a_2_para);
        std::function<complex<double> (const double &)> j_3_para = std::bind(&HardScattering::j3_massless, s, _1, m_B, a_1_para, a_2_para);
        std::function<complex<double> (const double &)> j_4_para = std::bind(&HardScattering::j4, s, _1, 0.0, m_B, mu, a_1_para, a_2_para);
        results.j0_parallel = integrate<1>(j_0_para, u_min, u_max, cub_conf);
        results.j1_parallel = integrate<1, 1, complex<double>>(j_1_para, u_min, u_max, cub_conf);
        results.j3_parallel = integrate<1, 1, complex<double>>(j_3_para, u_min, u_max, cub_conf);
        results.j4_parallel = integrate<1, 1, complex<double>>(j_4_para, u_min, u_max, cub_conf);

        // composite results
        const double sh = s / m_B / m_B;
        const double eh = (1.0 + power_of<2>(m_V / m_B) - sh) / 2.0;
        results.jtilde1_perp = 2.0 / eh * results.j1_perp + sh * results.j2_perp / (eh * eh);
        results.jtilde2_parallel = 2.0 / eh * results.j1_parallel + results.j3_parallel / (eh * eh);

        return results;
    }

    // charm case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_charm_case(const double & s,
            const double & m_c, const double & m_B, const double & m_V, const double & mu,
            const double & a_1_perp, const double & a_2_perp,
            const double & a_1_para, const double & a_2_para)
    {
        QCDFIntegrals<BToKstarDilepton> results;

        // avoid NaN at u=1
        static const double u_min = 0.0 + 1e-5;
        static const double u_max = 1.0 - 1e-5;
        // We use the same regularising cut-off x ~= Lambda / m_B as in j7_zero as to ensure
        // a smooth transition B->K^*ll -> B->K^*gamma for s -> 0. (Lambda = 0.5 / GeV).
        // The relative error for j7  in the QCDF region 1 <= q^2 <= 6 is less than 25%.
        // Since j7 enters only via subleading terms, it amounts to a relative error of A_FB
        // in the SM of < 0.3%.
        static const double u_max_7 = 1.0 - 0.5 / m_B;

        cubature::Config cub_conf = cubature::Config().epsrel(1e-3);

        // perpendicular amplitude
        std::function<double (const double &)>          j_0_perp    = std::bind(&HardScattering::j0, s, _1, m_B, a_1_perp, a_2_perp);
        std::function<double (const double &)>          j_0bar_perp = std::bind(&HardScattering::j0, s, _1, m_B, -a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_1_perp    = std::bind(&HardScattering::j1, s, _1, m_c, m_B, a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_2_perp    = std::bind(&HardScattering::j2, s, _1, m_c, m_B, a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_4_perp    = std::bind(&HardScattering::j4, s, _1, m_c, m_B, mu, a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_5_perp    = std::bind(&HardScattering::j5, s, _1, m_c, m_B, mu, a_1_perp, a_2_perp);
        // This integral arises in perpendicular amplitudes, but depends on parallel Gegenbauer moments!
        std::function<complex<double> (const double &)> j_6_perp    = std::bind(&HardScattering::j6, s, _1, m_c, m_B, mu, a_1_para, a_2_para);
        std::function<double (const double &)>          j_7_perp    = std::bind(&HardScattering::j7, s, _1, m_B, a_1_perp, a_2_perp);
        results.j0_perp    = integrate<1>(j_0_perp,    u_min, u_max, cub_conf);
        results.j0bar_perp = integrate<1>(j_0bar_perp, u_min, u_max, cub_conf);
        results.j1_perp    = integrate<1, 1, complex<double>>(j_1_perp,    u_min, u_max, cub_conf);
        results.j2_perp    = integrate<1, 1, complex<double>>(j_2_perp,    u_min, u_max, cub_conf);
        results.j4_perp    = integrate<1, 1, complex<double>>(j_4_perp,    u_min, u_max, cub_conf);
        results.j5_perp    = integrate<1, 1, complex<double>>(j_5_perp,    u_min, u_max, cub_conf);
        results.j6_perp    = integrate<1, 1, complex<double>>(j_6_perp,    u_min, u_max, cub_conf);
        results.j7_perp    = integrate<1>(j_7_perp,    u_min, u_max_7, cub_conf);

        // parallel amplitude
        std::function<double (const double &)>          j_0_para = std::bind(&HardScattering::j0, s, _1, m_B, a_1_para, a_2_para);
        std::function<complex<double> (const double &)> j_1_para = std::bind(&HardScattering::j1, s, _1, m_c, m_B, a_1_para, a_2_para);
        std::function<complex<double> (const double &)> j_3_para = std::bind(&HardScattering::j3, s, _1, m_c, m_B, a_1_para, a_2_para);
        std::function<complex<double> (const double &)> j_4_para = std::bind(&HardScattering::j4, s, _1, m_c, m_B, mu, a_1_para, a_2_para);
        results.j0_parallel = integrate<1>(j_0_para, u_min, u_max, cub_conf);
        results.j1_parallel = integrate<1, 1, complex<double>>(j_1_para, u_min, u_max, cub_conf);
        results.j3_parallel = integrate<1, 1, complex<double>>(j_3_para, u_min, u_max, cub_conf);
        results.j4_parallel = integrate<1, 1, complex<double>>(j_4_para, u_min, u_max, cub_conf);

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
    QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_bottom_case(const double & s,
            const double & m_b, const double & m_B, const double & m_V, const double & mu,
            const double & a_1_perp, const double & a_2_perp,
            const double & a_1_para, const double & a_2_para)
    {
        QCDFIntegrals<BToKstarDilepton> results;

        // avoid NaN at u=1
        static const double u_min = 0.0 + 1e-5;
        static const double u_max = 1.0 - 1e-5;
        // We use the same regularising cut-off x ~= Lambda / m_B as in j7_zero as to ensure
        // a smooth transition B->K^*ll -> B->K^*gamma for s -> 0. (Lambda = 0.5 / GeV).
        // The relative error for j7  in the QCDF region 1 <= q^2 <= 6 is less than 25%.
        // Since j7 enters only via subleading terms, it amounts to a relative error of A_FB
        // in the SM of < 0.3%.
        static const double u_max_7 = 1.0 - 0.5 / m_B;

        cubature::Config cub_conf = cubature::Config().epsrel(1e-3);

        // perpendicular amplitude
        std::function<double (const double &)>          j_0_perp    = std::bind(&HardScattering::j0, s, _1, m_B, a_1_perp, a_2_perp);
        std::function<double (const double &)>          j_0bar_perp = std::bind(&HardScattering::j0, s, _1, m_B, -a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_1_perp    = std::bind(&HardScattering::j1, s, _1, m_b, m_B, a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_2_perp    = std::bind(&HardScattering::j2, s, _1, m_b, m_B, a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_4_perp    = std::bind(&HardScattering::j4, s, _1, m_b, m_B, mu, a_1_perp, a_2_perp);
        std::function<complex<double> (const double &)> j_5_perp    = std::bind(&HardScattering::j5, s, _1, m_b, m_B, mu, a_1_perp, a_2_perp);
        // This integral arises in perpendicular amplitudes, but depends on parallel Gegenbauer moments!
        std::function<complex<double> (const double &)> j_6_perp    = std::bind(&HardScattering::j6, s, _1, m_b, m_B, mu, a_1_para, a_2_para);
        std::function<double (const double &)>          j_7_perp    = std::bind(&HardScattering::j7, s, _1, m_B, a_1_perp, a_2_perp);
        results.j0_perp    = integrate<1>(j_0_perp,    u_min, u_max, cub_conf);
        results.j0bar_perp = integrate<1>(j_0bar_perp, u_min, u_max, cub_conf);
        results.j1_perp    = integrate<1, 1, complex<double>>(j_1_perp,    u_min, u_max, cub_conf);
        results.j2_perp    = integrate<1, 1, complex<double>>(j_2_perp,    u_min, u_max, cub_conf);
        results.j4_perp    = integrate<1, 1, complex<double>>(j_4_perp,    u_min, u_max, cub_conf);
        results.j5_perp    = integrate<1, 1, complex<double>>(j_5_perp,    u_min, u_max, cub_conf);
        results.j6_perp    = integrate<1, 1, complex<double>>(j_6_perp,    u_min, u_max, cub_conf);
        results.j7_perp    = integrate<1>(j_7_perp,    u_min, u_max_7, cub_conf);

        // parallel amplitude
        std::function<double (const double &)>          j_0_para = std::bind(&HardScattering::j0, s, _1, m_B, a_1_para, a_2_para);
        std::function<complex<double> (const double &)> j_1_para = std::bind(&HardScattering::j1, s, _1, m_b, m_B, a_1_para, a_2_para);
        std::function<complex<double> (const double &)> j_3_para = std::bind(&HardScattering::j3, s, _1, m_b, m_B, a_1_para, a_2_para);
        std::function<complex<double> (const double &)> j_4_para = std::bind(&HardScattering::j4, s, _1, m_b, m_B, mu, a_1_para, a_2_para);
        results.j0_parallel = integrate<1>(j_0_para, u_min, u_max, cub_conf);
        results.j1_parallel = integrate<1, 1, complex<double>>(j_1_para, u_min, u_max, cub_conf);
        results.j3_parallel = integrate<1, 1, complex<double>>(j_3_para, u_min, u_max, cub_conf);
        results.j4_parallel = integrate<1, 1, complex<double>>(j_4_para, u_min, u_max, cub_conf);

        // composite results
        const double sh = s / m_B / m_B;
        const double eh = (1.0 + power_of<2>(m_V / m_B) - sh) / 2.0;
        results.jtilde1_perp = 2.0 / eh * results.j1_perp + sh * results.j2_perp / (eh * eh);
        results.jtilde2_parallel = 2.0 / eh * results.j1_parallel + results.j3_parallel / (eh * eh);

        return results;
    }
}
