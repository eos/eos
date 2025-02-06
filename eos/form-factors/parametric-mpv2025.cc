/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Danny van Dyk
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

#include <eos/form-factors/parametric-mpv2025.hh>
#include <eos/maths/derivative.hh>
#include <eos/utils/exception.hh>
#include <eos/maths/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/models/model.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

#include <functional>

namespace eos
{
    template <>
    struct Implementation<MPV2025FormFactors>
    {
        // hadronic parameters
        UsedParameter m_D;
        UsedParameter m_K;
        UsedParameter m_pi;

        // Ds pole parameters
        UsedParameter m_Ds1908;     // D_s^+(1908)        with J^P = 0^-
        UsedParameter m_Dsst2008;   // D_s^{*+}(2008)     with J^P = 1^-
        UsedParameter m_Dszst2317;  // D_{s,0}^{*+}(2317) with J^P = 0^+
        UsedParameter m_Dszst2460;  // D_{s,0}^{*+}(2460) with J^P = 1^+

        // K-pi resonance parameters
        // - S wave resonances
        UsedParameter m_Kzst700;
        UsedParameter gamma_Kzst700;
        UsedParameter m_Kzst1430;
        UsedParameter gamma_Kzst1430;
        // - P wave resonances
        UsedParameter m_Kzst892;
        UsedParameter gamma_Kzst892;
        UsedParameter m_Kzst1410;
        UsedParameter gamma_Kzst1410;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            m_D(p["mass::D_u"], u),
            m_K(p["mass::K_d"], u),
            m_pi(p["mass::pi^-"], u),
            m_Ds1908(p["mass::D_s(1908)@MPV2025"], u),
            m_Dsst2008(p["mass::D_s^*(2008)@MPV2025"], u),
            m_Dszst2317(p["mass::D_s0^*(2317)@MPV2025"], u),
            m_Dszst2460(p["mass::D_s0^*(2460)@MPV2025"], u),
            m_Kzst700(p["mass::K^*(700)@MPV2025"], u),
            gamma_Kzst700(p["decay-width::K^*(1410)@MPV2025"], u),
            m_Kzst1430(p["mass::K^*(1430)@MPV2025"], u),
            gamma_Kzst1430(p["decay-width::K^*(1430)@MPV2025"], u),
            m_Kzst892(p["mass::K^*(892)@MPV2025"], u),
            gamma_Kzst892(p["decay-width::K^*(892)@MPV2025"], u),
            m_Kzst1410(p["mass::K^*(1410)@MPV2025"], u),
            gamma_Kzst1410(p["decay-width::K^*(1410)@MPV2025"], u)
        {
        }

        // 1st z expansion, evaluation on the first Riemann sheet
        inline complex<double> zq2(const complex<double> & q2) const
        {
            static const complex<double> tplus = power_of<2>(1.9683 + 0.1396);
            static const complex<double> tzero = 0.0;

            return (std::sqrt(tplus - q2) - std::sqrt(tplus - tzero)) / (std::sqrt(tplus - q2) + std::sqrt(tplus - tzero));
        }

        // 2nd z expansion, evaluation on the first Riemann sheet
        inline complex<double> zk2(const complex<double> & k2) const
        {
            const complex<double> tplus = power_of<2>(m_K + m_pi);
            const complex<double> tzero = -2.0;

            return (std::sqrt(tplus - k2) - std::sqrt(tplus - tzero)) / (std::sqrt(tplus - k2) + std::sqrt(tplus - tzero));
        }

        inline complex<double> pole_factor_S(const complex<double> & k2) const
        {
            const auto z = zk2(k2);
            const complex<double> sr1 = power_of<2>(complex<double>(m_Kzst700, -0.5 * gamma_Kzst700));
            const complex<double> zr1 = 1.0 / zk2(sr1);
            const complex<double> sr2 = power_of<2>(complex<double>(m_Kzst1430, -0.5 * gamma_Kzst1430));
            const complex<double> zr2 = 1.0 / zk2(sr2);

            return (z - zr1) * (z - conj(zr1)) * (z - zr2) * (z - conj(zr2));
        }

        inline complex<double> pole_factor_P(const complex<double> & k2) const
        {
            const auto z = zk2(k2);
            const complex<double> sr1 = power_of<2>(complex<double>(m_Kzst892, -0.5 * gamma_Kzst892));
            const complex<double> zr1 = 1.0 / zk2(sr1);
            const complex<double> sr2 = power_of<2>(complex<double>(m_Kzst1410, -0.5 * gamma_Kzst1410));
            const complex<double> zr2 = 1.0 / zk2(sr2);

            return (z - zr1) * (z - conj(zr1)) * (z - zr2) * (z - conj(zr2));
        }

        complex<double> f_perp(const double & q2, const double & k2, const double & z) const
        {
            // no S-wave term
            // P-wave term
            auto f_perp_P = 1.0 / (q2 - power_of<2>(m_Dsst2008)) / pole_factor_P(k2);

            return sqrt(3.0 / 2.0) * f_perp_P;
        }

        complex<double> f_para(const double & q2, const double & k2, const double & z) const
        {
            // no S-wave term
            // P-wave term
            auto f_para_P = 1.0 / (q2 - power_of<2>(m_Dszst2460)) / pole_factor_P(k2);

            return sqrt(3.0 / 2.0) * f_para_P;
        }

        complex<double> f_long(const double & q2, const double & k2, const double & z) const
        {
            // S-wave term
            auto f_long_S = 1.0 / (q2 - power_of<2>(m_Dszst2317)) / pole_factor_S(k2);
            // P-wave term
            auto f_long_P = 1.0 / (q2 - power_of<2>(m_Dszst2317)) / pole_factor_P(k2);

            return f_long_S + sqrt(3.0) * f_long_P * z;
        }

        complex<double> f_time(const double & q2, const double & k2, const double & z) const
        {
            // S-wave term
            auto f_time_S = 1.0 / (q2 - power_of<2>(m_Ds1908)) / pole_factor_S(k2);
            // P-wave term
            auto f_time_P = 1.0 / (q2 - power_of<2>(m_Ds1908)) / pole_factor_P(k2);

            return f_time_S + sqrt(3.0) * f_time_P * z;
        }

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            return results;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<MPV2025FormFactors>::options
    {
    };

    MPV2025FormFactors::MPV2025FormFactors(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<MPV2025FormFactors>(new Implementation<MPV2025FormFactors>(p, o, *this))
    {
    }

    MPV2025FormFactors::~MPV2025FormFactors()
    {
    }

    FormFactors<PToPP> *
    MPV2025FormFactors::make(const Parameters & p, const Options & o)
    {
        return new MPV2025FormFactors(p, o);
    }

    complex<double>
    MPV2025FormFactors::f_perp(const double & q2, const double & k2, const double & z) const
    {
        return _imp->f_perp(q2, k2, z);
    }

    complex<double>
    MPV2025FormFactors::f_para(const double & q2, const double & k2, const double & z) const
    {
        return _imp->f_para(q2, k2, z);
    }

    complex<double>
    MPV2025FormFactors::f_long(const double & q2, const double & k2, const double & z) const
    {
        return _imp->f_long(q2, k2, z);
    }

    complex<double>
    MPV2025FormFactors::f_time(const double & q2, const double & k2, const double & z) const
    {
        return _imp->f_time(q2, k2, z);
    }

    Diagnostics
    MPV2025FormFactors::diagnostics() const
    {
        return _imp->diagnostics();
    }

    std::vector<OptionSpecification>::const_iterator
    MPV2025FormFactors::begin_options()
    {
        return Implementation<MPV2025FormFactors>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    MPV2025FormFactors::end_options()
    {
        return Implementation<MPV2025FormFactors>::options.cend();
    }
}
