/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Stephan Kürten
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

#include <eos/form-factors/parametric-kkvdz2022.hh>

using namespace std::literals::complex_literals;

namespace eos
{
    KKvDZ2022FormFactors::KKvDZ2022FormFactors(const Parameters & p, const Options &) :
        _a_omega{{{ UsedParameter(p[_par_name("omega_1_0")],  *this),
                   UsedParameter(p[_par_name("omega_1_1")],  *this),
                   UsedParameter(p[_par_name("omega_1_2")],  *this) },
                 { UsedParameter(p[_par_name("omega_2_0")],  *this),
                   UsedParameter(p[_par_name("omega_2_1")],  *this),
                   UsedParameter(p[_par_name("omega_2_2")],  *this) },
                 { UsedParameter(p[_par_name("omega_3_0")],  *this),
                   UsedParameter(p[_par_name("omega_3_1")],  *this),
                   UsedParameter(p[_par_name("omega_3_2")],  *this) },
                 { UsedParameter(p[_par_name("omega_4_0")],  *this),
                   UsedParameter(p[_par_name("omega_4_1")],  *this),
                   UsedParameter(p[_par_name("omega_4_2")],  *this) }}},
        _a_rho{{{ UsedParameter(p[_par_name("rho_1_0")],  *this),
                 UsedParameter(p[_par_name("rho_1_1")],  *this),
                 UsedParameter(p[_par_name("rho_1_2")],  *this) },
               { UsedParameter(p[_par_name("rho_2_0")],  *this),
                 UsedParameter(p[_par_name("rho_2_1")],  *this),
                 UsedParameter(p[_par_name("rho_2_2")],  *this) },
               { UsedParameter(p[_par_name("rho_3_0")],  *this),
                 UsedParameter(p[_par_name("rho_3_1")],  *this),
                 UsedParameter(p[_par_name("rho_3_2")],  *this) },
               { UsedParameter(p[_par_name("rho_4_0")],  *this),
                 UsedParameter(p[_par_name("rho_4_1")],  *this),
                 UsedParameter(p[_par_name("rho_4_2")],  *this) }}}
    {
    }

    FormFactors<PToGammaOffShell> *
    KKvDZ2022FormFactors::make(const Parameters & p, const Options & o)
    {
        return new KKvDZ2022FormFactors(p, o);
    }

    std::string
    KKvDZ2022FormFactors::_par_name(const std::string & ff_name)
    {
        return std::string("B->gamma^*") + std::string("::N^") + ff_name + std::string("@KKvDZ2022");
    }

    double
    KKvDZ2022FormFactors::_width_omega(const double & q2) const
    {
        const double constant_omega_width_approximation = 0.00868;
        const double m_pion = BToPi::m_P, m_pion_sq = m_pion * m_pion;

        if (q2 - 9.0 * m_pion_sq > 0.0)
        {
            return constant_omega_width_approximation;
        }
        if (q2 - 9.0 * m_pion_sq <= 0.0)
        {
            return 0.0;
        }
        else
        {
            throw InternalError("invalid value for q2");
        }
    }

    double
    KKvDZ2022FormFactors::_width_rho(const double & q2) const
    {
        const double constant_rho_width_approximation = 0.1474;
        const double m_pion = BToPi::m_P, m_pion_sq = m_pion * m_pion;
        const double m_rho = BToRho::m_V, m_rho_sq = m_rho * m_rho;

        if (q2 - 4.0 * m_pion_sq > 0.0)
        {
            return std::pow((q2 - 4.0 * m_pion_sq)/(m_rho_sq - 4.0 * m_pion_sq), 1.5) * m_rho_sq/q2 * constant_rho_width_approximation;
        }
        if (q2 - 4.0 * m_pion_sq <= 0.0)
        {
            return 0.0;
        }
        else
        {
            throw InternalError("invalid value for q2");
        }
    }

    double
    KKvDZ2022FormFactors::_z_omega(const double & k2) const
    {
        const double m_B = BToOmega::m_B, m_omega = BToOmega::m_V;
        const double t_p_omega = power_of<2>(m_B + m_omega), t_m_omega = power_of<2>(m_B - m_omega);
        const double t_0_omega = t_p_omega*(1.0 - sqrt(1.0 - t_m_omega/t_p_omega));

        return (sqrt(t_p_omega - k2) - sqrt(t_p_omega - t_0_omega))/(sqrt(t_p_omega - k2) + sqrt(t_p_omega - t_0_omega));
    }

    double
    KKvDZ2022FormFactors::_z_rho(const double & k2) const
    {
        const double m_B = BToOmega::m_B, m_rho = BToRho::m_V;
        const double t_p_rho = power_of<2>(m_B + m_rho), t_m_rho = power_of<2>(m_B - m_rho);
        const double t_0_rho = t_p_rho*(1.0 - sqrt(1.0 - t_m_rho/t_p_rho));

        return (sqrt(t_p_rho - k2) - sqrt(t_p_rho - t_0_rho))/(sqrt(t_p_rho - k2) + sqrt(t_p_rho - t_0_rho));
    }

    complex<double>
    KKvDZ2022FormFactors::_calc_ff_contribution_omega(const double & q2, const double & k2, const double & m_r_sq, const std::array<UsedParameter, 3> & a) const
    {
        double a_0(a[0]),a_1(a[1]),a_2(a[2]);
        const double diff_z = _z_omega(k2) - _z_omega(0.0);
        const double m_omega = BToOmega::m_V, m_omega_sq = m_omega * m_omega;

        return 1.0 / (1.0 - k2 / m_r_sq) * ( m_omega_sq / ( m_omega_sq - q2 - 1.0i * m_omega * _width_omega(q2) ) ) *
                (a_0 + a_1 * diff_z + a_2 * power_of<2>(diff_z) );
    }

    complex<double>
    KKvDZ2022FormFactors::_calc_ff_contribution_rho(const double & q2, const double & k2, const double & m_r_sq, const std::array<UsedParameter, 3> & a) const
    {
        double a_0(a[0]),a_1(a[1]),a_2(a[2]);
        const double diff_z = _z_rho(k2) - _z_rho(0.0);
        const double m_rho = BToRho::m_V, m_rho_sq = m_rho * m_rho;

        return 1.0 / (1.0 - k2 / m_r_sq) * ( m_rho_sq / ( m_rho_sq - q2 - 1.0i * sqrt(q2) * _width_rho(q2) ) ) *
                (a_0 + a_1 * diff_z + a_2 * power_of<2>(diff_z) );
    }

    complex<double>
    KKvDZ2022FormFactors::F_1(const double & q2, const double & k2) const
    {
        const double m_B1_sq = BToOmega::mR2_1p;

        return _calc_ff_contribution_omega(q2, k2, m_B1_sq, _a_omega[0]) + _calc_ff_contribution_rho(q2, k2, m_B1_sq, _a_rho[0]);
    }

    complex<double>
    KKvDZ2022FormFactors::F_2(const double & q2, const double & k2) const
    {
        const double m_B1_sq = BToOmega::mR2_1p;

        return _calc_ff_contribution_omega(q2, k2, m_B1_sq, _a_omega[1]) + _calc_ff_contribution_rho(q2, k2, m_B1_sq, _a_rho[1]);
    }

    complex<double>
    KKvDZ2022FormFactors::F_3(const double & q2, const double & k2) const
    {
        const double m_B = BToOmega::m_B;

        return _calc_ff_contribution_omega(q2, k2, m_B * m_B, _a_omega[2]) + _calc_ff_contribution_rho(q2, k2, m_B * m_B, _a_rho[2]);
    }

    complex<double>
    KKvDZ2022FormFactors::F_4(const double & q2, const double & k2) const
    {
        const double m_Bstar_sq = BToOmega::mR2_1m;

        return _calc_ff_contribution_omega(q2, k2, m_Bstar_sq, _a_omega[3]) + _calc_ff_contribution_rho(q2, k2, m_Bstar_sq, _a_rho[3]);
    }

    const std::set<ReferenceName>
    KKvDZ2022FormFactors::references
    {
        "KKvDZ:2022A"_rn
    };

    const std::vector<OptionSpecification> KKvDZ2022FormFactors::options
    {
    };

    std::vector<OptionSpecification>::const_iterator
    KKvDZ2022FormFactors::begin_options()
    {
        return options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    KKvDZ2022FormFactors::end_options()
    {
        return options.cend();
    }
}
