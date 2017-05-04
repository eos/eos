/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
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

#include <eos/form-factors/b-lcdas.hh>
#include <eos/utils/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

namespace eos
{
    template <>
    struct Implementation<BMesonLCDAs>
    {
        UsedParameter lambda_B_p;

        Implementation(const Parameters & p, ParameterUser & u) :
            lambda_B_p(p["lambda_B_p"], u)
        {
        }

        inline double phi_plus(const double & omega) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B_p;

            return omega / (omega_0 * omega_0) * std::exp(-omega / omega_0);
        }

        inline double phi_minus(const double & omega) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B_p;

            return 1.0 / omega_0 * std::exp(-omega / omega_0);
        }

        inline double Phibar(const double & omega) const
        {
            const double omega_0 = lambda_B_p;

            return -omega / omega_0 * std::exp(-omega / omega_0);
        }

        inline double psi_A(const double & omega, const double & xi) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B_p, omega_0_2 = omega_0 * omega_0, omega_0_4 = omega_0_2 * omega_0_2;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            return lambda_E_2 / (6.0 * omega_0_4) * xi * xi * std::exp(-(omega + xi) / omega_0);
        }

        inline double psi_V(const double & omega, const double & xi) const
        {
            return psi_A(omega, xi);
        }

        inline double X_A(const double & omega, const double & xi) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B_p, omega_0_2 = omega_0 * omega_0, omega_0_4 = omega_0_2 * omega_0_2;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            return lambda_E_2 / (6.0 * omega_0_4) * xi * (2.0 * omega - xi) * std::exp(-(omega + xi) / omega_0);
        }

        inline double Y_A(const double & omega, const double & xi) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B_p, omega_0_2 = omega_0 * omega_0, omega_0_4 = omega_0_2 * omega_0_2;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            return -lambda_E_2 / (24.0 * omega_0_4) * xi * (7.0 * omega_0 - 13.0 * omega + 3.0 * xi) * std::exp(-(omega + xi) / omega_0);
        }

        inline double Xbar_A(const double & omega, const double & xi) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B_p, omega_0_2 = omega_0 * omega_0, omega_0_3 = omega_0_2 * omega_0;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            // obtained by analytica integrating Y_A(tau, xi) over 0 <= tau <= omega.
            return lambda_E_2 / (6.0 * omega_0_3) * xi * std::exp(-(xi + omega) / omega_0)
                * (xi - 2.0 * (omega + omega_0) + std::exp(omega / omega_0) * (2.0 * omega_0 - xi));
        }

        inline double Ybar_A(const double & omega, const double & xi) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B_p, omega_0_2 = omega_0 * omega_0, omega_0_3 = omega_0_2 * omega_0;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            // obtained by analytica integrating Y_A(tau, xi) over 0 <= tau <= omega.
            return -lambda_E_2 / (24.0 * omega_0_3) * xi * std::exp(-(xi + omega) / omega_0)
                * (-3.0 * xi + 13.0 * omega + 6.0 * omega_0 + 3.0 * std::exp(omega / omega_0) * (xi - 2.0 * omega_0));
        }
    };

    BMesonLCDAs::BMesonLCDAs(const Parameters & p, const Options &) :
        PrivateImplementationPattern<BMesonLCDAs>(new Implementation<BMesonLCDAs>(p, *this))
    {
    }

    BMesonLCDAs::~BMesonLCDAs() = default;

    double
    BMesonLCDAs::phi_plus(const double & omega) const
    {
        return _imp->phi_plus(omega);
    }

    double
    BMesonLCDAs::phi_minus(const double & omega) const
    {
        return _imp->phi_minus(omega);
    }

    double
    BMesonLCDAs::Phibar(const double & omega) const
    {
        return _imp->Phibar(omega);
    }

    double
    BMesonLCDAs::inverse_lambda_plus() const
    {
        return 1.0 / _imp->lambda_B_p();
    }

    double
    BMesonLCDAs::psi_A(const double & omega, const double & xi) const
    {
        return _imp->psi_A(omega, xi);
    }

    double
    BMesonLCDAs::psi_V(const double & omega, const double & xi) const
    {
        return _imp->psi_V(omega, xi);
    }

    double
    BMesonLCDAs::X_A(const double & omega, const double & xi) const
    {
        return _imp->X_A(omega, xi);
    }

    double
    BMesonLCDAs::Y_A(const double & omega, const double & xi) const
    {
        return _imp->Y_A(omega, xi);
    }

    double
    BMesonLCDAs::Xbar_A(const double & omega, const double & xi) const
    {
        return _imp->Xbar_A(omega, xi);
    }

    double
    BMesonLCDAs::Ybar_A(const double & omega, const double & xi) const
    {
        return _imp->Ybar_A(omega, xi);
    }
}
