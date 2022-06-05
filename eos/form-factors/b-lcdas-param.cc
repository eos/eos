/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Philip Lüghausen
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

#include <eos/form-factors/b-lcdas-param.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/qualified-name.hh>

#include <gsl/gsl_sf_expint.h>

// TODO: EVERYTHING
#include <stdexcept>

class NotImplemented : public std::logic_error
{
public:
    NotImplemented() : std::logic_error("Function not yet implemented") { };
};

namespace eos
{
    template <>
    struct Implementation<b_lcdas::Param>
    {
        SwitchOption opt_q;

        UsedParameter lambda_B_inv;
        UsedParameter lambda_E2;
        UsedParameter lambda_H2;

        UsedParameter w0;
        UsedParameter a0;
        UsedParameter a1;
        UsedParameter a2;
        UsedParameter a3;
        UsedParameter a4;
        UsedParameter a5;
        UsedParameter a6;
        UsedParameter a7;
        UsedParameter a8;

        SwitchOption opt_gminus;

        double switch_gminus;

        inline
        QualifiedName parameter(const char * _name) const
        {
            qnp::Name name(_name);

            if (opt_q.value() == "s")
                return QualifiedName(qnp::Prefix("B_s"), name);

            return QualifiedName(qnp::Prefix("B"), name);
        }

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            w0(p[parameter("omega0").str()], u),
            a0(p[parameter("a0").str()], u),
            a1(p[parameter("a1").str()], u),
            a2(p[parameter("a2").str()], u),
            a3(p[parameter("a3").str()], u),
            a4(p[parameter("a4").str()], u),
            a5(p[parameter("a5").str()], u),
            a6(p[parameter("a6").str()], u),
            a7(p[parameter("a7").str()], u),
            a8(p[parameter("a8").str()], u),
            opt_q(o, "q", { "u", "d", "s" }, "u"),
            lambda_B_inv(p[parameter("1/lambda_B_p").str()], u),
            lambda_E2(p[parameter("lambda_E^2").str()], u),
            lambda_H2(p[parameter("lambda_H^2").str()], u),
            opt_gminus(o, "gminus", { "zero", "WW-limit" }, "WW-limit"),
            switch_gminus(1.0)
        {
            if (opt_gminus.value() == "zero")
            {
                switch_gminus = 0.0;
            }
        }

        // inline const std::array<const double, 9> get_a_vec()
        // {
        //     return std::array<const double, 9> {
        //         a0, a1, a2, a3, a4, a5, a6, a7, a8
        //     };
        // }

        inline double L0() const
        {
            throw NotImplemented();
        }

        inline double L0inc(const double & Omega) const
        {
            throw NotImplemented();
        }

        inline double Binc(const double & Omega, const double & sigma)
        {
            throw NotImplemented();
        }

        /* the inverse moment of phi_+ */
        inline double lambda_B() const {
            return 1.0 / lambda_B_inv();
        }

        /* Leading twist two-particle LCDAs */

        inline double phi_plus(const double & omega) const {
            throw NotImplemented();
        }

        inline double phi_minus(const double & omega) const {
            throw NotImplemented();
        }


        inline double phi_bar(const double & omega) const {
            throw NotImplemented();
        }

        inline double phi_bar_d1(const double & omega) const {
            throw NotImplemented();
        }


        /* Next-to-leading twist two-particle LCDAs */

        inline double g_minusWW(const double & omega) const {
            throw NotImplemented();
        }

       inline double g_minusWW_d1(const double & omega) const {
            throw NotImplemented();
        }

       inline double g_minusWW_d2(const double & omega) const {
            throw NotImplemented();
        }

        inline double g_plus(const double & omega) const {
            throw NotImplemented();
        }

        inline double g_plus_d1(const double & omega) const {
            throw NotImplemented();
        }

        inline double g_plus_d2(const double & omega) const {
            throw NotImplemented();
        }

        inline double g_bar(const double & omega) const {
            throw NotImplemented();
        }

        inline double g_bar_d1(const double & omega) const {
            throw NotImplemented();
        }

        inline double g_bar_d2(const double & omega) const {
            throw NotImplemented();
        }

        inline double g_bar_d3(const double & omega) const {
            throw NotImplemented();
        }

        /* Leading twist three-particle LCDAs */

        inline double phi_3(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double phi_4(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double phi_bar_3(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double phi_bar_4(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double phi_bar2_3(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double phi_bar2_4(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double phi_bar_bar_3(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double phi_bar_bar_4(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double psi_bar_4(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double psi_bar_bar_4(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }


        inline double chi_bar_4(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double chi_bar_bar_4(const double & omega_1, const double & omega_2) const {
            throw NotImplemented();
        }

        inline double psi_A(const double & omega, const double & xi) const {
            throw NotImplemented();
        }

        inline double psi_V(const double & omega, const double & xi) const {
            throw NotImplemented();
        }

        inline double X_A(const double & omega, const double & xi) const {
            throw NotImplemented();
        }

        inline double Y_A(const double & omega, const double & xi) const {
            throw NotImplemented();
        }

        inline double Xbar_A(const double & omega, const double & xi) const {
            throw NotImplemented();
        }

        inline double Ybar_A(const double & omega, const double & xi) const {
            throw NotImplemented();
        }
    };

namespace b_lcdas
{
    Param::~Param() = default;

    double
    Param::L0() const
    {
        return _imp->L0();
    }

    double
    Param::L0inc(const double & Omega) const
    {
        return _imp->L0inc(Omega);
    }

    double
    Param::Binc(const double & Omega, const double & sigma) const
    {
        return _imp->Binc(Omega, sigma);
    }

    Param::Param(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<Param>(new Implementation<Param>(p, o, *this))
    {
    }

    double
    Param::phi_plus(const double & omega) const
    {
        return _imp->phi_plus(omega);
    }

    double
    Param::phi_minus(const double & omega) const
    {
        return _imp->phi_minus(omega);
    }

    double
    Param::phi_bar(const double & omega) const
    {
        return _imp->phi_bar(omega);
    }

    double
    Param::phi_bar_d1(const double & omega) const
    {
        return _imp->phi_bar_d1(omega);
    }

    double
    Param::g_plus(const double & omega) const
    {
        return _imp->g_plus(omega);
    }

    double
    Param::g_plus_d1(const double & omega) const
    {
        return _imp->g_plus_d1(omega);
    }

    double
    Param::g_plus_d2(const double & omega) const
    {
        return _imp->g_plus_d2(omega);
    }

    double
    Param::g_minusWW(const double & omega) const
    {
        return _imp->g_minusWW(omega);
    }

    double
    Param::g_minusWW_d1(const double & omega) const
    {
        return _imp->g_minusWW_d1(omega);
    }

    double
    Param::g_minusWW_d2(const double & omega) const
    {
        return _imp->g_minusWW_d2(omega);
    }

    double
    Param::g_bar(const double & omega) const
    {
        return _imp->g_bar(omega);
    }

    double
    Param::g_bar_d1(const double & omega) const
    {
        return _imp->g_bar_d1(omega);
    }

    double
    Param::g_bar_d2(const double & omega) const
    {
        return _imp->g_bar_d2(omega);
    }

    double
    Param::g_bar_d3(const double & omega) const
    {
        return _imp->g_bar_d3(omega);
    }

    double
    Param::phi_3(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_3(omega_1, omega_2);
    }

    double
    Param::phi_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_4(omega_1, omega_2);
    }

    double
    Param::phi_bar_3(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar_3(omega_1, omega_2);
    }

    double
    Param::phi_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar_4(omega_1, omega_2);
    }

    double
    Param::phi_bar2_3(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar2_3(omega_1, omega_2);
    }

    double
    Param::phi_bar2_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar2_4(omega_1, omega_2);
    }

    double
    Param::phi_bar_bar_3(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar_bar_3(omega_1, omega_2);
    }

    double
    Param::phi_bar_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar_bar_4(omega_1, omega_2);
    }

    double
    Param::psi_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->psi_bar_4(omega_1, omega_2);
    }

    double
    Param::psi_bar_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->psi_bar_bar_4(omega_1, omega_2);
    }

    double
    Param::chi_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->chi_bar_4(omega_1, omega_2);
    }

    double
    Param::chi_bar_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->chi_bar_bar_4(omega_1, omega_2);
    }

    double
    Param::inverse_lambda_plus() const
    {
        return 1.0 / _imp->lambda_B();
    }

    double
    Param::psi_A(const double & omega, const double & xi) const
    {
        return _imp->psi_A(omega, xi);
    }

    double
    Param::psi_V(const double & omega, const double & xi) const
    {
        return _imp->psi_V(omega, xi);
    }

    double
    Param::X_A(const double & omega, const double & xi) const
    {
        return _imp->X_A(omega, xi);
    }

    double
    Param::Y_A(const double & omega, const double & xi) const
    {
        return _imp->Y_A(omega, xi);
    }

    double
    Param::Xbar_A(const double & omega, const double & xi) const
    {
        return _imp->Xbar_A(omega, xi);
    }

    double
    Param::Ybar_A(const double & omega, const double & xi) const
    {
        return _imp->Ybar_A(omega, xi);
    }
}
}