/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Stefan Meiser
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

#include <eos/form-factors/rho-lcdas.hh>
#include <eos/models/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    template <>
    struct Implementation<RhoLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 (even) para Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a2para_0;
        UsedParameter a4para_0;
        UsedParameter fpara;

        // twist 2 (even) perp Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a2perp_0;
        UsedParameter a4perp_0;
        UsedParameter fperp_0;

        // twist 3 LCDA parameters at mu = 1 GeV
        UsedParameter zeta3para_0;
        UsedParameter omega3paratilde_0;
        UsedParameter omega3para_0;
        UsedParameter omega3perp_0;

        // twist 4 LCDA parameters at mu = 1 GeV
        UsedParameter zeta4para_0;
        UsedParameter omega4paratilde_0;
        UsedParameter zeta4perp_0;
        UsedParameter zeta4perptilde_0;

        // mass and decay constant of the rho
        UsedParameter m_rho;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a2para_0(p["rho::a2para@1GeV"], u),
            a4para_0(p["rho::a4para@1GeV"], u),
            fpara(p["rho::fpara"], u),
            a2perp_0(p["rho::a2perp@1GeV"], u),
            a4perp_0(p["rho::a4perp@1GeV"], u),
            fperp_0(p["rho::fperp@1GeV"], u),
            zeta3para_0(p["rho::zeta3para@1GeV"], u),
            omega3paratilde_0(p["rho::omega3paratilde@1GeV"], u),
            omega3para_0(p["rho::omega3para@1GeV"], u),
            omega3perp_0(p["rho::omega3perp@1GeV"], u),
            zeta4para_0(p["rho::zeta4para@1GeV"], u),
            omega4paratilde_0(p["rho::omega4paratilde@1GeV"], u),
            zeta4perp_0(p["rho::zeta4perp@1GeV"], u),
            zeta4perptilde_0(p["rho::zeta4perptilde@1GeV"], u),
            m_rho(p["mass::rho^+"], u),
            _mu_c(p["QCD::mu_c"], u),
            _mu_b(p["QCD::mu_b"], u),
            _mu_t(p["QCD::mu_t"], u)
        {
        }

        inline double c_rge(const double & _mu) const
        {
            /*
             * RGE coefficient, basically
             *
             *     (alpha_s/alpha_s_0)^(1/beta_0),
             *
             * with matching between the individual n-flavor QCDs.
             */

            double mu = _mu, alpha_s_mu = model->alpha_s(mu);
            double mu_0 = 1.0, alpha_s_0 = model->alpha_s(mu_0);

            if (mu < _mu_c)
                return std::pow(alpha_s_mu / alpha_s_0, 1.0 / QCD::beta_function_nf_3[0]);

            double alpha_s_c = model->alpha_s(_mu_c);
            double result = std::pow(alpha_s_c / alpha_s_0, 1.0 / QCD::beta_function_nf_3[0]);

            if (mu < _mu_b)
                return result * std::pow(alpha_s_mu / alpha_s_c, 1.0 / QCD::beta_function_nf_4[0]);

            double alpha_s_b = model->alpha_s(_mu_b);
            result *= std::pow(alpha_s_b / alpha_s_c, 1.0 / QCD::beta_function_nf_4[0]);

            if (mu < _mu_t)
                return result * std::pow(alpha_s_mu / alpha_s_b, 1.0 / QCD::beta_function_nf_5[0]);

            throw InternalError("Implementation<RhoLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        /* running of twist 2 parameters */
        inline double a2para(const double & mu) const
        {
            return a2para_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double a4para(const double & mu) const
        {
            return a4para_0 * std::pow(c_rge(mu), 364.0 / 45.0);
        }

        inline double a2perp(const double & mu) const
        {
            return a2perp_0 * std::pow(c_rge(mu), 52.0 / 9.0);
        }

        inline double a4perp(const double & mu) const
        {
            return a4perp_0 * std::pow(c_rge(mu), 368.0 / 45.0);
        }

        inline double fperp(const double & mu) const
        {
            // [BBKT1998A], p. 23, eq. (3.59)
            return fperp_0 * std::pow(c_rge(mu), +4.0 / 3.0);
        }

        /* running of twist 3 parameters */
        inline double zeta3para(const double & mu) const
        {
            return zeta3para_0 * std::pow(c_rge(mu), +77.0 / 9.0);
        }
        inline double omega3paratilde(const double & mu) const
        {
            return
                (std::pow(c_rge(mu), (205 - std::sqrt(865)) / 18.0) * (6.0 * std::sqrt(865.0) * (-1.0 + std::pow(c_rge(mu), std::sqrt(865.0) / 9.0)) * omega3para_0 +
                (865.0 - 26.0 * std::sqrt(865.0) + (865.0 + 26.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), std::sqrt(865.0) / 9.0)) * omega3paratilde_0)) / 1730.0;
        }
        inline double omega3para(const double & mu) const
        {
            return
                (std::pow(c_rge(mu), (205.0 - std::sqrt(865.0)) / 18.0) * ((1730.0 + 52.0 * std::sqrt(865.0) + (1730.0 - 52.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), std::sqrt(865.0) / 9.0)) * omega3para_0 +
                63.0 * std::sqrt(865.0) * (-1.0 + std::pow(c_rge(mu), std::sqrt(865.0) / 9.0)) * omega3paratilde_0)) / 3460.0;
        }
        inline double omega3perp(const double & mu) const
        {
            return  omega3perp_0 * std::pow(c_rge(mu), +73.0 / 9.0) * fperp_0 / fperp(mu);
        }

        /* running of twist 4 parameters */
        inline double zeta4para(const double & mu) const
        {
            return zeta4para_0 * std::pow(c_rge(mu), +32.0 / 9.0);
        }
        inline double omega4paratilde(const double & mu) const
        {
            return omega4paratilde_0 * std::pow(c_rge(mu), 10.0);
        }
        inline double zeta4perp(const double & mu) const
        {
            return
                1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) + std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perp_0 +
                1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) - std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perptilde_0;
        }
        inline double zeta4perptilde(const double & mu) const
        {
            return
                1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) - std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perp_0 +
                1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) + std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perptilde_0;
        }

        // inline functions for two particle twist 3 LCDAs (only includes up to a2perp(para) while the leading twist LCDAs are implemented up to a4perp(para))
        inline double psi3para(const double & u, const double & mu) const
        {
            return 6.0 * (1.0 - u) * u * (1.0 + (a2perp(mu) / 6.0 + (5.0 * omega3perp(mu)) / 18.0) * (-1.5 + (15.0 * std::pow(-1.0 + 2.0 * u, 2)) / 2.0));
        }
        inline double phi3para(const double & u, const double & mu) const
        {
            return 3.0 * std::pow(-1.0 + 2.0 * u, 2) +
                (3.0 * a2perp(mu) * std::pow(-1.0 + 2.0 * u, 2) * (-3.0 + 5.0 * std::pow(-1.0 + 2.0 * u, 2))) / 2.0 +
                (5.0 * omega3perp(mu) * (3.0 - 30.0 * std::pow(-1.0 + 2.0 * u, 2) + 35.0 * std::pow(-1.0 + 2.0 * u, 4)))/8.0;
        }
        inline double psi3perp(const double & u, const double & mu) const
        {
            return 6.0 * (1.0 - u) * u * (1.0 + (-1.5 + (15.0 * std::pow(-1.0 + 2.0 * u, 2)) / 2.0) * (a2para(mu) / 6.0 +
                (5.0 *omega3para(mu)) / 12.0 - (5.0 * omega3paratilde(mu)) / 24.0 + (10.0 * zeta3para(mu)) / 9.0));
        }
        inline double phi3perp(const double & u, const double & mu) const
        {
            return (3.0 * (1.0 + std::pow(-1.0 + 2.0 * u, 2))) / 4.0 + ((9.0 * a2para(mu)) / 112.0 + (15.0 * omega3para(mu)) / 32.0 -
                (15.0 * omega3paratilde(mu)) / 64.0) * (3.0 - 30.0 * std::pow(-1.0 + 2.0 * u, 2) +
                35.0 * std::pow(-1.0 + 2.0 * u, 4)) + (-1.0 + 3.0 * std::pow(-1.0 + 2.0 * u, 2)) * ((3.0 * a2para(mu)) / 7.0 + 5.0 * zeta3para(mu));
        }
        // definition of chiral even parameters appearing in three-particle twist 4 LCDAs, [BBL:2007A], eq. 3.22 (renormalon model)
        inline double psi1para(const double & mu) const
        {
            return 7.0 / 12.0 * zeta4para(mu);
        }
        inline double psi0paratilde(const double & mu) const
        {
            return 1.0 / 3.0 * zeta4para(mu);
        }
        inline double psi2paratilde(const double & mu) const
        {
            return -7.0 / 12.0 * zeta4para(mu);
        }
        inline double phi0para(const double & mu) const
        {
            return 1.0 / 3.0 * zeta4para(mu);
        }
        inline double phi1para(const double & mu) const
        {
            return -7.0 / 18.0 * zeta4para(mu);
        }
        inline double phi2para(const double & mu) const
        {
            return -7.0 / 9.0 * zeta4para(mu);
        }

        // inline functions for three-particle twist 4 LCDAs
        inline double Psi4para(const double & u1, const double & u2, const double & u3, const double & mu) const
        {
            return 120.0 * u1 * u2 * u3 * psi1para(mu) * (u1 - u2);
        }
        inline double Psi4paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const
        {
            return 120.0 * u1 * u2 * u3 * (psi0paratilde(mu) + psi2paratilde(mu) * (3.0 * u3 - 1.0));
        }
        inline double Phi4para(const double & u1, const double & u2, const double & u3, const double & mu) const
        {
            return -30.0 * u3 * u3 * (u1 - u2) * (phi0para(mu) + u3 * phi1para(mu) + 1.0 / 2.0 * (5.0 * u3 - 3.0) * phi2para(mu));
        }
        inline double Phi4paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const
        {
            return 30.0 * u3 * u3 * (phi0para(mu) * (1.0 - u3) + phi1para(mu) * (u3 * (1.0 - u3) - 6.0 * u1 * u2) + phi2para(mu) * (u3 * (1.0 - u3) - 3.0 / 2.0 * (u1 * u1 + u2 * u2)));
        }

        // definition of chiral odd G conserving parameters appearing in three-particle twist 4 LCDAs, [BBL:2007A], eq. 4.13, setting zeta4perptilde = -zeta4perp (renormalon model)
        inline double psi0perp(const double & mu) const
        {
            return zeta4perp(mu);
        }
        inline double psi0perptilde(const double & mu) const
        {
            return -zeta4perp(mu);
        }

        // Twist-4 three particle chiral odd G conserving auxiliary parameters [BBL:2007A], eq. 4.17, setting twist 2, 3 and non-genuin twist 4 to zero
        inline double phi1perp(const double & mu) const
        {
            return 63.0 / 220.0 * Q1(mu) - 119.0 / 44.0 * Q3(mu);
        }
        inline double phi1perptilde(const double & mu) const
        {
            return -63.0 / 220.0 * Q1(mu) - 35.0 / 44.0 * Q3(mu);
        }
        inline double psi1perp(const double & mu) const
        {
            return 49.0 / 110.0 * Q1(mu) - 7.0 / 22.0 * Q3(mu);
        }
        inline double psi1perptilde(const double & mu) const
        {
            return -49.0 / 110.0 * Q1(mu) + 7.0 / 22.0 * Q3(mu);
        }
        inline double psi2perp(const double & mu) const
        {
            return 28.0 / 55.0 * Q1(mu) + 7.0 / 11.0 * Q3(mu);
        }
        inline double psi2perptilde(const double & mu) const
        {
            return -28.0 / 55.0 * Q1(mu) - 7.0 / 11.0 * Q3(mu);
        }

        // definition of chiral odd G conserving parameters appearing in three-particle twist 4 LCDAs, [BBL:2007A], eq. 4.19 (renormalon model)
        inline double Q1(const double & mu) const
        {
            return -10.0 / 3.0 * zeta4perp(mu);
        }
        inline double Q3(const double & mu) const
        {
            return -1.0 * zeta4perp(mu);
        }

        // inline functions for chriral odd three-particle twist 4 LCDAs
        inline double Psi4perp(const double & u1, const double & u2, const double & u3, const double & mu) const
        {
            return 30.0 * u3 * u3 * (psi0perp(mu) * (1.0 - u3) + psi1perp(mu) * (u3 * (1.0 - u3) - 6.0 * u1 * u2) +
                psi2perp(mu) * (u3 * (1.0 - u3) - 3.0 / 2.0 * (u1 * u1 + u2 * u2)));
        }
        inline double Psi4perptilde(const double & u1, const double & u2, const double & u3, const double & mu) const
        {
            return 30.0 * u3 * u3 * (psi0perptilde(mu) * (1.0 - u3) + psi1perptilde(mu) * (u3 * (1.0 - u3) - 6.0 * u1 * u2) +
                psi2perptilde(mu) * (u3 * (1.0 - u3) - 3.0 / 2.0 * (u1 * u1 + u2 * u2)));
        }
        inline double Phi4perp1(const double & u1, const double & u2, const double & u3, const double & mu) const
        {
            return 120.0 * u1 * u2 * u3 * phi1perp(mu) * (u1 - u2);
        }
        inline double Phi4perp2(const double & u1, const double & u2, const double & u3, const double & mu) const
        {
            return 30.0 * u3 * u3 * (u1 - u2) * (psi0perptilde(mu) + u3 * psi1perptilde(mu) + 1.0 / 2.0 * (5.0 * u3 - 3.0) * psi2perptilde(mu));
        }
        inline double Phi4perp3(const double & u1, const double & u2, const double & u3, const double & mu) const
        {
            return -120.0 * u1 * u2 * u3 * phi1perptilde(mu) * (u1 - u2);
        }
        inline double Phi4perp4(const double & u1, const double & u2, const double & u3, const double & mu) const
        {
            return -30.0 * u3 * u3 * (u1 - u2) * (psi0perp(mu) + u3 * psi1perp(mu) + 1.0 / 2.0 * (5.0 * u3 - 3.0) * psi2perp(mu));
        }

        // inline functions for chriral even two-particle twist 4 LCDAs
        inline double psi4paraT4(const double & u, const double & mu) const
        {
            static const GegenbauerPolynomial gp_2_1o2(2, 1.0 / 2.0);
            const double x = 2.0 * u - 1.0;
            const double c2 = gp_2_1o2.evaluate(x);

            return -20.0 / 3.0 * zeta4para(mu) * c2;
        }
        inline double psi4paraWW(const double & u, const double & mu) const
        {
            static const GegenbauerPolynomial gp_2_1o2(2, 1.0 / 2.0);
            static const GegenbauerPolynomial gp_4_1o2(4, 1.0 / 2.0);
            const double x = 2.0 * u - 1.0;
            const double c2 = gp_2_1o2.evaluate(x);
            const double c4 = gp_4_1o2.evaluate(x);

            return 1.0 +
                (-1.0 - 2.0 / 7.0 * a2para(mu) + 40.0 /3.0 * zeta3para(mu)) * c2 +
                (-27.0 / 28.0 * a2para(mu) + 5.0 / 4.0 * zeta3para(mu) - 15.0 / 8.0 * omega3para(mu) - 15.0 / 16.0 * omega3paratilde(mu)) * c4;
        }
        inline double phi4paraT4(const double & u, const double & mu) const
        {
            const double x = 2.0 * u - 1.0;

            return 30.0 * u * u * (1.0 - u) * (1.0 - u) * 20.0 / 9.0 * zeta4para(mu) -
                84.0 * omega4paratilde(mu) * ((1.0 / 8.0 * u * (1.0 - u)) * (21.0 - 13.0 * x * x) +
                power_of<3>(u) * (10.0 - 15.0 * u + 6.0 * u * u) * log(u) +
                power_of<3>(1.0 - u) * (10.0 - 15.0 * (1.0 - u) + 6.0 * power_of<2>(1.0 - u)) * log(1.0 - u));
        }
        inline double phi4paraWW(const double & u, const double & mu) const
        {
            static const GegenbauerPolynomial gp_2_5o2(2, 5.0 / 2.0);
            const double x = 2.0 * u - 1.0;
            const double c2 = gp_2_5o2.evaluate(x);

            return 30.0 * u * u * (1.0 - u) * (1.0 - u) * (4.0 / 5.0 * (1.0 + 1.0 / 21.0 * a2para(mu) + 10.0 / 9.0 * zeta3para(mu)) +
                1.0 / 10.0 * (9.0 / 7.0 * a2para(mu) + 1.0 / 9.0 * zeta3para(mu) + 7.0 / 6.0 * omega3para(mu) - 3.0 / 4.0 * omega3paratilde(mu)) * c2) +
                2.0 * (-2.0 * a2para(mu) - 14.0 / 3.0 * zeta3para(mu) + 3.0 * omega3para(mu)) * (1.0 / 8.0 * u * (1.0 - u) * (21.0 - 13.0 * x * x) +
                power_of<3>(u) * (10.0 - 15.0 * u + 6.0 * u * u) * log(u) + power_of<3>(1.0 - u) * (10.0 - 15.0 * (1.0 - u) + 6.0 * power_of<2>(1.0 - u)) * log(1.0 - u));
        }
    };

    RhoLCDAs::RhoLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<RhoLCDAs>(new Implementation<RhoLCDAs>(p, o, *this)),
        gp_2_3o2(2, 3.0 / 2.0),
        gp_4_3o2(4, 3.0 / 2.0)
    {
    }

    RhoLCDAs::~RhoLCDAs()
    {
    }

    VectorLCDAs *
    RhoLCDAs::make(const Parameters & p, const Options & o)
    {
        return new RhoLCDAs(p, o);
    }

    double
    RhoLCDAs::a2para(const double & mu) const
    {
        return _imp->a2para(mu);
    }

    double
    RhoLCDAs::a4para(const double & mu) const
    {
        return _imp->a4para(mu);
    }

    double
    RhoLCDAs::fpara() const
    {
        return _imp->fpara();
    }
    double
    RhoLCDAs::a2perp(const double & mu) const
    {
        return _imp->a2perp(mu);
    }

    double
    RhoLCDAs::a4perp(const double & mu) const
    {
        return _imp->a4perp(mu);
    }

    double
    RhoLCDAs::fperp(const double &mu) const
    {
        return _imp->fperp(mu);
    }

    double
    RhoLCDAs::zeta3para(const double &mu) const
    {
        return _imp->zeta3para(mu);
    }

    double
    RhoLCDAs::omega3paratilde(const double &mu) const
    {
        return _imp->omega3paratilde(mu);
    }

    double
    RhoLCDAs::omega3para(const double &mu) const
    {
        return _imp->omega3para(mu);
    }

    double
    RhoLCDAs::omega3perp(const double &mu) const
    {
        return _imp->omega3perp(mu);
    }

    double
    RhoLCDAs::zeta4para(const double &mu) const
    {
        return _imp->zeta4para(mu);
    }

    double
    RhoLCDAs::omega4paratilde(const double &mu) const
    {
        return _imp->omega4paratilde(mu);
    }

    double
    RhoLCDAs::zeta4perp(const double &mu) const
    {
        return _imp->zeta4perp(mu);
    }

    double
    RhoLCDAs::zeta4perptilde(const double &mu) const
    {
        return _imp->zeta4perptilde(mu);
    }

    double
    RhoLCDAs::phi2para(const double & u, const double & mu) const
    {
        const double x = 2.0 * u - 1.0;
        const double c2 = gp_2_3o2.evaluate(x);
        const double c4 = gp_4_3o2.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a2para(mu) * c2 + _imp->a4para(mu) * c4);
    }

    double
    RhoLCDAs::phi2perp(const double & u, const double & mu) const
    {
        const double x = 2.0 * u - 1.0;
        const double c2 = gp_2_3o2.evaluate(x);
        const double c4 = gp_4_3o2.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a2perp(mu) * c2 + _imp->a4perp(mu) * c4);
    }

    double
    RhoLCDAs::psi3para(const double & u, const double & mu) const
    {
        return _imp->psi3para(u, mu);
    }

    double
    RhoLCDAs::phi3para(const double & u, const double & mu) const
    {
        return _imp->phi3para(u, mu);
    }

    double
    RhoLCDAs::psi3perp(const double & u, const double & mu) const
    {
        return _imp->psi3perp(u, mu);
    }

    double
    RhoLCDAs::phi3perp(const double & u, const double & mu) const
    {
        return _imp->phi3perp(u, mu);
    }

    double
    RhoLCDAs::Phi3para(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return 360.0 * u1 * u2 * u3 * u3 * _imp->omega3para(mu) * (u1 - u2);
    }

    double
    RhoLCDAs::Phi3paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return 360.0 * u1 * u2 * u3 * u3 * (_imp->zeta3para(mu) + _imp->omega3paratilde(mu) * 1.0 / 2.0 * (7.0 * u3 - 3.0));
    }

    double
    RhoLCDAs::Phi3perp(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return 360.0 * u1 * u2 * u3 * u3 * _imp->omega3perp(mu) * (u1 - u2);
    }

    double
    RhoLCDAs::Psi4para(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return _imp->Psi4para(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Psi4paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return _imp->Psi4paratilde(u1, u2, u3, mu);
    }
    double
    RhoLCDAs::Phi4para(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return _imp->Phi4para(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Phi4paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return _imp->Phi4paratilde(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Psi4perp(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return _imp->Psi4perp(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Psi4perptilde(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return _imp->Psi4perptilde(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Phi4perp1(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return _imp->Phi4perp1(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Phi4perp2(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return _imp->Phi4perp2(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Phi4perp3(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return _imp->Phi4perp3(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Phi4perp4(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return _imp->Phi4perp4(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Phi4perptilde1(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return -Phi4perp3(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Phi4perptilde2(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return -Phi4perp4(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Phi4perptilde3(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return -Phi4perp1(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::Phi4perptilde4(const double & u1, const double & u2, const double & u3, const double & mu) const
    {
        return -Phi4perp2(u1, u2, u3, mu);
    }

    double
    RhoLCDAs::psi4para(const double & u, const double & mu) const
    {
        return _imp->psi4paraT4(u, mu) + _imp->psi4paraWW(u, mu);
    }

    double
    RhoLCDAs::phi4para(const double & u, const double & mu) const
    {
        return _imp->phi4paraT4(u, mu) + _imp->phi4paraWW(u, mu);
    }

    Diagnostics
    RhoLCDAs::diagnostics() const
    {
        Diagnostics results;

        results.add(Diagnostics::Entry{ _imp->c_rge(1.0), "RGE coefficient C(mu = 1.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(2.0), "RGE coefficient C(mu = 2.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(3.0), "RGE coefficient C(mu = 3.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(4.0), "RGE coefficient C(mu = 4.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(5.0), "RGE coefficient C(mu = 5.0 GeV)" });

        return results;
    }
}
