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

#include <eos/form-factors/analytic-b-to-kstar.hh>
#include <eos/form-factors/b-lcdas.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/polylog.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

#include <functional>

#include <gsl/gsl_sf_expint.h>

#include <iostream>

namespace eos
{
    template <>
    struct Implementation<AnalyticFormFactorBToKstarKMO2006>
    {
        std::shared_ptr<Model> model;

        // masses
        UsedParameter m_B;
        UsedParameter m_Kstar;

        // hadronic parameters
        UsedParameter f_B;
        UsedParameter f_Kstar;
        UsedParameter lambda_B_p;

        // sum rule parameters
        UsedParameter s0;
        UsedParameter M2;

        // renormalization scale
        UsedParameter mu;

        BMesonLCDAs b_lcdas;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            m_B(p["mass::B_d"], u),
            m_Kstar(p["mass::K^*_d"], u),
            f_B(p["decay-constant::B_d"], u),
            f_Kstar(p["K^*::f_para"], u),
            lambda_B_p(p["lambda_B_p"], u),
            s0(p["B->K^*::s_0@LCSR"], u),
            M2(p["B->K^*::M^2@LCSR"], u),
            mu(p["mu"], u),
            //mu(p["B->K^*::mu@KMO2006"], u),
            b_lcdas(p, o)
        {
            u.uses(b_lcdas);
        }

        double sigma0(const double & q2) const
        {
            const auto m_B2 = std::pow(m_B(), 2);

            return (1.0 + s0 / m_B2 - q2 / m_B2) / 2.0 - std::sqrt(std::pow(1 - s0 / m_B2 + q2 / m_B2, 2) / 4.0 - q2 / m_B2);
        }

        inline double etaf(const double & q2, const double & sigma) const
        {
            return 1.0 / (1.0 - q2 / pow((1.0 - sigma) * m_B(), 2));
        }

        inline double Detaf(const double & q2, const double & sigma) const
        {
            return 2.0 * q2 / (m_B() * m_B() * pow(1.0 - sigma, 3) * pow(1.0 - q2 / pow((1.0 - sigma) * m_B(), 2), 2));
        }

        /* V */

        inline double v_iota1(const double & /*q2*/, const double & /*sigma*/) const
        {
            return 0.0;
        }

        inline double v_iota2(const double & /*q2*/, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-8)
                return 0.0; // error ~ -5e-9

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return -(exp(-xi)*(-24*gamma*lambdaB2 +
                24*lambdaB*mB*sigma - 8*gamma*lambdaB*mB*sigma - 2*mB2*sigma +
                3*mB2*sigma2 + 24*lambdaB2*exp(xi)* Ei -
                16*lambdaB*mB*sigma*exp(xi)* Ei +
                4*mB2*sigma2*exp(xi)* Ei + 24*lambdaB2*log(lambdaB)
                + 8*lambdaB*mB*sigma* log(lambdaB) - 24*lambdaB2*log(mB*sigma) -
                8*lambdaB*mB*sigma* log(mB*sigma) + 12*lambdaB2*exp(xi)*log(xi) -
                8*lambdaB*mB*sigma*exp(xi)* log(xi) + 12*lambdaB2*exp(xi)*
                log(lambdaB*pow(mB,-1)* pow(sigma,-1)) -
                8*lambdaB*mB*sigma*exp(xi)* log(lambdaB*pow(mB,-1)*
                pow(sigma,-1)))*pow(mB2,-1)* pow(-1 + sigma,-3))/4.;
        }

        inline double v_iota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-5)
                return 0.0; // error ~ -7e-10

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return (mB2 - q2 - 2*mB2*sigma + mB2*sigma2)* exp(-xi)*
                (-24*gamma* lambdaB2 +
                24*lambdaB*mB*sigma - 8*gamma*lambdaB*mB*sigma + mB2*sigma2 +
                24*lambdaB2*exp(xi)* Ei -
                16*lambdaB*mB*sigma*exp(xi)* Ei +
                4*mB2*sigma2*exp(xi)* Ei + 24*lambdaB2*log(lambdaB)
                + 8*lambdaB*mB*sigma* log(lambdaB) + 24*lambdaB2*exp(xi)*
                log(lambdaB) - 16*lambdaB*mB*sigma*exp(xi)* log(lambdaB) -
                24*lambdaB2*log(mB*sigma) - 8*lambdaB*mB*sigma* log(mB*sigma) -
                24*lambdaB2*exp(xi)* log(mB*sigma) + 16*lambdaB*mB*sigma*exp(xi)*
                log(mB*sigma) + 24*lambdaB2*exp(xi)*log(xi) -
                16*lambdaB*mB*sigma*exp(xi)* log(xi)) / (4.0 * mB2 * pow(1 - sigma,4));
        }

        inline double v_Diota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-9)
                return 0.0; // error ~ 6e-8

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2), lambdaB3 = pow(lambdaB, 3);
            const auto mB = m_B(), mB2 = pow(mB, 2), mB3 = pow(mB, 3);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return (exp(-xi)*(-(mB3/lambdaB) -  (2*mB2*(2*(5 - 2*gamma)*
            lambdaB + mB - 4*lambdaB*log(xi)))/ (lambdaB*(-1 + sigma)) +
            (q2*(96*lambdaB2 - 48*gamma*lambdaB2 + 26*lambdaB*mB -
            8*gamma*lambdaB*mB + mB2 + 24*Ei*lambdaB* (-2*lambdaB +
            mB)*exp(xi)\
             - 12*lambdaB2*exp(xi)* log(lambdaB) - 12*lambdaB2*exp(xi)*
               log(1/(mB*sigma)) - 48*lambdaB2*log(xi) - 8*lambdaB*mB*log(xi) -
               12*lambdaB2*exp(xi)*log(xi)))/ (lambdaB*mB* pow(-1 + sigma,4))
               - (mB*(48*lambdaB2 - 32*gamma*lambdaB2 + 22*lambdaB*mB -
                 8*gamma*lambdaB*mB + mB2 - q2 + 8*Ei*lambdaB* (-2*lambdaB +
                 mB)*exp(xi)\
             - 4*lambdaB2*exp(xi)* log(lambdaB) - 4*lambdaB2*exp(xi)*
               log(1/(mB*sigma)) - 32*lambdaB2*log(xi) - 8*lambdaB*mB*log(xi) -
               4*lambdaB2*exp(xi)*log(xi)))/ (lambdaB*pow(-1 + sigma,2)) +
               (4*q2* (-24*gamma*lambdaB2 + 24*lambdaB*mB - 8*gamma*lambdaB*mB
               + mB2 + 4*Ei* (6*lambdaB2 - 4*lambdaB*mB + mB2)* exp(xi) -
               4*lambdaB* (-3*lambdaB + mB)*exp(xi)* log(lambdaB) +
               12*lambdaB2*exp(xi)* log(1/(mB*sigma)) - 4*lambdaB*mB*exp(xi)*
               log(1/(mB*sigma)) - 24*lambdaB2*log(xi) - 8*lambdaB*mB*log(xi) +
               12*lambdaB2*exp(xi)* log(xi) - 4*lambdaB*mB*exp(xi)* log(xi)))/
               (pow(mB,2)* pow(-1 + sigma,5)) + (2*(24*gamma*lambdaB3 -
               24*lambdaB2*mB + 8*gamma*lambdaB2*mB - lambdaB*mB2 +
               11*lambdaB*q2 - 4*gamma*lambdaB*q2 + mB*q2 - 4*Ei*lambdaB*
               (6*lambdaB2 - 4*lambdaB*mB + mB2 - q2) *exp(xi) + 4*lambdaB2*
               (-3*lambdaB + mB)*exp(xi)* log(lambdaB) - 12*lambdaB3*exp(xi)*
               log(1/(mB*sigma)) + 4*lambdaB2*mB*exp(xi)* log(1/(mB*sigma)) +
               24*lambdaB3*log(xi) + 8*lambdaB2*mB*log(xi) -
               4*lambdaB*q2*log(xi) - 12*lambdaB3*exp(xi)* log(xi) +
               4*lambdaB2*mB*exp(xi)* log(xi)))/ (lambdaB*pow(-1 +
               sigma,3)))) /4.;
        }

        double v_integrand(const double & q2, const double & sigma) const
        {
            const auto m_B2 = std::pow(m_B(), 2);

            // two-particle contribution
            const auto result_2p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * b_lcdas.phi_plus(sigma * m_B()) / (1.0 - sigma);

            // three-particle contributions
            const auto iota1 = v_iota1(q2, sigma);
            const auto iota2 = v_iota2(q2, sigma);
            const auto iota3 = v_iota3(q2, sigma);

            const auto result_3p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (-iota1 + iota2 / M2 - iota3 / (2.0 * M2 * M2));

            return result_2p + m_B() * result_3p;
        }

        inline double v(const double & q2) const
        {
            const auto m_B2 = pow(m_B(), 2);
            const auto m_Kstar2 = std::pow(m_Kstar(), 2);

            const std::function<double (const double &)> integrand = std::bind(&Implementation<AnalyticFormFactorBToKstarKMO2006>::v_integrand, this, q2, std::placeholders::_1);

            const auto sigma0  = this->sigma0(q2);
            const auto iota20  = v_iota2(q2, sigma0);
            const auto iota30  = v_iota3(q2, sigma0);
            const auto Diota30 = v_Diota3(q2, sigma0);
            const auto etaf0   = etaf(q2, sigma0);
            const auto Detaf0  = Detaf(q2, sigma0);

            const auto integral = integrate<GSL::QNG>(integrand, 0.0, sigma0);
            const auto delta = (iota20 - (1.0 / M2 + Detaf0 / m_B2) / 2.0 * iota30 - etaf0 / (2.0 * m_B2) * Diota30)
                * std::exp(-s0 / M2) / m_B2 * etaf0;

            return f_B * m_B * (m_B + m_Kstar) / (2.0 * f_Kstar * m_Kstar) * std::exp(m_Kstar2 / M2)
                * (integral + m_B * delta);
        }

        /* A_0 */

        inline double a_0_iota1(const double & /*q2*/, const double & /*sigma*/) const
        {
            return 0.0;
        }

        inline double a_0_iota2(const double & /*q2*/, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return exp(-xi)*((pow(mB,-1)*pow(-1 + sigma,-3)*
                (40*gamma*lambdaB2 + 40*gamma*lambdaB2*sigma - 40*lambdaB*mB*sigma + 18*gamma*lambdaB*mB*sigma - 40*lambdaB*mB*sigma2 +
                18*gamma*lambdaB*mB*sigma2 + 8*mB2*sigma2 - 40*lambdaB2*log(lambdaB) - 40*lambdaB2*sigma*log(lambdaB) - 18*lambdaB*mB*sigma*log(lambdaB)
                - 18*lambdaB*mB*sigma2*log(lambdaB) + 40*lambdaB2*log(mB*sigma) + 40*lambdaB2*sigma*log(mB*sigma) + 18*lambdaB*mB*sigma*log(mB*sigma) +
                18*lambdaB*mB*sigma2*log(mB*sigma) - 16*mB2*pow(sigma,3)))/8. - (exp(xi)*pow(mB,-1)*pow(-1 + sigma,-3)* (40*lambdaB2*Ei +
                40*lambdaB2*sigma*Ei - 22*lambdaB*mB*sigma*Ei - 22*lambdaB*mB*sigma2*Ei +
                2*mB2*sigma2*Ei + 2*mB2*sigma*sigma2*Ei + 40*lambdaB2*log(lambdaB) + 40*lambdaB2*sigma*log(lambdaB) -
                22*lambdaB*mB*sigma*log(lambdaB) - 22*lambdaB*mB*sigma2*log(lambdaB) - mB2*sigma*sigma2*log(lambdaB) - 40*lambdaB2*log(mB*sigma) -
                40*lambdaB2*sigma*log(mB*sigma) + 22*lambdaB*mB*sigma*log(mB*sigma) + 22*lambdaB*mB*sigma2*log(mB*sigma) + 40*lambdaB2*log(xi) +
                40*lambdaB2*sigma*log(xi) - 22*lambdaB*mB*sigma*log(xi) - 22*lambdaB*mB*sigma2*log(xi) + mB2*log(lambdaB)*pow(sigma,3)))/8.);
        }

        inline double a_0_iota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-4)
                return 0.0; // error ~ -1e-10

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return exp(-xi)*(-(sigma*(20*gamma*lambdaB2 - 20*lambdaB*mB*sigma + 9*gamma*lambdaB*mB*sigma - 2*mB2*sigma2 - 20*lambdaB2*log(lambdaB) - 9*lambdaB*mB*sigma*log(lambdaB) +
                20*lambdaB2*log(mB*sigma) + 9*lambdaB*mB*sigma*log(mB*sigma))* pow(mB,-1)*pow(-1 + sigma,-4)* (-q2 + mB2*pow(-1 + sigma,2)))/2. -
                (sigma*exp(xi)*(-20*lambdaB2* Ei + 11*lambdaB*mB*sigma*Ei - mB2*sigma2*Ei -
                20*lambdaB2*log(lambdaB) + 11*lambdaB*mB*sigma*log(lambdaB) + 20*lambdaB2*log(mB*sigma) -
                11*lambdaB*mB*sigma*log(mB*sigma) - 20*lambdaB2*log(xi) + 11*lambdaB*mB*sigma*log(xi))*pow(mB,-1)* pow(-1 +
                sigma,-4)* (-q2 + mB2*pow(-1 + sigma,2)))/2.);
        }

        inline double a_0_Diota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -4e-10

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);
            return exp(-xi)*(-(exp(xi)*pow(lambdaB,-1)*pow(mB,-1)*
                pow(-1 + sigma,-5)* (-20*lambdaB*lambdaB2*q2* Ei - 60*lambdaB*lambdaB2*q2*sigma* Ei +
                22*lambdaB2*mB*q2*sigma* Ei + 22*lambdaB2*mB*q2*sigma2* Ei - 3*lambdaB*mB2*q2*sigma2* Ei
                - lambdaB*mB2*q2*sigma*sigma2* Ei - 20*lambdaB*lambdaB2*q2*log(lambdaB) - 60*lambdaB*lambdaB2*q2*sigma* log(lambdaB) +
                22*lambdaB2*mB*q2*sigma*log(lambdaB) + 22*lambdaB2*mB*q2*sigma2*log(lambdaB) - 22*lambdaB2*mB*q2*sigma*log(mB*sigma) -
                22*lambdaB2*mB*q2*sigma2*log(mB*sigma) + 22*lambdaB2*mB*q2*sigma*log(xi) + 22*lambdaB2*mB*q2*sigma2*log(xi) -
                20*mB2*log(mB*sigma)*pow(lambdaB,3) + 20*q2*log(mB*sigma)*pow(lambdaB,3) + 20*mB2*sigma*log(mB*sigma)* pow(lambdaB,3) +
                60*q2*sigma*log(mB*sigma)* pow(lambdaB,3) + 20*mB2*sigma2*log(mB*sigma)* pow(lambdaB,3) + 20*mB2*log(xi)*pow(lambdaB,3) -
                20*q2*log(xi)*pow(lambdaB,3) - 20*mB2*sigma*log(xi)*pow(lambdaB,3) - 60*q2*sigma*log(xi)*pow(lambdaB,3) -
                20*mB2*sigma2*log(xi)*pow(lambdaB,3) + 22*lambdaB2*sigma*log(mB*sigma)* pow(mB,3) - 44*lambdaB2*sigma2*log(mB*sigma)* pow(mB,3) -
                22*lambdaB2*sigma*log(xi)*pow(mB,3) + 44*lambdaB2*sigma2*log(xi)*pow(mB,3) + 20*lambdaB*lambdaB2*mB2* Ei*pow(-1 +
                    sigma,2) + 20*lambdaB*lambdaB2*mB2*sigma* Ei*pow(-1 + sigma,2) - 22*lambdaB2*mB*mB2*sigma* Ei*pow(-1
                        + sigma,2) + 20*lambdaB*lambdaB2*mB2*log(lambdaB)* pow(-1 + sigma,2) + 20*lambdaB*lambdaB2*mB2*sigma* log(lambdaB)*pow(-1 +
                            sigma,2) - 22*lambdaB2*mB*mB2*sigma*log(lambdaB)* pow(-1 + sigma,2) + 3*lambdaB*sigma2*Ei* pow(mB2,2)*pow(-1
                                + sigma,2) - lambdaB*sigma*sigma2*Ei* pow(mB2,2)*pow(-1 + sigma,2) - 20*mB2*log(mB*sigma)*pow(lambdaB,3)*
                pow(sigma,3) + 20*mB2*log(xi)*pow(lambdaB,3)* pow(sigma,3) + 22*lambdaB2*log(mB*sigma)*pow(mB,3)* pow(sigma,3) -
                22*lambdaB2*log(xi)*pow(mB,3)* pow(sigma,3)))/2. + (pow(lambdaB,-1)*pow(mB,-1)*pow(-1 + sigma,-5)* (20*lambdaB2*mB*q2*sigma +
                2*gamma*lambdaB2*mB*q2*sigma + 60*lambdaB2*mB*q2*sigma2 - 38*gamma*lambdaB2*mB*q2*sigma2 - 13*lambdaB*mB2*q2*sigma2 +
                9*gamma*lambdaB*mB2*q2*sigma2 + 20*lambdaB*lambdaB2*q2*log(lambdaB) + 60*lambdaB*lambdaB2*q2*sigma* log(lambdaB) -
                2*lambdaB2*mB*q2*sigma*log(lambdaB) + 38*lambdaB2*mB*q2*sigma2*log(lambdaB) - 9*lambdaB*mB2*q2*sigma2*log(lambdaB) +
                9*lambdaB*mB2*q2*sigma*sigma2* log(lambdaB) + 2*lambdaB2*mB*q2*sigma*log(mB*sigma) - 38*lambdaB2*mB*q2*sigma2*log(mB*sigma) +
                9*lambdaB*mB2*q2*sigma2*log(mB*sigma) + 20*gamma*mB2*pow(lambdaB,3) - 20*gamma*q2*pow(lambdaB,3) -
                20*gamma*mB2*sigma*pow(lambdaB,3) - 60*gamma*q2*sigma*pow(lambdaB,3) - 20*gamma*mB2*sigma2*pow(lambdaB,3) +
                20*mB2*log(mB*sigma)*pow(lambdaB,3) - 20*q2*log(mB*sigma)*pow(lambdaB,3) - 20*mB2*sigma*log(mB*sigma)* pow(lambdaB,3) -
                60*q2*sigma*log(mB*sigma)* pow(lambdaB,3) - 20*mB2*sigma2*log(mB*sigma)* pow(lambdaB,3) - 20*lambdaB2*sigma*pow(mB,3) -
                2*gamma*lambdaB2*sigma*pow(mB,3) + 20*lambdaB2*sigma2*pow(mB,3) + 24*gamma*lambdaB2*sigma2*pow(mB,3) -
                2*lambdaB2*sigma*log(mB*sigma)* pow(mB,3) + 24*lambdaB2*sigma2*log(mB*sigma)* pow(mB,3) + 13*lambdaB*sigma2*pow(mB,4) -
                9*gamma*lambdaB*sigma2*pow(mB,4) - 9*lambdaB*sigma2*log(mB*sigma)* pow(mB,4) - 20*lambdaB*lambdaB2*mB2*log(lambdaB)* pow(-1 +
                    sigma,2) - 20*lambdaB*lambdaB2*mB2*sigma* log(lambdaB)*pow(-1 + sigma,2) + 2*lambdaB2*mB*mB2*sigma*log(lambdaB)* pow(-1 + sigma,2) -
                20*lambdaB2*mB*mB2*sigma2*log(lambdaB)* pow(-1 + sigma,2) + 9*lambdaB*sigma2*log(lambdaB)*pow(mB2,2)* pow(-1 + sigma,2) -
                9*lambdaB*sigma*sigma2*log(lambdaB)* pow(mB2,2)*pow(-1 + sigma,2) + 21*lambdaB*mB2*q2*pow(sigma,3) - 9*gamma*lambdaB*mB2*q2*
                pow(sigma,3) - 9*lambdaB*mB2*q2*log(mB*sigma)* pow(sigma,3) + 20*gamma*mB2*pow(lambdaB,3)* pow(sigma,3) +
                20*mB2*log(mB*sigma)*pow(lambdaB,3)* pow(sigma,3) + 20*lambdaB2*pow(mB,3)*pow(sigma,3) - 42*gamma*lambdaB2*pow(mB,3)* pow(sigma,3) -
                2*q2*pow(mB,3)*pow(sigma,3) - 42*lambdaB2*log(mB*sigma)*pow(mB,3)* pow(sigma,3) - 43*lambdaB*pow(mB,4)*pow(sigma,3) +
                27*gamma*lambdaB*pow(mB,4)* pow(sigma,3) + 27*lambdaB*log(mB*sigma)*pow(mB,4)* pow(sigma,3) + 2*pow(mB,5)*pow(sigma,3) -
                20*lambdaB2*pow(mB,3)*pow(sigma,4) + 20*gamma*lambdaB2*pow(mB,3)* pow(sigma,4) + 2*q2*pow(mB,3)*pow(sigma,4) +
                20*lambdaB2*log(mB*sigma)*pow(mB,3)* pow(sigma,4) + 47*lambdaB*pow(mB,4)*pow(sigma,4) - 27*gamma*lambdaB*pow(mB,4)* pow(sigma,4) -
                27*lambdaB*log(mB*sigma)*pow(mB,4)* pow(sigma,4) - 6*pow(mB,5)*pow(sigma,4) - 17*lambdaB*pow(mB,4)*pow(sigma,5) +
                9*gamma*lambdaB*pow(mB,4)* pow(sigma,5) + 9*lambdaB*log(mB*sigma)*pow(mB,4)* pow(sigma,5) + 6*pow(mB,5)*pow(sigma,5) -
                2*pow(mB,5)*pow(sigma,6)))/2.);
        }

        double a_0_integrand(const double & q2, const double & sigma) const
        {
            const auto m_B2 = std::pow(m_B(), 2);
            const auto phi_m = b_lcdas.phi_minus(sigma * m_B());
            const auto Phi   = b_lcdas.Phibar(sigma * m_B());

            // two-particle contribution
            const auto result_2p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (phi_m * sigma - Phi / m_B()) / (1.0 - sigma);

            // three-particle contributions
            const auto iota1 = a_0_iota1(q2, sigma);
            const auto iota2 = a_0_iota2(q2, sigma);
            const auto iota3 = a_0_iota3(q2, sigma);

            const auto result_3p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (-iota1 + iota2 / M2 - iota3 / (2.0 * M2 * M2));

            return result_2p + result_3p;
        }

        inline double a_0(const double & q2) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = std::pow(m_Kstar(), 2);
            const auto m_b      = model->m_b_msbar(mu());

            const std::function<double (const double &)> integrand = std::bind(&Implementation<AnalyticFormFactorBToKstarKMO2006>::a_0_integrand, this, q2, std::placeholders::_1);

            const auto sigma0  = this->sigma0(q2);
            const auto iota20  = a_0_iota2(q2, sigma0);
            const auto iota30  = a_0_iota3(q2, sigma0);
            const auto Diota30 = a_0_Diota3(q2, sigma0);
            const auto etaf0   = etaf(q2, sigma0);
            const auto Detaf0  = Detaf(q2, sigma0);

            const auto integral = integrate<GSL::QNG>(integrand, 0.0, sigma0);
            const auto delta = (iota20 - (1.0 / M2 + Detaf0 / m_B2) / 2.0 * iota30 - etaf0 / (2.0 * m_B2) * Diota30)
                * std::exp(-s0 / M2) / m_B2 * etaf0;

            return f_B * m_B2 * m_b / (2.0 * f_Kstar * m_Kstar2) * std::exp(m_Kstar2 / M2)
                * (integral + delta);
        }

        /* A_1 */

        inline double a_1_iota1(const double & /*q2*/, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-7)
                return 0.0; // error ~ -9e-9

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return exp(-xi)*(((24*gamma*lambdaB2 - 24*lambdaB*mB*sigma + 8*gamma*lambdaB*mB*sigma + 2*mB2*sigma - 3*mB2*sigma2 -
                            24*lambdaB2*log(lambdaB) - 8*lambdaB*mB*sigma*log(lambdaB) + 24*lambdaB2*log(mB*sigma) +
                            8*lambdaB*mB*sigma*log(mB*sigma))* pow(mB,-3)*pow(-1 + sigma,-3))/4. - exp(xi)*(6*lambdaB2*Ei - 4*lambdaB*mB*sigma*Ei +
                        mB2*sigma2*Ei + 3*lambdaB2*log(xi) - 2*lambdaB*mB*sigma*log(xi) + 3*lambdaB2*log(lambdaB*pow(mB,-1)* pow(sigma,-1)) -
                        2*lambdaB*mB*sigma* log(lambdaB*pow(mB,-1)*pow(sigma,-1)))* pow(mB,-3)*pow(-1 + sigma,-3));
        }

        inline double a_1_iota2(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return exp(-xi)*(-(exp(xi)*pow(mB,-3)*pow(-1 + sigma,-4)* (24*lambdaB2*q2*Ei - 16*lambdaB*mB*q2*sigma* Ei + 4*mB2*q2*sigma2*Ei +
                            24*lambdaB2*q2*log(lambdaB) - 16*lambdaB*mB*q2*sigma*log(lambdaB) - 20*lambdaB2*mB2*log(mB*sigma) -
                            24*lambdaB2*q2*log(mB*sigma) + 40*lambdaB2*mB2*sigma*log(mB*sigma) + 16*lambdaB*mB*q2*sigma*log(mB*sigma) -
                            20*lambdaB2*mB2*sigma2*log(mB*sigma) + 20*lambdaB2*mB2*log(xi) + 24*lambdaB2*q2*log(xi) - 40*lambdaB2*mB2*sigma*log(xi) -
                            16*lambdaB*mB*q2*sigma*log(xi) + 20*lambdaB2*mB2*sigma2*log(xi) + 11*lambdaB*sigma*log(mB*sigma)* pow(mB,3) -
                            22*lambdaB*sigma2*log(mB*sigma)* pow(mB,3) - 11*lambdaB*sigma*log(xi)*pow(mB,3) + 22*lambdaB*sigma2*log(xi)*pow(mB,3) +
                            20*lambdaB2*mB2*Ei* pow(-1 + sigma,2) - 11*lambdaB*mB*mB2*sigma* Ei*pow(-1 + sigma,2) + 20*lambdaB2*mB2*log(lambdaB)*
                            pow(-1 + sigma,2) - 11*lambdaB*mB*mB2*sigma*log(lambdaB)* pow(-1 + sigma,2) + sigma2*Ei*pow(mB2,2)* pow(-1 + sigma,2) +
                            11*lambdaB*log(mB*sigma)*pow(mB,3)* pow(sigma,3) - 11*lambdaB*log(xi)*pow(mB,3)*pow(sigma,3)))/2. - (pow(mB,-3)*pow(-1 +
                            sigma,-4)* (-20*gamma*lambdaB2*mB2 - 24*gamma*lambdaB2*q2 + 40*gamma*lambdaB2*mB2*sigma + 24*lambdaB*mB*q2*sigma -
                                8*gamma*lambdaB*mB*q2*sigma - mB2*q2*sigma - 20*gamma*lambdaB2*mB2*sigma2 + 2*mB2*q2*sigma2 +
                                24*lambdaB2*q2*log(lambdaB) + 8*lambdaB*mB*q2*sigma*log(lambdaB) - 20*lambdaB2*mB2*log(mB*sigma) -
                                24*lambdaB2*q2*log(mB*sigma) + 40*lambdaB2*mB2*sigma*log(mB*sigma) - 8*lambdaB*mB*q2*sigma*log(mB*sigma) -
                                20*lambdaB2*mB2*sigma2*log(mB*sigma) + 20*lambdaB*sigma*pow(mB,3) - 9*gamma*lambdaB*sigma*pow(mB,3) -
                                40*lambdaB*sigma2*pow(mB,3) + 18*gamma*lambdaB*sigma2*pow(mB,3) - 9*lambdaB*sigma*log(mB*sigma)*pow(mB,3) +
                                18*lambdaB*sigma2*log(mB*sigma)* pow(mB,3) + sigma*pow(mB,4) - sigma2*pow(mB,4) + 20*lambdaB2*mB2*log(lambdaB)* pow(-1
                                    + sigma,2) + 9*lambdaB*mB*mB2*sigma*log(lambdaB)* pow(-1 + sigma,2) + 20*lambdaB*pow(mB,3)*pow(sigma,3) -
                                9*gamma*lambdaB*pow(mB,3)* pow(sigma,3) - 9*lambdaB*log(mB*sigma)*pow(mB,3)* pow(sigma,3) - pow(mB,4)*pow(sigma,3) +
                                pow(mB,4)*pow(sigma,4)))/2.);
        }

        inline double a_1_iota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-4)
                return 0.0; // error ~ -1e-10

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return exp(-xi)*(((24*gamma*lambdaB2 - 24*lambdaB*mB*sigma + 8*gamma*lambdaB*mB*sigma - mB2*sigma2 - 24*lambdaB2*log(lambdaB) -
                            8*lambdaB*mB*sigma*log(lambdaB) + 24*lambdaB2*log(mB*sigma) + 8*lambdaB*mB*sigma*log(mB*sigma))* pow(mB,-3)*pow(-1 +
                            sigma,-5)* pow(q2 - mB2*pow(-1 + sigma,2),2))/4. - exp(xi)*(6*lambdaB2*Ei - 4*lambdaB*mB*sigma*Ei + mB2*sigma2*Ei +
                            6*lambdaB2*log(lambdaB) - 4*lambdaB*mB*sigma*log(lambdaB) - 6*lambdaB2*log(mB*sigma) + 4*lambdaB*mB*sigma*log(mB*sigma) +
                            6*lambdaB2*log(xi) - 4*lambdaB*mB*sigma*log(xi))*pow(mB,-3)* pow(-1 + sigma,-5)* pow(q2 - mB2*pow(-1 + sigma,2),2));
        }

        inline double a_1_Diota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -4e-10

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return exp(-xi)*(exp(xi)*pow(lambdaB,-1)*pow(mB,-3)* pow(-1 + sigma,-6)* (q2 - mB2*pow(-1 + sigma,2))* (30*lambdaB*lambdaB2*q2*Ei -
                        4*lambdaB2*mB*q2*Ei - 16*lambdaB2*mB*q2*sigma* Ei + 2*lambdaB*mB2*q2*sigma*Ei + 3*lambdaB*mB2*q2*sigma2* Ei +
                        30*lambdaB*lambdaB2*q2*log(lambdaB) - 4*lambdaB2*mB*q2*log(lambdaB) - 16*lambdaB2*mB*q2*sigma*log(lambdaB) +
                        4*lambdaB2*mB*q2*log(mB*sigma) + 16*lambdaB2*mB*q2*sigma*log(mB*sigma) - 4*lambdaB2*mB*q2*log(xi) -
                        16*lambdaB2*mB*q2*sigma*log(xi) + 6*mB2*log(mB*sigma)*pow(lambdaB,3) - 30*q2*log(mB*sigma)*pow(lambdaB,3) -
                        12*mB2*sigma*log(mB*sigma)*pow(lambdaB,3) + 6*mB2*sigma2*log(mB*sigma)*pow(lambdaB,3) - 6*mB2*log(xi)*pow(lambdaB,3) +
                        30*q2*log(xi)*pow(lambdaB,3) + 12*mB2*sigma*log(xi)*pow(lambdaB,3) - 6*mB2*sigma2*log(xi)*pow(lambdaB,3) -
                        4*lambdaB2*log(mB*sigma)*pow(mB,3) + 8*lambdaB2*sigma*log(mB*sigma)*pow(mB,3) - 4*lambdaB2*sigma2*log(mB*sigma)*pow(mB,3) +
                        4*lambdaB2*log(xi)*pow(mB,3) - 8*lambdaB2*sigma*log(xi)*pow(mB,3) + 4*lambdaB2*sigma2*log(xi)*pow(mB,3) -
                        6*lambdaB*lambdaB2*mB2*Ei* pow(-1 + sigma,2) + 4*lambdaB2*mB*mB2*Ei* pow(-1 + sigma,2) - 6*lambdaB*lambdaB2*mB2*log(lambdaB)*
                        pow(-1 + sigma,2) + 4*lambdaB2*mB*mB2*log(lambdaB)* pow(-1 + sigma,2) - 2*lambdaB*sigma*Ei* pow(mB2,2)*pow(-1 + sigma,2) +
                        lambdaB*sigma2*Ei* pow(mB2,2)*pow(-1 + sigma,2)) + (pow(lambdaB,-1)*pow(mB,-3)*pow(-1 + sigma,-6)* (q2 - mB2*pow(-1 +
                                sigma,2))* (16*gamma*lambdaB2*mB*q2 + 120*lambdaB2*mB*q2*sigma - 56*gamma*lambdaB2*mB*q2*sigma -
                                18*lambdaB*mB2*q2*sigma + 8*gamma*lambdaB*mB2*q2*sigma + 23*lambdaB*mB2*q2*sigma2 - 8*gamma*lambdaB*mB2*q2*sigma2 +
                                120*lambdaB*lambdaB2*q2*log(lambdaB) - 16*lambdaB2*mB*q2*log(lambdaB) + 56*lambdaB2*mB*q2*sigma*log(lambdaB) -
                                8*lambdaB*mB2*q2*sigma*log(lambdaB) + 8*lambdaB*mB2*q2*sigma2*log(lambdaB) + 16*lambdaB2*mB*q2*log(mB*sigma) -
                                56*lambdaB2*mB*q2*sigma*log(mB*sigma) + 8*lambdaB*mB2*q2*sigma*log(mB*sigma) - 8*lambdaB*mB2*q2*sigma2*log(mB*sigma) +
                                24*gamma*mB2*pow(lambdaB,3) - 120*gamma*q2*pow(lambdaB,3) - 48*gamma*mB2*sigma*pow(lambdaB,3) +
                                24*gamma*mB2*sigma2*pow(lambdaB,3) + 24*mB2*log(mB*sigma)*pow(lambdaB,3) - 120*q2*log(mB*sigma)*pow(lambdaB,3) -
                                48*mB2*sigma*log(mB*sigma)* pow(lambdaB,3) + 24*mB2*sigma2*log(mB*sigma)* pow(lambdaB,3) - 16*gamma*lambdaB2*pow(mB,3)
                                - 24*lambdaB2*sigma*pow(mB,3) + 56*gamma*lambdaB2*sigma*pow(mB,3) + 48*lambdaB2*sigma2*pow(mB,3) -
                                64*gamma*lambdaB2*sigma2*pow(mB,3) - q2*sigma2*pow(mB,3) - 16*lambdaB2*log(mB*sigma)*pow(mB,3) +
                                56*lambdaB2*sigma*log(mB*sigma)* pow(mB,3) - 64*lambdaB2*sigma2*log(mB*sigma)* pow(mB,3) + 18*lambdaB*sigma*pow(mB,4)
                                - 8*gamma*lambdaB*sigma*pow(mB,4) - 55*lambdaB*sigma2*pow(mB,4) + 24*gamma*lambdaB*sigma2*pow(mB,4) -
                                8*lambdaB*sigma*log(mB*sigma)*pow(mB,4) + 24*lambdaB*sigma2*log(mB*sigma)* pow(mB,4) + sigma2*pow(mB,5) -
                                24*lambdaB*lambdaB2*mB2*log(lambdaB)* pow(-1 + sigma,2) + 16*lambdaB2*mB*mB2*log(lambdaB)* pow(-1 + sigma,2) -
                                24*lambdaB2*mB*mB2*sigma*log(lambdaB)* pow(-1 + sigma,2) + 8*lambdaB*sigma*log(lambdaB)*pow(mB2,2)* pow(-1 + sigma,2)
                                - 8*lambdaB*sigma2*log(lambdaB)*pow(mB2,2)* pow(-1 + sigma,2) - 24*lambdaB2*pow(mB,3)*pow(sigma,3) +
                                24*gamma*lambdaB2*pow(mB,3)* pow(sigma,3) + q2*pow(mB,3)*pow(sigma,3) + 24*lambdaB2*log(mB*sigma)*pow(mB,3)*
                                pow(sigma,3) + 56*lambdaB*pow(mB,4)*pow(sigma,3) - 24*gamma*lambdaB*pow(mB,4)* pow(sigma,3) -
                                24*lambdaB*log(mB*sigma)*pow(mB,4)* pow(sigma,3) - 3*pow(mB,5)*pow(sigma,3) - 19*lambdaB*pow(mB,4)*pow(sigma,4) +
                                8*gamma*lambdaB*pow(mB,4)* pow(sigma,4) + 8*lambdaB*log(mB*sigma)*pow(mB,4)* pow(sigma,4) + 3*pow(mB,5)*pow(sigma,4) -
                                pow(mB,5)*pow(sigma,5)))/4.);
        }

        double a_1_integrand(const double & q2, const double & sigma) const
        {
            const auto m_B2 = std::pow(m_B(), 2);

            // two-particle contribution
            const auto result_2p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * b_lcdas.phi_plus(sigma * m_B()) * (m_B2 - q2 / pow(1.0 - sigma, 2));

            // three-particle contributions
            const auto iota1 = a_1_iota1(q2, sigma);
            const auto iota2 = a_1_iota2(q2, sigma);
            const auto iota3 = a_1_iota3(q2, sigma);

            const auto result_3p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (-iota1 + iota2 / M2 - iota3 / (2.0 * M2 * M2));

            return result_2p + m_B2 * result_3p;
        }

        inline double a_1(const double & q2) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = std::pow(m_Kstar(), 2);

            const std::function<double (const double &)> integrand = std::bind(&Implementation<AnalyticFormFactorBToKstarKMO2006>::a_1_integrand, this, q2, std::placeholders::_1);

            const auto sigma0  = this->sigma0(q2);
            const auto iota20  = a_1_iota2(q2, sigma0);
            const auto iota30  = a_1_iota3(q2, sigma0);
            const auto Diota30 = a_1_Diota3(q2, sigma0);
            const auto etaf0   = etaf(q2, sigma0);
            const auto Detaf0  = Detaf(q2, sigma0);

            const auto integral = integrate<GSL::QNG>(integrand, 0.0, sigma0);
            const auto delta = (iota20 - (1.0 / M2 + Detaf0 / m_B2) / 2.0 * iota30 - etaf0 / (2.0 * m_B2) * Diota30)
                * std::exp(-s0 / M2) / m_B2 * etaf0;

            return f_B * m_B / (m_B + m_Kstar) / (2.0 * f_Kstar * m_Kstar) * std::exp(m_Kstar2 / M2)
                * (integral + m_B2 * delta);
        }

        /* A_2 */

        inline double a_2_iota1(const double & /*q2*/, const double & /*sigma*/) const
        {
            return 0.0;
        }

        inline double a_2_iota2(const double & /*q2*/, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return exp(-xi)*(-(sigma*(24*gamma*lambdaB2 - 2*mB2 - 24*lambdaB*mB*sigma + 8*gamma*lambdaB*mB*sigma - 2*mB2*sigma + 3*mB2*sigma2 -
                            24*lambdaB2*log(lambdaB) - 8*lambdaB*mB*sigma*log(lambdaB) + 24*lambdaB2*log(mB*sigma) +
                            8*lambdaB*mB*sigma*log(mB*sigma))* pow(mB,-1)*pow(-1 + sigma,-3))/4. - sigma*exp(xi)*(-6*lambdaB2* Ei +
                        4*lambdaB*mB*sigma*Ei - mB2*sigma2*Ei - 6*lambdaB2*log(lambdaB) + 4*lambdaB*mB*sigma*log(lambdaB) + 6*lambdaB2*log(mB*sigma) -
                        4*lambdaB*mB*sigma*log(mB*sigma) - 6*lambdaB2*log(xi) + 4*lambdaB*mB*sigma*log(xi))*pow(mB,-1)* pow(-1 + sigma,-3));
        }

        inline double a_2_iota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-4)
                return 0.0; // error ~ -1e-10

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return exp(-xi)*(exp(xi)*Ei*pow(mB,-1)* pow(-1 + sigma,-4)* (-6*lambdaB2*q2 + 12*lambdaB2*q2*sigma + 4*lambdaB*mB*q2*sigma -
                        8*lambdaB*mB*q2*sigma2 - mB2*q2*sigma2 + 2*mB2*q2*sigma*sigma2 + 6*lambdaB2*mB2*pow(-1 + sigma,2) +
                        20*lambdaB2*mB2*sigma*pow(-1 + sigma,2) - 4*lambdaB*mB*mB2*sigma*pow(-1 + sigma,2) - 11*lambdaB*mB*mB2*sigma2* pow(-1 +
                            sigma,2) + sigma2*pow(mB2,2)*pow(-1 + sigma,2) + sigma*sigma2*pow(mB2,2)*pow(-1 + sigma,2)) - (pow(mB,-1)*pow(-1 +
                                sigma,-4)* (24*gamma*lambdaB2*mB2 - 24*gamma*lambdaB2*q2 + 32*gamma*lambdaB2*mB2*sigma + 48*gamma*lambdaB2*q2*sigma +
                                    24*lambdaB*mB*q2*sigma - 8*gamma*lambdaB*mB*q2*sigma - 136*gamma*lambdaB2*mB2*sigma2 - 48*lambdaB*mB*q2*sigma2 +
                                    16*gamma*lambdaB*mB*q2*sigma2 + mB2*q2*sigma2 + 24*lambdaB2*q2*log(lambdaB) - 48*lambdaB2*q2*sigma*log(lambdaB) +
                                    8*lambdaB*mB*q2*sigma*log(lambdaB) - 16*lambdaB*mB*q2*sigma2*log(lambdaB) + 24*lambdaB2*mB2*log(mB*sigma) -
                                    24*lambdaB2*q2*log(mB*sigma) + 32*lambdaB2*mB2*sigma*log(mB*sigma) + 48*lambdaB2*q2*sigma*log(mB*sigma) -
                                    8*lambdaB*mB*q2*sigma*log(mB*sigma) - 136*lambdaB2*mB2*sigma2*log(mB*sigma) +
                                    16*lambdaB*mB*q2*sigma2*log(mB*sigma) - 24*lambdaB*sigma*pow(mB,3) + 8*gamma*lambdaB*sigma*pow(mB,3) -
                                    32*lambdaB*sigma2*pow(mB,3) + 20*gamma*lambdaB*sigma2*pow(mB,3) + 8*lambdaB*sigma*log(mB*sigma)*pow(mB,3) +
                                    20*lambdaB*sigma2*log(mB*sigma)* pow(mB,3) - sigma2*pow(mB,4) - 24*lambdaB2*mB2*log(lambdaB)* pow(-1 + sigma,2) -
                                    80*lambdaB2*mB2*sigma*log(lambdaB)* pow(-1 + sigma,2) - 8*lambdaB*sigma*log(lambdaB)*pow(mB,3)* pow(-1 + sigma,2)
                                    - 36*lambdaB*sigma2*log(lambdaB)*pow(mB,3)* pow(-1 + sigma,2) + 80*gamma*lambdaB2*mB2*pow(sigma,3) -
                                    2*mB2*q2*pow(sigma,3) + 80*lambdaB2*mB2*log(mB*sigma)* pow(sigma,3) + 136*lambdaB*pow(mB,3)*pow(sigma,3) -
                                    64*gamma*lambdaB*pow(mB,3)* pow(sigma,3) - 64*lambdaB*log(mB*sigma)*pow(mB,3)* pow(sigma,3) -
                                    6*pow(mB,4)*pow(sigma,3) - 80*lambdaB*pow(mB,3)*pow(sigma,4) + 36*gamma*lambdaB*pow(mB,3)* pow(sigma,4) +
                                    36*lambdaB*log(mB*sigma)*pow(mB,3)* pow(sigma,4) + 15*pow(mB,4)*pow(sigma,4) - 8*pow(mB,4)*pow(sigma,5)))/4.);
        }

        inline double a_2_Diota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -4e-10

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return exp(-xi)*(lambdaB*exp(xi)*Ei* pow(lambdaB,-1)*pow(mB,-1)*pow(-1 + sigma,-5)* (12*lambdaB2*q2 - 4*lambdaB*mB*q2 -
                        36*lambdaB2*q2*sigma + 4*lambdaB*mB*q2*sigma + 2*mB2*q2*sigma + 16*lambdaB*mB*q2*sigma2 - 4*mB2*q2*sigma2 -
                        2*mB2*q2*sigma*sigma2 - 32*lambdaB2*mB2*pow(-1 + sigma,2) + 4*lambdaB*mB*mB2*pow(-1 + sigma,2) - 20*lambdaB2*mB2*sigma*pow(-1
                            + sigma,2) + 26*lambdaB*mB*mB2*sigma*pow(-1 + sigma,2) - 2*sigma*pow(mB2,2)*pow(-1 + sigma,2) - 3*sigma2*pow(mB2,2)*pow(-1
                                + sigma,2) + sigma*sigma2*pow(mB2,2)*pow(-1 + sigma,2)) - (pow(lambdaB,-1)*pow(mB,-1)* pow(-1 + sigma,-5)*
                                (-16*gamma*lambdaB2*mB*q2 - 48*lambdaB2*mB*q2*sigma + 64*gamma*lambdaB2*mB*q2*sigma + 18*lambdaB*mB2*q2*sigma -
                                 8*gamma*lambdaB*mB2*q2*sigma + 144*lambdaB2*mB*q2*sigma2 - 80*gamma*lambdaB2*mB*q2*sigma2 - 56*lambdaB*mB2*q2*sigma2
                                 + 24*gamma*lambdaB*mB2*q2*sigma2 - 48*lambdaB*lambdaB2*q2*log(lambdaB) + 16*lambdaB2*mB*q2*log(lambdaB) +
                                 144*lambdaB*lambdaB2*q2*sigma* log(lambdaB) - 64*lambdaB2*mB*q2*sigma*log(lambdaB) +
                                 8*lambdaB*mB2*q2*sigma*log(lambdaB) + 80*lambdaB2*mB*q2*sigma2*log(lambdaB) - 24*lambdaB*mB2*q2*sigma2*log(lambdaB) -
                                 16*lambdaB2*mB*q2*log(mB*sigma) + 64*lambdaB2*mB*q2*sigma*log(mB*sigma) - 8*lambdaB*mB2*q2*sigma*log(mB*sigma) -
                                 80*lambdaB2*mB*q2*sigma2*log(mB*sigma) + 24*lambdaB*mB2*q2*sigma2*log(mB*sigma) - 128*gamma*mB2*pow(lambdaB,3) +
                                 48*gamma*q2*pow(lambdaB,3) + 176*gamma*mB2*sigma*pow(lambdaB,3) - 144*gamma*q2*sigma*pow(lambdaB,3) +
                                 32*gamma*mB2*sigma2*pow(lambdaB,3) - 128*mB2*log(mB*sigma)*pow(lambdaB,3) + 48*q2*log(mB*sigma)*pow(lambdaB,3) +
                                 176*mB2*sigma*log(mB*sigma)* pow(lambdaB,3) - 144*q2*sigma*log(mB*sigma)* pow(lambdaB,3) +
                                 32*mB2*sigma2*log(mB*sigma)* pow(lambdaB,3) + 16*gamma*lambdaB2*pow(mB,3) + 128*lambdaB2*sigma*pow(mB,3) -
                                 56*gamma*lambdaB2*sigma*pow(mB,3) - 176*lambdaB2*sigma2*pow(mB,3) - 16*gamma*lambdaB2*sigma2*pow(mB,3) +
                                 q2*sigma2*pow(mB,3) + 16*lambdaB2*log(mB*sigma)*pow(mB,3) - 56*lambdaB2*sigma*log(mB*sigma)* pow(mB,3) -
                                 16*lambdaB2*sigma2*log(mB*sigma)* pow(mB,3) - 18*lambdaB*sigma*pow(mB,4) + 8*gamma*lambdaB*sigma*pow(mB,4) +
                                 4*lambdaB*sigma2*pow(mB,4) + 12*gamma*lambdaB*sigma2*pow(mB,4) + 8*lambdaB*sigma*log(mB*sigma)*pow(mB,4) +
                                 12*lambdaB*sigma2*log(mB*sigma)* pow(mB,4) - sigma2*pow(mB,5) + 128*lambdaB*lambdaB2*mB2*log(lambdaB)* pow(-1 +
                                     sigma,2) - 16*lambdaB2*mB*mB2*log(lambdaB)* pow(-1 + sigma,2) + 80*lambdaB*lambdaB2*mB2*sigma*
                                 log(lambdaB)*pow(-1 + sigma,2) + 24*lambdaB2*mB*mB2*sigma*log(lambdaB)* pow(-1 + sigma,2) +
                                 80*lambdaB2*mB*mB2*sigma2*log(lambdaB)* pow(-1 + sigma,2) - 8*lambdaB*sigma*log(lambdaB)*pow(mB2,2)* pow(-1 +
                                     sigma,2) - 28*lambdaB*sigma2*log(lambdaB)*pow(mB2,2)* pow(-1 + sigma,2) + 42*lambdaB*mB2*q2*pow(sigma,3) -
                                 16*gamma*lambdaB*mB2*q2* pow(sigma,3) + 16*lambdaB*mB2*q2*log(lambdaB)* pow(sigma,3) -
                                 16*lambdaB*mB2*q2*log(mB*sigma)* pow(sigma,3) - 80*gamma*mB2*pow(lambdaB,3)* pow(sigma,3) -
                                 80*mB2*log(mB*sigma)*pow(lambdaB,3)* pow(sigma,3) - 32*lambdaB2*pow(mB,3)*pow(sigma,3) +
                                 136*gamma*lambdaB2*pow(mB,3)* pow(sigma,3) - 3*q2*pow(mB,3)*pow(sigma,3) + 136*lambdaB2*log(mB*sigma)*pow(mB,3)*
                                 pow(sigma,3) + 114*lambdaB*pow(mB,4)*pow(sigma,3) - 84*gamma*lambdaB*pow(mB,4)* pow(sigma,3) -
                                 84*lambdaB*log(mB*sigma)*pow(mB,4)* pow(sigma,3) - 5*pow(mB,5)*pow(sigma,3) + 36*lambdaB*log(lambdaB)*pow(mB2,2)*
                                 pow(-1 + sigma,2)*pow(sigma,3) + 80*lambdaB2*pow(mB,3)*pow(sigma,4) - 80*gamma*lambdaB2*pow(mB,3)* pow(sigma,4) +
                                 2*q2*pow(mB,3)*pow(sigma,4) - 80*lambdaB2*log(mB*sigma)*pow(mB,3)* pow(sigma,4) - 168*lambdaB*pow(mB,4)*pow(sigma,4)
                                 + 100*gamma*lambdaB*pow(mB,4)* pow(sigma,4) + 100*lambdaB*log(mB*sigma)*pow(mB,4)* pow(sigma,4) +
                                 21*pow(mB,5)*pow(sigma,4) + 68*lambdaB*pow(mB,4)*pow(sigma,5) - 36*gamma*lambdaB*pow(mB,4)* pow(sigma,5) -
                                 36*lambdaB*log(mB*sigma)*pow(mB,4)* pow(sigma,5) - 23*pow(mB,5)*pow(sigma,5) + 8*pow(mB,5)*pow(sigma,6)))/4.);
        }

        double a_2_integrand(const double & q2, const double & sigma) const
        {
            const auto m_B2 = std::pow(m_B(), 2);
            const auto m_B3 = std::pow(m_B(), 3);

            // two-particle contributions
            const auto c_p     = 1.0 - sigma / (1.0 - sigma);
            const auto c_delta = 2.0 * sigma * (1.0 - sigma) * m_B2 / (pow(1.0 - sigma, 2) * m_B2 - q2);
            const auto c_bar   = 4.0 * sigma * pow(1.0 - sigma, 2) * m_B3
                               / pow(pow(1.0 - sigma, 2) * m_B2 - q2, 2)
                               + 2.0 * (1.0 - 2.0 * sigma) * m_B / (pow(1.0 - sigma, 2) * m_B2 - q2);

            const auto phi_p     = b_lcdas.phi_plus(sigma * m_B());
            const auto phi_m     = b_lcdas.phi_minus(sigma * m_B());
            const auto Phi_bar   = b_lcdas.Phibar(sigma * m_B());
            const auto phi_delta = phi_p - phi_m;

            const auto result_2p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (c_p * phi_p + c_delta * phi_delta + c_bar * Phi_bar);

            // three-particle contributions
            const auto iota1 = a_2_iota1(q2, sigma);
            const auto iota2 = a_2_iota2(q2, sigma);
            const auto iota3 = a_2_iota3(q2, sigma);

            const auto result_3p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (-iota1 + iota2 / M2 - iota3 / (2.0 * M2 * M2));

            return result_2p + result_3p;
        }

        inline double a_2(const double & q2) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = std::pow(m_Kstar(), 2);

            const std::function<double (const double &)> integrand = std::bind(&Implementation<AnalyticFormFactorBToKstarKMO2006>::a_2_integrand, this, q2, std::placeholders::_1);

            const auto sigma0  = this->sigma0(q2);
            const auto iota20  = a_2_iota2(q2, sigma0);
            const auto iota30  = a_2_iota3(q2, sigma0);
            const auto Diota30 = a_2_Diota3(q2, sigma0);
            const auto etaf0   = etaf(q2, sigma0);
            const auto Detaf0  = Detaf(q2, sigma0);

            const auto integral = integrate<GSL::QNG>(integrand, 0.0, sigma0);
            const auto delta = (iota20 - (1.0 / M2 + Detaf0 / m_B2) / 2.0 * iota30 - etaf0 / (2.0 * m_B2) * Diota30)
                * std::exp(-s0 / M2) / m_B2 * etaf0;

            return f_B * m_B * (m_B + m_Kstar) / (2.0 * f_Kstar * m_Kstar) * std::exp(m_Kstar2 / M2)
                * (integral + delta);
        }

        /* A_12 */

        inline double a_12(const double & q2) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = std::pow(m_Kstar(), 2);
            const auto lambda   = eos::lambda(m_B2, m_Kstar2, q2);

            return (pow(m_B + m_Kstar, 2) * (m_B2 - m_Kstar2 - q2) * this->a_1(q2)
                - lambda * this->a_2(q2)) / (16.0 * m_B * m_Kstar2 * (m_B + m_Kstar));
        }

        /* T_1 */

        inline double t_1_iota1(const double & /*q2*/, const double & /*sigma*/) const
        {
            return 0.0;
        }

        inline double t_1_iota2(const double & /*q2*/, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p();
            const auto mB = m_B();
            const auto xi = mB * sigma / lambdaB, xi2 = pow(xi, 2);
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return (lambdaB*exp(-xi)*(20*gamma*lambdaB - 20*lambdaB*xi +
            9*gamma*lambdaB*xi - 2*mB*xi - Ei*lambdaB*(20 - 11*xi +
            xi2)*exp(xi) + lambdaB*(20 + 9*xi)*log(xi)))/ (4.*mB*pow(-1 +
            sigma,2));
        }

        inline double t_1_iota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto xi = mB * sigma / lambdaB, xi2 = pow(xi, 2);
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return (lambdaB2*(-q2 + mB2*pow(-1 + sigma,2))*exp(-xi)*
            (24*gamma - 24*xi + 8*gamma*xi - xi2 - 4*Ei*(6 - 4*xi +
            xi2)*exp(xi) + 8*(3 + xi)*log(xi)))/(4.*mB*pow(-1 + sigma,3));
        }

        inline double t_1_Diota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(),
                lambdaB2 = pow(lambdaB, 2),
                lambdaB3 = pow(lambdaB, 3);
            const auto mB = m_B(),
                mB2 = pow(mB, 2),
                mB3 = pow(mB, 3),
                mB4 = pow(mB, 4);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return (exp(-xi)*((mB4*sigma)/lambdaB + (mB3* (19*lambdaB -
            8*gamma*lambdaB + mB - 4*Ei*lambdaB*exp(xi) - 8*lambdaB*log(xi)))
            /lambdaB + (mB2*(24*lambdaB2 - 24*gamma*lambdaB2 + 4*(5 -
            2*gamma)*lambdaB*mB + mB2 - q2 - 8*lambdaB*(3*lambdaB +
            mB)*log(xi)))/(lambdaB*(-1 + sigma)) - (3*q2*(-24*gamma*lambdaB2 +
            24*lambdaB*mB - 8*gamma*lambdaB*mB + mB2 + 4*Ei*(6*lambdaB2 -
            4*lambdaB*mB + mB2)*exp(xi) - 8*lambdaB*(3*lambdaB +
            mB)*log(xi)))/(mB*pow(-1 + sigma,4)) + (q2*(-72*lambdaB2 +
            40*gamma*lambdaB2 - 24*lambdaB*mB + 8*gamma*lambdaB*mB - mB2 -
            16*Ei*lambdaB*(-2*lambdaB + mB)*exp(xi) + 8*lambdaB*(5*lambdaB +
            mB)*log(xi)))/(lambdaB*pow(-1 + sigma,3)) +
            (mB*(-24*gamma*lambdaB3 + 24*lambdaB2*mB - 8*gamma*lambdaB2*mB +
            lambdaB*mB2 - 21*lambdaB*q2 + 8*gamma*lambdaB*q2 - 2*mB*q2 +
            4*Ei*lambdaB*(6*lambdaB2 - 4*lambdaB*mB + mB2 - q2)*exp(xi) +
            8*lambdaB*(-(lambdaB*(3*lambdaB + mB)) + q2)*log(xi)))/
            (lambdaB*pow(-1 + sigma,2))))/4.;
        }

        double t_1_integrand(const double & q2, const double & sigma) const
        {
            const auto m_B2 = std::pow(m_B(), 2);

            // two-particle contribution
            const auto result_2p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * b_lcdas.phi_plus(sigma * m_B());

            // three-particle contributions
            const auto iota1 = t_1_iota1(q2, sigma);
            const auto iota2 = t_1_iota2(q2, sigma);
            const auto iota3 = t_1_iota3(q2, sigma);

            const auto result_3p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (-iota1 + iota2 / M2 - iota3 / (2.0 * M2 * M2));

            return result_2p + result_3p;
        }

        inline double t_1(const double & q2) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = std::pow(m_Kstar(), 2);

            const std::function<double (const double &)> integrand = std::bind(&Implementation<AnalyticFormFactorBToKstarKMO2006>::t_1_integrand, this, q2, std::placeholders::_1);

            const auto sigma0  = this->sigma0(q2);
            const auto iota20  = t_1_iota2(q2, sigma0);
            const auto iota30  = t_1_iota3(q2, sigma0);
            const auto Diota30 = t_1_Diota3(q2, sigma0);
            const auto etaf0   = etaf(q2, sigma0);
            const auto Detaf0  = Detaf(q2, sigma0);

            const auto integral = integrate<GSL::QNG>(integrand, 0.0, sigma0);
            const auto delta = (iota20 - (1.0 / M2 + Detaf0 / m_B2) / 2.0 * iota30 - etaf0 / (2.0 * m_B2) * Diota30)
                * std::exp(-s0 / M2) / m_B2 * etaf0;

            return f_B * m_B2 / (2.0 * f_Kstar * m_Kstar) * std::exp(m_Kstar2 / M2)
                * (integral + delta);
        }

        /* T_23A */

        inline double t_23a_iota1(const double & /*q2*/, const double & /*sigma*/) const
        {
            return 0.0;
        }

        inline double t_23a_iota2(const double & /*q2*/, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto xi = mB * sigma / lambdaB, xi2 = pow(xi, 2);
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return (exp(-xi)*(20*gamma*lambdaB2 - 2*mB2*sigma - 20*lambdaB2*xi
            + 9*gamma*lambdaB2*xi - Ei*lambdaB2*(20 - 11*xi + xi2)* exp(xi) +
            lambdaB2*(20 + 9*xi)*log(xi)))/ (4.*mB*pow(-1 + sigma,2));
        }

        inline double t_23a_iota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2), lambdaB3 = pow(lambdaB, 3), lambdaB4 = pow(lambdaB, 4);
            const auto mB = m_B(), mB2 = pow(mB, 2), mB3 = pow(mB, 3), mB4 = pow(mB, 4);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB, xi2 = pow(xi, 2);
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return (exp(-xi)*(24*gamma*lambdaB2*mB2 - 24*gamma*lambdaB2*q2 -
            48*gamma*lambdaB2*mB2*sigma - 24*lambdaB*mB3*sigma +
            8*gamma*lambdaB*mB3*sigma + 128*gamma*lambdaB2*q2*sigma +
            2*mB4*pow(sigma,3) + 48*lambdaB*mB3*sigma2 -
            16*gamma*lambdaB*mB3*sigma2 - mB4*sigma2 - 128*lambdaB*mB*q2*sigma2
            + 52*gamma*lambdaB*mB*q2*sigma2 + 24*lambdaB2*q2*xi -
            8*gamma*lambdaB2*q2*xi - 24*lambdaB4*pow(xi,3) +
            8*gamma*lambdaB4*pow(xi,3) - (10*lambdaB3*q2*pow(xi,3))/ mB -
            lambdaB4*pow(xi,4) + 24*gamma*lambdaB4*xi2 + lambdaB2*q2*xi2 -
            4*Ei*(-4*lambdaB*mB3* pow(-1 + sigma,2)*sigma + mB4*pow(-1 +
            sigma,2)* sigma2 + mB2*(6*lambdaB2* pow(-1 + sigma,2) + q2*(-1 +
            3*sigma)*sigma2) + lambdaB2*q2* (-6 + 32*sigma + (4 -
            19*sigma)*xi))* exp(xi) + 4*lambdaB* (6*lambdaB*mB2* pow(-1 +
            sigma,2) + 2*mB3*pow(-1 + sigma,2)* sigma + lambdaB*q2* (-6 +
            32*sigma + (-2 + 13*sigma)*xi))* log(xi)))/ (4.*mB*pow(-1 +
            sigma,3));
        }

        inline double t_23a_Diota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(),
                lambdaB2 = pow(lambdaB, 2),
                lambdaB3 = pow(lambdaB, 3);
            const auto mB = m_B(),
                mB2 = pow(mB, 2),
                mB3 = pow(mB, 3),
                mB4 = pow(mB, 4);
            const auto xi = mB * sigma / lambdaB;
            const auto Ei = gsl_sf_expint_Ei(-xi);

            return (exp(-xi)*((mB4*sigma)/lambdaB + (mB2*(19*lambdaB*mB -
            8*gamma*lambdaB*mB + mB2 + 10*q2 - 4*Ei*lambdaB*mB*exp(xi) -
            8*lambdaB*mB*log(xi)))/ lambdaB + (mB*(4*(5 - 2*gamma)*lambdaB* mB2
            + mB3 + 4*(29 - 13*gamma)*lambdaB* q2 + mB*(-24*(-1 + gamma)*
            lambdaB2 + 29*q2) - 4*lambdaB* (6*lambdaB*mB + 2*mB2 +
            13*q2)*log(xi)))/ (lambdaB*(-1 + sigma)) +
            (3*q2*(-104*gamma*lambdaB2 + 104*lambdaB*mB - 44*gamma*lambdaB*mB +
            9*mB2 + 4*Ei* (26*lambdaB2 - 15*lambdaB*mB + 2*mB2)* exp(xi) -
            2*lambdaB* (-26*lambdaB + 11*mB)* exp(xi)*log(lambdaB) +
            52*lambdaB2*exp(xi)* log(1/(mB*sigma)) - 22*lambdaB*mB*exp(xi)*
            log(1/(mB*sigma)) - 104*lambdaB2*log(xi) - 44*lambdaB*mB*log(xi) +
            52*lambdaB2*exp(xi)* log(xi) - 22*lambdaB*mB*exp(xi)* log(xi)))/
            (mB*pow(-1 + sigma,4)) + (q2*(-256*gamma*lambdaB3 +
            568*lambdaB2*mB - 296*gamma*lambdaB2*mB + 152*lambdaB*mB2 -
            44*gamma*lambdaB*mB2 + 9*mB3 + 8*Ei*lambdaB* (32*lambdaB2 -
            34*lambdaB*mB + 7*mB2)* exp(xi) - 32*lambdaB2* (-4*lambdaB + 3*mB)*
            exp(xi)*log(lambdaB) + 128*lambdaB3*exp(xi)* log(1/(mB*sigma)) -
            96*lambdaB2*mB*exp(xi)* log(1/(mB*sigma)) - 256*lambdaB3*log(xi) -
            296*lambdaB2*mB*log(xi) - 44*lambdaB*mB2*log(xi) +
            128*lambdaB3*exp(xi)* log(xi) - 96*lambdaB2*mB*exp(xi)* log(xi)))/
            (lambdaB*mB* pow(-1 + sigma,3)) + (-24*gamma*lambdaB3*mB +
            24*lambdaB2*mB2 - 8*gamma*lambdaB2*mB2 + lambdaB*mB3 +
            256*lambdaB2*q2 - 180*gamma*lambdaB2*q2 + 241*lambdaB*mB*q2 -
            96*gamma*lambdaB*mB*q2 + 28*mB2*q2 + 4*Ei*lambdaB* (6*lambdaB2*mB -
            4*lambdaB*mB2 + mB3 - 19*lambdaB*q2 + 8*mB*q2)* exp(xi) +
            2*lambdaB2* (6*lambdaB*mB - 2*mB2 - 13*q2)*exp(xi)* log(lambdaB) +
            12*lambdaB3*mB*exp(xi)* log(1/(mB*sigma)) - 4*lambdaB2*mB2*exp(xi)*
            log(1/(mB*sigma)) - 26*lambdaB2*q2*exp(xi)* log(1/(mB*sigma)) -
            24*lambdaB3*mB*log(xi) - 8*lambdaB2*mB2*log(xi) -
            180*lambdaB2*q2*log(xi) - 96*lambdaB*mB*q2*log(xi) +
            12*lambdaB3*mB*exp(xi)* log(xi) - 4*lambdaB2*mB2*exp(xi)* log(xi) -
            26*lambdaB2*q2*exp(xi)* log(xi))/ (lambdaB*pow(-1 + sigma,2))))
            /4.;
        }

        inline double t_23a_integrand(const double & q2, const double & sigma) const
        {
            const auto m_B2 = pow(m_B(), 2);
            const auto sigma2 = pow(sigma, 2);

            // two-particle contributions
            const auto c_p     = 1.0;
            const auto c_delta = -2.0 * q2 * sigma / (-q2 + m_B2 * pow(1.0 - sigma, 2));
            const auto c_bar   = 2.0 * q2 * (q2 + m_B2 * (sigma2 - 1.0)) / (m_B * pow(q2 - m_B2 * pow(sigma - 1.0, 2), 2));

            const auto phi_p     = b_lcdas.phi_plus(sigma * m_B());
            const auto phi_m     = b_lcdas.phi_minus(sigma * m_B());
            const auto Phi_bar   = b_lcdas.Phibar(sigma * m_B());
            const auto phi_delta = phi_p - phi_m;

            const auto result_2p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (c_p * phi_p + c_delta * phi_delta + c_bar * Phi_bar);

            // three-particle contributions
            const auto iota1 = t_23a_iota1(q2, sigma);
            const auto iota2 = t_23a_iota2(q2, sigma);
            const auto iota3 = t_23a_iota3(q2, sigma);

            const auto result_3p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (-iota1 + iota2 / M2 - iota3 / (2.0 * M2 * M2));

            return result_2p + result_3p;
        }

        inline double t_23a(const double & q2) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = std::pow(m_Kstar(), 2);

            const std::function<double (const double &)> integrand = std::bind(&Implementation<AnalyticFormFactorBToKstarKMO2006>::t_23a_integrand, this, q2, std::placeholders::_1);

            const auto sigma0  = this->sigma0(q2);
            const auto iota20  = t_23a_iota2(q2, sigma0);
            const auto iota30  = t_23a_iota3(q2, sigma0);
            const auto Diota30 = t_23a_Diota3(q2, sigma0);
            const auto etaf0   = etaf(q2, sigma0);
            const auto Detaf0  = Detaf(q2, sigma0);

            const auto integral = integrate<GSL::QNG>(integrand, 0.0, sigma0);
            const auto delta = (iota20 - (1.0 / M2 + Detaf0 / m_B2) / 2.0 * iota30 - etaf0 / (2.0 * m_B2) * Diota30)
                * std::exp(-s0 / M2) / m_B2 * etaf0;

            return f_B * m_B2 / (2.0 * f_Kstar * m_Kstar) * std::exp(m_Kstar2 / M2)
                * (integral + delta);
        }

        /* T_23B */

        inline double t_23b_iota1(const double & /*q2*/, const double & /*sigma*/) const
        {
            return 0.0;
        }

        inline double t_23b_iota2(const double & /*q2*/, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2);
            const auto mB = m_B(), mB2 = pow(mB, 2);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB, xi2 = pow(xi, 2);
            const auto Ei = gsl_sf_expint_Ei(-xi);

            const auto result = (lambdaB*exp(-xi)* (-44*gamma*lambdaB*mB +
            84*gamma*lambdaB*mB*sigma - 84*mB2*sigma2 + 35*gamma*mB2*sigma2 +
            44*lambdaB*mB*xi - 17*gamma*lambdaB*mB*xi - 15*lambdaB2*pow(xi,3)
            + lambdaB*mB*xi2 - Ei*lambdaB*mB* (-44 + 84*sigma + (27 -
            49*sigma)*xi + (-5 + 7*sigma)*xi2)*exp(xi) + lambdaB*mB* (-44 +
            84*sigma - 17*xi + 35*sigma*xi)*log(xi)))/ (4.*pow(mB,2)*pow(-1
            + sigma,3));

            return result;
        }

        inline double t_23b_iota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(), lambdaB2 = pow(lambdaB, 2), lambdaB3 = pow(lambdaB, 3), lambdaB4 = pow(lambdaB, 4);
            const auto mB = m_B(), mB2 = pow(mB, 2), mB3 = pow(mB, 3), mB4 = pow(mB, 4);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB, xi2 = pow(xi, 2);
            const auto Ei = gsl_sf_expint_Ei(-xi);

            const auto result = (sigma*exp(-xi)* (44*gamma*lambdaB2*mB2 -
            44*gamma*lambdaB2*q2 - 56*gamma*lambdaB2*mB2*sigma -
            44*lambdaB*mB3*sigma + 17*gamma*lambdaB*mB3*sigma +
            64*gamma*lambdaB2*q2*sigma + mB4*pow(sigma,3) +
            56*lambdaB*mB3*sigma2 - 21*gamma*lambdaB*mB3*sigma2 - 3*mB4*sigma2
            - 64*lambdaB*mB*q2*sigma2 + 26*gamma*lambdaB*mB*q2*sigma2 +
            44*lambdaB2*q2*xi - 17*gamma*lambdaB2*q2*xi -
            12*lambdaB4*pow(xi,3) + 4*gamma*lambdaB4*pow(xi,3) -
            (5*lambdaB3*q2*pow(xi,3))/ mB + 2*lambdaB4*pow(xi,4) +
            12*gamma*lambdaB4*xi2 + 3*lambdaB2*q2*xi2 - Ei*(lambdaB*mB3*sigma*
            (-27 + 35*sigma - 8*sigma2) + mB4*sigma2* (5 - 7*sigma + 2*sigma2)
            + mB2*(q2*(-5 + 6*sigma)* sigma2 + 4*lambdaB2* (11 - 14*sigma +
            3*sigma2)) + lambdaB2*q2* (-44 + 64*sigma + (27 - 38*sigma)*xi))*
            exp(xi) + lambdaB*(4*lambdaB*mB2* (11 - 14*sigma + 3*sigma2) +
            mB3*sigma* (17 - 21*sigma + 4*sigma2) + lambdaB*q2* (-44 +
            64*sigma + (-17 + 26*sigma)*xi))* log(xi)))/ (2.*mB*pow(-1 +
            sigma,4));

            return result;
        }

        inline double t_23b_Diota3(const double & q2, const double & sigma) const
        {
            static const double gamma = 0.577215664901532861;

            if (sigma < 1e-6)
                return 0.0; // error ~ -1e-11

            // three-particle contribution
            const auto lambdaB = lambda_B_p(),
                lambdaB2 = pow(lambdaB, 2),
                lambdaB3 = pow(lambdaB, 3),
                lambdaB4 = pow(lambdaB, 4),
                lambdaB5 = pow(lambdaB, 5),
                lambdaB6 = pow(lambdaB, 6);
            const auto mB = m_B(),
                mB2 = pow(mB, 2),
                mB3 = pow(mB, 3),
                mB4 = pow(mB, 4),
                mB5 = pow(mB, 5);
            const auto sigma2 = pow(sigma, 2);
            const auto xi = mB * sigma / lambdaB,
                xi2 = pow(xi, 2);
            const auto Ei = gsl_sf_expint_Ei(-xi);

            const auto result = (exp(-xi)*(-44*gamma*lambdaB3*mB2 + 44*gamma*lambdaB3*q2 -
            20*gamma*lambdaB3*mB2*sigma + 44*lambdaB2*mB3*sigma +
            10*gamma*lambdaB2*mB3*sigma + 4*gamma*lambdaB3*q2*sigma +
            87*lambdaB*mB4*pow(sigma,3) - 38*gamma*lambdaB*mB4* pow(sigma,3) -
            3*mB5*pow(sigma,3) + 128*lambdaB2*mB*q2* pow(sigma,3) -
            90*gamma*lambdaB2*mB*q2* pow(sigma,3) + 4*mB5*pow(sigma,4) +
            20*lambdaB2*mB3*sigma2 - 71*gamma*lambdaB2*mB3*sigma2 -
            30*lambdaB*mB4*sigma2 + 17*gamma*lambdaB*mB4*sigma2 -
            128*gamma*lambdaB3*q2*sigma2 - 4*lambdaB2*mB*q2*sigma2 +
            64*gamma*lambdaB2*mB*q2* sigma2 - 44*lambdaB3*q2*xi -
            10*gamma*lambdaB3*q2*xi - 76*lambdaB5*pow(xi,3) +
            73*gamma*lambdaB5*pow(xi,3) - (12*gamma*lambdaB6*pow(xi,3))/ mB +
            3*lambdaB3*q2* pow(xi,3) - (80*lambdaB4*q2*pow(xi,3))/ mB +
            (43*gamma*lambdaB4*q2* pow(xi,3))/mB - 69*lambdaB5*pow(xi,4) +
            25*gamma*lambdaB5*pow(xi,4) + (12*lambdaB6*pow(xi,4))/mB -
            (12*gamma*lambdaB6*pow(xi,4))/ mB + (58*lambdaB5*q2*
            pow(xi,4))/pow(mB,2) - (26*gamma*lambdaB5*q2* pow(xi,4))/pow(mB,2)
            - (8*lambdaB4*q2*pow(xi,4))/ mB + lambdaB5*pow(xi,5) +
              (12*lambdaB6*pow(xi,5))/mB - (4*gamma*lambdaB6*pow(xi,5))/ mB +
              (5*lambdaB5*q2* pow(xi,5))/pow(mB,2) - (2*lambdaB6*pow(xi,6))/mB
              + 76*gamma*lambdaB5*xi2 + 30*lambdaB3*q2*xi2 -
              17*gamma*lambdaB3*q2*xi2 + Ei*lambdaB* (3*lambdaB*mB3*sigma* (-18
              + 17*sigma + sigma2) + mB4*sigma2* (15 - 23*sigma -
              2*pow(sigma,3) + 10*sigma2) + mB2*(4*lambdaB2* (11 + 5*sigma +
              3*pow(sigma,3) - 19*sigma2) + q2*(-15 + 19*sigma)*sigma2) -
              2*lambdaB2*q2* (22 + 2*sigma - 64*sigma2 + (-27 + 30*sigma +
              19*sigma2)*xi))*exp(xi)\
        - lambdaB* (lambdaB*mB3*sigma* (-10 + 71*sigma + 12*pow(sigma,3) -
          73*sigma2) + mB4*pow(-1 + sigma,2)* (-17 + 4*sigma)*sigma2 + mB2*(-1
          + sigma)* (q2*(-17 + 26*sigma)* sigma2 + 4*lambdaB2* (-11 - 16*sigma
          + 3*sigma2)) + 2*lambdaB2*q2* (-22 - 2*sigma + 64*sigma2 + (5 -
          32*sigma + 45*sigma2)*xi))*log(xi))) /(2.*lambdaB*mB* pow(-1 +
          sigma,5));

            return result;
        }

        inline double t_23b_integrand(const double & q2, const double & sigma) const
        {
            const auto m_B2 = pow(m_B(), 2);

            // two-particle contributions
            const auto c_p     = sigma / (1.0 - sigma);
            const auto c_delta = (m_B2 * pow(-1.0 + sigma, 2) + q2 * (2.0 * sigma - 1.0)) * sigma
                               / ((1.0 - sigma) * (m_B2 * pow(1.0 - sigma, 2) - q2));
            const auto c_bar   = (m_B2 * pow(-1.0 + sigma, 2) + q2 * (2.0 * sigma - 1.0))
                               * (q2 + m_B2 * (-1.0 + pow(sigma, 2)))
                               / (m_B * pow(q2 - m_B2 * pow(-1.0 + sigma, 2), 2) * (-1.0 + sigma));

            const auto phi_p     = b_lcdas.phi_plus(sigma * m_B());
            const auto phi_m     = b_lcdas.phi_minus(sigma * m_B());
            const auto Phi_bar   = b_lcdas.Phibar(sigma * m_B());
            const auto phi_delta = phi_p - phi_m;

            const auto result_2p = -std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (c_p * phi_p + c_delta * phi_delta + c_bar * Phi_bar);

            // three-particle contributions
            const auto iota1 = t_23b_iota1(q2, sigma);
            const auto iota2 = t_23b_iota2(q2, sigma);
            const auto iota3 = t_23b_iota3(q2, sigma);

            const auto result_3p = std::exp(-(m_B2 * sigma - q2 * sigma / (1.0 - sigma)) / M2)
                    * (-iota1 + iota2 / M2 - iota3 / (2.0 * M2 * M2));

            return result_2p + result_3p;
        }

        inline double t_23b(const double & q2) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = std::pow(m_Kstar(), 2);

            const std::function<double (const double &)> integrand = std::bind(&Implementation<AnalyticFormFactorBToKstarKMO2006>::t_23b_integrand, this, q2, std::placeholders::_1);

            const auto sigma0  = this->sigma0(q2);
            const auto iota20  = t_23b_iota2(q2, sigma0);
            const auto iota30  = t_23b_iota3(q2, sigma0);
            const auto Diota30 = t_23b_Diota3(q2, sigma0);
            const auto etaf0   = etaf(q2, sigma0);
            const auto Detaf0  = Detaf(q2, sigma0);

            const auto integral = integrate<GSL::QNG>(integrand, 0.0, sigma0);
            const auto delta = (iota20 - (1.0 / M2 + Detaf0 / m_B2) / 2.0 * iota30 - etaf0 / (2.0 * m_B2) * Diota30)
                * std::exp(-s0 / M2) / m_B2 * etaf0;

            return f_B * m_B2 / (f_Kstar * m_Kstar) * std::exp(m_Kstar2 / M2)
                * (integral + delta);
        }

        /* T_2 */

        inline double t_2_iota2(const double & q2, const double & sigma) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = pow(m_Kstar(), 2);

            const auto c_23a = (m_B2 - m_Kstar2 - q2) / (m_B2 - m_Kstar2);
            const auto c_23b = 2.0 * q2 / (m_B2 - m_Kstar2);

            return c_23a * t_23a_iota2(q2, sigma)  + c_23b * t_23b_iota2(q2, sigma);
        }

        inline double t_2_iota3(const double & q2, const double & sigma) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = pow(m_Kstar(), 2);

            const auto c_23a = (m_B2 - m_Kstar2 - q2) / (m_B2 - m_Kstar2);
            const auto c_23b = 2.0 * q2 / (m_B2 - m_Kstar2);

            return c_23a * t_23a_iota3(q2, sigma)  + c_23b * t_23b_iota3(q2, sigma);
        }

        inline double t_2_Diota3(const double & q2, const double & sigma) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = pow(m_Kstar(), 2);

            const auto c_23a = (m_B2 - m_Kstar2 - q2) / (m_B2 - m_Kstar2);
            const auto c_23b = 2.0 * q2 / (m_B2 - m_Kstar2);

            return c_23a * t_23a_Diota3(q2, sigma)  + c_23b * t_23b_Diota3(q2, sigma);
        }

        inline double t_2_integrand(const double & q2, const double & sigma) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = pow(m_Kstar(), 2);

            const auto c_23a = (m_B2 - m_Kstar2 - q2) / (m_B2 - m_Kstar2);
            const auto c_23b = 2.0 * q2 / (m_B2 - m_Kstar2);

            return c_23a * t_23a_integrand(q2, sigma) + c_23b * t_23b_integrand(q2, sigma);
        }

        inline double t_2(const double & q2) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = pow(m_Kstar(), 2);

            const std::function<double (const double &)> integrand = std::bind(&Implementation<AnalyticFormFactorBToKstarKMO2006>::t_2_integrand, this, q2, std::placeholders::_1);

            const auto sigma0  = this->sigma0(q2);
            const auto iota20  = t_2_iota2(q2, sigma0);
            const auto iota30  = t_2_iota3(q2, sigma0);
            const auto Diota30 = t_2_Diota3(q2, sigma0);
            const auto etaf0   = etaf(q2, sigma0);
            const auto Detaf0  = Detaf(q2, sigma0);

            const auto integral = integrate<GSL::QNG>(integrand, 0.0, sigma0);
            const auto delta = (iota20 - (1.0 / M2 + Detaf0 / m_B2) / 2.0 * iota30 - etaf0 / (2.0 * m_B2) * Diota30)
                * std::exp(-s0 / M2) / m_B2 * etaf0;

            return f_B * m_B2 / (2.0 * f_Kstar * m_Kstar) * std::exp(m_Kstar2 / M2)
                * (integral + delta);
        }

        /* T_3 */

        inline double t_3_iota2(const double & q2, const double & sigma) const
        {
            const auto c_23a = 1.0;
            const auto c_23b = -2.0;

            return c_23a * t_23a_iota2(q2, sigma)  + c_23b * t_23b_iota2(q2, sigma);
        }

        inline double t_3_iota3(const double & q2, const double & sigma) const
        {
            const auto c_23a = 1.0;
            const auto c_23b = -2.0;

            return c_23a * t_23a_iota3(q2, sigma)  + c_23b * t_23b_iota3(q2, sigma);
        }

        inline double t_3_Diota3(const double & q2, const double & sigma) const
        {
            const auto c_23a = 1.0;
            const auto c_23b = -2.0;

            return c_23a * t_23a_Diota3(q2, sigma)  + c_23b * t_23b_Diota3(q2, sigma);
        }

        inline double t_3_integrand(const double & q2, const double & sigma) const
        {
            const auto c_23a = 1.0;
            const auto c_23b = -2.0;

            return c_23a * t_23a_integrand(q2, sigma) + c_23b * t_23b_integrand(q2, sigma);
        }

        inline double t_3(const double & q2) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = pow(m_Kstar(), 2);

            const std::function<double (const double &)> integrand = std::bind(&Implementation<AnalyticFormFactorBToKstarKMO2006>::t_3_integrand, this, q2, std::placeholders::_1);

            const auto sigma0  = this->sigma0(q2);
            const auto iota20  = t_3_iota2(q2, sigma0);
            const auto iota30  = t_3_iota3(q2, sigma0);
            const auto Diota30 = t_3_Diota3(q2, sigma0);
            const auto etaf0   = etaf(q2, sigma0);
            const auto Detaf0  = Detaf(q2, sigma0);

            const auto integral = integrate<GSL::QNG>(integrand, 0.0, sigma0);
            const auto delta = (iota20 - (1.0 / M2 + Detaf0 / m_B2) / 2.0 * iota30 - etaf0 / (2.0 * m_B2) * Diota30)
                * std::exp(-s0 / M2) / m_B2 * etaf0;

            return f_B * m_B2 / (2.0 * f_Kstar * m_Kstar) * std::exp(m_Kstar2 / M2)
                * (integral + delta);
        }

        inline double t_23(const double & q2) const
        {
            const auto m_B2     = pow(m_B(), 2);
            const auto m_Kstar2 = pow(m_Kstar(), 2);
            const auto lambda   = eos::lambda(m_B2, m_Kstar2, q2);

            return (m_B + m_Kstar) / (8.0 * m_B * m_Kstar2) * ((m_B2 + 3.0 * m_Kstar2 - q2) * this->t_2(q2) - lambda / (m_B2 - m_Kstar2) * this->t_3(q2));
        }

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            // V
            // skipping iota_1, which is identically 0
            results.add({ v_iota2(-5.0, 0.04),    "Iota_2V(q^2 = -5, sigma = 0.04)"});
            results.add({ v_iota2(-1.0, 0.04),    "Iota_2V(q^2 = -1, sigma = 0.04)"});
            results.add({ v_iota2( 0.0, 0.04),    "Iota_2V(q^2 =  0, sigma = 0.04)"});

            results.add({ v_iota2(-5.0, 0.06),    "Iota_2V(q^2 = -5, sigma = 0.06)"});
            results.add({ v_iota2(-1.0, 0.06),    "Iota_2V(q^2 = -1, sigma = 0.06)"});
            results.add({ v_iota2( 0.0, 0.06),    "Iota_2V(q^2 =  0, sigma = 0.06)"});

            results.add({ v_iota3(-5.0, 0.04),    "Iota_3V(q^2 = -5, sigma = 0.04)"});
            results.add({ v_iota3(-1.0, 0.04),    "Iota_3V(q^2 = -1, sigma = 0.04)"});
            results.add({ v_iota3( 0.0, 0.04),    "Iota_3V(q^2 =  0, sigma = 0.04)"});

            results.add({ v_iota3(-5.0, 0.06),    "Iota_3V(q^2 = -5, sigma = 0.06)"});
            results.add({ v_iota3(-1.0, 0.06),    "Iota_3V(q^2 = -1, sigma = 0.06)"});
            results.add({ v_iota3( 0.0, 0.06),    "Iota_3V(q^2 =  0, sigma = 0.06)"});

            results.add({ v_Diota3(-5.0, 0.04),   "DIota_3V(q^2 = -5, sigma = 0.04)"});
            results.add({ v_Diota3(-1.0, 0.04),   "DIota_3V(q^2 = -1, sigma = 0.04)"});
            results.add({ v_Diota3( 0.0, 0.04),   "DIota_3V(q^2 =  0, sigma = 0.04)"});

            results.add({ v_Diota3(-5.0, 0.06),   "DIota_3V(q^2 = -5, sigma = 0.06)"});
            results.add({ v_Diota3(-1.0, 0.06),   "DIota_3V(q^2 = -1, sigma = 0.06)"});
            results.add({ v_Diota3( 0.0, 0.06),   "DIota_3V(q^2 =  0, sigma = 0.06)"});

            // A_0
            // skipping iota_1, which is identically 0
            results.add({ a_0_iota2(-5.0, 0.04),  "Iota_2A0(q^2 = -5, sigma = 0.04)"});
            results.add({ a_0_iota2(-1.0, 0.04),  "Iota_2A0(q^2 = -1, sigma = 0.04)"});
            results.add({ a_0_iota2( 0.0, 0.04),  "Iota_2A0(q^2 =  0, sigma = 0.04)"});

            results.add({ a_0_iota2(-5.0, 0.06),  "Iota_2A0(q^2 = -5, sigma = 0.06)"});
            results.add({ a_0_iota2(-1.0, 0.06),  "Iota_2A0(q^2 = -1, sigma = 0.06)"});
            results.add({ a_0_iota2( 0.0, 0.06),  "Iota_2A0(q^2 =  0, sigma = 0.06)"});

            results.add({ a_0_iota3(-5.0, 0.04),  "Iota_3A0(q^2 = -5, sigma = 0.04)"});
            results.add({ a_0_iota3(-1.0, 0.04),  "Iota_3A0(q^2 = -1, sigma = 0.04)"});
            results.add({ a_0_iota3( 0.0, 0.04),  "Iota_3A0(q^2 =  0, sigma = 0.04)"});

            results.add({ a_0_iota3(-5.0, 0.06),  "Iota_3A0(q^2 = -5, sigma = 0.06)"});
            results.add({ a_0_iota3(-1.0, 0.06),  "Iota_3A0(q^2 = -1, sigma = 0.06)"});
            results.add({ a_0_iota3( 0.0, 0.06),  "Iota_3A0(q^2 =  0, sigma = 0.06)"});

            results.add({ a_0_Diota3(-5.0, 0.04), "DIota_3A0(q^2 = -5, sigma = 0.04)"});
            results.add({ a_0_Diota3(-1.0, 0.04), "DIota_3A0(q^2 = -1, sigma = 0.04)"});
            results.add({ a_0_Diota3( 0.0, 0.04), "DIota_3A0(q^2 =  0, sigma = 0.04)"});

            results.add({ a_0_Diota3(-5.0, 0.06), "DIota_3A0(q^2 = -5, sigma = 0.06)"});
            results.add({ a_0_Diota3(-1.0, 0.06), "DIota_3A0(q^2 = -1, sigma = 0.06)"});
            results.add({ a_0_Diota3( 0.0, 0.06), "DIota_3A0(q^2 =  0, sigma = 0.06)"});

            // A_1
            results.add({ a_1_iota1(-5.0, 0.04),  "Iota_1A1(q^2 = -5, sigma = 0.04)"});
            results.add({ a_1_iota1(-1.0, 0.04),  "Iota_1A1(q^2 = -1, sigma = 0.04)"});
            results.add({ a_1_iota1( 0.0, 0.04),  "Iota_1A1(q^2 =  0, sigma = 0.04)"});

            results.add({ a_1_iota1(-5.0, 0.06),  "Iota_1A1(q^2 = -5, sigma = 0.06)"});
            results.add({ a_1_iota1(-1.0, 0.06),  "Iota_1A1(q^2 = -1, sigma = 0.06)"});
            results.add({ a_1_iota1( 0.0, 0.06),  "Iota_1A1(q^2 =  0, sigma = 0.06)"});

            results.add({ a_1_iota2(-5.0, 0.04),  "Iota_2A1(q^2 = -5, sigma = 0.04)"});
            results.add({ a_1_iota2(-1.0, 0.04),  "Iota_2A1(q^2 = -1, sigma = 0.04)"});
            results.add({ a_1_iota2( 0.0, 0.04),  "Iota_2A1(q^2 =  0, sigma = 0.04)"});

            results.add({ a_1_iota2(-5.0, 0.06),  "Iota_2A1(q^2 = -5, sigma = 0.06)"});
            results.add({ a_1_iota2(-1.0, 0.06),  "Iota_2A1(q^2 = -1, sigma = 0.06)"});
            results.add({ a_1_iota2( 0.0, 0.06),  "Iota_2A1(q^2 =  0, sigma = 0.06)"});

            results.add({ a_1_iota3(-5.0, 0.04),  "Iota_3A1(q^2 = -5, sigma = 0.04)"});
            results.add({ a_1_iota3(-1.0, 0.04),  "Iota_3A1(q^2 = -1, sigma = 0.04)"});
            results.add({ a_1_iota3( 0.0, 0.04),  "Iota_3A1(q^2 =  0, sigma = 0.04)"});

            results.add({ a_1_iota3(-5.0, 0.06),  "Iota_3A1(q^2 = -5, sigma = 0.06)"});
            results.add({ a_1_iota3(-1.0, 0.06),  "Iota_3A1(q^2 = -1, sigma = 0.06)"});
            results.add({ a_1_iota3( 0.0, 0.06),  "Iota_3A1(q^2 =  0, sigma = 0.06)"});

            results.add({ a_1_Diota3(-5.0, 0.04), "DIota_3A1(q^2 = -5, sigma = 0.04)"});
            results.add({ a_1_Diota3(-1.0, 0.04), "DIota_3A1(q^2 = -1, sigma = 0.04)"});
            results.add({ a_1_Diota3( 0.0, 0.04), "DIota_3A1(q^2 =  0, sigma = 0.04)"});

            results.add({ a_1_Diota3(-5.0, 0.06), "DIota_3A1(q^2 = -5, sigma = 0.06)"});
            results.add({ a_1_Diota3(-1.0, 0.06), "DIota_3A1(q^2 = -1, sigma = 0.06)"});
            results.add({ a_1_Diota3( 0.0, 0.06), "DIota_3A1(q^2 =  0, sigma = 0.06)"});

            // A_2
            // skipping iota_1, which is identically 0
            results.add({ a_2_iota2(-5.0, 0.04),  "Iota_2A2(q^2 = -5, sigma = 0.04)"});
            results.add({ a_2_iota2(-1.0, 0.04),  "Iota_2A2(q^2 = -1, sigma = 0.04)"});
            results.add({ a_2_iota2( 0.0, 0.04),  "Iota_2A2(q^2 =  0, sigma = 0.04)"});

            results.add({ a_2_iota2(-5.0, 0.06),  "Iota_2A2(q^2 = -5, sigma = 0.06)"});
            results.add({ a_2_iota2(-1.0, 0.06),  "Iota_2A2(q^2 = -1, sigma = 0.06)"});
            results.add({ a_2_iota2( 0.0, 0.06),  "Iota_2A2(q^2 =  0, sigma = 0.06)"});

            results.add({ a_2_iota3(-5.0, 0.04),  "Iota_3A2(q^2 = -5, sigma = 0.04)"});
            results.add({ a_2_iota3(-1.0, 0.04),  "Iota_3A2(q^2 = -1, sigma = 0.04)"});
            results.add({ a_2_iota3( 0.0, 0.04),  "Iota_3A2(q^2 =  0, sigma = 0.04)"});

            results.add({ a_2_iota3(-5.0, 0.06),  "Iota_3A2(q^2 = -5, sigma = 0.06)"});
            results.add({ a_2_iota3(-1.0, 0.06),  "Iota_3A2(q^2 = -1, sigma = 0.06)"});
            results.add({ a_2_iota3( 0.0, 0.06),  "Iota_3A2(q^2 =  0, sigma = 0.06)"});

            results.add({ a_2_Diota3(-5.0, 0.04), "DIota_3A2(q^2 = -5, sigma = 0.04)"});
            results.add({ a_2_Diota3(-1.0, 0.04), "DIota_3A2(q^2 = -1, sigma = 0.04)"});
            results.add({ a_2_Diota3( 0.0, 0.04), "DIota_3A2(q^2 =  0, sigma = 0.04)"});

            results.add({ a_2_Diota3(-5.0, 0.06), "DIota_3A2(q^2 = -5, sigma = 0.06)"});
            results.add({ a_2_Diota3(-1.0, 0.06), "DIota_3A2(q^2 = -1, sigma = 0.06)"});
            results.add({ a_2_Diota3( 0.0, 0.06), "DIota_3A2(q^2 =  0, sigma = 0.06)"});

            // T_1
            // skipping iota_1, which is identically 0
            results.add({ t_1_iota2(-5.0, 0.04),  "Iota_2T1(q^2 = -5, sigma = 0.04)"});
            results.add({ t_1_iota2(-1.0, 0.04),  "Iota_2T1(q^2 = -1, sigma = 0.04)"});
            results.add({ t_1_iota2( 0.0, 0.04),  "Iota_2T1(q^2 =  0, sigma = 0.04)"});

            results.add({ t_1_iota2(-5.0, 0.06),  "Iota_2T1(q^2 = -5, sigma = 0.06)"});
            results.add({ t_1_iota2(-1.0, 0.06),  "Iota_2T1(q^2 = -1, sigma = 0.06)"});
            results.add({ t_1_iota2( 0.0, 0.06),  "Iota_2T1(q^2 =  0, sigma = 0.06)"});

            results.add({ t_1_iota3(-5.0, 0.04),  "Iota_3T1(q^2 = -5, sigma = 0.04)"});
            results.add({ t_1_iota3(-1.0, 0.04),  "Iota_3T1(q^2 = -1, sigma = 0.04)"});
            results.add({ t_1_iota3( 0.0, 0.04),  "Iota_3T1(q^2 =  0, sigma = 0.04)"});

            results.add({ t_1_iota3(-5.0, 0.06),  "Iota_3T1(q^2 = -5, sigma = 0.06)"});
            results.add({ t_1_iota3(-1.0, 0.06),  "Iota_3T1(q^2 = -1, sigma = 0.06)"});
            results.add({ t_1_iota3( 0.0, 0.06),  "Iota_3T1(q^2 =  0, sigma = 0.06)"});

            results.add({ t_1_Diota3(-5.0, 0.04), "DIota_3T1(q^2 = -5, sigma = 0.04)"});
            results.add({ t_1_Diota3(-1.0, 0.04), "DIota_3T1(q^2 = -1, sigma = 0.04)"});
            results.add({ t_1_Diota3( 0.0, 0.04), "DIota_3T1(q^2 =  0, sigma = 0.04)"});

            results.add({ t_1_Diota3(-5.0, 0.06), "DIota_3T1(q^2 = -5, sigma = 0.06)"});
            results.add({ t_1_Diota3(-1.0, 0.06), "DIota_3T1(q^2 = -1, sigma = 0.06)"});
            results.add({ t_1_Diota3( 0.0, 0.06), "DIota_3T1(q^2 =  0, sigma = 0.06)"});

            // T_23A
            // skipping iota_1, which is identically 0
            results.add({ t_23a_iota2(-5.0, 0.04),  "Iota_2T23A(q^2 = -5, sigma = 0.04)"});
            results.add({ t_23a_iota2(-1.0, 0.04),  "Iota_2T23A(q^2 = -1, sigma = 0.04)"});
            results.add({ t_23a_iota2( 0.0, 0.04),  "Iota_2T23A(q^2 =  0, sigma = 0.04)"});

            results.add({ t_23a_iota2(-5.0, 0.06),  "Iota_2T23A(q^2 = -5, sigma = 0.06)"});
            results.add({ t_23a_iota2(-1.0, 0.06),  "Iota_2T23A(q^2 = -1, sigma = 0.06)"});
            results.add({ t_23a_iota2( 0.0, 0.06),  "Iota_2T23A(q^2 =  0, sigma = 0.06)"});

            results.add({ t_23a_iota3(-5.0, 0.04),  "Iota_3T23A(q^2 = -5, sigma = 0.04)"});
            results.add({ t_23a_iota3(-1.0, 0.04),  "Iota_3T23A(q^2 = -1, sigma = 0.04)"});
            results.add({ t_23a_iota3( 0.0, 0.04),  "Iota_3T23A(q^2 =  0, sigma = 0.04)"});

            results.add({ t_23a_iota3(-5.0, 0.06),  "Iota_3T23A(q^2 = -5, sigma = 0.06)"});
            results.add({ t_23a_iota3(-1.0, 0.06),  "Iota_3T23A(q^2 = -1, sigma = 0.06)"});
            results.add({ t_23a_iota3( 0.0, 0.06),  "Iota_3T23A(q^2 =  0, sigma = 0.06)"});

            results.add({ t_23a_Diota3(-5.0, 0.04), "DIota_3T23A(q^2 = -5, sigma = 0.04)"});
            results.add({ t_23a_Diota3(-1.0, 0.04), "DIota_3T23A(q^2 = -1, sigma = 0.04)"});
            results.add({ t_23a_Diota3( 0.0, 0.04), "DIota_3T23A(q^2 =  0, sigma = 0.04)"});

            results.add({ t_23a_Diota3(-5.0, 0.06), "DIota_3T23A(q^2 = -5, sigma = 0.06)"});
            results.add({ t_23a_Diota3(-1.0, 0.06), "DIota_3T23A(q^2 = -1, sigma = 0.06)"});
            results.add({ t_23a_Diota3( 0.0, 0.06), "DIota_3T23A(q^2 =  0, sigma = 0.06)"});

            results.add({ t_23a(-5.0), "T_23A(q^2 = -5)"});
            results.add({ t_23a(-1.0), "T_23A(q^2 = -1)"});
            results.add({ t_23a( 0.0), "T_23A(q^2 =  0)"});

            // T_23B
            // skipping iota_1, which is identically 0
            results.add({ t_23b_iota2(-5.0, 0.04),  "Iota_2T23B(q^2 = -5, sigma = 0.04)"});
            results.add({ t_23b_iota2(-1.0, 0.04),  "Iota_2T23B(q^2 = -1, sigma = 0.04)"});
            results.add({ t_23b_iota2( 0.0, 0.04),  "Iota_2T23B(q^2 =  0, sigma = 0.04)"});

            results.add({ t_23b_iota2(-5.0, 0.06),  "Iota_2T23B(q^2 = -5, sigma = 0.06)"});
            results.add({ t_23b_iota2(-1.0, 0.06),  "Iota_2T23B(q^2 = -1, sigma = 0.06)"});
            results.add({ t_23b_iota2( 0.0, 0.06),  "Iota_2T23B(q^2 =  0, sigma = 0.06)"});

            results.add({ t_23b_iota3(-5.0, 0.04),  "Iota_3T23B(q^2 = -5, sigma = 0.04)"});
            results.add({ t_23b_iota3(-1.0, 0.04),  "Iota_3T23B(q^2 = -1, sigma = 0.04)"});
            results.add({ t_23b_iota3( 0.0, 0.04),  "Iota_3T23B(q^2 =  0, sigma = 0.04)"});

            results.add({ t_23b_iota3(-5.0, 0.06),  "Iota_3T23B(q^2 = -5, sigma = 0.06)"});
            results.add({ t_23b_iota3(-1.0, 0.06),  "Iota_3T23B(q^2 = -1, sigma = 0.06)"});
            results.add({ t_23b_iota3( 0.0, 0.06),  "Iota_3T23B(q^2 =  0, sigma = 0.06)"});

            results.add({ t_23b_Diota3(-5.0, 0.04), "DIota_3T23B(q^2 = -5, sigma = 0.04)"});
            results.add({ t_23b_Diota3(-1.0, 0.04), "DIota_3T23B(q^2 = -1, sigma = 0.04)"});
            results.add({ t_23b_Diota3( 0.0, 0.04), "DIota_3T23B(q^2 =  0, sigma = 0.04)"});

            results.add({ t_23b_Diota3(-5.0, 0.06), "DIota_3T23B(q^2 = -5, sigma = 0.06)"});
            results.add({ t_23b_Diota3(-1.0, 0.06), "DIota_3T23B(q^2 = -1, sigma = 0.06)"});
            results.add({ t_23b_Diota3( 0.0, 0.06), "DIota_3T23B(q^2 =  0, sigma = 0.06)"});

            results.add({ t_23b(-5.0), "T_23B(q^2 = -5)"});
            results.add({ t_23b(-1.0), "T_23B(q^2 = -1)"});
            results.add({ t_23b( 0.0), "T_23B(q^2 =  0)"});

            // T_2
            // skipping iota_1, which is identically 0
            results.add({ t_2_iota2(-5.0, 0.04),  "Iota_2T2(q^2 = -5, sigma = 0.04)"});
            results.add({ t_2_iota2(-1.0, 0.04),  "Iota_2T2(q^2 = -1, sigma = 0.04)"});
            results.add({ t_2_iota2( 0.0, 0.04),  "Iota_2T2(q^2 =  0, sigma = 0.04)"});

            results.add({ t_2_iota2(-5.0, 0.06),  "Iota_2T2(q^2 = -5, sigma = 0.06)"});
            results.add({ t_2_iota2(-1.0, 0.06),  "Iota_2T2(q^2 = -1, sigma = 0.06)"});
            results.add({ t_2_iota2( 0.0, 0.06),  "Iota_2T2(q^2 =  0, sigma = 0.06)"});

            results.add({ t_2_iota3(-5.0, 0.04),  "Iota_3T2(q^2 = -5, sigma = 0.04)"});
            results.add({ t_2_iota3(-1.0, 0.04),  "Iota_3T2(q^2 = -1, sigma = 0.04)"});
            results.add({ t_2_iota3( 0.0, 0.04),  "Iota_3T2(q^2 =  0, sigma = 0.04)"});

            results.add({ t_2_iota3(-5.0, 0.06),  "Iota_3T2(q^2 = -5, sigma = 0.06)"});
            results.add({ t_2_iota3(-1.0, 0.06),  "Iota_3T2(q^2 = -1, sigma = 0.06)"});
            results.add({ t_2_iota3( 0.0, 0.06),  "Iota_3T2(q^2 =  0, sigma = 0.06)"});

            results.add({ t_2_Diota3(-5.0, 0.04), "DIota_3T2(q^2 = -5, sigma = 0.04)"});
            results.add({ t_2_Diota3(-1.0, 0.04), "DIota_3T2(q^2 = -1, sigma = 0.04)"});
            results.add({ t_2_Diota3( 0.0, 0.04), "DIota_3T2(q^2 =  0, sigma = 0.04)"});

            results.add({ t_2_Diota3(-5.0, 0.06), "DIota_3T2(q^2 = -5, sigma = 0.06)"});
            results.add({ t_2_Diota3(-1.0, 0.06), "DIota_3T2(q^2 = -1, sigma = 0.06)"});
            results.add({ t_2_Diota3( 0.0, 0.06), "DIota_3T2(q^2 =  0, sigma = 0.06)"});

            // T_3
            // skipping iota_1, which is identically 0
            results.add({ t_3_iota2(-5.0, 0.04),  "Iota_2T3(q^2 = -5, sigma = 0.04)"});
            results.add({ t_3_iota2(-1.0, 0.04),  "Iota_2T3(q^2 = -1, sigma = 0.04)"});
            results.add({ t_3_iota2( 0.0, 0.04),  "Iota_2T3(q^2 =  0, sigma = 0.04)"});

            results.add({ t_3_iota2(-5.0, 0.06),  "Iota_2T3(q^2 = -5, sigma = 0.06)"});
            results.add({ t_3_iota2(-1.0, 0.06),  "Iota_2T3(q^2 = -1, sigma = 0.06)"});
            results.add({ t_3_iota2( 0.0, 0.06),  "Iota_2T3(q^2 =  0, sigma = 0.06)"});

            results.add({ t_3_iota3(-5.0, 0.04),  "Iota_3T3(q^2 = -5, sigma = 0.04)"});
            results.add({ t_3_iota3(-1.0, 0.04),  "Iota_3T3(q^2 = -1, sigma = 0.04)"});
            results.add({ t_3_iota3( 0.0, 0.04),  "Iota_3T3(q^2 =  0, sigma = 0.04)"});

            results.add({ t_3_iota3(-5.0, 0.06),  "Iota_3T3(q^2 = -5, sigma = 0.06)"});
            results.add({ t_3_iota3(-1.0, 0.06),  "Iota_3T3(q^2 = -1, sigma = 0.06)"});
            results.add({ t_3_iota3( 0.0, 0.06),  "Iota_3T3(q^2 =  0, sigma = 0.06)"});

            results.add({ t_3_Diota3(-5.0, 0.04), "DIota_3T3(q^2 = -5, sigma = 0.04)"});
            results.add({ t_3_Diota3(-1.0, 0.04), "DIota_3T3(q^2 = -1, sigma = 0.04)"});
            results.add({ t_3_Diota3( 0.0, 0.04), "DIota_3T3(q^2 =  0, sigma = 0.04)"});

            results.add({ t_3_Diota3(-5.0, 0.06), "DIota_3T3(q^2 = -5, sigma = 0.06)"});
            results.add({ t_3_Diota3(-1.0, 0.06), "DIota_3T3(q^2 = -1, sigma = 0.06)"});
            results.add({ t_3_Diota3( 0.0, 0.06), "DIota_3T3(q^2 =  0, sigma = 0.06)"});

            return results;
        }
    };

    AnalyticFormFactorBToKstarKMO2006::AnalyticFormFactorBToKstarKMO2006(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<AnalyticFormFactorBToKstarKMO2006>(new Implementation<AnalyticFormFactorBToKstarKMO2006>(p, o, *this))
    {
    }

    AnalyticFormFactorBToKstarKMO2006::~AnalyticFormFactorBToKstarKMO2006()
    {
    }

    FormFactors<PToV> *
    AnalyticFormFactorBToKstarKMO2006::make(const Parameters & p, unsigned)
    {
        return new AnalyticFormFactorBToKstarKMO2006(p, Options{});
    }

    double
    AnalyticFormFactorBToKstarKMO2006::v(const double & q2) const
    {
        return _imp->v(q2);
    }

    double
    AnalyticFormFactorBToKstarKMO2006::a_0(const double & q2) const
    {
        return _imp->a_0(q2);
    }

    double
    AnalyticFormFactorBToKstarKMO2006::a_1(const double & q2) const
    {
        return _imp->a_1(q2);
    }

    double
    AnalyticFormFactorBToKstarKMO2006::a_2(const double & q2) const
    {
        return _imp->a_2(q2);
    }

    double
    AnalyticFormFactorBToKstarKMO2006::a_12(const double & q2) const
    {
        return _imp->a_12(q2);
    }

    double
    AnalyticFormFactorBToKstarKMO2006::t_1(const double & q2) const
    {
        return _imp->t_1(q2);
    }

    double
    AnalyticFormFactorBToKstarKMO2006::t_2(const double & q2) const
    {
        return _imp->t_2(q2);
    }

    double
    AnalyticFormFactorBToKstarKMO2006::t_3(const double & q2) const
    {
        return _imp->t_3(q2);
    }

    double
    AnalyticFormFactorBToKstarKMO2006::t_23(const double & q2) const
    {
        return _imp->t_23(q2);
    }

    Diagnostics
    AnalyticFormFactorBToKstarKMO2006::diagnostics() const
    {
        return _imp->diagnostics();
    }
}
