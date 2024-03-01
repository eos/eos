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

#include <eos/form-factors/k-star-lcdas.hh>
#include <eos/models/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    template <>
    struct Implementation<AntiKStarLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 (even) para Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a1para_0;
        UsedParameter a2para_0;
        UsedParameter a3para_0;
        UsedParameter a4para_0;
        UsedParameter fpara;

        // twist 2 (tensor) Gegenbauer coefficients and normalization at mu = 1 GeV
        UsedParameter a1perp_0;
        UsedParameter a2perp_0;
        UsedParameter a3perp_0;
        UsedParameter a4perp_0;
        UsedParameter fperp_0;

        // twist 3 LCDA parameters at mu = 1 Gev
        UsedParameter zeta3para_0;
        UsedParameter lambda3paratilde_0;
        UsedParameter omega3paratilde_0;
        UsedParameter kappa3para_0;
        UsedParameter omega3para_0;
        UsedParameter lambda3para_0;
        UsedParameter kappa3perp_0;
        UsedParameter omega3perp_0;
        UsedParameter lambda3perp_0;

        // twist 4 LCDA parameters at mu = 1 Gev
        UsedParameter zeta4para_0;
        UsedParameter omega4paratilde_0;
        UsedParameter zeta4perp_0;
        UsedParameter zeta4perptilde_0;

        // Kstar mass
        UsedParameter M_V;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a1para_0(p["K^*::a1para@1GeV"], u),
            a2para_0(p["K^*::a2para@1GeV"], u),
            a3para_0(p["K^*::a3para@1GeV"], u),
            a4para_0(p["K^*::a4para@1GeV"], u),
            fpara(p["K^*::fpara"], u),
            a1perp_0(p["K^*::a1perp@1GeV"], u),
            a2perp_0(p["K^*::a2perp@1GeV"], u),
            a3perp_0(p["K^*::a3perp@1GeV"], u),
            a4perp_0(p["K^*::a4perp@1GeV"], u),
            fperp_0(p["K^*::fperp@1GeV"], u),
            zeta3para_0(p["K^*::zeta3para@1GeV"], u),
            lambda3paratilde_0(p["K^*::lambda3paratilde@1GeV"], u),
            omega3paratilde_0(p["K^*::omega3paratilde@1GeV"], u),
            kappa3para_0(p["K^*::kappa3para@1GeV"], u),
            omega3para_0(p["K^*::omega3para@1GeV"], u),
            lambda3para_0(p["K^*::lambda3para@1GeV"], u),
            kappa3perp_0(p["K^*::kappa3perp@1GeV"], u),
            omega3perp_0(p["K^*::omega3perp@1GeV"], u),
            lambda3perp_0(p["K^*::lambda3perp@1GeV"], u),
            zeta4para_0(p["K^*::zeta4para@1GeV"], u),
            omega4paratilde_0(p["K^*::omega4paratilde@1GeV"], u),
            zeta4perp_0(p["K^*::zeta4perp@1GeV"], u),
            zeta4perptilde_0(p["K^*::zeta4perptilde@1GeV"], u),
            M_V(p["mass::K_u^*"], u),
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

            throw InternalError("Implementation<AntiKStarLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        inline double a1para(const double & mu) const
        {
            return a1para_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }

        inline double a2para(const double & mu) const
        {
            return a2para_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double a3para(const double & mu) const
        {
            return a3para_0 * std::pow(c_rge(mu), 314.0 / 45.0);
        }

        inline double a4para(const double & mu) const
        {
            return a4para_0 * std::pow(c_rge(mu), 364.0 / 45.0);
        }

        inline double a1perp(const double & mu) const
        {
            return a1perp_0 * std::pow(c_rge(mu), 36.0 / 9.0);
        }

        inline double a2perp(const double & mu) const
        {
            return a2perp_0 * std::pow(c_rge(mu), 52.0 / 9.0);
        }

        inline double a3perp(const double & mu) const
        {
            return a3perp_0 * std::pow(c_rge(mu), 64.0 / 9.0);
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

        // running of twist 3 parameters
        const double ms0 = model->m_s_msbar(1.0);
        const double mq0 = model->m_ud_msbar(1.0) / 2.0;

        inline double zeta3para(const double & mu) const
        {
            return (zeta3para_0 * std::pow(c_rge(mu), 77.0 / 9.0) * fpara +
                (6.0 * a1perp_0 * (-1.0 + std::pow(c_rge(mu), 5.0 / 9.0)) * std::pow(c_rge(mu), 8.0) * fperp_0 * (mq0 - ms0)) / (25.0 * M_V) +
                (2.0 * (std::pow(c_rge(mu), 16.0 / 3.0) - std::pow(c_rge(mu), 77.0 / 9.0)) * fperp_0 * (mq0 + ms0)) / (29.0 * M_V)) / fpara;
        }

        inline double kappa3para(const double & mu) const
        {
            return (std::pow(c_rge(mu), 77.0 / 9.0) * fpara * kappa3para_0 -
                (2.0 * std::pow(c_rge(mu), 16.0 / 3.0) * (-1.0 + std::pow(c_rge(mu), 29.0 / 9.0)) * fperp_0 * (mq0 - ms0)) / (29.0 * M_V) +
                (6.0 * a1perp_0 * (-1.0 + std::pow(c_rge(mu), 5.0 / 9.0)) * std::pow(c_rge(mu), 8.0) * fperp_0 * (mq0 + ms0)) / (25.0 * M_V)) / fpara;
        }

        inline double kappa3perp(const double & mu) const
        {
            return (std::pow(c_rge(mu), 55.0 / 9.0) * fperp_0 * kappa3perp_0 -
                (4.0 * std::pow(c_rge(mu), 4.0) * (-1.0 + std::pow(c_rge(mu), 19.0 / 9.0)) * fpara * (mq0 - ms0)) / (19.0 * M_V) +
                (12.0 * a1para_0 * std::pow(c_rge(mu), 55.0 / 9.0) * (-1.0 + std::pow(c_rge(mu), 13.0 / 9.0)) * fpara * (mq0 + ms0)) / (65.0 * M_V)) / fperp(mu);
        }

        inline double omega3perp(const double & mu) const
        {
            return ((-42.0 * a1para_0 * (-1.0 + std::pow(c_rge(mu), 5.0 / 9.0)) * std::pow(c_rge(mu), 68.0 / 9.0) * fpara * (mq0 - ms0)) / (25.0 * M_V) +
                (12.0 * a2para_0 * std::pow(c_rge(mu), 73.0 / 9.0) * (-1 + std::pow(c_rge(mu), 13.0 / 9.0)) * fpara * (mq0 + ms0)) / (13.0 * M_V) +
                (14.0 * (std::pow(c_rge(mu), 4.0) - std::pow(c_rge(mu), 73.0 / 9.0)) * fpara *(mq0 + ms0)) / (37.0 * M_V) +
                std::pow(c_rge(mu), 73.0 / 9.0) * fperp_0 * omega3perp_0) / fperp(mu);
        }

        inline double lambda3perp(const double & mu) const
        {
            return (std::pow(c_rge(mu), 4.0) * (3.0 * fpara * (mq0 - ms0) +
                68.0 * a2para_0 * std::pow(c_rge(mu), 5.555555555555555) * fpara * (mq0 - ms0) -
                51.0 * a1para_0 * std::pow(c_rge(mu), 3.5555555555555554) * fpara * (mq0 + ms0) +
                std::pow(c_rge(mu), 7.555555555555555) * (255.0 * fperp_0 * lambda3perp_0 * M_V +
                (-3.0 + 51.0 * a1para_0 - 68.0 * a2para_0) * fpara * mq0 +
                (3.0 + 51.0 * a1para_0 + 68.0 * a2para_0) * fpara * ms0))) / (255.0 * fperp(mu) * M_V);
        }

        inline double omega3para(const double & mu) const
        {
            return -3.7780044580452606e-7 * (887490.0 * a1perp_0 * std::pow(c_rge(mu), 8.0 + std::sqrt(865) / 18.0) * fperp_0 * (ms0 - mq0) -
                448070.0 * std::pow(c_rge(mu), (96.0 + std::sqrt(865.0)) / 18.0) * fperp_0 * (ms0 + mq0) -
                89465220.0 * a2perp_0 * std::pow(c_rge(mu), (176.0 + std::sqrt(865)) / 18.0) * fperp_0 * (ms0 + mq0) +
                std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0) * (fperp_0 * ((224035.0 - 6811.0 * std::sqrt(865.0) +
                27.0 * (-16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 -
                306.0 * (-146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * ms0 +
                (224035.0 - 6811.0 * std::sqrt(865.0) -
                27.0 * (-16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 -
                306.0 * (-146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * mq0) +
                765.0 * fpara * M_V * (2.0 * (-865.0 + 26.0 * std::sqrt(865.0)) * omega3para_0 -
                63.0 * std::sqrt(865.0) * omega3paratilde_0)) +
                std::pow(c_rge(mu), 11.38888888888889) * (fperp_0 * ((7.0 * (32005.0 + 973.0 * std::sqrt(865.0)) -
                27.0 * (16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 +
                306.0 * (146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * ms0 +
                (7.0 * (32005.0 + 973.0 * std::sqrt(865.0)) +
                27.0 * (16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 +
                306.0 * (146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * mq0) -
                765.0 * fpara * M_V * (2.0 * (865.0 + 26.0 * std::sqrt(865.0)) * omega3para_0 -
                63.0 * std::sqrt(865.0) * omega3paratilde_0))) / (std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V);
        }

        inline double omega3paratilde(const double & mu) const
        {
            return (9.0 * a1perp_0 * ((865.0 - 7.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889) - 1730.0 * std::pow(c_rge(mu), 8.0 + std::sqrt(865.0) / 18.0) +
                (865.0 + 7.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0)) * fperp_0 * (mq0 - ms0)) / (147050.0 * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V) -
                ((-((9515.0 + 6707.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889)) + (-9515.0 + 6707.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0) +
                19030.0 * std::pow(c_rge(mu), (96.0 + std::sqrt(865.0)) / 18.0)) * fperp_0 * (mq0 + ms0)) / (1.191105e7 * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V) -
                (3.0 * a2perp_0 * (-((2595.0 + 91.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889)) + (-2595.0 + 91.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0) +
                5190.0 * std::pow(c_rge(mu), (176.0 + std::sqrt(865.0)) / 18.0)) * fperp_0 * (mq0 + ms0)) / (4325.0 * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V) +
                (3.0 * (-std::pow(c_rge(mu), 11.38888888888889) + std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0)) * omega3para_0) / (std::sqrt(865.0) * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0)) +
                (((865.0 - 26.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889) + (865.0 + 26.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0)) * omega3paratilde_0) /
                (1730.0 * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0));
        }

        inline double lambda3para(const double & mu) const
        {
            return (19030.0 * std::pow(c_rge(mu), (96.0 + std::sqrt(865.0)) / 18.0) * fperp_0 * (ms0 - mq0) +
                42879780.0 * a2perp_0 * std::pow(c_rge(mu), (176 + std::sqrt(865)) / 18.0) * fperp_0 *(ms0 - mq0) -
                1261170.0 * a1perp_0 * std::pow(c_rge(mu), 8.0 + std::sqrt(865.0) / 18.0) * fperp_0 * (ms0 + mq0) +
                std::pow(c_rge(mu), 11.38888888888889) * (fperp_0 * (-((9515.0 + 6707.0 * std::sqrt(865.0) +
                729.0 * (-865.0 + 7.0 * std::sqrt(865.0)) * a1perp_0 + 8262.0 * (2595.0 + 91.0 * std::sqrt(865.0)) * a2perp_0) * ms0) +
                (9515.0 + 6707.0 * std::sqrt(865.0) - 729.0 * (-865.0 + 7.0 * std::sqrt(865.0)) * a1perp_0 +
                8262.0 * (2595.0 + 91.0 * std::sqrt(865.0)) * a2perp_0) * mq0) -
                6885.0 * fpara * ((-865.0 + 26.0 * std::sqrt(865.0)) * lambda3para_0 +
                6.0 * std::sqrt(865.0) * lambda3paratilde_0) * M_V) +
                std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0) * (fperp_0 * ((-9515.0 + 6707.0 * std::sqrt(865.0) + 729.0 * (865.0 + 7.0 * std::sqrt(865.0)) * a1perp_0 +
                8262.0 * (-2595.0 + 91.0 * std::sqrt(865.0)) * a2perp_0) * ms0 +
                (9515.0 - 6707.0 * std::sqrt(865.0) + 729.0 * (865.0 + 7.0 * std::sqrt(865.0)) * a1perp_0 - 8262.0 *(-2595.0 + 91.0 * std::sqrt(865.0)) * a2perp_0) * mq0) +
                6885.0 * fpara * ((865.0 + 26.0 * std::sqrt(865.0)) * lambda3para_0 +
                6.0 *std::sqrt(865.0) * lambda3paratilde_0) * M_V)) / (1.191105e7 * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V);
        }

        inline double lambda3paratilde(const double & mu) const
        {
            return -3.7780044580452606e-7 * (448070.0 * std::pow(c_rge(mu), (96.0 + std::sqrt(865.0)) / 18.0) * fperp_0 * (ms0 - mq0) +
                89465220.0 * a2perp_0 * std::pow(c_rge(mu), (176.0 + std::sqrt(865.0)) / 18.0) * fperp_0 * (ms0 - mq0) -
                887490.0 * a1perp_0 * std::pow(c_rge(mu), 8.0 + std::sqrt(865.0) / 18.0) * fperp_0 * (ms0 + mq0) +
                std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0) * (fperp_0 * ((7.0 * (-32005.0 + 973.0 * std::sqrt(865.0)) -
                27.0 * (-16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 +
                306.0 * (-146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * ms0 + (224035.0 - 6811.0 * std::sqrt(865.0) - 27.0 *(-16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 -
                306.0 * (-146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * mq0) - 765.0 * fpara * (63.0 * std::sqrt(865.0) * lambda3para_0 +
                2.0 * (865.0 - 26.0 * std::sqrt(865.0)) * lambda3paratilde_0) * M_V) +
                std::pow(c_rge(mu), 11.38888888888889) * (fperp_0 * ((-7.0 * (32005.0 + 973.0 * std::sqrt(865.0)) + 27.0 * (16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 -
                306.0 * (146185.0 + 4961.0 * std::sqrt(865.0)) *a2perp_0) * ms0 +
                (7.0 * (32005.0 + 973.0 * std::sqrt(865.0)) + 27.0 * (16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 +
                306.0 * (146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * mq0) +
                765.0 * fpara * (63.0 * std::sqrt(865.0) * lambda3para_0 -
                2.0 * (865.0 + 26.0 * std::sqrt(865.0)) * lambda3paratilde_0) * M_V)) / (std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V);
        }

        // running of twist 4 parameters
        inline double zeta4para(const double & mu) const
        {
            return zeta4para_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }
        inline double omega4paratilde(const double & mu) const
        {
            return omega4paratilde_0 * std::pow(c_rge(mu), 10.0);
        }
        // mass corrections for running of zeta4perp, zeta4perptilde unknown
        inline double zeta4perp(const double & mu) const
        {
            return  1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) + std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perp_0 +
                1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) - std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perptilde_0;
        }
        inline double zeta4perptilde(const double & mu) const
        {
            return 1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) - std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perp_0 +
                1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) + std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perptilde_0;
        }
        inline double kappa4para(const double & mu) const
        {
            return -3.0 / 20.0 * std::pow(c_rge(mu), 32.0 / 9.0) * a1para_0 -
                std::pow(c_rge(mu), 16.0 / 3.0) * fperp_0 / fpara * (ms0 - mq0) / (4.0 * M_V) +
                std::pow(c_rge(mu), 8.0) * (std::pow(ms0, 2.0) - std::pow(mq0, 2.0)) / (2.0 * std::pow(M_V, 2.0));
        }
        inline double kappa4perp(const double & mu) const
        {
            return 1.0 / 10.0 * std::pow(c_rge(mu), 8.0 / 3.0) * a1perp_0 +
                std::pow(c_rge(mu), 8.0 / 3.0) * fpara / fperp_0 * (ms0 - mq0) / (12.0 * M_V) -
                std::pow(c_rge(mu), 8.0) * (std::pow(ms0, 2.0) - std::pow(mq0, 2.0)) / (4.0 * std::pow(M_V, 2.0));
        }
    };

    AntiKStarLCDAs::AntiKStarLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<AntiKStarLCDAs>(new Implementation<AntiKStarLCDAs>(p, o, *this))
    {
    }

    AntiKStarLCDAs::~AntiKStarLCDAs()
    {
    }

    VectorLCDAs *
    AntiKStarLCDAs::make(const Parameters & p, const Options & o)
    {
        return new AntiKStarLCDAs(p, o);
    }

    double
    AntiKStarLCDAs::a1para(const double & mu) const
    {
        return _imp->a1para(mu);
    }

    double
    AntiKStarLCDAs::a2para(const double & mu) const
    {
        return _imp->a2para(mu);
    }

    double
    AntiKStarLCDAs::a3para(const double & mu) const
    {
        return _imp->a3para(mu);
    }

    double
    AntiKStarLCDAs::a4para(const double & mu) const
    {
        return _imp->a4para(mu);
    }

    double
    AntiKStarLCDAs::fpara() const
    {
        return _imp->fpara();
    }

    double
    AntiKStarLCDAs::a1perp(const double & mu) const
    {
        return _imp->a1perp(mu);
    }

    double
    AntiKStarLCDAs::a2perp(const double & mu) const
    {
        return _imp->a2perp(mu);
    }

    double
    AntiKStarLCDAs::a3perp(const double & mu) const
    {
        return _imp->a3perp(mu);
    }

    double
    AntiKStarLCDAs::a4perp(const double & mu) const
    {
        return _imp->a4perp(mu);
    }

    double
    AntiKStarLCDAs::fperp(const double & mu) const
    {
        return _imp->fperp(mu);
    }

    double
    AntiKStarLCDAs::zeta3para(const double & mu) const
    {
        return _imp->zeta3para(mu);
    }

    double
    AntiKStarLCDAs::lambda3paratilde(const double & mu) const
    {
        return _imp->lambda3paratilde(mu);
    }

    double
    AntiKStarLCDAs::omega3paratilde(const double & mu) const
    {
        return _imp->omega3paratilde(mu);
    }

    double
    AntiKStarLCDAs::kappa3para(const double & mu) const
    {
        return _imp->kappa3para(mu);
    }

    double
    AntiKStarLCDAs::omega3para(const double & mu) const
    {
        return _imp->omega3para(mu);
    }

    double
    AntiKStarLCDAs::lambda3para(const double & mu) const
    {
        return _imp->lambda3para(mu);
    }

    double
    AntiKStarLCDAs::kappa3perp(const double & mu) const
    {
        return _imp->kappa3perp(mu);
    }

    double
    AntiKStarLCDAs::omega3perp(const double & mu) const
    {
        return _imp->omega3perp(mu);
    }

    double
    AntiKStarLCDAs::lambda3perp(const double & mu) const
    {
        return _imp->lambda3perp(mu);
    }

    double
    AntiKStarLCDAs::zeta4para(const double & mu) const
    {
        return _imp->zeta4para(mu);
    }

    double
    AntiKStarLCDAs::omega4paratilde(const double & mu) const
    {
        return _imp->omega4paratilde(mu);
    }

    double
    AntiKStarLCDAs::zeta4perp(const double & mu) const
    {
        return _imp->zeta4perp(mu);
    }

    double
    AntiKStarLCDAs::zeta4perptilde(const double & mu) const
    {
        return _imp->zeta4perptilde(mu);
    }

    double
    AntiKStarLCDAs::kappa4para(const double & mu) const
    {
        return _imp->kappa4para(mu);
    }

    double
    AntiKStarLCDAs::kappa4perp(const double & mu) const
    {
        return _imp->kappa4perp(mu);
    }

    double
    AntiKStarLCDAs::phipara(const double & u, const double & mu) const
    {
         // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1(1, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2(2, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_3(3, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_4(4, 3.0 / 2.0);

        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1.evaluate(x);
        const double c2 = gp_2.evaluate(x);
        const double c3 = gp_3.evaluate(x);
        const double c4 = gp_4.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1para(mu) * c1 + _imp->a2para(mu) * c2 + _imp->a3para(mu) * c3 + _imp->a4para(mu) * c4);
    }

    double
    AntiKStarLCDAs::phiperp(const double & u, const double & mu) const
    {
         // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1(1, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2(2, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_3(3, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_4(4, 3.0 / 2.0);

        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1.evaluate(x);
        const double c2 = gp_2.evaluate(x);
        const double c3 = gp_3.evaluate(x);
        const double c4 = gp_4.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1perp(mu) * c1 + _imp->a2perp(mu) * c2 + _imp->a3perp(mu) * c3 + _imp->a4perp(mu) * c4);
    }

    Diagnostics
    AntiKStarLCDAs::diagnostics() const
    {
        Diagnostics results;

        results.add(Diagnostics::Entry{ _imp->c_rge(1.0), "RGE coefficient C(mu = 1.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(2.0), "RGE coefficient C(mu = 2.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(3.0), "RGE coefficient C(mu = 3.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(4.0), "RGE coefficient C(mu = 4.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(5.0), "RGE coefficient C(mu = 5.0 GeV)" });

        return results;
    }

    template <>
    struct Implementation<KStarLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 (even) para Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a1para_0;
        UsedParameter a2para_0;
        UsedParameter a3para_0;
        UsedParameter a4para_0;
        UsedParameter fpara;

        // twist 2 (tensor) Gegenbauer coefficients and normalization at mu = 1 GeV
        UsedParameter a1perp_0;
        UsedParameter a2perp_0;
        UsedParameter a3perp_0;
        UsedParameter a4perp_0;
        UsedParameter fperp_0;

        // twist 3 LCDA parameters at mu = 1 GeV
        UsedParameter zeta3para_0;
        UsedParameter lambda3paratilde_0;
        UsedParameter omega3paratilde_0;
        UsedParameter kappa3para_0;
        UsedParameter omega3para_0;
        UsedParameter lambda3para_0;
        UsedParameter kappa3perp_0;
        UsedParameter omega3perp_0;
        UsedParameter lambda3perp_0;

        // twist 4 LCDA parameters at mu = 1 Gev
        UsedParameter zeta4para_0;
        UsedParameter omega4paratilde_0;
        UsedParameter zeta4perp_0;
        UsedParameter zeta4perptilde_0;

        // Kstar mass
        UsedParameter M_V;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a1para_0(p["K^*::a1para@1GeV"], u),
            a2para_0(p["K^*::a2para@1GeV"], u),
            a3para_0(p["K^*::a3para@1GeV"], u),
            a4para_0(p["K^*::a4para@1GeV"], u),
            fpara(p["K^*::fpara"], u),
            a1perp_0(p["K^*::a1perp@1GeV"], u),
            a2perp_0(p["K^*::a2perp@1GeV"], u),
            a3perp_0(p["K^*::a3perp@1GeV"], u),
            a4perp_0(p["K^*::a4perp@1GeV"], u),
            fperp_0(p["K^*::fperp@1GeV"], u),
            zeta3para_0(p["K^*::zeta3para@1GeV"], u),
            lambda3paratilde_0(p["K^*::lambda3paratilde@1GeV"], u),
            omega3paratilde_0(p["K^*::omega3paratilde@1GeV"], u),
            kappa3para_0(p["K^*::kappa3para@1GeV"], u),
            omega3para_0(p["K^*::omega3para@1GeV"], u),
            lambda3para_0(p["K^*::lambda3para@1GeV"], u),
            kappa3perp_0(p["K^*::kappa3perp@1GeV"], u),
            omega3perp_0(p["K^*::omega3perp@1GeV"], u),
            lambda3perp_0(p["K^*::lambda3perp@1GeV"], u),
            zeta4para_0(p["K^*::zeta4para@1GeV"], u),
            omega4paratilde_0(p["K^*::omega4paratilde@1GeV"], u),
            zeta4perp_0(p["K^*::zeta4perp@1GeV"], u),
            zeta4perptilde_0(p["K^*::zeta4perptilde@1GeV"], u),
            M_V(p["mass::K_u^*"], u),
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

            throw InternalError("Implementation<KStarLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        // running of twist 2 parameters
        inline double a1para(const double & mu) const
        {
            return -1.0 * a1para_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }

        inline double a2para(const double & mu) const
        {
            return a2para_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double a3para(const double & mu) const
        {
            return -1.0 * a3para_0 * std::pow(c_rge(mu), 314.0 / 45.0);
        }

        inline double a4para(const double & mu) const
        {
            return a4para_0 * std::pow(c_rge(mu), 364.0 / 45.0);
        }

        inline double a1perp(const double & mu) const
        {
            return -1.0 * a1perp_0 * std::pow(c_rge(mu), 36.0 / 9.0);
        }

        inline double a2perp(const double & mu) const
        {
            return a2perp_0 * std::pow(c_rge(mu), 52.0 / 9.0);
        }

         inline double a3perp(const double & mu) const
        {
            return -1.0 * a3perp_0 * std::pow(c_rge(mu), 64.0 / 9.0);
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

        // running of twist 3 parameters
        const double ms0 = model->m_s_msbar(1.0);
        const double mq0 = model->m_ud_msbar(1.0) / 2.0;

        inline double zeta3para(const double & mu) const
        {
            return (zeta3para_0 * std::pow(c_rge(mu), 77.0 / 9.0) * fpara +
                (6.0 * a1perp_0 * (-1.0 + std::pow(c_rge(mu), 5.0 / 9.0)) * std::pow(c_rge(mu), 8.0) * fperp_0 * (mq0 - ms0)) / (25.0 * M_V) +
                (2.0 * (std::pow(c_rge(mu), 16.0 / 3.0) - std::pow(c_rge(mu), 77.0 / 9.0)) * fperp_0 * (mq0 + ms0)) / (29.0 * M_V)) / fpara;
        }

        inline double kappa3para(const double & mu) const
        {
            return -1.0 * (std::pow(c_rge(mu), 77.0 / 9.0) * fpara * kappa3para_0 -
                (2.0 * std::pow(c_rge(mu), 16.0 / 3.0) * (-1.0 + std::pow(c_rge(mu), 29.0 / 9.0)) * fperp_0 * (mq0 - ms0)) / (29.0 * M_V) +
                (6.0 * a1perp_0 * (-1.0 + std::pow(c_rge(mu), 5.0 / 9.0)) * std::pow(c_rge(mu), 8.0) * fperp_0 * (mq0 + ms0)) / (25.0 * M_V)) / fpara;
        }

        inline double kappa3perp(const double & mu) const
        {
            return -1.0 * (std::pow(c_rge(mu), 55.0 / 9.0) * fperp_0 * kappa3perp_0 -
                (4.0 * std::pow(c_rge(mu), 4.0) * (-1.0 + std::pow(c_rge(mu), 19.0 / 9.0)) * fpara * (mq0 - ms0)) / (19.0 * M_V) +
                (12.0 * a1para_0 * std::pow(c_rge(mu), 55.0 / 9.0) * (-1.0 + std::pow(c_rge(mu), 13.0 / 9.0)) * fpara * (mq0 + ms0)) / (65.0 * M_V)) / fperp(mu);
        }

        inline double omega3perp(const double & mu) const
        {
            return ((-42.0 * a1para_0 * (-1.0 + std::pow(c_rge(mu), 5.0 / 9.0)) * std::pow(c_rge(mu), 68.0 / 9.0) * fpara * (mq0 - ms0)) / (25.0 * M_V) +
                (12.0 * a2para_0 * std::pow(c_rge(mu), 73.0 / 9.0) * (-1 + std::pow(c_rge(mu), 13.0 / 9.0)) * fpara * (mq0 + ms0)) / (13.0 * M_V) +
                (14.0 * (std::pow(c_rge(mu), 4.0) - std::pow(c_rge(mu), 73.0 / 9.0)) * fpara *(mq0 + ms0)) / (37.0 * M_V) +
                std::pow(c_rge(mu), 73.0 / 9.0) * fperp_0 * omega3perp_0) / fperp(mu);
        }

        inline double lambda3perp(const double & mu) const
        {
            return -1.0 * (std::pow(c_rge(mu), 4.0) * (3.0 * fpara * (mq0 - ms0) +
                68.0 * a2para_0 * std::pow(c_rge(mu), 5.555555555555555) * fpara * (mq0 - ms0) -
                51.0 * a1para_0 * std::pow(c_rge(mu), 3.5555555555555554) * fpara * (mq0 + ms0) +
                std::pow(c_rge(mu), 7.555555555555555) * (255.0 * fperp_0 * lambda3perp_0 * M_V +
                (-3.0 + 51.0 * a1para_0 - 68.0 * a2para_0) * fpara * mq0 +
                (3.0 + 51.0 * a1para_0 + 68.0 * a2para_0) * fpara * ms0))) / (255.0 * fperp(mu) * M_V);
        }

        inline double omega3para(const double & mu) const
        {
            return -3.7780044580452606e-7 * (887490.0 * a1perp_0 * std::pow(c_rge(mu), 8.0 + std::sqrt(865) / 18.0) * fperp_0 * (ms0 - mq0) -
                448070.0 * std::pow(c_rge(mu), (96.0 + std::sqrt(865.0)) / 18.0) * fperp_0 * (ms0 + mq0) -
                89465220.0 * a2perp_0 * std::pow(c_rge(mu), (176.0 + std::sqrt(865)) / 18.0) * fperp_0 * (ms0 + mq0) +
                std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0) * (fperp_0 * ((224035.0 - 6811.0 * std::sqrt(865.0) +
                27.0 * (-16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 -
                306.0 * (-146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * ms0 +
                (224035.0 - 6811.0 * std::sqrt(865.0) -
                27.0 * (-16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 -
                306.0 * (-146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * mq0) +
                765.0 * fpara * M_V * (2.0 * (-865.0 + 26.0 * std::sqrt(865.0)) * omega3para_0 -
                63.0 * std::sqrt(865.0) * omega3paratilde_0)) +
                std::pow(c_rge(mu), 11.38888888888889) * (fperp_0 * ((7.0 * (32005.0 + 973.0 * std::sqrt(865.0)) -
                27.0 * (16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 +
                306.0 * (146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * ms0 +
                (7.0 * (32005.0 + 973.0 * std::sqrt(865.0)) +
                27.0 * (16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 +
                306.0 * (146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * mq0) -
                765.0 * fpara * M_V * (2.0 * (865.0 + 26.0 * std::sqrt(865.0)) * omega3para_0 -
                63.0 * std::sqrt(865.0) * omega3paratilde_0))) / (std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V);
        }

        inline double omega3paratilde(const double & mu) const
        {
            return (9.0 * a1perp_0 * ((865.0 - 7.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889) - 1730.0 * std::pow(c_rge(mu), 8.0 + std::sqrt(865.0) / 18.0) +
                (865.0 + 7.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0)) * fperp_0 * (mq0 - ms0)) / (147050.0 * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V) -
                ((-((9515.0 + 6707.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889)) + (-9515.0 + 6707.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0) +
                19030.0 * std::pow(c_rge(mu), (96.0 + std::sqrt(865.0)) / 18.0)) * fperp_0 * (mq0 + ms0)) / (1.191105e7 * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V) -
                (3.0 * a2perp_0 * (-((2595.0 + 91.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889)) + (-2595.0 + 91.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0) +
                5190.0 * std::pow(c_rge(mu), (176.0 + std::sqrt(865.0)) / 18.0)) * fperp_0 * (mq0 + ms0)) / (4325.0 * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V) +
                (3.0 * (-std::pow(c_rge(mu), 11.38888888888889) + std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0)) * omega3para_0) / (std::sqrt(865.0) * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0)) +
                (((865.0 - 26.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889) + (865.0 + 26.0 * std::sqrt(865.0)) * std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0)) * omega3paratilde_0) /
                (1730.0 * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0));
        }

        inline double lambda3para(const double & mu) const
        {
            return -1.0 * (19030.0 * std::pow(c_rge(mu), (96.0 + std::sqrt(865.0)) / 18.0) * fperp_0 * (ms0 - mq0) +
                42879780.0 * a2perp_0 * std::pow(c_rge(mu), (176 + std::sqrt(865)) / 18.0) * fperp_0 *(ms0 - mq0) -
                1261170.0 * a1perp_0 * std::pow(c_rge(mu), 8.0 + std::sqrt(865.0) / 18.0) * fperp_0 * (ms0 + mq0) +
                std::pow(c_rge(mu), 11.38888888888889) * (fperp_0 * (-((9515.0 + 6707.0 * std::sqrt(865.0) +
                729.0 * (-865.0 + 7.0 * std::sqrt(865.0)) * a1perp_0 + 8262.0 * (2595.0 + 91.0 * std::sqrt(865.0)) * a2perp_0) * ms0) +
                (9515.0 + 6707.0 * std::sqrt(865.0) - 729.0 * (-865.0 + 7.0 * std::sqrt(865.0)) * a1perp_0 +
                8262.0 * (2595.0 + 91.0 * std::sqrt(865.0)) * a2perp_0) * mq0) -
                6885.0 * fpara * ((-865.0 + 26.0 * std::sqrt(865.0)) * lambda3para_0 +
                6.0 * std::sqrt(865.0) * lambda3paratilde_0) * M_V) +
                std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0) * (fperp_0 * ((-9515.0 + 6707.0 * std::sqrt(865.0) + 729.0 * (865.0 + 7.0 * std::sqrt(865.0)) * a1perp_0 +
                8262.0 * (-2595.0 + 91.0 * std::sqrt(865.0)) * a2perp_0) * ms0 +
                (9515.0 - 6707.0 * std::sqrt(865.0) + 729.0 * (865.0 + 7.0 * std::sqrt(865.0)) * a1perp_0 - 8262.0 *(-2595.0 + 91.0 * std::sqrt(865.0)) * a2perp_0) * mq0) +
                6885.0 * fpara * ((865.0 + 26.0 * std::sqrt(865.0)) * lambda3para_0 +
                6.0 *std::sqrt(865.0) * lambda3paratilde_0) * M_V)) / (1.191105e7 * std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V);
        }

        inline double lambda3paratilde(const double & mu) const
        {
            return -1.0 * -3.7780044580452606e-7 * (448070.0 * std::pow(c_rge(mu), (96.0 + std::sqrt(865.0)) / 18.0) * fperp_0 * (ms0 - mq0) +
                89465220.0 * a2perp_0 * std::pow(c_rge(mu), (176.0 + std::sqrt(865.0)) / 18.0) * fperp_0 * (ms0 - mq0) -
                887490.0 * a1perp_0 * std::pow(c_rge(mu), 8.0 + std::sqrt(865.0) / 18.0) * fperp_0 * (ms0 + mq0) +
                std::pow(c_rge(mu), 11.38888888888889 + std::sqrt(865.0) / 9.0) * (fperp_0 * ((7.0 * (-32005.0 + 973.0 * std::sqrt(865.0)) -
                27.0 * (-16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 +
                306.0 * (-146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * ms0 + (224035.0 - 6811.0 * std::sqrt(865.0) - 27.0 *(-16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 -
                306.0 * (-146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * mq0) - 765.0 * fpara * (63.0 * std::sqrt(865.0) * lambda3para_0 +
                2.0 * (865.0 - 26.0 * std::sqrt(865.0)) * lambda3paratilde_0) * M_V) +
                std::pow(c_rge(mu), 11.38888888888889) * (fperp_0 * ((-7.0 * (32005.0 + 973.0 * std::sqrt(865.0)) + 27.0 * (16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 -
                306.0 * (146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * ms0 +
                (7.0 * (32005.0 + 973.0 * std::sqrt(865.0)) + 27.0 * (16435.0 + 683.0 * std::sqrt(865.0)) * a1perp_0 +
                306.0 * (146185.0 + 4961.0 * std::sqrt(865.0)) * a2perp_0) * mq0) +
                765.0 * fpara * (63.0 * std::sqrt(865.0) * lambda3para_0 -
                2.0 * (865.0 + 26.0 * std::sqrt(865.0)) * lambda3paratilde_0) * M_V)) / (std::pow(c_rge(mu), std::sqrt(865.0) / 18.0) * fpara * M_V);
        }

        // running of twist 4 parameters
        inline double zeta4para(const double & mu) const
        {
            return zeta4para_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }
        inline double omega4paratilde(const double & mu) const
        {
            return omega4paratilde_0 * std::pow(c_rge(mu), 10.0);
        }
        // mass corrections for running of zeta4perp, zeta4perptilde unknown
        inline double zeta4perp(const double & mu) const
        {
            return 1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) + std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perp_0 +
                1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) - std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perptilde_0;
        }
        inline double zeta4perptilde(const double & mu) const
        {
            return 1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) - std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perp_0 +
                1.0 / 2.0 * (std::pow(c_rge(mu), 49.0 / 9.0) + std::pow(c_rge(mu), 20.0 / 3.0)) * zeta4perptilde_0;
        }
        inline double kappa4para(const double & mu) const
        {
            return -1.0 * (-3.0 / 20.0 * std::pow(c_rge(mu), 32.0 / 9.0) * a1para_0 -
                std::pow(c_rge(mu), 16.0 / 3.0) * fperp_0 / fpara * (ms0 - mq0) / (4.0 * M_V) +
                std::pow(c_rge(mu), 8.0) * (std::pow(ms0, 2.0) - std::pow(mq0, 2.0)) / (2.0 * std::pow(M_V, 2.0)));
        }
        inline double kappa4perp(const double & mu) const
        {
            return -1.0 * (1.0 / 10.0 * std::pow(c_rge(mu), 8.0 / 3.0) * a1perp_0 +
                std::pow(c_rge(mu), 8.0 / 3.0) * fpara / fperp_0 * (ms0 - mq0) / (12.0 * M_V) -
                std::pow(c_rge(mu), 8.0) * (std::pow(ms0, 2.0) - std::pow(mq0, 2.0)) / (4.0 * std::pow(M_V, 2.0)));
        }
    };

    KStarLCDAs::KStarLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<KStarLCDAs>(new Implementation<KStarLCDAs>(p, o, *this))
    {
    }

    KStarLCDAs::~KStarLCDAs()
    {
    }

    VectorLCDAs *
    KStarLCDAs::make(const Parameters & p, const Options & o)
    {
        return new KStarLCDAs(p, o);
    }

    double
    KStarLCDAs::a1para(const double & mu) const
    {
        return _imp->a1para(mu);
    }

    double
    KStarLCDAs::a2para(const double & mu) const
    {
        return _imp->a2para(mu);
    }

    double
    KStarLCDAs::a3para(const double & mu) const
    {
        return _imp->a3para(mu);
    }

    double
    KStarLCDAs::a4para(const double & mu) const
    {
        return _imp->a4para(mu);
    }

    double
    KStarLCDAs::fpara() const
    {
        return _imp->fpara();
    }

    double
    KStarLCDAs::a1perp(const double & mu) const
    {
        return _imp->a1perp(mu);
    }

    double
    KStarLCDAs::a2perp(const double & mu) const
    {
        return _imp->a2perp(mu);
    }

    double
    KStarLCDAs::a3perp(const double & mu) const
    {
        return _imp->a3perp(mu);
    }

    double
    KStarLCDAs::a4perp(const double & mu) const
    {
        return _imp->a4perp(mu);
    }

    double
    KStarLCDAs::fperp(const double & mu) const
    {
        return _imp->fperp(mu);
    }

    double
    KStarLCDAs::zeta3para(const double & mu) const
    {
        return _imp->zeta3para(mu);
    }

    double
    KStarLCDAs::lambda3paratilde(const double & mu) const
    {
        return _imp->lambda3paratilde(mu);
    }

    double
    KStarLCDAs::omega3paratilde(const double & mu) const
    {
        return _imp->omega3paratilde(mu);
    }

    double
    KStarLCDAs::kappa3para(const double & mu) const
    {
        return _imp->kappa3para(mu);
    }

    double
    KStarLCDAs::omega3para(const double & mu) const
    {
        return _imp->omega3para(mu);
    }

    double
    KStarLCDAs::lambda3para(const double & mu) const
    {
        return _imp->lambda3para(mu);
    }

    double
    KStarLCDAs::kappa3perp(const double & mu) const
    {
        return _imp->kappa3perp(mu);
    }

    double
    KStarLCDAs::omega3perp(const double & mu) const
    {
        return _imp->omega3perp(mu);
    }

    double
    KStarLCDAs::lambda3perp(const double & mu) const
    {
        return _imp->lambda3perp(mu);
    }

    double
    KStarLCDAs::zeta4para(const double & mu) const
    {
        return _imp->zeta4para(mu);
    }

    double
    KStarLCDAs::omega4paratilde(const double & mu) const
    {
        return _imp->omega4paratilde(mu);
    }

    double
    KStarLCDAs::zeta4perp(const double & mu) const
    {
        return _imp->zeta4perp(mu);
    }

    double
    KStarLCDAs::zeta4perptilde(const double & mu) const
    {
        return _imp->zeta4perptilde(mu);
    }

    double
    KStarLCDAs::kappa4para(const double & mu) const
    {
        return _imp->kappa4para(mu);
    }

    double
    KStarLCDAs::kappa4perp(const double & mu) const
    {
        return _imp->kappa4perp(mu);
    }

    double
    KStarLCDAs::phipara(const double & u, const double & mu) const
    {
         // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1(1, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2(2, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_3(3, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_4(4, 3.0 / 2.0);

        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1.evaluate(x);
        const double c2 = gp_2.evaluate(x);
        const double c3 = gp_3.evaluate(x);
        const double c4 = gp_4.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1para(mu) * c1 + _imp->a2para(mu) * c2 + _imp->a3para(mu) * c3 + _imp->a4para(mu) * c4);
    }

    double
    KStarLCDAs::phiperp(const double & u, const double & mu) const
    {
         // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1(1, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2(2, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_3(3, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_4(4, 3.0 / 2.0);

        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1.evaluate(x);
        const double c2 = gp_2.evaluate(x);
        const double c3 = gp_3.evaluate(x);
        const double c4 = gp_4.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1perp(mu) * c1 + _imp->a2perp(mu) * c2 + _imp->a3perp(mu) * c3 + _imp->a4perp(mu) * c4);
    }

    Diagnostics
    KStarLCDAs::diagnostics() const
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
