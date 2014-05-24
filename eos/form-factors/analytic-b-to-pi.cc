/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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

#include <eos/form-factors/analytic-b-to-pi.hh>
#include <eos/form-factors/pi-lcdas.hh>
#include <eos/utils/derivative.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/model.hh>
#include <eos/utils/polylog.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

#include <functional>

#include <gsl/gsl_sf_gamma.h>

namespace eos
{
    template <>
    struct Implementation<AnalyticFormFactorBToPiDKMMO2008>
    {
        std::shared_ptr<Model> model;

        // hadronic parameters
        UsedParameter MB;
        UsedParameter mpi;
        UsedParameter fpi;

        // Borel parameters, thresholds and renormalization scale
        UsedParameter M2;
        UsedParameter Mprime2;
        UsedParameter s0B;
        UsedParameter sprime0B;
        UsedParameter mu;

        // Parameter for the estimation of NNLO corrections
        UsedParameter zeta_nnlo;

        // QCD parameters
        UsedParameter m0;
        UsedParameter cond_GG;
        UsedParameter r_vac;

        PionLCDAs pi;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            MB(p["mass::B_d"], u),
            mpi(p["mass::pi^+"], u),
            fpi(p["decay-constant::pi"], u),
            M2(p["B->pi::M^2@DKMMO2008"], u),
            Mprime2(p["B->pi::Mp^2@DKMMO2008"], u),
            s0B(p["B->pi::s_0^B@DKMMO2008"], u),
            sprime0B(p["B->pi::sp_0^B@DKMMO2008"], u),
            mu(p["B->pi::mu@DKMMO2008"], u),
            zeta_nnlo(p["B->pi::zeta(NNLO)@DKMMO2008"], u),
            m0(p["QCD::m_0"], u),
            cond_GG(p["QCD::cond_GG"], u),
            r_vac(p["QCD::r_vac"], u),
            pi(p, o)
        {
        }

        inline double m_b_msbar(const double & mu) const
        {
            return model->m_b_msbar(mu);
        }

        static double rho_1(const double & s, const double & mb, const double & mu)
        {
            const double mb2 = mb * mb, x = mb2 / s;
            const double lnx = std::log(x), ln1mx = std::log(1.0 - x), re_li2_x = real(dilog(complex<double>(x, 0.0)));
            const double lnmumb = std::log(mu / mb);

            double result = s / 2 * (1.0 - x) * (
                    (1.0 - x) * (4.0 * re_li2_x + 2.0 * lnx * ln1mx - (5.0 - 2.0 * x) * ln1mx)
                  + (1.0 - 2.0 * x) * (3.0 - x) * lnx + 3.0 * (1.0 - 3.0 * x) * 2.0 * lnmumb
                  + (17.0 - 33.0 * x) / 2.0
                );

            return result;
        }

        static double delta_1(const double & mb, const double & mu, const double & Mprime2)
        {
            const double mb2 = mb * mb, mu2 = mu * mu;
            const double gamma = gsl_sf_gamma_inc(0.0, mb2 / Mprime2);

            return -3.0 / 2.0 * (gamma * std::exp(mb2 / Mprime2) - 1.0 - (1.0 - mb2 / Mprime2) * (std::log(mu2 / mb2) + 4.0 / 3.0));
        }

        double decay_constant() const
        {
            static const double pi = M_PI, pi2 = pi * pi;
            static const double eps = 1.0e-10;

            const double MB2 = MB * MB, MB4 = MB2 * MB2;
            const double mb = this->m_b_msbar(mu), mb2 = mb * mb, mb4 = mb2 * mb2;
            const double Mprime4 = Mprime2 * Mprime2;
            const double m02 = m0 * m0;

            const double cond_qq_mu = -fpi * fpi * this->pi.mupi(mu) / 2.0; // <qq>@mu
            const double cond_qq_1 = -fpi * fpi * this->pi.mupi(1.0) / 2.0; // <qq>@1GeV

            const double alpha_s_mu = model->alpha_s(mu());
            const double alpha_s_1 = model->alpha_s(1.0); // alpha_s@1GeV

            using namespace std::placeholders;

            std::function<double (const double &)> integrand(
                [&] (const double & s) -> double
                {
                    return std::exp(-s / Mprime2) * ((s - mb2) * (s - mb2) / s + 4.0 * alpha_s_mu / (3.0 * pi) * rho_1(s, mb, mu));
                }
            );
            const double integral = integrate(integrand, 64, mb2 + eps, sprime0B);

            double result = std::exp(MB2 / Mprime2) / MB4 * (3.0 * mb2 / (8.0 * pi2) * integral
                + mb2 * std::exp(-mb2 / Mprime2) * (
                    - mb * cond_qq_mu * (1.0 + 4.0 * alpha_s_mu / (3.0 * pi) * delta_1(mb, mu, Mprime2))
                    - mb * cond_qq_1  * m02 / (2.0 * Mprime2) * (1.0 - mb2 / (2 * Mprime2))
                    + cond_GG / 12.0
                    - 16.0 * pi * alpha_s_1 * cond_qq_1 * cond_qq_1 * r_vac / (27.0 * Mprime2) * (1.0 - mb2 / (4.0 * Mprime2) - mb4 / (12.0 * Mprime4))
                )
            );

            return std::sqrt(result);
        }

        static double delta_1_M2_deriv(const double & mb, const double & mu, const double & Mprime2)
        {
            const double mb2 = mb * mb, mu2 = mu * mu;
            const double gamma = gsl_sf_gamma_inc(0.0, mb2 / Mprime2);

            return -3.0 / 2.0 * (Mprime2 - mb2 * gamma * std::exp(mb2 / Mprime2) - mb2 * (std::log(mu2 / mb2) + 4.0 / 3.0));
        }

        double MB_svz() const
        {
            static const double pi = M_PI, pi2 = pi * pi;
            static const double eps = 1.0e-10;

            const double mb = this->m_b_msbar(mu), mb2 = mb * mb, mb4 = mb2 * mb2;
            const double Mprime4 = Mprime2 * Mprime2;
            const double m02 = m0 * m0;

            const double cond_qq_mu = -fpi * fpi * this->pi.mupi(mu) / 2.0; // <qq>@mu
            const double cond_qq_1 = -fpi * fpi * this->pi.mupi(1.0) / 2.0; // <qq>@1GeV

            const double alpha_s_mu = model->alpha_s(mu());
            const double alpha_s_1 = model->alpha_s(1.0); // alpha_s@1GeV

            using namespace std::placeholders;

            std::function<double (const double &)> integrand_numerator(
                [&] (const double & s) -> double
                {
                    return std::exp(-s / Mprime2) * ((s - mb2) * (s - mb2) + 4.0 * s * alpha_s_mu / (3.0 * pi) * rho_1(s, mb, mu));
                }
            );
            const double integral_numerator = integrate(integrand_numerator, 64, mb2 + eps, sprime0B);
            std::function<double (const double &)> integrand_denominator(
                [&] (const double & s) -> double
                {
                    return std::exp(-s / Mprime2) * ((s - mb2) * (s - mb2) / s + 4.0 * alpha_s_mu / (3.0 * pi) * rho_1(s, mb, mu));
                }
            );
            const double integral_denominator = integrate(integrand_denominator, 64, mb2 + eps, sprime0B);

            double numerator = 3.0 * mb2 / (8.0 * pi2) * integral_numerator
                + mb4 * std::exp(-mb2 / Mprime2) * (
                    - mb * cond_qq_mu * (1.0 + 4.0 * alpha_s_mu / (3.0 * pi) * delta_1(mb, mu, Mprime2))
                    - mb * cond_qq_1  * m02 / (2.0 * Mprime2) * (1.0 - mb2 / (2 * Mprime2))
                    + cond_GG / 12.0
                    - 16.0 * pi * alpha_s_1 * cond_qq_1 * cond_qq_1 / (27.0 * Mprime2) * (1.0 - mb2 / (4.0 * Mprime2) - mb4 / (12.0 * Mprime4))
                )
                + mb2 * std::exp(-mb2 / Mprime2) * (
                    - mb * cond_qq_mu * 4.0 * alpha_s_mu / (3.0 * pi) * delta_1(mb, mu, Mprime2)
                    - mb * cond_qq_1  * m02 / (2.0 * Mprime2) * (mb2 - Mprime2)
                    + 16.0 * pi * alpha_s_1 * cond_qq_1 * cond_qq_1 / (27.0 * 4.0 * Mprime4) * (4.0 * Mprime4 - 2.0 * Mprime2 * mb2 - mb4)
                );
            double denominator = 3.0 * mb2 / (8.0 * pi2) * integral_denominator
                + mb2 * std::exp(-mb2 / Mprime2) * (
                    - mb * cond_qq_mu * (1.0 + 4.0 * alpha_s_mu / (3.0 * pi) * delta_1(mb, mu, Mprime2))
                    - mb * cond_qq_1  * m02 / (2.0 * Mprime2) * (1.0 - mb2 / (2 * Mprime2))
                    + cond_GG / 12.0
                    - 16.0 * pi * alpha_s_1 * cond_qq_1 * cond_qq_1 / (27.0 * Mprime2) * (1.0 - mb2 / (4.0 * Mprime2) - mb4 / (12.0 * Mprime4))
                );

            return std::sqrt(numerator / denominator);
        }

        double F_lo_tw2_integrand(const double & u, const double & q2, const double _M2) const
        {
            const double mb = this->m_b_msbar(mu), mb2 = mb * mb, mpi2 = mpi * mpi;

            return std::exp(-(mb2 - q2 * (1.0 - u) + mpi2 * u * (1.0 - u)) / (u * _M2)) / u * this->pi.phi(u, mu);
        }

        double F_lo_tw2(const double & q2, const double & _M2) const
        {
            const double mb = this->m_b_msbar(mu), mb2 = mb * mb;
            const double u0 = std::max(1e-10, (mb2 - q2) / (s0B - q2));

            std::function<double (const double &)> integrand(std::bind(&Implementation<AnalyticFormFactorBToPiDKMMO2008>::F_lo_tw2_integrand, this, std::placeholders::_1, q2, _M2));

            return mb2 * fpi * integrate(integrand, 64, u0, 1.000);
        }

        double F_lo_tw3_integrand(const double & u, const double & q2, const double & _M2) const
        {
            const double mb = this->m_b_msbar(mu), mb2 = mb * mb, mpi2 = mpi * mpi;
            const double mupi = pi.mupi(mu);
            const double omega3pi = pi.omega3pi(mu);

            // auxilliary functions and their first derivatives
            auto I3 = [&] (const double & u) -> double
            {
                const double u3 = u * u * u, ubar2 = (1.0 - u) * (1.0 - u);

                return 5.0 / 2.0 * u3 * ubar2 * (12.0 + (7.0 * u - 4) * omega3pi);
            };
            auto I3_d1 = [&] (const double & u) -> double
            {
                const double u2 = u * u, ubar = 1.0 - u;

                return 15.0 * u2 * ubar * (6.0 - 10.0 * u - (2.0 - 8.0 * u + 7.0 * u2) * omega3pi);
            };
            auto I3bar = [&] (const double & u) -> double
            {
                const double u3 = u * u * u, ubar2 = (1.0 - u) * (1.0 - u);

                return 5.0 / 2.0 * u3 * ubar2 * (24.0 * u + 6.0 * u * omega3pi - 3.0 * (omega3pi + 4.0));
            };
            auto I3bar_d1 = [&] (const double & u) -> double
            {
                const double u2 = u * u, u3 = u2 * u;

                return 15.0 / 2.0 * u2 * (12.0 * u3 - 25.0 * u2 + 16.0 * u - 3.0) * (omega3pi + 4.0);
            };

            const double u2 = u * u;
            const double tw3a = pi.phi3p(u, mu)
                + (
                    pi.phi3s(u, mu) / u
                    - (mb2 + q2 - u2 * mpi2) / (2 * (mb2 - q2 + u2 * mpi2)) * pi.phi3s_d1(u, mu)
                    - (2 * u * mpi2 * mb2) / power_of<2>(mb2 - q2 + u2 * mpi2) * pi.phi3s(u, mu)
                ) / 3.0;
            const double tw3b = 2.0 / u * (mb2 - q2 - u2 * mpi2) / (mb2 - q2 + u2 * mpi2)
                * (I3_d1(u) - (2.0 * u * mpi2) / (mb2 - q2 + u2 * mpi2) * I3(u));
            const double tw3c = 3.0 * mpi2 / (mb2 - q2 + u2 * mpi2)
                * (I3bar_d1(u) - (2.0 * u * mpi2) / (mb2 - q2 + u2 * mpi2) * I3bar(u));

            return std::exp(-(mb2 - q2 * (1.0 - u) + mpi2 * u * (1.0 - u)) / (u * _M2))
                * (mupi / mb * tw3a - pi.f3pi(mu) / (mb * fpi) * (tw3b + tw3c));
        }

        double F_lo_tw3(const double & q2, const double & _M2) const
        {
            const double mb = this->m_b_msbar(mu), mb2 = mb * mb;
            const double u0 = std::max(1e-10, (mb2 - q2) / (s0B - q2));

            std::function<double (const double &)> integrand(std::bind(&Implementation<AnalyticFormFactorBToPiDKMMO2008>::F_lo_tw3_integrand, this, std::placeholders::_1, q2, _M2));

            return mb2 * fpi * integrate(integrand, 64, u0, 1.000);
        }

        double F_lo_tw4(const double & q2, const double & _M2) const
        {
            const double mb = this->m_b_msbar(mu), mb2 = mb * mb, mpi2 = mpi * mpi, mpi4 = mpi2 * mpi2;
            const double u0 = std::max(1e-10, (mb2 - q2) / (s0B - q2));
            const double a2pi = pi.a2pi(mu);
            const double deltapipi = pi.deltapipi(mu);
            const double omega4pi = pi.omega4pi(mu);

            // auxilliary functions and their first derivatives
            auto I4 = [&] (const double & u) -> double
            {
                const double u2 = u * u, u3 = u2 * u;
                const double ubar = 1.0 - u;

                return -1.0 / 24.0 * u * ubar * (
                        mpi2 * (54.0 * u3 - 81.0 * u2 + 27.0 * ubar + 27.0 * a2pi * (16.0 * u3 - 29.0 * u2 + 13.0 * u - 1.0))
                        + 16.0 * u * (20.0 * u - 30.0) * deltapipi
                    );
            };
            auto I4_d1 = [&] (const double & u) -> double
            {
                const double u2 = u * u, u3 = u2 * u, u4 = u2 * u2;

                return 1.0 / 24 * (
                        27.0 * mpi2 * (
                            (10.0 * u4 - 20.0 * u3 + 6.0 * u2 + 4.0 * u - 1.0)
                            + a2pi * (80.0 * u4 - 180.0 * u3 + 126.0 * u2 - 28.0 * u + 1)
                        )
                        + 160.0 * u * (6.0 - 15.0 * u + 8.0 * u2) * deltapipi
                    );
            };
            auto I4bar = [&] (const double & u) -> double
            {
                const double u2 = u * u, u3 = u2 * u;
                const double ubar = 1.0 - u;

                return 1.0 / 48.0 * u * ubar * (
                        mpi2 * (
                            -(54.0 * u3 - 81.0 * u2 - 27.0 * u + 27.0)
                            + 27.0 * a2pi * (32.0 * u3 - 43.0 * u2 + 11.0 * u + 1.0)
                        )
                        - 20.0 * u * (
                            (12.0 - 20.0 * u)
                            + (378.0 * u2 - 567.0 * u + 189.0) * omega4pi
                        ) * deltapipi
                    );
            };
            auto I4barI = [&] (const double & u) -> double
            {
                const double u2 = u * u;
                const double ubar = 1.0 - u, ubar2 = ubar * ubar;

                return 1.0 / 96.0 * u2 * ubar2 * (
                        mpi2 * (
                            9.0 * (3.0 + 2.0 * ubar * u)
                            + 9.0 * a2pi * (32.0 * u2 - 26.0 * u - 3.0)
                        )
                        + 40.0 * u * (4.0 + 63.0 * ubar * omega4pi) * deltapipi
                    );
            };
            auto I4bar_d1 = [&] (const double & u) -> double
            {
                const double u2 = u * u, u3 = u2 * u, u4 = u2 * u2;

                return 1.0 / 48.0 * (
                        27.0 * mpi2 * (
                            (10.0 * u4 - 20.0 * u3 + 6.0 * u2 + 4.0 * u - 1.0)
                            - a2pi * (160.0 * u4 - 300.0 * u3 + 162.0 * u2 - 20.0 * u - 1.0)
                        )
                        + 40.0 * u * (
                            (-40.0 * u2 + 48.0 * u - 12.0)
                            + 189.0 * (5.0 * u3 - 10.0 *  u2 + 6.0 * u - 1.0) * omega4pi
                        ) * deltapipi
                    );
            };
            std::function<double (const double &)> integrand(
                [&] (const double & u) -> double
                {
                    const double u2 = u * u;

                    const double tw4psi = u * pi.psi4(u, mu) + (mb2 - q2 - u2 * mpi2) / (mb2 - q2 + u2 * mpi2) * pi.psi4_i(u, mu);
                    const double tw4phi = (
                            pi.phi4_d2(u, mu)
                            - 6.0 * u * mpi2 / (mb2 - q2 + u2 * mpi2) * pi.phi4_d1(u, mu)
                            + 12.0 * u * mpi4 / power_of<2>(mb2 - q2 + u2 * mpi2) * pi.phi4(u, mu)
                        ) * mb2 * u / (4 * (mb2 - q2 + u2 * mpi2));
                    const double tw4I4 = I4_d1(u) - 2.0 * u * mpi2 / (mb2 - q2 + u2 * mpi2) * I4(u);
                    const double tw4I4bar1 = (u * I4bar_d1(u) + (mb2 - q2 - 3.0 * u2 * mpi2) / (mb2 - q2 + u2 * mpi2) * I4bar(u)) * 2.0 * u * mpi2 / (mb2 - q2 + u2 * mpi2);
                    const double tw4I4bar2 = (I4bar(u) + 6.0 * u * mpi2 / (mb2 - q2 + u2 * mpi2) * I4barI(u)) * 2.0 * u * mpi2 * (mb2 - q2 - u2 * mpi2) / (mb2 - q2 + u2 * mpi2);

                    return std::exp(-(mb2 - q2 * (1.0 - u) + mpi2 * u * (1.0 - u)) / (u * _M2)) *
                            (tw4psi - tw4phi - tw4I4 - tw4I4bar1 - tw4I4bar2) / (mb2 - q2 + u2 * mpi2);
                }
            );

            return mb2 * fpi * integrate(integrand, 64, u0, 1 - 1e-10);
        }

        double F_nlo_tw2(const double & q2, const double & _M2) const
        {
            // Reminder: q2 is the kinematic variable associated with the momentum
            // transfer, while s is the kinematic variable in which the function is
            // analytically continued. See also comment at beginning of Appendix B
            // of [DKMMO2008], p. 21.
            const double mb = this->m_b_msbar(mu), mb2 = mb * mb;
            const double a2pi = pi.a2pi(mu), a4pi = pi.a4pi(mu);
            const double r1 = q2 / mb2;

            // imaginary parts of the hard scattering kernel, integrated over rho.
            auto T1tw2theta1mrho = [&] (const double & r1, const double & r2) -> double
            {
                const double r12 = r1 * r1, r13 = r12 * r1, r14 = r12 * r12, r15 = r14 * r1;
                const double r22 = r2 * r2, r23 = r22 * r2, r24 = r22 * r22, r25 = r24 * r2;
                const double L = std::log(power_of<2>(r2 - 1.0) * mb2 / (mu * mu * r2));

                const double ca0 = power_of<4>(r1 - r2) * (-3.0 + r1 + r2 * 2.0);
                const double ca2 = power_of<2>(r1 - r2) * (
                               (-125.0 + r1  * 155.0  - r12 * 43.0 + r13)
                        + r2 * (220.0 - r1 * 224.0 + r12 * 40.0)
                        + r22 * (-108.0 + 72.0 * r1)
                        + r23 * 12.0);
                const double ca4 = (-3087.0 + r1 * 6804.0 - r12 * 5096.0 + r13 * 1484.0 - r14 * 136.0 + r15)
                        + r2  * (8631.0 - 17024.0 * r1 + 10836.0 * r12 - 2424.0 * r13 + 131.0 * r14)
                        + r22 * (-8750.0 + 14700.0 * r1 - 7200.0 * r12 + 950.0 * r13)
                        + r23 * (3850.0 - r1 * 5000.0 + r12 * 1450.0)
                        + r24 * (-675.0 + r1 * 525.0 )
                        + r25 * 30.0;

                const double cb0 = power_of<4>(r1 - r2);
                const double cb2 = power_of<2>(r1 - r2) * (15.0 - r1 * 10.0 + r12 + r2 * (-20.0 + r1 * 8.0) + r22 * 6.0);
                const double cb4 = (210.0 - r1  * 336.0 + r12 * 168.0 - r13 * 28.0 + r14)
                        + r2  * (-504.0 + r1 * 672.0 - r12 * 252.0 + r13 * 24.0)
                        + r22 * (420.0 - r1 * 420.0 + r12 * 90.0)
                        + r23 * (-140.0 + r1 * 80.0)
                        + r24 * 15.0;

                return (
                        (r1 - r2) * (L - 1.0 / r2) * (ca0 + ca2 * a2pi + ca4 * a4pi)
                        + (r1 - 1.0) * (1.0 / r2 - 1.0) * (r2 - r1) * (cb0 + cb2 * a2pi + cb4 * a4pi)
                        + (1.0 - r1) * (r1 - 1.0) * (L - 1.0) * (cb0 + cb2 * a2pi + cb4 * a4pi)
                    ) * (r1 - 1.0) * 3.0 / power_of<8>(r1 - r2);
            };
            auto T1tw2thetarhom1 = [&] (const double & r1, const double & r2) -> double
            {
                const double r12 = r1 * r1, r13 = r12 * r1, r14 = r12 * r12, r15 = r14 * r1, r16 = r13 * r13;
                const double r22 = r2 * r2, r23 = r22 * r2, r24 = r22 * r22, r25 = r24 * r2, r26 = r23 * r23, r27 = r24 * r23,
                      r28 = r24 * r24;
                const double Lr2 = std::log(r2), Lr2m1 = std::log(r2 - 1.0), Lmu = std::log(mb2 / (mu * mu));

                const double ca00 = (-r1 * 4. + r12 * 4.)
                    + r2  * (3. + r1 * 12. - r12 * 12.)
                    + r22 * (-13. - r1 * 4. + r12 * 8.)
                    + r23 * (13. - r1 * 4.)
                    - r24 * 3.;
                const double ca0mu = r2 * (1. - r1 * 3. + r12 * 2.)
                    + r22 * (r1 * 2. - r12 * 2.)
                    + r23 * (-1. + r1);
                const double ca0r2 = r2 * (-1. + r12)
                    + r22 * (3. - r1 * 4. + r12);
                const double ca0r2m1 = 2.0 * ca0mu;

                const double ca20 = (r1 * 1680. - r12 * 3120 + r13 * 1728 - r14 * 288.)
                    + r2 * (-1500. - r1 * 8675. + r12 * 17308. - r13 * 8208. + r14 * 864.)
                    + r22 * (10895. + r1 * 2160. - r12 * 21084. + r13 * 10080. - r14 * 576.)
                    + r23 * (-19396. + r1 * 15264. + r12 * 5412. - r13 * 3600.)
                    + r24 * (12516. - r1 * 12880. + r12 * 1484.)
                    + r25 * (-2576. + r1 * 2451.)
                    + r26 * 61.;
                const double ca2mu = r2 * (-180. + r1 * 1740. - r12 * 2712. + r13 * 1296. - r14 * 144.)
                    + r22 * (-840. - r1 * 1536. + r12 * 4248. - r13 * 2016. + r14 * 144.)
                    + r23 * (2448. - r1 * 1944. - r12 * 1224. + r13 * 720.)
                    + r24 * (-1800. + r1 * 2112. - r12 * 312.)
                    + r25 * (372. - r1 * 372.);
                const double ca2r2 = r2 * (180. + r1 * 840. - r12 * 1728. + r13 * 720. - r14 * 72.)
                    + r22 * (-1740. + r1 * 1536. + r12 * 144. + r13 * 432. - r14 * 72.)
                    + r23 * (1992. - r1 * 2448. + r12 * 1512. - r13 * 576.)
                    + r24 * (-216. - r1 * 672. + r12 * 168.)
                    + r25 * (-300. + r1 * 300.);
                const double ca2r2m1 = 2.0 * ca2mu;

                const double ca40 = r1 * 98910. - r12 * 281610. + r13 * 294000. - r14 * 136500. + r15 * 27000. - r16 * 1800.
                    + r2  * (-92610. - r1 * 628467. + r12 * 2091411. - r13 * 2110325. + r14 * 869950. - r15 * 136800. + r16 * 5400.)
                    + r22 * (865977. - r1 * 51660. - r12 * 3323460. + r13 * 3765400. - r14 * 1417650. + r15 * 181800. - r16 * 3600.)
                    + r23 * (-2201451. + r1 * 2911860. + r12 * 894420. - r13 * 2358600. + r14 * 840450. - r15 * 72000.)
                    + r24 * (2437925. - r1 * 4042510. + r12 * 1372230. + r13 * 345800. - r14 * 156250.)
                    + r25 * (-1293760. + r1 * 2102595. - r12 * 890655. + r13 * 63725.)
                    + r26 * (307725. - r1 * 414708. + r12 * 137664.)
                    + r27 * (-23987. + r1 * 23980)
                    + r28 * 181.;
                const double ca4mu = r2 * (-6300. + r1 * 107730. - r12 * 271530. + r13 * 266700. - r14 * 115950. + r15 * 20250. - r16 * 900.)
                    + r22 * (-63630. - r1 * 103320. + r12 * 557550. - r13 * 603000. + r14 * 246600. - r15 * 35100. + r16 * 900.)
                    + r23 * (242550. - r1 * 299250. - r12 * 210600. + r13 * 411300. - r14 * 158850. + r15 * 14850.)
                    + r24 * (-304500. + r1 * 539400. - r12 * 200700. - r13 * 62400. + r14 * 28200.)
                    + r25 * (169650. - r1 * 304200. + r12 * 147150. - r13 * 12600.)
                    + r26 * (-40950. + r1 * 62820. - r12 * 21870.)
                    + r27 * (3180. - r1 * 3180.);
                const double ca4r2 = r2 * (6300. + r1 * 63630. - r12 * 204750. + r13 * 210000. - r14 * 87750. + r15 * 12600. - r16 * 450.)
                    + r22 * (-107730. + r1 * 103320. + r12 * 166950. - r13 * 237000. + r14 * 74250. + r15 * 3600. - r16 * 450.)
                    + r23 * (233730. - r1 * 425250. + r12 * 210600. - r13 * 45000. + r14 * 65700 - r15 * 10800.)
                    + r24 * (-172200. + r1 * 300600. - r12 * 165600. + r13 * 71400. - r14 * 23700.)
                    + r25 * (34050. - r1 * 16650. - r12 * 54900. + r13 * 8100.)
                    + r26 * (8100. - r1 * 38520. + r12 * 17820.)
                    + r27 * (-2730. + r1 * 2730.);
                const double ca4r2m1 = 2.0 * ca4mu;

                return -3.0 / (r2 * power_of<4>(r1 - r2)) * (ca00 + ca0mu * Lmu + ca0r2 * Lr2 + ca0r2m1 * Lr2m1)
                    + 1.0 / (4.0 * r2 * power_of<6>(r1 - r2)) * (ca20 + ca2mu * Lmu + ca2r2 * Lr2 + ca2r2m1 * Lr2m1) * a2pi
                    + 1.0 / (10.0 * r2 * power_of<8>(r1 - r2)) * (ca40 + ca4mu * Lmu + ca4r2 * Lr2 + ca4r2m1 * Lr2m1) * a4pi;
            };
            auto T1tw2delta = [&] (const double & r1, const double & r2) -> double
            {
                static const double pi = M_PI, pi2 = pi * pi;

                const double r12 = r1 * r1, r13 = r12 * r1, r14 = r12 * r12, r15 = r13 * r12;
                const double r22 = r2 * r2, r23 = r22 * r2, r24 = r22 * r22, r25 = r23 * r22, r26 = r23 * r23;
                const double L1mr1 = std::log(1.0 - r1), Lr2 = std::log(r2), Lr2m1 = std::log(r2 - 1.0), Lmu = std::log(mb2 / (mu * mu));
                const double L1mr12 = L1mr1 * L1mr1, Lr2m12 = Lr2m1 * Lr2m1;
                const double dilogr1 = real(dilog(complex<double>(r1, 0.0)));
                const double dilog1mr2 = real(dilog(complex<double>(1.0 - r2, 0.0)));

                const double ca00 = r2 * (18.0 + pi2 - r1 * (10.0 + pi2)) + r22 * (-10.0 - pi2 + r1 * (2.0 + pi2));
                const double ca0mu = r2 * (-15.0 + r1 * 9.0) + r22 * (9.0 - r1 * 3.0);
                const double ca0r1 = -2.0 + r1 * 2.0 + r2 * (4.0 - r1 * 4.0) + r22 * (-2.0 + r1 * 2.0);
                const double ca0r12 = r2 * (-2.0 + r1 * 2.0) + r22 * (2.0 - r1 * 2.0);

                const double ca20 = r2 * (5.0 * (34.0 + pi2) - r1 * 10.0 * (26.0 + pi2) + r12 * 6.0 * (18.0 + pi2) + r13 * (-10.0 - pi2))
                    + r22 * (-10.0 * (26.0 + pi2) + r1 * 18.0 * (18.0 + pi2) - r12 * 9.0 * (10.0 + pi2) + r13 * (2.0 + pi2))
                    + r23 * (6.0 * (18.0 + pi2) - r1 * 9.0 * (10.0 + pi2) + r12 * 3.0 * (2.0 + pi2))
                    + r24 * (-10.0 - pi2 + r1 * (2.0 + pi2));
                const double ca2mu = r2 * (-135.0 + r1 * 210.0 - r12 * 90.0 + r13 * 9.0)
                    + r22 * (210.0 - r1 * 270.0 + r12 * 81.0 - r13 * 3.0)
                    + r23 * (-90.0 + r1 * 81.0 - r12 * 9.0)
                    + r24 * (9.0 - r1 * 3.0);
                const double ca2r1 = -10.0 + r1 * 20.0 - r12 * 12.0 + r13 * 2.0
                    + r2  * (30.0 - r1 * 56.0 + r12 * 30.0 - r13 * 4.0)
                    + r22 * (-32.0 + r1 * 54.0 - r12 * 24.0 + r13 * 2.0)
                    + r23 * (14.0 - r1 * 20.0 + r12 * 6.0)
                    + r24 * (-2.0 + r1 * 2.0);
                const double ca2r12 = r2 * (-10.0 + r1 * 20.0 - r12 * 12.0 + r13 * 2.0)
                    + r22 * (20.0 - r1 * 36.0 + r12 * 18.0 - r13 * 2.0)
                    + r23 * (-12.0 + r1 * 18.0 - r12 * 6.0)
                    + r24 * (2.0 - r1 * 2.0);

                const double ca40 = r2 * (42.0 * (50.0 + pi2) - r1 * 126.0 * (42.0 + pi2) + r12 * 140.0 * (34.0 + pi2) - r13 * 70.0 * (26.0 + pi2) + r14 * 15.0 * (18.0 + pi2) + r15 * (-10.0 - pi2))
                    + r22 * (-126.0 * (42.0 + pi2) + r1 * 350.0 * (34.0 + pi2) - r12 * 350.0 * (26.0 + pi2) + r13 * 150.0 * (18.0 + pi2) - r14 * 25.0 * (10.0 + pi2) + r15 * (2.0 + pi2))
                    + r23 * (140.0 * (34.0 + pi2) - r1 * 350.0 * (26.0 + pi2) + r12 * 300.0 * (18.0 + pi2) - r13 * 100.0 * (10.0 + pi2) + r14 * 10.0 * (2.0 + pi2))
                    + r24 * (-70.0 * (26.0 + pi2) + r1 * 150.0 * (18.0 + pi2) - r12 * 100.0 * (10.0 + pi2) + r13 * 20.0 * (2.0 + pi2))
                    + r25 * (15.0 * (18.0 + pi2) - r1 * 25.0 * (10.0 + pi2) + r12 * 10.0 * (2.0 + pi2))
                    + r26 * (-10.0 - pi2 + r1 * (2.0 + pi2));
                const double ca4mu = r2 * (-1638.0 + r1 * 4158.0 - r12 * 3780.0 + r13 * 1470.0 - r14 * 225.0 + r15 * 9.0)
                    + r22 * (4158.0 - r1 * 9450.0 + r12 * 7350.0 - r13 * 2250.0 + r14 * 225.0 - r15 * 3.0)
                    + r23 * (-3780.0 + r1 * 7350.0 - r12 * 4500.0 + r13 * 900.0 - r14 * 30.0)
                    + r24 * (1470.0 - r1 * 2250.0 + r12 * 900.0 - r13 * 60.0)
                    + r25 * (-225.0 + r1 * 225.0 - r12 * 30.0)
                    + r26 * (9.0 - r1 * 3.0);
                const double ca4r1 = -84.0 + r1 * 252.0 - r12 * 280.0 + r13 * 140.0 - r14 * 30.0 + r15 * 2.0
                    + r2  * (336.0 - r1 * 952.0 + r12 * 980.0 - r13 * 440.0 + r14 * 80.0 - r15 * 4.0)
                    + r22 * (-532.0 + r1 * 1400.0 - r12 * 1300.0 + r13 * 500.0 - r14 * 70.0 + r15 * 2.0)
                    + r23 * (420.0 - r1 * 1000.0 + r12 * 800.0 - r13 * 240.0 + r14 * 20.0)
                    + r24 * (-170.0 + r1 * 350.0 - r12 * 220.0 + r13 * 40.0)
                    + r25 * (32.0 - r1 * 52.0 + r12 * 20.0)
                    + r26 * (-2.0 + r1 * 2.0);
                const double ca4r12 = r2 * (-84.0 + r1 * 252.0 - r12 * 280.0 + r13 * 140.0 - r14 * 30.0 + r15 * 2.0)
                    + r22 * (252.0 - r1 * 700.0 + r12 * 700.0 - r13 * 300.0 + r14 * 50.0 - r15 * 2.0)
                    + r23 * (-280.0 + r1 * 700.0 - r12 * 600.0 + r13 * 200.0 - r14 * 20.0)
                    + r24 * (140.0 - r1 * 300.0 + r12 * 200.0 - r13 * 40.0)
                    + r25 * (-30.0 + r1 * 50.0 - r12 * 20.0)
                    + r26 * (2.0 - r1 * 2.0);

                return -3.0 / (r2 * power_of<7>(r1 - r2)) * (
                            power_of<4>(r1 - r2) * (ca00 + ca0mu * Lmu + ca0r1 * (L1mr1 - 2.0 * Lr2m1) + ca0r12 * (L1mr12 + Lr2m12 - 2.0 * Lr2 * Lr2m1 + L1mr1 * (Lr2 - 2.0 * Lr2m1) + dilogr1 - 3.0 * dilog1mr2))
                            + 6.0 * power_of<2>(r1 - r2) * (ca20 + ca2mu * Lmu + ca2r1 * (L1mr1 - 2.0 * Lr2m1) + ca2r12 * (L1mr12 + Lr2m12 - 2.0 * Lr2 * Lr2m1 + L1mr1 * (Lr2 - 2.0 * Lr2m1) + dilogr1 - 3.0 * dilog1mr2)) * a2pi
                            + 15.0 * (ca40 + ca4mu * Lmu + ca4r1 * (L1mr1 - 2.0 * Lr2m1) + ca4r12 * (L1mr12 + Lr2m12 - 2.0 * Lr2 * Lr2m1 + L1mr1 * (Lr2 - 2.0 * Lr2m1) + dilogr1 - 3.0 * dilog1mr2)) * a4pi
                        );
            };
            std::function<double (const double &)> integrand(
                [&] (const double & r2) -> double
                {
                    return -2.0 * (T1tw2thetarhom1(r1, r2) + T1tw2theta1mrho(r1, r2) + T1tw2delta(r1, r2))
                        * std::exp(-mb2 * r2 / _M2);
                }
            );

            static const double eps = 1e-12;

            return mb2 * fpi * integrate(integrand, 64, 1.0 + eps, s0B / mb2);
        }

        double F_nlo_tw3(const double & q2, const double _M2) const
        {
            // Reminder: q2 is the kinematic variable associated with the momentum
            // transfer, while s is the kinematic variable in which the function is
            // analytically continued. See also comment at beginning of Appendix B
            // of [DKMMO2008], p. 21.

            static const double pi2 = M_PI * M_PI;

            const double mb = this->m_b_msbar(mu), mb2 = mb * mb;
            const double r1 = q2 / mb2;
            const double lmu = 2.0 * std::log(mb / mu());

            const double mupi = pi.mupi(mu);

            auto T1tw3ptheta1mrho = [&] (const double & r1, const double & r2) -> double
            {
                const double l1 = std::log((r2 - r1) / (r2 - 1.0)), l2 = lmu + std::log((r2 - 1.0) * (r2 - 1.0) / r2);

                return (r1 - r2 * (1.0 + r1 + r2) * l2) * l1 / (r2 * (r1 - r2));
            };
            auto T1tw3pthetarhom1 = [&] (const double & r1, const double & r2) -> double
            {
                const double logr2 = std::log(r2);
                const double l1 = std::log((1.0 - r1) / (r2 - r1));
                const double dl1 = pi2 / 6.0 + std::real(dilog(1.0 / r2)) + logr2 * (logr2 - std::log(r2 - 1.0));
                const double dl2 = std::real(-dilog(r1 / r2) + dilog(r1) - 2.0 * dilog((r2 - 1.0) / (r1 - 1.0)))
                    - logr2 * logr2 / 2.0 + logr2 * std::log(r2 - r1) - 2.0 * std::log((r2 - r1) / (1.0 - r1)) * std::log(r2 - 1.0);

                return (
                        dl1 * (1.0 + r1 + r2) + dl2 * (4.0 * r1 - 1.0)
                        + ((r1 + r2) * (r2 - 1.0) + (r1 * (2.0 - 3.0 * r2) + r2) * logr2) / (2.0 * r2)
                        + l1 * (1.0 - 2.0 * r1 + lmu * (4.0 * r1 - 1.0))
                    ) / (r2 - r1);
            };

            auto T1tw3pdeltarhom1 = [&] (const double & r1, const double & r2) -> double
            {
                const double l1mr1 = std::log(1.0 - r1);
                const double lr2 = std::log(r2), lr2m1 = std::log(r2 - 1.0);
                const double dlr1 = std::real(dilog(complex<double>(r1, 0.0)));
                const double dl1mr2 = std::real(dilog(complex<double>(1.0 - r2, 0.0)));

                return (
                        6.0 - 2.0 * r1 - pi2 / 6.0 * (1.0 + 4.0 * r1)
                        + lr2 * (l1mr1 * r1 - lr2m1 * 2.0 * r1)
                        + lr2m1 * (lr2m1 * (1.0 + 2.0 * r1) - 4.0 + 2.0 * r1 * (r2 - 1.0) / r2 - l1mr1 * 2.0 * r1 + lmu * (1.0 + r1))
                        + lmu * 3.0 / 2.0 * (r1 - 3.0)
                        + l1mr1 * (-l1mr1 + 2.0 + r1 + r1 / r2 - (1.0 + r1) * lmu)

                        - dlr1 + (1.0 - 2.0 * r1) * dl1mr2
                    ) / (r2 - r1);
            };
            auto T1tw3sigmatheta1mrho = [&] (const double & r1, const double & r2) -> double
            {
                const double lr2 = std::log(r2), lr2m1 = std::log(r2 - 1.0);
                const double lr2mr1 = std::log(r2 - r1);

                return (
                        - 6.0 * (r1 * r1 + 2.0 * (r2 - 1.0) * r2 + r1 * (-1.0 + 2.0 * r2 - 2.0 * r2 * r2))
                            / (r2 * (r1 - r2) * (r1 - r2))
                        + lr2mr1 * ((lmu - lr2 + 2.0 * lr2m1) * 6.0 * (1.0 + r1 + r2) / (r1 - r2) - 6.0 * r1 / (r2 * (r1 - r2)))
                        + lr2m1 * ((-2.0 * lr2m1 - lmu + lr2) * 6.0 * (1.0 + r1 + r2) / (r1 - r2)
                            + 6.0 * (-2.0 * (r2 - 1.0) * r2 + r1 * r2 * (2.0 * r2 - 5.0) + r1 * r1 * (1.0 + 2.0 * r2))
                                / ((r2 - r1) * (r2 - r1) * r2)
                        )
                        + (lmu - lr2) * 6.0 * (r1 - 1.0) * (-1.0 + r1 + r2) / ((r2 - r1) * (r2 - r1))
                    ) / (r2 - r1);
            };
            auto T1tw3sigmathetarhom1 = [&] (const double & r1, const double & r2) -> double
            {
                const double l1mr1 = std::log(1.0 - r1);
                const double lr2 = std::log(r2), lr2m1 = std::log(r2 - 1.0);
                const double lr2mr1 = std::log(r2 - r1);
                const double l1 = 2.0 * lr2m1 + lmu - lr2;
                const double dl1 = std::real(dilog(r1) - dilog(r1 / r2) - 2.0 * dilog((r2 - 1.0) / (r1 - 1.0)));
                const double dl2 = std::real(dilog(1.0 / r2)) - l1 * l1;

                return 3.0 * (
                        - dl1 * 2.0 * (4.0 * r1 - 1.0) * (r1 - r2) * r2
                        - dl2 * 2.0 * (r1 - r2) * r2 * (1.0 + r1 + r2)
                        + l1 * (
                            - l1 * (r1 - r2) * r2 * (5.0 + 4.0 * r2)
                            + lr2mr1 * 2.0 * (4.0 * r1 - 1.0) * (r1 - r2) * r2
                            - lr2m1 * 2.0 * (-5.0 + 5.0 * r1 - 3.0 * r2) * (r1 - r2) * r2
                            - lmu * 2.0 * (-3.0 + 2.0 * r1 - 2.0 * r2) * (r1 - r2) * r2
                            + r1 * (r2 - 1.0) * r2 - 5.0 * r2 * r2 + r1 * r1 * (2.0 + r2 - 2.0 * r2 * r2)
                        )
                        + lr2mr1 * (
                            - 2.0 * (-1.0 + 2.0 * r1) * (r1 - r2) * r2
                        )
                        + lr2m1 * (
                            lr2m1 * 4.0 * (r1 - r2) * (-2.0 + 3.0 * r1 - r2) * r2
                            - l1mr1 * 4.0 * (4.0 * r1 - 1.0) * (r1 - r2) * r2
                            + lmu * 2.0 * (-5.0 + 5.0 * r1 - 3.0 * r2) * (r1 - r2) * r2
                            - 2.0 * r1 * (-1.0 + r2) * r2 + 2.0 * r2 * (2.0 + 3.0 * r2) + r1 * r1 * (-4.0 - 2.0 * r2 + 4.0 * r2 * r2)
                        )
                        + l1mr1 * (
                            - lmu * 2.0 * (4.0 * r1 - 1.0) * (r1 - r2) * r2
                            + 2.0 * (-1.0 + 2.0 * r1) * (r1 - r2) * r2
                        )
                        + lmu * (
                            lmu * (-3.0 + 2.0 * r1 - 2.0 * r2) * (r1 - r2) * r2
                            -r1 * (r2 - 1.0) * r2 + r2 * (2.0 + 3.0 * r2) + r1 * r1 * (-2.0 + r2 * (-1.0 + 2.0 * r2))
                        )
                        + (
                            r2 * r2 * (pi2 - 3.0 + (3.0 + pi2) * r2)
                            + r1 * (6.0 - (6.0 + pi2) * r2)
                            - r1 * r1 * (3.0 + r2 * (pi2 - 9.0 + 6.0 * r2))
                        ) / 3.0
                    ) / (power_of<3>(r1 - r2) * r2);

            };
            auto T1tw3sigmadeltarhom1 = [&] (const double & r1, const double & r2) -> double
            {
                const double l1mr1 = std::log(1.0 - r1);
                const double lr2 = std::log(r2), lr2m1 = std::log(r2 - 1.0);
                const double l1 = 2.0 * lr2m1 + lmu - lr2;
                const double l2 = l1mr1 - 2.0 * lr2m1;
                const double dl1 = std::real(dilog(r1)) + l1mr1 * (l1mr1 + lmu);
                const double dl2 = std::real(dilog(1.0 - r2)) + lr2m1 * lr2m1;

                return (
                        dl1 * 6.0 * (r1 * (3.0 - 4.0 * r2) + r2)
                        + dl2 * (-30.0 * r2 + 6.0 * r1 * (-7.0 + 2.0 * r1 + 10.0 * r2))
                        + l1 * l2 * (-12.0 * r2 + 6.0 * r1 * (-2.0 + r1 + 3.0 * r2))
                        + lr2m1 * (
                            lmu * (-18.0 * r2 + 6.0 * r1 * (-5.0 + r1 + 7.0 * r2))
                            - 12.0 * (r2 + r1 * (2.0 - r1 - 3.0 * r2 + r2 * r2)) / r2
                        )
                        - l1mr1 * 6.0 * ((-2.0 + r1) * r1 - 2.0 * r2 + r1 * (5.0 + r1) * r2 + (2.0 - 5.0 * r1) * r2 * r2) / r2
                        + lmu * (-3.0 * r1 * (-17.0 + r1 - 5.0 * r2) + 9.0 * r2)
                        + r1 * (-72.0 + pi2 * (-5.0 + 4.0 * r1)) + r2 * (6.0 * (-1.0 + r1) * r1 + pi2 * (-7.0 + 8.0 * r1))
                        - 6.0 * (1.0 + 3.0 * r2)
                    ) / ((r1 - r2) * (r1 - r2) * (r1 - r2));
            };
            std::function<double (const double &)> integrand(
                [&] (const double & r2) -> double
                {
                    return (
                            2.0 / (r2 - r1) * (T1tw3pthetarhom1(r1, r2) + T1tw3ptheta1mrho(r1, r2) + T1tw3pdeltarhom1(r1, r2))
                            + 1.0 / 3.0 * (T1tw3sigmathetarhom1(r1, r2) + T1tw3sigmatheta1mrho(r1, r2) + T1tw3sigmadeltarhom1(r1, r2))
                        ) * std::exp(-mb2 * r2 / _M2);
                }
            );

// for diagnostics purposes only
#if 0
            std::function<double (const double &)> iT1tw3ptheta1mrho(
                [&] (const double & r2) -> double
                { return T1tw3ptheta1mrho(0.0, r2) * std::exp(-mb2 * r2 / _M2); }
            );
            std::function<double (const double &)> iT1tw3pthetarhom1(
                [&] (const double & r2) -> double
                { return T1tw3pthetarhom1(0.0, r2) * std::exp(-mb2 * r2 / _M2); }
            );
            std::function<double (const double &)> iT1tw3pdeltarhom1(
                [&] (const double & r2) -> double
                { return T1tw3pdeltarhom1(0.0, r2) * std::exp(-mb2 * r2 / _M2); }
            );
            std::function<double (const double &)> iT1tw3sigmatheta1mrho(
                [&] (const double & r2) -> double
                { return T1tw3sigmatheta1mrho(0.0, r2) * std::exp(-mb2 * r2 / _M2); }
            );
            std::function<double (const double &)> iT1tw3sigmathetarhom1(
                [&] (const double & r2) -> double
                { return T1tw3sigmathetarhom1(0.0, r2) * std::exp(-mb2 * r2 / _M2); }
            );
            std::function<double (const double &)> iT1tw3sigmadeltarhom1(
                [&] (const double & r2) -> double
                { return T1tw3sigmadeltarhom1(0.0, r2) * std::exp(-mb2 * r2 / _M2); }
            );
#endif

            static const double eps = 1e-12;

            return fpi * mupi * mb * (
                    integrate(integrand, 64, 1.0 + eps, s0B / mb2)
                    - (
                        2.0 / (1.0 - r1) * (4.0 - 3.0 * lmu)
                        + 2.0 * (1.0 + r1) / power_of<2>(1.0 - r1) * (4.0 - 3.0 * lmu)
                    ) * std::exp(-mb2 / _M2)
                );
        }

        double rescale_factor(const double & q2) const
        {
            const double mb = this->m_b_msbar(mu), mb2 = mb * mb;
            const double u0_q2 = std::max(1e-10, (mb2 - q2) / (s0B - q2));
            const double u0_zero = std::max(1e-10, mb2 / s0B);

            std::function<double (const double &)> integrand_numerator_q2(
                [&] (const double & u) -> double
                {
                    return u * (F_lo_tw2_integrand(u, q2, this->M2()) + F_lo_tw3_integrand(u, q2, this->M2));
                }
            );
            std::function<double (const double &)> integrand_denominator_q2(
                [&] (const double & u) -> double
                {
                    return (F_lo_tw2_integrand(u, q2, this->M2()) + F_lo_tw3_integrand(u, q2, this->M2()));
                }
            );
            std::function<double (const double &)> integrand_numerator_zero(
                [&] (const double & u) -> double
                {
                    return u * (F_lo_tw2_integrand(u, 0.0, this->M2()) + F_lo_tw3_integrand(u, 0.0, this->M2()));
                }
            );
            std::function<double (const double &)> integrand_denominator_zero(
                [&] (const double & u) -> double
                {
                    return (F_lo_tw2_integrand(u, 0.0, this->M2()) + F_lo_tw3_integrand(u, 0.0, this->M2()));
                }
            );

            double result = integrate(integrand_numerator_zero, 64, u0_zero, 1.000) / integrate(integrand_numerator_q2, 64, u0_q2, 1.000)
                / integrate(integrand_denominator_zero, 64, u0_zero, 1.000) * integrate(integrand_denominator_q2, 64, u0_q2, 1.000);

            return result;
        }

        double MB_lcsr(const double & q2) const
        {
            const double M2_rescaled = this->M2() * this->rescale_factor(q2);
            const double alpha_s = model->alpha_s(mu);

            std::function<double (const double &)> F(
                [&] (const double & _M2) -> double
                {
                    const double F_lo  = F_lo_tw2(q2, _M2) + F_lo_tw3(q2, _M2) + F_lo_tw4(q2, _M2);
                    const double F_nlo = F_nlo_tw2(q2, _M2) + F_lo_tw3(q2, _M2);

                    return F_lo + alpha_s / (3.0 * M_PI) * F_nlo;
                }
            );


            double MB2 = M2_rescaled * M2_rescaled * derivative<1, deriv::TwoSided>(F, M2_rescaled) / F(M2_rescaled);

            if (MB2 < 0.0)
                return 0.0;

            return std::sqrt(MB2);
        }

        double f_p(const double & q2) const
        {
            const double MB2 = MB * MB;
            const double M2_rescaled = this->M2() * this->rescale_factor(q2);
            const double fB = decay_constant();
            const double F_lo = F_lo_tw2(q2, M2_rescaled) + F_lo_tw3(q2, M2_rescaled) + F_lo_tw4(q2, M2_rescaled);
            const double F_nlo = F_nlo_tw2(q2, M2_rescaled) + F_nlo_tw3(q2, M2_rescaled);
            /*
             * we estimate the NNLO corrections to obey the relation |F_nnlo / F_nlo| = |F_nlo / F_lo|.
             * Therefore we set F_nnlo = F_nlo^2 / F_lo * zeta_nnlo, where zeta ranges between -1 and +1.
             */
            const double F_nnlo = F_nlo * F_nlo / F_lo * zeta_nnlo;
            const double alpha_s = model->alpha_s(mu);

            return std::exp(MB2 / M2_rescaled) / (2.0 * MB2 * fB) * (F_lo + alpha_s / (3.0 * M_PI) * F_nlo + alpha_s * alpha_s / (9.0 * M_PI * M_PI) * F_nnlo);
        }

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            // Function rho_1, cf. [DKMMO2008], eq. (C.2)
            {
                results.add(Diagnostics::Entry{ rho_1(19.60, 4.16, 4.16), "rho_1(s = 19.60, m_b = 4.16, mu = 4.16), [DKMMO2008]" });
                results.add(Diagnostics::Entry{ rho_1(22.05, 4.16, 4.16), "rho_1(s = 22.05, m_b = 4.16, mu = 4.16), [DKMMO2008]" });
                results.add(Diagnostics::Entry{ rho_1(25.20, 4.16, 4.16), "rho_1(s = 25.20, m_b = 4.16, mu = 4.16), [DKMMO2008]" });
            }

            results.add(Diagnostics::Entry{ this->decay_constant(), "f_B, [DKMM02008]" });

            results.add(Diagnostics::Entry{ this->rescale_factor( 0.0), "rescale_factor(s =  0.0), [DKMMO2008]" });
            results.add(Diagnostics::Entry{ this->rescale_factor(10.0), "rescale_factor(s = 10.0), [DKMMO2008]" });

            return results;
        }
    };

    AnalyticFormFactorBToPiDKMMO2008::AnalyticFormFactorBToPiDKMMO2008(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<AnalyticFormFactorBToPiDKMMO2008>(new Implementation<AnalyticFormFactorBToPiDKMMO2008>(p, o, *this))
    {
    }

    AnalyticFormFactorBToPiDKMMO2008::~AnalyticFormFactorBToPiDKMMO2008()
    {
    }

    FormFactors<PToP> *
    AnalyticFormFactorBToPiDKMMO2008::make(const Parameters & p, unsigned)
    {
        return new AnalyticFormFactorBToPiDKMMO2008(p, Options{ });
    }

    double
    AnalyticFormFactorBToPiDKMMO2008::F_lo_tw2(const double & q2) const
    {
        return _imp->F_lo_tw2(q2, _imp->M2() * _imp->rescale_factor(q2));
    }

    double
    AnalyticFormFactorBToPiDKMMO2008::F_lo_tw3(const double & q2) const
    {
        return _imp->F_lo_tw3(q2, _imp->M2() * _imp->rescale_factor(q2));
    }

    double
    AnalyticFormFactorBToPiDKMMO2008::F_lo_tw4(const double & q2) const
    {
        return _imp->F_lo_tw4(q2, _imp->M2() * _imp->rescale_factor(q2));
    }

    double
    AnalyticFormFactorBToPiDKMMO2008::F_nlo_tw2(const double & q2) const
    {
        return _imp->F_nlo_tw2(q2, _imp->M2() * _imp->rescale_factor(q2));
    }

    double
    AnalyticFormFactorBToPiDKMMO2008::F_nlo_tw3(const double & q2) const
    {
        return _imp->F_nlo_tw3(q2, _imp->M2() * _imp->rescale_factor(q2));
    }

    double
    AnalyticFormFactorBToPiDKMMO2008::f_p(const double & q2) const
    {
        return _imp->f_p(q2);
    }

    double
    AnalyticFormFactorBToPiDKMMO2008::f_0(const double &) const
    {
        throw InternalError("AnalyticFormFactorBToPiDKMMO2008::f_0: Evaluation of time-like form factor not yet implemented");
    }

    double
    AnalyticFormFactorBToPiDKMMO2008::f_t(const double &) const
    {
        throw InternalError("AnalyticFormFactorBToPiDKMMO2008::f_t: Evaluation of tensor form factor not yet implemented");
    }

    double
    AnalyticFormFactorBToPiDKMMO2008::MB_lcsr(const double & q2) const
    {
        return _imp->MB_lcsr(q2);
    }

    double
    AnalyticFormFactorBToPiDKMMO2008::MB_svz() const
    {
        return _imp->MB_svz();
    }

    Diagnostics
    AnalyticFormFactorBToPiDKMMO2008::diagnostics() const
    {
        return _imp->diagnostics();
    }
}

