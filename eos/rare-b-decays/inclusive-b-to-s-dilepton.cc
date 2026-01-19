/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
 *
 * Credit goes to Christoph Bobeth for proofreading and
 * extensive checks.
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
#include <eos/models/model.hh>
#include <eos/rare-b-decays/bremsstrahlung.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/rare-b-decays/em-contributions.hh>
#include <eos/rare-b-decays/inclusive-b-to-s-dilepton.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

#include <cmath>
#include <functional>
#include <utility>
#include <map>
#include <vector>

#include <gsl/gsl_sf_dilog.h>

namespace eos
{
    using namespace std::literals::string_literals;

    /* HLMW2005 */

    template <>
    struct Implementation<BToXsDilepton<HLMW2005>>
    {
        std::shared_ptr<Model> model;

        LeptonFlavorOption opt_l;

        UsedParameter gfermi;

        UsedParameter hbar;

        UsedParameter tau_B;

        UsedParameter m_b_MSbar;

        UsedParameter m_c_MSbar;

        UsedParameter m_tau;

        UsedParameter m_l;

        UsedParameter m_Z;

        UsedParameter mu2_g;

        UsedParameter mu2_pi;

        UsedParameter mu;

        UsedParameter alpha_e;

        cubature::Config cub_conf;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            opt_l(o, options, "l"_ok),
            gfermi(p["WET::G_Fermi"], u),
            hbar(p["QM::hbar"], u),
            tau_B(p["life_time::B" + (destringify<bool>(o.get("admixture"_ok, "true")) ? ("@Y(4S)") : ("_" + o.get("q"_ok, "d")))], u),
            m_b_MSbar(p["mass::b(MSbar)"], u),
            m_c_MSbar(p["mass::c"], u),
            m_tau(p["mass::tau"], u),
            m_l(p["mass::" + opt_l.str()], u),
            m_Z(p["mass::Z"], u),
            mu2_g(p["B->B::mu_G^2@1GeV"], u),
            mu2_pi(p["B->B::mu_pi^2@1GeV"], u),
            mu(p["sb" + opt_l.str() + opt_l.str() + "::mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            cub_conf(cubature::Config().epsrel(1e-4))
        {
            Context ctx("When constructing B->X_sll observables");

            u.uses(*model);
        }

        double m_b_pole() const
        {
            return model->m_b_pole();
        }

        double m_b_msbar() const
        {
            return model->m_b_msbar(mu());
        }

        double m_c_pole() const
        {
            return model->m_c_pole();
        }

        double s_hat(const double & s) const
        {
            double m_b = m_b_pole();

            return s / m_b / m_b;
        }

        /* NLO functions */
        // cf. [BMU:1999A], Eq. (34), p. 9 and [HLMW:2005A], Eq. (127), p. 29
        double omega1_99(const double & s_hat) const
        {
            double li2 = gsl_sf_dilog(s_hat);
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);
            double s_hat2 = s_hat * s_hat;

            return -4.0/3.0 * li2 - 2.0/3.0 * ln1 * ln - 2.0/9.0 * M_PI * M_PI
                - (5.0 + 4.0 * s_hat) / (3.0 * (1.0 + 2.0 * s_hat)) * ln1
                - (2.0 * s_hat * (1.0 + s_hat) * (1.0 - 2.0 * s_hat)) / (3.0 * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat)) * ln
                + (5.0 + 9.0 * s_hat - 6.0 * s_hat2) / (6.0 * (1.0 - s_hat) * (1.0 + 2.0 * s_hat));
        }

        // cf. [HLMW:2005A], Eq. (128), p. 29
        // only valid for 0 < s_hat < 0.4
        double omega2_99(const double & s_hat) const
        {
            double ln = log(s_hat);;
            double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;

            return (-19.2 + 6.1 * s_hat + (37.9 + 17.2 * ln) * s_hat2 - 18.7 * s_hat3) / ((1.0 - s_hat) * (1.0 - s_hat) * (1.0 + 2.0 * s_hat));
        }

        // cf. [HLMW:2005A], Eq. (130), p. 29
        double omega1_77(const double & s_hat) const
        {
            double li2 = gsl_sf_dilog(s_hat);
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);
            double s_hat2 = s_hat * s_hat;

            return -4.0/3.0 * li2 - 2.0/3.0 * ln1 * ln - 2.0/9.0 * M_PI * M_PI
                - (8.0 + s_hat) / (3.0 * (2.0 + s_hat)) * ln1
                - (2.0 * s_hat * (2.0 - 2.0 * s_hat - s_hat2)) / (3.0 * power_of<2>(1.0 - s_hat) * (2.0 + s_hat)) * ln
                - (16.0 - 11.0 * s_hat - 17.0 * s_hat2) / (18.0 * (1.0 - s_hat) * (2.0 + s_hat))
                // We use mu_b in MSbar scheme globally, so use m_b_MSbar here instead of m_b_pole
                - 8.0/3.0 * log(mu / m_b_MSbar());
        }

        // cf.[HLMW:2005A], Eq. (131), p. 29
        double omega1_79(const double & s_hat) const
        {
            double li2 = gsl_sf_dilog(s_hat);
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);

            return -4.0/3.0 * li2 - 2.0/3.0 * ln1 * ln - 2.0/9.0 * M_PI * M_PI
                - (2.0 + 7.0 * s_hat) / (9.0 * s_hat) * ln1
                - (2.0 * s_hat * (3.0 - 2.0 * s_hat)) / (9.0 * power_of<2>(1.0 - s_hat)) * ln
                + (5.0 - 9.0 * s_hat) / (18.0 * (1.0 - s_hat))
                // We use mu_b in MSbar scheme globally, so use m_b_MSbar here instead of m_b_pole
                - 4.0/3.0 * log(mu / m_b_MSbar());
        }

        // cf. [HLMW:2005A], Eq. (126), p. 28
        complex<double> g(const double & y) const
        {
            complex<double> x;

            if (y < 1)
                x = complex<double>(log(1.0 + sqrt(1.0 - y) / (1.0 - sqrt(1.0 - y))), -M_PI);
            else
                x = 2.0 * atan(1.0 / sqrt(y - 1.0));

            return 20.0 / 27.0 + 4.0 / 9.0 * y - 2.0 / 9.0 * (2.0 + y) * sqrt(abs(y - 1)) * x;
        }

        // cf. [HLMW:2005A], Eq. (72), p. 17
        // i = { 1 .. 10, Q3 .. Q6, b }
        complex<double> f(unsigned i, const double & s_hat) const
        {
            if ((7 == i) || (8 == i))
                throw InternalError("[HLMW:2005A] f_i not defined for i = 7,8!");

            static const std::vector<double> rho_b = {
                /* 1 .. 6 */
                0.0, 0.0, -7.0 / 2.0, -2.0 / 3.0, -38.0, -32.0 / 3.0,
                /* 7 .. 10 */
                0.0, 0.0, 0.0, 0.0,
                /* Q3 .. Q6 */
                7.0 / 6.0, 2.0 / 9.0, 38.0 / 3.0, 32.0 / 9.0, -2.0
            };

            static const std::vector<double> rho_c = {
                /* 1 .. 6 */
                4.0 / 3.0, 1.0, 6.0, 0.0, 60.0, 0.0,
                /* 7 .. 10 */
                0.0, 0.0, 0.0, 0.0,
                /* Q3 .. Q6 */
                4.0, 0.0, 40.0, 0.0, 0.0
            };

            static const std::vector<double> rho_0 = {
                /* 1 .. 6 */
                0.0, 0.0, 2.0/9.0, 8.0 / 27.0, 32.0 / 9.0, 128.0 / 27.0,
                /* 7 .. 10 */
                0.0, 0.0, 0.0, 0.0,
                /* Q3 .. Q6 */
                -74.0 / 27.0, -8.0 / 81.0, -752.0 / 27.0, -128.0 / 81.0, 0.0
            };

            static const std::vector<double> rho_sharp = {
                /* 1 .. 6 */
                -16.0 / 27.0, -4.0 / 9.0, 2.0 / 27.0, 8.0 / 81.0, -136.0 / 27.0, 320.0 / 81.0,
                /* 7 .. 10 */
                0.0, 0.0, 0.0, 0.0,
                /* Q3 .. Q6 */
                358.0 / 81.0, -8.0 / 243.0, 1144.0 / 81.0, -320.0 / 243.0, 26.0 / 27.0
            };

            static const std::vector<double> gamma9 = {
                /* 1 .. 6 */
                -32.0 / 27.0, -8.0 / 9.0, -16.0 / 9.0, 32.0 / 27.0, -112.0 / 9.0, 512.0 / 27.0,
                /* 7 .. 10 */
                0.0, 0.0, 8.0, -4.0,
                /* Q3 .. Q6 */
                -272.0 / 27.0, -32.0 / 81.0, 2768.0 / 27.0, -512.0 / 81.0, 16.0 / 9.0
            };

            double m_b = m_b_pole(), m_c = m_c_pole();
            double s = s_hat * power_of<2>(m_b);

            complex<double> g_b = g(4.0 / s_hat);
            complex<double> g_c = g(4.0 * power_of<2>(m_c) / s);

            /* mu == mu in MSbar scheme, so use m_b_MSbar here */
            return gamma9[i-1] * log(m_b_MSbar / mu)
                + rho_c[i-1] * (g_c + 8.0 / 9.0 * log(m_b / m_c))
                + rho_b[i-1] * g_b
                + rho_0[i-1] * complex<double>(log(s_hat), -M_PI)
                + rho_sharp[i-1];
        }

        complex<double> f9pen(const double & s_hat) const
        {
            complex<double> g_tau = g(4.0 * power_of<2>(m_tau / m_b_pole()) / s_hat);

            return 8.0 * log(m_b_MSbar / mu)
                - 3.0 * (g_tau + 8.0 / 9.0 * log(m_b_MSbar / m_tau))
                + 8.0 / 3.0 * std::complex<double>(log(s_hat), -M_PI)
                - 40.0 / 9.0;
        }

        // cf. [HLMW:2005A], Eq. (132), p. 29
        complex<double> F(const double & s_hat) const
        {
            double r = s_hat * power_of<2>(m_b_pole() / m_c_pole()) / 4.0;
            double x = sqrt(1.0 - 1.0 / r);
            complex<double> result;

            if ((0 < r) && (r < 1))
                result = 3.0 / (2.0 * r) * (std::atan(std::sqrt(r / (1.0 - r))) / std::sqrt(r * (1.0 - r)) - 1.0);
            else
                result = 3.0 / (2.0 * r) * (
                    complex<double>(log((1.0 - x) / (1.0 + x)), M_PI) / (2.0 * sqrt(r * (r - 1.0)))
                    - 1.0);

            return result;
        }

        // cf. [HLMW:2005A], Eq. (6), p. 4
        // see also comments on removing the factor phi_u from the ratio phi_ll / phi_u below.
        double phi_ll(const double & s) const
        {
            double m_c = m_c_pole(), m_b_msbar = this->m_b_msbar();
            double m_b_pole = this->m_b_pole(), m_b_kin = model->m_b_kin(1.0);
            double log_m_l_hat = std::log(m_l / m_b_msbar);
            double m_s_hat = model->m_s_msbar(mu()) / m_b_msbar;
            double s_hat = s / power_of<2>(m_b_pole), s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;
            // We express lambda_2 as mu^2_G / 3.0 and neglect terms of
            // higher order in 1/m_b for that relation.
            double lambda_1_hat = mu2_pi / power_of<2>(m_b_kin), lambda_2_hat = mu2_g / (3.0 * power_of<2>(m_b_kin));
            double alpha_s = model->alpha_s(mu);
            double kappa = alpha_e / alpha_s, alpha_s_tilde = alpha_s / (4.0 * M_PI);

            // We do not use the ratio phi_ll / phi_u anymore. Instead, we express
            // the rate in terms of the m_b mass.
            // For this we need to eliminate corrections from BR[B->X_u l nu], both perturbative (u1,u2,uem)
            // and from the HQE (see below).
            static const double u1 = 0.0, u2 = 0.0, uem = 0.0;
            //static const double u1 = (4.0 * M_PI * M_PI - 25.0) / 12.0;
            //double u2 = 27.1 + 23.0 / 3.0 * u1 * log(mu / m_b);
            //double uem = 12.0 / 23.0 * (model->alpha_s(m_Z) / alpha_s - 1.0);

            WilsonCoefficients<BToS> w = model->wilson_coefficients_b_to_s(mu(), LeptonFlavor::muon /* fake lepton flavor */);

            // cf. [HLMW:2005A], Eq. (69), p. 16
            complex<double> c7eff = w.c7() - w.c3() / 3.0 - 4.0 * w.c4() / 9.0 - 20.0 * w.c5() / 3.0 - 80.0 * w.c6() / 9.0;
            complex<double> c8eff = w.c8() + w.c3() - w.c4() / 6.0 + 20.0 * w.c5() - 10.0 * w.c6() / 3.0;

            /* S_{AB} */
            // cf. [HLMW:2005A], Eqs. (112)-(115), p. 26
            // The HQE contributions proportional to lambda_{1,2}_hat have been adjusted to remove
            // the B->X_ulnu contributions. See also [LT:2007A].
            double s77 = power_of<2>(1.0 - s_hat) * (4.0 + 8.0 / s_hat) * (
                    1.0
                    + 8.0 * alpha_s_tilde * (omega1_77(s_hat) + u1)
                    + kappa * uem
                    + 8.0 * alpha_s_tilde * kappa * EMContributions::omegaem_77(s_hat, log_m_l_hat)
                )
                + lambda_1_hat * power_of<2>(1.0 - s_hat) * (2.0 + 4.0 / s_hat)
                + lambda_2_hat * (30.0 * s_hat2 - 18.0 - 36.0 / s_hat);

            double s79 = 12.0 * power_of<2>(1.0 - s_hat) * (
                    1.0
                    + 8.0 * alpha_s_tilde * (omega1_79(s_hat) + u1)
                    + kappa * uem
                    + 8.0 * alpha_s_tilde * kappa * EMContributions::omegaem_79(s_hat, log_m_l_hat)
                )
                + 6.0 * lambda_1_hat * power_of<2>(1.0 - s_hat)
                + 6.0 * lambda_2_hat * (-5.0 - 6.0 * s_hat + 7.0 * s_hat2);

            double s99 = power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat) * (
                    1.0
                    + 8.0 * alpha_s_tilde * (omega1_99(s_hat) + u1)
                    + kappa * uem
                    + 8.0 * alpha_s_tilde * kappa * EMContributions::omegaem_99(s_hat, log_m_l_hat)
                    + 16.0 * power_of<2>(alpha_s_tilde) * (omega2_99(s_hat) + u2 + 4.0 * u1 * omega1_99(s_hat))
                )
                + lambda_1_hat * 0.5 * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat)
                + lambda_2_hat * 1.5 * (1.0 - 15.0 * s_hat2 + 10.0 * s_hat3);

            double s1010 = s99
                + 8.0 * alpha_s_tilde * kappa * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat) * (EMContributions::omegaem_1010(s_hat, log_m_l_hat) - EMContributions::omegaem_99(s_hat, log_m_l_hat));

            /* Wilson coefficients */
            std::vector<complex<double>> wc = {
                w.c1(), w.c2(), w.c3(), w.c4(), w.c5(), w.c6(), c7eff, c8eff,
                // We use a different basis of operators: O_{9,10} = alpha_e_tilde * P_{9,10} */
                alpha_s_tilde * kappa * w.c9(),
                alpha_s_tilde * kappa * w.c10(),
                // cf. [HLMW:2005A], Table 3, p. 17. Using values at mu = 5.0 GeV
                kappa * -3.72e-2,
                kappa * -1.04e-2,
                kappa * -1.71e-6,
                kappa * -1.03e-3,
                0.0
            };

            /* Corrections, cf. [HLMW:2005A], Table 6, p. 18 */
            std::vector<complex<double>> m7 = {
                -power_of<2>(alpha_s_tilde) * kappa * memoise(CharmLoops::F17_massive, mu(), s, m_b_msbar, m_c),
                -power_of<2>(alpha_s_tilde) * kappa * memoise(CharmLoops::F27_massive, mu(), s, m_b_msbar, m_c),
                0.0,
                0.0,
                0.0,
                0.0,
                alpha_s_tilde * kappa,
                -power_of<2>(alpha_s_tilde) * kappa * CharmLoops::F87_massless(mu, s, m_b_msbar),
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            };

            std::vector<complex<double>> m9 = {
                alpha_s_tilde * kappa * f(1, s_hat) - power_of<2>(alpha_s_tilde) * kappa * memoise(CharmLoops::F19_massive, mu(), s, m_b_msbar, m_c),
                alpha_s_tilde * kappa * f(2, s_hat) - power_of<2>(alpha_s_tilde) * kappa * memoise(CharmLoops::F29_massive, mu(), s, m_b_msbar, m_c),
                alpha_s_tilde * kappa * f(3, s_hat),
                alpha_s_tilde * kappa * f(4, s_hat),
                alpha_s_tilde * kappa * f(5, s_hat),
                alpha_s_tilde * kappa * f(6, s_hat),
                0.0,
                -power_of<2>(alpha_s_tilde) * kappa * CharmLoops::F89_massless(s, m_b_msbar),
                1.0 + alpha_s_tilde * kappa * f9pen(s_hat),
                0.0,
                alpha_s_tilde * kappa * f(11, s_hat),
                alpha_s_tilde * kappa * f(12, s_hat),
                alpha_s_tilde * kappa * f(13, s_hat),
                alpha_s_tilde * kappa * f(14, s_hat),
            };

            std::vector<complex<double>> m10(14, 0.0);
            m10[9] = 1.0; // M^10_i = delta_{10,i}

            // cf. [HLMW:2005A], Eq. (111)
            double phi_ll = 0.0;
            for (unsigned i(0) ; i < 14 ; ++i)
            {
                /* diagonal */
                phi_ll += norm(wc[i]) * real(
                        s77 * norm(m7[i])
                        + s99 * norm(m9[i])
                        + s1010 * norm(m10[i])
                        + s79 * m7[i] * conj(m9[i])
                    );

                /* upper */
                for (unsigned j(i + 1) ; j < 14 ; ++j)
                {
                    phi_ll += real(wc[i] * conj(wc[j]) * (
                            2.0 * s77 * m7[i] * conj(m7[j])
                            + 2.0 * s99 * m9[i] * conj(m9[j])
                            + 2.0 * s1010 * m10[i] * conj(m10[j])
                            + s79 * (m7[i] * conj(m9[j]) + m9[i] * conj(m7[j])))
                        );
                }
            }

            /*
             * We consider also contributions from chirality-flipped operators.
             * For m_s != 0, EOS provides c7',c8' != 0 in the "SM" model.
             */
            {
                /*
                 * diagonal (i = j) and interference terms are read off the
                 * matching between Eq. (111) from [HLMW:2005A] and Eq. (3.9)
                 * from [GN:1997A].
                 */

                // We use a different basis of operators: O_{9,10} = alpha_e_tilde * P_{9,10} */
                complex<double> wc7prime  = w.c7prime();
                complex<double> wc8prime  = w.c8prime();
                complex<double> wc9prime  = alpha_s_tilde * kappa * w.c9prime();
                complex<double> wc10prime = alpha_s_tilde * kappa * w.c10prime();

                /* only chirality-flipped */
                phi_ll += norm(wc7prime) * s77   * norm(m7[7 - 1])
                    + norm(wc9prime)     * s99   * norm(m9[9 - 1])
                    + norm(wc10prime)    * s1010 * norm(m10[10 - 1])
                    + real(wc7prime * m7[7 - 1] * conj(wc9prime * m9[9 - 1]) * s79);

                double s77p =  16.0 * m_s_hat / s_hat * (1.0 - s_hat) * (2.0 - 4.0 * s_hat - s_hat2);
                double s79p = -24.0 * m_s_hat * (s_hat - s_hat2), s97p = s79p;
                double s99p =  0.5 * s79p, s1010p = s99p;
                /* interference between chirality-flipped and SM-like*/
                phi_ll += real(wc[7 - 1] * m7[7 - 1]   * conj(wc7prime  * m7[7 - 1])   * s77p)
                    + real(wc[7 - 1]  * m7[7 - 1]   * conj(wc9prime  * m9[9 - 1])   * s79p)
                    + real(wc[9 - 1]  * m9[9 - 1]   * conj(wc7prime  * m7[7 - 1])   * s97p)
                    + real(wc[9 - 1]  * m9[9 - 1]   * conj(wc9prime  * m9[9 - 1])   * s99p)
                    + real(wc[10 - 1] * m10[10 - 1] * conj(wc10prime * m10[10 - 1]) * s1010p);

                /* we also include c8' contributions */
                phi_ll += norm(wc8prime) * (
                          norm(m7[8 - 1]) * s77
                        + norm(m9[8 - 1]) * s99
                        + real(m7[8 - 1] * conj(m9[8 - 1])) * s79);
            }

            /* bremsstrahlung */
            static const double c_tau1 = 1.0 / 27.0;
            static const double c_tau2 = - 2.0 / 9.0;
            double z = power_of<2>(m_c / m_b_msbar);
            double itau_22 = real(memoise(Bremsstrahlung::itau_22, s_hat, z));
            double itau_27 = real(memoise(Bremsstrahlung::itau_27, s_hat, z));
            double itau_28 = real(memoise(Bremsstrahlung::itau_28, s_hat, z));
            double itau_29 = real(memoise(Bremsstrahlung::itau_29, s_hat, z));
            double tau_78 = memoise(Bremsstrahlung::tau_78, s_hat);
            double tau_88 = memoise(Bremsstrahlung::tau_88, s_hat);
            double tau_89 = memoise(Bremsstrahlung::tau_89, s_hat);
            double b11 = power_of<3>(alpha_s_tilde) * power_of<2>(kappa) * itau_22 * c_tau1;
            double b12 = power_of<3>(alpha_s_tilde) * power_of<2>(kappa) * itau_22 * c_tau2 * 2.0;
            double b22 = power_of<3>(alpha_s_tilde) * power_of<2>(kappa) * itau_22 * QCD::casimir_f;
            double b17 = power_of<3>(alpha_s_tilde) * power_of<2>(kappa) * itau_27 * c_tau2 * 2.0;
            double b27 = power_of<3>(alpha_s_tilde) * power_of<2>(kappa) * itau_27 * QCD::casimir_f * 2.0;
            double b18 = power_of<3>(alpha_s_tilde) * power_of<2>(kappa) * itau_28 * c_tau2 * 2.0;
            double b28 = power_of<3>(alpha_s_tilde) * power_of<2>(kappa) * itau_28 * QCD::casimir_f * 2.0;
            double b19 = power_of<2>(alpha_s_tilde) * kappa * itau_29 * c_tau2 * 2.0;
            double b29 = power_of<2>(alpha_s_tilde) * kappa * itau_29 * QCD::casimir_f * 2.0;
            double b78 = power_of<3>(alpha_s_tilde) * power_of<2>(kappa) * tau_78 * 2.0 * QCD::casimir_f;
            double b88 = power_of<3>(alpha_s_tilde) * power_of<2>(kappa) * tau_88 * QCD::casimir_f;
            double b89 = power_of<3>(alpha_s_tilde) * power_of<2>(kappa) * tau_89 * 2.0 * QCD::casimir_f;
            phi_ll += norm(wc[1 - 1]) * b11 + real(wc[1 - 1] * conj(wc[2 - 1])) * b12 + norm(wc[2 - 1]) * b22;
            phi_ll += real(conj(wc[7 - 1]) * (wc[1 - 1] * b17 + wc[2 - 1] * b27));
            phi_ll += real(conj(wc[8 - 1]) * (wc[1 - 1] * b18 + wc[2 - 1] * b28));
            phi_ll += real(conj(wc[9 - 1]) * (wc[1 - 1] * b19 + wc[2 - 1] * b29));
            phi_ll += real(conj(wc[8 - 1]) * wc[7 - 1]) * b78;
            phi_ll += real(conj(wc[9 - 1]) * wc[8 - 1]) * b89;
            phi_ll += norm(wc[8 - 1]) * b88;

            /* non-perturbative 1/m_c^2 */
            complex<double> cF = F(s_hat);
            // We use lambda_2 = mu2_g / 3.0 and neglect higher orders in 1/m_b in that
            // relation.
            double c27 = - power_of<2>(alpha_s_tilde * kappa) * 8.0 * mu2_g / (27.0 * power_of<2>(m_c)) * power_of<2>(1.0 - s_hat)
                * (1.0 + 6.0 * s_hat - s_hat2) / s_hat * real(cF);
            double c29 = - alpha_s_tilde * kappa * 8.0 * mu2_g / (27.0 * power_of<2>(m_c)) * power_of<2>(1.0 - s_hat) * (2.0 + s_hat) * real(cF);
            double c22 = - alpha_s_tilde * kappa * 8.0 * mu2_g / (27.0 * power_of<2>(m_c)) * power_of<2>(1.0 - s_hat) * (2.0 + s_hat) * real(cF * conj(m9[1]));
            phi_ll += c22 * (-2.0 / 9.0 * norm(wc[0]) + 7.0 / 6.0 * real(wc[0] * conj(wc[1])) + norm(wc[1]));
            phi_ll += c27 * real((-1.0 / 6.0 * wc[0] + wc[1]) * conj(wc[6]));
            phi_ll += c29 * real((-1.0 / 6.0 * wc[0] + wc[1]) * conj(wc[8]));

            /* log enhanced em */
            double e22 = 8.0 * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat) * power_of<3>(alpha_s_tilde * kappa) * EMContributions::omegaem_22(s_hat, log_m_l_hat, mu);
            complex<double> e27 = 96.0 * power_of<2>(1.0 - s_hat) * power_of<3>(alpha_s_tilde * kappa) * EMContributions::omegaem_27(s_hat, log_m_l_hat, mu);
            complex<double> e29 = 8.0 * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat) * power_of<2>(alpha_s_tilde * kappa) * EMContributions::omegaem_29(s_hat, log_m_l_hat, mu);
            phi_ll += e22 * (16.0 / 9.0 * norm(wc[0]) + 8.0 / 3.0 * real(wc[0] * conj(wc[1])) + norm(wc[1]));
            phi_ll += real(e27 * (4.0 / 3.0 * wc[0] + wc[1]) * conj(wc[6]));
            phi_ll += real(e29 * (4.0 / 3.0 * wc[0] + wc[1]) * conj(wc[8]));

            return phi_ll;
        }

        // cf. [HLMW:2005A], Eq. (4), p. 4
        double branching_ratio(const double & s) const
        {
            static const double pi3 = power_of<3>(M_PI);

            return power_of<2>(gfermi()) * power_of<5>(m_b_pole())
                * std::norm(model->ckm_tb() * std::conj(model->ckm_ts())) * phi_ll(s) * tau_B()
                / (48.0 * pi3 * hbar());
        }

        // diagnostic values
        Diagnostics diagnostics() const
        {
            Diagnostics results;

            // Function phi_ll, cf. [HLMW:2005A]
            {
                results.add(Diagnostics::Entry{ phi_ll(1.0), "phi_ll(s = 1.0Gev^2), [HLMW:2005A]" });
                results.add(Diagnostics::Entry{ phi_ll(6.0), "phi_ll(s = 6.0Gev^2), [HLMW:2005A]" });
            }

            return results;
        }
    };

    BToXsDilepton<HLMW2005>::BToXsDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToXsDilepton<HLMW2005>>(new Implementation<BToXsDilepton<HLMW2005>>(parameters, options, *this))
    {
    }

    BToXsDilepton<HLMW2005>::~BToXsDilepton()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<BToXsDilepton<HLMW2005>>::options
    {
        { "l"_ok, { "e"s, "mu"s, "tau"s }, "mu"s },
        { "q"_ok, { "d"s, "u"s }, "d"s }
    };

    double
    BToXsDilepton<HLMW2005>::differential_branching_ratio(const double & s) const
    {
        return _imp->branching_ratio(s) / power_of<2>(_imp->m_b_pole());
    }

    double
    BToXsDilepton<HLMW2005>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrate<1>(std::function<double (const double &)>(
                    std::bind(&BToXsDilepton<HLMW2005>::differential_branching_ratio, this, std::placeholders::_1)),
                    s_min, s_max, _imp->cub_conf);
    }

    Diagnostics
    BToXsDilepton<HLMW2005>::diagnostics() const
    {
        return _imp->diagnostics();
    }

    const std::set<ReferenceName>
    BToXsDilepton<HLMW2005>::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BToXsDilepton<HLMW2005>::begin_options()
    {
        return Implementation<BToXsDilepton<HLMW2005>>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToXsDilepton<HLMW2005>::end_options()
    {
        return Implementation<BToXsDilepton<HLMW2005>>::options.cend();
    }
}
