/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2025 Danny van Dyk
 * Copyright (c) 2018      Nico Gubernari
 * Copyright (c) 2018      Ahmet Kokulu
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_P_LCSR_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_P_LCSR_IMPL_HH 1

#include <eos/form-factors/analytic-b-to-p-lcsr.hh>
#include <eos/form-factors/heavy-meson-lcdas.hh>
#include <eos/utils/exception.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/kinematic.hh>
#include <eos/models/model.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

#include <functional>
#include <memory>

#include <boost/predef.h>

#if BOOST_COMP_GNUC
#  pragma GCC optimize("no-var-tracking")
#endif

namespace eos
{
    using namespace std::literals::string_literals;

    template <typename Transition_>
    struct Implementation<AnalyticFormFactorBToPLCSR<Transition_>>
    {
        std::shared_ptr<Model> model;

        // B-meson parameters
        UsedParameter m_B;
        UsedParameter f_B;

        // final state meson parameters
        UsedParameter m_P;
        UsedParameter f_P;

        // sum rule parameters
        UsedParameter s0_0_p;
        UsedParameter s0_1_p;
        UsedParameter s0_0_pm;
        UsedParameter s0_1_pm;
        UsedParameter s0_0_t;
        UsedParameter s0_1_t;
        UsedParameter M2;

        // properties of the virtual quark in the sum rule's correlator
        std::function<double ()> m_v;

        // renormalization scale
        UsedParameter mu;

        std::shared_ptr<HeavyMesonLCDAs> b_lcdas;

        // switches to enable/disable certain contributions
        SwitchOption opt_2pt;
        SwitchOption opt_3pt;
        double switch_2pt_phi;
        double switch_2pt_g;
        double switch_3pt;

        // switch to select the QHD matching method (fully dispersive representation vs Borel transformation)
        SwitchOption opt_method;
        std::function<double (const Implementation *, const double &, const double &)> integrand_fp_2pt;
        std::function<double (const Implementation *, const double &, const double &)> integrand_fpm_2pt;
        std::function<double (const Implementation *, const double &, const double &)> integrand_fT_2pt;
        bool switch_borel;

        static const std::vector<OptionSpecification> options;

        using Traits = AnalyticFormFactorBToPLCSRProcessTraits<Transition_>;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            m_B(p[Traits::name_B], u),
            f_B(p[Traits::f_B], u),
            m_P(p[Traits::name_P], u),
            f_P(p[Traits::f_P], u),
            s0_0_p(p[stringify(Traits::label) + "::s_0^+,0@B-LCSR"], u),
            s0_1_p(p[stringify(Traits::label) + "::s_0^+,1@B-LCSR"], u),
            s0_0_pm(p[stringify(Traits::label) + "::s_0^+/-,0@B-LCSR"], u),
            s0_1_pm(p[stringify(Traits::label) + "::s_0^+/-,1@B-LCSR"], u),
            s0_0_t(p[stringify(Traits::label) + "::s_0^T,0@B-LCSR"], u),
            s0_1_t(p[stringify(Traits::label) + "::s_0^T,1@B-LCSR"], u),
            M2(p[stringify(Traits::label) + "::M^2@B-LCSR"], u),
            mu(p[stringify(Traits::label) + "::mu@B-LCSR"], u),
            b_lcdas(HeavyMesonLCDAs::make("exponential", p, o + Options{ { "q"_ok, stringify(Traits::spectator_flavor) } })), // operator+ is ordered!
            opt_2pt(o, "2pt"_ok, { "tw2+3", "all", "off" }, "all"),
            opt_3pt(o, "3pt"_ok, { "tw3+4", "all", "off" }, "all"),
            switch_2pt_phi(1.0),
            switch_2pt_g(1.0),
            switch_3pt(1.0),
            opt_method(o, "method"_ok, { "borel", "dispersive" }, "borel"),
            switch_borel(opt_method.value() == "borel")
        {
            Context ctx("When creating a B->P LCSR form factor with B-meson LCDAs");

            u.uses(*b_lcdas);

            // quark masses for the propagating quark
            const QuarkFlavor q_v = std::get<1>(Traits::partonic_transition);
            switch (q_v)
            {
                case QuarkFlavor::up:
                case QuarkFlavor::down:
                    m_v = [this]() -> double { return this->model->m_ud_msbar(this->mu()) / 2.0; };
                    break;

                case QuarkFlavor::strange:
                    m_v = [this]() -> double { return this->model->m_s_msbar(this->mu()); };
                    break;

                case QuarkFlavor::charm:
                    m_v = [this]() -> double { return this->model->m_c_msbar(this->mu()); };
                    break;

                default:
                    throw InternalError("Unknown valence quark flavor: '" + stringify(q_v) + "'");
            }

            // selectively enable/disable two-particle contributions
            if ("off" == opt_2pt.value())
            {
                switch_2pt_phi = 0.0;
                switch_2pt_g   = 0.0;
            }
            else if ("tw2+3" == opt_2pt.value())
            {
                switch_2pt_phi = 1.0;
                switch_2pt_g   = 0.0;
            }

            // selectively enable/disable three-particle contributions
            if ("off" == opt_3pt.value())
            {
                switch_3pt = 0.0;
            }

            // select the apropriate integrand based on method for the QHD matching
            if ("borel" == opt_method.value())
            {
                integrand_fp_2pt  = &Implementation::integrand_fp_2pt_borel;
                integrand_fpm_2pt = &Implementation::integrand_fpm_2pt_borel;
                integrand_fT_2pt  = &Implementation::integrand_fT_2pt_borel;
            }
            else
            {
                integrand_fp_2pt  = &Implementation::integrand_fp_2pt_disp;
                integrand_fpm_2pt = &Implementation::integrand_fpm_2pt_disp;
                integrand_fT_2pt  = &Implementation::integrand_fT_2pt_disp;
            }

        }

        ~Implementation() = default;

        /* forwarding the LCDAs */
        inline
        double phi_plus(const double & omega) const
        {
            return switch_2pt_phi * b_lcdas->phi_plus(omega);
        }

        inline
        double phi_bar(const double & omega) const
        {
            return switch_2pt_phi * b_lcdas->phi_bar(omega);
        }

        inline
        double phi_bar_d1(const double & omega) const
        {
            return switch_2pt_phi * b_lcdas->phi_bar_d1(omega);
        }

        inline
        double g_plus(const double & omega) const
        {
            return switch_2pt_g * b_lcdas->g_plus(omega);
        }

        inline
        double g_plus_d1(const double & omega) const
        {
            return switch_2pt_g * b_lcdas->g_plus_d1(omega);
        }

        inline
        double g_plus_d2(const double & omega) const
        {
            return switch_2pt_g * b_lcdas->g_plus_d2(omega);
        }

        inline
        double g_bar(const double & omega) const
        {
            return switch_2pt_g * b_lcdas->g_bar(omega);
        }

        inline
        double g_bar_d1(const double & omega) const
        {
            return switch_2pt_g * b_lcdas->g_bar_d1(omega);
        }

        inline
        double g_bar_d2(const double & omega) const
        {
            return switch_2pt_g * b_lcdas->g_bar_d2(omega);
        }

        inline
        double g_bar_d3(const double & omega) const
        {
            return switch_2pt_g * b_lcdas->g_bar_d3(omega);
        }

        inline
        double phi_3(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->phi_3(omega_1, omega_2);
        }

        inline
        double phi_bar_3(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->phi_bar_3(omega_1, omega_2);
        }

        inline
        double phi_bar2_3(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->phi_bar2_3(omega_1, omega_2);
        }

        inline
        double phi_bar_bar_3(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->phi_bar_bar_3(omega_1, omega_2);
        }

        inline
        double phi_4(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->phi_4(omega_1, omega_2);
        }

        inline
        double phi_bar_4(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->phi_bar_4(omega_1, omega_2);
        }

        inline
        double phi_bar2_4(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->phi_bar2_4(omega_1, omega_2);
        }

        inline
        double phi_bar_bar_4(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->phi_bar_bar_4(omega_1, omega_2);
        }

        inline
        double psi_bar_4(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->psi_bar_4(omega_1, omega_2);
        }

        inline
        double psi_bar_bar_4(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->psi_bar_bar_4(omega_1, omega_2);
        }

        inline
        double chi_bar_4(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->chi_bar_4(omega_1, omega_2);
        }

        inline
        double chi_bar_bar_4(const double omega_1, const double omega_2) const
        {
            return switch_3pt * b_lcdas->chi_bar_bar_4(omega_1, omega_2);
        }

        /* auxilliary functions */

        double s(const double & sigma, const double & q2) const
        {
            const double sigmabar = 1.0 - sigma;

            return sigma * power_of<2>(m_B()) + (power_of<2>(m_v()) - sigma * q2) / sigmabar;
        }

        double sigma(const double & s, const double & q2) const
        {
            const double m_B2 = power_of<2>(m_B());
            const double m_v2 = power_of<2>(m_v());

            return (m_B2 - q2 + s - std::sqrt(4.0 * (m_v2 - s) * m_B2 + power_of<2>(m_B2 - q2 + s))) / (2.0 * m_B2);
        }

        double sigma_0(const double & q2, const double & s0_0, const double & s0_1) const
        {
            const double s0 = s0_0 + s0_1 * q2;

            return sigma(s0, q2);
        }

        /* f_+ : 2-particle functions */

        inline
        double I1_fp_2pt_phi_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to phi_+

            //const double sigmabar = 1.0 - sigma;
            const double C_1 = -1.0; /* * sigmabar */

            const double phi_plus  = this->phi_plus(sigma * m_B());

            return C_1 * phi_plus; /* / power_of<1>(sigmabar) */
        }

        inline
        double I2_fp_2pt_phi_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to phibar

            // const double sigmabar = 1.0 - sigma;

            const double phi_bar  = this->phi_bar(sigma * m_B);

            const double C_2 = -m_B(); /* * power_of<2>(sigmabar)*/

            return C_2 * phi_bar; /* / power_of<2>(sigmabar) */
        }

        inline
        double I2d1_fp_2pt_phi_bar(const double & sigma, const double & /*q2*/) const
        {
            // first derivative of two-particle contribution to f_+ proportional to phibar

            // const double sigmabar = 1.0 - sigma;

            const double phi_bar_d1  = this->phi_bar_d1(sigma * m_B);

            const double C_2 = -m_B(); /* * power_of<2>(sigmabar)*/

            return C_2 * (m_B * phi_bar_d1); /* / power_of<2>(sigmabar) */
        }

        inline
        double I2_fp_2pt_g_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_+
            const double sigmabar = 1.0 - sigma;

            const double g_plus   = this->g_plus(sigma * m_B);

            return -4.0 * g_plus / sigmabar;
        }

        inline
        double I2d1_fp_2pt_g_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_+
            const double sigmabar = 1.0 - sigma;

            const double g_plus    = this->g_plus(sigma * m_B);
            const double g_plus_d1 = this->g_plus_d1(sigma * m_B()) * m_B();

            return -4.0 * (sigmabar * g_plus_d1 + g_plus) / power_of<2>(sigmabar);
        }

        inline
        double I3_fp_2pt_g_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_+
            const double m_v2     = power_of<2>(m_v());
            const double sigmabar = 1.0 - sigma;

            const double g_plus    = this->g_plus(sigma * m_B());

            return +8.0 * m_v2 * g_plus / power_of<2>(sigmabar);
        }

        inline
        double I3d1_fp_2pt_g_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_+
            const double m_v2     = power_of<2>(m_v());
            const double sigmabar = 1.0 - sigma;

            const double g_plus    = this->g_plus(sigma * m_B());
            const double g_plus_d1 = this->g_plus_d1(sigma * m_B()) * m_B();

            return +8.0 * m_v2 * (g_plus_d1 * sigmabar + 2.0 * g_plus) / power_of<3>(sigmabar);
        }

        inline
        double I3d2_fp_2pt_g_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_+
            const double m_v2     = power_of<2>(m_v());
            const double sigmabar = 1.0 - sigma;

            const double g_plus    = this->g_plus(sigma * m_B());
            const double g_plus_d1 = this->g_plus_d1(sigma * m_B()) * m_B();
            const double g_plus_d2 = this->g_plus_d2(sigma * m_B()) * power_of<2>(m_B());

            return +8.0 * m_v2 * (g_plus_d2 * power_of<2>(sigmabar) + 4.0 * g_plus_d1 * sigmabar + 6.0 * g_plus) / power_of<4>(sigmabar);
        }

        inline
        double I3_fp_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_bar
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B());

            return -8.0 * m_B() * g_bar / sigmabar;
        }

        inline
        double I3d1_fp_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_bar
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B());
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B());

            return -8.0 * m_B() * (g_bar_d1 * sigmabar * m_B + g_bar) / power_of<2>(sigmabar);
        }

        inline
        double I3d2_fp_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_bar
            const double m_B2     = power_of<2>(m_B());
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B());
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B());
            const double g_bar_d2 = this->g_bar_d2(sigma * m_B());

            return -8.0 * m_B() * (g_bar_d2 * power_of<2>(sigmabar) * m_B2 + 2.0 * g_bar_d1 * sigmabar * m_B() + 2.0 * g_bar) / power_of<3>(sigmabar);
        }

        inline
        double I4_fp_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_bar
            const double m_v2     = power_of<2>(m_v());
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B());

            return 24.0 * m_B() * m_v2 * g_bar / power_of<2>(sigmabar);
        }

        inline
        double I4d1_fp_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_bar
            const double m_v2     = power_of<2>(m_v());
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B());
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B()) * m_B();

            return 24.0 * m_B() * m_v2 * (g_bar_d1 * sigmabar + 2.0 * g_bar) / power_of<3>(sigmabar);
        }

        inline
        double I4d2_fp_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_bar
            const double m_v2     = power_of<2>(m_v());
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B());
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B()) * m_B();
            const double g_bar_d2 = this->g_bar_d2(sigma * m_B()) * power_of<2>(m_B());

            return 24.0 * m_B() * m_v2 * (g_bar_d2 * power_of<2>(sigmabar) + 4.0 * g_bar_d1 * sigmabar + 6.0 * g_bar) / power_of<4>(sigmabar);
        }

        inline
        double I4d3_fp_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_+ proportional to g_bar
            const double m_v2     = power_of<2>(m_v());
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B());
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B()) * m_B();
            const double g_bar_d2 = this->g_bar_d2(sigma * m_B()) * power_of<2>(m_B());
            const double g_bar_d3 = this->g_bar_d3(sigma * m_B()) * power_of<3>(m_B());

            return 24.0 * m_B() * m_v2 * (g_bar_d3 * power_of<3>(sigmabar) + 6.0 * g_bar_d2 * power_of<2>(sigmabar) + 18.0 * g_bar_d1 * sigmabar + 24.0 * g_bar) / power_of<5>(sigmabar);
        }

        /* f_+ : 3-particle functions */
        // {{{
        double I2_fp_3pt_phi_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_3
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_3 = this->phi_3(omega_1, omega_2);

            const double C_2 = -(m_B * sigmabar * u + 2.0 * m_v) / (m_B * power_of<2>(sigmabar));

            return C_2 * phi_3;
        }

        double I2_fp_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_3
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_2 = u / (m_B * power_of<2>(sigmabar));

            return C_2 * phi_bar_3;
        }

        double I3_fp_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_3
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = -2.0 * (u * (m_B2 * sigmabar2 + q2) + 4.0 * m_B * m_v * sigmabar + m_v2 * u)
                             / (m_B * power_of<3>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I3d1A_fp_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_3
            const double m_B2     = power_of<2>(m_B),   m_B3     = power_of<3>(m_B);
            const double m_v      = this->m_v(),   m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = -((4.0 * m_B * m_v * omega_2 * sigmabar * (3.0 + sigmabar) + sigma3 * (-(4.0 * m_B * q2) + 6.0 * m_B3 * sigmabar) +
                               2.0 * sigma2 * (2.0 * omega_1 * q2 + 2.0 * m_B * (m_v2 + q2) - 3.0 * m_B2 * omega_1 * sigmabar -
                               3.0 * m_B * q2 * sigmabar + 3.0 * m_B3 * sigmabar * (-2.0 + sigmabar)) +
                               sigmabar * (2.0 * m_B * (m_B2 + q2) * sigmabar + omega_1 * q2 * (-7.0 + sigmabar) +
                               m_B2 * omega_1 * (-6.0 + 4.0 * sigmabar)) +
                               m_v2 * (-(2.0 * omega_2 * (-1 + sigmabar)) + m_B * sigmabar * (-1 + 3.0 * sigmabar) +
                               2.0 * omega_1 * (2.0 - 5.0 * sigmabar)) -
                               sigma * (2.0 * m_v2 * (2.0 * omega_1 + omega_2) + 4.0 * m_B2 * omega_1 * sigmabar * (-3.0 + sigmabar) +
                               2.0 * m_B3 * sigmabar * (-3.0 + 4.0 * sigmabar) + omega_1 * q2 * (4.0 - 5.0 * sigmabar) +
                               m_B * (12.0 * m_v * omega_2 * sigmabar + 2.0 * q2 * sigmabar * (-4.0 + sigmabar) + m_v2 * (4.0 - 11.0 * sigmabar)))))
                               / (m_B * omega_2 * power_of<5>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I3d1B_fp_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_3
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = (2.0 * (sigma - 2.0) * sigmabar * (2.0 * m_B2 * sigma - q2) + 4.0 * m_B2 * sigmabar + 8.0 * m_B * m_v * sigmabar
                             * (-sigma + sigmabar + 1.0) + m_v2 * (sigma + 5.0 * sigmabar - 1.0) - 2.0 * q2 * (sigma - 1.0) * sigma)
                             / (2.0 * (-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I3d1C_fp_3pt_phi_bar_3(const double & sigma, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_3
            const double omega_1  = m_B * sigma;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = -8.0 * m_B * m_v / (omega_2 * power_of<2>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I4_fp_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_3
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = -6.0 * m_v * u * (2.0 * m_B * sigmabar + m_v * (2.0 * u - 1.0))
                             / (power_of<3>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d1A_fp_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_3
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * (m_B * ((-(2.0 * m_B) + m_v) * omega_2 + 6.0 * m_B * (m_v - omega_2) * sigma2 -
                               2.0 * (m_v * omega_2 + 2.0 * m_B * (m_v - 2.0 * omega_2)) * sigma) * sigmabar -
                               4.0 * m_B * sigma * sigmabar * (2.0 * m_B * m_v * sigma - omega_2 * (m_v - 2.0 * m_B * sigmabar)) +
                               omega_1 * (4.0 * sigmabar * (4.0 * m_B * m_v * sigma - omega_2 * (m_v - 2.0 * m_B * sigmabar)) +
                               sigmabar * (m_v * omega_2 + 4.0 * m_B * (m_v - 2.0 * m_v * sigma - omega_2 * sigmabar))) +
                               2.0 * m_v * (-4.0 + 4.0 * sigma + sigmabar) * power_of<2>(omega_1))
                             / (power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d1B_fp_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_3
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * (m_B - m_B * sigma) * (2.0 * m_B * sigmabar + m_v)
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * phi_bar_bar_3;
        }

        double I4d1C_fp_3pt_phi_bar_bar_3(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_3

            return 0.0;
        }
        double I4d2A_fp_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_3
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = 12.0 * m_v * (m_B *
                              (10.0 * (-1 + sigma) * sigma * (-(omega_2 * (m_v + 2.0 * m_B * (-1 + sigma))) + 2.0 * m_B * m_v * sigma) +
                               sigmabar2 * (-(m_v * omega_2) + m_B * (4.0 * omega_2 - 6.0 * omega_2 * sigma + m_v * (-2.0 + 6.0 * sigma))) +
                               4.0 * ((-(2.0 * m_B) + m_v) * omega_2 + 6.0 * m_B * (m_v - omega_2) * sigma2 -
                               2.0 * (m_v * omega_2 + 2.0 * m_B * (m_v - 2.0 * omega_2)) * sigma) * sigmabar) +
                               2.0 * omega_1 * (m_B * (-(2.0 * m_v) + omega_2) * sigmabar2 -
                               5.0 * (-1 + sigma) * (-(omega_2 * (m_v + 2.0 * m_B * (-1 + sigma))) + 4.0 * m_B * m_v * sigma) +
                               2.0 * (m_v * omega_2 + 4.0 * m_B * (m_v + omega_2 * (-1 + sigma) - 2.0 * m_v * sigma)) * sigmabar) +
                               4.0 * m_v * (-5.0 + 5.0 * sigma + 2.0 * sigmabar) * power_of<2>(omega_1))
                             / (power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d2B_fp_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_3
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3     = this->phi_bar_3(omega_1, omega_2);
            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            return  6.0 * m_B * m_v *
                   (2.0 * (4.0 * m_B * (-m_v + 2.0 * m_B * (-1 + sigma)) * (-1 + sigma) * sigma +
                    m_B * (m_v - 2.0 * m_v * sigma + 4.0 * m_B * (-1 + sigma) * sigma) * sigmabar +
                    omega_1 * (-(4.0 * m_B * (-1 + sigma) * (-2.0 + 2.0 * sigma + sigmabar)) + m_v * (-4.0 + 4.0 * sigma + sigmabar))) *
                    phi_bar_bar_3 + m_B * (-m_v + 2.0 * m_B * (-1 + sigma)) * (-1 + sigma) * (-omega_1 + m_B * sigma) * sigmabar * phi_bar_3)
                    / (power_of<2>(omega_2) * power_of<5>(sigmabar));
        }

        double I4d2C_fp_3pt_phi_bar_bar_3(const double & sigma, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_3
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = 6.0 * m_B2 * m_v * (m_v - 2.0 * m_B * sigmabar) / (power_of<2>(omega_2) * power_of<3>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d2D_fp_3pt_phi_bar_bar_3(const double & sigma, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_3
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3     = this->phi_bar_3(omega_1, omega_2);
            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            return     (6.0 * m_v * ((-4.0 * m_B * sigmabar - m_v) * phi_bar_bar_3
                      - m_B * sigmabar * (-2.0 * m_B * sigmabar - m_v) * phi_bar_3)) / (power_of<4>(sigmabar));
        }

        double I2_fp_3pt_phi_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_4
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_4 = this->phi_4(omega_1, omega_2);

            const double C_2 = -(u - 1.0) / sigmabar;

            return C_2 * phi_4;
        }

        double I2_fp_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_4
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_2 = (u - 1.0) / (m_B * power_of<2>(sigmabar));

            return C_2 * phi_bar_4;
        }

        double I3_fp_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(),  m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma,  sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (m_B2 * sigmabar2 * u + 2.0 * m_B * m_v * sigmabar + m_v2 * (-(u - 1.0)) - q2 * u + q2)
                             / (m_B * power_of<3>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3d1A_fp_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = (sigma3 * (4.0 * m_B * q2 + 6.0 * m_B3 * sigmabar) +
                               2.0 * sigma2 * (-(2.0 * (omega_1 + omega_2) * q2) - 3.0 * m_B2 * omega_1 * sigmabar +
                               3.0 * m_B3 * sigmabar * (-2.0 + sigmabar) + m_B * q2 * (-2.0 + 3.0 * sigmabar)) +
                               sigmabar * (8.0 * m_B * m_v * omega_2 * sigmabar + m_v2 * (6.0 * omega_1 + 6.0 * omega_2 - 2.0 * m_B * sigmabar) -
                               q2 * (2.0 * m_B * sigmabar + omega_1 * (-7.0 + sigmabar) + omega_2 * (-7.0 + sigmabar)) +
                               2.0 * m_B2 * (m_B * sigmabar + omega_1 * (-3.0 + 2.0 * sigmabar))) +
                               sigma * (-(2.0 * m_B * sigmabar * (3.0 * m_v2 + 2.0 * m_B * omega_1 * (-3.0 + sigmabar) +
                               m_B2 * (-3.0 + 4.0 * sigmabar))) +
                               q2 * (4.0 * omega_1 - 5.0 * omega_1 * sigmabar + 2.0 * m_B * sigmabar * (-4.0 + sigmabar) +
                               omega_2 * (4.0 - 5.0 * sigmabar))))
                             / (m_B * omega_2 * power_of<5>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3d1B_fp_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * m_B * (m_B * sigmabar2 + 2.0 * m_v * sigmabar)
                             / ((-omega_1 + m_B * sigma) * power_of<3>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3d1C_fp_3pt_phi_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_4
            const double omega_1  = m_B * sigma;

            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;
            const double sigma2   = power_of<2>(sigma);

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = (2.0 * sigmabar * (2.0 * m_B * m_v * sigmabar + m_v2 + q2) - q2 * sigma2 + sigma * (q2 - q2 * sigmabar))
                             / (omega_2 * power_of<4>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3_fp_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * u * (m_B * sigmabar * (2.0 * u - 1.0) + 2.0 * m_v)
                             / (m_B * power_of<3>(sigmabar));

            return C_3 * phi_bar_bar_4;
        }

        double I3d1A_fp_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B),   m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2   = power_of<2>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (4.0 * m_B3 * sigma2 * sigmabar - m_B2 * sigmabar2 * (omega_2 + 4.0 * omega_1) +
                               2.0 * m_B * sigma * (m_B * sigmabar *  (2.0 * m_B * sigmabar - omega_2 - 4.0 * omega_1) + 3.0 * m_v * omega_2) +
                               2.0 * m_B * sigmabar * (m_v * omega_2 + omega_1 * (omega_2 + 2.0 * omega_1)) - 6.0 * m_v * omega_2 * omega_1)
                             / (m_B * power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_3 * phi_bar_bar_4;
        }

        double I3d1B_fp_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (m_B * sigmabar + 2.0 * m_v)
                             / ((-omega_1 + m_B * sigma) * power_of<3>(sigmabar));

            return C_3 * phi_bar_bar_4;
        }

        double I3d1C_fp_3pt_phi_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4

            return 0.0;
        }

        double I4_fp_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * u * (m_B2 * sigmabar2 - q2) * (m_B * sigmabar * (2.0 * u - 1.0) + 2.0 * m_v)
                             / (m_B * power_of<4>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d1A_fp_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * (-(m_B4 * (4.0 * omega_1 + omega_2) * sigmabar2) + 8.0 * m_v * omega_1 * omega_2 * q2 +
                               6.0 * m_B5 * sigma4 * sigmabar - m_B * (2.0 * m_v * omega_2 + 3.0 * omega_1 * (2.0 * omega_1 + omega_2)) * q2 *
                               sigmabar + m_B2 * ((4.0 * omega_1 + omega_2) * sigmabar2 * q2 +
                               4.0 * m_v * omega_1 * omega_2 * (-2.0 + sigmabar)) +
                               m_B3 * sigmabar * (2.0 * m_v * omega_2 - omega_1 * (2.0 * omega_1 + omega_2) * (-3.0 + 2.0 * sigmabar)) +
                               m_B3 * sigma3 * (8.0 * m_v * omega_2 + m_B * sigmabar *
                               (-(3.0 * (4.0 * omega_1 + omega_2)) + 4.0 * m_B * (-3.0 + 2.0 * sigmabar))) +
                               m_B2 * sigma2 * (-(8.0 * m_v * omega_1 * omega_2) -
                               3.0 * m_B2 * (4.0 * omega_1 + omega_2) * sigmabar * (-2.0 + sigmabar) +
                               6.0 * m_B3 * sigmabar * (1 - 2.0 * sigmabar) +
                               m_B * (3.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar - 6.0 * q2 * sigmabar +
                               2.0 * m_v * omega_2 * (-8.0 + 3.0 * sigmabar))) +
                               m_B * sigma * (4.0 * m_B4 * sigmabar2 - 8.0 * m_v * omega_2 * q2 +
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar * (-3.0 + 4.0 * sigmabar) +
                               m_B * (3.0 * (4.0 * omega_1 + omega_2) * q2 * sigmabar - 4.0 * m_v * omega_1 * omega_2 * (-4.0 + sigmabar)) -
                               2.0 * m_B2 * (4.0 * m_v * omega_2 * (-1 + sigmabar) +
                               sigmabar * (2.0 * q2 * sigmabar - omega_1 * (2.0 * omega_1 + omega_2) * (-3.0 + sigmabar)))))
                             / (m_B * power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d1B_fp_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * (m_B * sigmabar + 2.0 * m_v) * (m_B2 * sigmabar2 -q2)
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * phi_bar_bar_4;
        }

        double I4d1C_fp_3pt_phi_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4

            return 0.0;
        }

        double I4d2A_fp_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 =12.0 * (2.0 * m_B5 * sigmabar3 + 20.0 * m_v * omega_1 * omega_2 * q2 + 12.0 * m_B5 * sigma4 * sigmabar -
                              2.0 * m_B * (4.0 * m_v * omega_2 + 3.0 * omega_1 * (2.0 * omega_1 + omega_2)) * q2 * sigmabar +
                              m_B4 * (4.0 * omega_1 + omega_2) * sigmabar2 * (-3.0 + 2.0 * sigmabar) +
                              2.0 * m_B3 * sigma3 * (10.0 * m_v * omega_2 +
                              3.0 * m_B * (-(4.0 * omega_1) - omega_2 + 4.0 * m_B * (-1 + sigmabar)) * sigmabar) +
                              m_B3 * sigmabar * (-(2.0 * sigmabar2 * q2) - 4.0 * m_v * omega_2 * (-2.0 + sigmabar) +
                              omega_1 * (2.0 * omega_1 + omega_2) * (6.0 + sigmabar * (-6.0 + sigmabar))) +
                              m_B2 * (3.0 * (4.0 * omega_1 + omega_2) * sigmabar2 * q2 -
                              2.0 * m_v * omega_1 * omega_2 * (10.0 + sigmabar * (-8.0 + sigmabar))) +
                              m_B2 * sigma2 * (-(20.0 * m_v * omega_1 * omega_2) -
                              3.0 * m_B2 * (4.0 * omega_1 + omega_2) * sigmabar * (-4.0 + 3.0 * sigmabar) +
                              12.0 * m_B3 * sigmabar * (1 + sigmabar * (-3.0 + sigmabar)) +
                              2.0 * m_B * (3.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar - 6.0 * q2 * sigmabar +
                              4.0 * m_v * omega_2 * (-5.0 + 3.0 * sigmabar))) +
                              m_B * sigma * (-(20.0 * m_v * omega_2 * q2) - 12.0 * m_B4 * sigmabar2 * (-1 + sigmabar) +
                              6.0 * m_B * omega_2 * q2 * sigmabar - 3.0 * m_B3 * (4.0 * omega_1 + omega_2) * sigmabar *
                             (2.0 + sigmabar * (-4.0 + sigmabar)) +
                              8.0 * m_B * omega_1 * (3.0 * q2 * sigmabar + m_v * omega_2 * (5.0 - 2.0 * sigmabar)) +
                              2.0 * m_B2 * (3.0 * sigmabar * (-(2.0 * q2 * sigmabar) + omega_1 * (2.0 * omega_1 + omega_2) * (-2.0 + sigmabar)) +
                              m_v * omega_2 * (10.0 + sigmabar * (-16.0 + 3.0 * sigmabar)))))
                             / (m_B * power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d2B_fp_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_4     = this->phi_bar_4(omega_1, omega_2);
            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            return  -(6.0 * pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-5.0) *
                     (2.0 * (omega_1 * (8.0 * m_v * q2 + 3.0 * m_B * q2 * sigmabar +
                      4.0 * m_B2 * m_v * sigmabar * (-2.0 + 2.0 * sigma + sigmabar) +
                      m_B3 * sigmabar2 * (-3.0 + 3.0 * sigma + 2.0 * sigmabar)) +
                      m_B * (8.0 * m_v * (m_B2 * sigmabar2 - q2) * sigma +
                      m_B * sigmabar2 * (-q2 - m_B2 * (-1 + 3.0 * sigma) * sigmabar) +
                      m_B * sigma * sigmabar * (-(3.0 * q2) - m_B * sigmabar * (4.0 * m_v - 3.0 * m_B * sigmabar)))) * phi_bar_bar_4 +
                      m_B * (m_B2 * sigmabar2 - q2) * (-omega_1 + m_B * sigma) * sigmabar * (2.0 * m_v + m_B * sigmabar) * phi_bar_4));
        }

        double I4d2C_fp_3pt_phi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * m_B * (m_B * sigmabar - 2.0 * m_v) * (m_B2 * sigmabar2 - q2)
                             / (power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d2D_fp_3pt_phi_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to phi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_4     = this->phi_bar_4(omega_1, omega_2);
            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            return    -(6.0 * pow(sigmabar,-4.0) * ((q2 - m_B * sigmabar * (m_B + 4.0 * m_v - m_B * sigma + 2.0 * m_B * sigmabar)) *
                        phi_bar_bar_4 + (m_B2 * sigmabar2 - q2) * (2.0 * m_v + m_B * sigmabar) * phi_bar_4));
        }

        double I2_fp_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_4
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_2 = (1.0 - 2.0 * u) / (m_B * power_of<2>(sigmabar));

            return C_2 * psi_bar_4;
        }

        double I3_fp_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (2.0 * u - 1.0) * (-m_B2 * sigmabar2 + m_v2 + q2)
                             / (m_B * power_of<3>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I3d1A_fp_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = -(2.0 * (6.0 * m_B3 * sigma3 + 3.0 * (2.0 * omega_1 + omega_2) * (m_v2 + q2) + 2.0 * m_B3 * sigmabar -
                                 2.0 * m_B * (m_v2 + q2) * sigmabar + m_B2 * (2.0 * omega_1 + omega_2) * (-3.0 + 2.0 * sigmabar) +
                                 3.0 * m_B2 * sigma2 * (-(2.0 * omega_1) - omega_2 + 2.0 * m_B * (-2.0 + sigmabar)) -
                                 2.0 * m_B * sigma * (3.0 * (m_v2 + q2) + m_B * (2.0 * omega_1 + omega_2) * (-3.0 + sigmabar) +
                                 m_B2 * (-3.0 + 4.0 * sigmabar))))
                             / (m_B * omega_2 * power_of<4>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I3d1B_fp_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (m_B2 * sigmabar2 - m_v2 - q2) / ((-omega_1 + m_B * sigma) * power_of<3>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I3d1C_fp_3pt_psi_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (m_B2 * sigmabar2 - m_v2 - q2) / (omega_2 * power_of<3>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I4_fp_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * m_v * u * (2.0 * m_B * sigmabar + m_v * (2.0 * u - 1.0))
                             / (power_of<3>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1A_fp_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * (m_B * ((-(2.0 * m_B) + m_v) * omega_2 + 6.0 * m_B * (m_v - omega_2) * sigma2 -
                               2.0 * (m_v * omega_2 + 2.0 * m_B * (m_v - 2.0 * omega_2)) * sigma) * sigmabar -
                               4.0 * m_B * sigma * sigmabar * (2.0 * m_B * m_v * sigma - omega_2 * (m_v - 2.0 * m_B * sigmabar)) +
                               omega_1 * (4.0 * sigmabar * (4.0 * m_B * m_v * sigma - omega_2 * (m_v - 2.0 * m_B * sigmabar)) +
                               sigmabar * (m_v * omega_2 + 4.0 * m_B * (m_v - 2.0 * m_v * sigma - omega_2 * sigmabar))) +
                               2.0 * m_v * (-4.0 + 4.0 * sigma + sigmabar) * power_of<2>(omega_1))
                             / (power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1B_fp_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * (m_B - m_B * sigma) * (2.0 * m_B * sigmabar + m_v)
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * psi_bar_bar_4;
        }

        double I4d1C_fp_3pt_psiA_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4

            return 0.0;
        }
        double I4d2A_fp_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 12.0 * m_v * (m_B *
                              (10.0 * (-1 + sigma) * sigma * (-(omega_2 * (m_v + 2.0 * m_B * (-1 + sigma))) + 2.0 * m_B * m_v * sigma) +
                               sigmabar2 * (-(m_v * omega_2) + m_B * (4.0 * omega_2 - 6.0 * omega_2 * sigma + m_v * (-2.0 + 6.0 * sigma))) +
                               4.0 * ((-(2.0 * m_B) + m_v) * omega_2 + 6.0 * m_B * (m_v - omega_2) * sigma2 -
                               2.0 * (m_v * omega_2 + 2.0 * m_B * (m_v - 2.0 * omega_2)) * sigma) * sigmabar) +
                               2.0 * omega_1 * (m_B * (-(2.0 * m_v) + omega_2) * sigmabar2 -
                               5.0 * (-1 + sigma) * (-(omega_2 * (m_v + 2.0 * m_B * (-1 + sigma))) + 4.0 * m_B * m_v * sigma) +
                               2.0 * (m_v * omega_2 + 4.0 * m_B * (m_v + omega_2 * (-1 + sigma) - 2.0 * m_v * sigma)) * sigmabar) +
                               4.0 * m_v * (-5.0 + 5.0 * sigma + 2.0 * sigmabar) * power_of<2>(omega_1))
                             / (power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2B_fp_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return  6.0 * m_B * m_v *
                   (2.0 * (4.0 * m_B * (-m_v + 2.0 * m_B * (-1 + sigma)) * (-1 + sigma) * sigma +
                    m_B * (m_v - 2.0 * m_v * sigma + 4.0 * m_B * (-1 + sigma) * sigma) * sigmabar +
                    omega_1 * (-(4.0 * m_B * (-1 + sigma) * (-2.0 + 2.0 * sigma + sigmabar)) + m_v * (-4.0 + 4.0 * sigma + sigmabar))) *
                    psi_bar_bar_4 + m_B * (-m_v + 2.0 * m_B * (-1 + sigma)) * (-1 + sigma) * (-omega_1 + m_B * sigma) * sigmabar * psi_bar_4)
                    / (power_of<2>(omega_2) * power_of<5>(sigmabar));
        }

        double I4d2C_fp_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * m_B2 * m_v * (m_v - 2.0 * m_B * sigmabar) / (power_of<2>(omega_2) * power_of<3>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2D_fp_3pt_psiA_bar_bar_4(const double & sigma, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return     (6.0 * m_v * ((-4.0 * m_B * sigmabar - m_v) * psi_bar_bar_4
                      - m_B * sigmabar * (-2.0 * m_B * sigmabar - m_v) * psi_bar_4)) / (power_of<4>(sigmabar));
        }

        double I3_fp_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * u * (m_B * sigmabar * (2.0 * u - 1.0) + 2.0 * m_v)
                             / (m_B * power_of<3>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1A_fp_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B),   m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2   = power_of<2>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (4.0 * m_B3 * sigma2 * sigmabar - m_B2 * sigmabar2 * (omega_2 + 4.0 * omega_1) +
                               2.0 * m_B * sigma * (m_B * sigmabar *  (2.0 * m_B * sigmabar - omega_2 - 4.0 * omega_1) + 3.0 * m_v * omega_2) +
                               2.0 * m_B * sigmabar * (m_v * omega_2 + omega_1 * (omega_2 + 2.0 * omega_1)) - 6.0 * m_v * omega_2 * omega_1)
                             / (m_B * power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1B_fp_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (m_B * sigmabar + 2.0 * m_v)
                             / ((-omega_1 + m_B * sigma) * power_of<3>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1C_fp_3pt_psiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4

            return 0.0;
        }

        double I4_fp_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * u * (m_B2 * sigmabar2 - q2) * (m_B * sigmabar * (2.0 * u - 1.0) + 2.0 * m_v)
                             / (m_B * power_of<4>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1A_fp_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * (-(m_B4 * (4.0 * omega_1 + omega_2) * sigmabar2) + 8.0 * m_v * omega_1 * omega_2 * q2 +
                               6.0 * m_B5 * sigma4 * sigmabar - m_B * (2.0 * m_v * omega_2 + 3.0 * omega_1 * (2.0 * omega_1 + omega_2)) * q2 *
                               sigmabar + m_B2 * ((4.0 * omega_1 + omega_2) * sigmabar2 * q2 +
                               4.0 * m_v * omega_1 * omega_2 * (-2.0 + sigmabar)) +
                               m_B3 * sigmabar * (2.0 * m_v * omega_2 - omega_1 * (2.0 * omega_1 + omega_2) * (-3.0 + 2.0 * sigmabar)) +
                               m_B3 * sigma3 * (8.0 * m_v * omega_2 + m_B * sigmabar *
                               (-(3.0 * (4.0 * omega_1 + omega_2)) + 4.0 * m_B * (-3.0 + 2.0 * sigmabar))) +
                               m_B2 * sigma2 * (-(8.0 * m_v * omega_1 * omega_2) -
                               3.0 * m_B2 * (4.0 * omega_1 + omega_2) * sigmabar * (-2.0 + sigmabar) +
                               6.0 * m_B3 * sigmabar * (1 - 2.0 * sigmabar) +
                               m_B * (3.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar - 6.0 * q2 * sigmabar +
                               2.0 * m_v * omega_2 * (-8.0 + 3.0 * sigmabar))) +
                               m_B * sigma * (4.0 * m_B4 * sigmabar2 - 8.0 * m_v * omega_2 * q2 +
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar * (-3.0 + 4.0 * sigmabar) +
                               m_B * (3.0 * (4.0 * omega_1 + omega_2) * q2 * sigmabar - 4.0 * m_v * omega_1 * omega_2 * (-4.0 + sigmabar)) -
                               2.0 * m_B2 * (4.0 * m_v * omega_2 * (-1 + sigmabar) +
                               sigmabar * (2.0 * q2 * sigmabar - omega_1 * (2.0 * omega_1 + omega_2) * (-3.0 + sigmabar)))))
                             / (m_B * power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1B_fp_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * (m_B * sigmabar + 2.0 * m_v) * (m_B2 * sigmabar2 -q2)
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * psi_bar_bar_4;
        }

        double I4d1C_fp_3pt_psiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4

            return 0.0;
        }

        double I4d2A_fp_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 =12.0 * (2.0 * m_B5 * sigmabar3 + 20.0 * m_v * omega_1 * omega_2 * q2 + 12.0 * m_B5 * sigma4 * sigmabar -
                              2.0 * m_B * (4.0 * m_v * omega_2 + 3.0 * omega_1 * (2.0 * omega_1 + omega_2)) * q2 * sigmabar +
                              m_B4 * (4.0 * omega_1 + omega_2) * sigmabar2 * (-3.0 + 2.0 * sigmabar) +
                              2.0 * m_B3 * sigma3 * (10.0 * m_v * omega_2 +
                              3.0 * m_B * (-(4.0 * omega_1) - omega_2 + 4.0 * m_B * (-1 + sigmabar)) * sigmabar) +
                              m_B3 * sigmabar * (-(2.0 * sigmabar2 * q2) - 4.0 * m_v * omega_2 * (-2.0 + sigmabar) +
                              omega_1 * (2.0 * omega_1 + omega_2) * (6.0 + sigmabar * (-6.0 + sigmabar))) +
                              m_B2 * (3.0 * (4.0 * omega_1 + omega_2) * sigmabar2 * q2 -
                              2.0 * m_v * omega_1 * omega_2 * (10.0 + sigmabar * (-8.0 + sigmabar))) +
                              m_B2 * sigma2 * (-(20.0 * m_v * omega_1 * omega_2) -
                              3.0 * m_B2 * (4.0 * omega_1 + omega_2) * sigmabar * (-4.0 + 3.0 * sigmabar) +
                              12.0 * m_B3 * sigmabar * (1 + sigmabar * (-3.0 + sigmabar)) +
                              2.0 * m_B * (3.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar - 6.0 * q2 * sigmabar +
                              4.0 * m_v * omega_2 * (-5.0 + 3.0 * sigmabar))) +
                              m_B * sigma * (-(20.0 * m_v * omega_2 * q2) - 12.0 * m_B4 * sigmabar2 * (-1 + sigmabar) +
                              6.0 * m_B * omega_2 * q2 * sigmabar - 3.0 * m_B3 * (4.0 * omega_1 + omega_2) * sigmabar *
                             (2.0 + sigmabar * (-4.0 + sigmabar)) +
                              8.0 * m_B * omega_1 * (3.0 * q2 * sigmabar + m_v * omega_2 * (5.0 - 2.0 * sigmabar)) +
                              2.0 * m_B2 * (3.0 * sigmabar * (-(2.0 * q2 * sigmabar) + omega_1 * (2.0 * omega_1 + omega_2) * (-2.0 + sigmabar)) +
                              m_v * omega_2 * (10.0 + sigmabar * (-16.0 + 3.0 * sigmabar)))))
                             / (m_B * power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2B_fp_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return  -(6.0 * pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-5.0) *
                     (2.0 * (omega_1 * (8.0 * m_v * q2 + 3.0 * m_B * q2 * sigmabar +
                      4.0 * m_B2 * m_v * sigmabar * (-2.0 + 2.0 * sigma + sigmabar) +
                      m_B3 * sigmabar2 * (-3.0 + 3.0 * sigma + 2.0 * sigmabar)) +
                      m_B * (8.0 * m_v * (m_B2 * sigmabar2 - q2) * sigma +
                      m_B * sigmabar2 * (-q2 - m_B2 * (-1 + 3.0 * sigma) * sigmabar) +
                      m_B * sigma * sigmabar * (-(3.0 * q2) - m_B * sigmabar * (4.0 * m_v - 3.0 * m_B * sigmabar)))) * psi_bar_bar_4 +
                      m_B * (m_B2 * sigmabar2 - q2) * (-omega_1 + m_B * sigma) * sigmabar * (2.0 * m_v + m_B * sigmabar) * psi_bar_4));
        }

        double I4d2C_fp_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * m_B * (m_B * sigmabar - 2.0 * m_v) * (m_B2 * sigmabar2 - q2)
                             / (power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2D_fp_3pt_psiB_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return    -(6.0 * pow(sigmabar,-4.0) * ((q2 - m_B * sigmabar * (m_B + 4.0 * m_v - m_B * sigma + 2.0 * m_B * sigmabar)) *
                        psi_bar_bar_4 + (m_B2 * sigmabar2 - q2) * (2.0 * m_v + m_B * sigmabar) * psi_bar_4));
        }

        double I3_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - 0.0                                                       - I3_fp_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1A_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - 0.0                                                       - I3d1A_fp_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1B_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return - 0.0                                                       - I3d1B_fp_3pt_psiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I3d1C_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return - 0.0                                                       - I3d1C_fp_3pt_psiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I4_fp_3pt_psiA_bar_bar_4( sigma, omega_1, omega_2, q2)    - I4_fp_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1A_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I4d1A_fp_3pt_psiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d1A_fp_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1B_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return - I4d1B_fp_3pt_psiA_bar_bar_4( sigma, omega_1, q2)          - I4d1B_fp_3pt_psiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d1C_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return - I4d1C_fp_3pt_psiA_bar_bar_4( sigma, omega_2, q2)          - I4d1C_fp_3pt_psiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2A_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I4d2A_fp_3pt_psiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d2A_fp_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d2B_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return - I4d2B_fp_3pt_psiA_bar_bar_4( sigma, omega_1, q2)          - I4d2B_fp_3pt_psiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d2C_fp_3pt_psi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return - I4d2C_fp_3pt_psiA_bar_bar_4( sigma, omega_2, q2)          - I4d2C_fp_3pt_psiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2D_fp_3pt_psi_bar_bar_4(const double & sigma, const double & q2) const
        {
            return - I4d2D_fp_3pt_psiA_bar_bar_4( sigma, q2)                   - I4d2D_fp_3pt_psiB_bar_bar_4(sigma, q2);
        }

        double I2_fp_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_4
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_2 = 1.0 / (m_B * power_of<2>(sigmabar));

            return C_2 * chi_bar_4;
        }

        double I3_fp_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (m_B2 * sigmabar2 * (2.0 * u - 1.0) + 4.0 * m_B * m_v * sigmabar + m_v2 + q2)
                             / (m_B * power_of<3>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I3d1A_fp_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (2.0 * omega_2 * sigma2 * q2 - 2.0 * m_B3 * sigmabar3 * sigma + 2.0 * omega_2 * q2 * sigma * (-1 + sigmabar) +
                               sigmabar * (-(3.0 * m_v2 * omega_2) - 3.0 * omega_2 * q2 - 8.0 * m_B * m_v * omega_2 * sigmabar +
                               m_B2 * sigmabar2 * (2.0 * omega_1 + omega_2 - 2.0 * m_B * sigmabar)))
                             / (m_B * omega_2 * power_of<5>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I3d1B_fp_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = (2.0 * sigmabar * (m_B2 * sigmabar2 + 4.0 * m_B * m_v * sigmabar + m_v2 + q2) - q2 * sigma2 + sigma * (q2 - q2 * sigmabar))
                             / ((-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I3d1C_fp_3pt_chi_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = (-2.0 * sigmabar * (-m_B2 * sigmabar2 + 4.0 * m_B * m_v * sigmabar + m_v2 + q2) + q2 * sigma * sigma + q2 * sigma * (sigmabar - 1.0))
                             / (omega_2 * power_of<4>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I4_fp_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * m_v * u * (2.0 * m_B * sigmabar + m_v * (2.0 * u - 1.0))
                             / (power_of<3>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1A_fp_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * (m_B * ((-(2.0 * m_B) + m_v) * omega_2 + 6.0 * m_B * (m_v - omega_2) * sigma2 -
                               2.0 * (m_v * omega_2 + 2.0 * m_B * (m_v - 2.0 * omega_2)) * sigma) * sigmabar -
                               4.0 * m_B * sigma * sigmabar * (2.0 * m_B * m_v * sigma - omega_2 * (m_v - 2.0 * m_B * sigmabar)) +
                               omega_1 * (4.0 * sigmabar * (4.0 * m_B * m_v * sigma - omega_2 * (m_v - 2.0 * m_B * sigmabar)) +
                               sigmabar * (m_v * omega_2 + 4.0 * m_B * (m_v - 2.0 * m_v * sigma - omega_2 * sigmabar))) +
                               2.0 * m_v * (-4.0 + 4.0 * sigma + sigmabar) * power_of<2>(omega_1))
                             / (power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1B_fp_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * (m_B - m_B * sigma) * (2.0 * m_B * sigmabar + m_v)
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * chi_bar_bar_4;
        }

        double I4d1C_fp_3pt_chiA_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4

            return 0.0;
        }
        double I4d2A_fp_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 12.0 * m_v * (m_B *
                              (10.0 * (-1 + sigma) * sigma * (-(omega_2 * (m_v + 2.0 * m_B * (-1 + sigma))) + 2.0 * m_B * m_v * sigma) +
                               sigmabar2 * (-(m_v * omega_2) + m_B * (4.0 * omega_2 - 6.0 * omega_2 * sigma + m_v * (-2.0 + 6.0 * sigma))) +
                               4.0 * ((-(2.0 * m_B) + m_v) * omega_2 + 6.0 * m_B * (m_v - omega_2) * sigma2 -
                               2.0 * (m_v * omega_2 + 2.0 * m_B * (m_v - 2.0 * omega_2)) * sigma) * sigmabar) +
                               2.0 * omega_1 * (m_B * (-(2.0 * m_v) + omega_2) * sigmabar2 -
                               5.0 * (-1 + sigma) * (-(omega_2 * (m_v + 2.0 * m_B * (-1 + sigma))) + 4.0 * m_B * m_v * sigma) +
                               2.0 * (m_v * omega_2 + 4.0 * m_B * (m_v + omega_2 * (-1 + sigma) - 2.0 * m_v * sigma)) * sigmabar) +
                               4.0 * m_v * (-5.0 + 5.0 * sigma + 2.0 * sigmabar) * power_of<2>(omega_1))
                             / (power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2B_fp_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return  6.0 * m_B * m_v *
                   (2.0 * (4.0 * m_B * (-m_v + 2.0 * m_B * (-1 + sigma)) * (-1 + sigma) * sigma +
                    m_B * (m_v - 2.0 * m_v * sigma + 4.0 * m_B * (-1 + sigma) * sigma) * sigmabar +
                    omega_1 * (-(4.0 * m_B * (-1 + sigma) * (-2.0 + 2.0 * sigma + sigmabar)) + m_v * (-4.0 + 4.0 * sigma + sigmabar))) *
                    chi_bar_bar_4 + m_B * (-m_v + 2.0 * m_B * (-1 + sigma)) * (-1 + sigma) * (-omega_1 + m_B * sigma) * sigmabar * chi_bar_4)
                    / (power_of<2>(omega_2) * power_of<5>(sigmabar));
        }

        double I4d2C_fp_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * m_B2 * m_v * (m_v - 2.0 * m_B * sigmabar) / (power_of<2>(omega_2) * power_of<3>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2D_fp_3pt_chiA_bar_bar_4(const double & sigma, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return     (6.0 * m_v * ((-4.0 * m_B * sigmabar - m_v) * chi_bar_bar_4
                      - m_B * sigmabar * (-2.0 * m_B * sigmabar - m_v) * chi_bar_4)) / (power_of<4>(sigmabar));
        }

        double I3_fp_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * u * (m_B * sigmabar * (2.0 * u - 1.0) + 2.0 * m_v)
                             / (m_B * power_of<3>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1A_fp_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B),   m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2   = power_of<2>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (4.0 * m_B3 * sigma2 * sigmabar - m_B2 * sigmabar2 * (omega_2 + 4.0 * omega_1) +
                               2.0 * m_B * sigma * (m_B * sigmabar *  (2.0 * m_B * sigmabar - omega_2 - 4.0 * omega_1) + 3.0 * m_v * omega_2) +
                               2.0 * m_B * sigmabar * (m_v * omega_2 + omega_1 * (omega_2 + 2.0 * omega_1)) - 6.0 * m_v * omega_2 * omega_1)
                             / (m_B * power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1B_fp_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (m_B * sigmabar + 2.0 * m_v)
                             / ((-omega_1 + m_B * sigma) * power_of<3>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1C_fp_3pt_chiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4

            return 0.0;
        }

        double I4_fp_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * u * (m_B2 * sigmabar2 - q2) * (m_B * sigmabar * (2.0 * u - 1.0) + 2.0 * m_v)
                             / (m_B * power_of<4>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1A_fp_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * (-(m_B4 * (4.0 * omega_1 + omega_2) * sigmabar2) + 8.0 * m_v * omega_1 * omega_2 * q2 +
                               6.0 * m_B5 * sigma4 * sigmabar - m_B * (2.0 * m_v * omega_2 + 3.0 * omega_1 * (2.0 * omega_1 + omega_2)) * q2 *
                               sigmabar + m_B2 * ((4.0 * omega_1 + omega_2) * sigmabar2 * q2 +
                               4.0 * m_v * omega_1 * omega_2 * (-2.0 + sigmabar)) +
                               m_B3 * sigmabar * (2.0 * m_v * omega_2 - omega_1 * (2.0 * omega_1 + omega_2) * (-3.0 + 2.0 * sigmabar)) +
                               m_B3 * sigma3 * (8.0 * m_v * omega_2 + m_B * sigmabar *
                               (-(3.0 * (4.0 * omega_1 + omega_2)) + 4.0 * m_B * (-3.0 + 2.0 * sigmabar))) +
                               m_B2 * sigma2 * (-(8.0 * m_v * omega_1 * omega_2) -
                               3.0 * m_B2 * (4.0 * omega_1 + omega_2) * sigmabar * (-2.0 + sigmabar) +
                               6.0 * m_B3 * sigmabar * (1 - 2.0 * sigmabar) +
                               m_B * (3.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar - 6.0 * q2 * sigmabar +
                               2.0 * m_v * omega_2 * (-8.0 + 3.0 * sigmabar))) +
                               m_B * sigma * (4.0 * m_B4 * sigmabar2 - 8.0 * m_v * omega_2 * q2 +
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar * (-3.0 + 4.0 * sigmabar) +
                               m_B * (3.0 * (4.0 * omega_1 + omega_2) * q2 * sigmabar - 4.0 * m_v * omega_1 * omega_2 * (-4.0 + sigmabar)) -
                               2.0 * m_B2 * (4.0 * m_v * omega_2 * (-1 + sigmabar) +
                               sigmabar * (2.0 * q2 * sigmabar - omega_1 * (2.0 * omega_1 + omega_2) * (-3.0 + sigmabar)))))
                             / (m_B * power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1B_fp_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * (m_B * sigmabar + 2.0 * m_v) * (m_B2 * sigmabar2 -q2)
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * chi_bar_bar_4;
        }

        double I4d1C_fp_3pt_chiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4

            return 0.0;
        }

        double I4d2A_fp_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 =12.0 * (2.0 * m_B5 * sigmabar3 + 20.0 * m_v * omega_1 * omega_2 * q2 + 12.0 * m_B5 * sigma4 * sigmabar -
                              2.0 * m_B * (4.0 * m_v * omega_2 + 3.0 * omega_1 * (2.0 * omega_1 + omega_2)) * q2 * sigmabar +
                              m_B4 * (4.0 * omega_1 + omega_2) * sigmabar2 * (-3.0 + 2.0 * sigmabar) +
                              2.0 * m_B3 * sigma3 * (10.0 * m_v * omega_2 +
                              3.0 * m_B * (-(4.0 * omega_1) - omega_2 + 4.0 * m_B * (-1 + sigmabar)) * sigmabar) +
                              m_B3 * sigmabar * (-(2.0 * sigmabar2 * q2) - 4.0 * m_v * omega_2 * (-2.0 + sigmabar) +
                              omega_1 * (2.0 * omega_1 + omega_2) * (6.0 + sigmabar * (-6.0 + sigmabar))) +
                              m_B2 * (3.0 * (4.0 * omega_1 + omega_2) * sigmabar2 * q2 -
                              2.0 * m_v * omega_1 * omega_2 * (10.0 + sigmabar * (-8.0 + sigmabar))) +
                              m_B2 * sigma2 * (-(20.0 * m_v * omega_1 * omega_2) -
                              3.0 * m_B2 * (4.0 * omega_1 + omega_2) * sigmabar * (-4.0 + 3.0 * sigmabar) +
                              12.0 * m_B3 * sigmabar * (1 + sigmabar * (-3.0 + sigmabar)) +
                              2.0 * m_B * (3.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar - 6.0 * q2 * sigmabar +
                              4.0 * m_v * omega_2 * (-5.0 + 3.0 * sigmabar))) +
                              m_B * sigma * (-(20.0 * m_v * omega_2 * q2) - 12.0 * m_B4 * sigmabar2 * (-1 + sigmabar) +
                              6.0 * m_B * omega_2 * q2 * sigmabar - 3.0 * m_B3 * (4.0 * omega_1 + omega_2) * sigmabar *
                             (2.0 + sigmabar * (-4.0 + sigmabar)) +
                              8.0 * m_B * omega_1 * (3.0 * q2 * sigmabar + m_v * omega_2 * (5.0 - 2.0 * sigmabar)) +
                              2.0 * m_B2 * (3.0 * sigmabar * (-(2.0 * q2 * sigmabar) + omega_1 * (2.0 * omega_1 + omega_2) * (-2.0 + sigmabar)) +
                              m_v * omega_2 * (10.0 + sigmabar * (-16.0 + 3.0 * sigmabar)))))
                             / (m_B * power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2B_fp_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return  -(6.0 * pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-5.0) *
                     (2.0 * (omega_1 * (8.0 * m_v * q2 + 3.0 * m_B * q2 * sigmabar +
                      4.0 * m_B2 * m_v * sigmabar * (-2.0 + 2.0 * sigma + sigmabar) +
                      m_B3 * sigmabar2 * (-3.0 + 3.0 * sigma + 2.0 * sigmabar)) +
                      m_B * (8.0 * m_v * (m_B2 * sigmabar2 - q2) * sigma +
                      m_B * sigmabar2 * (-q2 - m_B2 * (-1 + 3.0 * sigma) * sigmabar) +
                      m_B * sigma * sigmabar * (-(3.0 * q2) - m_B * sigmabar * (4.0 * m_v - 3.0 * m_B * sigmabar)))) * chi_bar_bar_4 +
                      m_B * (m_B2 * sigmabar2 - q2) * (-omega_1 + m_B * sigma) * sigmabar * (2.0 * m_v + m_B * sigmabar) * chi_bar_4));
        }

        double I4d2C_fp_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * m_B * (m_B * sigmabar - 2.0 * m_v) * (m_B2 * sigmabar2 - q2)
                             / (power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2D_fp_3pt_chiB_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_+ proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return    -(6.0 * pow(sigmabar,-4.0) * ((q2 - m_B * sigmabar * (m_B + 4.0 * m_v - m_B * sigma + 2.0 * m_B * sigmabar)) *
                        chi_bar_bar_4 + (m_B2 * sigmabar2 - q2) * (2.0 * m_v + m_B * sigmabar) * chi_bar_4));
        }

        double I3_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + 0.0                                                       - I3_fp_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1A_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + 0.0                                                       - I3d1A_fp_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1B_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return + 0.0                                                       - I3d1B_fp_3pt_chiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I3d1C_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return + 0.0                                                       - I3d1C_fp_3pt_chiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I4_fp_3pt_chiA_bar_bar_4( sigma, omega_1, omega_2, q2)    - I4_fp_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1A_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I4d1A_fp_3pt_chiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d1A_fp_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1B_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return + I4d1B_fp_3pt_chiA_bar_bar_4( sigma, omega_1, q2)          - I4d1B_fp_3pt_chiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d1C_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return + I4d1C_fp_3pt_chiA_bar_bar_4( sigma, omega_2, q2)          - I4d1C_fp_3pt_chiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2A_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I4d2A_fp_3pt_chiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d2A_fp_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d2B_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return + I4d2B_fp_3pt_chiA_bar_bar_4( sigma, omega_1, q2)          - I4d2B_fp_3pt_chiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d2C_fp_3pt_chi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return + I4d2C_fp_3pt_chiA_bar_bar_4( sigma, omega_2, q2)          - I4d2C_fp_3pt_chiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2D_fp_3pt_chi_bar_bar_4(const double & sigma, const double & q2) const
        {
            return + I4d2D_fp_3pt_chiA_bar_bar_4( sigma, q2)                   - I4d2D_fp_3pt_chiB_bar_bar_4(sigma, q2);
        }
        // }}}

        /* f_+ : integrands and surface terms */
        // {{{
        double integrand_fp_2pt_disp(const double & sigma, const double & q2) const
        {
            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);
            const double etad3    = 24.0 * (eta - 1.0) * power_of<2>(eta) * (2.0 * eta - 1.0) / power_of<3>(sigmabar);

            const double I1   = I1_fp_2pt_phi_p(sigma, q2);
            const double I2   = I2_fp_2pt_phi_bar(sigma, q2)   + I2_fp_2pt_g_p(sigma, q2);
            const double I2d1 = I2d1_fp_2pt_phi_bar(sigma, q2) + I2d1_fp_2pt_g_p(sigma, q2);
            const double I3   = I3_fp_2pt_g_p(sigma, q2)       + I3_fp_2pt_g_bar(sigma, q2);
            const double I3d1 = I3d1_fp_2pt_g_p(sigma, q2)     + I3d1_fp_2pt_g_bar(sigma, q2);
            const double I3d2 = I3d2_fp_2pt_g_p(sigma, q2)     + I3d2_fp_2pt_g_bar(sigma, q2);
            const double I4   = I4_fp_2pt_g_bar(sigma, q2);
            const double I4d1 = I4d1_fp_2pt_g_bar(sigma, q2);
            const double I4d2 = I4d2_fp_2pt_g_bar(sigma, q2);
            const double I4d3 = I4d3_fp_2pt_g_bar(sigma, q2);

            double result = 0.0;
            result += -1.0 * I1;
            result += (etad1 * I2 + eta * I2d1) / m_B2;
            result += -1.0 * (I3 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I3d1 * eta * etad1 + I3d2 * power_of<2>(eta)) / (2.0 * m_B4);
            result += I4 * (power_of<2>(eta) * etad3 + 4.0 * eta * etad1 * etad2 + power_of<3>(etad1)) / (6.0 * m_B6);
            result += I4d1 * eta * (4.0 * eta * etad2 + 7.0 * power_of<2>(etad1)) / (6.0 * m_B6);
            result += I4d2 * 6.0 * power_of<2>(eta) * etad1 / (6.0 * m_B6);
            result += I4d3 * power_of<3>(eta) / (6.0 * m_B6);
            result *= exp;
            return result;
        }

        double integrand_fp_2pt_borel(const double & sigma, const double & q2) const
        {
            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double I1   = I1_fp_2pt_phi_p(sigma, q2);
            const double I2   = I2_fp_2pt_phi_bar(sigma, q2)   + I2_fp_2pt_g_p(sigma, q2);
            const double I3   = I3_fp_2pt_g_p(sigma, q2)       + I3_fp_2pt_g_bar(sigma, q2);
            const double I4   = I4_fp_2pt_g_bar(sigma, q2);

            double result = 0.0;
            result += - I1;
            result +=   I2 / M2;
            result += - I3 / (2 * M4);
            result +=   I4 / (6 * M6);
            result *= exp;

            return result;
        }

        double surface_fp_2pt(const double & sigma, const double & q2) const
        {
            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fp_2pt_phi_bar(sigma, q2)   + I2_fp_2pt_g_p(sigma, q2);
            const double I3   = I3_fp_2pt_g_p(sigma, q2)       + I3_fp_2pt_g_bar(sigma, q2);
            const double I3d1 = I3d1_fp_2pt_g_p(sigma, q2)     + I3d1_fp_2pt_g_bar(sigma, q2);
            const double I4   = I4_fp_2pt_g_bar(sigma, q2);
            const double I4d1 = I4d1_fp_2pt_g_bar(sigma, q2);
            const double I4d2 = I4d2_fp_2pt_g_bar(sigma, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / (power_of<2>( M2())));
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;

            return result;
        }

        /*
         * rewrite integration ranges such that:
         * 1.)
         *    0 <= x_1 <= 1,     and    0 <= x_2 <= 1,
         * 2.)
         *    x_1 and x_2 integration boundaries do not depend on the other variables
         *
         * We obtain the integrand
         *
         *    sigma m_B f(sigma m_B x_1, sigma m_B (xbar_1 xbar_2 + x_2) / xbar_2) / (xbar_1 xbar_2^2 + x_2 xbar_2),
         *
         * where
         *
         *    xbar_1 = 1.0 - x_1,    and    xbar_2 = 1.0 - x_2.
         */
        double integrand_fp_3pt(const std::array<double, 3> & args, const double & q2) const
        {
            const double sigma  = args[0];
            const double x_1    = args[1];
            const double x_2    = args[2];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            constexpr double I1 = 0.0;
            const     double I2 = I2_fp_3pt_phi_3(sigma, omega_1, omega_2, q2)         + I2_fp_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                                + I2_fp_3pt_phi_4(sigma, omega_1, omega_2, q2)         + I2_fp_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)
                                + I2_fp_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I2_fp_3pt_chi_bar_4(sigma, omega_1, omega_2, q2);
            const     double I3 = I3_fp_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                                + I3_fp_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                                + I3_fp_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fp_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                                + I3_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const     double I4 = I4_fp_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                                + I4_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result = 0.0;
            result += - I1;
            result +=   I2 / M2;
            result += - I3 / (2.0 * M4);
            result +=   I4 / (6.0 * M6);
            result *=   prefactor * exp;

            return result;
        }

        double surface_fp_3pt_A(const std::array<double, 2> & args, const double & sigma, const double & q2) const
        {
            const double x_1    = args[0];
            const double x_2    = args[1];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fp_3pt_phi_3(sigma, omega_1, omega_2, q2)            + I2_fp_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I2_fp_3pt_phi_4(sigma, omega_1, omega_2, q2)            + I2_fp_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)
                              + I2_fp_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I2_fp_3pt_chi_bar_4(sigma, omega_1, omega_2, q2);
            const double I3   = I3_fp_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I3_fp_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fp_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fp_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I3_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I3d1 = I3d1A_fp_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I3d1A_fp_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fp_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fp_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3d1A_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4   = I4_fp_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)    + I4_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I4_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d1 = I4d1A_fp_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d1A_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d1A_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d1A_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d2 = I4d2A_fp_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d2A_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d2A_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d2A_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        double surface_fp_3pt_B(const double & x_1, const double & sigma, const double & q2) const
        {
            // this ONLY includes the Jacobian from the transformation (omega_1 -> x_1).
            const double prefactor = sigma * m_B();

            const double omega_1 = sigma * m_B() * x_1;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1B_fp_3pt_phi_bar_3(sigma, omega_1, q2)
                                  + I3d1B_fp_3pt_phi_bar_4(sigma, omega_1, q2)     + I3d1B_fp_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fp_3pt_psi_bar_4(sigma, omega_1, q2)     + I3d1B_fp_3pt_chi_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fp_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I3d1B_fp_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1B_fp_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d1B_fp_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d1B_fp_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d1B_fp_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            const     double I4d2 = I4d2B_fp_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d2B_fp_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d2B_fp_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d2B_fp_3pt_chi_bar_bar_4(sigma, omega_1, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        double surface_fp_3pt_C(const double & x_2, const double & sigma, const double & q2) const
        {
            const double xbar_2 = 1.0 - x_2;

            // this ONLY includes the Jacobian from the transformation (omega_2 -> x_2).
            const double prefactor = sigma * m_B() / (xbar_2 * xbar_2);

            const double omega_2 = sigma * m_B() * (x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1C_fp_3pt_phi_bar_3(sigma, omega_2, q2)
                                  + I3d1C_fp_3pt_phi_bar_4(sigma, omega_2, q2)     + I3d1C_fp_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fp_3pt_psi_bar_4(sigma, omega_2, q2)     + I3d1C_fp_3pt_chi_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fp_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I3d1C_fp_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1C_fp_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d1C_fp_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d1C_fp_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d1C_fp_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            const     double I4d2 = I4d2C_fp_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d2C_fp_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d2C_fp_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d2C_fp_3pt_chi_bar_bar_4(sigma, omega_2, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        double surface_fp_3pt_D(const double & sigma, const double & q2) const
        {
            // this does NOT includes the original factor of 1 / omega_2
            const double prefactor = 1.0;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);


            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            constexpr double I3d1 = 0.0;
            constexpr double I4   = 0.0;
            constexpr double I4d1 = 0.0;
            const     double I4d2 = I4d2D_fp_3pt_phi_bar_bar_3(sigma, q2) + I4d2D_fp_3pt_phi_bar_bar_4(sigma, q2)
                                  + I4d2D_fp_3pt_psi_bar_bar_4(sigma, q2) + I4d2D_fp_3pt_chi_bar_bar_4(sigma, q2);


            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        /*
         * Integrands for the first moments. Only the borel method is implemented
         */

        double integrand_fp_2pt_borel_m1(const double & sigma, const double & q2) const
        {
            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double I1   = I1_fp_2pt_phi_p(sigma, q2);
            const double I2   = I2_fp_2pt_phi_bar(sigma, q2)   + I2_fp_2pt_g_p(sigma, q2);
            const double I3   = I3_fp_2pt_g_p(sigma, q2)       + I3_fp_2pt_g_bar(sigma, q2);
            const double I4   = I4_fp_2pt_g_bar(sigma, q2);

            double result1 = 0.0;
            result1 += - I1;
            result1 +=   I2 / M2;
            result1 += - I3 / (2.0 * M4);
            result1 +=   I4 / (6.0 * M6);
            result1 *=   exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - I2;
            result2 +=   I3 / M2;
            result2 += - I4 / (2.0 * M4);
            result2 *=   exp;

            return result1 + result2;
        }

        double surface_fp_2pt_m1(const double & sigma, const double & q2) const
        {
            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double m_P2 = power_of<2>(m_P());
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fp_2pt_phi_bar(sigma, q2)   + I2_fp_2pt_g_p(sigma, q2);
            const double I3   = I3_fp_2pt_g_p(sigma, q2)       + I3_fp_2pt_g_bar(sigma, q2);
            const double I3d1 = I3d1_fp_2pt_g_p(sigma, q2)     + I3d1_fp_2pt_g_bar(sigma, q2);
            const double I4   = I4_fp_2pt_g_bar(sigma, q2);
            const double I4d1 = I4d1_fp_2pt_g_bar(sigma, q2);
            const double I4d2 = I4d2_fp_2pt_g_bar(sigma, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return result1 + result2;
        }
        double integrand_fp_3pt_m1(const std::array<double, 3> & args, const double & q2) const
        {
            const double sigma  = args[0];
            const double x_1    = args[1];
            const double x_2    = args[2];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double I1 = 0.0;
            const double I2 = I2_fp_3pt_phi_3(sigma, omega_1, omega_2, q2)         + I2_fp_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                            + I2_fp_3pt_phi_4(sigma, omega_1, omega_2, q2)         + I2_fp_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)
                            + I2_fp_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I2_fp_3pt_chi_bar_4(sigma, omega_1, omega_2, q2);
            const double I3 = I3_fp_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                            + I3_fp_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                            + I3_fp_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fp_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                            + I3_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4 = I4_fp_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                            + I4_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result1 = 0.0;
            result1 += - I1;
            result1 +=   I2 / M2;
            result1 += - I3 / (2.0 * M4);
            result1 +=   I4 / (6.0 * M6);
            result1 *=   exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - I2;
            result2 +=   I3 / M2;
            result2 += - I4 / (2.0 * M4);
            result2 *=   exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fp_3pt_A_m1(const std::array<double, 2> & args, const double & sigma, const double & q2) const
        {
            const double x_1    = args[0];
            const double x_2    = args[1];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fp_3pt_phi_3(sigma, omega_1, omega_2, q2)            + I2_fp_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I2_fp_3pt_phi_4(sigma, omega_1, omega_2, q2)            + I2_fp_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)
                              + I2_fp_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I2_fp_3pt_chi_bar_4(sigma, omega_1, omega_2, q2);
            const double I3   = I3_fp_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I3_fp_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fp_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fp_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I3_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I3d1 = I3d1A_fp_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I3d1A_fp_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fp_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fp_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3d1A_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4   = I4_fp_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)    + I4_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I4_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d1 = I4d1A_fp_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d1A_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d1A_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d1A_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d2 = I4d2A_fp_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d2A_fp_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d2A_fp_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d2A_fp_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fp_3pt_B_m1(const double & x_1, const double & sigma, const double & q2) const
        {
            // this ONLY includes the Jacobian from the transformation (omega_1 -> x_1).
            const double prefactor = sigma * m_B();

            const double omega_1 = sigma * m_B() * x_1;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1B_fp_3pt_phi_bar_3(sigma, omega_1, q2)
                                  + I3d1B_fp_3pt_phi_bar_4(sigma, omega_1, q2)     + I3d1B_fp_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fp_3pt_psi_bar_4(sigma, omega_1, q2)     + I3d1B_fp_3pt_chi_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fp_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I3d1B_fp_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1B_fp_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d1B_fp_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d1B_fp_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d1B_fp_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            const     double I4d2 = I4d2B_fp_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d2B_fp_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d2B_fp_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d2B_fp_3pt_chi_bar_bar_4(sigma, omega_1, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fp_3pt_C_m1(const double & x_2, const double & sigma, const double & q2) const
        {
            const double xbar_2 = 1.0 - x_2;

            // this ONLY includes the Jacobian from the transformation (omega_2 -> x_2).
            const double prefactor = sigma * m_B() / (xbar_2 * xbar_2);

            const double omega_2 = sigma * m_B() * (x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1C_fp_3pt_phi_bar_3(sigma, omega_2, q2)
                                  + I3d1C_fp_3pt_phi_bar_4(sigma, omega_2, q2)     + I3d1C_fp_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fp_3pt_psi_bar_4(sigma, omega_2, q2)     + I3d1C_fp_3pt_chi_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fp_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I3d1C_fp_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1C_fp_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d1C_fp_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d1C_fp_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d1C_fp_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            const     double I4d2 = I4d2C_fp_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d2C_fp_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d2C_fp_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d2C_fp_3pt_chi_bar_bar_4(sigma, omega_2, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fp_3pt_D_m1(const double & sigma, const double & q2) const
        {
            // this does NOT includes the original factor of 1 / omega_2
            const double prefactor = 1.0;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);


            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            constexpr double I3d1 = 0.0;
            constexpr double I4   = 0.0;
            constexpr double I4d1 = 0.0;
            const     double I4d2 = I4d2D_fp_3pt_phi_bar_bar_3(sigma, q2) + I4d2D_fp_3pt_phi_bar_bar_4(sigma, q2)
                                  + I4d2D_fp_3pt_psi_bar_bar_4(sigma, q2) + I4d2D_fp_3pt_chi_bar_bar_4(sigma, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }
        // }}}

        /* f_+ : form factor and moments */
        // {{{
        double f_p(const double & q2) const
        {
            const double sigma_0 = this->sigma_0(q2, s0_0_p(), s0_1_p());

            const std::function<double (const double &)> integrand_2pt = std::bind(integrand_fp_2pt, this, std::placeholders::_1, q2);

            const double integral_2pt = integrate<GSL::QAGS>(integrand_2pt, 0.0, sigma_0);
            const double surface_2pt  = 0.0 - surface_fp_2pt(switch_borel ? sigma_0 : 0.0, q2);

            double integral_3pt = 0.0;
            double surface_3pt  = 0.0;

            if (switch_3pt != 0.0)
            {
                const std::function<double (const double &)> surface_3pt_B = std::bind(&Implementation::surface_fp_3pt_B, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const double &)> surface_3pt_C = std::bind(&Implementation::surface_fp_3pt_C, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const std::array<double, 3> &)> integrand_3pt = std::bind(&Implementation::integrand_fp_3pt, this, std::placeholders::_1, q2);
                const std::function<double (const std::array<double, 2> &)> surface_3pt_A = std::bind(&Implementation::surface_fp_3pt_A, this, std::placeholders::_1, sigma_0, q2);

                integral_3pt = integrate(integrand_3pt, { 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0 }, cubature::Config());
                surface_3pt  = 0.0
                             - integrate(surface_3pt_A, { 0.0, 0.0 }, { 1.0, 1.0 }, cubature::Config()) // integrate over x_1 and x_2
                             - integrate<GSL::QAGS>(surface_3pt_B, 0.0, 1.0)                            // integrate over x_1
                             - integrate<GSL::QAGS>(surface_3pt_C, 0.0, 1.0)                            // integrate over x_2
                             - surface_fp_3pt_D(sigma_0, q2);
            }

            return f_B() * m_B() / f_P() * (integral_2pt + surface_2pt + integral_3pt + surface_3pt) / ( Traits::chi2);
        }

        double normalized_moment_1_f_p(const double & q2) const
        {
            const double sigma_0 = this->sigma_0(q2, s0_0_p(), s0_1_p());

            const std::function<double (const double &)> integrand_2pt_m1 = std::bind(&Implementation::integrand_fp_2pt_borel_m1, this, std::placeholders::_1, q2);


            const std::function<double (const double &)> integrand_2pt    = std::bind(&Implementation::integrand_fp_2pt_borel, this, std::placeholders::_1, q2);

            const double integral_2pt_m1 = integrate<GSL::QAGS>(integrand_2pt_m1, 0.0, sigma_0);
            const double surface_2pt_m1  = 0.0 - surface_fp_2pt_m1(sigma_0, q2);

            double integral_3pt_m1 = 0.0;
            double surface_3pt_m1  = 0.0;

            if (switch_3pt != 0.0)
            {
                const std::function<double (const double &)> surface_3pt_B_m1 = std::bind(&Implementation::surface_fp_3pt_B_m1, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const double &)> surface_3pt_C_m1 = std::bind(&Implementation::surface_fp_3pt_C_m1, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const std::array<double, 3> &)> integrand_3pt_m1 = std::bind(&Implementation::integrand_fp_3pt_m1, this, std::placeholders::_1, q2);
                const std::function<double (const std::array<double, 2> &)> surface_3pt_A_m1 = std::bind(&Implementation::surface_fp_3pt_A_m1, this, std::placeholders::_1, sigma_0, q2);

                integral_3pt_m1 = integrate(integrand_3pt_m1, { 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0 }, cubature::Config());
                surface_3pt_m1  = 0.0
                                - integrate(surface_3pt_A_m1, { 0.0, 0.0 }, { 1.0, 1.0 }, cubature::Config()) // integrate over x_1 and x_2
                                - integrate<GSL::QAGS>(surface_3pt_B_m1, 0.0, 1.0)                            // integrate over x_1
                                - integrate<GSL::QAGS>(surface_3pt_C_m1, 0.0, 1.0)                            // integrate over x_2
                                - surface_fp_3pt_D_m1(sigma_0, q2);
            }
            const double numerator       = integral_2pt_m1 + surface_2pt_m1 + integral_3pt_m1 + surface_3pt_m1;

            const double integral_2pt    = integrate<GSL::QAGS>(integrand_2pt, 0.0, sigma_0);
            const double surface_2pt     = 0.0 - surface_fp_2pt(sigma_0, q2);

            double integral_3pt    = 0.0;
            double surface_3pt     = 0.0;

            if (switch_3pt != 0.0)
            {
                const std::function<double (const double &)> surface_3pt_B    = std::bind(&Implementation::surface_fp_3pt_B, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const double &)> surface_3pt_C    = std::bind(&Implementation::surface_fp_3pt_C, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const std::array<double, 3> &)> integrand_3pt = std::bind(&Implementation::integrand_fp_3pt, this, std::placeholders::_1, q2);
                const std::function<double (const std::array<double, 2> &)> surface_3pt_A = std::bind(&Implementation::surface_fp_3pt_A, this, std::placeholders::_1, sigma_0, q2);

                integral_3pt    = integrate(integrand_3pt, { 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0 }, cubature::Config());
                surface_3pt     = 0.0
                                - integrate(surface_3pt_A, { 0.0, 0.0 }, { 1.0, 1.0 }, cubature::Config()) // integrate over x_1 and x_2
                                - integrate<GSL::QAGS>(surface_3pt_B, 0.0, 1.0)                            // integrate over x_1
                                - integrate<GSL::QAGS>(surface_3pt_C, 0.0, 1.0)                            // integrate over x_2
                                - surface_fp_3pt_D(sigma_0, q2);
            }
            const double denominator     = integral_2pt + surface_2pt + integral_3pt + surface_3pt;

            return numerator / denominator;
        }
        // }}}

        /* f_± : 2-particle functions */

        inline
        double I1_fpm_2pt_phi_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to phi_+

            const double sigmabar = 1.0 - sigma;

            const double phi_plus  = this->phi_plus(sigma * m_B);

            const double C_1 =  1.0 / sigmabar - 2.0;

            return C_1 * phi_plus;
        }

        inline
        double I2_fpm_2pt_phi_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to phibar
            const double sigmabar = 1.0 - sigma;
            const double m_v = this->m_v();

            const double phi_bar  = this->phi_bar(sigma * m_B);

            const double C_2 = (2.0 * m_B * sigma * sigmabar - m_v) / (power_of<2>(sigmabar));

            return C_2 * phi_bar;
        }

        inline
        double I2d1_fpm_2pt_phi_bar(const double & sigma, const double & /*q2*/) const
        {
            // first derivative of two-particle contribution to f_± proportional to phibar
            const double sigmabar = 1.0 - sigma;
            const double m_v = this->m_v();

            const double phi_bar     = this->phi_bar(sigma * m_B);
            const double phi_bar_d1  = this->phi_bar_d1(sigma * m_B);

            const double C_2   = (2.0 * m_B * sigma * sigmabar - m_v) / (power_of<2>(sigmabar)) * m_B;
            const double C_2d1 = 2.0 * (m_B * sigmabar - m_v) / (power_of<3>(sigmabar));

            return C_2 * phi_bar_d1 + C_2d1 * phi_bar;
        }

        inline
        double I2_fpm_2pt_g_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to g_+
            const double sigmabar = 1.0 - sigma;

            const double g_plus   = this->g_plus(sigma * m_B);

            const double C_2 = (4.0 - 8.0 * sigmabar) / (power_of<2>(sigmabar));

            return C_2 * g_plus;
        }

        inline
        double I2d1_fpm_2pt_g_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to g_+
            const double sigmabar = 1.0 - sigma;

            const double g_plus    = this->g_plus(sigma * m_B);
            const double g_plus_d1 = this->g_plus_d1(sigma * m_B);

            const double C_2   = (4.0 - 8.0 * sigmabar) / (power_of<2>(sigmabar)) * m_B;
            const double C_2d1 = (8.0 - 8.0 * sigmabar) / (power_of<3>(sigmabar));

            return C_2  * g_plus_d1 + C_2d1 * g_plus;
        }

        inline
        double I3_fpm_2pt_g_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to g_+
            const double sigmabar = 1.0 - sigma;
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_plus    = this->g_plus(sigma * m_B);

            const double C_3 = 8.0 * m_v2 * (2.0 * sigmabar - 1.0) / (power_of<3>(sigmabar));

            return C_3 * g_plus;
        }

        inline
        double I3d1_fpm_2pt_g_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to g_+
            const double sigmabar = 1.0 - sigma;
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_plus    = this->g_plus(sigma * m_B);
            const double g_plus_d1 = this->g_plus_d1(sigma * m_B);

            const double C_3   = 8.0 * m_v2 * (2.0 * sigmabar - 1.0) / (power_of<3>(sigmabar)) * m_B;
            const double C_3d1 = 8.0 * m_v2 * (4.0 * sigmabar - 3.0) / (power_of<4>(sigmabar));

            return C_3 * g_plus_d1 + C_3d1 * g_plus;
        }

        inline
        double I3d2_fpm_2pt_g_p(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to g_+
            const double sigmabar = 1.0 - sigma;
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_plus    = this->g_plus(sigma * m_B);
            const double g_plus_d1 = this->g_plus_d1(sigma * m_B);
            const double g_plus_d2 = this->g_plus_d2(sigma * m_B);

            const double C_3   = 8.0 * m_v2 * (2.0 * sigmabar - 1.0) / (power_of<3>(sigmabar)) * m_B2;
            const double C_3d1 = 16.0 * m_v2 * (4.0 * sigmabar - 3.0) / (power_of<4>(sigmabar)) * m_B;
            const double C_3d2 = - 96.0 * m_v2 * sigma / (power_of<5>(sigmabar));

            return C_3 * g_plus_d2 + C_3d1 * g_plus_d1 + C_3d2 * g_plus;
        }

        inline
        double I3_fpm_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to gbar
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B);

            const double C_3 = 16.0 * sigma * m_B / (power_of<2>(sigmabar));

            return C_3 * g_bar;
        }

        inline
        double I3d1_fpm_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to gbar
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);

            const double C_3   = 16.0 * sigma * m_B / (power_of<2>(sigmabar)) * m_B;
            const double C_3d1 = - 16.0 * (sigmabar - 2.0) * m_B / (power_of<3>(sigmabar));

            return C_3 * g_bar_d1 + C_3d1 * g_bar;
        }

        inline
        double I3d2_fpm_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to gbar
            const double sigmabar = 1.0 - sigma;
            const double m_B2 = power_of<2>(m_B);

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);
            const double g_bar_d2 = this->g_bar_d2(sigma * m_B);

            const double C_3   = 16.0 * sigma * m_B / (power_of<2>(sigmabar)) * m_B2;
            const double C_3d1 = - 32.0 * (sigmabar - 2.0) * m_B2 / (power_of<3>(sigmabar));
            const double C_3d2 = - 32.0 * (sigmabar - 3.0) * m_B / (power_of<4>(sigmabar));

            return C_3 * g_bar_d2 + C_3d1 * g_bar_d1 + C_3d2 * g_bar;
        }

        inline
        double I4_fpm_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to gbar
            const double sigmabar = 1.0 - sigma;
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);

            const double C_4 = 24.0 * m_v2 * (m_v - 2.0 * m_B * sigma * sigmabar) / (power_of<4>(sigmabar));

            return C_4 * g_bar;
        }

        inline
        double I4d1_fpm_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to gbar
            const double sigmabar = 1.0 - sigma;
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);

            const double C_4   = 24.0 * m_v2 * (m_v - 2.0 * m_B * sigma * sigmabar) / (power_of<4>(sigmabar)) * m_B;
            const double C_4d1 = 48.0 * m_v2 * (m_B * sigmabar * (2.0 * sigmabar - 3.0) + 2.0 * m_v) / (power_of<5>(sigmabar));

            return C_4 * g_bar_d1 + C_4d1 * g_bar;
        }

        inline
        double I4d2_fpm_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to gbar
            const double sigmabar = 1.0 - sigma;
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);
            const double g_bar_d2 = this->g_bar_d2(sigma * m_B);

            const double C_4   = 24.0 * m_v2 * (m_v - 2.0 * m_B * sigma * sigmabar) / (power_of<4>(sigmabar)) * m_B2;
            const double C_4d1 = 96.0 * m_v2 * (m_B * sigmabar * (2.0 * sigmabar - 3.0) + 2.0 * m_v) / (power_of<5>(sigmabar)) * m_B;
            const double C_4d2 = 96.0 * m_v2 * (3.0 * m_B * sigmabar * (sigmabar - 2.0) + 5.0 * m_v) / (power_of<6>(sigmabar));

            return C_4 * g_bar_d2 + C_4d1 * g_bar_d1 + C_4d2 * g_bar;
        }

        inline
        double I4d3_fpm_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_± proportional to gbar
            const double sigmabar = 1.0 - sigma;
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);
            const double g_bar_d2 = this->g_bar_d2(sigma * m_B);
            const double g_bar_d3 = this->g_bar_d3(sigma * m_B);

            const double C_4   = 24.0 * m_v2 * (m_v - 2.0 * m_B * sigma * sigmabar) / (power_of<4>(sigmabar)) * m_B2 * m_B;
            const double C_4d1 = 144.0 * m_v2 * (m_B * sigmabar * (2.0 * sigmabar - 3.0) + 2.0 * m_v) / (power_of<5>(sigmabar)) * m_B2;
            const double C_4d2 = 288.0 * m_v2 * (3.0 * m_B * sigmabar * (sigmabar - 2.0) + 5.0 * m_v) / (power_of<6>(sigmabar)) * m_B;
            const double C_4d3 = 576.0 * m_v2 * (m_B * sigmabar * (2.0 * sigmabar - 5.0) + 5.0 * m_v) / (power_of<7>(sigmabar));

            return C_4 * g_bar_d3 + C_4d1 * g_bar_d2 + C_4d2 * g_bar_d1 + C_4d3 * g_bar;
        }

        /* f_± : 3-particle functions */
        // {{{
        double I2_fpm_3pt_phi_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_3
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_3 = this->phi_3(omega_1, omega_2);

            const double C_2 = -(m_B * (2.0 * sigmabar - 3.0) * u + 4.0 * m_v) / (m_B * power_of<2>(sigmabar));

            return C_2 * phi_3;
        }

        double I2_fpm_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_3
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_2 = - 2.0 * sigma * u / (m_B * power_of<3>(sigmabar));

            return C_2 * phi_bar_3;
        }

        double I3_fpm_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_3
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = - 2.0 * (m_B2 * sigmabar2 * (2.0 * sigmabar - 3.0) * u + 4.0 * m_B * m_v * sigmabar * (2.0 * sigmabar - 1.0)
                             + m_v2 * (2.0 * sigmabar * u + u) + q2 * (2.0 * sigmabar - 1.0) * u)/ (m_B * power_of<4>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I3d1A_fpm_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_3
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(),   m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = -(2.0 * (omega_1 * (4.0 * sigma * (-m_v2 + q2 * sigma) + sigmabar2 * (q2 + m_B2 * (1 - 4.0 * sigma)) +
                                (-(10.0 * m_v2) + 3.0 * m_B2 * (1 - 2.0 * sigma2 + sigma) + q2 * (-3.0 + 5.0 * sigma)) * sigmabar) +
                                 m_B * (4.0 * sigma2 * (m_v2 - q2 * sigma) +
                                 sigmabar2 * (3.0 * m_v2 - 8.0 * m_v * omega_2 + q2 - 2.0 * q2 * sigma +
                                 m_B2 * (-1 + 6.0 * sigma2 - 2.0 * sigma)) +
                                 sigmabar * (12.0 * m_v * omega_2 * (1 - 2.0 * sigma) +
                                 sigma * (11.0 * m_v2 + 3.0 * q2 - 6.0 * q2 * sigma - 3.0 * m_B2 * (1 + 2.0 * sigma) * sigmabar)))))
                               / (m_B * omega_2 * power_of<5>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I3d1B_fpm_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_3
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;
            const double sigma2   = power_of<2>(sigma);
            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = -2.0 * (sigmabar * (m_B2 * (-2.0 * sigma2 + sigma + 1.0) + 4.0 * m_B * m_v * (2.0 * sigma - 1.0) -
                                3.0 * m_v2 - q2 * sigmabar) + sigma * (q2 * sigma - m_v2))
                             / ((-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I3d1C_fpm_3pt_phi_bar_3(const double & sigma, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_3
            const double omega_1  = m_B * sigma;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = 8.0 * m_B * m_v * (2.0 * sigma - 1.0) / (omega_2 * power_of<3>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I4_fpm_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_3
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = -6.0 * m_v * u * (m_v * (2.0 * sigmabar + 1.0) * (2.0 * u - 1.0) - 4.0 * m_B * sigma * sigmabar)
                             / (power_of<4>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d1A_fpm_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_3
            const double m_v      = this->m_v(), m_B2 = m_B * m_B;
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * (16.0 * m_B2 * (m_v - omega_2) * sigma3 +
                               m_B * (-(4.0 * omega_1 * omega_2) + 3.0 * m_v * (4.0 * omega_1 + omega_2)) * sigmabar +
                               2.0 * m_v * omega_1 * (2.0 * omega_1 + omega_2) * (-6.0 + sigmabar) -
                               4.0 * sigma * (-(2.0 * m_v * omega_1 * (2.0 * omega_1 + omega_2)) +
                               m_B2 * (3.0 * m_v - 2.0 * omega_2) * sigmabar - 2.0 * m_B * omega_1 * omega_2 * (-2.0 + sigmabar) +
                               m_B * m_v * (4.0 * omega_1 + omega_2) * (-3.0 + sigmabar)) +
                               4.0 * m_B * sigma2 * (4.0 * omega_1 * omega_2 - 2.0 * m_v * (4.0 * omega_1 + omega_2) +
                               3.0 * m_B * m_v * (-2.0 + sigmabar) + m_B * omega_2 * (4.0 - 3.0 * sigmabar)))
                             / (power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d1B_fpm_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_3
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * m_B * (-4.0 * m_B * sigma * sigmabar - 2.0 * m_v * sigma + 3.0 * m_v)
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * phi_bar_bar_3;
        }

        double I4d1C_fpm_3pt_phi_bar_bar_3(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_3

            return 0.0;
        }
        double I4d2A_fpm_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_3
            const double m_v      = this->m_v(), m_B2 = power_of<2>(m_B);
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = 24.0 * m_v * (20.0 * m_B2 * (m_v - omega_2) * sigma3 + m_B2 * (-(3.0 * m_v) + 2.0 * omega_2) * sigmabar2 +
                               m_v * omega_1 * (2.0 * omega_1 + omega_2) * (-15.0 + 4.0 * sigmabar) +
                               m_B * sigmabar * (2.0 * omega_1 * omega_2 * (-4.0 + sigmabar) -
                               m_v * (4.0 * omega_1 + omega_2) * (-6.0 + sigmabar)) +
                               2.0 * m_B * sigma2 * (10.0 * omega_1 * omega_2 - 5.0 * m_v * (4.0 * omega_1 + omega_2) +
                               3.0 * m_B * m_v * (-5.0 + 4.0 * sigmabar) + 2.0 * m_B * omega_2 * (5.0 - 6.0 * sigmabar)) +
                               sigma * (10.0 * m_v * omega_1 * (2.0 * omega_1 + omega_2) +
                               2.0 * m_B2 * sigmabar * (3.0 * m_v * (-4.0 + sigmabar) + omega_2 * (8.0 - 3.0 * sigmabar)) +
                               m_B * (4.0 * omega_1 * omega_2 * (-5.0 + 4.0 * sigmabar) -
                               m_v * (4.0 * omega_1 + omega_2) * (-15.0 + 8.0 * sigmabar))))
                             / (power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d2B_fpm_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_3
            const double omega_2  = m_B * sigma - omega_1;
            const double sigma2   = power_of<2>(sigma);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3     = this->phi_bar_3(omega_1, omega_2);
            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            return  6.0 * m_B * m_v * pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-5.0) *
                   (2.0 * (m_B * (3.0 * m_v + 8.0 * m_B * sigma2 - 4.0 * (m_B + m_v) * sigma) * sigmabar +
                    4.0 * m_B * sigma * (3.0 * m_v - 2.0 * m_v * sigma - 4.0 * m_B * sigma * sigmabar) +
                    2.0 * omega_1 * (-(8.0 * m_B * sigma2) + 2.0 * m_B * sigmabar + m_v * (-6.0 + sigmabar) +
                    4.0 * sigma * (2.0 * m_B + m_v - m_B * sigmabar))) * phi_bar_bar_3 +
                    m_B * (-omega_1 + m_B * sigma) * sigmabar * (3.0 * m_v - 2.0 * m_v * sigma - 4.0 * m_B * sigma * sigmabar) * phi_bar_3);
        }

        double I4d2C_fpm_3pt_phi_bar_bar_3(const double & sigma, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_3
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = -6.0 * m_B2 * m_v * (2.0 * sigma * (m_v - 2.0 * m_B * sigmabar) - 3.0 * m_v) / (power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d2D_fpm_3pt_phi_bar_bar_3(const double & sigma, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_3
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3     = this->phi_bar_3(omega_1, omega_2);
            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            return     (6.0 * m_v * pow(sigmabar,-4.0) * (-(2.0 * (m_v + m_B * (2.0 - 4.0 * sigma)) * phi_bar_bar_3) +
                        m_B * (3.0 * m_v - 2.0 * m_v * sigma - 4.0 * m_B * sigma * sigmabar) * phi_bar_3));
        }

        double I2_fpm_3pt_phi_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_4
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_4 = this->phi_4(omega_1, omega_2);

            const double C_2 = -(2.0 * sigmabar + 1.0) * (u - 1.0) / power_of<2>(sigmabar);

            return C_2 * phi_4;
        }

        double I2_fpm_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_4
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_2 = -2.0 * sigma * (u - 1.0) / (m_B * power_of<3>(sigmabar));

            return C_2 * phi_bar_4;
        }

        double I3_fpm_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(),  m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma,  sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (m_B2 * sigmabar2 * (-2.0 * sigmabar * u + u + 1.0) + m_B * m_v * (1.0 - 4.0 * sigmabar)
                             * sigmabar + m_v2 * (2.0 * sigmabar + 1.0) * (u - 1.0) + q2 * (2.0 * sigmabar - 1.0) * (u - 1.0))
                             / (m_B * power_of<4>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3d1A_fpm_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (sigma3 * (4.0 * m_B * q2 + 6.0 * m_B3 * sigmabar) +
                               sigma2 * (-(4.0 * (omega_1 + omega_2) * q2) - 6.0 * m_B2 * omega_1 * sigmabar +
                               3.0 * m_B3 * sigmabar * (-3.0 + 2.0 * sigmabar) + m_B * (-(4.0 * m_v2) + 6.0 * q2 * sigmabar)) +
                               sigmabar * (9.0 * m_B * m_v * omega_2 + m_B3 * sigmabar -
                               m_B * (3.0 * m_v2 + 4.0 * m_v * omega_2 + q2) * sigmabar +
                              (omega_1 + omega_2) * (10.0 * m_v2 + 3.0 * q2 - q2 * sigmabar) +
                               m_B2 * (3.0 * omega_1 * (-1 + sigmabar) + omega_2 * (-3.0 + sigmabar))) +
                               sigma * (-(12.0 * m_B * m_v * omega_2 * sigmabar) +
                               m_v2 * (4.0 * omega_1 + 4.0 * omega_2 - 11.0 * m_B * sigmabar) +
                               sigmabar * (-(5.0 * (omega_1 + omega_2) * q2) + m_B * q2 * (-3.0 + 2.0 * sigmabar) +
                               m_B3 * (3.0 - 6.0 * sigmabar) + m_B2 * (9.0 * omega_1 + 3.0 * omega_2 - 4.0 * omega_1 * sigmabar))))
                             / (m_B * omega_2 * power_of<5>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3d1B_fpm_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * m_B * (2.0 * sigma * (-m_B * sigma + m_B + 2.0 * m_v) - 3.0 * m_v)
                             / ((-omega_1 + m_B * sigma) * power_of<3>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3d1C_fpm_3pt_phi_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (sigmabar * (m_B2 * (-sigmabar) + m_B * m_v * (3.0 - 4.0 * sigma) + 3.0 * m_v2 - q2 * sigma + q2) + sigma * (m_v2 - q2 * sigma))
                             / (omega_2 * power_of<4>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3_fpm_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * u * (2.0 * m_B * (sigmabar - 2.0) * sigmabar * (2.0 * u - 1.0) + m_v * (4.0 * sigmabar - 3.0))
                             / (m_B * power_of<4>(sigmabar));

            return C_3 * phi_bar_bar_4;
        }

        double I3d1A_fpm_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B),   m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma), sigma4   = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (-(8.0 * m_B * omega_1 * (2.0 * omega_1 + omega_2)) + 16.0 * m_B3 * sigma4 +
                               4.0 * m_v * omega_1 * omega_2 * (-1 + sigmabar) + m_B * m_v * omega_2 * sigmabar +
                               2.0 * m_B2 * (4.0 * omega_1 + omega_2) * sigmabar +
                               8.0 * m_B2 * sigma3 * (-(4.0 * omega_1) - omega_2 + 2.0 * m_B * sigmabar) -
                               2.0 * m_B * sigma2 * (8.0 * m_B2 + 8.0 * m_v * omega_2 - 4.0 * omega_1 * (2.0 * omega_1 + omega_2) +
                               3.0 * m_B * (4.0 * omega_1 + omega_2) * sigmabar) +
                               4.0 * sigma * (4.0 * m_v * omega_1 * omega_2 + 2.0 * m_B2 * (4.0 * omega_1 + omega_2) - 2.0 * m_B3 * sigmabar +
                               m_B * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar + m_B * m_v * (omega_2 - 2.0 * omega_2 * sigmabar)))
                             / (m_B * power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_3 * phi_bar_bar_4;
        }

        double I3d1B_fpm_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (2.0 * m_B * (sigma2 - 1.0) - 4.0 * m_v * sigma + m_v)
                             / ((-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * phi_bar_bar_4;
        }

        double I3d1C_fpm_3pt_phi_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4

            return 0.0;
        }

        double I4_fpm_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2     = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * u * (m_v * (4.0 * sigmabar - 1.0) * (m_B2 * sigmabar2 - q2) - 2.0 * m_B * sigma * sigmabar * (2.0 * u - 1.0) * (m_B2 * sigmabar2 - q2)
                             +m_B * m_v2 * sigmabar * (2.0 * u - 1.0) - m_v3)
                             / (m_B * power_of<5>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d1A_fpm_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma), sigma5 = power_of<5>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -(2.0 * (m_B3 * (-(9.0 * m_v * omega_2 * sigmabar2) + 6.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar4) +
                                 20.0 * m_B3 * sigma5 * q2 + m_v * omega_1 * omega_2 *
                                (m_v2 * (5.0 - 20.0 * sigmabar) + q2 * sigmabar * (-37 + 4.0 * sigmabar)) +
                                 m_B2 * m_v * sigmabar * (m_v * (4.0 * omega_1 + omega_2) * (-2.0 + 5.0 * sigmabar) +
                                 3.0 * omega_1 * omega_2 * (12.0 - 7.0 * sigmabar)) -
                                 2.0 * m_B2 * sigma4 * (5.0 * (4.0 * omega_1 + omega_2) * q2 + 2.0 * m_B * (5.0 * m_v2 - 8.0 * q2 * sigmabar)) +
                                 m_B * (9.0 * m_v * omega_2 * sigmabar2 * q2 -
                                 2.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                 m_v3 * omega_2 * sigmabar * (-1 + 4.0 * sigmabar) -
                                 2.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * (-5.0 + sigmabar * (10.0 + sigmabar))) +
                                 sigma * (-(12.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar4) +
                                 m_v * omega_1 * omega_2 * (-(5.0 * (4.0 * m_v2 + q2)) + 24 * q2 * sigmabar) +
                                 2.0 * m_B3 * sigmabar * (6.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar2 +
                                 m_v2 * (4.0 - 10.0 * sigmabar) + 3.0 * m_v * omega_2 * (-6.0 + 7.0 * sigmabar)) +
                                 2.0 * m_B2 * (2.0 * (4.0 * omega_1 + omega_2) * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                 6.0 * m_v * omega_1 * omega_2 * sigmabar * (-7.0 + 2.0 * sigmabar) +
                                 m_v2 * (4.0 * omega_1 + omega_2) * (-5.0 + 2.0 * sigmabar * (5.0 + sigmabar))) +
                                 m_B * (-(12.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) +
                                 m_v3 * omega_2 * (-5.0 + 24 * sigmabar) + 2.0 * m_v * omega_2 * q2 * sigmabar * (19.0 - 4.0 * sigmabar) +
                                 2.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * (-5.0 + sigmabar * (-7.0 + sigmabar)))) +
                                 2.0 * m_B * sigma3 * (12.0 * m_B4 * sigmabar3 +
                                 5.0 * (-(2.0 * m_v * omega_2) + omega_1 * (2.0 * omega_1 + omega_2)) * q2 +
                                 m_B * (4.0 * omega_1 + omega_2) * (5.0 * m_v2 - 7.0 * q2 * sigmabar) -
                                 2.0 * m_B2 * (4.0 * m_v * (2.0 * m_v + 3.0 * omega_2) * sigmabar +
                                 q2 * (5.0 + sigmabar * (7.0 - 3.0 * sigmabar)))) +
                                 sigma2 * (-(12.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar3) + 36 * m_B5 * sigmabar4 +
                                 20.0 * m_v * omega_1 * omega_2 * q2 + m_B *
                                (20.0 * m_v3 * omega_2 - 10.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) +
                                 m_v * omega_2 * q2 * (5.0 - 28 * sigmabar) + 12.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar) +
                                 2.0 * m_B2 * (m_v * (24 * omega_1 * omega_2 + 7.0 * m_v * (4.0 * omega_1 + omega_2)) * sigmabar -
                                (4.0 * omega_1 + omega_2) * q2 * (-5.0 + sigmabar * (-7.0 + 2.0 * sigmabar))) -
                                 4.0 * m_B3 * (3.0 * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                 3.0 * m_v * omega_2 * sigmabar * (-7.0 + 3.0 * sigmabar) +
                                 m_v2 * (-5.0 + sigmabar * (10.0 + 3.0 * sigmabar))))))
                             / (m_B * power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d1B_fpm_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 2.0 * (6.0 * m_B3 * sigmabar3 * sigma + (m_v + 2.0 * m_B * (-1 + sigma2) - 4.0 * m_v * sigma) *
                             (-m_v2 + q2 * sigma) - sigmabar * (-(4.0 * m_v3) - 2.0 * m_B * q2 * sigma * (-2.0 + sigma) +
                               m_B * m_v2 * (5.0 + 2.0 * sigma) + m_v * q2 * (-9.0 + 4.0 * sigma) -
                               3.0 * m_B2 * m_v * (-3.0 + 4.0 * sigma) * sigmabar))
                             / (power_of<5>(sigmabar) * omega_2);

            return C_4 * phi_bar_bar_4;
        }

        double I4d1C_fpm_3pt_phi_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4

            return 0.0;
        }

        double I4d2A_fpm_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma), sigma5 = power_of<5>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar), sigmabar5 = power_of<5>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 =-(4.0 * (-(6.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar5) + 60 * m_B3 * sigma5 * q2 +
                                5.0 * m_v * omega_1 * omega_2 * (m_v2 * (3.0 - 12.0 * sigmabar) + q2 * sigmabar * (-19.0 + 4.0 * sigmabar)) +
                                m_B3 * sigmabar2 * (12.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar2 +
                                m_v2 * (4.0 - 10.0 * sigmabar) + 3.0 * m_v * omega_2 * (-12.0 + 7.0 * sigmabar)) -
                                30 * m_B2 * sigma4 * ((4.0 * omega_1 + omega_2) * q2 + 2.0 * m_B * (m_v2 - 2.0 * q2 * sigmabar)) +
                                m_B * (m_v * omega_2 * sigmabar2 * q2 * (37 - 4.0 * sigmabar) +
                                5.0 * m_v3 * omega_2 * sigmabar * (-1 + 4.0 * sigmabar) -
                                2.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar * (5.0 + 7.0 * sigmabar) -
                                10.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * (-3.0 + sigmabar * (5.0 + sigmabar))) +
                                2.0 * m_B2 * sigmabar * ((4.0 * omega_1 + omega_2) * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                m_v2 * (4.0 * omega_1 + omega_2) * (-5.0 + sigmabar * (10.0 + sigmabar)) +
                                3.0 * m_v * omega_1 * omega_2 * (15.0 + 2.0 * sigmabar * (-7.0 + sigmabar))) +
                                sigma2 * (-(18.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar3) + 72 * m_B5 * sigmabar4 +
                                60 * m_v * omega_1 * omega_2 * q2 + 5.0 * m_B *
                               (12.0 * m_v3 * omega_2 - 6.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) +
                                8.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar + m_v * omega_2 * q2 * (3.0 - 20.0 * sigmabar)) -
                                2.0 * m_B3 * (3.0 * m_v * omega_2 * sigmabar * (-35 + 24 * sigmabar) +
                                m_v2 * (-30 + 36 * sigmabar2 + 50 * sigmabar) - 6.0 * q2 * sigmabar * (-5.0 + sigmabar * (-7.0 + sigmabar))) +
                                2.0 * m_B2 * (5.0 * m_v * (12.0 * omega_1 * omega_2 + 5.0 * m_v * (4.0 * omega_1 + omega_2)) * sigmabar -
                               (4.0 * omega_1 + omega_2) * q2 * (-15.0 + sigmabar * (-15.0 + 11.0 * sigmabar)))) +
                                sigma * (-24 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar4 + 36 * m_B5 * sigmabar5 -
                                5.0 * m_v * omega_1 * omega_2 * (12.0 * m_v2 + q2 * (3.0 - 16.0 * sigmabar)) +
                                m_B * (-40 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar +
                                4.0 * m_v * omega_2 * q2 * sigmabar * (25 - 11.0 * sigmabar) +
                                5.0 * m_v3 * omega_2 * (-3.0 + 16.0 * sigmabar) +
                                10.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * (-3.0 + sigmabar * (-3.0 + sigmabar))) +
                                2.0 * m_B2 * (3.0 * m_v * omega_1 * omega_2 * sigmabar * (-35 + 16.0 * sigmabar) -
                               (4.0 * omega_1 + omega_2) * q2 * sigmabar * (-10.0 + sigmabar * (-14.0 + sigmabar)) +
                                m_v2 * (4.0 * omega_1 + omega_2) * (-15.0 + sigmabar * (25 + 11.0 * sigmabar))) -
                                2.0 * m_B3 * sigmabar * (3.0 * m_v * omega_2 * (15.0 + 6.0 * sigmabar2 - 28 * sigmabar) +
                                m_v2 * (-20.0 + 6.0 * sigmabar2 + 40 * sigmabar) +
                                3.0 * sigmabar * (-(3.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) + q2 * (2.0 + 4.0 * sigmabar)))) +
                                2.0 * m_B * sigma3 * (18.0 * m_B4 * sigmabar3 +
                                15.0 * (-(2.0 * m_v * omega_2) + omega_1 * (2.0 * omega_1 + omega_2)) * q2 +
                                5.0 * m_B * (4.0 * omega_1 + omega_2) * (3.0 * m_v2 - 5.0 * q2 * sigmabar) -
                                6.0 * m_B2 * (10.0 * m_v * (m_v + omega_2) * sigmabar + q2 * (5.0 + sigmabar * (5.0 - 6.0 * sigmabar))))))
                             / (m_B * power_of<2>(omega_2) * power_of<7>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d2B_fpm_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);

            const double phi_bar_4     = this->phi_bar_4(omega_1, omega_2);
            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            return  pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-6.0) *
                   (4.0 * (-(6.0 * m_B3 * omega_1 * sigmabar4) + 10.0 * m_B2 * sigma4 * q2 +
                    m_v * omega_1 * (m_v2 * (5.0 - 20.0 * sigmabar) + q2 * sigmabar * (-37 + 4.0 * sigmabar)) +
                    m_B2 * m_v * sigmabar * (m_v * (2.0 - 5.0 * sigmabar) + 3.0 * omega_1 * (12.0 - 7.0 * sigmabar)) -
                    2.0 * m_B * sigma3 * (5.0 * (2.0 * m_v + omega_1) * q2 + 24 * m_B2 * m_v * sigmabar +
                    m_B * (5.0 * m_v2 - 7.0 * q2 * sigmabar)) +
                    2.0 * m_B * omega_1 * (q2 * sigmabar * (1 + 2.0 * sigmabar) + m_v2 * (-5.0 + sigmabar * (10.0 + sigmabar))) +
                    sigma * (12.0 * m_B4 * sigmabar4 + m_v * omega_1 * (-(5.0 * (4.0 * m_v2 + q2)) + 24 * q2 * sigmabar) -
                    3.0 * m_B3 * sigmabar * (4.0 * omega_1 * sigmabar2 + m_v * (12.0 - 7.0 * sigmabar)) -
                    2.0 * m_B2 * (6.0 * m_v * omega_1 * sigmabar * (7.0 - 2.0 * sigmabar) +
                    2.0 * q2 * sigmabar * (1 + 2.0 * sigmabar) + m_v2 * (-5.0 + 2.0 * sigmabar * (5.0 + sigmabar))) +
                    m_B * (12.0 * m_v2 * omega_1 * sigmabar + m_v * q2 * sigmabar * (37 - 4.0 * sigmabar) +
                    5.0 * m_v3 * (-1 + 4.0 * sigmabar) - 2.0 * omega_1 * q2 * (-5.0 + sigmabar * (-7.0 + sigmabar)))) +
                    sigma2 * (12.0 * m_B4 * sigmabar3 + 20.0 * m_v * omega_1 * q2 +
                    12.0 * m_B3 * m_v * sigmabar * (7.0 - 2.0 * sigmabar) +
                    m_B * (20.0 * m_v3 + 10.0 * m_v2 * omega_1 + m_v * q2 * (5.0 - 24 * sigmabar) -
                    12.0 * omega_1 * q2 * sigmabar) + 2.0 * m_B2 *
                   (m_v * (-(7.0 * m_v) + 24 * omega_1) * sigmabar + q2 * (-5.0 + sigmabar * (-7.0 + 2.0 * sigmabar))))) * phi_bar_bar_4 +
                    2.0 * m_B * (-omega_1 + m_B * sigma) * sigmabar *
                   (6.0 * m_B3 * sigmabar3 * sigma + 3.0 * m_B2 * m_v * sigmabar2 * (-3.0 + 4.0 * sigma) +
                    m_v3 * (-1 + 4.0 * sigma + 4.0 * sigmabar) +
                    m_v * q2 * (-(4.0 * sigma2) + sigma + 9.0 * sigmabar - 4.0 * sigma * sigmabar) +
                    2.0 * m_B * q2 * sigma * (-1 + sigma2 + (-2.0 + sigma) * sigmabar) -
                    m_B * m_v2 * (-2.0 + 5.0 * sigmabar + 2.0 * sigma * (sigma + sigmabar))) * phi_bar_4);
        }

        double I4d2C_fpm_3pt_phi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 2.0 * m_B * (6.0 * m_B3 * sigmabar3 * sigma +
                             (-m_v2 + q2 * sigma) * (2.0 * m_B * (-1 + sigma2) + m_v * (-1 + 4.0 * sigma)) +
                               sigmabar * (-(4.0 * m_v3) + 2.0 * m_B * q2 * sigma * (-2.0 + sigma) - m_B * m_v2 * (5.0 + 2.0 * sigma) +
                               m_v * q2 * (-9.0 + 4.0 * sigma) - 3.0 * m_B2 * m_v * (-3.0 + 4.0 * sigma) * sigmabar))
                             / (power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d2D_fpm_3pt_phi_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_± proportional to phi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double phi_bar_4     = this->phi_bar_4(omega_1, omega_2);
            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            return    2.0 * pow(m_B,-1) * pow(sigmabar,-6.0) * ((3.0 * m_B2 * m_v * sigmabar2 * (7.0 - 8.0 * sigma) +
                      6.0 * m_B3 * sigmabar3 * (-(2.0 * sigma) + sigmabar) +
                      2.0 * m_B * q2 * (sigma3 - sigma + sigmabar2 * (-2.0 + sigma) + (2.0 * sigma2 - sigmabar) * sigmabar) +
                      m_v * (-1 + 4.0 * sigma + 4.0 * sigmabar) * (m_v2 - q2 * (sigma + sigmabar))) * phi_bar_bar_4 +
                      m_B * sigmabar * (6.0 * m_B3 * sigmabar3 * sigma +
                     (m_v + 2.0 * m_B * (-1 + sigma2) - 4.0 * m_v * sigma) * (-m_v2 + q2 * sigma) -
                      sigmabar * (-(4.0 * m_v3) - 2.0 * m_B * q2 * sigma * (-2.0 + sigma) + m_B * m_v2 * (5.0 + 2.0 * sigma) +
                      m_v * q2 * (-9.0 + 4.0 * sigma) - 3.0 * m_B2 * m_v * (-3.0 + 4.0 * sigma) * sigmabar)) * phi_bar_4);
        }

        double I2_fpm_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_4
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_2 = 2.0 * sigma * (2.0 * u - 1.0) / (m_B * power_of<3>(sigmabar));

            return C_2 * psi_bar_4;
        }

        double I3_fpm_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * ((2.0 * sigmabar - 1.0) * (2.0 * u - 1.0) * (m_B2 * sigmabar2 - q2) + 2.0 * m_B * m_v * sigmabar +
                                m_v2 * (2.0 * sigmabar + 1.0) * (-(2.0 * u - 1.0)))
                             / (m_B * power_of<4>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I3d1A_fpm_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = -(2.0 * (4.0 * sigma3 * (2.0 * m_B * q2 + 3.0 * m_B3 * sigmabar) +
                                 sigmabar * (6.0 * m_B * m_v * omega_2 + 3.0 * m_B2 * (2.0 * omega_1 + omega_2) * (-1 + sigmabar) +
                                 2.0 * m_B3 * sigmabar - 2.0 * m_B * (3.0 * m_v2 + q2) * sigmabar +
                                (2.0 * omega_1 + omega_2) * (10.0 * m_v2 + 3.0 * q2 - q2 * sigmabar)) +
                                 2.0 * sigma2 * (-(2.0 * (2.0 * omega_1 + omega_2) * q2) - 3.0 * m_B2 * (2.0 * omega_1 + omega_2) * sigmabar +
                                 3.0 * m_B3 * sigmabar * (-3.0 + 2.0 * sigmabar) + m_B * (-(4.0 * m_v2) + 6.0 * q2 * sigmabar)) +
                                 sigma * (m_v2 * (8.0 * omega_1 + 4.0 * omega_2 - 22 * m_B * sigmabar) -
                                 sigmabar * (5.0 * (2.0 * omega_1 + omega_2) * q2 + 2.0 * m_B * q2 * (3.0 - 2.0 * sigmabar) +
                                 6.0 * m_B3 * (-1 + 2.0 * sigmabar) + m_B2 * (2.0 * omega_1 + omega_2) * (-9.0 + 4.0 * sigmabar)))))
                             /  (m_B * omega_2 * power_of<5>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I3d1B_fpm_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (sigmabar * (m_B2 * (-(2.0 * sigma - 1.0)) * sigmabar + 2.0 * m_B * m_v - 3.0 * m_v2 - q2 * sigmabar) + sigma * (q2 * sigma - m_v2))
                             / ((-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I3d1C_fpm_3pt_psi_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (sigmabar * (m_B2 * (-(2.0 * sigma - 1.0)) * sigmabar - 2.0 * m_B * m_v - 3.0 * m_v2 - q2 * sigmabar) + sigma * (q2 * sigma - m_v2))
                             / (omega_2 * power_of<4>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I4_fpm_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * m_v * u * (m_v * (2.0 * sigmabar + 1.0) * (2.0 * u - 1.0) - 4.0 * m_B * sigma * sigmabar)
                             / (power_of<4>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1A_fpm_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double m_v      = this->m_v(), m_B2 = m_B * m_B;
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * (16.0 * m_B2 * (m_v - omega_2) * sigma3 +
                               m_B * (-(4.0 * omega_1 * omega_2) + 3.0 * m_v * (4.0 * omega_1 + omega_2)) * sigmabar +
                               2.0 * m_v * omega_1 * (2.0 * omega_1 + omega_2) * (-6.0 + sigmabar) -
                               4.0 * sigma * (-(2.0 * m_v * omega_1 * (2.0 * omega_1 + omega_2)) +
                               m_B2 * (3.0 * m_v - 2.0 * omega_2) * sigmabar - 2.0 * m_B * omega_1 * omega_2 * (-2.0 + sigmabar) +
                               m_B * m_v * (4.0 * omega_1 + omega_2) * (-3.0 + sigmabar)) +
                               4.0 * m_B * sigma2 * (4.0 * omega_1 * omega_2 - 2.0 * m_v * (4.0 * omega_1 + omega_2) +
                               3.0 * m_B * m_v * (-2.0 + sigmabar) + m_B * omega_2 * (4.0 - 3.0 * sigmabar)))
                             / (power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1B_fpm_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * m_B * (-4.0 * m_B * sigma * sigmabar - 2.0 * m_v * sigma + 3.0 * m_v)
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * psi_bar_bar_4;
        }

        double I4d1C_fpm_3pt_psiA_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4

            return 0.0;
        }
        double I4d2A_fpm_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double m_v      = this->m_v(), m_B2 = power_of<2>(m_B);
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 24.0 * m_v * (20.0 * m_B2 * (m_v - omega_2) * sigma3 + m_B2 * (-(3.0 * m_v) + 2.0 * omega_2) * sigmabar2 +
                               m_v * omega_1 * (2.0 * omega_1 + omega_2) * (-15.0 + 4.0 * sigmabar) +
                               m_B * sigmabar * (2.0 * omega_1 * omega_2 * (-4.0 + sigmabar) -
                               m_v * (4.0 * omega_1 + omega_2) * (-6.0 + sigmabar)) +
                               2.0 * m_B * sigma2 * (10.0 * omega_1 * omega_2 - 5.0 * m_v * (4.0 * omega_1 + omega_2) +
                               3.0 * m_B * m_v * (-5.0 + 4.0 * sigmabar) + 2.0 * m_B * omega_2 * (5.0 - 6.0 * sigmabar)) +
                               sigma * (10.0 * m_v * omega_1 * (2.0 * omega_1 + omega_2) +
                               2.0 * m_B2 * sigmabar * (3.0 * m_v * (-4.0 + sigmabar) + omega_2 * (8.0 - 3.0 * sigmabar)) +
                               m_B * (4.0 * omega_1 * omega_2 * (-5.0 + 4.0 * sigmabar) -
                               m_v * (4.0 * omega_1 + omega_2) * (-15.0 + 8.0 * sigmabar))))
                             / (power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2B_fpm_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;
            const double sigma2   = power_of<2>(sigma);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return  6.0 * m_B * m_v * pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-5.0) *
                   (2.0 * (m_B * (3.0 * m_v + 8.0 * m_B * sigma2 - 4.0 * (m_B + m_v) * sigma) * sigmabar +
                    4.0 * m_B * sigma * (3.0 * m_v - 2.0 * m_v * sigma - 4.0 * m_B * sigma * sigmabar) +
                    2.0 * omega_1 * (-(8.0 * m_B * sigma2) + 2.0 * m_B * sigmabar + m_v * (-6.0 + sigmabar) +
                    4.0 * sigma * (2.0 * m_B + m_v - m_B * sigmabar))) * psi_bar_bar_4 +
                    m_B * (-omega_1 + m_B * sigma) * sigmabar * (3.0 * m_v - 2.0 * m_v * sigma - 4.0 * m_B * sigma * sigmabar) * psi_bar_4);
        }

        double I4d2C_fpm_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * m_B2 * m_v * (2.0 * sigma * (m_v - 2.0 * m_B * sigmabar) - 3.0 * m_v) / (power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2D_fpm_3pt_psiA_bar_bar_4(const double & sigma, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return     (6.0 * m_v * pow(sigmabar,-4.0) * (-(2.0 * (m_v + m_B * (2.0 - 4.0 * sigma)) * psi_bar_bar_4) +
                        m_B * (3.0 * m_v - 2.0 * m_v * sigma - 4.0 * m_B * sigma * sigmabar) * psi_bar_4));
        }

        double I3_fpm_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * u * (2.0 * m_B * (sigmabar - 2.0) * sigmabar * (2.0 * u - 1.0) + m_v * (4.0 * sigmabar - 3.0))
                             / (m_B * power_of<4>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1A_fpm_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B),   m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma), sigma4   = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (-(8.0 * m_B * omega_1 * (2.0 * omega_1 + omega_2)) + 16.0 * m_B3 * sigma4 +
                               4.0 * m_v * omega_1 * omega_2 * (-1 + sigmabar) + m_B * m_v * omega_2 * sigmabar +
                               2.0 * m_B2 * (4.0 * omega_1 + omega_2) * sigmabar +
                               8.0 * m_B2 * sigma3 * (-(4.0 * omega_1) - omega_2 + 2.0 * m_B * sigmabar) -
                               2.0 * m_B * sigma2 * (8.0 * m_B2 + 8.0 * m_v * omega_2 - 4.0 * omega_1 * (2.0 * omega_1 + omega_2) +
                               3.0 * m_B * (4.0 * omega_1 + omega_2) * sigmabar) +
                               4.0 * sigma * (4.0 * m_v * omega_1 * omega_2 + 2.0 * m_B2 * (4.0 * omega_1 + omega_2) - 2.0 * m_B3 * sigmabar +
                               m_B * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar + m_B * m_v * (omega_2 - 2.0 * omega_2 * sigmabar)))
                             / (m_B * power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1B_fpm_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (2.0 * m_B * (sigma2 - 1.0) - 4.0 * m_v * sigma + m_v)
                             / ((-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1C_fpm_3pt_psiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4

            return 0.0;
        }

        double I4_fpm_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2     = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * u * (m_v * (4.0 * sigmabar - 1.0) * (m_B2 * sigmabar2 - q2) - 2.0 * m_B * sigma * sigmabar * (2.0 * u - 1.0) * (m_B2 * sigmabar2 - q2)
                             +m_B * m_v2 * sigmabar * (2.0 * u - 1.0) - m_v3)
                             / (m_B * power_of<5>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1A_fpm_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma), sigma5 = power_of<5>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -(2.0 * (m_B3 * (-(9.0 * m_v * omega_2 * sigmabar2) + 6.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar4) +
                                 20.0 * m_B3 * sigma5 * q2 + m_v * omega_1 * omega_2 *
                                (m_v2 * (5.0 - 20.0 * sigmabar) + q2 * sigmabar * (-37 + 4.0 * sigmabar)) +
                                 m_B2 * m_v * sigmabar * (m_v * (4.0 * omega_1 + omega_2) * (-2.0 + 5.0 * sigmabar) +
                                 3.0 * omega_1 * omega_2 * (12.0 - 7.0 * sigmabar)) -
                                 2.0 * m_B2 * sigma4 * (5.0 * (4.0 * omega_1 + omega_2) * q2 + 2.0 * m_B * (5.0 * m_v2 - 8.0 * q2 * sigmabar)) +
                                 m_B * (9.0 * m_v * omega_2 * sigmabar2 * q2 -
                                 2.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                 m_v3 * omega_2 * sigmabar * (-1 + 4.0 * sigmabar) -
                                 2.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * (-5.0 + sigmabar * (10.0 + sigmabar))) +
                                 sigma * (-(12.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar4) +
                                 m_v * omega_1 * omega_2 * (-(5.0 * (4.0 * m_v2 + q2)) + 24 * q2 * sigmabar) +
                                 2.0 * m_B3 * sigmabar * (6.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar2 +
                                 m_v2 * (4.0 - 10.0 * sigmabar) + 3.0 * m_v * omega_2 * (-6.0 + 7.0 * sigmabar)) +
                                 2.0 * m_B2 * (2.0 * (4.0 * omega_1 + omega_2) * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                 6.0 * m_v * omega_1 * omega_2 * sigmabar * (-7.0 + 2.0 * sigmabar) +
                                 m_v2 * (4.0 * omega_1 + omega_2) * (-5.0 + 2.0 * sigmabar * (5.0 + sigmabar))) +
                                 m_B * (-(12.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) +
                                 m_v3 * omega_2 * (-5.0 + 24 * sigmabar) + 2.0 * m_v * omega_2 * q2 * sigmabar * (19.0 - 4.0 * sigmabar) +
                                 2.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * (-5.0 + sigmabar * (-7.0 + sigmabar)))) +
                                 2.0 * m_B * sigma3 * (12.0 * m_B4 * sigmabar3 +
                                 5.0 * (-(2.0 * m_v * omega_2) + omega_1 * (2.0 * omega_1 + omega_2)) * q2 +
                                 m_B * (4.0 * omega_1 + omega_2) * (5.0 * m_v2 - 7.0 * q2 * sigmabar) -
                                 2.0 * m_B2 * (4.0 * m_v * (2.0 * m_v + 3.0 * omega_2) * sigmabar +
                                 q2 * (5.0 + sigmabar * (7.0 - 3.0 * sigmabar)))) +
                                 sigma2 * (-(12.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar3) + 36 * m_B5 * sigmabar4 +
                                 20.0 * m_v * omega_1 * omega_2 * q2 + m_B *
                                (20.0 * m_v3 * omega_2 - 10.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) +
                                 m_v * omega_2 * q2 * (5.0 - 28 * sigmabar) + 12.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar) +
                                 2.0 * m_B2 * (m_v * (24 * omega_1 * omega_2 + 7.0 * m_v * (4.0 * omega_1 + omega_2)) * sigmabar -
                                (4.0 * omega_1 + omega_2) * q2 * (-5.0 + sigmabar * (-7.0 + 2.0 * sigmabar))) -
                                 4.0 * m_B3 * (3.0 * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                 3.0 * m_v * omega_2 * sigmabar * (-7.0 + 3.0 * sigmabar) +
                                 m_v2 * (-5.0 + sigmabar * (10.0 + 3.0 * sigmabar))))))
                             / (m_B * power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1B_fpm_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 2.0 * (6.0 * m_B3 * sigmabar3 * sigma + (m_v + 2.0 * m_B * (-1 + sigma2) - 4.0 * m_v * sigma) *
                             (-m_v2 + q2 * sigma) - sigmabar * (-(4.0 * m_v3) - 2.0 * m_B * q2 * sigma * (-2.0 + sigma) +
                               m_B * m_v2 * (5.0 + 2.0 * sigma) + m_v * q2 * (-9.0 + 4.0 * sigma) -
                               3.0 * m_B2 * m_v * (-3.0 + 4.0 * sigma) * sigmabar))
                             / (power_of<5>(sigmabar) * omega_2);

            return C_4 * psi_bar_bar_4;
        }

        double I4d1C_fpm_3pt_psiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4

            return 0.0;
        }

        double I4d2A_fpm_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma), sigma5 = power_of<5>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar), sigmabar5 = power_of<5>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 =-(4.0 * (-(6.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar5) + 60 * m_B3 * sigma5 * q2 +
                                5.0 * m_v * omega_1 * omega_2 * (m_v2 * (3.0 - 12.0 * sigmabar) + q2 * sigmabar * (-19.0 + 4.0 * sigmabar)) +
                                m_B3 * sigmabar2 * (12.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar2 +
                                m_v2 * (4.0 - 10.0 * sigmabar) + 3.0 * m_v * omega_2 * (-12.0 + 7.0 * sigmabar)) -
                                30 * m_B2 * sigma4 * ((4.0 * omega_1 + omega_2) * q2 + 2.0 * m_B * (m_v2 - 2.0 * q2 * sigmabar)) +
                                m_B * (m_v * omega_2 * sigmabar2 * q2 * (37 - 4.0 * sigmabar) +
                                5.0 * m_v3 * omega_2 * sigmabar * (-1 + 4.0 * sigmabar) -
                                2.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar * (5.0 + 7.0 * sigmabar) -
                                10.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * (-3.0 + sigmabar * (5.0 + sigmabar))) +
                                2.0 * m_B2 * sigmabar * ((4.0 * omega_1 + omega_2) * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                m_v2 * (4.0 * omega_1 + omega_2) * (-5.0 + sigmabar * (10.0 + sigmabar)) +
                                3.0 * m_v * omega_1 * omega_2 * (15.0 + 2.0 * sigmabar * (-7.0 + sigmabar))) +
                                sigma2 * (-(18.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar3) + 72 * m_B5 * sigmabar4 +
                                60 * m_v * omega_1 * omega_2 * q2 + 5.0 * m_B *
                               (12.0 * m_v3 * omega_2 - 6.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) +
                                8.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar + m_v * omega_2 * q2 * (3.0 - 20.0 * sigmabar)) -
                                2.0 * m_B3 * (3.0 * m_v * omega_2 * sigmabar * (-35 + 24 * sigmabar) +
                                m_v2 * (-30 + 36 * sigmabar2 + 50 * sigmabar) - 6.0 * q2 * sigmabar * (-5.0 + sigmabar * (-7.0 + sigmabar))) +
                                2.0 * m_B2 * (5.0 * m_v * (12.0 * omega_1 * omega_2 + 5.0 * m_v * (4.0 * omega_1 + omega_2)) * sigmabar -
                               (4.0 * omega_1 + omega_2) * q2 * (-15.0 + sigmabar * (-15.0 + 11.0 * sigmabar)))) +
                                sigma * (-24 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar4 + 36 * m_B5 * sigmabar5 -
                                5.0 * m_v * omega_1 * omega_2 * (12.0 * m_v2 + q2 * (3.0 - 16.0 * sigmabar)) +
                                m_B * (-40 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar +
                                4.0 * m_v * omega_2 * q2 * sigmabar * (25 - 11.0 * sigmabar) +
                                5.0 * m_v3 * omega_2 * (-3.0 + 16.0 * sigmabar) +
                                10.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * (-3.0 + sigmabar * (-3.0 + sigmabar))) +
                                2.0 * m_B2 * (3.0 * m_v * omega_1 * omega_2 * sigmabar * (-35 + 16.0 * sigmabar) -
                               (4.0 * omega_1 + omega_2) * q2 * sigmabar * (-10.0 + sigmabar * (-14.0 + sigmabar)) +
                                m_v2 * (4.0 * omega_1 + omega_2) * (-15.0 + sigmabar * (25 + 11.0 * sigmabar))) -
                                2.0 * m_B3 * sigmabar * (3.0 * m_v * omega_2 * (15.0 + 6.0 * sigmabar2 - 28 * sigmabar) +
                                m_v2 * (-20.0 + 6.0 * sigmabar2 + 40 * sigmabar) +
                                3.0 * sigmabar * (-(3.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) + q2 * (2.0 + 4.0 * sigmabar)))) +
                                2.0 * m_B * sigma3 * (18.0 * m_B4 * sigmabar3 +
                                15.0 * (-(2.0 * m_v * omega_2) + omega_1 * (2.0 * omega_1 + omega_2)) * q2 +
                                5.0 * m_B * (4.0 * omega_1 + omega_2) * (3.0 * m_v2 - 5.0 * q2 * sigmabar) -
                                6.0 * m_B2 * (10.0 * m_v * (m_v + omega_2) * sigmabar + q2 * (5.0 + sigmabar * (5.0 - 6.0 * sigmabar))))))
                             / (m_B * power_of<2>(omega_2) * power_of<7>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2B_fpm_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return  pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-6.0) *
                   (4.0 * (-(6.0 * m_B3 * omega_1 * sigmabar4) + 10.0 * m_B2 * sigma4 * q2 +
                    m_v * omega_1 * (m_v2 * (5.0 - 20.0 * sigmabar) + q2 * sigmabar * (-37 + 4.0 * sigmabar)) +
                    m_B2 * m_v * sigmabar * (m_v * (2.0 - 5.0 * sigmabar) + 3.0 * omega_1 * (12.0 - 7.0 * sigmabar)) -
                    2.0 * m_B * sigma3 * (5.0 * (2.0 * m_v + omega_1) * q2 + 24 * m_B2 * m_v * sigmabar +
                    m_B * (5.0 * m_v2 - 7.0 * q2 * sigmabar)) +
                    2.0 * m_B * omega_1 * (q2 * sigmabar * (1 + 2.0 * sigmabar) + m_v2 * (-5.0 + sigmabar * (10.0 + sigmabar))) +
                    sigma * (12.0 * m_B4 * sigmabar4 + m_v * omega_1 * (-(5.0 * (4.0 * m_v2 + q2)) + 24 * q2 * sigmabar) -
                    3.0 * m_B3 * sigmabar * (4.0 * omega_1 * sigmabar2 + m_v * (12.0 - 7.0 * sigmabar)) -
                    2.0 * m_B2 * (6.0 * m_v * omega_1 * sigmabar * (7.0 - 2.0 * sigmabar) +
                    2.0 * q2 * sigmabar * (1 + 2.0 * sigmabar) + m_v2 * (-5.0 + 2.0 * sigmabar * (5.0 + sigmabar))) +
                    m_B * (12.0 * m_v2 * omega_1 * sigmabar + m_v * q2 * sigmabar * (37 - 4.0 * sigmabar) +
                    5.0 * m_v3 * (-1 + 4.0 * sigmabar) - 2.0 * omega_1 * q2 * (-5.0 + sigmabar * (-7.0 + sigmabar)))) +
                    sigma2 * (12.0 * m_B4 * sigmabar3 + 20.0 * m_v * omega_1 * q2 +
                    12.0 * m_B3 * m_v * sigmabar * (7.0 - 2.0 * sigmabar) +
                    m_B * (20.0 * m_v3 + 10.0 * m_v2 * omega_1 + m_v * q2 * (5.0 - 24 * sigmabar) -
                    12.0 * omega_1 * q2 * sigmabar) + 2.0 * m_B2 *
                   (m_v * (-(7.0 * m_v) + 24 * omega_1) * sigmabar + q2 * (-5.0 + sigmabar * (-7.0 + 2.0 * sigmabar))))) * psi_bar_bar_4 +
                    2.0 * m_B * (-omega_1 + m_B * sigma) * sigmabar *
                   (6.0 * m_B3 * sigmabar3 * sigma + 3.0 * m_B2 * m_v * sigmabar2 * (-3.0 + 4.0 * sigma) +
                    m_v3 * (-1 + 4.0 * sigma + 4.0 * sigmabar) +
                    m_v * q2 * (-(4.0 * sigma2) + sigma + 9.0 * sigmabar - 4.0 * sigma * sigmabar) +
                    2.0 * m_B * q2 * sigma * (-1 + sigma2 + (-2.0 + sigma) * sigmabar) -
                    m_B * m_v2 * (-2.0 + 5.0 * sigmabar + 2.0 * sigma * (sigma + sigmabar))) * psi_bar_4);
        }

        double I4d2C_fpm_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 2.0 * m_B * (6.0 * m_B3 * sigmabar3 * sigma +
                             (-m_v2 + q2 * sigma) * (2.0 * m_B * (-1 + sigma2) + m_v * (-1 + 4.0 * sigma)) +
                               sigmabar * (-(4.0 * m_v3) + 2.0 * m_B * q2 * sigma * (-2.0 + sigma) - m_B * m_v2 * (5.0 + 2.0 * sigma) +
                               m_v * q2 * (-9.0 + 4.0 * sigma) - 3.0 * m_B2 * m_v * (-3.0 + 4.0 * sigma) * sigmabar))
                             / (power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2D_fpm_3pt_psiB_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_± proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return    2.0 * pow(m_B,-1) * pow(sigmabar,-6.0) * ((3.0 * m_B2 * m_v * sigmabar2 * (7.0 - 8.0 * sigma) +
                      6.0 * m_B3 * sigmabar3 * (-(2.0 * sigma) + sigmabar) +
                      2.0 * m_B * q2 * (sigma3 - sigma + sigmabar2 * (-2.0 + sigma) + (2.0 * sigma2 - sigmabar) * sigmabar) +
                      m_v * (-1 + 4.0 * sigma + 4.0 * sigmabar) * (m_v2 - q2 * (sigma + sigmabar))) * psi_bar_bar_4 +
                      m_B * sigmabar * (6.0 * m_B3 * sigmabar3 * sigma +
                     (m_v + 2.0 * m_B * (-1 + sigma2) - 4.0 * m_v * sigma) * (-m_v2 + q2 * sigma) -
                      sigmabar * (-(4.0 * m_v3) - 2.0 * m_B * q2 * sigma * (-2.0 + sigma) + m_B * m_v2 * (5.0 + 2.0 * sigma) +
                      m_v * q2 * (-9.0 + 4.0 * sigma) - 3.0 * m_B2 * m_v * (-3.0 + 4.0 * sigma) * sigmabar)) * psi_bar_4);
        }

        double I3_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - 0.0                                                       - I3_fpm_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1A_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - 0.0                                                       - I3d1A_fpm_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1B_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return - 0.0                                                       - I3d1B_fpm_3pt_psiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I3d1C_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return - 0.0                                                       - I3d1C_fpm_3pt_psiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I4_fpm_3pt_psiA_bar_bar_4( sigma, omega_1, omega_2, q2)    - I4_fpm_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1A_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I4d1A_fpm_3pt_psiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d1A_fpm_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1B_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return - I4d1B_fpm_3pt_psiA_bar_bar_4( sigma, omega_1, q2)          - I4d1B_fpm_3pt_psiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d1C_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return - I4d1C_fpm_3pt_psiA_bar_bar_4( sigma, omega_2, q2)          - I4d1C_fpm_3pt_psiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2A_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I4d2A_fpm_3pt_psiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d2A_fpm_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d2B_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return - I4d2B_fpm_3pt_psiA_bar_bar_4( sigma, omega_1, q2)          - I4d2B_fpm_3pt_psiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d2C_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return - I4d2C_fpm_3pt_psiA_bar_bar_4( sigma, omega_2, q2)          - I4d2C_fpm_3pt_psiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2D_fpm_3pt_psi_bar_bar_4(const double & sigma, const double & q2) const
        {
            return - I4d2D_fpm_3pt_psiA_bar_bar_4( sigma, q2)                   - I4d2D_fpm_3pt_psiB_bar_bar_4(sigma, q2);
        }

        double I2_fpm_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_4
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_2 = -2.0 * sigma / (m_B * power_of<3>(sigmabar));

            return C_2 * chi_bar_4;
        }

        double I3_fpm_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (m_B2 * sigmabar2 * (-2.0 * sigmabar - 4.0 * sigma * u + 1.0) + 4.0 * m_B * m_v * sigmabar * (2.0 * sigmabar - 1.0)
                             + m_v2 * (2.0 * sigmabar + 1.0) + q2 * (2.0 * sigmabar - 1.0))
                             / (m_B * power_of<4>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I3d1A_fpm_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = -(2.0 * (12.0 * m_B3 * sigma3 * sigmabar - 2.0 * sigma2 *
                                (2.0 * omega_2 * q2 + 3.0 * m_B2 * sigmabar * (2.0 * m_B + 2.0 * omega_1 + omega_2 - 2.0 * m_B * sigmabar)) +
                                 sigmabar * (4.0 * m_B2 * omega_1 * sigmabar +
                                 omega_2 * (10.0 * m_v2 + 3.0 * q2 + 3.0 * m_B2 * (-1 + sigmabar) - q2 * sigmabar +
                                 4.0 * m_B * m_v * (3.0 - 2.0 * sigmabar))) +
                                 sigma * (4.0 * m_v2 * omega_2 - 24 * m_B * m_v * omega_2 * sigmabar -
                                 sigmabar * (5.0 * omega_2 * q2 + 8.0 * m_B3 * sigmabar +
                                 m_B2 * (-(12.0 * omega_1) - 9.0 * omega_2 + 8.0 * omega_1 * sigmabar + 4.0 * omega_2 * sigmabar)))))
                             / (m_B * omega_2 * power_of<5>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I3d1B_fpm_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (sigmabar * (m_B2 * (-2.0 * sigma2 + sigma + 1.0) + 4.0 * m_B * m_v
                             * (2.0 * sigma - 1.0) - 3.0 * m_v2 - q2 * sigmabar) + sigma * (q2 * sigma - m_v2))
                             / ((-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I3d1C_fpm_3pt_chi_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (sigmabar * (m_B2 * (-(2.0 * sigma - 1.0)) * sigmabar + 4.0 * m_B * m_v
                             * (2.0 * sigma - 1.0) - 3.0 * m_v2 - q2 * sigmabar) + sigma * (q2 * sigma - m_v2))
                             / (omega_2 * power_of<4>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I4_fpm_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * m_v * u * (m_v * (2.0 * sigmabar + 1.0) * (2.0 * u - 1.0) - 4.0 * m_B * sigma * sigmabar)
                             / (power_of<4>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1A_fpm_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double m_v      = this->m_v(), m_B2 = m_B * m_B;
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * (16.0 * m_B2 * (m_v - omega_2) * sigma3 +
                               m_B * (-(4.0 * omega_1 * omega_2) + 3.0 * m_v * (4.0 * omega_1 + omega_2)) * sigmabar +
                               2.0 * m_v * omega_1 * (2.0 * omega_1 + omega_2) * (-6.0 + sigmabar) -
                               4.0 * sigma * (-(2.0 * m_v * omega_1 * (2.0 * omega_1 + omega_2)) +
                               m_B2 * (3.0 * m_v - 2.0 * omega_2) * sigmabar - 2.0 * m_B * omega_1 * omega_2 * (-2.0 + sigmabar) +
                               m_B * m_v * (4.0 * omega_1 + omega_2) * (-3.0 + sigmabar)) +
                               4.0 * m_B * sigma2 * (4.0 * omega_1 * omega_2 - 2.0 * m_v * (4.0 * omega_1 + omega_2) +
                               3.0 * m_B * m_v * (-2.0 + sigmabar) + m_B * omega_2 * (4.0 - 3.0 * sigmabar)))
                             / (power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1B_fpm_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * m_v * m_B * (-4.0 * m_B * sigma * sigmabar - 2.0 * m_v * sigma + 3.0 * m_v)
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * chi_bar_bar_4;
        }

        double I4d1C_fpm_3pt_chiA_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4

            return 0.0;
        }
        double I4d2A_fpm_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double m_v      = this->m_v(), m_B2 = power_of<2>(m_B);
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 24.0 * m_v * (20.0 * m_B2 * (m_v - omega_2) * sigma3 + m_B2 * (-(3.0 * m_v) + 2.0 * omega_2) * sigmabar2 +
                               m_v * omega_1 * (2.0 * omega_1 + omega_2) * (-15.0 + 4.0 * sigmabar) +
                               m_B * sigmabar * (2.0 * omega_1 * omega_2 * (-4.0 + sigmabar) -
                               m_v * (4.0 * omega_1 + omega_2) * (-6.0 + sigmabar)) +
                               2.0 * m_B * sigma2 * (10.0 * omega_1 * omega_2 - 5.0 * m_v * (4.0 * omega_1 + omega_2) +
                               3.0 * m_B * m_v * (-5.0 + 4.0 * sigmabar) + 2.0 * m_B * omega_2 * (5.0 - 6.0 * sigmabar)) +
                               sigma * (10.0 * m_v * omega_1 * (2.0 * omega_1 + omega_2) +
                               2.0 * m_B2 * sigmabar * (3.0 * m_v * (-4.0 + sigmabar) + omega_2 * (8.0 - 3.0 * sigmabar)) +
                               m_B * (4.0 * omega_1 * omega_2 * (-5.0 + 4.0 * sigmabar) -
                               m_v * (4.0 * omega_1 + omega_2) * (-15.0 + 8.0 * sigmabar))))
                             / (power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2B_fpm_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;
            const double sigma2   = power_of<2>(sigma);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return  6.0 * m_B * m_v * pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-5.0) *
                   (2.0 * (m_B * (3.0 * m_v + 8.0 * m_B * sigma2 - 4.0 * (m_B + m_v) * sigma) * sigmabar +
                    4.0 * m_B * sigma * (3.0 * m_v - 2.0 * m_v * sigma - 4.0 * m_B * sigma * sigmabar) +
                    2.0 * omega_1 * (-(8.0 * m_B * sigma2) + 2.0 * m_B * sigmabar + m_v * (-6.0 + sigmabar) +
                    4.0 * sigma * (2.0 * m_B + m_v - m_B * sigmabar))) * chi_bar_bar_4 +
                    m_B * (-omega_1 + m_B * sigma) * sigmabar * (3.0 * m_v - 2.0 * m_v * sigma - 4.0 * m_B * sigma * sigmabar) * chi_bar_4);
        }

        double I4d2C_fpm_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -6.0 * m_B2 * m_v * (2.0 * sigma * (m_v - 2.0 * m_B * sigmabar) - 3.0 * m_v) / (power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2D_fpm_3pt_chiA_bar_bar_4(const double & sigma, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return     (6.0 * m_v * pow(sigmabar,-4.0) * (-(2.0 * (m_v + m_B * (2.0 - 4.0 * sigma)) * chi_bar_bar_4) +
                        m_B * (3.0 * m_v - 2.0 * m_v * sigma - 4.0 * m_B * sigma * sigmabar) * chi_bar_4));
        }

        double I3_fpm_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * u * (2.0 * m_B * (sigmabar - 2.0) * sigmabar * (2.0 * u - 1.0) + m_v * (4.0 * sigmabar - 3.0))
                             / (m_B * power_of<4>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1A_fpm_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B),   m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma), sigma4   = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (-(8.0 * m_B * omega_1 * (2.0 * omega_1 + omega_2)) + 16.0 * m_B3 * sigma4 +
                               4.0 * m_v * omega_1 * omega_2 * (-1 + sigmabar) + m_B * m_v * omega_2 * sigmabar +
                               2.0 * m_B2 * (4.0 * omega_1 + omega_2) * sigmabar +
                               8.0 * m_B2 * sigma3 * (-(4.0 * omega_1) - omega_2 + 2.0 * m_B * sigmabar) -
                               2.0 * m_B * sigma2 * (8.0 * m_B2 + 8.0 * m_v * omega_2 - 4.0 * omega_1 * (2.0 * omega_1 + omega_2) +
                               3.0 * m_B * (4.0 * omega_1 + omega_2) * sigmabar) +
                               4.0 * sigma * (4.0 * m_v * omega_1 * omega_2 + 2.0 * m_B2 * (4.0 * omega_1 + omega_2) - 2.0 * m_B3 * sigmabar +
                               m_B * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar + m_B * m_v * (omega_2 - 2.0 * omega_2 * sigmabar)))
                             / (m_B * power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1B_fpm_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (2.0 * m_B * (sigma2 - 1.0) - 4.0 * m_v * sigma + m_v)
                             / ((-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1C_fpm_3pt_chiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4

            return 0.0;
        }

        double I4_fpm_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2     = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * u * (m_v * (4.0 * sigmabar - 1.0) * (m_B2 * sigmabar2 - q2) - 2.0 * m_B * sigma * sigmabar * (2.0 * u - 1.0) * (m_B2 * sigmabar2 - q2)
                             +m_B * m_v2 * sigmabar * (2.0 * u - 1.0) - m_v3)
                             / (m_B * power_of<5>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1A_fpm_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma), sigma5 = power_of<5>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -(2.0 * (m_B3 * (-(9.0 * m_v * omega_2 * sigmabar2) + 6.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar4) +
                                 20.0 * m_B3 * sigma5 * q2 + m_v * omega_1 * omega_2 *
                                (m_v2 * (5.0 - 20.0 * sigmabar) + q2 * sigmabar * (-37 + 4.0 * sigmabar)) +
                                 m_B2 * m_v * sigmabar * (m_v * (4.0 * omega_1 + omega_2) * (-2.0 + 5.0 * sigmabar) +
                                 3.0 * omega_1 * omega_2 * (12.0 - 7.0 * sigmabar)) -
                                 2.0 * m_B2 * sigma4 * (5.0 * (4.0 * omega_1 + omega_2) * q2 + 2.0 * m_B * (5.0 * m_v2 - 8.0 * q2 * sigmabar)) +
                                 m_B * (9.0 * m_v * omega_2 * sigmabar2 * q2 -
                                 2.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                 m_v3 * omega_2 * sigmabar * (-1 + 4.0 * sigmabar) -
                                 2.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * (-5.0 + sigmabar * (10.0 + sigmabar))) +
                                 sigma * (-(12.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar4) +
                                 m_v * omega_1 * omega_2 * (-(5.0 * (4.0 * m_v2 + q2)) + 24 * q2 * sigmabar) +
                                 2.0 * m_B3 * sigmabar * (6.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar2 +
                                 m_v2 * (4.0 - 10.0 * sigmabar) + 3.0 * m_v * omega_2 * (-6.0 + 7.0 * sigmabar)) +
                                 2.0 * m_B2 * (2.0 * (4.0 * omega_1 + omega_2) * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                 6.0 * m_v * omega_1 * omega_2 * sigmabar * (-7.0 + 2.0 * sigmabar) +
                                 m_v2 * (4.0 * omega_1 + omega_2) * (-5.0 + 2.0 * sigmabar * (5.0 + sigmabar))) +
                                 m_B * (-(12.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) +
                                 m_v3 * omega_2 * (-5.0 + 24 * sigmabar) + 2.0 * m_v * omega_2 * q2 * sigmabar * (19.0 - 4.0 * sigmabar) +
                                 2.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * (-5.0 + sigmabar * (-7.0 + sigmabar)))) +
                                 2.0 * m_B * sigma3 * (12.0 * m_B4 * sigmabar3 +
                                 5.0 * (-(2.0 * m_v * omega_2) + omega_1 * (2.0 * omega_1 + omega_2)) * q2 +
                                 m_B * (4.0 * omega_1 + omega_2) * (5.0 * m_v2 - 7.0 * q2 * sigmabar) -
                                 2.0 * m_B2 * (4.0 * m_v * (2.0 * m_v + 3.0 * omega_2) * sigmabar +
                                 q2 * (5.0 + sigmabar * (7.0 - 3.0 * sigmabar)))) +
                                 sigma2 * (-(12.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar3) + 36 * m_B5 * sigmabar4 +
                                 20.0 * m_v * omega_1 * omega_2 * q2 + m_B *
                                (20.0 * m_v3 * omega_2 - 10.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) +
                                 m_v * omega_2 * q2 * (5.0 - 28 * sigmabar) + 12.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar) +
                                 2.0 * m_B2 * (m_v * (24 * omega_1 * omega_2 + 7.0 * m_v * (4.0 * omega_1 + omega_2)) * sigmabar -
                                (4.0 * omega_1 + omega_2) * q2 * (-5.0 + sigmabar * (-7.0 + 2.0 * sigmabar))) -
                                 4.0 * m_B3 * (3.0 * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                 3.0 * m_v * omega_2 * sigmabar * (-7.0 + 3.0 * sigmabar) +
                                 m_v2 * (-5.0 + sigmabar * (10.0 + 3.0 * sigmabar))))))
                             / (m_B * power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1B_fpm_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 2.0 * (6.0 * m_B3 * sigmabar3 * sigma + (m_v + 2.0 * m_B * (-1 + sigma2) - 4.0 * m_v * sigma) *
                             (-m_v2 + q2 * sigma) - sigmabar * (-(4.0 * m_v3) - 2.0 * m_B * q2 * sigma * (-2.0 + sigma) +
                               m_B * m_v2 * (5.0 + 2.0 * sigma) + m_v * q2 * (-9.0 + 4.0 * sigma) -
                               3.0 * m_B2 * m_v * (-3.0 + 4.0 * sigma) * sigmabar))
                             / (power_of<5>(sigmabar) * omega_2);

            return C_4 * chi_bar_bar_4;
        }

        double I4d1C_fpm_3pt_chiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4

            return 0.0;
        }

        double I4d2A_fpm_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma), sigma5 = power_of<5>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar), sigmabar5 = power_of<5>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 =-(4.0 * (-(6.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar5) + 60 * m_B3 * sigma5 * q2 +
                                5.0 * m_v * omega_1 * omega_2 * (m_v2 * (3.0 - 12.0 * sigmabar) + q2 * sigmabar * (-19.0 + 4.0 * sigmabar)) +
                                m_B3 * sigmabar2 * (12.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar2 +
                                m_v2 * (4.0 - 10.0 * sigmabar) + 3.0 * m_v * omega_2 * (-12.0 + 7.0 * sigmabar)) -
                                30 * m_B2 * sigma4 * ((4.0 * omega_1 + omega_2) * q2 + 2.0 * m_B * (m_v2 - 2.0 * q2 * sigmabar)) +
                                m_B * (m_v * omega_2 * sigmabar2 * q2 * (37 - 4.0 * sigmabar) +
                                5.0 * m_v3 * omega_2 * sigmabar * (-1 + 4.0 * sigmabar) -
                                2.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar * (5.0 + 7.0 * sigmabar) -
                                10.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * (-3.0 + sigmabar * (5.0 + sigmabar))) +
                                2.0 * m_B2 * sigmabar * ((4.0 * omega_1 + omega_2) * q2 * sigmabar * (1 + 2.0 * sigmabar) +
                                m_v2 * (4.0 * omega_1 + omega_2) * (-5.0 + sigmabar * (10.0 + sigmabar)) +
                                3.0 * m_v * omega_1 * omega_2 * (15.0 + 2.0 * sigmabar * (-7.0 + sigmabar))) +
                                sigma2 * (-(18.0 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar3) + 72 * m_B5 * sigmabar4 +
                                60 * m_v * omega_1 * omega_2 * q2 + 5.0 * m_B *
                               (12.0 * m_v3 * omega_2 - 6.0 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) +
                                8.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * sigmabar + m_v * omega_2 * q2 * (3.0 - 20.0 * sigmabar)) -
                                2.0 * m_B3 * (3.0 * m_v * omega_2 * sigmabar * (-35 + 24 * sigmabar) +
                                m_v2 * (-30 + 36 * sigmabar2 + 50 * sigmabar) - 6.0 * q2 * sigmabar * (-5.0 + sigmabar * (-7.0 + sigmabar))) +
                                2.0 * m_B2 * (5.0 * m_v * (12.0 * omega_1 * omega_2 + 5.0 * m_v * (4.0 * omega_1 + omega_2)) * sigmabar -
                               (4.0 * omega_1 + omega_2) * q2 * (-15.0 + sigmabar * (-15.0 + 11.0 * sigmabar)))) +
                                sigma * (-24 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar4 + 36 * m_B5 * sigmabar5 -
                                5.0 * m_v * omega_1 * omega_2 * (12.0 * m_v2 + q2 * (3.0 - 16.0 * sigmabar)) +
                                m_B * (-40 * m_v2 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar +
                                4.0 * m_v * omega_2 * q2 * sigmabar * (25 - 11.0 * sigmabar) +
                                5.0 * m_v3 * omega_2 * (-3.0 + 16.0 * sigmabar) +
                                10.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * (-3.0 + sigmabar * (-3.0 + sigmabar))) +
                                2.0 * m_B2 * (3.0 * m_v * omega_1 * omega_2 * sigmabar * (-35 + 16.0 * sigmabar) -
                               (4.0 * omega_1 + omega_2) * q2 * sigmabar * (-10.0 + sigmabar * (-14.0 + sigmabar)) +
                                m_v2 * (4.0 * omega_1 + omega_2) * (-15.0 + sigmabar * (25 + 11.0 * sigmabar))) -
                                2.0 * m_B3 * sigmabar * (3.0 * m_v * omega_2 * (15.0 + 6.0 * sigmabar2 - 28 * sigmabar) +
                                m_v2 * (-20.0 + 6.0 * sigmabar2 + 40 * sigmabar) +
                                3.0 * sigmabar * (-(3.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) + q2 * (2.0 + 4.0 * sigmabar)))) +
                                2.0 * m_B * sigma3 * (18.0 * m_B4 * sigmabar3 +
                                15.0 * (-(2.0 * m_v * omega_2) + omega_1 * (2.0 * omega_1 + omega_2)) * q2 +
                                5.0 * m_B * (4.0 * omega_1 + omega_2) * (3.0 * m_v2 - 5.0 * q2 * sigmabar) -
                                6.0 * m_B2 * (10.0 * m_v * (m_v + omega_2) * sigmabar + q2 * (5.0 + sigmabar * (5.0 - 6.0 * sigmabar))))))
                             / (m_B * power_of<2>(omega_2) * power_of<7>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2B_fpm_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return  pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-6.0) *
                   (4.0 * (-(6.0 * m_B3 * omega_1 * sigmabar4) + 10.0 * m_B2 * sigma4 * q2 +
                    m_v * omega_1 * (m_v2 * (5.0 - 20.0 * sigmabar) + q2 * sigmabar * (-37 + 4.0 * sigmabar)) +
                    m_B2 * m_v * sigmabar * (m_v * (2.0 - 5.0 * sigmabar) + 3.0 * omega_1 * (12.0 - 7.0 * sigmabar)) -
                    2.0 * m_B * sigma3 * (5.0 * (2.0 * m_v + omega_1) * q2 + 24 * m_B2 * m_v * sigmabar +
                    m_B * (5.0 * m_v2 - 7.0 * q2 * sigmabar)) +
                    2.0 * m_B * omega_1 * (q2 * sigmabar * (1 + 2.0 * sigmabar) + m_v2 * (-5.0 + sigmabar * (10.0 + sigmabar))) +
                    sigma * (12.0 * m_B4 * sigmabar4 + m_v * omega_1 * (-(5.0 * (4.0 * m_v2 + q2)) + 24 * q2 * sigmabar) -
                    3.0 * m_B3 * sigmabar * (4.0 * omega_1 * sigmabar2 + m_v * (12.0 - 7.0 * sigmabar)) -
                    2.0 * m_B2 * (6.0 * m_v * omega_1 * sigmabar * (7.0 - 2.0 * sigmabar) +
                    2.0 * q2 * sigmabar * (1 + 2.0 * sigmabar) + m_v2 * (-5.0 + 2.0 * sigmabar * (5.0 + sigmabar))) +
                    m_B * (12.0 * m_v2 * omega_1 * sigmabar + m_v * q2 * sigmabar * (37 - 4.0 * sigmabar) +
                    5.0 * m_v3 * (-1 + 4.0 * sigmabar) - 2.0 * omega_1 * q2 * (-5.0 + sigmabar * (-7.0 + sigmabar)))) +
                    sigma2 * (12.0 * m_B4 * sigmabar3 + 20.0 * m_v * omega_1 * q2 +
                    12.0 * m_B3 * m_v * sigmabar * (7.0 - 2.0 * sigmabar) +
                    m_B * (20.0 * m_v3 + 10.0 * m_v2 * omega_1 + m_v * q2 * (5.0 - 24 * sigmabar) -
                    12.0 * omega_1 * q2 * sigmabar) + 2.0 * m_B2 *
                   (m_v * (-(7.0 * m_v) + 24 * omega_1) * sigmabar + q2 * (-5.0 + sigmabar * (-7.0 + 2.0 * sigmabar))))) * chi_bar_bar_4 +
                    2.0 * m_B * (-omega_1 + m_B * sigma) * sigmabar *
                   (6.0 * m_B3 * sigmabar3 * sigma + 3.0 * m_B2 * m_v * sigmabar2 * (-3.0 + 4.0 * sigma) +
                    m_v3 * (-1 + 4.0 * sigma + 4.0 * sigmabar) +
                    m_v * q2 * (-(4.0 * sigma2) + sigma + 9.0 * sigmabar - 4.0 * sigma * sigmabar) +
                    2.0 * m_B * q2 * sigma * (-1 + sigma2 + (-2.0 + sigma) * sigmabar) -
                    m_B * m_v2 * (-2.0 + 5.0 * sigmabar + 2.0 * sigma * (sigma + sigmabar))) * chi_bar_4);
        }

        double I4d2C_fpm_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 2.0 * m_B * (6.0 * m_B3 * sigmabar3 * sigma +
                             (-m_v2 + q2 * sigma) * (2.0 * m_B * (-1 + sigma2) + m_v * (-1 + 4.0 * sigma)) +
                               sigmabar * (-(4.0 * m_v3) + 2.0 * m_B * q2 * sigma * (-2.0 + sigma) - m_B * m_v2 * (5.0 + 2.0 * sigma) +
                               m_v * q2 * (-9.0 + 4.0 * sigma) - 3.0 * m_B2 * m_v * (-3.0 + 4.0 * sigma) * sigmabar))
                             / (power_of<2>(omega_2) * power_of<5>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2D_fpm_3pt_chiB_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_± proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return    2.0 * pow(m_B,-1) * pow(sigmabar,-6.0) * ((3.0 * m_B2 * m_v * sigmabar2 * (7.0 - 8.0 * sigma) +
                      6.0 * m_B3 * sigmabar3 * (-(2.0 * sigma) + sigmabar) +
                      2.0 * m_B * q2 * (sigma3 - sigma + sigmabar2 * (-2.0 + sigma) + (2.0 * sigma2 - sigmabar) * sigmabar) +
                      m_v * (-1 + 4.0 * sigma + 4.0 * sigmabar) * (m_v2 - q2 * (sigma + sigmabar))) * chi_bar_bar_4 +
                      m_B * sigmabar * (6.0 * m_B3 * sigmabar3 * sigma +
                     (m_v + 2.0 * m_B * (-1 + sigma2) - 4.0 * m_v * sigma) * (-m_v2 + q2 * sigma) -
                      sigmabar * (-(4.0 * m_v3) - 2.0 * m_B * q2 * sigma * (-2.0 + sigma) + m_B * m_v2 * (5.0 + 2.0 * sigma) +
                      m_v * q2 * (-9.0 + 4.0 * sigma) - 3.0 * m_B2 * m_v * (-3.0 + 4.0 * sigma) * sigmabar)) * chi_bar_4);
        }

        double I3_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + 0.0                                                       - I3_fpm_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1A_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + 0.0                                                       - I3d1A_fpm_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1B_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return + 0.0                                                       - I3d1B_fpm_3pt_chiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I3d1C_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return + 0.0                                                       - I3d1C_fpm_3pt_chiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I4_fpm_3pt_chiA_bar_bar_4( sigma, omega_1, omega_2, q2)    - I4_fpm_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1A_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I4d1A_fpm_3pt_chiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d1A_fpm_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1B_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return + I4d1B_fpm_3pt_chiA_bar_bar_4( sigma, omega_1, q2)          - I4d1B_fpm_3pt_chiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d1C_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return + I4d1C_fpm_3pt_chiA_bar_bar_4( sigma, omega_2, q2)          - I4d1C_fpm_3pt_chiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2A_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I4d2A_fpm_3pt_chiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d2A_fpm_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d2B_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return + I4d2B_fpm_3pt_chiA_bar_bar_4( sigma, omega_1, q2)          - I4d2B_fpm_3pt_chiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d2C_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return + I4d2C_fpm_3pt_chiA_bar_bar_4( sigma, omega_2, q2)          - I4d2C_fpm_3pt_chiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2D_fpm_3pt_chi_bar_bar_4(const double & sigma, const double & q2) const
        {
            return + I4d2D_fpm_3pt_chiA_bar_bar_4( sigma, q2)                   - I4d2D_fpm_3pt_chiB_bar_bar_4(sigma, q2);
        }
        // }}}

        /* f_± : integrands and surface terms */
        // {{{
        double integrand_fpm_2pt_disp(const double & sigma, const double & q2) const
        {
            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);
            const double etad3    = 24.0 * (eta - 1.0) * power_of<2>(eta) * (2.0 * eta - 1.0) / power_of<3>(sigmabar);

            const double I1   = I1_fpm_2pt_phi_p(sigma, q2);
            const double I2   = I2_fpm_2pt_phi_bar(sigma, q2)   + I2_fpm_2pt_g_p(sigma, q2);
            const double I2d1 = I2d1_fpm_2pt_phi_bar(sigma, q2) + I2d1_fpm_2pt_g_p(sigma, q2);
            const double I3   = I3_fpm_2pt_g_p(sigma, q2)       + I3_fpm_2pt_g_bar(sigma, q2);
            const double I3d1 = I3d1_fpm_2pt_g_p(sigma, q2)     + I3d1_fpm_2pt_g_bar(sigma, q2);
            const double I3d2 = I3d2_fpm_2pt_g_p(sigma, q2)     + I3d2_fpm_2pt_g_bar(sigma, q2);
            const double I4   = I4_fpm_2pt_g_bar(sigma, q2);
            const double I4d1 = I4d1_fpm_2pt_g_bar(sigma, q2);
            const double I4d2 = I4d2_fpm_2pt_g_bar(sigma, q2);
            const double I4d3 = I4d3_fpm_2pt_g_bar(sigma, q2);

            double result = 0.0;
            result += -1.0 * I1;
            result += (etad1 * I2 + eta * I2d1) / m_B2;
            result += -1.0 * (I3 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I3d1 * eta * etad1 + I3d2 * power_of<2>(eta)) / (2.0 * m_B4);
            result += I4 * (power_of<2>(eta) * etad3 + 4.0 * eta * etad1 * etad2 + power_of<3>(etad1)) / (6.0 * m_B6);
            result += I4d1 * eta * (4.0 * eta * etad2 + 7.0 * power_of<2>(etad1)) / (6.0 * m_B6);
            result += I4d2 * 6.0 * power_of<2>(eta) * etad1 / (6.0 * m_B6);
            result += I4d3 * power_of<3>(eta) / (6.0 * m_B6);
            result *= exp;

            return result;
        }

        double integrand_fpm_2pt_borel(const double & sigma, const double & q2) const
        {
            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double I1   = I1_fpm_2pt_phi_p(sigma, q2);
            const double I2   = I2_fpm_2pt_phi_bar(sigma, q2)   + I2_fpm_2pt_g_p(sigma, q2);
            const double I3   = I3_fpm_2pt_g_p(sigma, q2)       + I3_fpm_2pt_g_bar(sigma, q2);
            const double I4   = I4_fpm_2pt_g_bar(sigma, q2);

            double result = 0.0;
            result += - I1;
            result +=   I2 / M2;
            result += - I3 / (2 * M4);
            result +=   I4 / (6 * M6);
            result *= exp;

            return result;
        }

        double surface_fpm_2pt(const double & sigma, const double & q2) const
        {
            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fpm_2pt_phi_bar(sigma, q2)   + I2_fpm_2pt_g_p(sigma, q2);
            const double I3   = I3_fpm_2pt_g_p(sigma, q2)       + I3_fpm_2pt_g_bar(sigma, q2);
            const double I3d1 = I3d1_fpm_2pt_g_p(sigma, q2)     + I3d1_fpm_2pt_g_bar(sigma, q2);
            const double I4   = I4_fpm_2pt_g_bar(sigma, q2);
            const double I4d1 = I4d1_fpm_2pt_g_bar(sigma, q2);
            const double I4d2 = I4d2_fpm_2pt_g_bar(sigma, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / (power_of<2>( M2())));
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;

            return result;
        }

        /*
         * rewrite integration ranges such that:
         * 1.)
         *    0 <= x_1 <= 1,     and    0 <= x_2 <= 1,
         * 2.)
         *    x_1 and x_2 integration boundaries do not depend on the other variables
         *
         * We obtain the integrand
         *
         *    sigma m_B f(sigma m_B x_1, sigma m_B (xbar_1 xbar_2 + x_2) / xbar_2) / (xbar_1 xbar_2^2 + x_2 xbar_2),
         *
         * where
         *
         *    xbar_1 = 1.0 - x_1,    and    xbar_2 = 1.0 - x_2.
         */
        double integrand_fpm_3pt(const std::array<double, 3> & args, const double & q2) const
        {
            const double sigma  = args[0];
            const double x_1    = args[1];
            const double x_2    = args[2];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            constexpr double I1 = 0.0;
            const     double I2 = I2_fpm_3pt_phi_3(sigma, omega_1, omega_2, q2)         + I2_fpm_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                                + I2_fpm_3pt_phi_4(sigma, omega_1, omega_2, q2)         + I2_fpm_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)
                                + I2_fpm_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I2_fpm_3pt_chi_bar_4(sigma, omega_1, omega_2, q2);
            const     double I3 = I3_fpm_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                                + I3_fpm_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                                + I3_fpm_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fpm_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                                + I3_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const     double I4 = I4_fpm_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                                + I4_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result = 0.0;
            result += - I1;
            result +=   I2 / M2;
            result += - I3 / (2.0 * M4);
            result +=   I4 / (6.0 * M6);
            result *=   prefactor * exp;

            return result;
        }

        double surface_fpm_3pt_A(const std::array<double, 2> & args, const double & sigma, const double & q2) const
        {
            const double x_1    = args[0];
            const double x_2    = args[1];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fpm_3pt_phi_3(sigma, omega_1, omega_2, q2)            + I2_fpm_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I2_fpm_3pt_phi_4(sigma, omega_1, omega_2, q2)            + I2_fpm_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)
                              + I2_fpm_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I2_fpm_3pt_chi_bar_4(sigma, omega_1, omega_2, q2);
            const double I3   = I3_fpm_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I3_fpm_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fpm_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fpm_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I3_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I3d1 = I3d1A_fpm_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I3d1A_fpm_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fpm_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fpm_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3d1A_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4   = I4_fpm_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)    + I4_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I4_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d1 = I4d1A_fpm_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d1A_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d1A_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d1A_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d2 = I4d2A_fpm_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d2A_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d2A_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d2A_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        double surface_fpm_3pt_B(const double & x_1, const double & sigma, const double & q2) const
        {
            // this ONLY includes the Jacobian from the transformation (omega_1 -> x_1).
            const double prefactor = sigma * m_B();

            const double omega_1 = sigma * m_B() * x_1;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1B_fpm_3pt_phi_bar_3(sigma, omega_1, q2)
                                  + I3d1B_fpm_3pt_phi_bar_4(sigma, omega_1, q2)     + I3d1B_fpm_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fpm_3pt_psi_bar_4(sigma, omega_1, q2)     + I3d1B_fpm_3pt_chi_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fpm_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I3d1B_fpm_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1B_fpm_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d1B_fpm_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d1B_fpm_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d1B_fpm_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            const     double I4d2 = I4d2B_fpm_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d2B_fpm_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d2B_fpm_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d2B_fpm_3pt_chi_bar_bar_4(sigma, omega_1, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        double surface_fpm_3pt_C(const double & x_2, const double & sigma, const double & q2) const
        {
            const double xbar_2 = 1.0 - x_2;

            // this ONLY includes the Jacobian from the transformation (omega_2 -> x_2).
            const double prefactor = sigma * m_B() / (xbar_2 * xbar_2);

            const double omega_2 = sigma * m_B() * (x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1C_fpm_3pt_phi_bar_3(sigma, omega_2, q2)
                                  + I3d1C_fpm_3pt_phi_bar_4(sigma, omega_2, q2)     + I3d1C_fpm_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fpm_3pt_psi_bar_4(sigma, omega_2, q2)     + I3d1C_fpm_3pt_chi_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fpm_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I3d1C_fpm_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1C_fpm_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d1C_fpm_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d1C_fpm_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d1C_fpm_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            const     double I4d2 = I4d2C_fpm_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d2C_fpm_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d2C_fpm_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d2C_fpm_3pt_chi_bar_bar_4(sigma, omega_2, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        double surface_fpm_3pt_D(const double & sigma, const double & q2) const
        {
            // this does NOT includes the original factor of 1 / omega_2
            const double prefactor = 1.0;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);


            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            constexpr double I3d1 = 0.0;
            constexpr double I4   = 0.0;
            constexpr double I4d1 = 0.0;
            const     double I4d2 = I4d2D_fpm_3pt_phi_bar_bar_3(sigma, q2) + I4d2D_fpm_3pt_phi_bar_bar_4(sigma, q2)
                                  + I4d2D_fpm_3pt_psi_bar_bar_4(sigma, q2) + I4d2D_fpm_3pt_chi_bar_bar_4(sigma, q2);


            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        /*
         * Integrands for the first moments. Only the borel method is implemented
         */

        double integrand_fpm_2pt_borel_m1(const double & sigma, const double & q2) const
        {
            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double I1   = I1_fpm_2pt_phi_p(sigma, q2);
            const double I2   = I2_fpm_2pt_phi_bar(sigma, q2)   + I2_fpm_2pt_g_p(sigma, q2);
            const double I3   = I3_fpm_2pt_g_p(sigma, q2)       + I3_fpm_2pt_g_bar(sigma, q2);
            const double I4   = I4_fpm_2pt_g_bar(sigma, q2);

            double result1 = 0.0;
            result1 += - I1;
            result1 +=   I2 / M2;
            result1 += - I3 / (2.0 * M4);
            result1 +=   I4 / (6.0 * M6);
            result1 *=   exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - I2;
            result2 +=   I3 / M2;
            result2 += - I4 / (2.0 * M4);
            result2 *=   exp;

            return result1 + result2;
        }

        double surface_fpm_2pt_m1(const double & sigma, const double & q2) const
        {
            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double m_P2 = power_of<2>(m_P());
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fpm_2pt_phi_bar(sigma, q2)   + I2_fpm_2pt_g_p(sigma, q2);
            const double I3   = I3_fpm_2pt_g_p(sigma, q2)       + I3_fpm_2pt_g_bar(sigma, q2);
            const double I3d1 = I3d1_fpm_2pt_g_p(sigma, q2)     + I3d1_fpm_2pt_g_bar(sigma, q2);
            const double I4   = I4_fpm_2pt_g_bar(sigma, q2);
            const double I4d1 = I4d1_fpm_2pt_g_bar(sigma, q2);
            const double I4d2 = I4d2_fpm_2pt_g_bar(sigma, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return result1 + result2;
        }
        double integrand_fpm_3pt_m1(const std::array<double, 3> & args, const double & q2) const
        {
            const double sigma  = args[0];
            const double x_1    = args[1];
            const double x_2    = args[2];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double I1 = 0.0;
            const double I2 = I2_fpm_3pt_phi_3(sigma, omega_1, omega_2, q2)         + I2_fpm_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                            + I2_fpm_3pt_phi_4(sigma, omega_1, omega_2, q2)         + I2_fpm_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)
                            + I2_fpm_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I2_fpm_3pt_chi_bar_4(sigma, omega_1, omega_2, q2);
            const double I3 = I3_fpm_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                            + I3_fpm_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                            + I3_fpm_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fpm_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                            + I3_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4 = I4_fpm_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                            + I4_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result1 = 0.0;
            result1 += - I1;
            result1 +=   I2 / M2;
            result1 += - I3 / (2.0 * M4);
            result1 +=   I4 / (6.0 * M6);
            result1 *=   exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - I2;
            result2 +=   I3 / M2;
            result2 += - I4 / (2.0 * M4);
            result2 *=   exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fpm_3pt_A_m1(const std::array<double, 2> & args, const double & sigma, const double & q2) const
        {
            const double x_1    = args[0];
            const double x_2    = args[1];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fpm_3pt_phi_3(sigma, omega_1, omega_2, q2)            + I2_fpm_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I2_fpm_3pt_phi_4(sigma, omega_1, omega_2, q2)            + I2_fpm_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)
                              + I2_fpm_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I2_fpm_3pt_chi_bar_4(sigma, omega_1, omega_2, q2);
            const double I3   = I3_fpm_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I3_fpm_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fpm_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fpm_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I3_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I3d1 = I3d1A_fpm_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I3d1A_fpm_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fpm_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fpm_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3d1A_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4   = I4_fpm_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)    + I4_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I4_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d1 = I4d1A_fpm_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d1A_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d1A_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d1A_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d2 = I4d2A_fpm_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d2A_fpm_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d2A_fpm_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d2A_fpm_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fpm_3pt_B_m1(const double & x_1, const double & sigma, const double & q2) const
        {
            // this ONLY includes the Jacobian from the transformation (omega_1 -> x_1).
            const double prefactor = sigma * m_B();

            const double omega_1 = sigma * m_B() * x_1;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1B_fpm_3pt_phi_bar_3(sigma, omega_1, q2)
                                  + I3d1B_fpm_3pt_phi_bar_4(sigma, omega_1, q2)     + I3d1B_fpm_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fpm_3pt_psi_bar_4(sigma, omega_1, q2)     + I3d1B_fpm_3pt_chi_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fpm_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I3d1B_fpm_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1B_fpm_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d1B_fpm_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d1B_fpm_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d1B_fpm_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            const     double I4d2 = I4d2B_fpm_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d2B_fpm_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d2B_fpm_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d2B_fpm_3pt_chi_bar_bar_4(sigma, omega_1, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fpm_3pt_C_m1(const double & x_2, const double & sigma, const double & q2) const
        {
            const double xbar_2 = 1.0 - x_2;

            // this ONLY includes the Jacobian from the transformation (omega_2 -> x_2).
            const double prefactor = sigma * m_B() / (xbar_2 * xbar_2);

            const double omega_2 = sigma * m_B() * (x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1C_fpm_3pt_phi_bar_3(sigma, omega_2, q2)
                                  + I3d1C_fpm_3pt_phi_bar_4(sigma, omega_2, q2)     + I3d1C_fpm_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fpm_3pt_psi_bar_4(sigma, omega_2, q2)     + I3d1C_fpm_3pt_chi_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fpm_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I3d1C_fpm_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1C_fpm_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d1C_fpm_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d1C_fpm_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d1C_fpm_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            const     double I4d2 = I4d2C_fpm_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d2C_fpm_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d2C_fpm_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d2C_fpm_3pt_chi_bar_bar_4(sigma, omega_2, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fpm_3pt_D_m1(const double & sigma, const double & q2) const
        {
            // this does NOT includes the original factor of 1 / omega_2
            const double prefactor = 1.0;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);


            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            constexpr double I3d1 = 0.0;
            constexpr double I4   = 0.0;
            constexpr double I4d1 = 0.0;
            const     double I4d2 = I4d2D_fpm_3pt_phi_bar_bar_3(sigma, q2) + I4d2D_fpm_3pt_phi_bar_bar_4(sigma, q2)
                                  + I4d2D_fpm_3pt_psi_bar_bar_4(sigma, q2) + I4d2D_fpm_3pt_chi_bar_bar_4(sigma, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }
        // }}}

        /* f_+ : form factor and moments */
        // {{{
        double f_pm(const double & q2) const
        {
            const double sigma_0 = this->sigma_0(q2, s0_0_pm(), s0_1_pm());

            const std::function<double (const double &)> integrand_2pt = std::bind(integrand_fpm_2pt, this, std::placeholders::_1, q2);

            const double integral_2pt = integrate<GSL::QAGS>(integrand_2pt, 0.0, sigma_0);
            const double surface_2pt  = 0.0 - surface_fpm_2pt(switch_borel ? sigma_0 : 0.0, q2);

            double integral_3pt = 0.0;
            double surface_3pt  = 0.0;

            if (switch_3pt != 0.0)
            {
                const std::function<double (const double &)> surface_3pt_B = std::bind(&Implementation::surface_fpm_3pt_B, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const double &)> surface_3pt_C = std::bind(&Implementation::surface_fpm_3pt_C, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const std::array<double, 3> &)> integrand_3pt = std::bind(&Implementation::integrand_fpm_3pt, this, std::placeholders::_1, q2);
                const std::function<double (const std::array<double, 2> &)> surface_3pt_A = std::bind(&Implementation::surface_fpm_3pt_A, this, std::placeholders::_1, sigma_0, q2);

                integral_3pt = integrate(integrand_3pt, { 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0 }, cubature::Config());
                surface_3pt  = 0.0
                             - integrate(surface_3pt_A, { 0.0, 0.0 }, { 1.0, 1.0 }, cubature::Config()) // integrate over x_1 and x_2
                             - integrate<GSL::QAGS>(surface_3pt_B, 0.0, 1.0)                            // integrate over x_1
                             - integrate<GSL::QAGS>(surface_3pt_C, 0.0, 1.0)                            // integrate over x_2
                             - surface_fpm_3pt_D(sigma_0, q2);
            }

            return f_B() * m_B() / f_P() * (integral_2pt + surface_2pt + integral_3pt + surface_3pt) / ( Traits::chi2);
        }

        double normalized_moment_1_f_pm(const double & q2) const
        {
            const double sigma_0 = this->sigma_0(q2, s0_0_pm(), s0_1_pm());

            const std::function<double (const double &)> integrand_2pt_m1 = std::bind(&Implementation::integrand_fpm_2pt_borel_m1, this, std::placeholders::_1, q2);


            const std::function<double (const double &)> integrand_2pt    = std::bind(&Implementation::integrand_fpm_2pt_borel, this, std::placeholders::_1, q2);

            const double integral_2pt_m1 = integrate<GSL::QAGS>(integrand_2pt_m1, 0.0, sigma_0);
            const double surface_2pt_m1  = 0.0 - surface_fpm_2pt_m1(sigma_0, q2);

            double integral_3pt_m1 = 0.0;
            double surface_3pt_m1  = 0.0;

            if (switch_3pt != 0.0)
            {
                const std::function<double (const double &)> surface_3pt_B_m1 = std::bind(&Implementation::surface_fpm_3pt_B_m1, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const double &)> surface_3pt_C_m1 = std::bind(&Implementation::surface_fpm_3pt_C_m1, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const std::array<double, 3> &)> integrand_3pt_m1 = std::bind(&Implementation::integrand_fpm_3pt_m1, this, std::placeholders::_1, q2);
                const std::function<double (const std::array<double, 2> &)> surface_3pt_A_m1 = std::bind(&Implementation::surface_fpm_3pt_A_m1, this, std::placeholders::_1, sigma_0, q2);

                integral_3pt_m1 = integrate(integrand_3pt_m1, { 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0 }, cubature::Config());
                surface_3pt_m1  = 0.0
                                - integrate(surface_3pt_A_m1, { 0.0, 0.0 }, { 1.0, 1.0 }, cubature::Config()) // integrate over x_1 and x_2
                                - integrate<GSL::QAGS>(surface_3pt_B_m1, 0.0, 1.0)                            // integrate over x_1
                                - integrate<GSL::QAGS>(surface_3pt_C_m1, 0.0, 1.0)                            // integrate over x_2
                                - surface_fpm_3pt_D_m1(sigma_0, q2);
            }
            const double numerator       = integral_2pt_m1 + surface_2pt_m1 + integral_3pt_m1 + surface_3pt_m1;

            const double integral_2pt    = integrate<GSL::QAGS>(integrand_2pt, 0.0, sigma_0);
            const double surface_2pt     = 0.0 - surface_fpm_2pt(sigma_0, q2);

            double integral_3pt    = 0.0;
            double surface_3pt     = 0.0;

            if (switch_3pt != 0.0)
            {
                const std::function<double (const double &)> surface_3pt_B    = std::bind(&Implementation::surface_fpm_3pt_B, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const double &)> surface_3pt_C    = std::bind(&Implementation::surface_fpm_3pt_C, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const std::array<double, 3> &)> integrand_3pt = std::bind(&Implementation::integrand_fpm_3pt, this, std::placeholders::_1, q2);
                const std::function<double (const std::array<double, 2> &)> surface_3pt_A = std::bind(&Implementation::surface_fpm_3pt_A, this, std::placeholders::_1, sigma_0, q2);

                integral_3pt    = integrate(integrand_3pt, { 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0 }, cubature::Config());
                surface_3pt     = 0.0
                                - integrate(surface_3pt_A, { 0.0, 0.0 }, { 1.0, 1.0 }, cubature::Config()) // integrate over x_1 and x_2
                                - integrate<GSL::QAGS>(surface_3pt_B, 0.0, 1.0)                            // integrate over x_1
                                - integrate<GSL::QAGS>(surface_3pt_C, 0.0, 1.0)                            // integrate over x_2
                                - surface_fpm_3pt_D(sigma_0, q2);
            }
            const double denominator     = integral_2pt + surface_2pt + integral_3pt + surface_3pt;

            return numerator / denominator;
        }
        // }}}

        /* f_T */

        inline
        double I1_fT_2pt_phi_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_T proportional to phibar
            const double sigmabar = 1.0 - sigma;

            const double phi_bar  = this->phi_bar(sigma * m_B);

            const double C_1 = 1.0 / (sigmabar * m_B);

            return C_1 * phi_bar;
        }

        inline
        double I2_fT_2pt_phi_bar(const double & sigma, const double & q2) const
        {
            // two-particle contribution to f_T proportional to phibar
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double phi_bar  = this->phi_bar(sigma * m_B);

            const double C_2 = (-m_B2 * sigmabar2 + m_v2 + 2.0 * q2 * sigmabar - q2) / (power_of<2>(sigmabar) * m_B);

            return C_2 * phi_bar;
        }

        inline
        double I2d1_fT_2pt_phi_bar(const double & sigma, const double & q2) const
        {
            // first derivative of two-particle contribution to f_T proportional to phibar
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double phi_bar     = this->phi_bar(sigma * m_B);
            const double phi_bar_d1  = this->phi_bar_d1(sigma * m_B);

            const double C_2   = (-m_B2 * sigmabar2 + m_v2 + 2.0 * q2 * sigmabar - q2) / power_of<2>(sigmabar);
            const double C_2d1 = 2.0 * (m_v2 - q2 * sigma) / (power_of<3>(sigmabar) * m_B);

            return C_2 * phi_bar_d1 + C_2d1 * phi_bar;
        }

        inline
        double I2_fT_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_T proportional to gbar
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B);

            const double C_2 = 8.0 / (power_of<2>(sigmabar) * m_B);

            return C_2 * g_bar;
        }

        inline
        double I2d1_fT_2pt_g_bar(const double & sigma, const double & /*q2*/) const
        {
            // two-particle contribution to f_T proportional to gbar
            const double sigmabar = 1.0 - sigma;

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);

            const double C_2   =  8.0 / power_of<2>(sigmabar);
            const double C_2d1 =  16.0 / (power_of<3>(sigmabar) * m_B);

            return C_2 * g_bar_d1 + C_2d1 * g_bar;
        }

        inline
        double I3_fT_2pt_g_bar(const double & sigma, const double & q2) const
        {
            // two-particle contribution to f_T proportional to gbar
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);

            const double C_3 = -8.0 * (m_B2 * sigmabar2 + 2.0 * m_v2 - 2.0 * q2 * sigmabar + q2) / (sigmabar3 * m_B);

            return C_3 * g_bar;
        }

        inline
        double I3d1_fT_2pt_g_bar(const double & sigma, const double & q2) const
        {
            // two-particle contribution to f_T proportional to gbar
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);

            const double C_3   = -8.0 * (m_B2 * sigmabar2 + 2.0 * m_v2 - 2.0 * q2 * sigmabar + q2) / sigmabar3;
            const double C_3d1 = -8.0 * (m_B2 * sigmabar2 + 6.0 * m_v2 + q2 * (3.0 - 4.0 * sigmabar)) / (power_of<4>(sigmabar) * m_B);

            return C_3 * g_bar_d1 + C_3d1 * g_bar;
        }

        inline
        double I3d2_fT_2pt_g_bar(const double & sigma, const double & q2) const
        {
            // two-particle contribution to f_T proportional to gbar
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);
            const double g_bar_d2 = this->g_bar_d2(sigma * m_B);

            const double C_3   = -8.0 * (m_B2 * sigmabar2 + 2.0 * m_v2 - 2.0 * q2 * sigmabar + q2) / sigmabar3 * m_B;
            const double C_3d1 = -16.0 * (m_B2 * sigmabar2 + 6.0 * m_v2 + q2 * (3.0 - 4.0 * sigmabar)) / (power_of<4>(sigmabar));
            const double C_3d2 = -16.0 * (m_B2 * sigmabar2 + 12.0 * m_v2 + 6.0 * q2 * sigma) / (power_of<5>(sigmabar) * m_B);

            return C_3 * g_bar_d2 + C_3d1 * g_bar_d1 + C_3d2 * g_bar;
        }

        inline
        double I4_fT_2pt_g_bar(const double & sigma, const double & q2) const
        {
            // two-particle contribution to f_T proportional to gbar
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar4 = power_of<4>(sigmabar);
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);

            const double C_4 = 24.0 * m_v2 * (m_B2 * sigmabar2 - m_v2 - 2.0 * q2 * sigmabar + q2) / (sigmabar4 * m_B);

            return C_4 * g_bar;
        }

        inline
        double I4d1_fT_2pt_g_bar(const double & sigma, const double & q2) const
        {
            // two-particle contribution to f_T proportional to gbar
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar4 = power_of<4>(sigmabar), sigmabar5 = power_of<5>(sigmabar);
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);

            const double C_4   = 24.0 * m_v2 * (m_B2 * sigmabar2 - m_v2 - 2.0 * q2 * sigmabar + q2) / sigmabar4;
            const double C_4d1 = 48.0 * m_v2 * (m_B2 * sigmabar2 - 2.0 * m_v2 + q2 * (2.0 - 3.0 * sigmabar)) / (sigmabar5 * m_B);
            return C_4 * g_bar_d1 + C_4d1 * g_bar;
        }

        inline
        double I4d2_fT_2pt_g_bar(const double & sigma, const double & q2) const
        {
            // two-particle contribution to f_T proportional to gbar
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar4 = power_of<4>(sigmabar), sigmabar5 = power_of<5>(sigmabar);
            const double sigmabar6 = power_of<6>(sigmabar);
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);
            const double g_bar_d2 = this->g_bar_d2(sigma * m_B);

            const double C_4   = 24.0 * m_v2 * (m_B2 * sigmabar2 - m_v2 - 2.0 * q2 * sigmabar + q2) / sigmabar4 * m_B;
            const double C_4d1 = 96.0 * m_v2 * (m_B2 * sigmabar2 - 2.0 * m_v2 + q2 * (2.0 - 3.0 * sigmabar)) / sigmabar5;
            const double C_4d2 = 48.0 * m_v2 * (3.0 * m_B2 * sigmabar2 - 10.0 * m_v2 - 2.0 * q2 * (6.0 * sigmabar - 5.0)) / (sigmabar6 * m_B);

            return C_4 * g_bar_d2 + C_4d1 * g_bar_d1 + C_4d2 * g_bar;
        }

        inline
        double I4d3_fT_2pt_g_bar(const double & sigma, const double & q2) const
        {
            // two-particle contribution to f_T proportional to gbar
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar4 = power_of<4>(sigmabar), sigmabar5 = power_of<5>(sigmabar);
            const double sigmabar6 = power_of<6>(sigmabar);
            const double m_B2 = power_of<2>(m_B);
            const double m_v = this->m_v(), m_v2 = power_of<2>(m_v);

            const double g_bar    = this->g_bar(sigma * m_B);
            const double g_bar_d1 = this->g_bar_d1(sigma * m_B);
            const double g_bar_d2 = this->g_bar_d2(sigma * m_B);
            const double g_bar_d3 = this->g_bar_d3(sigma * m_B);

            const double C_4   = 24.0 * m_v2 * (m_B2 * sigmabar2 - m_v2 - 2.0 * q2 * sigmabar + q2) / sigmabar4 * m_B2;
            const double C_4d1 = 144.0 * m_v2 * (m_B2 * sigmabar2 - 2.0 * m_v2 + q2 * (2.0 - 3.0 * sigmabar)) / sigmabar5 * m_B;
            const double C_4d2 = 24.0 * m_v2 * (18.0 * m_B2 * sigmabar2 - 60.0 * m_v2 + q2 * (60.0 - 72.0 * sigmabar)) / sigmabar6;
            const double C_4d3 = 576.0 * m_v2 * (m_B2 * sigmabar2 - 5.0 * m_v2 + 5.0 * q2 * sigma) / (power_of<7>(sigmabar) * m_B);

            return C_4 * g_bar_d3 + C_4d1 * g_bar_d2 + C_4d2 * g_bar_d1 + C_4d3 * g_bar;
        }

        /* f_T : 3-particle functions */
        // {{{
        double I1_fT_3pt_phi_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_3
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_3 = this->phi_3(omega_1, omega_2);

            const double C_1 = 2.0 * u / (m_B * m_B * power_of<2>(sigmabar));

            return C_1 * phi_3;
        }

        double I2_fT_3pt_phi_3(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_3
            const double m_B2 = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_3 = this->phi_3(omega_1, omega_2);

            const double C_2 = -2.0 * u * (-m_v2 + q2 - 2.0 * q2 * sigmabar + m_B2 * sigmabar2) / (m_B2 * power_of<3>(sigmabar));

            return C_2 * phi_3;
        }

        double I2_fT_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_bar_3
            const double m_B2 = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_2 = 4.0 * (m_v + m_B * u * sigmabar) / (m_B2 * power_of<3>(sigmabar));

            return C_2 * phi_bar_3;
        }

        double I3_fT_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_3
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = - 4.0 * (m_v + m_B * u * sigmabar) * (-m_v2 + q2 - 2.0 * q2 * sigmabar + m_B2 * sigmabar2) / (m_B2 * power_of<4>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I3d1A_fT_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_3
            const double m_B2     = power_of<2>(m_B),   m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(),   m_v2 = power_of<2>(m_v),   m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = 4.0 * (4.0 * m_v3 * omega_2 + m_B * m_v2 * sigmabar * (-(3.0 * omega_1) + 3.0 * m_B * sigma + m_B * sigmabar) +
                               m_v * omega_2 * (2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar) +
                               sigma * (-(4.0 * q2) + 3.0 * m_B2 * sigmabar)) +
                               m_B * sigmabar * (-(3.0 * m_B * sigma2 * q2) - 2.0 * m_B3 * sigmabar2 * sigma +
                               m_B * sigmabar2 * (q2 + m_B2 * (-1 + 2.0 * sigma)) +
                               omega_1 * (sigma * (3.0 * q2 - 2.0 * m_B2 * sigmabar) - sigmabar * (q2 + m_B2 * (-2.0 + sigmabar)))))
                             / (m_B2 * omega_2 * power_of<5>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I3d1B_fT_3pt_phi_bar_3(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_3
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = -4.0 * (m_v + m_B * sigmabar) * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B * (-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I3d1C_fT_3pt_phi_bar_3(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_3
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

            const double C_3 = 4.0 * m_v * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar)) / (m_B * omega_2 * power_of<4>(sigmabar));

            return C_3 * phi_bar_3;
        }

        double I3_fT_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_3 = 12.0 * m_v * u / (m_B * power_of<3>(sigmabar));

            return C_3 * phi_bar_bar_3;
        }

        double I3d1A_fT_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_3 = 12.0 * m_v * (3.0 * sigma * m_B + m_B * sigmabar - 3.0 * omega_1)
                               / (m_B * omega_2 * power_of<4>(sigmabar));

            return C_3 * phi_bar_bar_3;
        }

        double I3d1B_fT_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_3 = -12.0 * m_v / ((-omega_1 + m_B * sigma) * power_of<3>(sigmabar));

            return C_3 * phi_bar_bar_3;
        }

        double I3d1C_fT_3pt_phi_bar_bar_3(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3

            return 0.0;
        }

        double I4_fT_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = -12.0 * m_v * u * (-m_v2 + q2 - 2.0 * q2 * sigmabar + m_B2 * sigmabar2)
                             / (m_B * power_of<4>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d1A_fT_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = 12.0 * m_v * (-(omega_1 * (4.0 * m_v2 + 2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar) +
                               sigma * (-(4.0 * q2) + 3.0 * m_B2 * sigmabar))) +
                               m_B * (4.0 * sigma * (m_v2 - q2 * sigma) + sigmabar2 * (q2 + m_B2 * (-1 + 2.0 * sigma)) +
                               sigmabar * (m_v2 + sigma * (q2 - 3.0 * m_B2 * sigmabar))))
                             / (m_B * omega_2 * power_of<5>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d1B_fT_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = -12.0 * m_v * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar))
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * phi_bar_bar_3;
        }

        double I4d1C_fT_3pt_phi_bar_bar_3(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3

            return 0.0;
        }
        double I4d2A_fT_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_B2 = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = 24.0 * m_v * (2.0 * sigma2 * (-(5.0 * m_B * q2) + 3.0 * m_B3 * sigmabar) -
                               omega_1 * (10.0 * m_v2 + 2.0 * q2 * sigmabar + 3.0 * m_B2 * sigmabar * (-2.0 + sigmabar)) +
                               m_B * sigmabar * (4.0 * m_v2 + 2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar)) +
                               2.0 * sigma * (5.0 * omega_1 * q2 - 3.0 * m_B2 * omega_1 * sigmabar + 3.0 * m_B3 * (-1 + sigmabar) * sigmabar +
                               m_B * (5.0 * m_v2 - q2 * sigmabar)))
                             / (m_B * omega_2 * power_of<6>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d2B_fT_3pt_phi_bar_bar_3(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3
            const double omega_2  = m_B * sigma - omega_1;
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_B2 = power_of<2>(m_B);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3     = this->phi_bar_3(omega_1, omega_2);
            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            return  pow(-omega_1 + m_B * sigma,-1) * pow(sigmabar,-5.0) *
                    (-24.0 * m_v * (4.0 * m_v2 + 2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar) +
                    sigma * (-(4.0 * q2) + 3.0 * m_B2 * sigmabar)) * phi_bar_bar_3 -
                    12.0 * m_B * m_v * sigmabar * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * phi_bar_3);
        }

        double I4d2C_fT_3pt_phi_bar_bar_3(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2     = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            const double C_4 = 12.0 * m_B * m_v * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar)) / (power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_4 * phi_bar_bar_3;
        }

        double I4d2D_fT_3pt_phi_bar_bar_3(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_3
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_B2 = power_of<2>(m_B);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_3     = this->phi_bar_3(omega_1, omega_2);
            const double phi_bar_bar_3 = this->phi_bar_bar_3(omega_1, omega_2);

            return     12.0 * m_v * pow(m_B,-1) * pow(sigmabar,-4.0) *
                       ((2.0 * q2 - 2.0 * m_B2 * sigmabar) * phi_bar_bar_3 -
                        m_B * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * phi_bar_3);
        }

        double I2_fT_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_bar_4
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_2 = -2.0 / (m_B * power_of<2>(sigmabar));

            return C_2 * phi_bar_4;
        }

        double I3_fT_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(),  m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma,  sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (-m_v2 + q2 - 2.0 * q2 * sigmabar + m_B2 * sigmabar2)
                             / (m_B * power_of<3>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3d1A_fT_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (3.0 * m_v2 + (q2 + m_B2 * (-2.0 + sigmabar)) * sigmabar + sigma * (-3.0 * q2 + 2.0 * m_B2 * sigmabar))
                             / (m_B * power_of<4>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3d1B_fT_3pt_phi_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar))
                             / ((-omega_1 + m_B * sigma) * power_of<3>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I3d1C_fT_3pt_phi_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

            const double C_3 = -2.0 * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar))
                             / (omega_2 * power_of<3>(sigmabar));

            return C_3 * phi_bar_4;
        }

        double I2_fT_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double m_B2 = power_of<2>(m_B);
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_2 = -4.0 * u * (-1.0 + 2.0 * u) / (m_B2 * power_of<3>(sigmabar));

            return C_2 * phi_bar_bar_4;
        }

        double I3_fT_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * u * (-1.0 + 2.0 * u) * (m_v2 + q2 * (5.0 - 4.0 * sigmabar) - m_B2 * sigmabar2)
                             / (m_B2 * power_of<4>(sigmabar));

            return C_3 * phi_bar_bar_4;
        }

        double I3d1A_fT_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B),   m_B3 = power_of<3>(m_B),   m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(),   m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma), sigma4   = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (m_B3 * (4.0 * omega_1 + omega_2) * sigmabar2 +
                               2.0 * m_B2 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar * (-2.0 + sigmabar) -
                               8.0 * sigma4 * (5.0 * m_B2 * q2 + m_B4 * sigmabar) +
                               4.0 * m_B * sigma3 * (5.0 * (4.0 * omega_1 + omega_2) * q2 + 10.0 * m_B * (m_v2 + q2) +
                               m_B2 * (4.0 * omega_1 + omega_2) * sigmabar - 8.0 * m_B * q2 * sigmabar - 2.0 * m_B3 * sigmabar * (-2.0 + sigmabar)) -
                               m_B * (4.0 * omega_1 + omega_2) * sigmabar * (q2 * sigmabar + m_v2 * (-4.0 + 5.0 * sigmabar)) +
                               4.0 * omega_1 * (2.0 * omega_1 + omega_2) * (2.0 * q2 * sigmabar + m_v2 * (-5.0 + 6.0 * sigmabar)) -
                               2.0 * sigma * (2.0 * m_B4 * sigmabar2 + 2.0 * m_B3 * (4.0 * omega_1 + omega_2) * (-1 + sigmabar) * sigmabar -
                               2.0 * omega_1 * (2.0 * omega_1 + omega_2) * (5.0 * (m_v2 + q2) - 2.0 * q2 * sigmabar) +
                               m_B2 * sigmabar * (-(2.0 * q2 * sigmabar) + omega_1 * (2.0 * omega_1 + omega_2) * (-4.0 + sigmabar) +
                               m_v2 * (8.0 - 10.0 * sigmabar)) + 2.0 * m_B * (4.0 * omega_1 + omega_2) *
                              (3.0 * q2 * sigmabar + m_v2 * (-5.0 + 7.0 * sigmabar))) +
                               sigma2 * (-(20.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2) + 4.0 * m_B4 * sigmabar * (-2.0 + 3.0 * sigmabar) +
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar * (-8.0 + 3.0 * sigmabar) -
                               4.0 * m_B * (4.0 * omega_1 + omega_2) * (5.0 * (m_v2 + q2) - 3.0 * q2 * sigmabar) +
                               4.0 * m_B2 * (-(omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) + 8.0 * q2 * sigmabar +
                               2.0 * m_v2 * (-5.0 + 8.0 * sigmabar))))
                             / (m_B2 * power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_3 * phi_bar_bar_4;
        }

        double I3d1B_fT_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(),   m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (m_v2 * (4.0 - 4.0 * sigma - 5.0 * sigmabar) - 4.0 * sigma * q2 * sigmabar + sigmabar * (-q2 + m_B2 * sigmabar2))
                             / (m_B * (-omega_1 + m_B * sigma) * power_of<5>(sigmabar));

            return C_3 * phi_bar_bar_4;
        }

        double I3d1C_fT_3pt_phi_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4

            return 0.0;
        }

        double I4_fT_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B4     = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2     = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * u * (-1.0 + 2.0 * u) * (m_v4 - 2.0 * m_B2 * q2 * sigmabar3 + m_B4 * sigmabar4 + q2 * q2 * (-1.0 + 2.0 * sigmabar)
                             + 2.0 * m_v2 * sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B2 * power_of<5>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d1A_fT_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -(2.0 * (6.0 * (-omega_1 + m_B * sigma) * (2.0 * omega_1 + omega_2 - 2.0 * m_B * sigma) *
                                (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) + (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) *
                                sigmabar) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) -
                                2.0 * m_B * (-omega_1 + m_B * sigma) * sigmabar *
                                (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) + (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) *
                                sigmabar) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) +
                                m_B * (2.0 * omega_1 + omega_2 - 2.0 * m_B * sigma) * sigmabar *
                                (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) + (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) *
                                sigmabar) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) +
                                (-omega_1 + m_B * sigma) * (-(2.0 * omega_1) - omega_2 + 2.0 * m_B * sigma) * sigmabar *
                                (3.0 * m_v4 + 12.0 * m_v2 * q2 * sigmabar - 3.0 * m_B4 * sigmabar3 * (-2.0 + 2.0 * sigma - 3.0 * sigmabar) +
                                m_B2 * (sigmabar2 * (-m_v2 + q2 * sigma) + sigmabar2 * (-(5.0 * m_v2) + q2 * (-11.0 + 10.0 * sigma)) -
                                sigmabar2 * (12.0 * m_v2 + 13.0 * q2 * sigmabar)) +
                                (-(6.0 * sigma2) + sigmabar * (7.0 + 2.0 * sigmabar) + sigma * (3.0 - 4.0 * sigmabar)) * power_of<2>(q2))))
                             / (m_B2 * power_of<2>(omega_2) * power_of<7>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d1B_fT_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 2.0 * (3.0 * m_B2 * sigmabar3 - 2.0 * q2 * sigma * sigmabar + q2 * (-3.0 + 2.0 * sigma) * sigmabar +
                               m_v2 * (2.0 - 2.0 * sigma - 5.0 * sigmabar)) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B * power_of<6>(sigmabar) * omega_2);

            return C_4 * phi_bar_bar_4;
        }

        double I4d1C_fT_3pt_phi_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4

            return 0.0;
        }

        double I4d2A_fT_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B), m_B6 = power_of<6>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma), sigma5 = power_of<5>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 4.0 * (6.0 * m_B6 * sigmabar4 + 3.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar3 * (-4.0 + 3.0 * sigmabar) +
                               m_B4 * sigmabar2 * (-(12.0 * sigmabar2 * q2) - 2.0 * m_v2 * sigmabar * (1 + 5.0 * sigmabar) +
                               3.0 * omega_1 * (2.0 * omega_1 + omega_2) * (13.0 + 3.0 * sigmabar * (-5.0 + sigmabar))) -
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar2 *
                               (q2 * sigmabar * (-23 + 11.0 * sigmabar) + m_v2 * (-5.0 + sigmabar * (-18.0 + 5.0 * sigmabar))) +
                               m_B2 * sigmabar * (2.0 * m_v4 * sigmabar * (-2.0 + 5.0 * sigmabar) +
                               m_v2 * (2.0 * sigmabar2 * q2 * (1 + 5.0 * sigmabar) +
                               omega_1 * (2.0 * omega_1 + omega_2) * (-15.0 + sigmabar * (-40 + 19.0 * sigmabar))) +
                               q2 * sigmabar * (6.0 * sigmabar2 * q2 - omega_1 * (2.0 * omega_1 + omega_2) * (66 + sigmabar * (-53 + 5.0 * sigmabar)))) +
                               sigma5 * (-60 * m_B6 * sigmabar2 + 30 * m_B4 * q2 * sigmabar + 84 * m_B2 * power_of<2>(q2)) +
                               sigma3 * (6.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar2 * (-15.0 + 8.0 * sigmabar) +
                               15.0 * m_B3 * (4.0 * omega_1 + omega_2) * sigmabar * (m_v2 + 2.0 * q2 * (1 + sigmabar)) -
                               6.0 * m_B6 * sigmabar2 * (33 + 2.0 * sigmabar * (-24 + 5.0 * sigmabar)) +
                               6.0 * m_B * (4.0 * omega_1 + omega_2) * q2 * (14.0 * m_v2 + q2 * (7.0 - 8.0 * sigmabar)) +
                               2.0 * m_B4 * sigmabar * (-(15.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) + 30 * m_v2 * (1 + sigmabar) +
                               q2 * (15.0 + (81 - 70 * sigmabar) * sigmabar)) +
                               m_B2 * (84 * m_v4 + 6.0 * m_v2 * q2 * (28 - 57 * sigmabar) +
                               q2 * sigmabar * (15.0 * omega_1 * (2.0 * omega_1 + omega_2) - 2.0 * q2 * (63 + 2.0 * sigmabar))) +
                               42 * omega_1 * (2.0 * omega_1 + omega_2) * power_of<2>(q2)) -
                               m_B * sigma4 * (-30 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar2 +
                               15.0 * m_B2 * (4.0 * omega_1 + omega_2) * q2 * sigmabar + 60 * m_B5 * sigmabar2 * (-3.0 + 2.0 * sigmabar) +
                               12.0 * m_B * q2 * (14.0 * m_v2 + q2 * (7.0 - 10.0 * sigmabar)) +
                               10.0 * m_B3 * sigmabar * (3.0 * m_v2 + q2 * (6.0 + 5.0 * sigmabar)) + 42 * (4.0 * omega_1 + omega_2) * power_of<2>(q2)
                               ) + omega_1 * (2.0 * omega_1 + omega_2) *
                               (m_v4 * (-42 + 87 * sigmabar) + 3.0 * m_v2 * q2 * sigmabar * (13.0 + 7.0 * sigmabar) +
                               5.0 * sigmabar2 * (5.0 - 2.0 * sigmabar) * power_of<2>(q2)) +
                               m_B * (4.0 * omega_1 + omega_2) * sigmabar *
                               (-(3.0 * m_v2 * q2 * sigmabar * (3.0 + 5.0 * sigmabar)) + 3.0 * m_v4 * (4.0 - 9.0 * sigmabar) +
                               sigmabar2 * (-11.0 + 2.0 * sigmabar) * power_of<2>(q2)) +
                               sigma2 * (6.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * (-(7.0 * (2.0 * m_v2 + q2)) + 6.0 * q2 * sigmabar) +
                               9.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar2 * (11.0 + 2.0 * sigmabar * (-6.0 + sigmabar)) +
                               6.0 * m_B6 * sigmabar2 * (13.0 + 3.0 * sigmabar * (-13.0 + 6.0 * sigmabar)) -
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar *
                               (5.0 * m_v2 * (6.0 + 7.0 * sigmabar) + q2 * (15.0 + (91 - 54 * sigmabar) * sigmabar)) -
                               2.0 * m_B4 * sigmabar * (9.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar * (-5.0 + 2.0 * sigmabar) +
                               m_v2 * (15.0 - 54 * sigmabar2 + 20.0 * sigmabar) + q2 * sigmabar * (56 + sigmabar * (-131 + 30 * sigmabar))) +
                               m_B2 * (6.0 * m_v4 * (-14.0 + 37 * sigmabar) -
                               m_v2 * sigmabar * (15.0 * omega_1 * (2.0 * omega_1 + omega_2) + 2.0 * q2 * (-87 + 49 * sigmabar)) -
                               q2 * sigmabar * (2.0 * q2 * sigmabar * (-5.0 + 26 * sigmabar) +
                               5.0 * omega_1 * (2.0 * omega_1 + omega_2) * (6.0 + 7.0 * sigmabar))) -
                               m_B * (4.0 * omega_1 + omega_2) * (42 * m_v4 + 21 * m_v2 * q2 * (4.0 - 7.0 * sigmabar) -
                               sigmabar * (51 + 10.0 * sigmabar) * power_of<2>(q2))) +
                               sigma * (6.0 * m_B6 * sigmabar3 * (8.0 - 9.0 * sigmabar) -
                               3.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar2 * (13.0 + 9.0 * sigmabar * (-3.0 + sigmabar)) +
                               m_B4 * sigmabar2 * (m_v2 * (-20.0 + 30 * sigmabar2 - 68 * sigmabar) +
                               6.0 * q2 * sigmabar * (-15.0 + 11.0 * sigmabar) -
                               9.0 * omega_1 * (2.0 * omega_1 + omega_2) * (11.0 + sigmabar * (-8.0 + sigmabar))) +
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar *
                               (m_v2 * (15.0 + (30 - 37 * sigmabar) * sigmabar) + q2 * sigmabar * (61 + 3.0 * sigmabar * (-31 + 5.0 * sigmabar))) +
                               m_B2 * sigmabar * (16.0 * m_v4 * (-3.0 + 7.0 * sigmabar) +
                               m_v2 * (2.0 * q2 * sigmabar * (22 + 25 * sigmabar) +
                               10.0 * omega_1 * (2.0 * omega_1 + omega_2) * (3.0 + 4.0 * sigmabar)) +
                               q2 * (6.0 * sigmabar2 * q2 * (7.0 - 2.0 * sigmabar) -
                               omega_1 * (2.0 * omega_1 + omega_2) * (-15.0 + sigmabar * (-101 + 37 * sigmabar)))) +
                               omega_1 * (2.0 * omega_1 + omega_2) * (42 * m_v4 + 3.0 * m_v2 * q2 * (28 - 41 * sigmabar) -
                               sigmabar * (39 + 16.0 * sigmabar) * power_of<2>(q2)) -
                               m_B * (4.0 * omega_1 + omega_2) * (m_v4 * (-42 + 99 * sigmabar) + 3.0 * m_v2 * q2 * sigmabar * (21 - 4.0 * sigmabar) +
                               2.0 * sigmabar2 * (8.0 - 9.0 * sigmabar) * power_of<2>(q2))))
                             / (m_B2 * power_of<2>(omega_2) * power_of<8>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d2B_fT_3pt_phi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double phi_bar_4     = this->phi_bar_4(omega_1, omega_2);
            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            return  pow(m_B,-1) * pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-7.0) *
                    (4.0 * phi_bar_bar_4 * (-(m_B * sigmabar * (-m_v2 + (m_B2 - q2) * sigmabar) *
                    (3.0 * (m_B2 - q2) * sigmabar + m_v2 * (2.0 - 5.0 * sigmabar))) +
                    m_B * sigma * (m_v4 * (12.0 - 29 * sigmabar) +
                    2.0 * sigmabar2 * (m_B2 - q2) * (q2 * (5.0 - 2.0 * sigmabar) + m_B2 * (-6.0 + 9.0 * sigmabar)) +
                    m_v2 * sigmabar * (-(q2 * (13.0 + 10.0 * sigmabar)) + m_B2 * (5.0 + 2.0 * sigmabar * (8.0 - 5.0 * sigmabar)))) +
                    omega_1 * (3.0 * m_v4 * (-4.0 + 9.0 * sigmabar) -
                    sigmabar2 * (m_B2 - q2) * (q2 * (11.0 - 2.0 * sigmabar) + 3.0 * m_B2 * (-4.0 + 3.0 * sigmabar)) +
                    m_v2 * sigmabar * (3.0 * q2 * (3.0 + 5.0 * sigmabar) + m_B2 * (-5.0 + sigmabar * (-18.0 + 5.0 * sigmabar)))) +
                    sigma4 * (12.0 * m_B5 * sigmabar2 - 5.0 * m_B3 * q2 * sigmabar - 12.0 * m_B * power_of<2>(q2)) +
                    sigma3 * (-(12.0 * m_B4 * omega_1 * sigmabar2) + 5.0 * m_B2 * omega_1 * q2 * sigmabar +
                    12.0 * m_B5 * sigmabar2 * (-3.0 + sigmabar) + 4.0 * m_B * q2 * (6.0 * m_v2 + 3.0 * q2 - 2.0 * q2 * sigmabar) +
                    m_B3 * sigmabar * (5.0 * m_v2 + 2.0 * q2 * (5.0 + 8.0 * sigmabar)) + 12.0 * omega_1 * power_of<2>(q2)) +
                    sigma2 * (-(9.0 * m_B4 * omega_1 * sigmabar2 * (-4.0 + sigmabar)) +
                    9.0 * m_B5 * sigmabar2 * (4.0 - 3.0 * sigmabar) +
                    6.0 * omega_1 * q2 * (-(4.0 * m_v2) + q2 * (-2.0 + sigmabar)) -
                    m_B2 * omega_1 * sigmabar * (5.0 * m_v2 + q2 * (10.0 + 17.0 * sigmabar)) +
                    m_B3 * sigmabar * (-(m_v2 * (10.0 + 17.0 * sigmabar)) + q2 * (-5.0 + sigmabar * (-38 + 15.0 * sigmabar))) +
                    m_B * (-(12.0 * m_v4) + m_v2 * q2 * (-24 + 37 * sigmabar) + sigmabar * (11.0 + 8.0 * sigmabar) * power_of<2>(q2))) +
                    omega_1 * sigma * (12.0 * m_v4 + m_v2 * (q2 * (24 - 33 * sigmabar) + 2.0 * m_B2 * sigmabar * (5.0 + 9.0 * sigmabar)) +
                    sigmabar * (18.0 * m_B4 * sigmabar * (-2.0 + sigmabar) +
                    5.0 * m_B2 * q2 * (1 - 2.0 * sigmabar * (-4.0 + sigmabar)) - (9.0 + 8.0 * sigmabar) * power_of<2>(q2)))) +
                    2.0 * m_B * (-omega_1 + m_B * sigma) * sigmabar *
                    (3.0 * m_B2 * sigmabar3 - 2.0 * q2 * sigma * sigmabar + q2 * (-3.0 + 2.0 * sigma) * sigmabar +
                    m_v2 * (2.0 - 2.0 * sigma - 5.0 * sigmabar)) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * phi_bar_4);
        }

        double I4d2C_fT_3pt_phi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -(2.0 * (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) +
                                (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) * sigmabar) *
                                (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)))
                             /  (power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * phi_bar_bar_4;
        }

        double I4d2D_fT_3pt_phi_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_T proportional to phi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double phi_bar_4     = this->phi_bar_4(omega_1, omega_2);
            const double phi_bar_bar_4 = this->phi_bar_bar_4(omega_1, omega_2);

            return   2.0 * pow(m_B,-2.0) * pow(sigmabar,-7.0) * (phi_bar_bar_4 *
                     ((m_v2 - q2 * sigma) * (-(7.0 * m_v2) + m_B2 * sigmabar2 - 5.0 * q2 + 6.0 * q2 * sigma) * sigmabar +
                     sigmabar2 * (q2 * (5.0 * m_v2 + q2) - m_B2 * (2.0 * m_v2 + q2 - 3.0 * q2 * sigma) * sigmabar) +
                     sigmabar3 * (9.0 * m_B4 * sigmabar2 + m_B2 * (-(5.0 * m_v2) - 11.0 * q2 + 10.0 * q2 * sigma) +
                     2.0 * power_of<2>(q2)) + 4.0 * sigmabar * power_of<2>(m_v2 - q2 * sigma)) +
                     m_B * sigmabar * (3.0 * m_B2 * sigmabar3 - 2.0 * q2 * sigma * sigmabar + q2 * (-3.0 + 2.0 * sigma) * sigmabar +
                     m_v2 * (2.0 - 2.0 * sigma - 5.0 * sigmabar)) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * phi_bar_4);
        }

        double I2_fT_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to psi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double sigmabar = 1.0 - sigma;
            const double m_v      = this->m_v();

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_2 = -4.0 * m_v / (m_B2 * power_of<3>(sigmabar));

            return C_2 * psi_bar_4;
        }

        double I3_fT_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = 4.0 * m_v * (-m_v2 + q2 - 2.0 * q2 * sigmabar + m_B2 * sigmabar2)
                             / (m_B2 * power_of<4>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I3d1A_fT_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = -4.0 * m_v * (4.0 * m_v2 + 2.0 * q2 * sigmabar + m_B2 * (-3.0 + sigmabar) * sigmabar + sigma * (-4.0 * q2 + 3.0 * m_B2 * sigmabar))
                             /  (m_B2 * power_of<5>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I3d1B_fT_3pt_psi_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = 4.0 * m_v * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B * (-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I3d1C_fT_3pt_psi_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

            const double C_3 = -4.0 * m_v * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B * omega_2 * power_of<4>(sigmabar));

            return C_3 * psi_bar_4;
        }

        double I3_fT_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 12.0 * m_v * u / (m_B * power_of<3>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1A_fT_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 12.0 * m_v * (3.0 * sigma * m_B + m_B * sigmabar - 3.0 * omega_1)
                               / (m_B * omega_2 * power_of<4>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1B_fT_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = -12.0 * m_v / ((-omega_1 + m_B * sigma) * power_of<3>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1C_fT_3pt_psiA_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4

            return 0.0;
        }

        double I4_fT_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -12.0 * m_v * u * (-m_v2 + q2 - 2.0 * q2 * sigmabar + m_B2 * sigmabar2)
                             / (m_B * power_of<4>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1A_fT_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 12.0 * m_v * (-(omega_1 * (4.0 * m_v2 + 2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar) +
                               sigma * (-(4.0 * q2) + 3.0 * m_B2 * sigmabar))) +
                               m_B * (4.0 * sigma * (m_v2 - q2 * sigma) + sigmabar2 * (q2 + m_B2 * (-1 + 2.0 * sigma)) +
                               sigmabar * (m_v2 + sigma * (q2 - 3.0 * m_B2 * sigmabar))))
                             / (m_B * omega_2 * power_of<5>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1B_fT_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -12.0 * m_v * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar))
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * psi_bar_bar_4;
        }

        double I4d1C_fT_3pt_psiA_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4

            return 0.0;
        }
        double I4d2A_fT_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_B2 = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 24.0 * m_v * (2.0 * sigma2 * (-(5.0 * m_B * q2) + 3.0 * m_B3 * sigmabar) -
                               omega_1 * (10.0 * m_v2 + 2.0 * q2 * sigmabar + 3.0 * m_B2 * sigmabar * (-2.0 + sigmabar)) +
                               m_B * sigmabar * (4.0 * m_v2 + 2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar)) +
                               2.0 * sigma * (5.0 * omega_1 * q2 - 3.0 * m_B2 * omega_1 * sigmabar + 3.0 * m_B3 * (-1 + sigmabar) * sigmabar +
                               m_B * (5.0 * m_v2 - q2 * sigmabar)))
                             / (m_B * omega_2 * power_of<6>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2B_fT_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_B2 = power_of<2>(m_B);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return  pow(-omega_1 + m_B * sigma,-1) * pow(sigmabar,-5.0) *
                    (-24.0 * m_v * (4.0 * m_v2 + 2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar) +
                    sigma * (-(4.0 * q2) + 3.0 * m_B2 * sigmabar)) * psi_bar_bar_4 -
                    12.0 * m_B * m_v * sigmabar * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * psi_bar_4);
        }

        double I4d2C_fT_3pt_psiA_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2     = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 12.0 * m_B * m_v * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar)) / (power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2D_fT_3pt_psiA_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_B2 = power_of<2>(m_B);
            const double sigmabar = 1.0 - sigma;

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return     12.0 * m_v * pow(m_B,-1) * pow(sigmabar,-4.0) *
                       ((2.0 * q2 - 2.0 * m_B2 * sigmabar) * psi_bar_bar_4 -
                        m_B * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * psi_bar_4);
        }

        double I2_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_B2 = power_of<2>(m_B);
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_2 = -4.0 * u * (-1.0 + 2.0 * u) / (m_B2 * power_of<3>(sigmabar));

            return C_2 * psi_bar_bar_4;
        }

        double I3_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * u * (-1.0 + 2.0 * u) * (m_v2 + q2 * (5.0 - 4.0 * sigmabar) - m_B2 * sigmabar2)
                             / (m_B2 * power_of<4>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1A_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B),   m_B3 = power_of<3>(m_B),   m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(),   m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma), sigma4   = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (m_B3 * (4.0 * omega_1 + omega_2) * sigmabar2 +
                               2.0 * m_B2 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar * (-2.0 + sigmabar) -
                               8.0 * sigma4 * (5.0 * m_B2 * q2 + m_B4 * sigmabar) +
                               4.0 * m_B * sigma3 * (5.0 * (4.0 * omega_1 + omega_2) * q2 + 10.0 * m_B * (m_v2 + q2) +
                               m_B2 * (4.0 * omega_1 + omega_2) * sigmabar - 8.0 * m_B * q2 * sigmabar - 2.0 * m_B3 * sigmabar * (-2.0 + sigmabar)) -
                               m_B * (4.0 * omega_1 + omega_2) * sigmabar * (q2 * sigmabar + m_v2 * (-4.0 + 5.0 * sigmabar)) +
                               4.0 * omega_1 * (2.0 * omega_1 + omega_2) * (2.0 * q2 * sigmabar + m_v2 * (-5.0 + 6.0 * sigmabar)) -
                               2.0 * sigma * (2.0 * m_B4 * sigmabar2 + 2.0 * m_B3 * (4.0 * omega_1 + omega_2) * (-1 + sigmabar) * sigmabar -
                               2.0 * omega_1 * (2.0 * omega_1 + omega_2) * (5.0 * (m_v2 + q2) - 2.0 * q2 * sigmabar) +
                               m_B2 * sigmabar * (-(2.0 * q2 * sigmabar) + omega_1 * (2.0 * omega_1 + omega_2) * (-4.0 + sigmabar) +
                               m_v2 * (8.0 - 10.0 * sigmabar)) + 2.0 * m_B * (4.0 * omega_1 + omega_2) *
                              (3.0 * q2 * sigmabar + m_v2 * (-5.0 + 7.0 * sigmabar))) +
                               sigma2 * (-(20.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2) + 4.0 * m_B4 * sigmabar * (-2.0 + 3.0 * sigmabar) +
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar * (-8.0 + 3.0 * sigmabar) -
                               4.0 * m_B * (4.0 * omega_1 + omega_2) * (5.0 * (m_v2 + q2) - 3.0 * q2 * sigmabar) +
                               4.0 * m_B2 * (-(omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) + 8.0 * q2 * sigmabar +
                               2.0 * m_v2 * (-5.0 + 8.0 * sigmabar))))
                             / (m_B2 * power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1B_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(),   m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (m_v2 * (4.0 - 4.0 * sigma - 5.0 * sigmabar) - 4.0 * sigma * q2 * sigmabar + sigmabar * (-q2 + m_B2 * sigmabar2))
                             / (m_B * (-omega_1 + m_B * sigma) * power_of<5>(sigmabar));

            return C_3 * psi_bar_bar_4;
        }

        double I3d1C_fT_3pt_psiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4

            return 0.0;
        }

        double I4_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B4     = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2     = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * u * (-1.0 + 2.0 * u) * (m_v4 - 2.0 * m_B2 * q2 * sigmabar3 + m_B4 * sigmabar4 + q2 * q2 * (-1.0 + 2.0 * sigmabar)
                             + 2.0 * m_v2 * sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B2 * power_of<5>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1A_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -(2.0 * (6.0 * (-omega_1 + m_B * sigma) * (2.0 * omega_1 + omega_2 - 2.0 * m_B * sigma) *
                                (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) + (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) *
                                sigmabar) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) -
                                2.0 * m_B * (-omega_1 + m_B * sigma) * sigmabar *
                                (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) + (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) *
                                sigmabar) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) +
                                m_B * (2.0 * omega_1 + omega_2 - 2.0 * m_B * sigma) * sigmabar *
                                (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) + (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) *
                                sigmabar) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) +
                                (-omega_1 + m_B * sigma) * (-(2.0 * omega_1) - omega_2 + 2.0 * m_B * sigma) * sigmabar *
                                (3.0 * m_v4 + 12.0 * m_v2 * q2 * sigmabar - 3.0 * m_B4 * sigmabar3 * (-2.0 + 2.0 * sigma - 3.0 * sigmabar) +
                                m_B2 * (sigmabar2 * (-m_v2 + q2 * sigma) + sigmabar2 * (-(5.0 * m_v2) + q2 * (-11.0 + 10.0 * sigma)) -
                                sigmabar2 * (12.0 * m_v2 + 13.0 * q2 * sigmabar)) +
                                (-(6.0 * sigma2) + sigmabar * (7.0 + 2.0 * sigmabar) + sigma * (3.0 - 4.0 * sigmabar)) * power_of<2>(q2))))
                             / (m_B2 * power_of<2>(omega_2) * power_of<7>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d1B_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 2.0 * (3.0 * m_B2 * sigmabar3 - 2.0 * q2 * sigma * sigmabar + q2 * (-3.0 + 2.0 * sigma) * sigmabar +
                               m_v2 * (2.0 - 2.0 * sigma - 5.0 * sigmabar)) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B * power_of<6>(sigmabar) * omega_2);

            return C_4 * psi_bar_bar_4;
        }

        double I4d1C_fT_3pt_psiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4

            return 0.0;
        }

        double I4d2A_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B), m_B6 = power_of<6>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma), sigma5 = power_of<5>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 4.0 * (6.0 * m_B6 * sigmabar4 + 3.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar3 * (-4.0 + 3.0 * sigmabar) +
                               m_B4 * sigmabar2 * (-(12.0 * sigmabar2 * q2) - 2.0 * m_v2 * sigmabar * (1 + 5.0 * sigmabar) +
                               3.0 * omega_1 * (2.0 * omega_1 + omega_2) * (13.0 + 3.0 * sigmabar * (-5.0 + sigmabar))) -
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar2 *
                               (q2 * sigmabar * (-23 + 11.0 * sigmabar) + m_v2 * (-5.0 + sigmabar * (-18.0 + 5.0 * sigmabar))) +
                               m_B2 * sigmabar * (2.0 * m_v4 * sigmabar * (-2.0 + 5.0 * sigmabar) +
                               m_v2 * (2.0 * sigmabar2 * q2 * (1 + 5.0 * sigmabar) +
                               omega_1 * (2.0 * omega_1 + omega_2) * (-15.0 + sigmabar * (-40 + 19.0 * sigmabar))) +
                               q2 * sigmabar * (6.0 * sigmabar2 * q2 - omega_1 * (2.0 * omega_1 + omega_2) * (66 + sigmabar * (-53 + 5.0 * sigmabar)))) +
                               sigma5 * (-60 * m_B6 * sigmabar2 + 30 * m_B4 * q2 * sigmabar + 84 * m_B2 * power_of<2>(q2)) +
                               sigma3 * (6.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar2 * (-15.0 + 8.0 * sigmabar) +
                               15.0 * m_B3 * (4.0 * omega_1 + omega_2) * sigmabar * (m_v2 + 2.0 * q2 * (1 + sigmabar)) -
                               6.0 * m_B6 * sigmabar2 * (33 + 2.0 * sigmabar * (-24 + 5.0 * sigmabar)) +
                               6.0 * m_B * (4.0 * omega_1 + omega_2) * q2 * (14.0 * m_v2 + q2 * (7.0 - 8.0 * sigmabar)) +
                               2.0 * m_B4 * sigmabar * (-(15.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) + 30 * m_v2 * (1 + sigmabar) +
                               q2 * (15.0 + (81 - 70 * sigmabar) * sigmabar)) +
                               m_B2 * (84 * m_v4 + 6.0 * m_v2 * q2 * (28 - 57 * sigmabar) +
                               q2 * sigmabar * (15.0 * omega_1 * (2.0 * omega_1 + omega_2) - 2.0 * q2 * (63 + 2.0 * sigmabar))) +
                               42 * omega_1 * (2.0 * omega_1 + omega_2) * power_of<2>(q2)) -
                               m_B * sigma4 * (-30 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar2 +
                               15.0 * m_B2 * (4.0 * omega_1 + omega_2) * q2 * sigmabar + 60 * m_B5 * sigmabar2 * (-3.0 + 2.0 * sigmabar) +
                               12.0 * m_B * q2 * (14.0 * m_v2 + q2 * (7.0 - 10.0 * sigmabar)) +
                               10.0 * m_B3 * sigmabar * (3.0 * m_v2 + q2 * (6.0 + 5.0 * sigmabar)) + 42 * (4.0 * omega_1 + omega_2) * power_of<2>(q2)
                               ) + omega_1 * (2.0 * omega_1 + omega_2) *
                               (m_v4 * (-42 + 87 * sigmabar) + 3.0 * m_v2 * q2 * sigmabar * (13.0 + 7.0 * sigmabar) +
                               5.0 * sigmabar2 * (5.0 - 2.0 * sigmabar) * power_of<2>(q2)) +
                               m_B * (4.0 * omega_1 + omega_2) * sigmabar *
                               (-(3.0 * m_v2 * q2 * sigmabar * (3.0 + 5.0 * sigmabar)) + 3.0 * m_v4 * (4.0 - 9.0 * sigmabar) +
                               sigmabar2 * (-11.0 + 2.0 * sigmabar) * power_of<2>(q2)) +
                               sigma2 * (6.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * (-(7.0 * (2.0 * m_v2 + q2)) + 6.0 * q2 * sigmabar) +
                               9.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar2 * (11.0 + 2.0 * sigmabar * (-6.0 + sigmabar)) +
                               6.0 * m_B6 * sigmabar2 * (13.0 + 3.0 * sigmabar * (-13.0 + 6.0 * sigmabar)) -
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar *
                               (5.0 * m_v2 * (6.0 + 7.0 * sigmabar) + q2 * (15.0 + (91 - 54 * sigmabar) * sigmabar)) -
                               2.0 * m_B4 * sigmabar * (9.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar * (-5.0 + 2.0 * sigmabar) +
                               m_v2 * (15.0 - 54 * sigmabar2 + 20.0 * sigmabar) + q2 * sigmabar * (56 + sigmabar * (-131 + 30 * sigmabar))) +
                               m_B2 * (6.0 * m_v4 * (-14.0 + 37 * sigmabar) -
                               m_v2 * sigmabar * (15.0 * omega_1 * (2.0 * omega_1 + omega_2) + 2.0 * q2 * (-87 + 49 * sigmabar)) -
                               q2 * sigmabar * (2.0 * q2 * sigmabar * (-5.0 + 26 * sigmabar) +
                               5.0 * omega_1 * (2.0 * omega_1 + omega_2) * (6.0 + 7.0 * sigmabar))) -
                               m_B * (4.0 * omega_1 + omega_2) * (42 * m_v4 + 21 * m_v2 * q2 * (4.0 - 7.0 * sigmabar) -
                               sigmabar * (51 + 10.0 * sigmabar) * power_of<2>(q2))) +
                               sigma * (6.0 * m_B6 * sigmabar3 * (8.0 - 9.0 * sigmabar) -
                               3.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar2 * (13.0 + 9.0 * sigmabar * (-3.0 + sigmabar)) +
                               m_B4 * sigmabar2 * (m_v2 * (-20.0 + 30 * sigmabar2 - 68 * sigmabar) +
                               6.0 * q2 * sigmabar * (-15.0 + 11.0 * sigmabar) -
                               9.0 * omega_1 * (2.0 * omega_1 + omega_2) * (11.0 + sigmabar * (-8.0 + sigmabar))) +
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar *
                               (m_v2 * (15.0 + (30 - 37 * sigmabar) * sigmabar) + q2 * sigmabar * (61 + 3.0 * sigmabar * (-31 + 5.0 * sigmabar))) +
                               m_B2 * sigmabar * (16.0 * m_v4 * (-3.0 + 7.0 * sigmabar) +
                               m_v2 * (2.0 * q2 * sigmabar * (22 + 25 * sigmabar) +
                               10.0 * omega_1 * (2.0 * omega_1 + omega_2) * (3.0 + 4.0 * sigmabar)) +
                               q2 * (6.0 * sigmabar2 * q2 * (7.0 - 2.0 * sigmabar) -
                               omega_1 * (2.0 * omega_1 + omega_2) * (-15.0 + sigmabar * (-101 + 37 * sigmabar)))) +
                               omega_1 * (2.0 * omega_1 + omega_2) * (42 * m_v4 + 3.0 * m_v2 * q2 * (28 - 41 * sigmabar) -
                               sigmabar * (39 + 16.0 * sigmabar) * power_of<2>(q2)) -
                               m_B * (4.0 * omega_1 + omega_2) * (m_v4 * (-42 + 99 * sigmabar) + 3.0 * m_v2 * q2 * sigmabar * (21 - 4.0 * sigmabar) +
                               2.0 * sigmabar2 * (8.0 - 9.0 * sigmabar) * power_of<2>(q2))))
                             / (m_B2 * power_of<2>(omega_2) * power_of<8>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2B_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return  pow(m_B,-1) * pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-7.0) *
                    (4.0 * psi_bar_bar_4 * (-(m_B * sigmabar * (-m_v2 + (m_B2 - q2) * sigmabar) *
                    (3.0 * (m_B2 - q2) * sigmabar + m_v2 * (2.0 - 5.0 * sigmabar))) +
                    m_B * sigma * (m_v4 * (12.0 - 29 * sigmabar) +
                    2.0 * sigmabar2 * (m_B2 - q2) * (q2 * (5.0 - 2.0 * sigmabar) + m_B2 * (-6.0 + 9.0 * sigmabar)) +
                    m_v2 * sigmabar * (-(q2 * (13.0 + 10.0 * sigmabar)) + m_B2 * (5.0 + 2.0 * sigmabar * (8.0 - 5.0 * sigmabar)))) +
                    omega_1 * (3.0 * m_v4 * (-4.0 + 9.0 * sigmabar) -
                    sigmabar2 * (m_B2 - q2) * (q2 * (11.0 - 2.0 * sigmabar) + 3.0 * m_B2 * (-4.0 + 3.0 * sigmabar)) +
                    m_v2 * sigmabar * (3.0 * q2 * (3.0 + 5.0 * sigmabar) + m_B2 * (-5.0 + sigmabar * (-18.0 + 5.0 * sigmabar)))) +
                    sigma4 * (12.0 * m_B5 * sigmabar2 - 5.0 * m_B3 * q2 * sigmabar - 12.0 * m_B * power_of<2>(q2)) +
                    sigma3 * (-(12.0 * m_B4 * omega_1 * sigmabar2) + 5.0 * m_B2 * omega_1 * q2 * sigmabar +
                    12.0 * m_B5 * sigmabar2 * (-3.0 + sigmabar) + 4.0 * m_B * q2 * (6.0 * m_v2 + 3.0 * q2 - 2.0 * q2 * sigmabar) +
                    m_B3 * sigmabar * (5.0 * m_v2 + 2.0 * q2 * (5.0 + 8.0 * sigmabar)) + 12.0 * omega_1 * power_of<2>(q2)) +
                    sigma2 * (-(9.0 * m_B4 * omega_1 * sigmabar2 * (-4.0 + sigmabar)) +
                    9.0 * m_B5 * sigmabar2 * (4.0 - 3.0 * sigmabar) +
                    6.0 * omega_1 * q2 * (-(4.0 * m_v2) + q2 * (-2.0 + sigmabar)) -
                    m_B2 * omega_1 * sigmabar * (5.0 * m_v2 + q2 * (10.0 + 17.0 * sigmabar)) +
                    m_B3 * sigmabar * (-(m_v2 * (10.0 + 17.0 * sigmabar)) + q2 * (-5.0 + sigmabar * (-38 + 15.0 * sigmabar))) +
                    m_B * (-(12.0 * m_v4) + m_v2 * q2 * (-24 + 37 * sigmabar) + sigmabar * (11.0 + 8.0 * sigmabar) * power_of<2>(q2))) +
                    omega_1 * sigma * (12.0 * m_v4 + m_v2 * (q2 * (24 - 33 * sigmabar) + 2.0 * m_B2 * sigmabar * (5.0 + 9.0 * sigmabar)) +
                    sigmabar * (18.0 * m_B4 * sigmabar * (-2.0 + sigmabar) +
                    5.0 * m_B2 * q2 * (1 - 2.0 * sigmabar * (-4.0 + sigmabar)) - (9.0 + 8.0 * sigmabar) * power_of<2>(q2)))) +
                    2.0 * m_B * (-omega_1 + m_B * sigma) * sigmabar *
                    (3.0 * m_B2 * sigmabar3 - 2.0 * q2 * sigma * sigmabar + q2 * (-3.0 + 2.0 * sigma) * sigmabar +
                    m_v2 * (2.0 - 2.0 * sigma - 5.0 * sigmabar)) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * psi_bar_4);
        }

        double I4d2C_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -(2.0 * (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) +
                                (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) * sigmabar) *
                                (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)))
                             /  (power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * psi_bar_bar_4;
        }

        double I4d2D_fT_3pt_psiB_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_T proportional to psi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double psi_bar_4     = this->psi_bar_4(omega_1, omega_2);
            const double psi_bar_bar_4 = this->psi_bar_bar_4(omega_1, omega_2);

            return   2.0 * pow(m_B,-2.0) * pow(sigmabar,-7.0) * (psi_bar_bar_4 *
                     ((m_v2 - q2 * sigma) * (-(7.0 * m_v2) + m_B2 * sigmabar2 - 5.0 * q2 + 6.0 * q2 * sigma) * sigmabar +
                     sigmabar2 * (q2 * (5.0 * m_v2 + q2) - m_B2 * (2.0 * m_v2 + q2 - 3.0 * q2 * sigma) * sigmabar) +
                     sigmabar3 * (9.0 * m_B4 * sigmabar2 + m_B2 * (-(5.0 * m_v2) - 11.0 * q2 + 10.0 * q2 * sigma) +
                     2.0 * power_of<2>(q2)) + 4.0 * sigmabar * power_of<2>(m_v2 - q2 * sigma)) +
                     m_B * sigmabar * (3.0 * m_B2 * sigmabar3 - 2.0 * q2 * sigma * sigmabar + q2 * (-3.0 + 2.0 * sigma) * sigmabar +
                     m_v2 * (2.0 - 2.0 * sigma - 5.0 * sigmabar)) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * psi_bar_4);
        }

        double I2_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - 0.0                                                        - I2_fT_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I3_fT_3pt_psiA_bar_bar_4(sigma, omega_1, omega_2, q2)      - I3_fT_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1A_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I3d1A_fT_3pt_psiA_bar_bar_4(sigma, omega_1, omega_2, q2)  - I3d1A_fT_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1B_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return - I3d1B_fT_3pt_psiA_bar_bar_4(sigma, omega_1, q2)           - I3d1B_fT_3pt_psiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I3d1C_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return - I3d1C_fT_3pt_psiA_bar_bar_4(sigma, omega_2, q2)           - I3d1C_fT_3pt_psiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I4_fT_3pt_psiA_bar_bar_4( sigma, omega_1, omega_2, q2)    - I4_fT_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1A_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I4d1A_fT_3pt_psiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d1A_fT_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1B_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return - I4d1B_fT_3pt_psiA_bar_bar_4( sigma, omega_1, q2)          - I4d1B_fT_3pt_psiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d1C_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return - I4d1C_fT_3pt_psiA_bar_bar_4( sigma, omega_2, q2)          - I4d1C_fT_3pt_psiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2A_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return - I4d2A_fT_3pt_psiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d2A_fT_3pt_psiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d2B_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return - I4d2B_fT_3pt_psiA_bar_bar_4( sigma, omega_1, q2)          - I4d2B_fT_3pt_psiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d2C_fT_3pt_psi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return - I4d2C_fT_3pt_psiA_bar_bar_4( sigma, omega_2, q2)          - I4d2C_fT_3pt_psiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2D_fT_3pt_psi_bar_bar_4(const double & sigma, const double & q2) const
        {
            return - I4d2D_fT_3pt_psiA_bar_bar_4( sigma, q2)                   - I4d2D_fT_3pt_psiB_bar_bar_4(sigma, q2);
        }

        double I2_fT_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to chi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_2 = 4.0 * (m_v + m_B * u * sigmabar) / (m_B2 * power_of<3>(sigmabar));

            return C_2 * chi_bar_4;
        }

        double I3_fT_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = -4.0 * (m_v + m_B * u * sigmabar) * (-m_v2 + q2 - 2.0 * q2 * sigmabar + m_B2 * sigmabar2)
                             / (m_B2 * power_of<4>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I3d1A_fT_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v3 = power_of<3>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = 4.0 * (4.0 * m_v3 * omega_2 + m_B * m_v2 * sigmabar * (-(3.0 * omega_1) + 3.0 * m_B * sigma + m_B * sigmabar) +
                               m_v * omega_2 * (2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar) +
                               sigma * (-(4.0 * q2) + 3.0 * m_B2 * sigmabar)) +
                               m_B * sigmabar * (-(3.0 * m_B * sigma2 * q2) - 2.0 * m_B3 * sigmabar2 * sigma +
                               m_B * sigmabar2 * (q2 + m_B2 * (-1 + 2.0 * sigma)) +
                               omega_1 * (sigma * (3.0 * q2 - 2.0 * m_B2 * sigmabar) - sigmabar * (q2 + m_B2 * (-2.0 + sigmabar)))))
                             / (m_B2 * omega_2 * power_of<5>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I3d1B_fT_3pt_chi_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = -4.0 * (m_v + m_B * sigmabar) * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B * (-omega_1 + m_B * sigma) * power_of<4>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I3d1C_fT_3pt_chi_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

            const double C_3 = 4.0 * m_v * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B * omega_2 * power_of<4>(sigmabar));

            return C_3 * chi_bar_4;
        }

        double I3_fT_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 12.0 * m_v * u / (m_B * power_of<3>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1A_fT_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 12.0 * m_v * (3.0 * sigma * m_B + m_B * sigmabar - 3.0 * omega_1)
                               / (m_B * omega_2 * power_of<4>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1B_fT_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v();
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = -12.0 * m_v / ((-omega_1 + m_B * sigma) * power_of<3>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1C_fT_3pt_chiA_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4

            return 0.0;
        }

        double I4_fT_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -12.0 * m_v * u * (-m_v2 + q2 - 2.0 * q2 * sigmabar + m_B2 * sigmabar2)
                             / (m_B * power_of<4>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1A_fT_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 12.0 * m_v * (-(omega_1 * (4.0 * m_v2 + 2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar) +
                               sigma * (-(4.0 * q2) + 3.0 * m_B2 * sigmabar))) +
                               m_B * (4.0 * sigma * (m_v2 - q2 * sigma) + sigmabar2 * (q2 + m_B2 * (-1 + 2.0 * sigma)) +
                               sigmabar * (m_v2 + sigma * (q2 - 3.0 * m_B2 * sigmabar))))
                             / (m_B * omega_2 * power_of<5>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1B_fT_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -12.0 * m_v * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar))
                             / (power_of<4>(sigmabar) * omega_2);

            return C_4 * chi_bar_bar_4;
        }

        double I4d1C_fT_3pt_chiA_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4

            return 0.0;
        }
        double I4d2A_fT_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_B2 = power_of<2>(m_B), m_B3 = power_of<3>(m_B);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 24.0 * m_v * (2.0 * sigma2 * (-(5.0 * m_B * q2) + 3.0 * m_B3 * sigmabar) -
                               omega_1 * (10.0 * m_v2 + 2.0 * q2 * sigmabar + 3.0 * m_B2 * sigmabar * (-2.0 + sigmabar)) +
                               m_B * sigmabar * (4.0 * m_v2 + 2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar)) +
                               2.0 * sigma * (5.0 * omega_1 * q2 - 3.0 * m_B2 * omega_1 * sigmabar + 3.0 * m_B3 * (-1 + sigmabar) * sigmabar +
                               m_B * (5.0 * m_v2 - q2 * sigmabar)))
                             / (m_B * omega_2 * power_of<6>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2B_fT_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_B2 = power_of<2>(m_B);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return  pow(-omega_1 + m_B * sigma,-1) * pow(sigmabar,-5.0) *
                    (-24.0 * m_v * (4.0 * m_v2 + 2.0 * q2 * sigmabar + m_B2 * sigmabar * (-3.0 + sigmabar) +
                    sigma * (-(4.0 * q2) + 3.0 * m_B2 * sigmabar)) * chi_bar_bar_4 -
                    12.0 * m_B * m_v * sigmabar * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * chi_bar_4);
        }

        double I4d2C_fT_3pt_chiA_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2     = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 12.0 * m_B * m_v * (m_v2 - sigma * q2 + sigmabar * (q2 - m_B2 * sigmabar)) / (power_of<2>(omega_2) * power_of<4>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2D_fT_3pt_chiA_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_B2 = power_of<2>(m_B);
            const double sigmabar = 1.0 - sigma;

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return     12.0 * m_v * pow(m_B,-1) * pow(sigmabar,-4.0) *
                       ((2.0 * q2 - 2.0 * m_B2 * sigmabar) * chi_bar_bar_4 -
                        m_B * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * chi_bar_4);
        }

        double I2_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_B2 = power_of<2>(m_B);
            const double sigmabar = 1.0 - sigma;
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_2 = -4.0 * u * (-1.0 + 2.0 * u) / (m_B2 * power_of<3>(sigmabar));

            return C_2 * chi_bar_bar_4;
        }

        double I3_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * u * (-1.0 + 2.0 * u) * (m_v2 + q2 * (5.0 - 4.0 * sigmabar) - m_B2 * sigmabar2)
                             / (m_B2 * power_of<4>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1A_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B),   m_B3 = power_of<3>(m_B),   m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(),   m_v2 = power_of<2>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3   = power_of<3>(sigma), sigma4   = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (m_B3 * (4.0 * omega_1 + omega_2) * sigmabar2 +
                               2.0 * m_B2 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar * (-2.0 + sigmabar) -
                               8.0 * sigma4 * (5.0 * m_B2 * q2 + m_B4 * sigmabar) +
                               4.0 * m_B * sigma3 * (5.0 * (4.0 * omega_1 + omega_2) * q2 + 10.0 * m_B * (m_v2 + q2) +
                               m_B2 * (4.0 * omega_1 + omega_2) * sigmabar - 8.0 * m_B * q2 * sigmabar - 2.0 * m_B3 * sigmabar * (-2.0 + sigmabar)) -
                               m_B * (4.0 * omega_1 + omega_2) * sigmabar * (q2 * sigmabar + m_v2 * (-4.0 + 5.0 * sigmabar)) +
                               4.0 * omega_1 * (2.0 * omega_1 + omega_2) * (2.0 * q2 * sigmabar + m_v2 * (-5.0 + 6.0 * sigmabar)) -
                               2.0 * sigma * (2.0 * m_B4 * sigmabar2 + 2.0 * m_B3 * (4.0 * omega_1 + omega_2) * (-1 + sigmabar) * sigmabar -
                               2.0 * omega_1 * (2.0 * omega_1 + omega_2) * (5.0 * (m_v2 + q2) - 2.0 * q2 * sigmabar) +
                               m_B2 * sigmabar * (-(2.0 * q2 * sigmabar) + omega_1 * (2.0 * omega_1 + omega_2) * (-4.0 + sigmabar) +
                               m_v2 * (8.0 - 10.0 * sigmabar)) + 2.0 * m_B * (4.0 * omega_1 + omega_2) *
                              (3.0 * q2 * sigmabar + m_v2 * (-5.0 + 7.0 * sigmabar))) +
                               sigma2 * (-(20.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2) + 4.0 * m_B4 * sigmabar * (-2.0 + 3.0 * sigmabar) +
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar * (-8.0 + 3.0 * sigmabar) -
                               4.0 * m_B * (4.0 * omega_1 + omega_2) * (5.0 * (m_v2 + q2) - 3.0 * q2 * sigmabar) +
                               4.0 * m_B2 * (-(omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) + 8.0 * q2 * sigmabar +
                               2.0 * m_v2 * (-5.0 + 8.0 * sigmabar))))
                             / (m_B2 * power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1B_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(),   m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_3 = 2.0 * (m_v2 * (4.0 - 4.0 * sigma - 5.0 * sigmabar) - 4.0 * sigma * q2 * sigmabar + sigmabar * (-q2 + m_B2 * sigmabar2))
                             / (m_B * (-omega_1 + m_B * sigma) * power_of<5>(sigmabar));

            return C_3 * chi_bar_bar_4;
        }

        double I3d1C_fT_3pt_chiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4

            return 0.0;
        }

        double I4_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B4     = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2     = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);
            const double u        = (sigma * m_B() - omega_1) / omega_2;

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 6.0 * u * (-1.0 + 2.0 * u) * (m_v4 - 2.0 * m_B2 * q2 * sigmabar3 + m_B4 * sigmabar4 + q2 * q2 * (-1.0 + 2.0 * sigmabar)
                             + 2.0 * m_v2 * sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B2 * power_of<5>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1A_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigma2   = power_of<2>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -(2.0 * (6.0 * (-omega_1 + m_B * sigma) * (2.0 * omega_1 + omega_2 - 2.0 * m_B * sigma) *
                                (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) + (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) *
                                sigmabar) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) -
                                2.0 * m_B * (-omega_1 + m_B * sigma) * sigmabar *
                                (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) + (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) *
                                sigmabar) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) +
                                m_B * (2.0 * omega_1 + omega_2 - 2.0 * m_B * sigma) * sigmabar *
                                (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) + (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) *
                                sigmabar) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) +
                                (-omega_1 + m_B * sigma) * (-(2.0 * omega_1) - omega_2 + 2.0 * m_B * sigma) * sigmabar *
                                (3.0 * m_v4 + 12.0 * m_v2 * q2 * sigmabar - 3.0 * m_B4 * sigmabar3 * (-2.0 + 2.0 * sigma - 3.0 * sigmabar) +
                                m_B2 * (sigmabar2 * (-m_v2 + q2 * sigma) + sigmabar2 * (-(5.0 * m_v2) + q2 * (-11.0 + 10.0 * sigma)) -
                                sigmabar2 * (12.0 * m_v2 + 13.0 * q2 * sigmabar)) +
                                (-(6.0 * sigma2) + sigmabar * (7.0 + 2.0 * sigmabar) + sigma * (3.0 - 4.0 * sigmabar)) * power_of<2>(q2))))
                             / (m_B2 * power_of<2>(omega_2) * power_of<7>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d1B_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar3 = power_of<3>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 2.0 * (3.0 * m_B2 * sigmabar3 - 2.0 * q2 * sigma * sigmabar + q2 * (-3.0 + 2.0 * sigma) * sigmabar +
                               m_v2 * (2.0 - 2.0 * sigma - 5.0 * sigmabar)) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar))
                             / (m_B * power_of<6>(sigmabar) * omega_2);

            return C_4 * chi_bar_bar_4;
        }

        double I4d1C_fT_3pt_chiB_bar_bar_4(const double & /*sigma*/, const double & /*omega_2*/, const double & /*q2*/) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4

            return 0.0;
        }

        double I4d2A_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B), m_B6 = power_of<6>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma), sigma5 = power_of<5>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar), sigmabar4 = power_of<4>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = 4.0 * (6.0 * m_B6 * sigmabar4 + 3.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar3 * (-4.0 + 3.0 * sigmabar) +
                               m_B4 * sigmabar2 * (-(12.0 * sigmabar2 * q2) - 2.0 * m_v2 * sigmabar * (1 + 5.0 * sigmabar) +
                               3.0 * omega_1 * (2.0 * omega_1 + omega_2) * (13.0 + 3.0 * sigmabar * (-5.0 + sigmabar))) -
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar2 *
                               (q2 * sigmabar * (-23 + 11.0 * sigmabar) + m_v2 * (-5.0 + sigmabar * (-18.0 + 5.0 * sigmabar))) +
                               m_B2 * sigmabar * (2.0 * m_v4 * sigmabar * (-2.0 + 5.0 * sigmabar) +
                               m_v2 * (2.0 * sigmabar2 * q2 * (1 + 5.0 * sigmabar) +
                               omega_1 * (2.0 * omega_1 + omega_2) * (-15.0 + sigmabar * (-40 + 19.0 * sigmabar))) +
                               q2 * sigmabar * (6.0 * sigmabar2 * q2 - omega_1 * (2.0 * omega_1 + omega_2) * (66 + sigmabar * (-53 + 5.0 * sigmabar)))) +
                               sigma5 * (-60 * m_B6 * sigmabar2 + 30 * m_B4 * q2 * sigmabar + 84 * m_B2 * power_of<2>(q2)) +
                               sigma3 * (6.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar2 * (-15.0 + 8.0 * sigmabar) +
                               15.0 * m_B3 * (4.0 * omega_1 + omega_2) * sigmabar * (m_v2 + 2.0 * q2 * (1 + sigmabar)) -
                               6.0 * m_B6 * sigmabar2 * (33 + 2.0 * sigmabar * (-24 + 5.0 * sigmabar)) +
                               6.0 * m_B * (4.0 * omega_1 + omega_2) * q2 * (14.0 * m_v2 + q2 * (7.0 - 8.0 * sigmabar)) +
                               2.0 * m_B4 * sigmabar * (-(15.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar) + 30 * m_v2 * (1 + sigmabar) +
                               q2 * (15.0 + (81 - 70 * sigmabar) * sigmabar)) +
                               m_B2 * (84 * m_v4 + 6.0 * m_v2 * q2 * (28 - 57 * sigmabar) +
                               q2 * sigmabar * (15.0 * omega_1 * (2.0 * omega_1 + omega_2) - 2.0 * q2 * (63 + 2.0 * sigmabar))) +
                               42 * omega_1 * (2.0 * omega_1 + omega_2) * power_of<2>(q2)) -
                               m_B * sigma4 * (-30 * m_B4 * (4.0 * omega_1 + omega_2) * sigmabar2 +
                               15.0 * m_B2 * (4.0 * omega_1 + omega_2) * q2 * sigmabar + 60 * m_B5 * sigmabar2 * (-3.0 + 2.0 * sigmabar) +
                               12.0 * m_B * q2 * (14.0 * m_v2 + q2 * (7.0 - 10.0 * sigmabar)) +
                               10.0 * m_B3 * sigmabar * (3.0 * m_v2 + q2 * (6.0 + 5.0 * sigmabar)) + 42 * (4.0 * omega_1 + omega_2) * power_of<2>(q2)
                               ) + omega_1 * (2.0 * omega_1 + omega_2) *
                               (m_v4 * (-42 + 87 * sigmabar) + 3.0 * m_v2 * q2 * sigmabar * (13.0 + 7.0 * sigmabar) +
                               5.0 * sigmabar2 * (5.0 - 2.0 * sigmabar) * power_of<2>(q2)) +
                               m_B * (4.0 * omega_1 + omega_2) * sigmabar *
                               (-(3.0 * m_v2 * q2 * sigmabar * (3.0 + 5.0 * sigmabar)) + 3.0 * m_v4 * (4.0 - 9.0 * sigmabar) +
                               sigmabar2 * (-11.0 + 2.0 * sigmabar) * power_of<2>(q2)) +
                               sigma2 * (6.0 * omega_1 * (2.0 * omega_1 + omega_2) * q2 * (-(7.0 * (2.0 * m_v2 + q2)) + 6.0 * q2 * sigmabar) +
                               9.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar2 * (11.0 + 2.0 * sigmabar * (-6.0 + sigmabar)) +
                               6.0 * m_B6 * sigmabar2 * (13.0 + 3.0 * sigmabar * (-13.0 + 6.0 * sigmabar)) -
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar *
                               (5.0 * m_v2 * (6.0 + 7.0 * sigmabar) + q2 * (15.0 + (91 - 54 * sigmabar) * sigmabar)) -
                               2.0 * m_B4 * sigmabar * (9.0 * omega_1 * (2.0 * omega_1 + omega_2) * sigmabar * (-5.0 + 2.0 * sigmabar) +
                               m_v2 * (15.0 - 54 * sigmabar2 + 20.0 * sigmabar) + q2 * sigmabar * (56 + sigmabar * (-131 + 30 * sigmabar))) +
                               m_B2 * (6.0 * m_v4 * (-14.0 + 37 * sigmabar) -
                               m_v2 * sigmabar * (15.0 * omega_1 * (2.0 * omega_1 + omega_2) + 2.0 * q2 * (-87 + 49 * sigmabar)) -
                               q2 * sigmabar * (2.0 * q2 * sigmabar * (-5.0 + 26 * sigmabar) +
                               5.0 * omega_1 * (2.0 * omega_1 + omega_2) * (6.0 + 7.0 * sigmabar))) -
                               m_B * (4.0 * omega_1 + omega_2) * (42 * m_v4 + 21 * m_v2 * q2 * (4.0 - 7.0 * sigmabar) -
                               sigmabar * (51 + 10.0 * sigmabar) * power_of<2>(q2))) +
                               sigma * (6.0 * m_B6 * sigmabar3 * (8.0 - 9.0 * sigmabar) -
                               3.0 * m_B5 * (4.0 * omega_1 + omega_2) * sigmabar2 * (13.0 + 9.0 * sigmabar * (-3.0 + sigmabar)) +
                               m_B4 * sigmabar2 * (m_v2 * (-20.0 + 30 * sigmabar2 - 68 * sigmabar) +
                               6.0 * q2 * sigmabar * (-15.0 + 11.0 * sigmabar) -
                               9.0 * omega_1 * (2.0 * omega_1 + omega_2) * (11.0 + sigmabar * (-8.0 + sigmabar))) +
                               m_B3 * (4.0 * omega_1 + omega_2) * sigmabar *
                               (m_v2 * (15.0 + (30 - 37 * sigmabar) * sigmabar) + q2 * sigmabar * (61 + 3.0 * sigmabar * (-31 + 5.0 * sigmabar))) +
                               m_B2 * sigmabar * (16.0 * m_v4 * (-3.0 + 7.0 * sigmabar) +
                               m_v2 * (2.0 * q2 * sigmabar * (22 + 25 * sigmabar) +
                               10.0 * omega_1 * (2.0 * omega_1 + omega_2) * (3.0 + 4.0 * sigmabar)) +
                               q2 * (6.0 * sigmabar2 * q2 * (7.0 - 2.0 * sigmabar) -
                               omega_1 * (2.0 * omega_1 + omega_2) * (-15.0 + sigmabar * (-101 + 37 * sigmabar)))) +
                               omega_1 * (2.0 * omega_1 + omega_2) * (42 * m_v4 + 3.0 * m_v2 * q2 * (28 - 41 * sigmabar) -
                               sigmabar * (39 + 16.0 * sigmabar) * power_of<2>(q2)) -
                               m_B * (4.0 * omega_1 + omega_2) * (m_v4 * (-42 + 99 * sigmabar) + 3.0 * m_v2 * q2 * sigmabar * (21 - 4.0 * sigmabar) +
                               2.0 * sigmabar2 * (8.0 - 9.0 * sigmabar) * power_of<2>(q2))))
                             / (m_B2 * power_of<2>(omega_2) * power_of<8>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2B_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B3 = power_of<3>(m_B), m_B4 = power_of<4>(m_B), m_B5 = power_of<5>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v), m_v4 = power_of<4>(m_v);
            const double sigma2   = power_of<2>(sigma), sigma3 = power_of<3>(sigma), sigma4 = power_of<4>(sigma);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return  pow(m_B,-1) * pow(omega_1 - m_B * sigma,-2.0) * pow(sigmabar,-7.0) *
                    (4.0 * chi_bar_bar_4 * (-(m_B * sigmabar * (-m_v2 + (m_B2 - q2) * sigmabar) *
                    (3.0 * (m_B2 - q2) * sigmabar + m_v2 * (2.0 - 5.0 * sigmabar))) +
                    m_B * sigma * (m_v4 * (12.0 - 29 * sigmabar) +
                    2.0 * sigmabar2 * (m_B2 - q2) * (q2 * (5.0 - 2.0 * sigmabar) + m_B2 * (-6.0 + 9.0 * sigmabar)) +
                    m_v2 * sigmabar * (-(q2 * (13.0 + 10.0 * sigmabar)) + m_B2 * (5.0 + 2.0 * sigmabar * (8.0 - 5.0 * sigmabar)))) +
                    omega_1 * (3.0 * m_v4 * (-4.0 + 9.0 * sigmabar) -
                    sigmabar2 * (m_B2 - q2) * (q2 * (11.0 - 2.0 * sigmabar) + 3.0 * m_B2 * (-4.0 + 3.0 * sigmabar)) +
                    m_v2 * sigmabar * (3.0 * q2 * (3.0 + 5.0 * sigmabar) + m_B2 * (-5.0 + sigmabar * (-18.0 + 5.0 * sigmabar)))) +
                    sigma4 * (12.0 * m_B5 * sigmabar2 - 5.0 * m_B3 * q2 * sigmabar - 12.0 * m_B * power_of<2>(q2)) +
                    sigma3 * (-(12.0 * m_B4 * omega_1 * sigmabar2) + 5.0 * m_B2 * omega_1 * q2 * sigmabar +
                    12.0 * m_B5 * sigmabar2 * (-3.0 + sigmabar) + 4.0 * m_B * q2 * (6.0 * m_v2 + 3.0 * q2 - 2.0 * q2 * sigmabar) +
                    m_B3 * sigmabar * (5.0 * m_v2 + 2.0 * q2 * (5.0 + 8.0 * sigmabar)) + 12.0 * omega_1 * power_of<2>(q2)) +
                    sigma2 * (-(9.0 * m_B4 * omega_1 * sigmabar2 * (-4.0 + sigmabar)) +
                    9.0 * m_B5 * sigmabar2 * (4.0 - 3.0 * sigmabar) +
                    6.0 * omega_1 * q2 * (-(4.0 * m_v2) + q2 * (-2.0 + sigmabar)) -
                    m_B2 * omega_1 * sigmabar * (5.0 * m_v2 + q2 * (10.0 + 17.0 * sigmabar)) +
                    m_B3 * sigmabar * (-(m_v2 * (10.0 + 17.0 * sigmabar)) + q2 * (-5.0 + sigmabar * (-38 + 15.0 * sigmabar))) +
                    m_B * (-(12.0 * m_v4) + m_v2 * q2 * (-24 + 37 * sigmabar) + sigmabar * (11.0 + 8.0 * sigmabar) * power_of<2>(q2))) +
                    omega_1 * sigma * (12.0 * m_v4 + m_v2 * (q2 * (24 - 33 * sigmabar) + 2.0 * m_B2 * sigmabar * (5.0 + 9.0 * sigmabar)) +
                    sigmabar * (18.0 * m_B4 * sigmabar * (-2.0 + sigmabar) +
                    5.0 * m_B2 * q2 * (1 - 2.0 * sigmabar * (-4.0 + sigmabar)) - (9.0 + 8.0 * sigmabar) * power_of<2>(q2)))) +
                    2.0 * m_B * (-omega_1 + m_B * sigma) * sigmabar *
                    (3.0 * m_B2 * sigmabar3 - 2.0 * q2 * sigma * sigmabar + q2 * (-3.0 + 2.0 * sigma) * sigmabar +
                    m_v2 * (2.0 - 2.0 * sigma - 5.0 * sigmabar)) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * chi_bar_4);
        }

        double I4d2C_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;

            const double m_B2     = power_of<2>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);

            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            const double C_4 = -(2.0 * (-(2.0 * (m_v2 - q2 * sigma) * sigmabar) +
                                (5.0 * m_v2 - 3.0 * m_B2 * sigmabar2 + q2 * (3.0 - 2.0 * sigma)) * sigmabar) *
                                (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)))
                             /  (power_of<2>(omega_2) * power_of<6>(sigmabar));

            return C_4 * chi_bar_bar_4;
        }

        double I4d2D_fT_3pt_chiB_bar_bar_4(const double & sigma, const double & q2) const
        {
            // three-particle contribution to f_T proportional to chi_bar_bar_4
            const double omega_1  = m_B * sigma;
            const double omega_2  = m_B * sigma - omega_1;

            const double m_B2     = power_of<2>(m_B), m_B4 = power_of<4>(m_B);
            const double m_v      = this->m_v(), m_v2 = power_of<2>(m_v);
            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar), sigmabar3 = power_of<3>(sigmabar);

            const double chi_bar_4     = this->chi_bar_4(omega_1, omega_2);
            const double chi_bar_bar_4 = this->chi_bar_bar_4(omega_1, omega_2);

            return   2.0 * pow(m_B,-2.0) * pow(sigmabar,-7.0) * (chi_bar_bar_4 *
                     ((m_v2 - q2 * sigma) * (-(7.0 * m_v2) + m_B2 * sigmabar2 - 5.0 * q2 + 6.0 * q2 * sigma) * sigmabar +
                     sigmabar2 * (q2 * (5.0 * m_v2 + q2) - m_B2 * (2.0 * m_v2 + q2 - 3.0 * q2 * sigma) * sigmabar) +
                     sigmabar3 * (9.0 * m_B4 * sigmabar2 + m_B2 * (-(5.0 * m_v2) - 11.0 * q2 + 10.0 * q2 * sigma) +
                     2.0 * power_of<2>(q2)) + 4.0 * sigmabar * power_of<2>(m_v2 - q2 * sigma)) +
                     m_B * sigmabar * (3.0 * m_B2 * sigmabar3 - 2.0 * q2 * sigma * sigmabar + q2 * (-3.0 + 2.0 * sigma) * sigmabar +
                     m_v2 * (2.0 - 2.0 * sigma - 5.0 * sigmabar)) * (m_v2 - q2 * sigma + sigmabar * (q2 - m_B2 * sigmabar)) * chi_bar_4);
        }

        double I2_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + 0.0                                                        - I2_fT_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I3_fT_3pt_chiA_bar_bar_4(sigma, omega_1, omega_2, q2)     - I3_fT_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1A_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I3d1A_fT_3pt_chiA_bar_bar_4(sigma, omega_1, omega_2, q2)  - I3d1A_fT_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I3d1B_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return + I3d1B_fT_3pt_chiA_bar_bar_4(sigma, omega_1, q2)           - I3d1B_fT_3pt_chiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I3d1C_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return + I3d1C_fT_3pt_chiA_bar_bar_4(sigma, omega_2, q2)           - I3d1C_fT_3pt_chiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I4_fT_3pt_chiA_bar_bar_4( sigma, omega_1, omega_2, q2)    - I4_fT_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1A_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I4d1A_fT_3pt_chiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d1A_fT_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d1B_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return + I4d1B_fT_3pt_chiA_bar_bar_4( sigma, omega_1, q2)          - I4d1B_fT_3pt_chiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d1C_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return + I4d1C_fT_3pt_chiA_bar_bar_4( sigma, omega_2, q2)          - I4d1C_fT_3pt_chiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2A_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & omega_2, const double & q2) const
        {
            return + I4d2A_fT_3pt_chiA_bar_bar_4( sigma, omega_1, omega_2, q2) - I4d2A_fT_3pt_chiB_bar_bar_4(sigma, omega_1, omega_2, q2);
        }

        double I4d2B_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_1, const double & q2) const
        {
            return + I4d2B_fT_3pt_chiA_bar_bar_4( sigma, omega_1, q2)          - I4d2B_fT_3pt_chiB_bar_bar_4(sigma, omega_1, q2);
        }

        double I4d2C_fT_3pt_chi_bar_bar_4(const double & sigma, const double & omega_2, const double & q2) const
        {
            return + I4d2C_fT_3pt_chiA_bar_bar_4( sigma, omega_2, q2)          - I4d2C_fT_3pt_chiB_bar_bar_4(sigma, omega_2, q2);
        }

        double I4d2D_fT_3pt_chi_bar_bar_4(const double & sigma, const double & q2) const
        {
            return + I4d2D_fT_3pt_chiA_bar_bar_4( sigma, q2)                   - I4d2D_fT_3pt_chiB_bar_bar_4(sigma, q2);
        }
        // }}}

        /* fT : integrands and surface terms */
        // {{{
        double integrand_fT_2pt_disp(const double & sigma, const double & q2) const
        {
            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);
            const double etad3    = 24.0 * (eta - 1.0) * power_of<2>(eta) * (2.0 * eta - 1.0) / power_of<3>(sigmabar);

            const double I1   = I1_fT_2pt_phi_bar(sigma, q2);
            const double I2   = I2_fT_2pt_phi_bar(sigma, q2)   + I2_fT_2pt_g_bar(sigma, q2);
            const double I2d1 = I2d1_fT_2pt_phi_bar(sigma, q2) + I2d1_fT_2pt_g_bar(sigma, q2);
            const double I3   = I3_fT_2pt_g_bar(sigma, q2);
            const double I3d1 = I3d1_fT_2pt_g_bar(sigma, q2);
            const double I3d2 = I3d2_fT_2pt_g_bar(sigma, q2);
            const double I4   = I4_fT_2pt_g_bar(sigma, q2);
            const double I4d1 = I4d1_fT_2pt_g_bar(sigma, q2);
            const double I4d2 = I4d2_fT_2pt_g_bar(sigma, q2);
            const double I4d3 = I4d3_fT_2pt_g_bar(sigma, q2);

            double result = 0.0;
            result += -1.0 * I1;
            result += (etad1 * I2 + eta * I2d1) / m_B2;
            result += -1.0 * (I3 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I3d1 * eta * etad1 + I3d2 * power_of<2>(eta)) / (2.0 * m_B4);
            result += I4 * (power_of<2>(eta) * etad3 + 4.0 * eta * etad1 * etad2 + power_of<3>(etad1)) / (6.0 * m_B6);
            result += I4d1 * eta * (4.0 * eta * etad2 + 7.0 * power_of<2>(etad1)) / (6.0 * m_B6);
            result += I4d2 * 6.0 * power_of<2>(eta) * etad1 / (6.0 * m_B6);
            result += I4d3 * power_of<3>(eta) / (6.0 * m_B6);
            result *= exp;

            return result;
        }

        double integrand_fT_2pt_borel(const double & sigma, const double & q2) const
        {
            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double I1   = I1_fT_2pt_phi_bar(sigma, q2);
            const double I2   = I2_fT_2pt_phi_bar(sigma, q2)   + I2_fT_2pt_g_bar(sigma, q2);
            const double I3   = I3_fT_2pt_g_bar(sigma, q2);
            const double I4   = I4_fT_2pt_g_bar(sigma, q2);

            double result = 0.0;
            result += - I1;
            result +=   I2 / M2;
            result += - I3 / (2 * M4);
            result +=   I4 / (6 * M6);
            result *= exp;

            return result;
        }

        double surface_fT_2pt(const double & sigma, const double & q2) const
        {
            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fT_2pt_phi_bar(sigma, q2) + I2_fT_2pt_g_bar(sigma, q2);
            const double I3   = I3_fT_2pt_g_bar(sigma, q2);
            const double I3d1 = I3d1_fT_2pt_g_bar(sigma, q2);
            const double I4   = I4_fT_2pt_g_bar(sigma, q2);
            const double I4d1 = I4d1_fT_2pt_g_bar(sigma, q2);
            const double I4d2 = I4d2_fT_2pt_g_bar(sigma, q2);


            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / (power_of<2>( M2())));
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;

            return result;
        }

        /*
         * rewrite integration ranges such that:
         * 1.)
         *    0 <= x_1 <= 1,     and    0 <= x_2 <= 1,
         * 2.)
         *    x_1 and x_2 integration boundaries do not depend on the other variables
         *
         * We obtain the integrand
         *
         *    sigma m_B f(sigma m_B x_1, sigma m_B (xbar_1 xbar_2 + x_2) / xbar_2) / (xbar_1 xbar_2^2 + x_2 xbar_2),
         *
         * where
         *
         *    xbar_1 = 1.0 - x_1,    and    xbar_2 = 1.0 - x_2.
         */
        double integrand_fT_3pt(const std::array<double, 3> & args, const double & q2) const
        {
            const double sigma  = args[0];
            const double x_1    = args[1];
            const double x_2    = args[2];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double I1 = I1_fT_3pt_phi_3(sigma, omega_1, omega_2, q2);
            const double I2 = I2_fT_3pt_phi_3(sigma, omega_1, omega_2, q2)         + I2_fT_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                            + I2_fT_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I2_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                            + I2_fT_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I2_fT_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                            + I2_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I2_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I3 = I3_fT_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)     + I3_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)
                            + I3_fT_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                            + I3_fT_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fT_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                            + I3_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4 = I4_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                            + I4_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result = 0.0;
            result += - I1;
            result +=   I2 / M2;
            result += - I3 / (2.0 * M4);
            result +=   I4 / (6.0 * M6);
            result *=   prefactor * exp;

            return result;
        }

        double surface_fT_3pt_A(const std::array<double, 2> & args, const double & sigma, const double & q2) const
        {
            const double x_1    = args[0];
            const double x_2    = args[1];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fT_3pt_phi_3(sigma, omega_1, omega_2, q2)            + I2_fT_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I2_fT_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)        + I2_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I2_fT_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I2_fT_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I2_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I2_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I3   = I3_fT_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)        + I3_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)
                              + I3_fT_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fT_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fT_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I3_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I3d1 = I3d1A_fT_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)     + I3d1A_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)
                              + I3d1A_fT_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fT_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fT_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3d1A_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4   = I4_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)    + I4_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I4_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d1 = I4d1A_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d1A_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d1A_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d1A_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d2 = I4d2A_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d2A_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d2A_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d2A_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        double surface_fT_3pt_B(const double & x_1, const double & sigma, const double & q2) const
        {
            // this ONLY includes the Jacobian from the transformation (omega_1 -> x_1).
            const double prefactor = sigma * m_B();

            const double omega_1 = sigma * m_B() * x_1;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1B_fT_3pt_phi_bar_3(sigma, omega_1, q2)     + I3d1B_fT_3pt_phi_bar_bar_3(sigma, omega_1, q2)
                                  + I3d1B_fT_3pt_phi_bar_4(sigma, omega_1, q2)     + I3d1B_fT_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fT_3pt_psi_bar_4(sigma, omega_1, q2)     + I3d1B_fT_3pt_chi_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fT_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I3d1B_fT_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1B_fT_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d1B_fT_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d1B_fT_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d1B_fT_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            const     double I4d2 = I4d2B_fT_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d2B_fT_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d2B_fT_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d2B_fT_3pt_chi_bar_bar_4(sigma, omega_1, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        double surface_fT_3pt_C(const double & x_2, const double & sigma, const double & q2) const
        {
            const double xbar_2 = 1.0 - x_2;

            // this ONLY includes the Jacobian from the transformation (omega_2 -> x_2).
            const double prefactor = sigma * m_B() / (xbar_2 * xbar_2);

            const double omega_2 = sigma * m_B() * (x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1C_fT_3pt_phi_bar_3(sigma, omega_2, q2)     + I3d1C_fT_3pt_phi_bar_bar_3(sigma, omega_2, q2)
                                  + I3d1C_fT_3pt_phi_bar_4(sigma, omega_2, q2)     + I3d1C_fT_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fT_3pt_psi_bar_4(sigma, omega_2, q2)     + I3d1C_fT_3pt_chi_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fT_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I3d1C_fT_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1C_fT_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d1C_fT_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d1C_fT_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d1C_fT_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            const     double I4d2 = I4d2C_fT_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d2C_fT_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d2C_fT_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d2C_fT_3pt_chi_bar_bar_4(sigma, omega_2, q2);

            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        double surface_fT_3pt_D(const double & sigma, const double & q2) const
        {
            // this does NOT includes the original factor of 1 / omega_2
            const double prefactor = 1.0;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);


            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            constexpr double I3d1 = 0.0;
            constexpr double I4   = 0.0;
            constexpr double I4d1 = 0.0;
            const     double I4d2 = I4d2D_fT_3pt_phi_bar_bar_3(sigma, q2) + I4d2D_fT_3pt_phi_bar_bar_4(sigma, q2)
                                  + I4d2D_fT_3pt_psi_bar_bar_4(sigma, q2) + I4d2D_fT_3pt_chi_bar_bar_4(sigma, q2);


            double result = 0.0;
            result += -1.0 * eta * I2 / m_B2;
            result +=  0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result += -1.0 / 6.0 * eta / m_B2 * (I4 / M4);
            result += -1.0 / 6.0 * eta / (m_B4 * M2() ) * (eta * I4d1 + I4 * etad1);
            result += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result *= exp;
            result *= prefactor;

            return result;
        }

        /*
         * Integrands for the first moments. Only the borel method is implemented
         */

        double integrand_fT_2pt_borel_m1(const double & sigma, const double & q2) const
        {
            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double I1   = I1_fT_2pt_phi_bar(sigma, q2);
            const double I2   = I2_fT_2pt_phi_bar(sigma, q2)   + I2_fT_2pt_g_bar(sigma, q2);
            const double I3   = I3_fT_2pt_g_bar(sigma, q2);
            const double I4   = I4_fT_2pt_g_bar(sigma, q2);

            double result1 = 0.0;
            result1 += - I1;
            result1 +=   I2 / M2;
            result1 += - I3 / (2.0 * M4);
            result1 +=   I4 / (6.0 * M6);
            result1 *=   exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - I2;
            result2 +=   I3 / M2;
            result2 += - I4 / (2.0 * M4);
            result2 *=   exp;

            return result1 + result2;
        }

        double surface_fT_2pt_m1(const double & sigma, const double & q2) const
        {
            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double m_P2 = power_of<2>(m_P());
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fT_2pt_phi_bar(sigma, q2)   + I2_fT_2pt_g_bar(sigma, q2);
            const double I3   = I3_fT_2pt_g_bar(sigma, q2);
            const double I3d1 = I3d1_fT_2pt_g_bar(sigma, q2);
            const double I4   = I4_fT_2pt_g_bar(sigma, q2);
            const double I4d1 = I4d1_fT_2pt_g_bar(sigma, q2);
            const double I4d2 = I4d2_fT_2pt_g_bar(sigma, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return result1 + result2;
        }

        double integrand_fT_3pt_m1(const std::array<double, 3> & args, const double & q2) const
        {
            const double sigma  = args[0];
            const double x_1    = args[1];
            const double x_2    = args[2];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_P2 = power_of<2>(m_P());
            const double M4   = power_of<2>(M2), M6 = power_of<3>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double I1 = I1_fT_3pt_phi_3(sigma, omega_1, omega_2, q2);
            const double I2 = I2_fT_3pt_phi_3(sigma, omega_1, omega_2, q2)         + I2_fT_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                            + I2_fT_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)
                            + I2_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                            + I2_fT_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I2_fT_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                            + I2_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I2_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I3 = I3_fT_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)     + I3_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)
                            + I3_fT_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                            + I3_fT_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3_fT_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                            + I3_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4 = I4_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                            + I4_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result1 = 0.0;
            result1 += - I1;
            result1 +=   I2 / M2;
            result1 += - I3 / (2.0 * M4);
            result1 +=   I4 / (6.0 * M6);
            result1 *=   exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - I2;
            result2 +=   I3 / M2;
            result2 += - I4 / (2.0 * M4);
            result2 *=   exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fT_3pt_A_m1(const std::array<double, 2> & args, const double & sigma, const double & q2) const
        {
            const double x_1    = args[0];
            const double x_2    = args[1];
            const double xbar_1 = 1.0 - x_1;
            const double xbar_2 = 1.0 - x_2;

            // this includes the original factor of 1 / omega_2 (which corresponds to 1 / xi in the notation of
            // [KMO:2006A]), as well as the Jacobian from the transformation (omega_1, omega_2 -> x_1, x_2).
            const double prefactor = sigma * m_B() / ((xbar_1 * xbar_2 + x_2) * xbar_2);

            const double omega_1 = sigma * m_B() * x_1;
            const double omega_2 = sigma * m_B() * (xbar_1 + x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            const double I2   = I2_fT_3pt_phi_3(sigma, omega_1, omega_2, q2)            + I2_fT_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)
                              + I2_fT_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)
                              + I2_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I2_fT_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I2_fT_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I2_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I2_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I3   = I3_fT_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)        + I3_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)
                              + I3_fT_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fT_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)        + I3_fT_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I3_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I3d1 = I3d1A_fT_3pt_phi_bar_3(sigma, omega_1, omega_2, q2)     + I3d1A_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)
                              + I3d1A_fT_3pt_phi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fT_3pt_psi_bar_4(sigma, omega_1, omega_2, q2)     + I3d1A_fT_3pt_chi_bar_4(sigma, omega_1, omega_2, q2)
                              + I3d1A_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I3d1A_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4   = I4_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2)    + I4_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2)    + I4_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d1 = I4d1A_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d1A_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d1A_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d1A_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);
            const double I4d2 = I4d2A_fT_3pt_phi_bar_bar_3(sigma, omega_1, omega_2, q2) + I4d2A_fT_3pt_phi_bar_bar_4(sigma, omega_1, omega_2, q2)
                              + I4d2A_fT_3pt_psi_bar_bar_4(sigma, omega_1, omega_2, q2) + I4d2A_fT_3pt_chi_bar_bar_4(sigma, omega_1, omega_2, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fT_3pt_B_m1(const double & x_1, const double & sigma, const double & q2) const
        {
            // this ONLY includes the Jacobian from the transformation (omega_1 -> x_1).
            const double prefactor = sigma * m_B();

            const double omega_1 = sigma * m_B() * x_1;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1B_fT_3pt_phi_bar_3(sigma, omega_1, q2)     + I3d1B_fT_3pt_phi_bar_bar_3(sigma, omega_1, q2)
                                  + I3d1B_fT_3pt_phi_bar_4(sigma, omega_1, q2)     + I3d1B_fT_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fT_3pt_psi_bar_4(sigma, omega_1, q2)     + I3d1B_fT_3pt_chi_bar_4(sigma, omega_1, q2)
                                  + I3d1B_fT_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I3d1B_fT_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1B_fT_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d1B_fT_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d1B_fT_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d1B_fT_3pt_chi_bar_bar_4(sigma, omega_1, q2);
            const     double I4d2 = I4d2B_fT_3pt_phi_bar_bar_3(sigma, omega_1, q2) + I4d2B_fT_3pt_phi_bar_bar_4(sigma, omega_1, q2)
                                  + I4d2B_fT_3pt_psi_bar_bar_4(sigma, omega_1, q2) + I4d2B_fT_3pt_chi_bar_bar_4(sigma, omega_1, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fT_3pt_C_m1(const double & x_2, const double & sigma, const double & q2) const
        {
            const double xbar_2 = 1.0 - x_2;

            // this ONLY includes the Jacobian from the transformation (omega_2 -> x_2).
            const double prefactor = sigma * m_B() / (xbar_2 * xbar_2);

            const double omega_2 = sigma * m_B() * (x_2 / xbar_2);

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);

            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            const     double I3d1 = I3d1C_fT_3pt_phi_bar_3(sigma, omega_2, q2)     + I3d1C_fT_3pt_phi_bar_bar_3(sigma, omega_2, q2)
                                  + I3d1C_fT_3pt_phi_bar_4(sigma, omega_2, q2)     + I3d1C_fT_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fT_3pt_psi_bar_4(sigma, omega_2, q2)     + I3d1C_fT_3pt_chi_bar_4(sigma, omega_2, q2)
                                  + I3d1C_fT_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I3d1C_fT_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            constexpr double I4   = 0.0;
            const     double I4d1 = I4d1C_fT_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d1C_fT_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d1C_fT_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d1C_fT_3pt_chi_bar_bar_4(sigma, omega_2, q2);
            const     double I4d2 = I4d2C_fT_3pt_phi_bar_bar_3(sigma, omega_2, q2) + I4d2C_fT_3pt_phi_bar_bar_4(sigma, omega_2, q2)
                                  + I4d2C_fT_3pt_psi_bar_bar_4(sigma, omega_2, q2) + I4d2C_fT_3pt_chi_bar_bar_4(sigma, omega_2, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }

        double surface_fT_3pt_D_m1(const double & sigma, const double & q2) const
        {
            // this does NOT includes the original factor of 1 / omega_2
            const double prefactor = 1.0;

            const double m_B2 = power_of<2>(m_B()), m_B4 = power_of<4>(m_B()), m_B6 = power_of<6>(m_B());
            const double m_P2 = power_of<2>(m_P());
            const double m_v2 = power_of<2>(m_v());
            const double M4   = power_of<2>(M2);
            const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

            const double sigmabar = 1.0 - sigma, sigmabar2 = power_of<2>(sigmabar);
            const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
            const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;
            const double etad2    = 2.0 * (eta - 1.0) * eta * (4.0 * eta - 1.0) / power_of<2>(sigmabar);


            constexpr double I2   = 0.0;
            constexpr double I3   = 0.0;
            constexpr double I3d1 = 0.0;
            constexpr double I4   = 0.0;
            constexpr double I4d1 = 0.0;
            const     double I4d2 = I4d2D_fT_3pt_phi_bar_bar_3(sigma, q2) + I4d2D_fT_3pt_phi_bar_bar_4(sigma, q2)
                                  + I4d2D_fT_3pt_psi_bar_bar_4(sigma, q2) + I4d2D_fT_3pt_chi_bar_bar_4(sigma, q2);

            double result1 = 0.0;
            result1 += -1.0 * eta * I2 / m_B2;
            result1 +=  0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
            result1 += -1.0 / 6.0 * eta / m_B2 * (I4 / (M4));
            result1 += -1.0 / 6.0 * eta / (m_B4 * M2 ) * (eta * I4d1 + I4 * etad1);
            result1 += -1.0 / 6.0 * eta / m_B6 * (I4 * (power_of<2>(etad1) + eta * etad2) + 3.0 * I4d1 * eta * etad1 + I4d2 * power_of<2>(eta));
            result1 *=  exp * s(sigma, q2);

            double result2 = 0.0;
            result2 += - 0.5 * eta * I3 / m_B2;
            result2 += + eta * I4 / (3.0 * M2 * m_B2) ;
            result2 += + eta * (eta * I4d1 + I4 * etad1) / (6.0 * m_B4);
            result2 *= exp;

            return (result1 + result2) * prefactor;
        }
        // }}}

        /* fT : form factor and moments */
        // {{{
        double f_t(const double & q2) const
        {
            const double sigma_0 = this->sigma_0(q2, s0_0_t(), s0_1_t());

            const std::function<double (const double &)> integrand_2pt = std::bind(integrand_fT_2pt, this, std::placeholders::_1, q2);

            const double integral_2pt = integrate<GSL::QAGS>(integrand_2pt, 0.0, sigma_0);
            const double surface_2pt  = 0.0 - surface_fT_2pt(switch_borel ? sigma_0 : 0.0, q2);

            double integral_3pt = 0.0;
            double surface_3pt  = 0.0;

            if (switch_3pt != 0.0)
            {
                const std::function<double (const double &)> surface_3pt_B = std::bind(&Implementation::surface_fT_3pt_B, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const double &)> surface_3pt_C = std::bind(&Implementation::surface_fT_3pt_C, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const std::array<double, 3> &)> integrand_3pt = std::bind(&Implementation::integrand_fT_3pt, this, std::placeholders::_1, q2);
                const std::function<double (const std::array<double, 2> &)> surface_3pt_A = std::bind(&Implementation::surface_fT_3pt_A, this, std::placeholders::_1, sigma_0, q2);

                integral_3pt = integrate(integrand_3pt, { 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0 }, cubature::Config());
                surface_3pt  = 0.0
                             - integrate(surface_3pt_A, { 0.0, 0.0 }, { 1.0, 1.0 }, cubature::Config()) // integrate over x_1 and x_2
                             - integrate<GSL::QAGS>(surface_3pt_B, 0.0, 1.0)                            // integrate over x_1
                             - integrate<GSL::QAGS>(surface_3pt_C, 0.0, 1.0)                            // integrate over x_2
                             - surface_fT_3pt_D(sigma_0, q2);
            }

            return f_B() * power_of<2>(m_B()) * (m_B() + m_P()) / (f_P() * (power_of<2>(m_B()) - power_of<2>(m_P()) - q2)) * (integral_2pt + surface_2pt + integral_3pt + surface_3pt) / ( Traits::chi2);
        }

        double normalized_moment_1_f_t(const double & q2) const
        {
            const double sigma_0 = this->sigma_0(q2, s0_0_t(), s0_1_t());

            const std::function<double (const double &)> integrand_2pt_m1 = std::bind(&Implementation::integrand_fT_2pt_borel_m1, this, std::placeholders::_1, q2);


            const std::function<double (const double &)> integrand_2pt    = std::bind(&Implementation::integrand_fT_2pt_borel, this, std::placeholders::_1, q2);

            const double integral_2pt_m1 = integrate<GSL::QAGS>(integrand_2pt_m1, 0.0, sigma_0);
            const double surface_2pt_m1  = 0.0 - surface_fT_2pt_m1(sigma_0, q2);

            double integral_3pt_m1 = 0.0;
            double surface_3pt_m1  = 0.0;

            if (switch_3pt != 0.0)
            {
                const std::function<double (const double &)> surface_3pt_B_m1 = std::bind(&Implementation::surface_fT_3pt_B_m1, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const double &)> surface_3pt_C_m1 = std::bind(&Implementation::surface_fT_3pt_C_m1, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const std::array<double, 3> &)> integrand_3pt_m1 = std::bind(&Implementation::integrand_fT_3pt_m1, this, std::placeholders::_1, q2);
                const std::function<double (const std::array<double, 2> &)> surface_3pt_A_m1 = std::bind(&Implementation::surface_fT_3pt_A_m1, this, std::placeholders::_1, sigma_0, q2);

                integral_3pt_m1 = integrate(integrand_3pt_m1, { 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0 }, cubature::Config());
                surface_3pt_m1  = 0.0
                                - integrate(surface_3pt_A_m1, { 0.0, 0.0 }, { 1.0, 1.0 }, cubature::Config()) // integrate over x_1 and x_2
                                - integrate<GSL::QAGS>(surface_3pt_B_m1, 0.0, 1.0)                            // integrate over x_1
                                - integrate<GSL::QAGS>(surface_3pt_C_m1, 0.0, 1.0)                            // integrate over x_2
                                - surface_fT_3pt_D_m1(sigma_0, q2);
            }
            const double numerator       = integral_2pt_m1 + surface_2pt_m1 + integral_3pt_m1 + surface_3pt_m1;

            const double integral_2pt    = integrate<GSL::QAGS>(integrand_2pt, 0.0, sigma_0);
            const double surface_2pt     = 0.0 - surface_fT_2pt(sigma_0, q2);

            double integral_3pt    = 0.0;
            double surface_3pt     = 0.0;

            if (switch_3pt != 0.0)
            {
                const std::function<double (const double &)> surface_3pt_B    = std::bind(&Implementation::surface_fT_3pt_B, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const double &)> surface_3pt_C    = std::bind(&Implementation::surface_fT_3pt_C, this, std::placeholders::_1, sigma_0, q2);
                const std::function<double (const std::array<double, 3> &)> integrand_3pt = std::bind(&Implementation::integrand_fT_3pt, this, std::placeholders::_1, q2);
                const std::function<double (const std::array<double, 2> &)> surface_3pt_A = std::bind(&Implementation::surface_fT_3pt_A, this, std::placeholders::_1, sigma_0, q2);

                integral_3pt    = integrate(integrand_3pt, { 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0 }, cubature::Config());
                surface_3pt     = 0.0
                                - integrate(surface_3pt_A, { 0.0, 0.0 }, { 1.0, 1.0 }, cubature::Config()) // integrate over x_1 and x_2
                                - integrate<GSL::QAGS>(surface_3pt_B, 0.0, 1.0)                            // integrate over x_1
                                - integrate<GSL::QAGS>(surface_3pt_C, 0.0, 1.0)                            // integrate over x_2
                                - surface_fT_3pt_D(sigma_0, q2);
            }
            const double denominator     = integral_2pt + surface_2pt + integral_3pt + surface_3pt;

            return numerator / denominator;
        }
        // }}}

        /* Diagnostics */

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            /* dependent variables */
            results.add({ this->m_v(),                       "m_v(mu) in the MSbar scheme"                         });

            results.add({this->s0_0_t(),                      "s_0 value for fT"                                    });

            /* f_+ */

            /* 2 particle */

            /* I_1 phi_+ */
            results.add({ this->I1_fp_2pt_phi_p(0.04, -5.0), "f_+: I_1^{2pt,phi_+}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I1_fp_2pt_phi_p(0.04,  0.0), "f_+: I_1^{2pt,phi_+}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I1_fp_2pt_phi_p(0.04, +5.0), "f_+: I_1^{2pt,phi_+}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I1_fp_2pt_phi_p(0.08, -5.0), "f_+: I_1^{2pt,phi_+}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I1_fp_2pt_phi_p(0.08,  0.0), "f_+: I_1^{2pt,phi_+}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I1_fp_2pt_phi_p(0.08, +5.0), "f_+: I_1^{2pt,phi_+}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2 phi_bar */
            results.add({ this->I2_fp_2pt_phi_bar(0.04, -5.0), "f_+: I_2^{2pt,phi_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fp_2pt_phi_bar(0.04,  0.0), "f_+: I_2^{2pt,phi_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fp_2pt_phi_bar(0.04, +5.0), "f_+: I_2^{2pt,phi_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2_fp_2pt_phi_bar(0.08, -5.0), "f_+: I_2^{2pt,phi_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fp_2pt_phi_bar(0.08,  0.0), "f_+: I_2^{2pt,phi_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fp_2pt_phi_bar(0.08, +5.0), "f_+: I_2^{2pt,phi_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2d1 phi_bar */
            results.add({ this->I2d1_fp_2pt_phi_bar(0.04, -5.0), "f_+: I_2d1^{2pt,phi_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fp_2pt_phi_bar(0.04,  0.0), "f_+: I_2d1^{2pt,phi_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fp_2pt_phi_bar(0.04, +5.0), "f_+: I_2d1^{2pt,phi_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2d1_fp_2pt_phi_bar(0.08, -5.0), "f_+: I_2d1^{2pt,phi_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fp_2pt_phi_bar(0.08,  0.0), "f_+: I_2d1^{2pt,phi_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fp_2pt_phi_bar(0.08, +5.0), "f_+: I_2d1^{2pt,phi_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2 g_+ */
            results.add({ this->I2_fp_2pt_g_p(0.04, -5.0), "f_+: I_2^{2pt,g_p}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fp_2pt_g_p(0.04,  0.0), "f_+: I_2^{2pt,g_p}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fp_2pt_g_p(0.04, +5.0), "f_+: I_2^{2pt,g_p}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2_fp_2pt_g_p(0.08, -5.0), "f_+: I_2^{2pt,g_p}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fp_2pt_g_p(0.08,  0.0), "f_+: I_2^{2pt,g_p}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fp_2pt_g_p(0.08, +5.0), "f_+: I_2^{2pt,g_p}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2d1 g_+ */
            results.add({ this->I2d1_fp_2pt_g_p(0.04, -5.0), "f_+: I_2d1^{2pt,g_p}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fp_2pt_g_p(0.04,  0.0), "f_+: I_2d1^{2pt,g_p}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fp_2pt_g_p(0.04, +5.0), "f_+: I_2d1^{2pt,g_p}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2d1_fp_2pt_g_p(0.08, -5.0), "f_+: I_2d1^{2pt,g_p}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fp_2pt_g_p(0.08,  0.0), "f_+: I_2d1^{2pt,g_p}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fp_2pt_g_p(0.08, +5.0), "f_+: I_2d1^{2pt,g_p}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3 g_+ */
            results.add({ this->I3_fp_2pt_g_p(0.04, -5.0), "f_+: I_3^{2pt,g_p}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3_fp_2pt_g_p(0.04,  0.0), "f_+: I_3^{2pt,g_p}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3_fp_2pt_g_p(0.04, +5.0), "f_+: I_3^{2pt,g_p}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3_fp_2pt_g_p(0.08, -5.0), "f_+: I_3^{2pt,g_p}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3_fp_2pt_g_p(0.08,  0.0), "f_+: I_3^{2pt,g_p}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3_fp_2pt_g_p(0.08, +5.0), "f_+: I_3^{2pt,g_p}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3d1 g_+ */
            results.add({ this->I3d1_fp_2pt_g_p(0.04, -5.0), "f_+: I_3d1^{2pt,g_p}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d1_fp_2pt_g_p(0.04,  0.0), "f_+: I_3d1^{2pt,g_p}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d1_fp_2pt_g_p(0.04, +5.0), "f_+: I_3d1^{2pt,g_p}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3d1_fp_2pt_g_p(0.08, -5.0), "f_+: I_3d1^{2pt,g_p}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d1_fp_2pt_g_p(0.08,  0.0), "f_+: I_3d1^{2pt,g_p}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d1_fp_2pt_g_p(0.08, +5.0), "f_+: I_3d1^{2pt,g_p}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3d2 g_+ */
            results.add({ this->I3d2_fp_2pt_g_p(0.04, -5.0), "f_+: I_3d2^{2pt,g_p}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d2_fp_2pt_g_p(0.04,  0.0), "f_+: I_3d2^{2pt,g_p}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d2_fp_2pt_g_p(0.04, +5.0), "f_+: I_3d2^{2pt,g_p}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3d2_fp_2pt_g_p(0.08, -5.0), "f_+: I_3d2^{2pt,g_p}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d2_fp_2pt_g_p(0.08,  0.0), "f_+: I_3d2^{2pt,g_p}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d2_fp_2pt_g_p(0.08, +5.0), "f_+: I_3d2^{2pt,g_p}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3 g_bar */
            results.add({ this->I3_fp_2pt_g_bar(0.04, -5.0), "f_+: I_3^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3_fp_2pt_g_bar(0.04,  0.0), "f_+: I_3^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3_fp_2pt_g_bar(0.04, +5.0), "f_+: I_3^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3_fp_2pt_g_bar(0.08, -5.0), "f_+: I_3^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3_fp_2pt_g_bar(0.08,  0.0), "f_+: I_3^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3_fp_2pt_g_bar(0.08, +5.0), "f_+: I_3^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3d1 g_bar */
            results.add({ this->I3d1_fp_2pt_g_bar(0.04, -5.0), "f_+: I_3d1^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d1_fp_2pt_g_bar(0.04,  0.0), "f_+: I_3d1^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d1_fp_2pt_g_bar(0.04, +5.0), "f_+: I_3d1^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3d1_fp_2pt_g_bar(0.08, -5.0), "f_+: I_3d1^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d1_fp_2pt_g_bar(0.08,  0.0), "f_+: I_3d1^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d1_fp_2pt_g_bar(0.08, +5.0), "f_+: I_3d1^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3d2 g_bar */
            results.add({ this->I3d2_fp_2pt_g_bar(0.04, -5.0), "f_+: I_3d2^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d2_fp_2pt_g_bar(0.04,  0.0), "f_+: I_3d2^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d2_fp_2pt_g_bar(0.04, +5.0), "f_+: I_3d2^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3d2_fp_2pt_g_bar(0.08, -5.0), "f_+: I_3d2^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d2_fp_2pt_g_bar(0.08,  0.0), "f_+: I_3d2^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d2_fp_2pt_g_bar(0.08, +5.0), "f_+: I_3d2^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4 g_bar */
            results.add({ this->I4_fp_2pt_g_bar(0.04, -5.0), "f_+: I_4^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4_fp_2pt_g_bar(0.04,  0.0), "f_+: I_4^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4_fp_2pt_g_bar(0.04, +5.0), "f_+: I_4^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4_fp_2pt_g_bar(0.08, -5.0), "f_+: I_4^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4_fp_2pt_g_bar(0.08,  0.0), "f_+: I_4^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4_fp_2pt_g_bar(0.08, +5.0), "f_+: I_4^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4d1 g_bar */
            results.add({ this->I4d1_fp_2pt_g_bar(0.04, -5.0), "f_+: I_4d1^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d1_fp_2pt_g_bar(0.04,  0.0), "f_+: I_4d1^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d1_fp_2pt_g_bar(0.04, +5.0), "f_+: I_4d1^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4d1_fp_2pt_g_bar(0.08, -5.0), "f_+: I_4d1^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d1_fp_2pt_g_bar(0.08,  0.0), "f_+: I_4d1^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d1_fp_2pt_g_bar(0.08, +5.0), "f_+: I_4d1^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4d2 g_bar */
            results.add({ this->I4d2_fp_2pt_g_bar(0.04, -5.0), "f_+: I_4d2^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d2_fp_2pt_g_bar(0.04,  0.0), "f_+: I_4d2^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d2_fp_2pt_g_bar(0.04, +5.0), "f_+: I_4d2^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4d2_fp_2pt_g_bar(0.08, -5.0), "f_+: I_4d2^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d2_fp_2pt_g_bar(0.08,  0.0), "f_+: I_4d2^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d2_fp_2pt_g_bar(0.08, +5.0), "f_+: I_4d2^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4d3 g_bar */
            results.add({ this->I4d3_fp_2pt_g_bar(0.04, -5.0), "f_+: I_4d3^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d3_fp_2pt_g_bar(0.04,  0.0), "f_+: I_4d3^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d3_fp_2pt_g_bar(0.04, +5.0), "f_+: I_4d3^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4d3_fp_2pt_g_bar(0.08, -5.0), "f_+: I_4d3^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d3_fp_2pt_g_bar(0.08,  0.0), "f_+: I_4d3^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d3_fp_2pt_g_bar(0.08, +5.0), "f_+: I_4d3^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* 3 particle */

            /* I_2 phi_3 */
            results.add({ this->I2_fp_3pt_phi_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),            "f_+: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),            "f_+: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),            "f_+: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),            "f_+: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_2 phi_bar_3 */
            results.add({ this->I2_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),        "f_+: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),        "f_+: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),        "f_+: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),        "f_+: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 phi_bar_3 */
            results.add({ this->I3_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),        "f_+: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),        "f_+: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),        "f_+: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),        "f_+: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A phi_bar_3 */
            results.add({ this->I3d1A_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),     "f_+: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),     "f_+: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),     "f_+: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),     "f_+: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B phi_bar_3 */
            results.add({ this->I3d1B_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 5.0),          "f_+: I_3d1B^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 5.0),          "f_+: I_3d1B^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C phi_bar_3 */
            results.add({ this->I3d1C_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 0.1, 5.0),          "f_+: I_3d1C^{3pt,phi_bar_3}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fp_3pt_phi_bar_3(this->sigma(s0_0_p(), 5.0), 0.5, 5.0),          "f_+: I_3d1C^{3pt,phi_bar_3}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4 phi_bar_bar_3 */
            results.add({ this->I4_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),    "f_+: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),    "f_+: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),    "f_+: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),    "f_+: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1A phi_bar_bar_3 */
            results.add({ this->I4d1A_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0), "f_+: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0), "f_+: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0), "f_+: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0), "f_+: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1B phi_bar_bar_3 */
            results.add({ this->I4d1B_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 5.0),      "f_+: I_4d1B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d1B_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 5.0),      "f_+: I_4d1B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d1C phi_bar_bar_3 */
            results.add({ this->I4d1C_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 0.1, 5.0),      "f_+: I_4d1C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1C_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 0.5, 5.0),      "f_+: I_4d1C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2A phi_bar_bar_3 */
            results.add({ this->I4d2A_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0), "f_+: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0), "f_+: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0), "f_+: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0), "f_+: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2B phi_bar_bar_3 */
            results.add({ this->I4d2B_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 1.0, 5.0),      "f_+: I_4d2B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d2B_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 3.0, 5.0),      "f_+: I_4d2B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d2C phi_bar_bar_3 */
            results.add({ this->I4d2C_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 0.1, 5.0),      "f_+: I_4d2C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2C_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 0.5, 5.0),      "f_+: I_4d2C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2D phi_bar_bar_3 */
            results.add({ this->I4d2D_fp_3pt_phi_bar_bar_3(this->sigma(s0_0_p(), 5.0), 5.0),           "f_+: I_4d2D^{3pt,phi_bar_bar_3}(sigma=sigma_0, q2=5.0 GeV^2)"});

            /* I_2 phi_4 */
            results.add({ this->I2_fp_3pt_phi_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),            "f_+: I_2^{3pt,phi_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),            "f_+: I_2^{3pt,phi_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),            "f_+: I_2^{3pt,phi_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),            "f_+: I_2^{3pt,phi_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_2 phi_bar_4 */
            results.add({ this->I2_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),        "f_+: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),        "f_+: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),        "f_+: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),        "f_+: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 phi_bar_4 */
            results.add({ this->I3_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),        "f_+: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),        "f_+: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),        "f_+: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),        "f_+: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A phi_bar_4 */
            results.add({ this->I3d1A_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),     "f_+: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),     "f_+: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),     "f_+: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),     "f_+: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B phi_bar_4 */
            results.add({ this->I3d1B_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 5.0),          "f_+: I_3d1B^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 5.0),          "f_+: I_3d1B^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C phi_bar_4 */
            results.add({ this->I3d1C_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 0.1, 5.0),          "f_+: I_3d1C^{3pt,phi_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fp_3pt_phi_bar_4(this->sigma(s0_0_p(), 5.0), 0.5, 5.0),          "f_+: I_3d1C^{3pt,phi_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 phi_bar_bar_4 */
            results.add({ this->I3_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),    "f_+: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),    "f_+: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),    "f_+: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),    "f_+: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A phi_bar_bar_4 */
            results.add({ this->I3d1A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0), "f_+: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0), "f_+: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0), "f_+: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0), "f_+: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B phi_bar_bar_4 */
            results.add({ this->I3d1B_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 5.0),      "f_+: I_3d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 5.0),      "f_+: I_3d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C phi_bar_bar_4 */
            results.add({ this->I3d1C_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 0.1, 5.0),      "f_+: I_3d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 0.5, 5.0),      "f_+: I_3d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4 phi_bar_bar_4 */
            results.add({ this->I4_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),    "f_+: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),    "f_+: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),    "f_+: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),    "f_+: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1A phi_bar_bar_4 */
            results.add({ this->I4d1A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0), "f_+: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0), "f_+: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0), "f_+: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0), "f_+: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1B phi_bar_bar_4 */
            results.add({ this->I4d1B_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 5.0),      "f_+: I_4d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d1B_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 5.0),      "f_+: I_4d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d1C phi_bar_bar_4 */
            results.add({ this->I4d1C_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 0.1, 5.0),      "f_+: I_4d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1C_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 0.5, 5.0),      "f_+: I_4d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2A phi_bar_bar_4 */
            results.add({ this->I4d2A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0), "f_+: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0), "f_+: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0), "f_+: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0), "f_+: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2B phi_bar_bar_4 */
            results.add({ this->I4d2B_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 5.0),      "f_+: I_4d2B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d2B_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 5.0),      "f_+: I_4d2B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d2C phi_bar_bar_4 */
            results.add({ this->I4d2C_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 0.1, 5.0),      "f_+: I_4d2C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2C_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 0.5, 5.0),      "f_+: I_4d2C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2D phi_bar_bar_4 */
            results.add({ this->I4d2D_fp_3pt_phi_bar_bar_4(this->sigma(s0_0_p(), 5.0), 5.0),           "f_+: I_4d2D^{3pt,phi_bar_bar_4}(sigma=sigma_0, q2=5.0 GeV^2)"});

            /* I_2 psi_bar_4 */
            results.add({ this->I2_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),        "f_+: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),        "f_+: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),        "f_+: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),        "f_+: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 psi_bar_4 */
            results.add({ this->I3_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),        "f_+: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),        "f_+: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),        "f_+: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),        "f_+: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A psi_bar_4 */
            results.add({ this->I3d1A_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),     "f_+: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),     "f_+: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),     "f_+: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),     "f_+: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B psi_bar_4 */
            results.add({ this->I3d1B_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 5.0),          "f_+: I_3d1B^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 5.0),          "f_+: I_3d1B^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C psi_bar_4 */
            results.add({ this->I3d1C_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 0.1, 5.0),          "f_+: I_3d1C^{3pt,psi_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fp_3pt_psi_bar_4(this->sigma(s0_0_p(), 5.0), 0.5, 5.0),          "f_+: I_3d1C^{3pt,psi_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_2 chi_bar_4 */
            results.add({ this->I2_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),        "f_+: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),        "f_+: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),        "f_+: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),        "f_+: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 chi_bar_4 */
            results.add({ this->I3_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),        "f_+: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),        "f_+: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),        "f_+: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),        "f_+: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A chi_bar_4 */
            results.add({ this->I3d1A_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.1, 5.0),     "f_+: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 0.5, 5.0),     "f_+: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.1, 5.0),     "f_+: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 0.5, 5.0),     "f_+: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B chi_bar_4 */
            results.add({ this->I3d1B_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 1.0, 5.0),          "f_+: I_3d1B^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 3.0, 5.0),          "f_+: I_3d1B^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C chi_bar_4 */
            results.add({ this->I3d1C_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 0.1, 5.0),          "f_+: I_3d1C^{3pt,chi_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fp_3pt_chi_bar_4(this->sigma(s0_0_p(), 5.0), 0.5, 5.0),          "f_+: I_3d1C^{3pt,chi_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* f_± */

            /* 2 particle */

            /* I_1 phi_+ */
            results.add({ this->I1_fpm_2pt_phi_p(0.04, -5.0), "f_±: I_1^{2pt,phi_+}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I1_fpm_2pt_phi_p(0.04,  0.0), "f_±: I_1^{2pt,phi_+}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I1_fpm_2pt_phi_p(0.04, +5.0), "f_±: I_1^{2pt,phi_+}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I1_fpm_2pt_phi_p(0.08, -5.0), "f_±: I_1^{2pt,phi_+}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I1_fpm_2pt_phi_p(0.08,  0.0), "f_±: I_1^{2pt,phi_+}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I1_fpm_2pt_phi_p(0.08, +5.0), "f_±: I_1^{2pt,phi_+}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2 phi_bar */
            results.add({ this->I2_fpm_2pt_phi_bar(0.04, -5.0), "f_±: I_2^{2pt,phi_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fpm_2pt_phi_bar(0.04,  0.0), "f_±: I_2^{2pt,phi_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fpm_2pt_phi_bar(0.04, +5.0), "f_±: I_2^{2pt,phi_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2_fpm_2pt_phi_bar(0.08, -5.0), "f_±: I_2^{2pt,phi_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fpm_2pt_phi_bar(0.08,  0.0), "f_±: I_2^{2pt,phi_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fpm_2pt_phi_bar(0.08, +5.0), "f_±: I_2^{2pt,phi_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2d1 phi_bar */
            results.add({ this->I2d1_fpm_2pt_phi_bar(0.04, -5.0), "f_±: I_2d1^{2pt,phi_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fpm_2pt_phi_bar(0.04,  0.0), "f_±: I_2d1^{2pt,phi_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fpm_2pt_phi_bar(0.04, +5.0), "f_±: I_2d1^{2pt,phi_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2d1_fpm_2pt_phi_bar(0.08, -5.0), "f_±: I_2d1^{2pt,phi_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fpm_2pt_phi_bar(0.08,  0.0), "f_±: I_2d1^{2pt,phi_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fpm_2pt_phi_bar(0.08, +5.0), "f_±: I_2d1^{2pt,phi_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2 g_+ */
            results.add({ this->I2_fpm_2pt_g_p(0.04, -5.0), "f_±: I_2^{2pt,g_p}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fpm_2pt_g_p(0.04,  0.0), "f_±: I_2^{2pt,g_p}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fpm_2pt_g_p(0.04, +5.0), "f_±: I_2^{2pt,g_p}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2_fpm_2pt_g_p(0.08, -5.0), "f_±: I_2^{2pt,g_p}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fpm_2pt_g_p(0.08,  0.0), "f_±: I_2^{2pt,g_p}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fpm_2pt_g_p(0.08, +5.0), "f_±: I_2^{2pt,g_p}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2d1 g_+ */
            results.add({ this->I2d1_fpm_2pt_g_p(0.04, -5.0), "f_±: I_2d1^{2pt,g_p}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fpm_2pt_g_p(0.04,  0.0), "f_±: I_2d1^{2pt,g_p}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fpm_2pt_g_p(0.04, +5.0), "f_±: I_2d1^{2pt,g_p}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2d1_fpm_2pt_g_p(0.08, -5.0), "f_±: I_2d1^{2pt,g_p}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fpm_2pt_g_p(0.08,  0.0), "f_±: I_2d1^{2pt,g_p}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fpm_2pt_g_p(0.08, +5.0), "f_±: I_2d1^{2pt,g_p}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3 g_+ */
            results.add({ this->I3_fpm_2pt_g_p(0.04, -5.0), "f_±: I_3^{2pt,g_p}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3_fpm_2pt_g_p(0.04,  0.0), "f_±: I_3^{2pt,g_p}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3_fpm_2pt_g_p(0.04, +5.0), "f_±: I_3^{2pt,g_p}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3_fpm_2pt_g_p(0.08, -5.0), "f_±: I_3^{2pt,g_p}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3_fpm_2pt_g_p(0.08,  0.0), "f_±: I_3^{2pt,g_p}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3_fpm_2pt_g_p(0.08, +5.0), "f_±: I_3^{2pt,g_p}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3d1 g_+ */
            results.add({ this->I3d1_fpm_2pt_g_p(0.04, -5.0), "f_±: I_3d1^{2pt,g_p}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d1_fpm_2pt_g_p(0.04,  0.0), "f_±: I_3d1^{2pt,g_p}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d1_fpm_2pt_g_p(0.04, +5.0), "f_±: I_3d1^{2pt,g_p}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3d1_fpm_2pt_g_p(0.08, -5.0), "f_±: I_3d1^{2pt,g_p}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d1_fpm_2pt_g_p(0.08,  0.0), "f_±: I_3d1^{2pt,g_p}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d1_fpm_2pt_g_p(0.08, +5.0), "f_±: I_3d1^{2pt,g_p}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3d2 g_+ */
            results.add({ this->I3d2_fpm_2pt_g_p(0.04, -5.0), "f_±: I_3d2^{2pt,g_p}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d2_fpm_2pt_g_p(0.04,  0.0), "f_±: I_3d2^{2pt,g_p}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d2_fpm_2pt_g_p(0.04, +5.0), "f_±: I_3d2^{2pt,g_p}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3d2_fpm_2pt_g_p(0.08, -5.0), "f_±: I_3d2^{2pt,g_p}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d2_fpm_2pt_g_p(0.08,  0.0), "f_±: I_3d2^{2pt,g_p}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d2_fpm_2pt_g_p(0.08, +5.0), "f_±: I_3d2^{2pt,g_p}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3 g_bar */
            results.add({ this->I3_fpm_2pt_g_bar(0.04, -5.0), "f_±: I_3^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3_fpm_2pt_g_bar(0.04,  0.0), "f_±: I_3^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3_fpm_2pt_g_bar(0.04, +5.0), "f_±: I_3^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3_fpm_2pt_g_bar(0.08, -5.0), "f_±: I_3^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3_fpm_2pt_g_bar(0.08,  0.0), "f_±: I_3^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3_fpm_2pt_g_bar(0.08, +5.0), "f_±: I_3^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3d1 g_bar */
            results.add({ this->I3d1_fpm_2pt_g_bar(0.04, -5.0), "f_±: I_3d1^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d1_fpm_2pt_g_bar(0.04,  0.0), "f_±: I_3d1^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d1_fpm_2pt_g_bar(0.04, +5.0), "f_±: I_3d1^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3d1_fpm_2pt_g_bar(0.08, -5.0), "f_±: I_3d1^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d1_fpm_2pt_g_bar(0.08,  0.0), "f_±: I_3d1^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d1_fpm_2pt_g_bar(0.08, +5.0), "f_±: I_3d1^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3d2 g_bar */
            results.add({ this->I3d2_fpm_2pt_g_bar(0.04, -5.0), "f_±: I_3d2^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d2_fpm_2pt_g_bar(0.04,  0.0), "f_±: I_3d2^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d2_fpm_2pt_g_bar(0.04, +5.0), "f_±: I_3d2^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3d2_fpm_2pt_g_bar(0.08, -5.0), "f_±: I_3d2^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d2_fpm_2pt_g_bar(0.08,  0.0), "f_±: I_3d2^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d2_fpm_2pt_g_bar(0.08, +5.0), "f_±: I_3d2^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4 g_bar */
            results.add({ this->I4_fpm_2pt_g_bar(0.04, -5.0), "f_±: I_4^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4_fpm_2pt_g_bar(0.04,  0.0), "f_±: I_4^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4_fpm_2pt_g_bar(0.04, +5.0), "f_±: I_4^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4_fpm_2pt_g_bar(0.08, -5.0), "f_±: I_4^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4_fpm_2pt_g_bar(0.08,  0.0), "f_±: I_4^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4_fpm_2pt_g_bar(0.08, +5.0), "f_±: I_4^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4d1 g_bar */
            results.add({ this->I4d1_fpm_2pt_g_bar(0.04, -5.0), "f_±: I_4d1^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d1_fpm_2pt_g_bar(0.04,  0.0), "f_±: I_4d1^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d1_fpm_2pt_g_bar(0.04, +5.0), "f_±: I_4d1^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4d1_fpm_2pt_g_bar(0.08, -5.0), "f_±: I_4d1^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d1_fpm_2pt_g_bar(0.08,  0.0), "f_±: I_4d1^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d1_fpm_2pt_g_bar(0.08, +5.0), "f_±: I_4d1^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4d2 g_bar */
            results.add({ this->I4d2_fpm_2pt_g_bar(0.04, -5.0), "f_±: I_4d2^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d2_fpm_2pt_g_bar(0.04,  0.0), "f_±: I_4d2^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d2_fpm_2pt_g_bar(0.04, +5.0), "f_±: I_4d2^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4d2_fpm_2pt_g_bar(0.08, -5.0), "f_±: I_4d2^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d2_fpm_2pt_g_bar(0.08,  0.0), "f_±: I_4d2^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d2_fpm_2pt_g_bar(0.08, +5.0), "f_±: I_4d2^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4d3 g_bar */
            results.add({ this->I4d3_fpm_2pt_g_bar(0.04, -5.0), "f_±: I_4d3^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d3_fpm_2pt_g_bar(0.04,  0.0), "f_±: I_4d3^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d3_fpm_2pt_g_bar(0.04, +5.0), "f_±: I_4d3^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4d3_fpm_2pt_g_bar(0.08, -5.0), "f_±: I_4d3^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d3_fpm_2pt_g_bar(0.08,  0.0), "f_±: I_4d3^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d3_fpm_2pt_g_bar(0.08, +5.0), "f_±: I_4d3^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* 3 particle */

            /* I_2 phi_3 */
            results.add({ this->I2_fpm_3pt_phi_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),            "f_±: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),            "f_±: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),            "f_±: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),            "f_±: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_2 phi_bar_3 */
            results.add({ this->I2_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),        "f_±: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),        "f_±: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),        "f_±: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),        "f_±: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 phi_bar_3 */
            results.add({ this->I3_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),        "f_±: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),        "f_±: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),        "f_±: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),        "f_±: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A phi_bar_3 */
            results.add({ this->I3d1A_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),     "f_±: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),     "f_±: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),     "f_±: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),     "f_±: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B phi_bar_3 */
            results.add({ this->I3d1B_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 5.0),          "f_±: I_3d1B^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 5.0),          "f_±: I_3d1B^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C phi_bar_3 */
            results.add({ this->I3d1C_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 0.1, 5.0),          "f_±: I_3d1C^{3pt,phi_bar_3}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fpm_3pt_phi_bar_3(this->sigma(s0_0_pm(), 5.0), 0.5, 5.0),          "f_±: I_3d1C^{3pt,phi_bar_3}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4 phi_bar_bar_3 */
            results.add({ this->I4_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),    "f_±: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),    "f_±: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),    "f_±: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),    "f_±: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1A phi_bar_bar_3 */
            results.add({ this->I4d1A_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0), "f_±: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0), "f_±: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0), "f_±: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0), "f_±: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1B phi_bar_bar_3 */
            results.add({ this->I4d1B_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 5.0),      "f_±: I_4d1B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d1B_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 5.0),      "f_±: I_4d1B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d1C phi_bar_bar_3 */
            results.add({ this->I4d1C_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 0.1, 5.0),      "f_±: I_4d1C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1C_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 0.5, 5.0),      "f_±: I_4d1C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2A phi_bar_bar_3 */
            results.add({ this->I4d2A_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0), "f_±: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0), "f_±: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0), "f_±: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0), "f_±: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2B phi_bar_bar_3 */
            results.add({ this->I4d2B_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 1.0, 5.0),      "f_±: I_4d2B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d2B_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 3.0, 5.0),      "f_±: I_4d2B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d2C phi_bar_bar_3 */
            results.add({ this->I4d2C_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 0.1, 5.0),      "f_±: I_4d2C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2C_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 0.5, 5.0),      "f_±: I_4d2C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2D phi_bar_bar_3 */
            results.add({ this->I4d2D_fpm_3pt_phi_bar_bar_3(this->sigma(s0_0_pm(), 5.0), 5.0),           "f_±: I_4d2D^{3pt,phi_bar_bar_3}(sigma=sigma_0, q2=5.0 GeV^2)"});

            /* I_2 phi_4 */
            results.add({ this->I2_fpm_3pt_phi_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),            "f_±: I_2^{3pt,phi_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),            "f_±: I_2^{3pt,phi_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),            "f_±: I_2^{3pt,phi_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),            "f_±: I_2^{3pt,phi_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_2 phi_bar_4 */
            results.add({ this->I2_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),        "f_±: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),        "f_±: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),        "f_±: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),        "f_±: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 phi_bar_4 */
            results.add({ this->I3_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),        "f_±: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),        "f_±: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),        "f_±: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),        "f_±: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A phi_bar_4 */
            results.add({ this->I3d1A_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),     "f_±: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),     "f_±: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),     "f_±: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),     "f_±: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B phi_bar_4 */
            results.add({ this->I3d1B_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 5.0),          "f_±: I_3d1B^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 5.0),          "f_±: I_3d1B^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C phi_bar_4 */
            results.add({ this->I3d1C_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 0.1, 5.0),          "f_±: I_3d1C^{3pt,phi_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fpm_3pt_phi_bar_4(this->sigma(s0_0_pm(), 5.0), 0.5, 5.0),          "f_±: I_3d1C^{3pt,phi_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 phi_bar_bar_4 */
            results.add({ this->I3_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),    "f_±: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),    "f_±: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),    "f_±: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),    "f_±: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A phi_bar_bar_4 */
            results.add({ this->I3d1A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0), "f_±: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0), "f_±: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0), "f_±: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0), "f_±: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B phi_bar_bar_4 */
            results.add({ this->I3d1B_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 5.0),      "f_±: I_3d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 5.0),      "f_±: I_3d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C phi_bar_bar_4 */
            results.add({ this->I3d1C_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 0.1, 5.0),      "f_±: I_3d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 0.5, 5.0),      "f_±: I_3d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4 phi_bar_bar_4 */
            results.add({ this->I4_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),    "f_±: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),    "f_±: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),    "f_±: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),    "f_±: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1A phi_bar_bar_4 */
            results.add({ this->I4d1A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0), "f_±: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0), "f_±: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0), "f_±: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0), "f_±: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1B phi_bar_bar_4 */
            results.add({ this->I4d1B_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 5.0),      "f_±: I_4d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d1B_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 5.0),      "f_±: I_4d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d1C phi_bar_bar_4 */
            results.add({ this->I4d1C_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 0.1, 5.0),      "f_±: I_4d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1C_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 0.5, 5.0),      "f_±: I_4d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2A phi_bar_bar_4 */
            results.add({ this->I4d2A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0), "f_±: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0), "f_±: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0), "f_±: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0), "f_±: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2B phi_bar_bar_4 */
            results.add({ this->I4d2B_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 5.0),      "f_±: I_4d2B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d2B_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 5.0),      "f_±: I_4d2B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d2C phi_bar_bar_4 */
            results.add({ this->I4d2C_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 0.1, 5.0),      "f_±: I_4d2C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2C_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 0.5, 5.0),      "f_±: I_4d2C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2D phi_bar_bar_4 */
            results.add({ this->I4d2D_fpm_3pt_phi_bar_bar_4(this->sigma(s0_0_pm(), 5.0), 5.0),           "f_±: I_4d2D^{3pt,phi_bar_bar_4}(sigma=sigma_0, q2=5.0 GeV^2)"});

            /* I_2 psi_bar_4 */
            results.add({ this->I2_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),        "f_±: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),        "f_±: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),        "f_±: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),        "f_±: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 psi_bar_4 */
            results.add({ this->I3_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),        "f_±: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),        "f_±: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),        "f_±: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),        "f_±: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A psi_bar_4 */
            results.add({ this->I3d1A_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),     "f_±: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),     "f_±: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),     "f_±: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),     "f_±: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B psi_bar_4 */
            results.add({ this->I3d1B_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 5.0),          "f_±: I_3d1B^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 5.0),          "f_±: I_3d1B^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C psi_bar_4 */
            results.add({ this->I3d1C_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 0.1, 5.0),          "f_±: I_3d1C^{3pt,psi_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fpm_3pt_psi_bar_4(this->sigma(s0_0_pm(), 5.0), 0.5, 5.0),          "f_±: I_3d1C^{3pt,psi_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_2 chi_bar_4 */
            results.add({ this->I2_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),        "f_±: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),        "f_±: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),        "f_±: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),        "f_±: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 chi_bar_4 */
            results.add({ this->I3_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),        "f_±: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),        "f_±: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),        "f_±: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),        "f_±: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A chi_bar_4 */
            results.add({ this->I3d1A_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.1, 5.0),     "f_±: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 0.5, 5.0),     "f_±: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.1, 5.0),     "f_±: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 0.5, 5.0),     "f_±: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B chi_bar_4 */
            results.add({ this->I3d1B_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 1.0, 5.0),          "f_±: I_3d1B^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 3.0, 5.0),          "f_±: I_3d1B^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C chi_bar_4 */
            results.add({ this->I3d1C_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 0.1, 5.0),          "f_±: I_3d1C^{3pt,chi_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fpm_3pt_chi_bar_4(this->sigma(s0_0_pm(), 5.0), 0.5, 5.0),          "f_±: I_3d1C^{3pt,chi_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* f_T */

            /* 2 particle */

            /* I_1 phi_bar */
            results.add({ this->I1_fT_2pt_phi_bar(0.04, -5.0), "f_T: I_1^{2pt,phi_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I1_fT_2pt_phi_bar(0.04,  0.0), "f_T: I_1^{2pt,phi_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I1_fT_2pt_phi_bar(0.04, +5.0), "f_T: I_1^{2pt,phi_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I1_fT_2pt_phi_bar(0.08, -5.0), "f_T: I_1^{2pt,phi_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I1_fT_2pt_phi_bar(0.08,  0.0), "f_T: I_1^{2pt,phi_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I1_fT_2pt_phi_bar(0.08, +5.0), "f_T: I_1^{2pt,phi_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2 phi_bar */
            results.add({ this->I2_fT_2pt_phi_bar(0.04, -5.0), "f_T: I_2^{2pt,phi_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fT_2pt_phi_bar(0.04,  0.0), "f_T: I_2^{2pt,phi_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fT_2pt_phi_bar(0.04, +5.0), "f_T: I_2^{2pt,phi_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2_fT_2pt_phi_bar(0.08, -5.0), "f_T: I_2^{2pt,phi_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fT_2pt_phi_bar(0.08,  0.0), "f_T: I_2^{2pt,phi_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fT_2pt_phi_bar(0.08, +5.0), "f_T: I_2^{2pt,phi_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2d1 phi_bar */
            results.add({ this->I2d1_fT_2pt_phi_bar(0.04, -5.0), "f_T: I_2d1^{2pt,phi_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fT_2pt_phi_bar(0.04,  0.0), "f_T: I_2d1^{2pt,phi_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fT_2pt_phi_bar(0.04, +5.0), "f_T: I_2d1^{2pt,phi_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2d1_fT_2pt_phi_bar(0.08, -5.0), "f_T: I_2d1^{2pt,phi_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fT_2pt_phi_bar(0.08,  0.0), "f_T: I_2d1^{2pt,phi_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fT_2pt_phi_bar(0.08, +5.0), "f_T: I_2d1^{2pt,phi_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2 g_bar */
            results.add({ this->I2_fT_2pt_g_bar(0.04, -5.0), "f_T: I_2^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fT_2pt_g_bar(0.04,  0.0), "f_T: I_2^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fT_2pt_g_bar(0.04, +5.0), "f_T: I_2^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2_fT_2pt_g_bar(0.08, -5.0), "f_T: I_2^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2_fT_2pt_g_bar(0.08,  0.0), "f_T: I_2^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2_fT_2pt_g_bar(0.08, +5.0), "f_T: I_2^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_2d1 g_bar */
            results.add({ this->I2d1_fT_2pt_g_bar(0.04, -5.0), "f_T: I_2d1^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fT_2pt_g_bar(0.04,  0.0), "f_T: I_2d1^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fT_2pt_g_bar(0.04, +5.0), "f_T: I_2d1^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I2d1_fT_2pt_g_bar(0.08, -5.0), "f_T: I_2d1^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I2d1_fT_2pt_g_bar(0.08,  0.0), "f_T: I_2d1^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I2d1_fT_2pt_g_bar(0.08, +5.0), "f_T: I_2d1^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3 g_bar */
            results.add({ this->I3_fT_2pt_g_bar(0.04, -5.0), "f_T: I_3^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3_fT_2pt_g_bar(0.04,  0.0), "f_T: I_3^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3_fT_2pt_g_bar(0.04, +5.0), "f_T: I_3^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3_fT_2pt_g_bar(0.08, -5.0), "f_T: I_3^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3_fT_2pt_g_bar(0.08,  0.0), "f_T: I_3^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3_fT_2pt_g_bar(0.08, +5.0), "f_T: I_3^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3d1 g_bar */
            results.add({ this->I3d1_fT_2pt_g_bar(0.04, -5.0), "f_T: I_3d1^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d1_fT_2pt_g_bar(0.04,  0.0), "f_T: I_3d1^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d1_fT_2pt_g_bar(0.04, +5.0), "f_T: I_3d1^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3d1_fT_2pt_g_bar(0.08, -5.0), "f_T: I_3d1^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d1_fT_2pt_g_bar(0.08,  0.0), "f_T: I_3d1^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d1_fT_2pt_g_bar(0.08, +5.0), "f_T: I_3d1^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_3d2 g_bar */
            results.add({ this->I3d2_fT_2pt_g_bar(0.04, -5.0), "f_T: I_3d2^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d2_fT_2pt_g_bar(0.04,  0.0), "f_T: I_3d2^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d2_fT_2pt_g_bar(0.04, +5.0), "f_T: I_3d2^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I3d2_fT_2pt_g_bar(0.08, -5.0), "f_T: I_3d2^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I3d2_fT_2pt_g_bar(0.08,  0.0), "f_T: I_3d2^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I3d2_fT_2pt_g_bar(0.08, +5.0), "f_T: I_3d2^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4 g_bar */
            results.add({ this->I4_fT_2pt_g_bar(0.04, -5.0), "f_T: I_4^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4_fT_2pt_g_bar(0.04,  0.0), "f_T: I_4^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4_fT_2pt_g_bar(0.04, +5.0), "f_T: I_4^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4_fT_2pt_g_bar(0.08, -5.0), "f_T: I_4^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4_fT_2pt_g_bar(0.08,  0.0), "f_T: I_4^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4_fT_2pt_g_bar(0.08, +5.0), "f_T: I_4^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4d1 g_bar */
            results.add({ this->I4d1_fT_2pt_g_bar(0.04, -5.0), "f_T: I_4d1^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d1_fT_2pt_g_bar(0.04,  0.0), "f_T: I_4d1^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d1_fT_2pt_g_bar(0.04, +5.0), "f_T: I_4d1^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4d1_fT_2pt_g_bar(0.08, -5.0), "f_T: I_4d1^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d1_fT_2pt_g_bar(0.08,  0.0), "f_T: I_4d1^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d1_fT_2pt_g_bar(0.08, +5.0), "f_T: I_4d1^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4d2 g_bar */
            results.add({ this->I4d2_fT_2pt_g_bar(0.04, -5.0), "f_T: I_4d2^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d2_fT_2pt_g_bar(0.04,  0.0), "f_T: I_4d2^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d2_fT_2pt_g_bar(0.04, +5.0), "f_T: I_4d2^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4d2_fT_2pt_g_bar(0.08, -5.0), "f_T: I_4d2^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d2_fT_2pt_g_bar(0.08,  0.0), "f_T: I_4d2^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d2_fT_2pt_g_bar(0.08, +5.0), "f_T: I_4d2^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* I_4d3 g_bar */
            results.add({ this->I4d3_fT_2pt_g_bar(0.04, -5.0), "f_T: I_4d3^{2pt,g_bar}(sigma = 0.04, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d3_fT_2pt_g_bar(0.04,  0.0), "f_T: I_4d3^{2pt,g_bar}(sigma = 0.04, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d3_fT_2pt_g_bar(0.04, +5.0), "f_T: I_4d3^{2pt,g_bar}(sigma = 0.04, q2 = +5.0 GeV^2)" });

            results.add({ this->I4d3_fT_2pt_g_bar(0.08, -5.0), "f_T: I_4d3^{2pt,g_bar}(sigma = 0.08, q2 = -5.0 GeV^2)" });
            results.add({ this->I4d3_fT_2pt_g_bar(0.08,  0.0), "f_T: I_4d3^{2pt,g_bar}(sigma = 0.08, q2 =  0.0 GeV^2)" });
            results.add({ this->I4d3_fT_2pt_g_bar(0.08, +5.0), "f_T: I_4d3^{2pt,g_bar}(sigma = 0.08, q2 = +5.0 GeV^2)" });

            /* 3 particle */

            /* I_1 phi_3 */
            results.add({ this->I1_fT_3pt_phi_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),            "f_T: I_1^{3pt,phi_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I1_fT_3pt_phi_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),            "f_T: I_1^{3pt,phi_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I1_fT_3pt_phi_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),            "f_T: I_1^{3pt,phi_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I1_fT_3pt_phi_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),            "f_T: I_1^{3pt,phi_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_2 phi_3 */
            results.add({ this->I2_fT_3pt_phi_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),            "f_T: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),            "f_T: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),            "f_T: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),            "f_T: I_2^{3pt,phi_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_2 phi_bar_3 */
            results.add({ this->I2_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),        "f_T: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),        "f_T: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),        "f_T: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),        "f_T: I_2^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 phi_bar_3 */
            results.add({ this->I3_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),        "f_T: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),        "f_T: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),        "f_T: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),        "f_T: I_3^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A phi_bar_3 */
            results.add({ this->I3d1A_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),     "f_T: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),     "f_T: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),     "f_T: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),     "f_T: I_3d1A^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B phi_bar_3 */
            results.add({ this->I3d1B_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 5.0),          "f_T: I_3d1B^{3pt,phi_bar_3}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 5.0),          "f_T: I_3d1B^{3pt,phi_bar_3}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C phi_bar_3 */
            results.add({ this->I3d1C_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 0.1, 5.0),          "f_T: I_3d1C^{3pt,phi_bar_3}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fT_3pt_phi_bar_3(this->sigma(s0_0_t(), 5.0), 0.5, 5.0),          "f_T: I_3d1C^{3pt,phi_bar_3}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 phi_bar_bar_3 */
            results.add({ this->I3_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),    "f_T: I_3^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),    "f_T: I_3^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),    "f_T: I_3^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),    "f_T: I_3^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A phi_bar_bar_3 */
            results.add({ this->I3d1A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0), "f_T: I_3d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0), "f_T: I_3d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0), "f_T: I_3d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0), "f_T: I_3d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B phi_bar_bar_3 */
            results.add({ this->I3d1B_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 5.0),      "f_T: I_3d1B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 5.0),      "f_T: I_3d1B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C phi_bar_bar_3 */
            results.add({ this->I3d1C_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 0.1, 5.0),      "f_T: I_3d1C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 0.5, 5.0),      "f_T: I_3d1C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4 phi_bar_bar_3 */
            results.add({ this->I4_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),    "f_T: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),    "f_T: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),    "f_T: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),    "f_T: I_4^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1A phi_bar_bar_3 */
            results.add({ this->I4d1A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0), "f_T: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0), "f_T: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0), "f_T: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0), "f_T: I_4d1A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1B phi_bar_bar_3 */
            results.add({ this->I4d1B_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 5.0),      "f_T: I_4d1B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d1B_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 5.0),      "f_T: I_4d1B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d1C phi_bar_bar_3 */
            results.add({ this->I4d1C_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 0.1, 5.0),      "f_T: I_4d1C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1C_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 0.5, 5.0),      "f_T: I_4d1C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2A phi_bar_bar_3 */
            results.add({ this->I4d2A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0), "f_T: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0), "f_T: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0), "f_T: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0), "f_T: I_4d2A^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2B phi_bar_bar_3 */
            results.add({ this->I4d2B_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 1.0, 5.0),      "f_T: I_4d2B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d2B_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 3.0, 5.0),      "f_T: I_4d2B^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d2C phi_bar_bar_3 */
            results.add({ this->I4d2C_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 0.1, 5.0),      "f_T: I_4d2C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2C_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 0.5, 5.0),      "f_T: I_4d2C^{3pt,phi_bar_bar_3}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2D phi_bar_bar_3 */
            results.add({ this->I4d2D_fT_3pt_phi_bar_bar_3(this->sigma(s0_0_t(), 5.0), 5.0),           "f_T: I_4d2D^{3pt,phi_bar_bar_3}(sigma=sigma_0, q2=5.0 GeV^2)"});

            /* I_2 phi_bar_4 */
            results.add({ this->I2_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),        "f_T: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),        "f_T: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),        "f_T: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),        "f_T: I_2^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 phi_bar_4 */
            results.add({ this->I3_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),        "f_T: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),        "f_T: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),        "f_T: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),        "f_T: I_3^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A phi_bar_4 */
            results.add({ this->I3d1A_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),     "f_T: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),     "f_T: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),     "f_T: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),     "f_T: I_3d1A^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B phi_bar_4 */
            results.add({ this->I3d1B_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 5.0),          "f_T: I_3d1B^{3pt,phi_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 5.0),          "f_T: I_3d1B^{3pt,phi_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C phi_bar_4 */
            results.add({ this->I3d1C_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 0.1, 5.0),          "f_T: I_3d1C^{3pt,phi_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fT_3pt_phi_bar_4(this->sigma(s0_0_t(), 5.0), 0.5, 5.0),          "f_T: I_3d1C^{3pt,phi_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_2 phi_bar_bar_4 */
            results.add({ this->I2_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),    "f_T: I_2^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),    "f_T: I_2^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),    "f_T: I_2^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),    "f_T: I_2^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 phi_bar_bar_4 */
            results.add({ this->I3_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),    "f_T: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),    "f_T: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),    "f_T: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),    "f_T: I_3^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A phi_bar_bar_4 */
            results.add({ this->I3d1A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0), "f_T: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0), "f_T: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0), "f_T: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0), "f_T: I_3d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B phi_bar_bar_4 */
            results.add({ this->I3d1B_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 5.0),      "f_T: I_3d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 5.0),      "f_T: I_3d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C phi_bar_bar_4 */
            results.add({ this->I3d1C_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 0.1, 5.0),      "f_T: I_3d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 0.5, 5.0),      "f_T: I_3d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4 phi_bar_bar_4 */
            results.add({ this->I4_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),    "f_T: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),    "f_T: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),    "f_T: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),    "f_T: I_4^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1A phi_bar_bar_4 */
            results.add({ this->I4d1A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0), "f_T: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0), "f_T: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0), "f_T: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0), "f_T: I_4d1A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d1B phi_bar_bar_4 */
            results.add({ this->I4d1B_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 5.0),      "f_T: I_4d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d1B_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 5.0),      "f_T: I_4d1B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d1C phi_bar_bar_4 */
            results.add({ this->I4d1C_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 0.1, 5.0),      "f_T: I_4d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d1C_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 0.5, 5.0),      "f_T: I_4d1C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2A phi_bar_bar_4 */
            results.add({ this->I4d2A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0), "f_T: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0), "f_T: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0), "f_T: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2A_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0), "f_T: I_4d2A^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2B phi_bar_bar_4 */
            results.add({ this->I4d2B_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 5.0),      "f_T: I_4d2B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I4d2B_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 5.0),      "f_T: I_4d2B^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_4d2C phi_bar_bar_4 */
            results.add({ this->I4d2C_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 0.1, 5.0),      "f_T: I_4d2C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I4d2C_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 0.5, 5.0),      "f_T: I_4d2C^{3pt,phi_bar_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_4d2D phi_bar_bar_4 */
            results.add({ this->I4d2D_fT_3pt_phi_bar_bar_4(this->sigma(s0_0_t(), 5.0), 5.0),           "f_T: I_4d2D^{3pt,phi_bar_bar_4}(sigma=sigma_0, q2=5.0 GeV^2)"});

            /* I_2 psi_bar_4 */
            results.add({ this->I2_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),        "f_T: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),        "f_T: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),        "f_T: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),        "f_T: I_2^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 psi_bar_4 */
            results.add({ this->I3_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),        "f_T: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),        "f_T: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),        "f_T: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),        "f_T: I_3^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A psi_bar_4 */
            results.add({ this->I3d1A_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),     "f_T: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),     "f_T: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),     "f_T: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),     "f_T: I_3d1A^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B psi_bar_4 */
            results.add({ this->I3d1B_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 5.0),          "f_T: I_3d1B^{3pt,psi_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 5.0),          "f_T: I_3d1B^{3pt,psi_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C psi_bar_4 */
            results.add({ this->I3d1C_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 0.1, 5.0),          "f_T: I_3d1C^{3pt,psi_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fT_3pt_psi_bar_4(this->sigma(s0_0_t(), 5.0), 0.5, 5.0),          "f_T: I_3d1C^{3pt,psi_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_2 chi_bar_4 */
            results.add({ this->I2_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),        "f_T: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),        "f_T: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),        "f_T: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I2_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),        "f_T: I_2^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3 chi_bar_4 */
            results.add({ this->I3_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),        "f_T: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),        "f_T: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),        "f_T: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),        "f_T: I_3^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1A chi_bar_4 */
            results.add({ this->I3d1A_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.1, 5.0),     "f_T: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 0.5, 5.0),     "f_T: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, w_2=0.5, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.1, 5.0),     "f_T: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1A_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 0.5, 5.0),     "f_T: I_3d1A^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, w_2=0.5, q2=5.0 GeV^2)"});

            /* I_3d1B chi_bar_4 */
            results.add({ this->I3d1B_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 1.0, 5.0),          "f_T: I_3d1B^{3pt,chi_bar_4}(sigma=sigma_0, w_1=1.0, q2=5.0 GeV^2)"});
            results.add({ this->I3d1B_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 3.0, 5.0),          "f_T: I_3d1B^{3pt,chi_bar_4}(sigma=sigma_0, w_1=3.0, q2=5.0 GeV^2)"});

            /* I_3d1C chi_bar_4 */
            results.add({ this->I3d1C_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 0.1, 5.0),          "f_T: I_3d1C^{3pt,chi_bar_4}(sigma=sigma_0, w_2=0.1, q2=5.0 GeV^2)"});
            results.add({ this->I3d1C_fT_3pt_chi_bar_4(this->sigma(s0_0_t(), 5.0), 0.5, 5.0),          "f_T: I_3d1C^{3pt,chi_bar_4}(sigma=sigma_0, w_2=0.5, q2=5.0 GeV^2)"});

            return results;
        }

    };

    template <typename Transition_>
    const std::vector<OptionSpecification>
    Implementation<AnalyticFormFactorBToPLCSR<Transition_>>::options
    {
        { "2pt"_ok,    { "tw2+3"s, "all"s, "off"s }, "all"s   },
        { "3pt"_ok,    { "tw3+4"s, "all"s, "off"s }, "all"s   },
        { "method"_ok, { "borel"s, "dispersive"s  }, "borel"s }
    };

    template <typename Transition_>
    AnalyticFormFactorBToPLCSR<Transition_>::AnalyticFormFactorBToPLCSR(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<AnalyticFormFactorBToPLCSR<Transition_>>(new Implementation<AnalyticFormFactorBToPLCSR<Transition_>>(p, o, *this))
    {
    }

    template <typename Transition_>
    AnalyticFormFactorBToPLCSR<Transition_>::~AnalyticFormFactorBToPLCSR()
    {
    }

    template <typename Transition_>
    FormFactors<PToP> *
    AnalyticFormFactorBToPLCSR<Transition_>::make(const Parameters & p, const Options & o)
    {
        return new AnalyticFormFactorBToPLCSR<Transition_>(p, o);
    }

    template <typename Transition_>
    double
    AnalyticFormFactorBToPLCSR<Transition_>::f_p(const double & q2) const
    {
        return this->_imp->f_p(q2);
    }

    template <typename Transition_>
    double
    AnalyticFormFactorBToPLCSR<Transition_>::f_0(const double & q2) const
    {
        const double m_B = this->_imp->m_B(), m_B2 = power_of<2>(m_B);
        const double m_P = this->_imp->m_P(), m_P2 = power_of<2>(m_P);

        return (this->_imp->f_pm(q2)-this->_imp->f_p(q2)) * q2 / (m_B2 - m_P2) + this->_imp->f_p(q2);
    }

    template <typename Transition_>
    double
    AnalyticFormFactorBToPLCSR<Transition_>::f_m(const double & q2) const
    {
        return this->_imp->f_pm(q2)-this->_imp->f_p(q2);
    }

    template <typename Transition_>
    double
    AnalyticFormFactorBToPLCSR<Transition_>::f_t(const double & q2) const
    {
        return this->_imp->f_t(q2);
    }

    template <typename Transition_>
    double
    AnalyticFormFactorBToPLCSR<Transition_>::f_plus_T(const double & q2) const
    {
        // Conventions of GvDV:2020A eq. (A.5)
        return this->_imp->f_t(q2) * q2 / this->_imp->m_B() / (this->_imp->m_B() + this->_imp->m_P());
    }

    template <typename Transition_>
    double
    AnalyticFormFactorBToPLCSR<Transition_>::normalized_moment_1_f_p(const double & q2) const
    {
        return this->_imp->normalized_moment_1_f_p(q2);
    }

    template <typename Transition_>
    double
    AnalyticFormFactorBToPLCSR<Transition_>::normalized_moment_1_f_pm(const double & q2) const
    {
        return this->_imp->normalized_moment_1_f_pm(q2);
    }

    template <typename Transition_>
    double
    AnalyticFormFactorBToPLCSR<Transition_>::normalized_moment_1_f_t(const double & q2) const
    {
        return this->_imp->normalized_moment_1_f_t(q2);
    }

    template <typename Transition_>
    Diagnostics
    AnalyticFormFactorBToPLCSR<Transition_>::diagnostics() const
    {
        return this->_imp->diagnostics();
    }

    template <typename Transition_>
    const std::set<ReferenceName>
    AnalyticFormFactorBToPLCSR<Transition_>::references
    {
        "KMO:2005A"_rn,
        "KMO:2006A"_rn,
        "FKKM:2008A"_rn,
        "GKvD:2018A"_rn
    };

    template <typename Transition_>
    std::vector<OptionSpecification>::const_iterator
    AnalyticFormFactorBToPLCSR<Transition_>::begin_options()
    {
        return Implementation<AnalyticFormFactorBToPLCSR<Transition_>>::options.cbegin();
    }

    template <typename Transition_>
    std::vector<OptionSpecification>::const_iterator
    AnalyticFormFactorBToPLCSR<Transition_>::end_options()
    {
        return Implementation<AnalyticFormFactorBToPLCSR<Transition_>>::options.cend();
    }
}

#endif
