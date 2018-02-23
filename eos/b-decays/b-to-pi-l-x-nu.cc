/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
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

#include <eos/form-factors/form-factors.hh>
#include <eos/b-decays/b-to-pi-l-x-nu.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/save.hh>

namespace eos
{
    template <>
    struct Implementation<BToPiLeptonInclusiveNeutrinos>
    {
        std::shared_ptr<FormFactors<PToP>> form_factors;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_pi;

        UsedParameter m_mu;

        UsedParameter m_tau;

        UsedParameter g_fermi;

        UsedParameter hbar;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            m_B(p["mass::B_" + o.get("q", "d")], u),
            tau_B(p["life_time::B_" + o.get("q", "d")], u),
            m_pi(p["mass::pi^" + std::string(o.get("q", "d") == "d" ? "+" : "0")], u),
            m_mu(p["mass::mu"], u),
            m_tau(p["mass::tau"], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u)
        {
            if ((o.get("q", "d") != "d") && (o.get("q", "d") != "u")) // q = d is the default
            {
                // only B_{d,u} mesons can decay in this channel
                throw InternalError("BToPiLeptonInclusiveNeutrinos: q = '" + o["q"] + "' is not a valid option for this decay channel");
            }

            form_factors = FormFactorFactory<PToP>::create("B->pi@" + o.get("form-factors", "BCL2008"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
        }

        // normalized to N_1 = |V_vb|^2 G_F^2 / (192 pi^3 MB^3)
        double differential_decay_width_1nu(const double & s) const
        {
            double fp = form_factors->f_p(s);
            double lam = lambda(m_B * m_B, m_pi * m_pi, s);

            return sqrt(lam) * (3.0 * fp * fp * lam);
        }

        double normalized_differential_decay_width_1nu(const double & s, const double & c_theta_mu) const
        {
            std::function<double (const double &)> integrand = std::bind(&Implementation<BToPiLeptonInclusiveNeutrinos>::differential_decay_width_1nu, this, std::placeholders::_1);
            const double s_min = 0.0, s_max = power_of<2>(m_B() - m_pi());
            double Gamma_1 = integrate<GSL::QNG>(integrand, s_min, s_max);

            double fp = form_factors->f_p(s);
            double lam = lambda(m_B * m_B, m_pi * m_pi, s);

            return 3.0 / 4.0 * fp * fp * lam * sqrt(lam) * (1.0 - c_theta_mu * c_theta_mu) / Gamma_1;
        }

        // normalized to N_3 = |V_ub|^2 G_F^2 / (384 pi^3 MB^3)
        //                   * tau_tau / hbar * G_F^2 m_tau^5 / (192 pi^3)
        double differential_decay_width_3nu(const double & s) const
        {
            double fp = form_factors->f_p(s);
            double f0 = form_factors->f_0(s);
            double lam = lambda(m_B * m_B, m_pi * m_pi, s);
            // make sure to return NaN if s < m_tau^2
            double sqrtv = sqrt(1.0 - m_tau() * m_tau() / s);
            double v = sqrtv * sqrtv, v2 = v * v;

            return sqrt(lam) * v2 * ((3.0 - v) * fp * fp * lam + 3.0 * (1.0 - v) * f0 * f0 * power_of<2>(m_B() * m_B() - m_pi() * m_pi())) * 4.0 / 3.0;
        }

        double normalized_differential_decay_width_3nu(const double & s, const double & snunubar,
                const double & z, const double & phi, const double & zst) const
        {
            std::function<double (const double &)> integrand = std::bind(&Implementation<BToPiLeptonInclusiveNeutrinos>::differential_decay_width_3nu, this, std::placeholders::_1);
            const double s_min = 3.16, s_max = power_of<2>(m_B() - m_pi());
            const double Gamma_3 = integrate<GSL::QNG>(integrand, s_min, s_max);

            const double fp = form_factors->f_p(s), fp2 = fp * fp;
            const double f0 = form_factors->f_0(s), f02 = f0 * f0;
            const double lam = lambda(m_B * m_B, m_pi * m_pi, s), sqrtlam = sqrt(lam);
            const double mtau = m_tau(), mtau2 = mtau * mtau, mtau4 = mtau2 * mtau2;
            const double mtau6 = mtau4 * mtau2, mtau8 = mtau4 * mtau4;
            const double mB2 = m_B() * m_B(), mD2 = m_pi() * m_pi();
            const double s2 = s * s, s3 = s2 * s, sqrts = sqrt(s);
            const double z2 = z * z;

            // constant in z
            const double a = power_of<2>((mtau2 - s) * (mtau2 - snunubar)) * sqrtlam / (mtau8 * M_PI * s3) * (
                        +(mtau2 + 2.0 * snunubar) * (f02 * power_of<2>(mB2 - mD2) * mtau2 + fp2 * s * lam)
                        -(mtau2 - 2.0 * snunubar) * (f02 * power_of<2>(mB2 - mD2) * mtau2 - fp2 * s * lam) * zst
                    );

            // multiplying z
            const double b = 2.0 * power_of<2>((mtau2 - s) * (mtau2 - snunubar)) * lam * f0 * fp * (mB2 - mD2) / (mtau6 * M_PI * s3) * (
                        +(mtau2 + 2.0 * snunubar)
                        -(mtau2 - 2.0 * snunubar) * zst
                    );

            // multiplying z^2
            const double c = power_of<2>((mtau2 - s) * (mtau2 - snunubar)) * lam * sqrtlam * fp2 / (mtau8 * M_PI * s3) * (
                        +(mtau2 + 2.0 * snunubar) * (mtau2 - s)
                        -(mtau2 - 2.0 * snunubar) * (mtau2 + s) * zst
                    );

            // multiplying sqrt(1 - z^2)
            const double d = 2.0 * mtau * sqrts * power_of<2>((mtau2 - snunubar) * (s - mtau2)) * (mtau2 - 2.0 * snunubar) * f0 * fp * (mB2 - mD2) * lam / (mtau8 * M_PI * s3)
                * sqrt(1.0 - zst * zst);

            // multiplying z sqrt(1 - z^2)
            const double e = 2.0 * mtau * sqrts * power_of<2>((mtau2 - snunubar) * (s - mtau2)) * (mtau2 - 2.0 * snunubar) * fp2 * sqrtlam * lam / (mtau8 * M_PI * s3)
                * sqrt(1.0 - zst * zst);

            return (a + b * z + c * z2 + (d + e * z) * sqrt(1.0 - z2) * cos(phi)) / Gamma_3;
        }

    };

    BToPiLeptonInclusiveNeutrinos::BToPiLeptonInclusiveNeutrinos(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToPiLeptonInclusiveNeutrinos>(new Implementation<BToPiLeptonInclusiveNeutrinos>(parameters, options, *this))
    {
    }

    BToPiLeptonInclusiveNeutrinos::~BToPiLeptonInclusiveNeutrinos()
    {
    }

    double
    BToPiLeptonInclusiveNeutrinos::normalized_differential_decay_width_1nu(const double & s, const double & c_theta_mu) const
    {
        return _imp->normalized_differential_decay_width_1nu(s, c_theta_mu);
    }

    double
    BToPiLeptonInclusiveNeutrinos::normalized_differential_decay_width_3nu(const double & s, const double & snunubar,
            const double & z, const double & phi, const double & zst) const
    {
        return _imp->normalized_differential_decay_width_3nu(s, snunubar, z, phi, zst);
    }

    const std::string
    BToPiLeptonInclusiveNeutrinos::description = "\
The neutrino-inclusive decay B->pi l X_nu, where l=e,mu is a light lepton, and \
X_nu = {nu, nu nubar nu} is either a one or three neutrino final state. The \
numerical implementations for this process follow [BIvD2016].";

    const std::string
    BToPiLeptonInclusiveNeutrinos::kinematics_description_s = "\
The invariant mass of the tau-nubar_tau pair in GeV^2.";

    const std::string
    BToPiLeptonInclusiveNeutrinos::kinematics_description_snunubar = "\
The invariant mass of the nu_tau-nubar_mu pair in GeV^2.";

    const std::string
    BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta = "\
The cosine of the charged lepton l's helicity angle theta in the l-nubar_l rest frame.";

    const std::string
    BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta_tau = "\
The cosine of the tau's helicity angle theta_tau in the tau-nubar_tau rest frame.";

    const std::string
    BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta_mu_star = "\
The cosine of the muon's helicity angle theta_mu^* in the tau's rest frame.";

    const std::string
    BToPiLeptonInclusiveNeutrinos::kinematics_description_phi = "\
The azimuthal angle between the B-D plane and the tau-nubar_tau plane.";
}
