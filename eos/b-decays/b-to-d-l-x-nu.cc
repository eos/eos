/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2025 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2018 Christoph Bobeth
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

#include <eos/b-decays/b-to-d-l-x-nu.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    using namespace std::literals::string_literals;

    template <>
    struct Implementation<BToDLeptonInclusiveNeutrinos>
    {
        std::shared_ptr<FormFactors<PToP>> form_factors;

        QuarkFlavorOption opt_q;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_D;

        UsedParameter m_mu;

        UsedParameter m_tau;

        UsedParameter g_fermi;

        UsedParameter hbar;

        SwitchOption opt_model;

        std::shared_ptr<Model> model;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            form_factors(FormFactorFactory<PToP>::create("B->D::" + o.get("form-factors"_ok, "BCL2008"), p, o)),
            opt_q(o, options, "q"_ok),
            m_B(p["mass::B_" + opt_q.str()], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            m_D(p["mass::D_" + opt_q.str()], u),
            m_mu(p["mass::mu"], u),
            m_tau(p["mass::tau"], u),
            g_fermi(p["WET::G_Fermi"], u),
            hbar(p["QM::hbar"], u),
            opt_model(o, "model"_ok, {"SM"}, "SM"),
            model(Model::make(opt_model.value(), p, o))
        {
            Context ctx("When constructing B->DlX observable");

            u.uses(*form_factors);
        }

        // normalized to N_1 = |V_vb|^2 G_F^2 / (192 pi^3 MB^3)
        double differential_decay_width_1nu_1var(const double & s) const
        {
            double fp = form_factors->f_p(s);
            double lam = lambda(m_B * m_B, m_D * m_D, s);

            return sqrt(lam) * (3.0 * fp * fp * lam);
        }

        double differential_decay_width_1nu(const double & s, const double & c_theta_mu) const
        {
            double fp = form_factors->f_p(s);
            double lam = lambda(m_B * m_B, m_D * m_D, s);

            return 3.0 / 4.0 * fp * fp * lam * sqrt(lam) * (1.0 - c_theta_mu * c_theta_mu);
        }

        // normalized to N_3 = |V_vb|^2 G_F^2 / (384 pi^3 MB^3)
        //                   * tau_tau / hbar * G_F^2 m_tau^5 / (192 pi^3)
        double differential_decay_width_3nu_1var(const double & s) const
        {
            double fp = form_factors->f_p(s);
            double f0 = form_factors->f_0(s);
            double lam = lambda(m_B * m_B, m_D * m_D, s);
            // make sure to return NaN if s < m_tau^2
            double sqrtv = sqrt(1.0 - m_tau() * m_tau() / s);
            double v = sqrtv * sqrtv, v2 = v * v;

            return sqrt(lam) * v2 * ((3.0 - v) * fp * fp * lam + 3.0 * (1.0 - v) * f0 * f0 * power_of<2>(m_B() * m_B() - m_D() * m_D())) * 4.0 / 3.0;
        }

        double differential_decay_width_3nu(const double & s, const double & snunubar,
                const double & z, const double & phi, const double & zst) const
        {
            const double fp = form_factors->f_p(s), fp2 = fp * fp;
            const double f0 = form_factors->f_0(s), f02 = f0 * f0;
            const double lam = lambda(m_B * m_B, m_D * m_D, s), sqrtlam = sqrt(lam);
            const double mtau = m_tau(), mtau2 = mtau * mtau, mtau4 = mtau2 * mtau2;
            const double mtau6 = mtau4 * mtau2, mtau8 = mtau4 * mtau4;
            const double mB2 = m_B() * m_B(), mD2 = m_D() * m_D();
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

            return (a + b * z + c * z2 + (d + e * z) * sqrt(1.0 - z2) * cos(phi));
        }

    };

    const std::vector<OptionSpecification>
    Implementation<BToDLeptonInclusiveNeutrinos>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToP>::option_specification(),
        { "q"_ok, { "d"s, "u"s }, "d"s }
    };

    BToDLeptonInclusiveNeutrinos::BToDLeptonInclusiveNeutrinos(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDLeptonInclusiveNeutrinos>(new Implementation<BToDLeptonInclusiveNeutrinos>(parameters, options, *this))
    {
    }

    BToDLeptonInclusiveNeutrinos::~BToDLeptonInclusiveNeutrinos()
    {
    }

    double
    BToDLeptonInclusiveNeutrinos::differential_decay_width_1nu(const double & s, const double & c_theta_mu) const
    {
        return _imp->differential_decay_width_1nu(s, c_theta_mu);
    }

    double
    BToDLeptonInclusiveNeutrinos::differential_decay_width_3nu(const double & s, const double & snunubar,
            const double & z, const double & phi, const double & zst) const
    {
        return _imp->differential_decay_width_3nu(s, snunubar, z, phi, zst);
    }

    double
    BToDLeptonInclusiveNeutrinos::integrated_decay_width_1nu(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonInclusiveNeutrinos>::differential_decay_width_1nu_1var,
                                                             _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max);
    }

    double
    BToDLeptonInclusiveNeutrinos::integrated_decay_width_3nu(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonInclusiveNeutrinos>::differential_decay_width_3nu_1var,
                                                             _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max);
    }

    const std::string
    BToDLeptonInclusiveNeutrinos::description = "\
The neutrino-inclusive decay B->D l X_nu, where l=e,mu is a light lepton, and \
X_nu = {nu, nu nubar nu} is either a one or three neutrino final state. The \
numerical implementations for this process follow [BIvD2016].";

    const std::string
    BToDLeptonInclusiveNeutrinos::kinematics_description_s = "\
The invariant mass of the tau-nubar_tau pair in GeV^2.";

    const std::string
    BToDLeptonInclusiveNeutrinos::kinematics_description_snunubar = "\
The invariant mass of the nu_tau-nubar_mu pair in GeV^2.";

    const std::string
    BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta = "\
The cosine of the charged lepton l's helicity angle theta in the l-nubar_l rest frame.";

    const std::string
    BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta_tau = "\
The cosine of the tau's helicity angle theta_tau in the tau-nubar_tau rest frame.";

    const std::string
    BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta_mu_star = "\
The cosine of the muon's helicity angle theta_mu^* in the tau's rest frame.";

    const std::string
    BToDLeptonInclusiveNeutrinos::kinematics_description_phi = "\
The azimuthal angle between the B-D plane and the tau-nubar_tau plane.";

    const std::set<ReferenceName>
    BToDLeptonInclusiveNeutrinos::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BToDLeptonInclusiveNeutrinos::begin_options()
    {
        return Implementation<BToDLeptonInclusiveNeutrinos>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToDLeptonInclusiveNeutrinos::end_options()
    {
        return Implementation<BToDLeptonInclusiveNeutrinos>::options.cend();
    }
}
