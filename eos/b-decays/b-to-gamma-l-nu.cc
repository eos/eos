/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#include <eos/b-decays/b-to-gamma-l-nu.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <tuple>

namespace eos
{
    using std::norm;

    /*
     * Decay: B_u -> gamma l nubar, cf. [BBJW:2018A], [BR:2011A]
     *
     * Lepton and neutrino are assumed to be massless
     */
    template <>
    struct Implementation<BToGammaLeptonNeutrino>
    {
        std::shared_ptr<Model> model;
        std::shared_ptr<FormFactors<PToGamma>> form_factors;

        UsedParameter alpha_qed;
        UsedParameter g_fermi;
        UsedParameter v_ub_abs;
        UsedParameter hbar;

        UsedParameter m_B;
        UsedParameter f_B;
        UsedParameter tau_B;

        static const constexpr double e_l = -1.0;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            form_factors(FormFactorFactory<PToGamma>::create("B->gamma::" + o.get("form-factors", "FLvD2022QCDF"), p, o)),
            alpha_qed(p["QED::alpha_e(m_b)"] ,u),
            g_fermi(p["WET::G_Fermi"], u),
            v_ub_abs(p["CKM::abs(V_ub)"], u),
            hbar(p["QM::hbar"], u),
            m_B(p["mass::B_u"], u),
            f_B(p["decay-constant::B_u"], u),
            tau_B(p["life_time::B_u"], u)
        {
            Context ctx("When constructing B->gammalnu observable");

            u.uses(*model);
            u.uses(*form_factors);
        }

        double fully_differential_decay_width(const double & E_gamma, const double & costheta) const
        {
            const double E_ell = 1.0 / 2.0 * ((-1.0 + costheta) * E_gamma + m_B);
            const double dEell_dcostheta = E_gamma / 2.0;

            const double x_gamma = 2.0 * E_gamma / m_B;
            const double x_ell   = 2.0 * E_ell   / m_B;
            const double x_nu    = 2.0 * (1.0 - (E_gamma + E_ell) / m_B);

            const double F_V = form_factors->F_V(E_gamma);
            const double F_A = form_factors->F_A(E_gamma) + e_l * f_B / E_gamma; // mind different definitions of F_A between [BBJW:2018A] and [BR:2011A]

            // cf. [BR:2011A], eq. (2.7)
            const double dGamma_dEgamma_dEell =
                alpha_qed * power_of<2>(g_fermi * v_ub_abs) / (16.0 * power_of<2>(M_PI))
                * power_of<3>(m_B) * (1.0 - x_gamma)
                * (
                      power_of<2>(1.0 - x_nu ) * power_of<2>(F_A + F_V)
                    + power_of<2>(1.0 - x_ell) * power_of<2>(F_A - F_V)
                  );
            return dGamma_dEgamma_dEell * dEell_dcostheta;
        }

        double differential_decay_width_dEgamma(const double & E_gamma) const
        {
            // Use analytic result of the angular integration
            // cf. [BBJW:2018A], eq. (2.6)
            const double prefactor = alpha_qed * power_of<2>(g_fermi * v_ub_abs) / (6.0 * power_of<2>(M_PI));
            const double F_V = form_factors->F_V(E_gamma);
            const double F_A = form_factors->F_A(E_gamma);

            return prefactor * m_B * power_of<3>(E_gamma) * (1.0 - 2.0 * E_gamma / m_B)
                * ( power_of<2>(std::abs(F_V)) + power_of<2>( F_A + e_l * f_B / E_gamma ) );
        }

        double integrated_decay_width(const double & E_gamma_min) const
        {
            return integrate<GSL::QAGS>([&](const double & E_gamma) { return differential_decay_width_dEgamma(E_gamma); }, E_gamma_min, m_B / 2.0);
        }

        double integrated_branching_ratio(const double & E_gamma_min) const
        {
            return integrated_decay_width(E_gamma_min) * tau_B / hbar;
        }

        std::tuple<double, double> forward_backward_decay_widths(const double & E_gamma_min) const
        {
            // Use analytic result of the angular integration

            const double m_B2 = m_B * m_B, m_B3 = m_B2 * m_B;

            auto dGamma_dEgamma_forward = [&](const double & E_gamma) -> double {
                const double F_V = form_factors->F_V(E_gamma);
                const double F_A = form_factors->F_A(E_gamma) + e_l * f_B / E_gamma; // mind different definitions of F_A between [BBJW:2018A] and [BR:2011A]

                return (2 * 1. / m_B3 * (2 * E_gamma - m_B) * (3 * F_A * F_V + 2 * power_of<2>(F_A) + 2 * power_of<2>(F_V)) * power_of<3>(E_gamma)) / 3.;
            };

            auto dGamma_dEgamma_backward = [&](const double & E_gamma) -> double {
                const double F_V = form_factors->F_V(E_gamma);
                const double F_A = form_factors->F_A(E_gamma) + e_l * f_B / E_gamma; // mind different definitions of F_A between [BBJW:2018A] and [BR:2011A]

                return (2 * 1. / m_B3 * (2 * E_gamma - m_B) * (-3 * F_A * F_V + 2 * power_of<2>(F_A) + 2 * power_of<2>(F_V)) * power_of<3>(E_gamma)) / 3.;
            };

            const double prefactor = alpha_qed * power_of<2>(g_fermi * v_ub_abs)
                    / (16.0 * power_of<2>(M_PI)) * power_of<3>(m_B);

            const double Gamma_forward  = prefactor * integrate<GSL::QAGS>(dGamma_dEgamma_forward,  E_gamma_min, m_B / 2.0);
            const double Gamma_backward = prefactor * integrate<GSL::QAGS>(dGamma_dEgamma_backward, E_gamma_min, m_B / 2.0);

            return { Gamma_forward, Gamma_backward };
        }

        double forward_backward_asymmetry(const double & E_gamma_min) const
        {
            const auto [Gamma_forward, Gamma_backward] = forward_backward_decay_widths(E_gamma_min);

            return (Gamma_forward - Gamma_backward) / (Gamma_forward + Gamma_backward);
        }

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            // check consitency of Gamma_F_Gamma_B with forward_backward_asymmetry
            {
                const constexpr double E_gamma_min = 2.0;
                auto [Gamma_F, Gamma_B] = this->forward_backward_decay_widths(E_gamma_min);
                const double Gamma = this->integrated_decay_width(E_gamma_min);

                results.add({ Gamma_F + Gamma_B - Gamma, "Gamma_F + Gamma_B - Gamma" });
            }

            return results;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToGammaLeptonNeutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToGamma>::option_specification()
    };

    BToGammaLeptonNeutrino::BToGammaLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToGammaLeptonNeutrino>(new Implementation<BToGammaLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToGammaLeptonNeutrino::~BToGammaLeptonNeutrino()
    {
    }

    double
    BToGammaLeptonNeutrino::integrated_branching_ratio(const double & E_gamma_min) const
    {
        return _imp->integrated_branching_ratio(E_gamma_min);
    }

    double
    BToGammaLeptonNeutrino::forward_backward_asymmetry(const double & E_gamma_min) const
    {
        return _imp->forward_backward_asymmetry(E_gamma_min);
    }

    double
    BToGammaLeptonNeutrino::fully_differential_decay_width(const double & E_gamma, const double & costheta) const
    {
        return _imp->fully_differential_decay_width(E_gamma, costheta);
    }

    double
    BToGammaLeptonNeutrino::differential_decay_width_dEgamma(const double & E_gamma) const
    {
        return _imp->differential_decay_width_dEgamma(E_gamma);
    }

    double
    BToGammaLeptonNeutrino::integrated_decay_width(const double & E_gamma_min) const
    {
        return _imp->integrated_decay_width(E_gamma_min);
    }

    Diagnostics
    BToGammaLeptonNeutrino::diagnostics() const
    {
        return _imp->diagnostics();
    }

    const std::string
    BToGammaLeptonNeutrino::description = "\
The decay B_u -> gamma l nu, where l=e, mu, tau is a lepton.";

    const std::string
    BToGammaLeptonNeutrino::kinematics_description_Egamma = "\
The energy of the photon in the B meson rest frame. The approach of Ref. [BBJW:2018A] is valid in the region Egamma > 1.5 GeV.";

    const std::string
    BToGammaLeptonNeutrino::kinematics_description_c_theta_l = "\
The cosine of the polar angle theta_l between the charged lepton and the direction opposite to the photon in the l-nubar rest frame.";

    const std::set<ReferenceName>
    BToGammaLeptonNeutrino::references
    {
        "BBJW:2018A"_rn,
        "BR:2011A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BToGammaLeptonNeutrino::begin_options()
    {
        return Implementation<BToGammaLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToGammaLeptonNeutrino::end_options()
    {
        return Implementation<BToGammaLeptonNeutrino>::options.cend();
    }
}
