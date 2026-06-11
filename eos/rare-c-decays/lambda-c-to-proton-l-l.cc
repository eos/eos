/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Ahmet Kokulu
 * Copyright (c) 2019,2021 Danny van Dyk
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2026 Dominik Suelmann
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

#include <eos/form-factors/baryonic.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/rare-c-decays/lambda-c-to-proton-l-l.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <functional>
#include <map>
#include <string>

namespace eos
{
    using std::norm;

    namespace lambda_c_to_proton_l_l
    {
        struct Amplitudes
        {
                double aU11;
                double aL11;
                double aU22;
                double aL22;
                double aP12;
                double aS22;

                double beta_l;
                double beta_l_squared;
                double four_m_l_squared;
                double q2;
                double norm_squared;
        };

        struct AngularObservables
        {
                std::array<double, 3> _k;

                AngularObservables(const Amplitudes & a)
                {
                    using std::conj;
                    using std::imag;
                    using std::norm;
                    using std::real;
                    using std::sqrt;

                    double beta_l           = a.beta_l;
                    double beta_l_squared   = a.beta_l_squared;
                    double four_m_l_squared = a.four_m_l_squared;
                    double q2               = a.q2;
                    double aU11U22          = a.aU11 + a.aU22;
                    double aULS             = a.aU11 + a.aL11 + a.aS22;

                    // cf. [GHM:2021A], eq. (20), p. 6.

                    // K_{1ss}
                    _k[0] = a.norm_squared * (q2 * beta_l_squared * (aU11U22 / 2.0 + a.aL11 + a.aL22) + four_m_l_squared * aULS);

                    // K_{1cc}
                    _k[1] = a.norm_squared * (q2 * beta_l_squared * aU11U22 + four_m_l_squared * aULS);

                    // K_{1c}
                    _k[2] = a.norm_squared * (-2.0 * q2 * beta_l * a.aP12);
                }

                AngularObservables(const std::array<double, 3> & k) :
                    _k(k)
                {
                }

                inline double
                k1ss() const
                {
                    return _k[0];
                }

                inline double
                k1cc() const
                {
                    return _k[1];
                }

                inline double
                k1c() const
                {
                    return _k[2];
                }

                inline double
                decay_width() const
                {
                    return 2.0 * k1ss() + k1cc();
                }

                inline double
                a_fb_leptonic() const
                {
                    return 3.0 / 2.0 * k1c() / decay_width();
                }

                inline double
                a_fb_leptonic_num() const
                {
                    return 3.0 / 2.0 * k1c();
                }

                inline double
                f_l() const
                {
                    return (2.0 * k1ss() - k1cc()) / decay_width();
                }

                inline double
                f_l_num() const
                {
                    return (2.0 * k1ss() - k1cc());
                }

                inline double
                d2gamma(const double & c_lep) const
                {
                    const double c2_lep = c_lep * c_lep;
                    const double s2_lep = 1.0 - c2_lep;

                    // cf. [GHM:2021A], p. 6, eq. (9)
                    return 3.0 / 2.0 * (k1ss() * s2_lep + k1cc() * c2_lep + k1c() * c_lep);
                }
        };
    } // namespace lambda_c_to_proton_l_l

    template <> struct Implementation<LambdaCToProtonLeptonLepton>
    {
            std::shared_ptr<Model> model;

            Parameters parameters;

            UsedParameter hbar;
            UsedParameter tau_Lambda_c;

            UsedParameter g_fermi;
            UsedParameter alpha_e;

            LeptonFlavorOption opt_l;
            UsedParameter      m_l;

            UsedParameter m_Lambda_c;
            UsedParameter m_proton;
            UsedParameter m_rho;
            UsedParameter m_omega;
            UsedParameter m_phi;
            UsedParameter tau_rho;
            UsedParameter tau_omega;
            UsedParameter tau_phi;

            // resonance parameters
            UsedParameter a_rho;
            UsedParameter a_omega;
            UsedParameter a_phi;
            UsedParameter delta_rho;
            UsedParameter delta_omega_m_rho;
            UsedParameter delta_phi_m_rho;

            BooleanOption opt_cp_conjugate;

            UsedParameter mu;

            std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

            static const std::vector<OptionSpecification> options;

            Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
                model(Model::make(o.get("model"_ok, "WET"_ov), p, o)),
                parameters(p),
                hbar(p["QM::hbar"], u),
                tau_Lambda_c(p["life_time::Lambda_c"], u),
                g_fermi(p["WET::G_Fermi"], u),
                alpha_e(p["QED::alpha_e(m_c)"], u),
                opt_l(o, options, "l"_ok),
                m_l(p["mass::" + opt_l.str()], u),
                m_Lambda_c(p["mass::Lambda_c"], u),
                m_proton(p["mass::proton"], u),
                m_rho(p["mass::rho^0"], u),
                m_omega(p["mass::omega"], u),
                m_phi(p["mass::phi"], u),
                tau_rho(p["life_time::rho^0"], u),
                tau_omega(p["life_time::omega"], u),
                tau_phi(p["life_time::phi"], u),
                a_rho(p["Lambda_c->proton::res_a_rho@GHM2021"], u),
                a_omega(p["Lambda_c->proton::res_a_omega@GHM2021"], u),
                a_phi(p["Lambda_c->proton::res_a_phi@GHM2021"], u),
                delta_rho(p["Lambda_c->proton::res_delta_rho@GHM2021"], u),
                delta_omega_m_rho(p["Lambda_c->proton::res_delta_omega_m_rho@GHM2021"], u),
                delta_phi_m_rho(p["Lambda_c->proton::res_delta_phi_m_rho@GHM2021"], u),
                opt_cp_conjugate(o, options, "cp-conjugate"_ok),
                mu(p["uc::mu"], u),
                form_factors(FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_c->neutron::" + o.get("form-factors"_ok, "BMRvD2022"_ov).str(), p, o)) // using isospin
            {
                u.uses(*form_factors);
                u.uses(*model);
            }

            const double
            normalization(const double & q2) const
            {
                double lam = lambda(m_Lambda_c * m_Lambda_c, m_proton * m_proton, q2);

                if ((lam <= 0) || (q2 <= 4.0 * m_l * m_l))
                {
                    return 0.0;
                }
                // there is no overall CKM factor in c -> u transitions, see [dBMS:2016A] eq. 4, instead the individual factors are applied in the amplitude calculation.
                return g_fermi() * alpha_e() * std::sqrt(std::sqrt(1.0 - 4.0 * m_l * m_l / q2))
                       * std::sqrt(std::sqrt(lam) / 3.0 / power_of<11>(2) / power_of<5>(M_PI) / power_of<3>(m_Lambda_c()));
            }

            lambda_c_to_proton_l_l::Amplitudes
            amplitudes(const double & s)
            {
                using std::conj;
                using std::imag;
                using std::norm;
                using std::real;
                using std::sqrt;

                lambda_c_to_proton_l_l::Amplitudes result;

                const auto wc = model->wilson_coefficients_uc(opt_l.value(), opt_cp_conjugate.value());

                complex<double> ckm_factor_ds = model->ckm_ud() * conj(model->ckm_cd()) + model->ckm_us() * conj(model->ckm_cs());
                if (opt_cp_conjugate.value())
                {
                    ckm_factor_ds = conj(ckm_factor_ds);
                };

                const complex<double> c7minus  = (wc.c7() - wc.c7prime()) * ckm_factor_ds;
                const complex<double> c7plus   = (wc.c7() + wc.c7prime()) * ckm_factor_ds;
                const complex<double> c9minus  = (wc.c9() - wc.c9prime()) * ckm_factor_ds;
                const complex<double> c9plus   = (wc.c9() + wc.c9prime()) * ckm_factor_ds;
                const complex<double> c10minus = (wc.c10() - wc.c10prime()) * ckm_factor_ds;
                const complex<double> c10plus  = (wc.c10() + wc.c10prime()) * ckm_factor_ds;

                const double res_delta_rho   = delta_rho;
                const double res_delta_omega = delta_omega_m_rho;
                const double res_delta_phi   = delta_phi_m_rho;

                const double hb = hbar;

                const complex<double> c9R =
                        (a_rho() / (complex<double>(s - m_rho * m_rho, m_rho * hb / tau_rho))
                         + a_omega() * complex<double>(std::cos(res_delta_omega), std::sin(res_delta_omega)) / (complex<double>(s - m_omega * m_omega, m_omega * hb / tau_omega))
                         + a_phi() * complex<double>(std::cos(res_delta_phi), std::sin(res_delta_phi)) / (complex<double>(s - m_phi * m_phi, m_phi * hb / tau_phi)));
                const complex<double> cPhase(std::cos(res_delta_rho), std::sin(res_delta_rho));

                // baryonic form factors (10)
                const double fff0         = form_factors->f_time_v(s);  // f0
                const double fffplus      = form_factors->f_long_v(s);  // fplus
                const double fffperp      = form_factors->f_perp_v(s);  // fperp
                const double ffg0         = form_factors->f_time_a(s);  // g0
                const double ffgplus      = form_factors->f_long_a(s);  // gplus
                const double ffgperp      = form_factors->f_perp_a(s);  // gperp
                const double ffhplus      = form_factors->f_long_t(s);  // hplus
                const double ffhtildeplus = form_factors->f_long_t5(s); // htildeplus
                const double ffhperp      = form_factors->f_perp_t(s);  // hperp
                const double ffhtildeperp = form_factors->f_perp_t5(s); // htildeperp
                // running quark masses
                const double mcatmu       = model->m_c_msbar(mu);

                // kinematics
                const double beta_l_squared   = (1.0 - 4.0 * m_l * m_l / s);
                const double beta_l           = std::sqrt(beta_l_squared);
                const double four_m_l_squared = 4.0 * m_l * m_l;
                const double mplus            = m_Lambda_c + m_proton;
                const double mminus           = m_Lambda_c - m_proton;
                const double mplus_squared    = power_of<2>(mplus);
                const double mminus_squared   = power_of<2>(mminus);
                const double splus            = mplus_squared - s;
                const double sminus           = mminus_squared - s;
                const double sqrtsminussplus  = std::sqrt(sminus) * std::sqrt(splus);

                // normalization
                const complex<double> N = normalization(s);

                // cf. [GHM:2021A], eq. (11), p. 6, excluding contributions from scalar and tensor operators.
                result.aU11 = +4.0
                              * ((norm(c7plus * 2.0 * mcatmu / s * mplus * ffhperp + c9plus * fffperp)
                                  + 2.0 * real(conj(c7plus * 2.0 * mcatmu / s * mplus * ffhperp + c9plus * fffperp) * c9R * cPhase * fffperp) + norm(c9R * fffperp))
                                         * sminus
                                 + (norm(c7minus * 2.0 * mcatmu / s * mminus * ffhtildeperp + c9minus * ffgperp)
                                    + 2.0 * real(conj(c7minus * 2.0 * mcatmu / s * mminus * ffhtildeperp + c9minus * ffgperp) * c9R * cPhase * ffgperp) + norm(c9R * ffgperp))
                                           * splus);

                result.aL11 = +2.0 / s
                              * ((norm(c7plus * 2.0 * mcatmu * ffhplus + c9plus * mplus * fffplus)
                                  + 2.0 * real(conj(c7plus * 2.0 * mcatmu * ffhplus + c9plus * mplus * fffplus) * c9R * cPhase * mplus * fffplus) + norm(c9R * mplus * fffplus))
                                         * sminus
                                 + (norm(c7minus * 2.0 * mcatmu * ffhtildeplus + c9minus * mminus * ffgplus)
                                    + 2.0 * real(conj(c7minus * 2.0 * mcatmu * ffhtildeplus + c9minus * mminus * ffgplus) * c9R * cPhase * mminus * ffgplus)
                                    + norm(c9R * mminus * ffgplus))
                                           * splus);
                result.aU22 = +4.0 * (norm(c10plus) * fffperp * fffperp * sminus + norm(c10minus) * ffgperp * ffgperp * splus);
                result.aL22 = +2.0 / s * (norm(c10plus) * mplus_squared * fffplus * fffplus * sminus + norm(c10minus) * mminus_squared * ffgplus * ffgplus * splus);
                result.aS22 = +2.0 / s * (norm(c10plus) * mminus_squared * fff0 * fff0 * splus + norm(c10minus) * mplus_squared * ffg0 * ffg0 * sminus);
                result.aP12 =
                        -8
                        * (real(c7minus * conj(c10plus)) * mcatmu / s * mminus * fffperp * ffhtildeperp + real(c7plus * conj(c10minus)) * mcatmu / s * mplus * ffgperp * ffhperp
                           + real((wc.c9() * ckm_factor_ds + c9R * cPhase) * conj(wc.c10() * ckm_factor_ds) - wc.c9prime() * ckm_factor_ds * conj(wc.c10prime() * ckm_factor_ds))
                                     * ffgperp * fffperp)
                        * sqrtsminussplus;

                result.beta_l           = beta_l;
                result.beta_l_squared   = beta_l_squared;
                result.four_m_l_squared = four_m_l_squared;
                result.q2               = s;
                result.norm_squared     = norm(N);

                return result;
            }

            std::array<double, 3>
            _differential_angular_observables(const double & q2)
            {
                return lambda_c_to_proton_l_l::AngularObservables(this->amplitudes(q2))._k;
            }

            std::array<double, 3>
            _integrated_angular_observables(const double & q2_min, const double & q2_max)
            {
                std::function<std::array<double, 3>(const double &)> integrand = [this](const double & q2) -> std::array<double, 3>
                { return this->_differential_angular_observables(q2); };

                return integrate<1, 3>(integrand, q2_min, q2_max, cubature::Config().epsrel(1e-5));
            }

            inline lambda_c_to_proton_l_l::AngularObservables
            differential_angular_observables(const double & q2)
            {
                return lambda_c_to_proton_l_l::AngularObservables{ _differential_angular_observables(q2) };
            }

            inline lambda_c_to_proton_l_l::AngularObservables
            integrated_angular_observables(const double & q2_min, const double & q2_max)
            {
                return lambda_c_to_proton_l_l::AngularObservables{ _integrated_angular_observables(q2_min, q2_max) };
            }
    };

    const std::vector<OptionSpecification> Implementation<LambdaCToProtonLeptonLepton>::options{
        Model::option_specification(),
        { "cp-conjugate"_ok,     { "true"_ov, "false"_ov }, "false"_ov },
        {            "l"_ok, { "e"_ov, "mu"_ov, "tau"_ov },    "mu"_ov },
    };

    LambdaCToProtonLeptonLepton::LambdaCToProtonLeptonLepton(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<LambdaCToProtonLeptonLepton>(new Implementation<LambdaCToProtonLeptonLepton>(p, o, *this))
    {
    }

    LambdaCToProtonLeptonLepton::~LambdaCToProtonLeptonLepton() {}

    /* for double-differential signal PDF */
    double
    LambdaCToProtonLeptonLepton::double_differential_decay_width(const double & q2, const double & c_lep) const
    {
        return _imp->differential_angular_observables(q2).d2gamma(c_lep);
    }

    double
    LambdaCToProtonLeptonLepton::integrated_decay_width(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).decay_width();
    }

    /* q^2-differential observables */

    double
    LambdaCToProtonLeptonLepton::differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).decay_width() * _imp->tau_Lambda_c / _imp->hbar;
    }

    double
    LambdaCToProtonLeptonLepton::differential_a_fb_leptonic(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).a_fb_leptonic();
    }

    double
    LambdaCToProtonLeptonLepton::differential_f_l(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).f_l();
    }

    /* q^2-integrated observables */

    double
    LambdaCToProtonLeptonLepton::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).decay_width() * _imp->tau_Lambda_c / _imp->hbar;
    }

    double
    LambdaCToProtonLeptonLepton::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_a_fb_leptonic_num(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic_num();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_f_l(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_l();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_f_l_num(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_l_num();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k1ss(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1ss() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k1cc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1cc() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k1c(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1c() / o.decay_width();
    }

    const std::string LambdaCToProtonLeptonLepton::description = "\
      The decay Lambda_c -> proton lbar l, where lbar=e^+,mu^+,tau^+ is a charged antilepton.";

    const std::string LambdaCToProtonLeptonLepton::kinematics_description_q2 = "\
      The invariant mass of the lbar-l pair in GeV^2.";

    const std::string LambdaCToProtonLeptonLepton::kinematics_description_c_theta_l = "\
      The cosine of the helicity angle between the direction of flight of the charged antilepton and of the Lambda_c in the lbar-l rest frame.";

    const std::set<ReferenceName> LambdaCToProtonLeptonLepton::references{ "GHM:2021A"_rn, "BMRvD:2022A"_rn, "dBMS:2016A"_rn };

    std::vector<OptionSpecification>::const_iterator
    LambdaCToProtonLeptonLepton::begin_options()
    {
        return Implementation<LambdaCToProtonLeptonLepton>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    LambdaCToProtonLeptonLepton::end_options()
    {
        return Implementation<LambdaCToProtonLeptonLepton>::options.cend();
    }
} // namespace eos
