/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018-2019 Ahmet Kokulu
 * Copyright (c) 2019-2025 Danny van Dyk
 * Copyright (c) 2021      Christoph Bobeth
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

#include <eos/b-decays/b-to-vec-l-nu-impl.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/models/model.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <functional>
#include <map>
#include <string>

namespace eos
{
    using namespace std::literals::string_literals;
    using std::norm;

    template <> struct Implementation<BToVectorLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        Parameters parameters;

        QuarkFlavorOption opt_q;
        RestrictedOption opt_V;

        UsedParameter hbar;

        UsedParameter tau_B;

        UsedParameter g_fermi;

        LeptonFlavorOption opt_l;

        UsedParameter m_l;

        UsedParameter m_B;

        UsedParameter m_V;

        const double isospin_factor;

        BooleanOption opt_cp_conjugate;

        UsedParameter mu;

        std::function<double (const double &)> m_U_msbar;
        std::function<complex<double> ()> v_Ub;
        std::function<WilsonCoefficients<ChargedCurrent> (LeptonFlavor, bool)> wc;

        SwitchOption opt_int_points;

        int int_points;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        using IntermediateResult = BToVectorLeptonNeutrino::IntermediateResult;

        IntermediateResult intermediate_result;

        static const std::vector<OptionSpecification> options;

        // { q, V } -> { process, U, B_name, V_name, c_I }
        // q: u, d, s: the spectar quark flavor
        // V: D, K, pi: the type of daughter meson
        // process: string that can be used to obtain the form factor
        // U: the quark flavor in the weak transition
        // B_name: name of the B meson
        // V_name: name of the daughter meson
        // c_I: isospin factor by which the amplitudes are multiplied
        static const std::map<std::tuple<QuarkFlavor, std::string>, std::tuple<std::string, QuarkFlavor, std::string, std::string, double>> process_map;

        inline std::string _process() const
        {
            const QuarkFlavor q = opt_q.value();
            const std::string V = opt_V.value();
            const auto p = process_map.find(std::make_tuple(q, V));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of q = " + stringify(q) + ", V = " + V);

            return std::get<0>(p->second);
        }

        inline QuarkFlavor _U() const
        {
            const QuarkFlavor q = opt_q.value();
            const std::string V = opt_V.value();
            const auto p = process_map.find(std::make_tuple(q, V));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of q = " + stringify(q) + ", V = " + V);

            return std::get<1>(p->second);
        }

        inline std::string _B() const
        {
            const QuarkFlavor q = opt_q.value();
            const std::string V = opt_V.value();
            const auto p = process_map.find(std::make_tuple(q, V));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of q = " + stringify(q) + ", V = " + V);

            return std::get<2>(p->second);
        }

        inline std::string _V() const
        {
            const QuarkFlavor q = opt_q.value();
            const std::string V = opt_V.value();
            const auto p = process_map.find(std::make_tuple(q, V));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of q = " + stringify(q) + ", V = " + V);

            return std::get<3>(p->second);
        }

        inline double _isospin_factor() const
        {
            const QuarkFlavor q = opt_q.value();
            const std::string V = opt_V.value();
            const auto p = process_map.find(std::make_tuple(q, V));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of q = " + stringify(q) + ", V = " + V);

            return std::get<4>(p->second);
        }

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            parameters(p),
            opt_q(o, options, "q"_ok),
            opt_V(o, options, "V"_ok),
            hbar(p["QM::hbar"], u),
            tau_B(p["life_time::" + _B()], u),
            g_fermi(p["WET::G_Fermi"], u),
            opt_l(o, options, "l"_ok),
            m_l(p["mass::" + opt_l.str()], u),
            m_B(p["mass::" + _B()], u),
            m_V(p["mass::" + _V()], u),
            isospin_factor(_isospin_factor()),
            opt_cp_conjugate(o, options, "cp-conjugate"_ok),
            mu(p[stringify(_U()) + "b" + opt_l.str() + "nu" + opt_l.str() + "::mu"], u),
            opt_int_points(o, "integration-points"_ok, {"256", "4096"}, "256"),
            int_points(destringify<int>(opt_int_points.value())),
            form_factors(FormFactorFactory<PToV>::create(_process() + "::" + o.get("form-factors"_ok, "BSZ2015"), p, o))
        {
            Context ctx("When constructing B->Vlnu observable");

            using std::placeholders::_1;
            using std::placeholders::_2;

            switch (_U())
            {
                case QuarkFlavor::up:
                    m_U_msbar = std::bind(&ModelComponent<components::QCD>::m_u_msbar, model.get(), _1);
                    v_Ub      = std::bind(&ModelComponent<components::CKM>::ckm_ub, model.get());
                    wc        = std::bind(&ModelComponent<components::WET::UBLNu>::wet_ublnu, model.get(), _1, _2);
                    break;
                case QuarkFlavor::charm:
                    m_U_msbar = std::bind(&ModelComponent<components::QCD>::m_c_msbar, model.get(), _1);
                    v_Ub      = std::bind(&ModelComponent<components::CKM>::ckm_cb, model.get());
                    wc        = std::bind(&ModelComponent<components::WET::CBLNu>::wet_cblnu, model.get(), _1, _2);
                    break;
                default:
                    throw InternalError("Invalid quark flavor: " + stringify(_U()));
            }

            u.uses(*form_factors);
            u.uses(*model);
        }

        // normalization cf. [DSD2014] eq. (7), p. 5
        double norm(const double & q2) const
        {
            const double lam  = lambda(m_B * m_B, m_V * m_V, q2);
            const double p    = (lam > 0.0) ? std::sqrt(lam) / (2.0 * m_B) : 0.0;

            // normalized prefactor (|Vcb|^2=1)
            return power_of<2>(g_fermi()) * p * q2 * power_of<2>(1.0 - m_l * m_l / q2) / (3.0 * 64.0 * power_of<3>(M_PI) * m_B * m_B);
        }

        b_to_vec_l_nu::Amplitudes amplitudes(const double & q2) const
        {
            b_to_vec_l_nu::Amplitudes result;

            // NP contributions in EFT including tensor operator cf. [DSD2014], p. 3
            const WilsonCoefficients<ChargedCurrent> wc = this->wc(opt_l.value(), opt_cp_conjugate.value());
            const complex<double> gV_pl = wc.cvl() + wc.cvr();  // gV_pl = 1 + gV = 1 + VL + VR = cVL + cVR
            const complex<double> gV_mi = wc.cvl() - wc.cvr();  // gV_mi = 1 - gA = 1 + VL - VR = cVL - cVR
            const complex<double> gP = wc.csr() - wc.csl();
            const complex<double> TL = wc.ct();

            // form factors
            const double aff0  = form_factors->a_0(q2);
            const double aff1  = form_factors->a_1(q2);
            const double aff12 = form_factors->a_12(q2);
            const double vff   = form_factors->v(q2);
            const double tff1  = form_factors->t_1(q2);
            const double tff2  = form_factors->t_2(q2);
            const double tff3  = form_factors->t_3(q2);
            // meson & lepton masses
            const double m_l = this->m_l();
            const double m_B = this->m_B();
            const double m_V = this->m_V();
            // running quark masses
            const double mbatmu = model->m_b_msbar(mu);
            const double mUatmu = this->m_U_msbar(mu);
            // kinematic variables
            const double lam      = lambda(m_B * m_B, m_V * m_V, q2);
            const double sqrt_lam = (lam > 0.0) ? std::sqrt(lam) : 0.0;
            const double sqrtq2   = std::sqrt(q2);
            // isospin factor
            const double isospin = this->isospin_factor;

            // transversity amplitudes A's. cf. [DSD2014], p.17
            if ((q2 >= power_of<2>(m_l)) && (q2 <= power_of<2>(m_B - m_V))) {
                result.a_0          = isospin * gV_mi * 8.0 * m_B * m_V / sqrtq2 * aff12;
                result.a_0_T        = isospin * TL / (2.0 * m_V) * ( (m_B * m_B + 3.0 * m_V * m_V - q2) * tff2 - lam * tff3 / (m_B * m_B - m_V * m_V) );
                result.a_plus       = isospin * ( (m_B + m_V) * aff1 * gV_mi - sqrt_lam * vff * gV_pl / (m_B + m_V) );
                result.a_minus      = isospin * ( (m_B + m_V) * aff1 * gV_mi + sqrt_lam * vff * gV_pl / (m_B + m_V) );
                result.a_plus_T     = isospin * TL / sqrtq2 * ( (m_B * m_B - m_V * m_V) * tff2 + sqrt_lam * tff1 );
                result.a_minus_T    = isospin * TL / sqrtq2 * ( (m_B * m_B - m_V * m_V) * tff2 - sqrt_lam * tff1 );
                result.a_t          = isospin * sqrt_lam * aff0 * gV_mi / sqrtq2;
                result.a_P          = isospin * sqrt_lam * aff0 * gP / (mbatmu + mUatmu);
                result.a_para       = (result.a_plus + result.a_minus) / std::sqrt(2.0);
                result.a_para_T     = (result.a_plus_T + result.a_minus_T) / std::sqrt(2.0);
                result.a_perp       = (result.a_plus - result.a_minus) / std::sqrt(2.0);
                result.a_perp_T     = (result.a_plus_T - result.a_minus_T) / std::sqrt(2.0);

                result.mlH = (m_l > 0.0) ? std::sqrt(m_l * m_l / q2) : 0.0;
                result.NF = norm(q2);
            }
            else
            {
                result.a_0          = 0.0;
                result.a_0_T        = 0.0;
                result.a_plus       = 0.0;
                result.a_minus      = 0.0;
                result.a_plus_T     = 0.0;
                result.a_minus_T    = 0.0;
                result.a_t          = 0.0;
                result.a_P          = 0.0;
                result.a_para       = 0.0;
                result.a_para_T     = 0.0;
                result.a_perp       = 0.0;
                result.a_perp_T     = 0.0;

                result.mlH = 0.0;
                result.NF = 0.0;
            }

            return result;

        }

        std::array<double, 12> _differential_angular_observables(const double & q2) const
        {
            return b_to_vec_l_nu::AngularObservables(this->amplitudes(q2))._vv;
        }

        // define below integrated observables in generic form
        std::array<double, 12> _integrated_angular_observables(const double & q2_min, const double & q2_max) const
        {
            std::function<std::array<double, 12> (const double &)> integrand(std::bind(&Implementation::_differential_angular_observables, this, std::placeholders::_1));
            // second argument of integrate1D is some power of 2
            return integrate1D(integrand, int_points, q2_min, q2_max);
        }

        inline b_to_vec_l_nu::AngularObservables differential_angular_observables(const double & q2) const
        {
            return b_to_vec_l_nu::AngularObservables{ _differential_angular_observables(q2) };
        }

        inline b_to_vec_l_nu::AngularObservables integrated_angular_observables(const double & q2_min, const double & q2_max) const
        {
            return b_to_vec_l_nu::AngularObservables{ _integrated_angular_observables(q2_min, q2_max) };
        }

        const IntermediateResult * prepare(const double & q2_min, const double & q2_max)
        {
            intermediate_result.ao = integrated_angular_observables(q2_min, q2_max);
            return &intermediate_result;
        }

        double normalized_decay_width(const double & q2) const
        {
            return this->differential_angular_observables(q2).normalized_decay_width();
        }

        double integrated_pdf_q2(const double & q2_min, const double & q2_max) const
        {
            const double q2_abs_min = power_of<2>(m_l());
            const double q2_abs_max = power_of<2>(m_B() - m_V());

            std::function<double (const double &)> f = std::bind(&Implementation<BToVectorLeptonNeutrino>::normalized_decay_width, this, std::placeholders::_1);
            const double num   = integrate<GSL::QAGS>(f, q2_min,     q2_max);
            const double denom = integrate<GSL::QAGS>(f, q2_abs_min, q2_abs_max);

            return num / denom / (q2_max - q2_min);
        }

        double integrated_pdf_w(const double & w_min, const double & w_max) const
        {
            const double m_B    = this->m_B(), m_B2 = m_B * m_B;
            const double m_V    = this->m_V(), m_V2 = m_V * m_V;
            const double q2_max = m_B2 + m_V2 - 2.0 * m_B * m_V * w_min;
            const double q2_min = m_B2 + m_V2 - 2.0 * m_B * m_V * w_max;

            return integrated_pdf_q2(q2_min, q2_max) * (q2_max - q2_min) / (w_max - w_min);
        }

        double differential_pdf_l(const double & c_theta_l, const double & q2_min, const double & q2_max)
        {
            b_to_vec_l_nu::AngularObservables ao = integrated_angular_observables(q2_min, q2_max);

            const double a_fb = ao.a_fb_leptonic();
            const double ftilde_L = ao.ftilde_L();

            return 0.5 + a_fb * c_theta_l + (1.0 - 3.0 * ftilde_L) / 8.0 * (3.0 * c_theta_l * c_theta_l - 1.0);
        }

        double differential_pdf_v(const double & c_theta_v, const double & q2_min, const double & q2_max)
        {
            b_to_vec_l_nu::AngularObservables ao = integrated_angular_observables(q2_min, q2_max);

            const double f_L = ao.f_L();

            return 3.0 / 4.0 * ((1.0 - f_L) - (1.0 - 3.0 * f_L) * c_theta_v * c_theta_v);
        }

        double differential_pdf_phi(const double & phi, const double & q2_min, const double & q2_max)
        {
            b_to_vec_l_nu::AngularObservables ao = integrated_angular_observables(q2_min, q2_max);

            const double J3 = 3.0 / 4.0 * ao.vv4T();
            const double J9 = 3.0 / 4.0 * ao.vv5T();

            return (0.5 + 2.0 / 3.0 * J3 * std::cos(2.0 * phi) + 2.0 / 3.0 * J9 * std::sin(2.0 * phi)) / M_PI;
        }

        double integrated_pdf_l(const double & c_theta_l_min, const double & c_theta_l_max, const double & q2_min, const double & q2_max)
        {
            b_to_vec_l_nu::AngularObservables ao = integrated_angular_observables(q2_min, q2_max);

            const double a_fb = ao.a_fb_leptonic();
            const double ftilde_L = ao.ftilde_L();

            return 0.5 * (c_theta_l_max - c_theta_l_min)
                + a_fb * 0.5 * (c_theta_l_max * c_theta_l_max - c_theta_l_min * c_theta_l_min)
                + (1.0 - 3.0 * ftilde_L) / 8.0 * ((c_theta_l_max * c_theta_l_max * c_theta_l_max - c_theta_l_min * c_theta_l_min * c_theta_l_min) - (c_theta_l_max - c_theta_l_min));
        }

        double integrated_pdf_v(const double & c_theta_v_min, const double & c_theta_v_max, const double & q2_min, const double & q2_max)
        {
            b_to_vec_l_nu::AngularObservables ao = integrated_angular_observables(q2_min, q2_max);

            const double f_L = ao.f_L();

            return 3.0 / 4.0 * (
                (1.0 - f_L) * (c_theta_v_max - c_theta_v_min)
                - (1.0 - 3.0 * f_L) * (c_theta_v_max * c_theta_v_max * c_theta_v_max - c_theta_v_min * c_theta_v_min * c_theta_v_min) / 3.0
            );
        }

        double integrated_pdf_phi(const double & phi_min, const double & phi_max, const double & q2_min, const double & q2_max)
        {
            b_to_vec_l_nu::AngularObservables ao = integrated_angular_observables(q2_min, q2_max);

            const double J3 = 3.0 / 4.0 * ao.vv4T();
            const double J9 = 3.0 / 4.0 * ao.vv5T();

            return (0.5 * (phi_max - phi_min)
                + 1.0 / 3.0 * J3 * (std::sin(2.0 * phi_max) - std::sin(2.0 * phi_min))
                - 1.0 / 3.0 * J9 * (std::cos(2.0 * phi_max) - std::cos(2.0 * phi_min))
            ) / M_PI;
        }
    };

    const std::map<std::tuple<QuarkFlavor, std::string>, std::tuple<std::string, QuarkFlavor, std::string, std::string, double>>
    Implementation<BToVectorLeptonNeutrino>::Implementation::process_map
    {
        { { QuarkFlavor::up,      "D^*"   }, { "B->D^*",     QuarkFlavor::charm, "B_u", "D_u^*", 1.0                  } },
        { { QuarkFlavor::down,    "D^*"   }, { "B->D^*",     QuarkFlavor::charm, "B_d", "D_d^*", 1.0                  } },
        { { QuarkFlavor::strange, "D_s^*" }, { "B_s->D_s^*", QuarkFlavor::charm, "B_s", "D_s^*", 1.0                  } },
        { { QuarkFlavor::up,      "rho"   }, { "B->rho",     QuarkFlavor::up,    "B_u", "rho^0", 1.0 / std::sqrt(2.0) } },
        { { QuarkFlavor::up,      "omega" }, { "B->omega",   QuarkFlavor::up,    "B_u", "omega", 1.0 / std::sqrt(2.0) } },
        { { QuarkFlavor::down,    "rho"   }, { "B->rho",     QuarkFlavor::up,    "B_d", "rho^+", 1.0                  } },
        { { QuarkFlavor::strange, "K^*"   }, { "B_s->K^*",   QuarkFlavor::up,    "B_s", "K_u^*", 1.0                  } },
    };

    const std::vector<OptionSpecification>
    Implementation<BToVectorLeptonNeutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToV>::option_specification(),
        { "V"_ok,            { "D^*"s, "D_s^*"s, "rho"s, "omega"s, "K^*"s }, ""s      },
        { "cp-conjugate"_ok, { "true"s, "false"s },                          "false"s },
        { "l"_ok,            { "e"s, "mu"s, "tau"s },                        "mu"s    },
        { "q"_ok,            { "u"s, "d"s, "s"s },                           "d"s     },
    };

    BToVectorLeptonNeutrino::BToVectorLeptonNeutrino(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BToVectorLeptonNeutrino>(new Implementation<BToVectorLeptonNeutrino>(p, o, *this))
    {
    }

    BToVectorLeptonNeutrino::~BToVectorLeptonNeutrino()
    {
    }

    /* q^2-differential observables */

    // |Vcb|=1
    double
    BToVectorLeptonNeutrino::normalized_differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).normalized_decay_width() * _imp->tau_B / _imp->hbar;
    }

    double
    BToVectorLeptonNeutrino::differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).normalized_decay_width() * std::norm(_imp->v_Ub()) * _imp->tau_B / _imp->hbar;
    }

    double
    BToVectorLeptonNeutrino::differential_a_fb_leptonic(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).a_fb_leptonic();
    }

    double
    BToVectorLeptonNeutrino::differential_J1c(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv10();
    }

    double
    BToVectorLeptonNeutrino::differential_J1s(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv1T();
    }

    double
    BToVectorLeptonNeutrino::differential_J2c(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv20();
    }

    double
    BToVectorLeptonNeutrino::differential_J2s(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv2T();
    }

    double
    BToVectorLeptonNeutrino::differential_J3(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv4T();
    }

    double
    BToVectorLeptonNeutrino::differential_J4(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv10T();
    }

    double
    BToVectorLeptonNeutrino::differential_J5(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv20T();
    }

    double
    BToVectorLeptonNeutrino::differential_J6c(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv30();
    }

    double
    BToVectorLeptonNeutrino::differential_J6s(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv3T();
    }

    double
    BToVectorLeptonNeutrino::differential_J7(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv30T();
    }

    double
    BToVectorLeptonNeutrino::differential_J8(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv40T();
    }

    double
    BToVectorLeptonNeutrino::differential_J9(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv5T();
    }

    /* q^2-integrated observables */

    const BToVectorLeptonNeutrino::IntermediateResult *
    BToVectorLeptonNeutrino::prepare(const double & q2_min, const double & q2_max) const
    {
        return _imp->prepare(q2_min, q2_max);
    }

    // |Vcb|=1
    double
    BToVectorLeptonNeutrino::normalized_decay_width(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).normalized_decay_width();
    }

    double
    BToVectorLeptonNeutrino::normalized_integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).normalized_decay_width() * _imp->tau_B / _imp->hbar;
    }

    double
    BToVectorLeptonNeutrino::integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).normalized_decay_width() * std::norm(_imp->v_Ub()) * _imp->tau_B / _imp->hbar;
    }

    double
    BToVectorLeptonNeutrino::integrated_branching_ratio_perp(const double & kperp_min, const double & kperp_max) const
    {
        std::function<double (const std::array<double, 2> &)> integrand = [this](const std::array<double, 2> & kinematics) -> double
        {
            const auto & [kperp, z_B] = kinematics;

            const double m_B = this->_imp->m_B(), m_V = this->_imp->m_V();
            const double kvec2 = kperp * kperp / (1.0 - z_B * z_B);
            const double q2 = m_B * m_B + m_V * m_V - 2 * m_B * sqrt(m_V * m_V + kvec2);

            const double jacobian1 = m_B / ((1 - z_B * z_B) * sqrt(m_V * m_V + kvec2));
            const double jacobian2 = 2.0 * kperp;
            const double jacobian  = jacobian1 * jacobian2;

            return this->differential_branching_ratio(q2) * jacobian / 2.0;
        };

        static cubature::Config config = cubature::Config().epsrel(0.5e-3).epsabs(1.0e-9);

        return integrate(integrand, std::array<double, 2>{kperp_min, -1.0}, std::array<double, 2>{kperp_max, 1.0},
                         config);
    }

    double
    BToVectorLeptonNeutrino::integrated_pdf_q2(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_pdf_q2(q2_min, q2_max);
    }

    double
    BToVectorLeptonNeutrino::integrated_pdf_w(const double & w_min, const double & w_max) const
    {
        return _imp->integrated_pdf_w(w_min, w_max);
    }

    double
    BToVectorLeptonNeutrino::differential_pdf_q2(const double & q2) const
    {
        const double q2_min = power_of<2>(_imp->m_l());
        const double q2_max = power_of<2>(_imp->m_B() - _imp->m_V());

        return _imp->differential_angular_observables(q2).normalized_decay_width() / _imp->integrated_angular_observables(q2_min, q2_max).normalized_decay_width();
    }

    double
    BToVectorLeptonNeutrino::differential_pdf_w(const double & w) const
    {
        const double m_B    = _imp->m_B();
        const double m_V    = _imp->m_V();
        const double q2     = m_B * m_B + m_V * m_V - 2.0 * m_B * m_V * w;
        const double q2_min = power_of<2>(_imp->m_l());
        const double q2_max = power_of<2>(m_B - m_V);

        // Jacobian |dq2/dw|
        const double jacobian = 2.0 * m_B * m_V;

        return jacobian * _imp->differential_angular_observables(q2).normalized_decay_width() / _imp->integrated_angular_observables(q2_min, q2_max).normalized_decay_width();
    }

    double
    BToVectorLeptonNeutrino::differential_pdf_l(const double & c_theta_l) const
    {
        const double q2_min = power_of<2>(_imp->m_l());
        const double q2_max = power_of<2>(_imp->m_B() - _imp->m_V());

        return _imp->differential_pdf_l(c_theta_l, q2_min, q2_max);
    }

    double
    BToVectorLeptonNeutrino::differential_pdf_v(const double & c_theta_v) const
    {
        const double q2_min = power_of<2>(_imp->m_l());
        const double q2_max = power_of<2>(_imp->m_B() - _imp->m_V());

        return _imp->differential_pdf_v(c_theta_v, q2_min, q2_max);
    }

    double
    BToVectorLeptonNeutrino::differential_pdf_phi(const double & phi) const
    {
        const double q2_min = power_of<2>(_imp->m_l());
        const double q2_max = power_of<2>(_imp->m_B() - _imp->m_V());

        return _imp->differential_pdf_phi(phi, q2_min, q2_max);
    }

    double
    BToVectorLeptonNeutrino::integrated_pdf_l(const double & c_theta_l_min, const double & c_theta_l_max) const
    {
        const double q2_min = power_of<2>(_imp->m_l());
        const double q2_max = power_of<2>(_imp->m_B() - _imp->m_V());

        return _imp->integrated_pdf_l(c_theta_l_min, c_theta_l_max, q2_min, q2_max);
    }

    double
    BToVectorLeptonNeutrino::integrated_pdf_v(const double & c_theta_v_min, const double & c_theta_v_max) const
    {
        const double q2_min = power_of<2>(_imp->m_l());
        const double q2_max = power_of<2>(_imp->m_B() - _imp->m_V());

        return _imp->integrated_pdf_v(c_theta_v_min, c_theta_v_max, q2_min, q2_max);
    }

    double
    BToVectorLeptonNeutrino::integrated_pdf_phi(const double & phi_min, const double & phi_max) const
    {
        const double q2_min = power_of<2>(_imp->m_l());
        const double q2_max = power_of<2>(_imp->m_B() - _imp->m_V());

        return _imp->integrated_pdf_phi(phi_min, phi_max, q2_min, q2_max);
    }

    double
    BToVectorLeptonNeutrino::integrated_a_fb_leptonic(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.a_fb_leptonic();
    }

    double
    BToVectorLeptonNeutrino::integrated_amplitude_polarization_L(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.normalized_amplitude_polarization_L() * std::norm(_imp->v_Ub());
    }

    double
    BToVectorLeptonNeutrino::integrated_amplitude_polarization_T(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.normalized_amplitude_polarization_T() * std::norm(_imp->v_Ub());
    }

    double
    BToVectorLeptonNeutrino::integrated_f_L(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.f_L();
    }

    double
    BToVectorLeptonNeutrino::integrated_ftilde_L(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.ftilde_L();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_c_1(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.a_c_1();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_c_2(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.a_c_2();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_c_3(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.a_c_3();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_t_1(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.a_t_1();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_t_2(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.a_t_2();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_t_3(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return o.a_t_3();
    }

    double
    BToVectorLeptonNeutrino::integrated_J1c(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv10();
    }

    double
    BToVectorLeptonNeutrino::integrated_J1s(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv1T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J2c(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv20();
    }

    double
    BToVectorLeptonNeutrino::integrated_J2s(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv2T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J3(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv4T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J4(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv10T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J5(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv20T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J6c(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv30();
    }

    double
    BToVectorLeptonNeutrino::integrated_J6s(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv3T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J7(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv30T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J8(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv40T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J9(const IntermediateResult * ir) const
    {
        const auto & o   = ir->ao;
        return 3.0 / 4.0 * o.vv5T();
    }


    //* cf. [DSD2014], eq. (6), p. 5 - normalized(|Vcb|=1)
    double
    BToVectorLeptonNeutrino::normalized_four_differential_decay_width(const double & q2, const double & c_theta_l, const double & c_theta_d, const double & phi) const
    {
        // compute d^4 Gamma, cf. [DSD2014], p. 5, eq. (6)
        // Trigonometric identities: Cosine squared of the angles
        double c_theta_d_2 = c_theta_d * c_theta_d;
        double c_theta_l_2 = c_theta_l * c_theta_l;
        double c_phi = std::cos(phi);
        // Sine squared of the angles
        double s_theta_d_2 = 1.0 - c_theta_d_2;
        double s_theta_l_2 = 1.0 - c_theta_l_2;
        // Sine of the angles
        double s_theta_d = std::sqrt(s_theta_d_2);
        double s_theta_l = std::sqrt(s_theta_l_2);
        double s_phi = std::sin(phi);
        // Cosine of twice the angle
        //double c_2_theta_d = 2.0 * c_theta_d_2 - 1.0;
        double c_2_theta_l = 2.0 * c_theta_l_2 - 1.0;
        double c_2_phi = std::cos(2.0 * phi);
        // Sine of twice the angle
        double s_2_theta_d = 2.0 * s_theta_d * c_theta_d;
        double s_2_theta_l = 2.0 * s_theta_l * c_theta_l;
        double s_2_phi = std::sin(2.0 * phi);

        b_to_vec_l_nu::AngularObservables a_o = _imp->differential_angular_observables(q2);

        double result = 9.0 / 32.0 / M_PI * (
                                            (a_o.vv10() + a_o.vv20() * c_2_theta_l + a_o.vv30() * c_theta_l) * c_theta_d_2
                                            +  (a_o.vv1T() + a_o.vv2T() * c_2_theta_l + a_o.vv3T() * c_theta_l) * s_theta_d_2
                                            +  a_o.vv4T() * s_theta_d_2 * s_theta_l_2 * c_2_phi + a_o.vv10T() * s_2_theta_d * s_2_theta_l * c_phi
                                            +  a_o.vv20T() * s_2_theta_d * s_theta_l * c_phi + a_o.vv5T() * s_theta_d_2 * s_theta_l_2 * s_2_phi
                                            +  a_o.vv30T() * s_2_theta_d * s_theta_l * s_phi +  a_o.vv40T() * s_2_theta_d * s_2_theta_l * s_phi
                                            );

        return result;
    }

    const std::string
    BToVectorLeptonNeutrino::description = "\
    The decay B->D^* l nu, where l=e,mu,tau is a lepton.";

    const std::string
    BToVectorLeptonNeutrino::kinematics_description_q2 = "\
    The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    BToVectorLeptonNeutrino::kinematics_description_c_theta_l = "\
    The cosine of the charged lepton's helicity angle theta_l in the l-nubar rest frame.";

    const std::string
    BToVectorLeptonNeutrino::kinematics_description_c_theta_d = "\
    The cosine of the D's helicity angle theta_d in the D-pi rest frame.";

    const std::string
    BToVectorLeptonNeutrino::kinematics_description_phi = "\
    The azimuthal angle between the D-pi plane and the l-nubar  plane.";

    const std::set<ReferenceName>
    BToVectorLeptonNeutrino::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BToVectorLeptonNeutrino::begin_options()
    {
        return Implementation<BToVectorLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToVectorLeptonNeutrino::end_options()
    {
        return Implementation<BToVectorLeptonNeutrino>::options.cend();
    }
}
