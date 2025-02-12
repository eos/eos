/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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

#include <eos/b-decays/b-to-psd-psd-l-nu-impl.hh>
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
#include <eos/utils/quantum-numbers.hh>

#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <iostream>

namespace eos
{
    using std::norm;

    template <> struct Implementation<BToPPLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        Parameters parameters;

        QuarkFlavorOption opt_U;
        QuarkFlavorOption opt_q;
        SwitchOption opt_I1, opt_I2, opt_C;

        UsedParameter hbar;

        UsedParameter tau_B;

        UsedParameter g_fermi;

        LeptonFlavorOption opt_l;

        UsedParameter m_l;

        UsedParameter m_B;

        UsedParameter m_P1, m_P2;

        const IsospinRepresentation Ip1, Ip2;

        BooleanOption opt_cp_conjugate;

        bool cp_conjugate;

        UsedParameter mu;

        std::function<complex<double> ()> v_Ub;

        IntegerOption opt_int_points;

        std::shared_ptr<FormFactors<PToPP>> form_factors;

        static const std::vector<OptionSpecification> options;

        // { U, q, I1, I2, C } -> { process, scattering, m_B, m_P1, m_P2, Ip1, Ip2, c_I1, c_I2 }
        static const std::map<std::tuple<QuarkFlavor, QuarkFlavor, std::string, std::string, std::string>,
                              std::tuple<std::string, std::string, std::string, std::string, std::string, IsospinRepresentation, IsospinRepresentation>> process_map;

        inline std::string _process() const
        {
            const QuarkFlavor U = opt_U.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I1 = opt_I1.value();
            const std::string I2 = opt_I2.value();
            const std::string C = opt_C.value();
            const auto p = process_map.find(std::make_tuple(U, q, I1, I2, C));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of U=" + opt_U.str() + ", q=" + opt_q.str() + ", I1=" + I1 + ", I2=" + I2 + ", C=" + C);

            return std::get<0>(p->second);
        }

        inline std::string _scattering_amps() const
        {
            const QuarkFlavor U = opt_U.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I1 = opt_I1.value();
            const std::string I2 = opt_I2.value();
            const std::string C = opt_C.value();
            const auto p = process_map.find(std::make_tuple(U, q, I1, I2, C));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of U=" + opt_U.str() + ", q=" + opt_q.str() + ", I1=" + I1 + ", I2=" + I2 + ", C=" + C);

            return std::get<1>(p->second);
        }

        inline std::string _m_B() const
        {
            const QuarkFlavor U = opt_U.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I1 = opt_I1.value();
            const std::string I2 = opt_I2.value();
            const std::string C = opt_C.value();
            const auto p = process_map.find(std::make_tuple(U, q, I1, I2, C));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of U=" + opt_U.str() + ", q=" + opt_q.str() + ", I1=" + I1 + ", I2=" + I2 + ", C=" + C);

            return std::get<2>(p->second);
        }

        inline std::string _m_P1() const
        {
            const QuarkFlavor U = opt_U.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I1 = opt_I1.value();
            const std::string I2 = opt_I2.value();
            const std::string C = opt_C.value();
            const auto p = process_map.find(std::make_tuple(U, q, I1, I2, C));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of U=" + opt_U.str() + ", q=" + opt_q.str() + ", I1=" + I1 + ", I2=" + I2 + ", C=" + C);

            return std::get<3>(p->second);
        }

        inline std::string _m_P2() const
        {
            const QuarkFlavor U = opt_U.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I1 = opt_I1.value();
            const std::string I2 = opt_I2.value();
            const std::string C = opt_C.value();
            const auto p = process_map.find(std::make_tuple(U, q, I1, I2, C));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of U=" + opt_U.str() + ", q=" + opt_q.str() + ", I1=" + I1 + ", I2=" + I2 + ", C=" + C);

            return std::get<4>(p->second);
        }

        inline IsospinRepresentation _isospin_label_1() const
        {
            const QuarkFlavor U = opt_U.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I1 = opt_I1.value();
            const std::string I2 = opt_I2.value();
            const std::string C = opt_C.value();
            const auto p = process_map.find(std::make_tuple(U, q, I1, I2, C));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of U=" + opt_U.str() + ", q=" + opt_q.str() + ", I1=" + I1 + ", I2=" + I2 + ", C=" + C);

            return std::get<5>(p->second);
        }

        inline IsospinRepresentation _isospin_label_2() const
        {
            const QuarkFlavor U = opt_U.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I1 = opt_I1.value();
            const std::string I2 = opt_I2.value();
            const std::string C = opt_C.value();
            const auto p = process_map.find(std::make_tuple(U, q, I1, I2, C));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of U=" + opt_U.str() + ", q=" + opt_q.str() + ", I1=" + I1 + ", I2=" + I2 + ", C=" + C);

            return std::get<6>(p->second);
        }

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            parameters(p),
            opt_U(o, options, "U"_ok),
            opt_q(o, options, "q"_ok),
            opt_I1(o, "I1"_ok, { "1", "0", "1/2" }),
            opt_I2(o, "I2"_ok, { "1", "0", "1/2" }),
            opt_C(o, "C"_ok, { "+-", "+0", "00" }),
            hbar(p["QM::hbar"], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            g_fermi(p["WET::G_Fermi"], u),
            opt_l(o, options, "l"_ok),
            m_l(p["mass::" + opt_l.str()], u),
            m_B(p["mass::" + _m_B()], u),
            m_P1(p["mass::" + _m_P1()], u),
            m_P2(p["mass::" + _m_P2()], u),
            Ip1(_isospin_label_1()),
            Ip2(_isospin_label_2()),
            opt_cp_conjugate(o, options, "cp-conjugate"_ok),
            cp_conjugate(opt_cp_conjugate.value()),
            mu(p[opt_U.str() + "b" + opt_l.str() + "nu" + opt_l.str() + "::mu"], u),
            opt_int_points(o, options, "integration-points"_ok),
            form_factors(FormFactorFactory<PToPP>::create(_process() + "::" + o.get("form-factors"_ok, "HKvT2025"), p, o))
        {
            Context ctx("When constructing B->PPlnu observable");

            using std::placeholders::_1;
            using std::placeholders::_2;

            switch (opt_U.value())
            {
                case QuarkFlavor::up:
                    v_Ub      = std::bind(&ModelComponent<components::CKM>::ckm_ub, model.get());
                    break;

                default:
                    throw InternalError("Unexpected quark flavor: '" + opt_U.str() + "'");
            }

            u.uses(*form_factors);
            u.uses(*model);
        }

    b_to_psd_psd_l_nu::Amplitudes amplitudes(const double & q2, const double & k2) const
    {
            b_to_psd_psd_l_nu::Amplitudes result;

            // meson & lepton masses
            const double m_l = this->m_l();
            const double m_B = this->m_B();
            const double m_P1 = this->m_P1();
            const double m_P2 = this->m_P2();
            // transversity amplitudes A's. cf. [DSD2014], p.17
            if ((q2 >= power_of<2>(m_l)) && (q2 <= power_of<2>(m_B - m_P1 - m_P2)) && (k2 >= power_of<2>(m_P1 + m_P2)) && (k2 <= power_of<2>(m_B - std::sqrt(q2)))) {
                double lamq3 = (q2 - power_of<2>(m_B + std::sqrt(k2))) * (q2 - power_of<2>(m_B - std::sqrt(k2)));
                double lams12 = (k2 - power_of<2>(m_P1 + m_P2)) * (k2 - power_of<2>(m_P1 - m_P2));
                result.f_perp = form_factors->f_perp(q2, k2);
                result.f_para = form_factors->f_para(q2, k2);
                result.f_long = form_factors->f_long(q2, k2);
                result.f_time = form_factors->f_time(q2, k2);
                result.q2 = q2;
                result.beta_l = (m_l > 0.0) ? (1.0 - m_l * m_l / q2) : 1.0;
                result.beta_pi = std::sqrt(lams12) / k2;
                result.pref = std::norm(v_Ub()) * power_of<2>(g_fermi()) / power_of<3>(m_B) * q2 * power_of<2>(result.beta_l) * result.beta_pi * std::sqrt(lamq3) / power_of<5>(4.0 * M_PI) / 4.0;
            }
            else
            {
                result.f_perp = { 0.0, 0.0, 0.0, 0.0 };
                result.f_para = { 0.0, 0.0, 0.0, 0.0 };
                result.f_long = { 0.0, 0.0, 0.0, 0.0 };
                result.f_time = { 0.0, 0.0, 0.0, 0.0 };

                result.q2 = 0.0;
                result.beta_l = 0.0;
                result.beta_pi = 0.0;
                result.pref = 0.0;
            }

            return result;

        }

        std::array<std::array<double, 5>, 9> _differential_angular_observables(const double & q2, const double & k2) const
        {
            return b_to_psd_psd_l_nu::AngularObservables(this->amplitudes(q2, k2))._M;
        }

        inline b_to_psd_psd_l_nu::AngularObservables differential_angular_observables(const double & q2, const double & k2) const
        {
            return b_to_psd_psd_l_nu::AngularObservables{ _differential_angular_observables(q2, k2) };
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToPPLeptonNeutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToPP>::option_specification(),
        { "cp-conjugate"_ok, { "true", "false" },  "false" },
        { "l"_ok,            { "e", "mu", "tau" }, "mu"    },
        { "U"_ok,            { "c", "u" },         "c"     },
        { "q"_ok,            { "u", "d" },    "d"     },
        { "I1"_ok,           { "1", "0", "1/2" },  "1"     },
        { "I2"_ok,           { "1", "0", "1/2" },  "1"     },
        { "C"_ok,            { "+-", "00", "+0" }, "+-"   },
        { "integration-points"_ok, {"256", "512", "1024", "2048", "4096", "8192", "16384"}, "4096" }
    };

    const std::map<std::tuple<QuarkFlavor, QuarkFlavor, std::string, std::string, std::string>,
                   std::tuple<std::string, std::string, std::string, std::string, std::string, IsospinRepresentation, IsospinRepresentation>>
    Implementation<BToPPLeptonNeutrino>::Implementation::process_map
    {
        // Here we already squared them....
        { { QuarkFlavor::up,    QuarkFlavor::up,    "1",    "1",    "00"   },
          { "B->pipi",  "pipi->pipi", "B_u",  "pi^0", "pi^0", IsospinRepresentation::zero, IsospinRepresentation::one } },
        { { QuarkFlavor::up,    QuarkFlavor::up,    "1",    "1",    "+-"   },
          { "B->pipi",  "pipi->pipi", "B_u",  "pi^+", "pi^+", IsospinRepresentation::zero, IsospinRepresentation::one } },
        { { QuarkFlavor::down,  QuarkFlavor::up,    "1",    "1",    "+0"   },
          { "B->pipi",  "pipi->pipi", "B_d",  "pi^+", "pi^0", IsospinRepresentation::zero, IsospinRepresentation::one } },
    };

    BToPPLeptonNeutrino::BToPPLeptonNeutrino(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BToPPLeptonNeutrino>(new Implementation<BToPPLeptonNeutrino>(p, o, *this))
    {
    }

    BToPPLeptonNeutrino::~BToPPLeptonNeutrino()
    {
    }

    double BToPPLeptonNeutrino::double_differential_decay_width(const double & q2, const double & k2) const
    {
        return _imp->differential_angular_observables(q2, k2).double_differential_decay_width();
    }

    double BToPPLeptonNeutrino::double_differential_mesonic_afb(const double & q2, const double & k2) const
    {
        return _imp->differential_angular_observables(q2, k2).double_differential_mesonic_afb();
    }

    double BToPPLeptonNeutrino::double_differential_branching_ratio(const double & q2, const double & k2) const
    {
        return _imp->differential_angular_observables(q2, k2).double_differential_decay_width() * _imp->tau_B / _imp->hbar;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio(const double & q2_min, const double & q2_max, const double & sqrt_k2_min, const double & sqrt_k2_max) const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobian
            return 2.0 * x[1] * this->double_differential_branching_ratio(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ q2_min, sqrt_k2_min };
        std::array<double, 2> x_max{ q2_max, sqrt_k2_max };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_mesonic_afb(const double & q2_min, const double & q2_max, const double & sqrt_k2_min, const double & sqrt_k2_max) const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobians
            return 2.0 * x[1] * this->double_differential_mesonic_afb(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ q2_min, sqrt_k2_min };
        std::array<double, 2> x_max{ q2_max, sqrt_k2_max };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::fully_integrated_branching_ratio() const
    {
        return integrated_branching_ratio(1e-4, power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), _imp->m_P1 + _imp->m_P2, _imp->m_B);
    }

    double BToPPLeptonNeutrino::fully_integrated_mesonic_afb() const
    {
        return integrated_mesonic_afb(1e-4, power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), _imp->m_P1 + _imp->m_P2, _imp->m_B);
    }

    double BToPPLeptonNeutrino::q2_integrated_branching_ratio(const double & sqrt_k2_min, const double & sqrt_k2_max) const
    {
        return integrated_branching_ratio(1e-4, power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), sqrt_k2_min, sqrt_k2_max);
    }

    double BToPPLeptonNeutrino::q2_integrated_mesonic_afb(const double & sqrt_k2_min, const double & sqrt_k2_max) const
    {
        return integrated_mesonic_afb(1e-4, power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), sqrt_k2_min, sqrt_k2_max);
    }

    double BToPPLeptonNeutrino::s_integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        return integrated_branching_ratio(q2_min, q2_max, _imp->m_P1 + _imp->m_P2, _imp->m_B);
    }

    double BToPPLeptonNeutrino::s_integrated_mesonic_afb(const double & q2_min, const double & q2_max) const
    {
        return integrated_mesonic_afb(q2_min, q2_max, _imp->m_P1 + _imp->m_P2, _imp->m_B);
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio_q2(const double & q2) const
    {
        double k2_min = power_of<2>(_imp->m_P1 + _imp->m_P2);
        double k2_max = power_of<2>(_imp->m_B - std::sqrt(q2));

        std::function<double(const double &)> integrand = [&] (const double & k2)
        {
            return this->double_differential_branching_ratio(q2, k2);
        };

        double res = integrate1D(integrand, _imp->opt_int_points.value(), k2_min, k2_max);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_mesonic_afb_q2(const double & q2) const
    {
        double k2_min = power_of<2>(_imp->m_P1 + _imp->m_P2);
        double k2_max = power_of<2>(_imp->m_B - std::sqrt(q2));

        std::function<double(const double &)> integrand = [&] (const double & k2)
        {
            return this->double_differential_mesonic_afb(q2, k2);
        };

        double res = integrate1D(integrand, _imp->opt_int_points.value(), k2_min, k2_max);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio_sqrt_k2(const double & sqrt_k2) const
    {
        double q2_min = 1e-4;
        double q2_max = power_of<2>(_imp->m_B - sqrt_k2);

        std::function<double(const double &)> integrand = [&] (const double & q2)
        {
            return 2 * sqrt_k2 * this->double_differential_branching_ratio(q2, sqrt_k2 * sqrt_k2);
        };

        double res = integrate1D(integrand, _imp->opt_int_points.value(), q2_min, q2_max);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_mesonic_afb_sqrt_k2(const double & sqrt_k2) const
    {
        double q2_min = 1e-4;
        double q2_max = power_of<2>(_imp->m_B - sqrt_k2);

        std::function<double(const double &)> integrand = [&] (const double & q2)
        {
            return 2 * sqrt_k2 * this->double_differential_mesonic_afb(q2, sqrt_k2 * sqrt_k2);
        };

        double res = integrate1D(integrand, _imp->opt_int_points.value(), q2_min, q2_max);
        return res;
    }

    const std::string
    BToPPLeptonNeutrino::description = "\
    The decay B->D^* l nu, where l=e,mu,tau is a lepton.";

    const std::string
    BToPPLeptonNeutrino::kinematics_description_q2 = "\
    The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    BToPPLeptonNeutrino::kinematics_description_k2 = "\
    The invariant mass of the p-p pair in GeV^2.";

    const std::string
    BToPPLeptonNeutrino::kinematics_description_c_theta_l = "\
    The cosine of the charged lepton's helicity angle theta_l in the l-nubar rest frame.";

    const std::string
    BToPPLeptonNeutrino::kinematics_description_c_theta_nu = "\
    The cosine of the D's helicity angle theta_d in the D-pi rest frame.";

    const std::string
    BToPPLeptonNeutrino::kinematics_description_phi = "\
    The azimuthal angle between the D-pi plane and the l-nubar  plane.";

    const std::set<ReferenceName>
    BToPPLeptonNeutrino::references
    {
        "HKvT:2025A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BToPPLeptonNeutrino::begin_options()
    {
        return Implementation<BToPPLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToPPLeptonNeutrino::end_options()
    {
        return Implementation<BToPPLeptonNeutrino>::options.cend();
    }
}
