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

#include <eos/b-decays/b-to-psd-psd-l-nu.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/scattering/scattering-amplitudes.hh>
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

        const std::array<double, 3> isospin_factor_1, isospin_factor_2;

        const IsospinRepresentation Ip1, Ip2;

        BooleanOption opt_cp_conjugate;

        bool cp_conjugate;

        UsedParameter mu;

        std::function<complex<double> ()> v_Ub;

        SwitchOption opt_int_points;

        int int_points;

        std::shared_ptr<FormFactors<PToPP>> form_factors;

        std::shared_ptr<ScatteringAmplitudes<PPToPP>> scattering_amplitudes;

        static const std::vector<OptionSpecification> options;

        // { U, q, I1, I2, C } -> { process, scattering, m_B, m_P1, m_P2, Ip1, Ip2, c_I1, c_I2 }
        static const std::map<std::tuple<QuarkFlavor, QuarkFlavor, std::string, std::string, std::string>,
                              std::tuple<std::string, std::string, std::string, std::string, std::string, IsospinRepresentation, IsospinRepresentation, std::array<double, 3>, std::array<double, 3>>> process_map;

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

        inline std::array<double, 3> _isospin_factor_1() const
        {
            const QuarkFlavor U = opt_U.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I1 = opt_I1.value();
            const std::string I2 = opt_I2.value();
            const std::string C = opt_C.value();
            const auto p = process_map.find(std::make_tuple(U, q, I1, I2, C));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of U=" + opt_U.str() + ", q=" + opt_q.str() + ", I1=" + I1 + ", I2=" + I2 + ", C=" + C);

            return std::get<7>(p->second);
        }

        inline std::array<double, 3> _isospin_factor_2() const
        {
            const QuarkFlavor U = opt_U.value();
            const QuarkFlavor q = opt_q.value();
            const std::string I1 = opt_I1.value();
            const std::string I2 = opt_I2.value();
            const std::string C = opt_C.value();
            const auto p = process_map.find(std::make_tuple(U, q, I1, I2, C));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of U=" + opt_U.str() + ", q=" + opt_q.str() + ", I1=" + I1 + ", I2=" + I2 + ", C=" + C);

            return std::get<8>(p->second);
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
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            opt_U(o, options, "U"),
            opt_q(o, options, "q"),
            opt_I1(o, "I1", { "1", "0", "1/2" }),
            opt_I2(o, "I2", { "1", "0", "1/2" }),
            opt_C(o, "C", { "+-", "+0", "00" }),
            hbar(p["QM::hbar"], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            g_fermi(p["WET::G_Fermi"], u),
            opt_l(o, options, "l"),
            m_l(p["mass::" + opt_l.str()], u),
            m_B(p["mass::" + _m_B()], u),
            m_P1(p["mass::" + _m_P1()], u),
            m_P2(p["mass::" + _m_P2()], u),
            isospin_factor_1(_isospin_factor_1()),
            isospin_factor_2(_isospin_factor_2()),
            Ip1(_isospin_label_1()),
            Ip2(_isospin_label_2()),
            opt_cp_conjugate(o, options, "cp-conjugate"),
            cp_conjugate(opt_cp_conjugate.value()),
            mu(p[opt_U.str() + "b" + opt_l.str() + "nu" + opt_l.str() + "::mu"], u),
            opt_int_points(o, "integration-points", {"256", "512", "1024", "2048", "4096", "8192", "16384"}, "4096"),
            int_points(destringify<int>(opt_int_points.value())),
            form_factors(FormFactorFactory<PToPP>::create(_process() + "::" + o.get("form-factors", "HKvT2025"), p, o)),
            scattering_amplitudes(ScatteringAmplitudeFactory<PPToPP>::create(_scattering_amps() + "::" + o.get("scattering-amplitudes", "HKvT2025"), p, o))
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
            u.uses(*scattering_amplitudes);
            u.uses(*model);
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToPPLeptonNeutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToPP>::option_specification(),
        { "cp-conjugate", { "true", "false" },  "false" },
        { "l",            { "e", "mu", "tau" }, "mu"    },
        { "U",            { "c", "u" },         "c"     },
        { "q",            { "u", "d", "s" },    "d"     },
        { "I1",           { "1", "0", "1/2" },  "1"     },
        { "I2",           { "1", "0", "1/2" },  "1"     },
        { "C",            { "+-", "00", "+0" }, "+-"   },
    };

    const std::map<std::tuple<QuarkFlavor, QuarkFlavor, std::string, std::string, std::string>,
                   std::tuple<std::string, std::string, std::string, std::string, std::string, IsospinRepresentation, IsospinRepresentation, std::array<double, 3>, std::array<double, 3>>>
    Implementation<BToPPLeptonNeutrino>::Implementation::process_map
    {
        // Here we already squared them....
        { { QuarkFlavor::up,    QuarkFlavor::up,    "1",    "1",    "00"   },
          { "B->pipi",  "pipi->pipi", "B_u",  "pi^0", "pi^0", IsospinRepresentation::zero, IsospinRepresentation::one,
            std::array<double, 3>{1.0 / 6.0,0.0,1.0 / 6.0}, std::array<double, 3>{0.0,0.0,0.0} } },
        { { QuarkFlavor::up,    QuarkFlavor::up,    "1",    "1",    "+-"   },
          { "B->pipi",  "pipi->pipi", "B_u",  "pi^+", "pi^+", IsospinRepresentation::zero, IsospinRepresentation::one,
            std::array<double, 3>{1.0 / 6.0,0.0,1.0 / 6.0}, std::array<double, 3>{0.0,1.0 / 4.0,0.0} } },
        { { QuarkFlavor::down,  QuarkFlavor::up,    "1",    "1",    "+0"   },
          { "B->pipi",  "pipi->pipi", "B_d",  "pi^+", "pi^0", IsospinRepresentation::zero, IsospinRepresentation::one,
            std::array<double, 3>{0.0,0.0,0.0}, std::array<double, 3>{0.0,1.0 / 2.0,0.0} } },
    };

    BToPPLeptonNeutrino::BToPPLeptonNeutrino(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BToPPLeptonNeutrino>(new Implementation<BToPPLeptonNeutrino>(p, o, *this))
    {
    }

    BToPPLeptonNeutrino::~BToPPLeptonNeutrino()
    {
    }

    double
    BToPPLeptonNeutrino::saturation_0_m() const
    {
        std::function<double(const double &)> integrand = [&] (const double & t)
        {
            double contrib = 0.0;
            for (unsigned l = 0 ; l < _imp->form_factors->numWaves() ; l++)
            {
                contrib += _imp->form_factors->unitarity_integrand_0m(1.0 / t, l, false) * std::norm(_imp->scattering_amplitudes->omnes_factor(1.0 / t, l, _imp->Ip1))
                          * std::norm(_imp->scattering_amplitudes->isospin_breaking(1.0 / t, l, _imp->Ip1));
                contrib += _imp->form_factors->unitarity_integrand_0m(1.0 / t, l, true) * std::norm(_imp->scattering_amplitudes->omnes_factor(1.0 / t, l, _imp->Ip2))
                          * std::norm(_imp->scattering_amplitudes->isospin_breaking(1.0 / t, l, _imp->Ip2));
            }
            return contrib / power_of<2>(t);
        };

        double res = integrate1D(integrand, _imp->int_points, 1e-5, 1.0 / power_of<2>(_imp->m_P1 + _imp->m_P2));

        return res;
    }

    double
    BToPPLeptonNeutrino::saturation_1_m() const
    {
        std::function<double(const double &)> integrand = [&] (const double & t)
        {
            double contrib = 0.0;
            for (unsigned l = 0 ; l < _imp->form_factors->numWaves() ; l++)
            {
                contrib += _imp->form_factors->unitarity_integrand_1m(1.0 / t, l, false) * std::norm(_imp->scattering_amplitudes->omnes_factor(1.0 / t, l, _imp->Ip1))
                          * std::norm(_imp->scattering_amplitudes->isospin_breaking(1.0 / t, l, _imp->Ip1));
                contrib += _imp->form_factors->unitarity_integrand_1m(1.0 / t, l, true) * std::norm(_imp->scattering_amplitudes->omnes_factor(1.0 / t, l, _imp->Ip2))
                          * std::norm(_imp->scattering_amplitudes->isospin_breaking(1.0 / t, l, _imp->Ip2));
            }
            return contrib / power_of<2>(t);
        };

        double res = integrate1D(integrand, _imp->int_points, 1e-5, 1.0 / power_of<2>(_imp->m_P1 + _imp->m_P2));

        return res;
    }

    double
    BToPPLeptonNeutrino::saturation_1_p() const
    {
        std::function<double(const double &)> integrand = [&] (const double & t)
        {
            double contrib = 0.0;
            for (unsigned l = 0 ; l < _imp->form_factors->numWaves() ; l++)
            {
                contrib += _imp->form_factors->unitarity_integrand_1p(1.0 / t, l, false) * std::norm(_imp->scattering_amplitudes->omnes_factor(1.0 / t, l, _imp->Ip1))
                          * std::norm(_imp->scattering_amplitudes->isospin_breaking(1.0 / t, l, _imp->Ip1));
                contrib += _imp->form_factors->unitarity_integrand_1p(1.0 / t, l, true) * std::norm(_imp->scattering_amplitudes->omnes_factor(1.0 / t, l, _imp->Ip2))
                          * std::norm(_imp->scattering_amplitudes->isospin_breaking(1.0 / t, l, _imp->Ip2));
            }
            return contrib / power_of<2>(t);
        };

        double res = integrate1D(integrand, _imp->int_points, 1e-5, 1.0 / power_of<2>(_imp->m_P1 + _imp->m_P2));

        return res;
    }

    double BToPPLeptonNeutrino::double_differential_decay_width(const double & q2, const double & s) const
    {
        const double lamq3 = (q2 - power_of<2>(_imp->m_B + std::sqrt(s))) * (q2 - power_of<2>(_imp->m_B - std::sqrt(s)));
        const double lams12 = (s - power_of<2>(_imp->m_P1 + _imp->m_P2)) * (s - power_of<2>(_imp->m_P1 - _imp->m_P2));
        if (lamq3 < 0.0 || lams12 < 0.0)
            return 0.0;
        const double kappa = std::sqrt(lams12) * std::sqrt(lamq3) / s;
        const double ml2q2 = power_of<2>(_imp->m_l) / q2;
        const double pref = std::norm(_imp->v_Ub()) * power_of<2>(_imp->g_fermi()) / power_of<3>(_imp->m_B) * q2 * power_of<2>(1.0 - ml2q2) * kappa / power_of<5>(4.0 * M_PI) / 4.0;
        double res1 = 0.0;
        double res2 = 0.0;

        const double kin1 = ml2q2 / q2;
        const double kin2 = (2.0 + ml2q2) / 12.0 / q2;
        const double kin3 = lams12 / s * (2.0 + ml2q2) / 3.0;

        // S-wave
        if (_imp->isospin_factor_1[0] != 0.0)
        {
            res1 += _imp->isospin_factor_1[0] * kin1 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, 0, _imp->Ip1)) * std::norm(_imp->form_factors->a_t(q2, s, 0, false));
            res1 += _imp->isospin_factor_1[0] * lamq3 * kin2 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, 0, _imp->Ip1)) * std::norm(_imp->form_factors->a_0(q2, s, 0, false));
        }

        if (_imp->isospin_factor_2[0] != 0.0)
        {
            res2 += _imp->isospin_factor_2[0] * kin1 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, 0, _imp->Ip2)) * std::norm(_imp->form_factors->a_t(q2, s, 0, true));
            res2 += _imp->isospin_factor_2[0] * lamq3 * kin2 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, 0, _imp->Ip2)) * std::norm(_imp->form_factors->a_0(q2, s, 0, true));
        }

        for (unsigned l = 1 ; l < _imp->form_factors->numWaves() ; l++)
        {
            const double twolplus1 = 1.0 / (2 * l + 1);
            if (_imp->isospin_factor_1[l] != 0.0)
            {
                res1 += _imp->isospin_factor_1[l] * kin1 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip1))
                                                  * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->a_t(q2, s, l, false)) * twolplus1;
                res1 += _imp->isospin_factor_1[l] * kin2 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip1))
                                                  * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->a_0(q2, s, l, false)) * twolplus1;
                res1 += _imp->isospin_factor_1[l] * kin3 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip1))
                                                  * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->a_par(q2, s, l, false)) * l * (l + 1) * twolplus1;
                res1 += _imp->isospin_factor_1[l] * lamq3 * kin3 / 16.0 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip1))
                                                  * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->v_perp(q2, s, l, false)) * l * (l + 1) * twolplus1;
            }
            if (_imp->isospin_factor_2[l] != 0.0)
            {
                res2 += _imp->isospin_factor_2[l] * kin1 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip2))
                                                  * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->a_t(q2, s, l, true)) * twolplus1;
                res2 += _imp->isospin_factor_2[l] * kin2 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip2))
                                                  * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->a_0(q2, s, l, true)) * twolplus1;
                res2 += _imp->isospin_factor_2[l] * kin3 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip2))
                                                  * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->a_par(q2, s, l, true)) * l * (l + 1) * twolplus1;
                res2 += _imp->isospin_factor_2[l] * lamq3 * kin3 / 16.0 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip2))
                                                  * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->v_perp(q2, s, l, true)) * l * (l + 1) * twolplus1;
            }
        }

        return 2.0 * pref * (res1 + res2);
    }

    double BToPPLeptonNeutrino::double_differential_decay_width_S(const double & q2, const double & s) const
    {
        const double lamq3 = (q2 - power_of<2>(_imp->m_B + std::sqrt(s))) * (q2 - power_of<2>(_imp->m_B - std::sqrt(s)));
        const double lams12 = (s - power_of<2>(_imp->m_P1 + _imp->m_P2)) * (s - power_of<2>(_imp->m_P1 - _imp->m_P2));
        if (lamq3 < 0.0 || lams12 < 0.0)
            return 0.0;
        const double kappa = std::sqrt(lams12) * std::sqrt(lamq3) / s;
        const double ml2q2 = power_of<2>(_imp->m_l) / q2;
        const double pref = std::norm(_imp->v_Ub()) * power_of<2>(_imp->g_fermi()) / power_of<3>(_imp->m_B) * q2 * power_of<2>(1.0 - ml2q2) * kappa / power_of<5>(4.0 * M_PI) / 4.0;
        double res1 = 0.0;
        double res2 = 0.0;

        const double kin1 = ml2q2 / q2;
        const double kin2 = (2.0 + ml2q2) / 12.0 / q2;

        // S-wave
        if (_imp->isospin_factor_1[0] != 0.0)
        {
            res1 += _imp->isospin_factor_1[0] * kin1 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, 0, _imp->Ip1)) * std::norm(_imp->form_factors->a_t(q2, s, 0, false));
            res1 += _imp->isospin_factor_1[0] * lamq3 * kin2 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, 0, _imp->Ip1)) * std::norm(_imp->form_factors->a_0(q2, s, 0, false));
        }

        if (_imp->isospin_factor_2[0] != 0.0)
        {
            res2 += _imp->isospin_factor_2[0] * kin1 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, 0, _imp->Ip2)) * std::norm(_imp->form_factors->a_t(q2, s, 0, true));
            res2 += _imp->isospin_factor_2[0] * lamq3 * kin2 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, 0, _imp->Ip2)) * std::norm(_imp->form_factors->a_0(q2, s, 0, true));
        }

        return 2.0 * pref * (res1 + res2);
    }

    double BToPPLeptonNeutrino::double_differential_decay_width_P(const double & q2, const double & s) const
    {
        const double lamq3 = (q2 - power_of<2>(_imp->m_B + std::sqrt(s))) * (q2 - power_of<2>(_imp->m_B - std::sqrt(s)));
        const double lams12 = (s - power_of<2>(_imp->m_P1 + _imp->m_P2)) * (s - power_of<2>(_imp->m_P1 - _imp->m_P2));
        if (lamq3 < 0.0 || lams12 < 0.0)
            return 0.0;
        const double kappa = std::sqrt(lams12) * std::sqrt(lamq3) / s;
        const double ml2q2 = power_of<2>(_imp->m_l) / q2;
        const double pref = std::norm(_imp->v_Ub()) * power_of<2>(_imp->g_fermi()) / power_of<3>(_imp->m_B) * q2 * power_of<2>(1.0 - ml2q2) * kappa / power_of<5>(4.0 * M_PI) / 4.0;
        double res1 = 0.0;
        double res2 = 0.0;

        const double kin1 = ml2q2 / q2;
        const double kin2 = (2.0 + ml2q2) / 12.0 / q2;
        const double kin3 = lams12 / s * (2.0 + ml2q2) / 3.0;

        unsigned l = 1;
        const double twolplus1 = 1.0 / (2 * l + 1);
        if (_imp->isospin_factor_1[l] != 0.0)
        {
            res1 += _imp->isospin_factor_1[l] * kin1 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip1))
                                              * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->a_t(q2, s, l, false)) * twolplus1;
            res1 += _imp->isospin_factor_1[l] * kin2 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip1))
                                              * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->a_0(q2, s, l, false)) * twolplus1;
            res1 += _imp->isospin_factor_1[l] * kin3 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip1))
                                              * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->a_par(q2, s, l, false)) * l * (l + 1) * twolplus1;
            res1 += _imp->isospin_factor_1[l] * lamq3 * kin3 / 16.0 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip1))
                                              * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->v_perp(q2, s, l, false)) * l * (l + 1) * twolplus1;
        }
        if (_imp->isospin_factor_2[l] != 0.0)
        {
            res2 += _imp->isospin_factor_2[l] * kin1 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip2))
                                              * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->a_t(q2, s, l, true)) * twolplus1;
            res2 += _imp->isospin_factor_2[l] * kin2 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip2))
                                              * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->a_0(q2, s, l, true)) * twolplus1;
            res2 += _imp->isospin_factor_2[l] * kin3 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip2))
                                              * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->a_par(q2, s, l, true)) * l * (l + 1) * twolplus1;
            res2 += _imp->isospin_factor_2[l] * lamq3 * kin3 / 16.0 * std::norm(_imp->scattering_amplitudes->isospin_breaking(s, l, _imp->Ip2))
                                              * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->v_perp(q2, s, l, true)) * l * (l + 1) * twolplus1;
        }

        return 2.0 * pref * (res1 + res2);
    }

    double BToPPLeptonNeutrino::double_differential_decay_width_D(const double & q2, const double & s) const
    {
        const double lamq3 = (q2 - power_of<2>(_imp->m_B + std::sqrt(s))) * (q2 - power_of<2>(_imp->m_B - std::sqrt(s)));
        const double lams12 = (s - power_of<2>(_imp->m_P1 + _imp->m_P2)) * (s - power_of<2>(_imp->m_P1 - _imp->m_P2));
        if (lamq3 < 0.0 || lams12 < 0.0)
            return 0.0;
        const double kappa = std::sqrt(lams12) * std::sqrt(lamq3) / s;
        const double ml2q2 = power_of<2>(_imp->m_l) / q2;
        const double pref = std::norm(_imp->v_Ub()) * power_of<2>(_imp->g_fermi()) / power_of<3>(_imp->m_B) * q2 * power_of<2>(1.0 - ml2q2) * kappa / power_of<5>(4.0 * M_PI) / 4.0;
        double res1 = 0.0;
        double res2 = 0.0;

        const double kin1 = ml2q2 / q2;
        const double kin2 = (2.0 + ml2q2) / 12.0 / q2;
        const double kin3 = lams12 / s * (2.0 + ml2q2) / 3.0;

        unsigned l = 2;
        const double twolplus1 = 1.0 / (2 * l + 1);
        if (_imp->isospin_factor_1[l] != 0.0)
        {
            res1 += _imp->isospin_factor_1[l] * kin1 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->a_t(q2, s, l, false)) * twolplus1;
            res1 += _imp->isospin_factor_1[l] * kin2 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->a_0(q2, s, l, false)) * twolplus1;
            res1 += _imp->isospin_factor_1[l] * kin3 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->a_par(q2, s, l, false)) * l * (l + 1) * twolplus1;
            res1 += _imp->isospin_factor_1[l] * lamq3 * kin3 / 16.0 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip1)) * std::norm(_imp->form_factors->v_perp(q2, s, l, false)) * l * (l + 1) * twolplus1;
        }
        if (_imp->isospin_factor_2[l] != 0.0)
        {
            res2 += _imp->isospin_factor_2[l] * kin1 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->a_t(q2, s, l, true)) * twolplus1;
            res2 += _imp->isospin_factor_2[l] * kin2 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->a_0(q2, s, l, true)) * twolplus1;
            res2 += _imp->isospin_factor_2[l] * kin3 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->a_par(q2, s, l, true)) * l * (l + 1) * twolplus1;
            res2 += _imp->isospin_factor_2[l] * lamq3 * kin3 / 16.0 * std::norm(_imp->scattering_amplitudes->omnes_factor(s, l, _imp->Ip2)) * std::norm(_imp->form_factors->v_perp(q2, s, l, true)) * l * (l + 1) * twolplus1;
        }

        return 2.0 * pref * (res1 + res2);
    }

    double BToPPLeptonNeutrino::double_differential_branching_ratio(const double & q2, const double & s) const
    {
        return double_differential_decay_width(q2, s) * _imp->tau_B / _imp->hbar;
    }

    double BToPPLeptonNeutrino::double_differential_branching_ratio_S(const double & q2, const double & s) const
    {
        return double_differential_decay_width_S(q2, s) * _imp->tau_B / _imp->hbar;
    }

    double BToPPLeptonNeutrino::double_differential_branching_ratio_P(const double & q2, const double & s) const
    {
        return double_differential_decay_width_P(q2, s) * _imp->tau_B / _imp->hbar;
    }

    double BToPPLeptonNeutrino::double_differential_branching_ratio_D(const double & q2, const double & s) const
    {
        return double_differential_decay_width_D(q2, s) * _imp->tau_B / _imp->hbar;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio(const double & q2_min, const double & q2_max, const double & sqrts_min, const double & sqrts_max) const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobian
            return 2.0 * x[1] * this->double_differential_branching_ratio(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ q2_min, sqrts_min };
        std::array<double, 2> x_max{ q2_max, sqrts_max };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::fully_integrated_branching_ratio() const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobian
            return 2.0 * x[1] * this->double_differential_branching_ratio(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ 1e-4, _imp->m_P1 + _imp->m_P2 };
        std::array<double, 2> x_max{ power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), _imp->m_B };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::fully_integrated_branching_ratio_S() const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobian
            return 2.0 * x[1] * this->double_differential_branching_ratio_S(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ 1e-4, _imp->m_P1 + _imp->m_P2 };
        std::array<double, 2> x_max{ power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), _imp->m_B };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::fully_integrated_branching_ratio_P() const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobian
            return 2.0 * x[1] * this->double_differential_branching_ratio_P(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ 1e-4, _imp->m_P1 + _imp->m_P2 };
        std::array<double, 2> x_max{ power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), _imp->m_B };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::fully_integrated_branching_ratio_D() const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobian
            return 2.0 * x[1] * this->double_differential_branching_ratio_D(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ 1e-4, _imp->m_P1 + _imp->m_P2 };
        std::array<double, 2> x_max{ power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), _imp->m_B };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::q2_integrated_branching_ratio(const double & sqrts_min, const double & sqrts_max) const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobian
            return 2.0 * x[1] * this->double_differential_branching_ratio(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ 1e-4, sqrts_min };
        std::array<double, 2> x_max{ power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), sqrts_max };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::q2_integrated_branching_ratio_S(const double & sqrts_min, const double & sqrts_max) const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobian
            return 2.0 * x[1] * this->double_differential_branching_ratio_S(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ 1e-4, sqrts_min };
        std::array<double, 2> x_max{ power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), sqrts_max };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::q2_integrated_branching_ratio_P(const double & sqrts_min, const double & sqrts_max) const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobian
            return 2.0 * x[1] * this->double_differential_branching_ratio_P(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ 1e-4, sqrts_min };
        std::array<double, 2> x_max{ power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), sqrts_max };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::q2_integrated_branching_ratio_D(const double & sqrts_min, const double & sqrts_max) const
    {
        std::function<double(const std::array<double, 2u> &)> integrand = [&] (const std::array<double, 2u> & x)
        {
            // Multiply by s -> sqrt(s) Jacobian
            return 2.0 * x[1] * this->double_differential_branching_ratio_D(x[0], x[1] * x[1]);
        };

        auto config_cubature = cubature::Config().epsrel(5e-3);

        std::array<double, 2> x_min{ 1e-4, sqrts_min };
        std::array<double, 2> x_max{ power_of<2>(_imp->m_B - _imp->m_P1 - _imp->m_P2), sqrts_max };

        double res = integrate(integrand, x_min, x_max, config_cubature);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio_q2(const double & q2) const
    {
        double s_min = power_of<2>(_imp->m_P1 + _imp->m_P2);
        double s_max = power_of<2>(_imp->m_B - std::sqrt(q2));

        std::function<double(const double &)> integrand = [&] (const double & s)
        {
            return this->double_differential_branching_ratio(q2, s);
        };

        double res = integrate1D(integrand, _imp->int_points, s_min, s_max);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio_sqrts(const double & sqrts) const
    {
        double q2_min = 1e-4;
        double q2_max = power_of<2>(_imp->m_B - sqrts);

        std::function<double(const double &)> integrand = [&] (const double & q2)
        {
            return 2 * sqrts * this->double_differential_branching_ratio(q2, sqrts * sqrts);
        };

        double res = integrate1D(integrand, _imp->int_points, q2_min, q2_max);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio_q2_S(const double & q2) const
    {
        double s_min = power_of<2>(_imp->m_P1 + _imp->m_P2);
        double s_max = power_of<2>(_imp->m_B - std::sqrt(q2));

        std::function<double(const double &)> integrand = [&] (const double & s)
        {
            return this->double_differential_branching_ratio_S(q2, s);
        };

        double res = integrate1D(integrand, _imp->int_points, s_min, s_max);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio_sqrts_S(const double & sqrts) const
    {
        double q2_min = 1e-4;
        double q2_max = power_of<2>(_imp->m_B - sqrts);

        std::function<double(const double &)> integrand = [&] (const double & q2)
        {
            return 2 * sqrts * this->double_differential_branching_ratio_S(q2, sqrts * sqrts);
        };

        double res = integrate1D(integrand, _imp->int_points, q2_min, q2_max);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio_q2_P(const double & q2) const
    {
        double s_min = power_of<2>(_imp->m_P1 + _imp->m_P2);
        double s_max = power_of<2>(_imp->m_B - std::sqrt(q2));

        std::function<double(const double &)> integrand = [&] (const double & s)
        {
            return this->double_differential_branching_ratio_P(q2, s);
        };

        double res = integrate1D(integrand, _imp->int_points, s_min, s_max);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio_sqrts_P(const double & sqrts) const
    {
        double q2_min = 1e-4;
        double q2_max = power_of<2>(_imp->m_B - sqrts);

        std::function<double(const double &)> integrand = [&] (const double & q2)
        {
            return 2 * sqrts * this->double_differential_branching_ratio_P(q2, sqrts * sqrts);
        };

        double res = integrate1D(integrand, _imp->int_points, q2_min, q2_max);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio_q2_D(const double & q2) const
    {
        double s_min = power_of<2>(_imp->m_P1 + _imp->m_P2);
        double s_max = power_of<2>(_imp->m_B - std::sqrt(q2));

        std::function<double(const double &)> integrand = [&] (const double & s)
        {
            return this->double_differential_branching_ratio_D(q2, s);
        };

        double res = integrate1D(integrand, _imp->int_points, s_min, s_max);
        return res;
    }

    double BToPPLeptonNeutrino::integrated_branching_ratio_sqrts_D(const double & sqrts) const
    {
        double q2_min = 1e-4;
        double q2_max = power_of<2>(_imp->m_B - sqrts);

        std::function<double(const double &)> integrand = [&] (const double & q2)
        {
            return 2 * sqrts * this->double_differential_branching_ratio_D(q2, sqrts * sqrts);
        };

        double res = integrate1D(integrand, _imp->int_points, q2_min, q2_max);
        return res;
    }

    const std::string
    BToPPLeptonNeutrino::description = "\
    The decay B->D^* l nu, where l=e,mu,tau is a lepton.";

    const std::string
    BToPPLeptonNeutrino::kinematics_description_q2 = "\
    The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    BToPPLeptonNeutrino::kinematics_description_s = "\
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
