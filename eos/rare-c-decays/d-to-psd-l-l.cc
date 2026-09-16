/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026    Carolina Bolognani
 * Copyright (c) 2026    Dominik Suelmann
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

#include <eos/maths/integrate-impl.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/rare-c-decays/d-to-psd-l-l-base.hh>
#include <eos/rare-c-decays/d-to-psd-l-l-bght2019.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    using std::abs;
    using std::norm;
    using std::sqrt;

    struct DToPseudoscalarLeptonLepton::AngularCoefficients
    {
            double a_l, b_l, c_l;

            AngularCoefficients() {}

            AngularCoefficients(const std::array<double, 3> & a) :
                a_l(a[0]),
                b_l(a[1]),
                c_l(a[2])
            {
            }
    };

    /*!
     * Implementation for the decay @f$D_{s} \to P \ell^+ \ell^-@f$.
     */
    template <> struct Implementation<DToPseudoscalarLeptonLepton>
    {
            std::shared_ptr<DToPseudoscalarLeptonLepton::AmplitudeGenerator> amplitude_generator;

            std::shared_ptr<Model> model;

            Parameters parameters;

            LeptonFlavorOption opt_l;
            QuarkFlavorOption  opt_q;
            SpecifiedOption    opt_P;

            UsedParameter hbar;
            UsedParameter m_D;
            UsedParameter m_P;
            UsedParameter m_l;
            UsedParameter tau;
            UsedParameter mu;

            static const std::vector<OptionSpecification> options;

            // { q, P } -> { process, D_name, P_name }
            // q: u, d, s: the spectator quark flavor
            // P: K, pi, eta: the type of daughter meson
            // process: string that can be used to obtain the form factor
            // D_name: name of the D meson
            // P_name: name of the daughter meson
            // c_I: isospin factor by which the amplitudes are multiplied
            static const std::map<std::tuple<QuarkFlavor, std::string>, std::tuple<std::string, std::string>> process_map;

            inline std::string
            _D() const
            {
                const QuarkFlavor q = opt_q.value();
                const std::string P = opt_P.value();
                const auto        p = process_map.find(std::make_tuple(q, P));

                if (p == process_map.end())
                {
                    throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + P);
                }

                return std::get<0>(p->second);
            }

            inline std::string
            _P() const
            {
                const QuarkFlavor q = opt_q.value();
                const std::string P = opt_P.value();
                const auto        p = process_map.find(std::make_tuple(q, P));

                if (p == process_map.end())
                {
                    throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + P);
                }

                return std::get<1>(p->second);
            }

            Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
                model(Model::make(o.get("model"_ok, "WET"_ov), p, o)),
                parameters(p),
                opt_l(o, options, "l"_ok),
                opt_q(o, options, "q"_ok),
                opt_P(o, options, "P"_ok),
                hbar(p["QM::hbar"], u),
                m_D(p["mass::" + _D()], u),
                m_P(p["mass::" + _P()], u),
                m_l(p["mass::" + opt_l.str()], u),
                tau(p["life_time::" + _D()], u),
                mu(p["uc::mu"], u)
            {
                Context ctx("When constructing D(s)->Pll observables");

                std::string tag = o.has("tag"_ok) ? o["tag"_ok].str() : "";

                if ("BGHT2019" == tag)
                {
                    amplitude_generator.reset(new DToPseudoscalarLeptonLeptonAmplitudes<tag::BGHT2019>(p, o));
                }
                else
                {
                    throw InternalError("DToPseudoscalarLeptonLepton: Unknown tag or no valid tag specified (tag = '" + tag + "')!");
                }

                u.uses(*amplitude_generator);
            }

            ~Implementation() {}

            inline std::array<double, 3>
            angular_coefficients_array(const DToPseudoscalarLeptonLepton::Amplitudes & A, const double & q2) const
            {
                // cf. [BHP:2007A], Eq. (4.2) - (4.4)
                std::array<double, 3> result;

                const double beta_l   = amplitude_generator->beta_l(q2);
                const double lambda_D = amplitude_generator->lambda(q2);

                // a_l
                result[0] = amplitude_generator->normalisation(q2)
                            * (q2 * (power_of<2>(beta_l) * norm(A.F_S) + norm(A.F_P)) + 0.25 * lambda_D * (norm(A.F_A) + norm(A.F_V))
                               + 2.0 * m_l * (m_D() * m_D() - m_P() * m_P() + q2) * std::real(A.F_P * std::conj(A.F_A)) + 4.0 * m_l * m_l * m_D() * m_D() * norm(A.F_A));

                // b_l
                result[1] =
                        2.0 * amplitude_generator->normalisation(q2)
                        * (q2 * (power_of<2>(beta_l) * std::real(A.F_S * std::conj(A.F_T)) + std::real(A.F_P * std::conj(A.F_T5)))
                           + m_l * (sqrt(lambda_D) * beta_l * std::real(A.F_S * std::conj(A.F_V)) + (m_D() * m_D() - m_P() * m_P() + q2) * std::real(A.F_T5 * std::conj(A.F_A))));

                // c_l
                result[2] = amplitude_generator->normalisation(q2)
                            * (q2 * (power_of<2>(beta_l) * norm(A.F_T) + norm(A.F_T5)) - 0.25 * lambda_D * power_of<2>(beta_l) * (norm(A.F_A) + norm(A.F_V))
                               + 2.0 * m_l * sqrt(lambda_D) * beta_l * std::real(A.F_T * std::conj(A.F_V)));

                return result;
            }

            inline std::array<double, 3>
            differential_angular_coefficients_array(const double & q2) const
            {
                return angular_coefficients_array(amplitude_generator->amplitudes(q2), q2);
            }

            inline DToPseudoscalarLeptonLepton::AngularCoefficients
            differential_angular_coefficients(const double & q2) const
            {
                return DToPseudoscalarLeptonLepton::AngularCoefficients(differential_angular_coefficients_array(q2));
            }

            // cf. [BHP:2007A], Eq. (4.8)
            inline double
            unnormalized_decay_width(const DToPseudoscalarLeptonLepton::AngularCoefficients & a) const
            {
                return 2.0 * (a.a_l + a.c_l / 3.0);
            }

            inline double
            differential_branching_ratio(const DToPseudoscalarLeptonLepton::AngularCoefficients & a) const
            {
                return unnormalized_decay_width(a) * tau() / hbar();
            }

            // cf. [BHP:2007A], Eq. (4.9)
            inline double
            differential_flat_term_numerator(const DToPseudoscalarLeptonLepton::AngularCoefficients & a) const
            {
                return 2.0 * (a.a_l + a.c_l);
            }

            inline double
            differential_forward_backward_asymmetry_numerator(const DToPseudoscalarLeptonLepton::AngularCoefficients & a) const
            {
                return a.b_l;
            }

            DToPseudoscalarLeptonLepton::AngularCoefficients
            integrated_angular_coefficients(const double & q2_min, const double & q2_max) const
            {
                std::function<std::array<double, 3>(const double &)> integrand =
                        std::bind(&Implementation<DToPseudoscalarLeptonLepton>::differential_angular_coefficients_array, this, std::placeholders::_1);
                std::array<double, 3> integrated_angular_coefficients_array = integrate<1, 3>(integrand, q2_min, q2_max, cubature::Config().epsrel(1e-5));

                return DToPseudoscalarLeptonLepton::AngularCoefficients(integrated_angular_coefficients_array);
            }
    };

    DToPseudoscalarLeptonLepton::DToPseudoscalarLeptonLepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<DToPseudoscalarLeptonLepton>(new Implementation<DToPseudoscalarLeptonLepton>(parameters, options, *this))
    {
    }

    DToPseudoscalarLeptonLepton::~DToPseudoscalarLeptonLepton() {}

    const std::map<std::tuple<QuarkFlavor, std::string>, std::tuple<std::string, std::string>> Implementation<DToPseudoscalarLeptonLepton>::Implementation::process_map{
        {     { QuarkFlavor::up, "pi" }, { "D_u", "pi^0" } },
        {   { QuarkFlavor::down, "pi" }, { "D_d", "pi^+" } },
        { { QuarkFlavor::strange, "K" },  { "D_s", "K_u" } },
    };

    const std::vector<OptionSpecification> Implementation<DToPseudoscalarLeptonLepton>::options{
        Model::option_specification(),
        { "l"_ok, { "e"_ov, "mu"_ov, "tau"_ov }, "mu"_ov },
        { "q"_ok,    { "u"_ov, "d"_ov, "s"_ov },  "d"_ov },
        { "P"_ok,           { "pi"_ov, "K"_ov }, "pi"_ov },
    };

    double
    DToPseudoscalarLeptonLepton::differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_branching_ratio(_imp->differential_angular_coefficients(q2));
    }

    double
    DToPseudoscalarLeptonLepton::differential_flat_term(const double & q2) const
    {
        AngularCoefficients a = _imp->differential_angular_coefficients(q2);

        return _imp->differential_flat_term_numerator(a) / _imp->unnormalized_decay_width(a);
    }

    double
    DToPseudoscalarLeptonLepton::differential_forward_backward_asymmetry(const double & q2) const
    {
        AngularCoefficients a = _imp->differential_angular_coefficients(q2);

        return _imp->differential_forward_backward_asymmetry_numerator(a) / _imp->unnormalized_decay_width(a);
    }

    double
    DToPseudoscalarLeptonLepton::double_differential_decay_width(const double & q2, const double & c_theta_l_LHCb) const
    {
        AngularCoefficients a = _imp->differential_angular_coefficients(q2);

        // using the angular convention of the LHCb experiment
        const double c_theta_l = -c_theta_l_LHCb;

        // cf. [BHP:2007A], Eq. (4.1)
        return a.a_l + a.b_l * c_theta_l + a.c_l * c_theta_l * c_theta_l;
    }

    // Integrated Observables
    double
    DToPseudoscalarLeptonLepton::integrated_decay_width(const double & q2_min, const double & q2_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(q2_min, q2_max);

        return _imp->unnormalized_decay_width(a);
    }

    double
    DToPseudoscalarLeptonLepton::integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(q2_min, q2_max);

        return _imp->differential_branching_ratio(a);
    }

    double
    DToPseudoscalarLeptonLepton::integrated_flat_term_numerator(const double & q2_min, const double & q2_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(q2_min, q2_max);

        return _imp->differential_flat_term_numerator(a);
    }

    double
    DToPseudoscalarLeptonLepton::integrated_forward_backward_asymmetry_numerator(const double & q2_min, const double & q2_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(q2_min, q2_max);

        return _imp->differential_forward_backward_asymmetry_numerator(a);
    }

    const std::string DToPseudoscalarLeptonLepton::description = "\
The decay D->Psd l^+ l^-, with l=e,mu,tau a charged lepton.";

    const std::string DToPseudoscalarLeptonLepton::kinematics_description_q2 = "\
The invariant mass of the charged lepton pair in GeV^2.";

    const std::string DToPseudoscalarLeptonLepton::kinematics_description_c_theta_l = "\
The cosine of the lepton's helicity angle theta_l in the l^+l^- rest frame using the LHCb convention.";

    /*
     * For diagnostic purposes only!
     */
    DToPseudoscalarLeptonLepton::Amplitudes
    DToPseudoscalarLeptonLepton::amplitudes(const double & q2) const
    {
        return _imp->amplitude_generator->amplitudes(q2);
    }

    std::array<double, 3>
    DToPseudoscalarLeptonLepton::angular_coefficients(const double & q2) const
    {
        return _imp->angular_coefficients_array(_imp->amplitude_generator->amplitudes(q2), q2);
    }

    const std::set<ReferenceName> DToPseudoscalarLeptonLepton::references{};

    std::vector<OptionSpecification>::const_iterator
    DToPseudoscalarLeptonLepton::begin_options()
    {
        return Implementation<DToPseudoscalarLeptonLepton>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    DToPseudoscalarLeptonLepton::end_options()
    {
        return Implementation<DToPseudoscalarLeptonLepton>::options.cend();
    }
} // namespace eos
