/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
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
#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-nu-nu.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-nu-nu-impl.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

using std::norm;
using std::sqrt;

namespace eos
{
    template <>
    struct Implementation<LambdaBToLambdaDineutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

        UsedParameter m_Lambda_b;
        UsedParameter tau_Lambda_b;
        UsedParameter m_Lambda;
        UsedParameter g_fermi;
        UsedParameter alpha_e;
        UsedParameter hbar;
        UsedParameter mu;

        using IntermediateResult = LambdaBToLambdaDineutrino::IntermediateResult;
        IntermediateResult intermediate_result;

        static const std::vector<OptionSpecification> options;

        std::function<complex<double> ()> lambda_t;
        std::function<WilsonCoefficients<wc::SBNuNu> ()> wc;

        BooleanOption opt_cp_conjugate;

        bool cp_conjugate;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            m_Lambda_b(p["mass::Lambda_b"], u),
            tau_Lambda_b(p["life_time::Lambda_b"], u),
            m_Lambda(p["mass::Lambda"], u),
            g_fermi(p["WET::G_Fermi"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            hbar(p["QM::hbar"], u),
            mu(p["sbnunu::mu"], u),
            opt_cp_conjugate(o, options, "cp-conjugate"_ok),
            cp_conjugate(opt_cp_conjugate.value())
        {
            Context ctx("When constructing Lb->Lnunu observables");

            form_factors = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda::" + o.get("form-factors"_ok, "BFvD2014"), p, o);

            lambda_t  = [*this] () { return model->ckm_tb() * std::conj(model->ckm_ts()); };
            wc        = [*this] () { return model->wet_sbnunu(cp_conjugate); };

            u.uses(*form_factors);
            u.uses(*model);
        }

        inline std::array<double, 2> angular_coefficients_array(const double & s) const
        {
            std::array<double, 2> result;

            const auto wc = this->wc();
            const double lambda = eos::lambda(power_of<2>(m_Lambda_b), power_of<2>(m_Lambda), s),
                         sqrt_lambda = std::sqrt(lambda),
                         s_minus = power_of<2>(m_Lambda_b - m_Lambda) - s,
                         s_plus = power_of<2>(m_Lambda_b + m_Lambda) - s;

            std::complex<double> N = g_fermi * alpha_e * lambda_t() * sqrt(s * sqrt_lambda /
                3.0 / 2048.0 / power_of<3>(m_Lambda_b) / power_of<5>(M_PI));

            std::complex<double>
                a_perp_plus =   2 * sqrt(2.0) * N * (wc.cVL() + wc.cVR()) * (- sqrt(2 * s_minus) * form_factors->f_perp_v(s)),
                a_para_plus = - 2 * sqrt(2.0) * N * (wc.cVL() - wc.cVR()) * (- sqrt(2 * s_plus) * form_factors->f_perp_a(s)),
                a_perp_long =   2 * sqrt(2.0) * N * (wc.cVL() + wc.cVR()) * ((m_Lambda_b + m_Lambda) * sqrt(s_minus / s) * form_factors->f_long_v(s)),
                a_para_long = - 2 * sqrt(2.0) * N * (wc.cVL() - wc.cVR()) * ((m_Lambda_b - m_Lambda) * sqrt(s_plus / s) * form_factors->f_long_a(s));

            // K1ss
            result[0] = 0.25 * (norm(a_perp_plus) + norm(a_para_plus) + 2 * norm(a_perp_long) + 2 * norm(a_para_long));
            // K1cc
            result[1] = 0.5 * (norm(a_perp_plus) + norm(a_para_plus));

            return result;
        }

        inline LambdaBToLambdaDineutrino::AngularCoefficients differential_angular_coefficients(const double & s) const
        {
            return LambdaBToLambdaDineutrino::AngularCoefficients(angular_coefficients_array(s));
        }

        LambdaBToLambdaDineutrino::AngularCoefficients integrated_angular_coefficients(const double & s_min, const double & s_max) const
        {
            std::function<std::array<double, 2> (const double &)> integrand =
                    std::bind(&Implementation<LambdaBToLambdaDineutrino>::angular_coefficients_array, this, std::placeholders::_1);
            std::array<double, 2> integrated_angular_coefficients_array = integrate1D(integrand, 64, s_min, s_max);

            return LambdaBToLambdaDineutrino::AngularCoefficients(integrated_angular_coefficients_array);
        }

        const IntermediateResult * prepare(const double & q2_min, const double & q2_max)
        {
            intermediate_result.ac = integrated_angular_coefficients(q2_min, q2_max);
            return &intermediate_result;
        }

        inline double decay_width(const LambdaBToLambdaDineutrino::AngularCoefficients & a_c)
        {
            // assume the production of 3 diagonal neutrino flavors (nu_i nubar_i)
            return 3.0 * (2.0 * a_c.K1ss + a_c.K1cc);
        }
    };

    LambdaBToLambdaDineutrino::LambdaBToLambdaDineutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<LambdaBToLambdaDineutrino>(new Implementation<LambdaBToLambdaDineutrino>(parameters, options, *this))
    {
    }

    LambdaBToLambdaDineutrino::~LambdaBToLambdaDineutrino()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<LambdaBToLambdaDineutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<OneHalfPlusToOneHalfPlus>::option_specification(),
        { "cp-conjugate"_ok, { "true", "false" },  "false" }
    };


    double
    LambdaBToLambdaDineutrino::differential_decay_width(const double & s) const
    {
        return _imp->decay_width(_imp->differential_angular_coefficients(s));
    }

    double
    LambdaBToLambdaDineutrino::differential_branching_ratio(const double & s) const
    {
        return differential_decay_width(s) * _imp->tau_Lambda_b() / _imp->hbar();
    }

    double
    LambdaBToLambdaDineutrino::differential_longitudinal_polarisation(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        // assume the production of 3 diagonal neutrino flavors (nu_i nubar_i)
        return 3.0 * (2.0 * a_c.K1ss - a_c.K1cc) / _imp->decay_width(a_c);
    }

    const LambdaBToLambdaDineutrino::IntermediateResult *
    LambdaBToLambdaDineutrino::prepare(const double & q2_min, const double & q2_max) const
    {
        return _imp->prepare(q2_min, q2_max);
    }

    double
    LambdaBToLambdaDineutrino::integrated_decay_width(const IntermediateResult * ir) const
    {
        return _imp->decay_width(ir->ac);
    }

    double
    LambdaBToLambdaDineutrino::integrated_branching_ratio(const IntermediateResult * ir) const
    {
        return integrated_decay_width(ir) * _imp->tau_Lambda_b() / _imp->hbar();
    }

    double
    LambdaBToLambdaDineutrino::integrated_longitudinal_polarisation(const IntermediateResult * ir) const
    {
        AngularCoefficients a_c = ir->ac;
        // assume the production of 3 diagonal neutrino flavors (nu_i nubar_i)
        return 3.0 * (2.0 * a_c.K1ss + a_c.K1cc) / _imp->decay_width(a_c);
    }


    const std::string
    LambdaBToLambdaDineutrino::description = "\
The decay Lambda_b->Lambda nu nu, assuming left-handed neutrinos and left-handed sb current";

    const std::string
    LambdaBToLambdaDineutrino::kinematics_description_q2 = "\
The invariant mass of the nu-nubar pair in GeV^2.";

    const std::set<ReferenceName>
    LambdaBToLambdaDineutrino::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambdaDineutrino::begin_options()
    {
        return Implementation<LambdaBToLambdaDineutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambdaDineutrino::end_options()
    {
        return Implementation<LambdaBToLambdaDineutrino>::options.cend();
    }
}
