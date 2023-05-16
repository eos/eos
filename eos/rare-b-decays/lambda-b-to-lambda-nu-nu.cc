/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
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

#include <eos/rare-b-decays/lambda-b-to-lambda-nu-nu.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    template <>
    struct Implementation<LambdaBToLambdaDineutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

        Parameters parameters;

        UsedParameter m_Lambda_b;

        UsedParameter tau_Lambda_b;

        UsedParameter m_Lambda;

        UsedParameter g_fermi;

        UsedParameter alpha_e;

        UsedParameter hbar;

        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        std::function<complex<double> ()> lambda_t;
        std::function<WilsonCoefficients<wc::SBNuNu> ()> wc;

        GSL::QAGS::Config int_config;

        bool cp_conjugate;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            m_Lambda_b(p["mass::Lambda_b"], u),
            tau_Lambda_b(p["life_time::Lambda_b"], u),
            m_Lambda(p["mass::Lambda"], u),
            g_fermi(p["WET::G_Fermi"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            hbar(p["QM::hbar"], u),
            mu(p["sbnunu::mu"], u),
            int_config(GSL::QAGS::Config().epsrel(0.5e-3)),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false")))
        {
            form_factors = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda::" + o.get("form-factors", "BFvD2014"), p, o);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            lambda_t  = [*this] () { return model->ckm_tb() * std::conj(model->ckm_ts()); };
            wc        = [*this] () { return model->wet_sbnunu(cp_conjugate); };

            u.uses(*form_factors);
            u.uses(*model);
        }

        // differential decay width
        double differential_decay_width(const double & q2) const
        {
            const double lambda = eos::lambda(power_of<2>(m_Lambda_b), power_of<2>(m_Lambda), q2),
                         sqrt_lambda = std::sqrt(lambda);

            const auto wc = this->wc();

            const double rho_perp_v = 32 * q2 * sqrt_lambda *
                (power_of<2>(m_Lambda_b - m_Lambda) - q2) / power_of<4>(m_Lambda_b) * std::norm(form_factors->f_perp_v(q2));
            const double rho_perp_a = 32 * q2 * sqrt_lambda *
                (power_of<2>(m_Lambda_b + m_Lambda) - q2) / power_of<4>(m_Lambda_b) * std::norm(form_factors->f_perp_a(q2));
            const double rho_long_v = 16 * sqrt_lambda * power_of<2>(m_Lambda_b + m_Lambda) *
                (power_of<2>(m_Lambda_b - m_Lambda) - q2) / power_of<4>(m_Lambda_b) * std::norm(form_factors->f_long_v(q2));
            const double rho_long_a = 16 * sqrt_lambda * power_of<2>(m_Lambda_b - m_Lambda) *
                (power_of<2>(m_Lambda_b + m_Lambda) - q2) / power_of<4>(m_Lambda_b) * std::norm(form_factors->f_long_a(q2));

            const double norm = std::norm(g_fermi * alpha_e * lambda_t() * wc.cVL() *
                std::sqrt(m_Lambda_b / 3.0 / M_PI) / 16.0 / power_of<2>(M_PI));

            // assume the production of 3 diagonal neutrino flavors (nu_i nubar_i)
            return 3.0 * norm * (rho_perp_v + rho_perp_a + rho_long_v + rho_long_a);
        }

        // differential branching_ratio
        double differential_branching_ratio(const double & s) const
        {
            return differential_decay_width(s) * tau_Lambda_b / hbar;
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
    };

    double
    LambdaBToLambdaDineutrino::differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_branching_ratio(q2);
    }

    double
    LambdaBToLambdaDineutrino::integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<LambdaBToLambdaDineutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, q2_min, q2_max, _imp->int_config);
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
