/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017      Elena Graverini
 * Copyright (c) 2017-2026 Danny van Dyk
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

#include <eos/b-decays/lambdab-to-lambdac2625-l-nu.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/baryonic.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <math.h>

namespace eos
{
    template <>
    struct Implementation<LambdaBToLambdaC2625LeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<OneHalfPlusToThreeHalfMinus>> form_factors;

        Parameters parameters;

        UsedParameter m_LambdaB;

        UsedParameter tau_LambdaB;

        UsedParameter m_LambdaC2625;

        LeptonFlavorOption opt_l;

        UsedParameter m_l;

        UsedParameter g_fermi;

        UsedParameter hbar;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"_ov), p, o)),
            form_factors(FormFactorFactory<OneHalfPlusToThreeHalfMinus>::create("Lambda_b->Lambda_c(2625)::" + o.get("form-factors"_ok, "HQET"_ov).str(), p)),
            parameters(p),
            m_LambdaB(p["mass::Lambda_b"], u),
            tau_LambdaB(p["life_time::Lambda_b"], u),
            m_LambdaC2625(p["mass::Lambda_c(2625)"], u),
            opt_l(o, options, "l"_ok),
            m_l(p["mass::" + opt_l.str()], u),
            g_fermi(p["WET::G_Fermi"], u),
            hbar(p["QM::hbar"], u)
        {
            Context ctx("When constructing L_b->L_c(2625) lnu observable");

            u.uses(*form_factors);
            u.uses(*model);
        }

        double F12T(const double & q2) const { return form_factors->f_time12_v(q2); };
        double F120(const double & q2) const { return form_factors->f_long12_v(q2); };
        double F12P(const double & q2) const { return form_factors->f_perp12_v(q2); };
        double F32P(const double & q2) const { return form_factors->f_perp32_v(q2); };
        double G12T(const double & q2) const { return form_factors->f_time12_a(q2); };
        double G120(const double & q2) const { return form_factors->f_long12_a(q2); };
        double G12P(const double & q2) const { return form_factors->f_perp12_a(q2); };
        double G32P(const double & q2) const { return form_factors->f_perp32_a(q2); };
        double s_plus(const double & q2) const { return power_of<2>(m_LambdaB + m_LambdaC2625) - q2; };
        double s_minus(const double & q2) const { return power_of<2>(m_LambdaB - m_LambdaC2625) - q2; };

        // [BBGIOvD] parametrization for the differential decay width
        double a_l(const double & q2) const
        {
            double val = power_of<2>(F12T(q2)) * (power_of<2>(m_l) / q2) * power_of<2>(m_LambdaB - m_LambdaC2625);
            val += power_of<2>(F120(q2)) * power_of<2>(m_LambdaB + m_LambdaC2625) + (power_of<2>(F12P(q2)) + 3. * power_of<2>(F32P(q2))) * (power_of<2>(m_l) + q2);
            val += power_of<2>(G12T(q2)) * power_of<2>(m_l) / q2 * power_of<2>(m_LambdaB + m_LambdaC2625);
            val += power_of<2>(G120(q2)) * power_of<2>(m_LambdaB - m_LambdaC2625) + (power_of<2>(G12P(q2)) + 3. * power_of<2>(G32P(q2))) * (power_of<2>(m_l) + q2);
            return val;
        }

        double b_l(const double & q2) const
        {
            double val = 2. * (F12T(q2) * F120(q2) + G12T(q2) * G120(q2)) * power_of<2>(m_l) / q2;
            val *= power_of<2>(m_LambdaB) - power_of<2>(m_LambdaC2625);
            val += (-4.) * q2 * (F12P(q2) * G12P(q2) + 3. * F32P(q2) * G32P(q2));
            return val;
        }

        double c_l(const double & q2) const
        {
            double val = power_of<2>(F120(q2)) * power_of<2>(m_LambdaB + m_LambdaC2625);
            val += (-1.) * q2 * (power_of<2>(F12P(q2)) + 3. * power_of<2>(F32P(q2)));
            val += power_of<2>(G120(q2)) * power_of<2>(m_LambdaB - m_LambdaC2625);
            val += (-1.) * q2 * (power_of<2>(G12P(q2)) + 3. * power_of<2>(G32P(q2)));
            val *= (-1.) * (1. - power_of<2>(m_l) / q2);
            return val;
        }

        double gamma_0(const double & q2) const
        {
            double val = power_of<2>(g_fermi) * sqrt(s_plus(q2) * s_minus(q2));
            val *= m_LambdaB * m_LambdaC2625;
            val *= (1. / (96. * power_of<3>(M_PI * m_LambdaB)));
            val *= power_of<2>(1. - power_of<2>(m_l) / q2);
            return val;
        }

        inline double lambda(const double & q2) const
        {
            return s_plus(q2) * s_minus(q2);
        }

        // normalized to V_cb = 1
        double normalized_differential_decay_width(const double & q2) const
        {
            if ((q2 < m_l * m_l) || (lambda(q2) < 0.0))
            {
                return 0.0;
            }

            return 2. * gamma_0(q2) * (a_l(q2) + c_l(q2) / 3.);
        }

        double normalized_differential_forward_backward_asymmetry(const double & q2) const
        {
            if ((q2 < m_l * m_l) || (lambda(q2) < 0.0))
            {
                return 0.0;
            }

            return gamma_0(q2) * b_l(q2);
        }

        double normalized_double_differential_decay_width(const double & q2, const double & cos_theta_l) const
        {
            if ((q2 < m_l * m_l) || (lambda(q2) < 0.0))
            {
                return 0.0;
            }

            return gamma_0(q2) * (a_l(q2) + b_l(q2) * cos_theta_l + c_l(q2) * power_of<2>(cos_theta_l));
        }

        double differential_decay_width(const double & q2) const
        {
            return normalized_differential_decay_width(q2) * std::norm(model->ckm_cb());
        }

        double double_differential_decay_width(const double & q2, const double & cos_theta_l) const
        {
            return normalized_double_differential_decay_width(q2, cos_theta_l) * std::norm(model->ckm_cb());
        }

        double differential_branching_ratio(const double & q2) const
        {
            return differential_decay_width(q2) * tau_LambdaB / hbar;
        }

        double differential_forward_backward_asymmetry(const double & q2) const
        {
            if ((q2 < m_l * m_l) || (lambda(q2) < 0.0))
            {
                return 0.0;
            }

            return b_l(q2) / (2. * (a_l(q2) + c_l(q2) / 3.));
        }

        double double_differential_branching_ratio(const double & q2, const double & cos_theta_l) const
        {
            return double_differential_decay_width(q2, cos_theta_l) * tau_LambdaB / hbar;
        }

        double integrated_branching_ratio(const double & q2_min, const double & q2_max) const
        {
            std::function<double (const double &)> f = std::bind(&Implementation<LambdaBToLambdaC2625LeptonNeutrino>::differential_branching_ratio,
                    *this, std::placeholders::_1);

            return integrate<GSL::QAGS>(f, q2_min, q2_max);
        }

        double integrated_forward_backward_asymmetry(const double & q2_min, const double & q2_max) const
        {
            std::function<double (const double &)> numerator   = std::bind(&Implementation<LambdaBToLambdaC2625LeptonNeutrino>::normalized_differential_forward_backward_asymmetry,
                    *this, std::placeholders::_1);
            std::function<double (const double &)> denominator = std::bind(&Implementation<LambdaBToLambdaC2625LeptonNeutrino>::normalized_differential_decay_width,
                    *this, std::placeholders::_1);

            return integrate<GSL::QAGS>(numerator, q2_min, q2_max) / integrate<GSL::QAGS>(denominator, q2_min, q2_max);
        }
    };

    const std::vector<OptionSpecification>
    Implementation<LambdaBToLambdaC2625LeptonNeutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<OneHalfPlusToThreeHalfMinus>::option_specification(),
        { "l"_ok, { "e"_ov, "mu"_ov, "tau"_ov }, "mu"_ov }
    };

    LambdaBToLambdaC2625LeptonNeutrino::LambdaBToLambdaC2625LeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<LambdaBToLambdaC2625LeptonNeutrino>(new Implementation<LambdaBToLambdaC2625LeptonNeutrino>(parameters, options, *this))
    {
    }

    LambdaBToLambdaC2625LeptonNeutrino::~LambdaBToLambdaC2625LeptonNeutrino()
    {
    }

    double
    LambdaBToLambdaC2625LeptonNeutrino::a_l(const double & q2) const
    {
        return _imp->a_l(q2);
    }

    double
    LambdaBToLambdaC2625LeptonNeutrino::b_l(const double & q2) const
    {
        return _imp->b_l(q2);
    }

    double
    LambdaBToLambdaC2625LeptonNeutrino::c_l(const double & q2) const
    {
        return _imp->c_l(q2);
    }

    double
    LambdaBToLambdaC2625LeptonNeutrino::differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_branching_ratio(q2);
    }

    double
    LambdaBToLambdaC2625LeptonNeutrino::differential_forward_backward_asymmetry(const double & q2) const
    {
        return _imp->differential_forward_backward_asymmetry(q2);
    }

    double
    LambdaBToLambdaC2625LeptonNeutrino::double_differential_branching_ratio(const double & q2, const double & cos_theta_l) const
    {
        return _imp->double_differential_branching_ratio(q2, cos_theta_l);
    }

    double
    LambdaBToLambdaC2625LeptonNeutrino::integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_branching_ratio(q2_min, q2_max);
    }

    double
    LambdaBToLambdaC2625LeptonNeutrino::integrated_forward_backward_asymmetry(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_forward_backward_asymmetry(q2_min, q2_max);
    }

    double
    LambdaBToLambdaC2625LeptonNeutrino::normalized_integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        const double abs_q2_min = power_of<2>(_imp->m_l);
        const double abs_q2_max = power_of<2>(_imp->m_LambdaB - _imp->m_LambdaC2625);

        return _imp->integrated_branching_ratio(q2_min, q2_max) / _imp->integrated_branching_ratio(abs_q2_min, abs_q2_max);
    }

    const std::string
    LambdaBToLambdaC2625LeptonNeutrino::description = "\
The decay Lambda_b -> Lambda_c(2625) l nu, where l=e,mu,tau is a lepton.";

    const std::string
    LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_q2 = "\
The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_c_theta_l = "\
The cosine of the helicity angle between the direction of flight of the muon and of the Lambda_c(2625) in the l-nubar rest frame.";

    const std::set<ReferenceName>
    LambdaBToLambdaC2625LeptonNeutrino::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambdaC2625LeptonNeutrino::begin_options()
    {
        return Implementation<LambdaBToLambdaC2625LeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambdaC2625LeptonNeutrino::end_options()
    {
        return Implementation<LambdaBToLambdaC2625LeptonNeutrino>::options.cend();
    }
}
