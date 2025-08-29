/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016-2025 Danny van Dyk
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

#include <eos/b-decays/b-to-pi-pi-l-nu.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>

namespace eos
{
    using namespace std::literals::string_literals;

    template <>
    struct Implementation<BToPiPiLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToPP>> form_factors;

        LeptonFlavorOption opt_l;
        QuarkFlavorOption opt_q;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_pi;

        UsedParameter m_l;

        UsedParameter g_fermi;

        UsedParameter hbar;

        // GSL elements for MC integration
        gsl_rng * rng;
        gsl_monte_function normalized_decay_width_integrand;
        gsl_monte_miser_state * state;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            form_factors(FormFactorFactory<PToPP>::create("B->pipi::" + o.get("form-factors"_ok, "BFvD2016"), p, o)),
            opt_l(o, options, "l"_ok),
            opt_q(o, options, "q"_ok),
            m_B(p["mass::B_" + opt_q.str()], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            m_pi(p["mass::pi^" + std::string(opt_q.value() == QuarkFlavor::down ? "+" : "0")], u),
            m_l(p["mass::" + opt_l.str()], u),
            g_fermi(p["WET::G_Fermi"], u),
            hbar(p["QM::hbar"], u),
            rng(gsl_rng_alloc(gsl_rng_mt19937)),
            state(gsl_monte_miser_alloc(3u))
        {
            Context ctx("When constructing B->pipilnu observable");

            u.uses(*form_factors);
            u.uses(*model);
        }

        ~Implementation()
        {
            gsl_monte_miser_free(state);
            gsl_rng_free(rng);
        }

        // normalized to V_ub = 1
        double normalized_differential_decay_width(const double & q2, const double & k2, const double & z) const
        {
            const double m_B = this->m_B(),
                  m_B2 = m_B * m_B,
                  m_B3 = m_B * m_B2;
            const double lambda = eos::lambda(m_B2, q2, k2),
                  sqrt_lambda = sqrt(lambda);
            const double m_l2 = m_l() * m_l();
            const double beta_l = 1.0 - m_l2 / q2;
            const double norm = power_of<2>(g_fermi()) * beta_l * q2 * sqrt_lambda
                / (3072.0 * power_of<5>(M_PI) * m_B3);

            const complex<double> F_perp = form_factors->f_perp(q2, k2, z);
            const complex<double> F_para = form_factors->f_para(q2, k2, z);
            const complex<double> F_long = form_factors->f_long(q2, k2, z);
            const complex<double> F_time = form_factors->f_time(q2, k2, z);

            return norm * beta_l / 4.0 *
                (
                    (3.0 - beta_l) * std::norm(F_long)
                    + (1.0 - z * z) * (3.0 - beta_l) * (std::norm(F_perp) + std::norm(F_para))
                    + 3.0 * m_l2 * std::norm(F_time)
                );
        }

        double differential_branching_ratio(const double & q2, const double & k2, const double & z) const
        {
            return differential_decay_width(q2, k2, z) * tau_B / hbar;
        }

        double differential_decay_width(const double & q2, const double & k2, const double & z) const
        {
            return normalized_differential_decay_width(q2, k2, z) * std::norm(model->ckm_ub());
        }

        static double normalized_differential_decay_width_gsl_adapter(double * parameters, size_t dim, void * _imp)
        {
            if (dim != 3u)
                throw InternalError("Implemenation<BToPiPiLeptonNeutrino::normalized_differentia_decay_width_gsl_adapter(): wrong number of parameters!");

            const double q2 = parameters[0];
            const double k2 = parameters[1];
            const double z  = parameters[2];

            auto imp = reinterpret_cast<const Implementation<BToPiPiLeptonNeutrino> *>(_imp);

            if ((lambda(q2, k2, imp->m_B() * imp->m_B()) <= 0) || (q2 <= power_of<2>(imp->m_l())))
                return 0.0;

            return imp->normalized_differential_decay_width(q2, k2, z);
        }

        double normalized_integrated_decay_width(const double & q2min, const double & q2max,
                const double k2min, const double & k2max,
                const double & zmin, const double & zmax) const
        {
            // Yields a numerical error of approximately 0.2%.
            static const size_t calls = 50000;

            const double x_min[3] = { q2min, k2min, zmin };
            const double x_max[3] = { q2max, k2max, zmax };

            gsl_monte_function integrand{ &normalized_differential_decay_width_gsl_adapter, 3u, const_cast<void *>(reinterpret_cast<const void *>(this)) };

            double result, error;

            gsl_monte_miser_integrate(&integrand, x_min, x_max, 3u, calls, rng, state, &result, &error);

            return result;
        }

        double normalized_integrated_forward_backward_asymmetry(const double & q2min, const double & q2max,
                const double k2min, const double & k2max) const
        {
            // Yields a numerical error of approximately 0.2%.
            static const size_t calls = 50000;

            const double x_forward_min[3]  = { q2min, k2min,  0.0 };
            const double x_forward_max[3]  = { q2max, k2max, +1.0 };
            const double x_backward_min[3] = { q2min, k2min, -1.0 };
            const double x_backward_max[3] = { q2max, k2max,  0.0 };

            gsl_monte_function integrand{ &normalized_differential_decay_width_gsl_adapter, 3u, const_cast<void *>(reinterpret_cast<const void *>(this)) };

            double forward,  forward_error;
            double backward, backward_error;

            gsl_monte_miser_integrate(&integrand, x_forward_min, x_forward_max, 3u, calls, rng, state, &forward, &forward_error);
            gsl_monte_miser_integrate(&integrand, x_backward_min, x_backward_max, 3u, calls, rng, state, &backward, &backward_error);

            return (forward - backward) / (forward + backward);
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToPiPiLeptonNeutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToPP>::option_specification(),
        { "l"_ok, { "e"s, "mu"s }, "mu"s },
        { "q"_ok, { "d"s, "u"s },  "d"s  }
    };

    BToPiPiLeptonNeutrino::BToPiPiLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToPiPiLeptonNeutrino>(new Implementation<BToPiPiLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToPiPiLeptonNeutrino::~BToPiPiLeptonNeutrino()
    {
    }

    double
    BToPiPiLeptonNeutrino::double_differential_branching_ratio(const double & q2, const double & k2) const
    {
        std::function<double (const double &)> integrand = std::bind(
                &Implementation<BToPiPiLeptonNeutrino>::differential_branching_ratio,
                _imp.get(), q2, k2, std::placeholders::_1);

        return integrate<GSL::QNG>(integrand, -1.0, +1.0);
    }

    double
    BToPiPiLeptonNeutrino::triple_differential_branching_ratio(const double & q2, const double & k2, const double & z) const
    {
        return _imp->differential_branching_ratio(q2, k2, z);
    }

    double
    BToPiPiLeptonNeutrino::double_differential_forward_backward_asymmetry(const double & q2, const double & k2) const
    {
        std::function<double (const double &)> integrand = std::bind(
                &Implementation<BToPiPiLeptonNeutrino>::normalized_differential_decay_width,
                _imp.get(), q2, k2, std::placeholders::_1);

        const double numerator = integrate<GSL::QNG>(integrand, 0.0, +1.0) - integrate<GSL::QNG>(integrand, -1.0, 0.0);
        const double denominator = integrate<GSL::QNG>(integrand, -1.0, +1.0);

        return numerator / denominator;
    }


    double
    BToPiPiLeptonNeutrino::partial_waves(const double & q2, const double & k2, const double & z) const
    {
        std::function<double (const double &)> integrand = std::bind(
                &Implementation<BToPiPiLeptonNeutrino>::normalized_differential_decay_width,
                _imp.get(), q2, k2, std::placeholders::_1);

        return integrand(z)
            / integrate<GSL::QNG>(integrand, -1.0, +1.0);
    }

    double
    BToPiPiLeptonNeutrino::integrated_branching_ratio(const double & q2min, const double & q2max,
                const double & k2min, const double & k2max,
                const double & zmin, const double & zmax) const
    {
        return _imp->normalized_integrated_decay_width(q2min, q2max, k2min, k2max, zmin, zmax)
            * std::norm(_imp->model->ckm_ub()) * _imp->tau_B / _imp->hbar;
    }

    double
    BToPiPiLeptonNeutrino::integrated_forward_backward_asymmetry(const double & q2min, const double & q2max,
                const double & k2min, const double & k2max) const
    {
        return _imp->normalized_integrated_forward_backward_asymmetry(q2min, q2max, k2min, k2max);
    }

    const std::string
    BToPiPiLeptonNeutrino::description = "\
The decay B->pi pi l nubar, where l=e,mu is a light lepton, see e.g. [FFKMvD:2013].";

    const std::string
    BToPiPiLeptonNeutrino::kinematics_description_q2 = "\
The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    BToPiPiLeptonNeutrino::kinematics_description_k2 = "\
The invariant mass of the pi-pi pair in GeV^2.";

    const std::string
    BToPiPiLeptonNeutrino::kinematics_description_z = "\
The cosine of the pion helicity angle in the pi-pi rest frame.";

    const std::set<ReferenceName>
    BToPiPiLeptonNeutrino::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BToPiPiLeptonNeutrino::begin_options()
    {
        return Implementation<BToPiPiLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToPiPiLeptonNeutrino::end_options()
    {
        return Implementation<BToPiPiLeptonNeutrino>::options.cend();
    }
}
