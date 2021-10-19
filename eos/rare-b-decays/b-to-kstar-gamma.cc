/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2015, 2016 Danny van Dyk
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

#include <eos/rare-b-decays/b-to-kstar-gamma-base.hh>
#include <eos/rare-b-decays/b-to-kstar-gamma-bfs2004.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/save.hh>

#include <cmath>

#include <gsl/gsl_sf.h>

namespace eos
{
    using namespace std::placeholders;
    using std::norm;

    /*!
     * Implementation for the decay @f$\bar{B} \to \bar{K}^* \gamma@f$.
     */
    template <>
    struct Implementation<BToKstarGamma>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;

        SwitchOption q;

        UsedParameter tau;

        SwitchOption tag;

        std::shared_ptr<BToKstarGamma::AmplitudeGenerator> amplitude_generator;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["QM::hbar"], u),
            q(o, "q", { "d", "u" }, "d"),
            tau(p["life_time::B_" + q.value()], u),
            tag(o, "tag", { "BFS2004"})
        {
            if ("BFS2004" == tag.value())
            {
                amplitude_generator.reset(new BToKstarGammaAmplitudes<tag::BFS2004>(p, o));
            }
            else
            {
                throw InternalError("BToKstarGamma: Unknown tag or no valid tag specified (tag = '" + tag.value() + "')!");
            }

            u.uses(*model);
            u.uses(*amplitude_generator);
        }

        double decay_rate()
        {
            auto amps = amplitude_generator->amplitudes();

            return std::norm(amps.a_perp) + std::norm(amps.a_para);
        }

        double s_kstar_gamma()
        {
            Save<bool> save(amplitude_generator->cp_conjugate, false);

            // S_K^*gamma is calculated for B as the first state, Bbar as the second.
            // opposite order than in B->K^*ll.
            BToKstarGamma::Amplitudes abar = amplitude_generator->amplitudes();
            amplitude_generator->cp_conjugate = true;
            BToKstarGamma::Amplitudes a = amplitude_generator->amplitudes();

            double phi_d = arg(pow(conj(model->ckm_td()) * model->ckm_tb(), 2));
            complex<double> q_over_p = std::polar(1.0, -phi_d);

            auto a_left     = (a.a_para    + a.a_perp)    / sqrt(2.0);
            auto a_right    = (a.a_para    - a.a_perp)    / sqrt(2.0);
            auto abar_left  = (abar.a_para + abar.a_perp) / sqrt(2.0);
            auto abar_right = (abar.a_para - abar.a_perp) / sqrt(2.0);

            double numerator = -2.0 * imag(q_over_p * (conj(a_left) * abar_right + conj(a_right) * abar_left));
            double denominator = std::norm(a_left) + std::norm(a_right) + std::norm(abar_left) + std::norm(abar_right);

            return numerator / denominator;
        }

        double c_kstar_gamma()
        {
            Save<bool> save(amplitude_generator->cp_conjugate, false);

            // S_K^*gamma is calculated for B as the first state, Bbar as the second.
            // opposite order than in B->K^*ll.
            BToKstarGamma::Amplitudes abar = amplitude_generator->amplitudes();
            amplitude_generator->cp_conjugate = true;
            BToKstarGamma::Amplitudes a = amplitude_generator->amplitudes();

            auto a_left     = (a.a_para    + a.a_perp)    / sqrt(2.0);
            auto a_right    = (a.a_para    - a.a_perp)    / sqrt(2.0);
            auto abar_left  = (abar.a_para + abar.a_perp) / sqrt(2.0);
            auto abar_right = (abar.a_para - abar.a_perp) / sqrt(2.0);

            double numerator = std::norm(a_left) + std::norm(a_right) - std::norm(abar_left) - std::norm(abar_right);
            double denominator = std::norm(a_left) + std::norm(a_right) + std::norm(abar_left) + std::norm(abar_right);

            return numerator / denominator;
        }

        double isospin_asymmetry()
        {
            Save<char> save_q(amplitude_generator->q, 'd');
            Save<double> save_eq(amplitude_generator->e_q, -1.0/3.0);

            double gamma_neutral = decay_rate();
            amplitude_generator->q = 'u';
            amplitude_generator->e_q = +2.0/3.0;
            double gamma_charged = decay_rate();

            return (gamma_neutral - gamma_charged) / (gamma_neutral + gamma_charged);
        }
    };

    BToKstarGamma::BToKstarGamma(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToKstarGamma>(new Implementation<BToKstarGamma>(parameters, options, *this))
    {
    }

    BToKstarGamma::~BToKstarGamma() = default;

    double
    BToKstarGamma::branching_ratio() const
    {
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        double Gamma_B = _imp->hbar() / _imp->tau;
        double gamma   = _imp->decay_rate();

        return gamma / Gamma_B;
    }

    double
    BToKstarGamma::branching_ratio_cp_averaged() const
    {
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        double Gamma_B = _imp->hbar() / _imp->tau;

        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);
        double gamma    = _imp->decay_rate();
        _imp->amplitude_generator->cp_conjugate = true;
        double gammabar = _imp->decay_rate();

        return (gamma + gammabar) / (2.0 * Gamma_B);
    }

    double
    BToKstarGamma::cp_asymmetry() const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);
        double gamma    = _imp->decay_rate();
        _imp->amplitude_generator->cp_conjugate = true;
        double gammabar = _imp->decay_rate();

        return (gamma - gammabar) / (gamma + gammabar);
    }

    double
    BToKstarGamma::s_kstar_gamma() const
    {
        return _imp->s_kstar_gamma();
    }

    double
    BToKstarGamma::c_kstar_gamma() const
    {
        return _imp->c_kstar_gamma();
    }

    double
    BToKstarGamma::isospin_asymmetry() const
    {
        return _imp->isospin_asymmetry();
    }

    const std::set<ReferenceName>
    BToKstarGamma::references
    {
    };
}
