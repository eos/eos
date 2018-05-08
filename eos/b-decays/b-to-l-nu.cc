/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2018 Christoph Bobeth
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

#include <eos/b-decays/b-to-l-nu.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    using std::norm;

    /*
     * Decay: B_q -> l nubar, cf. [FMvD2015]
     */
    template <>
    struct Implementation<BToLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter m_B;

        UsedParameter f_B;

        UsedParameter tau_B;

        SwitchOption opt_l;
        
        UsedParameter m_l;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["hbar"], u),
            g_fermi(p["G_Fermi"], u),
            m_B(p["mass::B_u"], u),
            f_B(p["decay-constant::B_u"], u),
            tau_B(p["life_time::B_u"], u),
            opt_l(o, "l", {"e", "mu", "tau"}, "tau"),
            m_l(p["mass::" + opt_l.value()], u)
        {
            u.uses(*model);
        }

        inline double beta_l() const
        {
            return (1.0 - power_of<2>(m_l() / m_B()));
        }

        double decay_width() const
        {
            const WilsonCoefficients<BToU> wc = model->wilson_coefficients_b_to_u(opt_l.value(), false);

            // cf. [DBG2013], eq. (5), p. 5
            const complex<double> ga = wc.cvl() - wc.cvr();

            return power_of<2>(g_fermi * std::abs(model->ckm_ub()) * f_B * m_l * beta_l())
                * m_B / (8.0 * M_PI)
                * norm(ga);
        }

        double branching_ratio() const
        {
            return decay_width() * tau_B / hbar;
        }
    };

    BToLeptonNeutrino::BToLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToLeptonNeutrino>(new Implementation<BToLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToLeptonNeutrino::~BToLeptonNeutrino()
    {
    }

    double
    BToLeptonNeutrino::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    BToLeptonNeutrino::decay_width() const
    {
        return _imp->decay_width();
    }
}
