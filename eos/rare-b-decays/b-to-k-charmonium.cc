/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Méril Reboud
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

#include <eos/form-factors/mesonic.hh>
#include <eos/rare-b-decays/b-to-k-charmonium.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <complex>

namespace eos
{
    using std::norm;
    using std::sqrt;

    /*!
     * Implementation for the decay @f$\bar{B} \to \bar{K} \psi@f$.
     */
    template <>
    struct Implementation<BToKCharmonium>
    {
        UsedParameter g_fermi;

        UsedParameter hbar;

        std::shared_ptr<Model> model;

        SwitchOption opt_q;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_K;

        SwitchOption opt_nonlocal_formfactor;

        NonlocalFormFactorPtr<nff::PToP> nonlocal_formfactor;

        SwitchOption opt_psi;

        UsedParameter m_psi;

        UsedParameter f_psi;

        std::function<complex<double> ()> residue_H_plus;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u),
            model(Model::make(o.get("model", "SM"), p, o)),
            opt_q(o, "q", { "d", "u" }, "d"),
            m_B(p["mass::B_" + opt_q.value()], u),
            tau_B(p["life_time::B_" + opt_q.value()], u),
            m_K(p["mass::K_" + opt_q.value()], u),
            opt_nonlocal_formfactor(o, "nonlocal-formfactor", { "GvDV2020", "GRvDV2021" }, "GvDV2020"),
            nonlocal_formfactor(NonlocalFormFactor<nff::PToP>::make("B->K::" + opt_nonlocal_formfactor.value(), p, o)),
            opt_psi(o, "psi", { "J/psi", "psi(2S)" }, "J/psi"),
            m_psi(p["mass::" + opt_psi.value()], u),
            f_psi(p["decay-constant::" + opt_psi.value()], u)
        {
            if (! nonlocal_formfactor.get())
                throw InternalError("Cannot construct the nonlocal formfactor");

            if ("J/psi" == opt_psi.value())
            {
                residue_H_plus = std::bind(&NonlocalFormFactor<nff::PToP>::H_plus_residue_jpsi, nonlocal_formfactor);
            }
            else
            {
                residue_H_plus = std::bind(&NonlocalFormFactor<nff::PToP>::H_plus_residue_psi2s, nonlocal_formfactor);
            }

            u.uses(*model);
            u.uses(*nonlocal_formfactor);
        }

        ~Implementation() = default;

        virtual double branching_ratio() const
        {
            const auto lambda = eos::lambda(pow(m_B, 2), pow(m_K, 2), pow(m_psi, 2));
            const auto prefactor = m_B * pow(g_fermi * abs(model->ckm_cb() * conj(model->ckm_cs())) / f_psi / pow(m_psi, 2), 2)
                    * tau_B() / hbar() * pow(lambda, 1.5) / ( 2.0 * M_PI );

            return prefactor * norm(this->residue_H_plus());
        }

    };

    BToKCharmonium::BToKCharmonium(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BToKCharmonium>(new Implementation<BToKCharmonium>(p, o, *this))
    {
    }

    BToKCharmonium::~BToKCharmonium() = default;

    double
    BToKCharmonium::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

}
