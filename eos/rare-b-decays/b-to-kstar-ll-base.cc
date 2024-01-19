/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016, 2017 Danny van Dyk
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

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>

namespace eos
{
    BToKstarDilepton::AmplitudeGenerator::AmplitudeGenerator(const Parameters & p, const Options & o) :
        model(Model::make(o.get("model", "SM"), p, o)),
        form_factors(FormFactorFactory<PToV>::create("B->K^*::" + o.get("form-factors", "BSZ2015"), p)),
        opt_l(o, options, "l"),
        opt_cp_conjugate(o, options, "cp-conjugate"),
        mu(p["sb" + opt_l.str() + opt_l.str() + "::mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["WET::G_Fermi"], *this),
        hbar(p["QM::hbar"], *this),
        tau(p["life_time::B_" + o.get("q", "d")], *this),
        m_B(p["mass::B_" + o.get("q", "d")], *this),
        m_Kstar(p["mass::K_d^*"], *this),
        m_l(p["mass::" + opt_l.str()], *this),
        cp_conjugate(opt_cp_conjugate.value()),
        lepton_flavor(opt_l.value())
    {
        Context ctx("When constructing B->K^*ll amplitudes");

        if (0.0 == m_l())
        {
            throw InternalError("Zero lepton mass leads to NaNs in timelike amplitudes. Use tiny lepton mass > 0!");
        }

        this->uses(*form_factors);
        this->uses(*model);
    }

    BToKstarDilepton::AmplitudeGenerator::~AmplitudeGenerator()
    {
    }

    const std::vector<OptionSpecification>
    BToKstarDilepton::AmplitudeGenerator::options
    {
        Model::option_specification(),
        FormFactorFactory<PToV>::option_specification(),
        { "cp-conjugate", { "true", "false" },  "false" },
        { "l", { "e", "mu", "tau" }, "mu" },
    };

    double
    BToKstarDilepton::AmplitudeGenerator::beta_l(const double & s) const
    {
        return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
    }

    double
    BToKstarDilepton::AmplitudeGenerator::lambda(const double & s) const
    {
        return eos::lambda(m_B() * m_B(), m_Kstar() * m_Kstar(), s);
    }

    double
    BToKstarDilepton::AmplitudeGenerator::energy(const double & s) const
    {
        return (m_B() * m_B() + m_Kstar() * m_Kstar() - s) / (2.0 * m_B());
    }

    double
    BToKstarDilepton::AmplitudeGenerator::s_hat(const double & s) const
    {
        return s / m_B() / m_B();
    }

}
