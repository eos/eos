/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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
#include <eos/utils/destringify.hh>

namespace eos
{

    BToKstarGamma::AmplitudeGenerator::AmplitudeGenerator(const Parameters & p, const Options & o) :
        model(Model::make(o.get("model"_ok, "SM"_ov), p, o)),
        form_factors(FormFactorFactory<PToV>::create("B->K^*::" + o.get("form-factors"_ok, "SSE"_ov).str(), p)),
        hbar(p["QM::hbar"], *this),
        mu(p["sb::mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["WET::G_Fermi"], *this),
        q(o, options, "q"_ok),
        tau(p["life_time::B_" + q.str()], *this),
        m_B(p["mass::B_" + q.str()], *this),
        m_Kstar(p["mass::K_d^*"], *this),
        l(o, options, "l"_ok),
        m_l(p["mass::" + l.str()], *this),
        opt_cp_conjugate(o, options, "cp-conjugate"_ok),
        cp_conjugate(opt_cp_conjugate.value())
    {
        Context ctx("When constructing B->K^*gamma amplitudes");

        switch (q.value())
        {
            case QuarkFlavor::down: e_q = -1.0 / 3.0; break;

            case QuarkFlavor::up: e_q = 2.0 / 3.0; break;

            default: throw InternalError("Unexpected quark flavor: '" + q.str() + "'");
        }

        this->uses(*form_factors);
        this->uses(*model);
    }

    BToKstarGamma::AmplitudeGenerator::~AmplitudeGenerator() = default;

    const std::vector<OptionSpecification> BToKstarGamma::AmplitudeGenerator::options{
        Model::option_specification(),
        FormFactorFactory<PToV>::option_specification(),
        {            "l"_ok,       { "e"_ov, "mu"_ov },    "mu"_ov },
        {            "q"_ok,        { "d"_ov, "u"_ov },     "d"_ov },
        { "cp-conjugate"_ok, { "true"_ov, "false"_ov }, "false"_ov }
    };
} // namespace eos
