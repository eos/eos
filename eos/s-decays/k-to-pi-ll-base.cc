/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Danny van Dyk
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

#include <eos/maths/power-of.hh>
#include <eos/s-decays/k-to-pi-ll-base.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>

namespace eos
{
    KToPiDilepton::AmplitudeGenerator::AmplitudeGenerator(const Parameters & p, const Options & o) :
        model(Model::make(o.get("model", "SM"), p, o)),
        form_factors(FormFactorFactory<PToP>::create("B->K::" + o.get("form-factors", "KMPW2010"), p)),
        opt_l(o, options, "l"),
        mu(p["ds" + opt_l.str() + opt_l.str() + "::mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["WET::G_Fermi"], *this),
        hbar(p["QM::hbar"], *this),
        opt_q(o, options, "q"),
        tau(p[opt_q.value() == QuarkFlavor::up ? "life_time::K_u" : "life_time::K_d" ], *this),
        m_K(p[opt_q.value() == QuarkFlavor::up ? "mass::K_u" : "mass::K_d"], *this),
        m_pi(p[opt_q.value() == QuarkFlavor::up ? "mass::pi^-" : "mass::pi^0"], *this),
        m_l(p["mass::" + opt_l.str()], *this),
        opt_cp_conjugate(o, options, "cp-conjugate"),
        cp_conjugate(opt_cp_conjugate.value()),
        lepton_flavor(opt_l.value())
    {
        Context ctx("When constructing K->pill amplitudes");

        if (0.0 == m_l())
        {
            throw InternalError("Zero lepton mass leads to NaNs in timelike amplitudes. Use tiny lepton mass > 0!");
        }

        this->uses(*form_factors);
        this->uses(*model);
    }

    KToPiDilepton::AmplitudeGenerator::~AmplitudeGenerator()
    {
    }

    const std::vector<OptionSpecification>
    KToPiDilepton::AmplitudeGenerator::options
    {
        Model::option_specification(),
        FormFactorFactory<PToP>::option_specification(),
        { "cp-conjugate", { "true", "false" },  "false" },
        { "q", { "u", "d" },                    "u"     },
        { "l", { "e", "mu" },                   "mu"    },
    };

    double
    KToPiDilepton::AmplitudeGenerator::beta_l(const double & s) const
    {
        return std::sqrt(1.0 - 4.0 * m_l() * m_l() / s);
    }

    double
    KToPiDilepton::AmplitudeGenerator::lambda(const double & q2) const
    {
        return eos::lambda(m_K() * m_K(), m_pi() * m_pi(), q2);
    }

    double
    KToPiDilepton::AmplitudeGenerator::energy(const double & q2) const
    {
        return (m_K() * m_K() + m_pi() * m_pi() - q2) / (2.0 * m_K());
    }

    double
    KToPiDilepton::AmplitudeGenerator::normalisation(const double & s) const
    {
        // cf. [BHP2007], Eq. (4.2) - (4.4)
        double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));

        return power_of<2>(g_fermi * alpha_e() * lambda_t) * sqrt(lambda(s)) * beta_l(s) /
                       (512.0 * power_of<5>(M_PI) * power_of<3>(m_K()));
    }
}
