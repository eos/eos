/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
 * Copyright (c) 2021      Méril Reboud
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
#include <eos/rare-b-decays/b-to-k-ll-base.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>

namespace eos
{

    BToKDilepton::AmplitudeGenerator::AmplitudeGenerator(const Parameters & p, const Options & o) :
        model(Model::make(o.get("model"_ok, "SM"_ov), p, o)),
        form_factors(FormFactorFactory<PToP>::create("B->K::" + o.get("form-factors"_ok, "KMPW2010"_ov).str(), p)),
        opt_l(o, options, "l"_ok),
        mu(p["sb" + opt_l.str() + opt_l.str() + "::mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["WET::G_Fermi"], *this),
        hbar(p["QM::hbar"], *this),
        tau(p["life_time::B_" + o.get("q"_ok, "d"_ov).str()], *this),
        m_B(p["mass::B_" + o.get("q"_ok, "d"_ov).str()], *this),
        m_K(p["mass::K_" + o.get("q"_ok, "d"_ov).str()], *this),
        m_l(p["mass::" + opt_l.str()], *this),
        opt_cp_conjugate(o, options, "cp-conjugate"_ok),
        cp_conjugate(opt_cp_conjugate.value()),
        lepton_flavor(opt_l.value())
    {
        Context ctx("When constructing B->Kll amplitudes");

        if (0.0 == m_l())
        {
            throw InternalError("Zero lepton mass leads to NaNs in timelike amplitudes. Use tiny lepton mass > 0!");
        }

        this->uses(*form_factors);
        this->uses(*model);
    }

    BToKDilepton::AmplitudeGenerator::~AmplitudeGenerator()
    {
    }

    const std::vector<OptionSpecification>
    BToKDilepton::AmplitudeGenerator::options
    {
        Model::option_specification(),
        FormFactorFactory<PToP>::option_specification(),
        { "cp-conjugate"_ok, { "true"_ov, "false"_ov },  "false"_ov },
        { "l"_ok, { "e"_ov, "mu"_ov, "tau"_ov }, "mu"_ov },
    };

    double
    BToKDilepton::AmplitudeGenerator::beta_l(const double & q2) const
    {
        return std::sqrt(1.0 - 4.0 * m_l() * m_l() / q2);
    }

    double
    BToKDilepton::AmplitudeGenerator::lambda(const double & q2) const
    {
        return eos::lambda(m_B() * m_B(), m_K() * m_K(), q2);
    }

    double
    BToKDilepton::AmplitudeGenerator::energy(const double & q2) const
    {
        return (m_B() * m_B() + m_K() * m_K() - q2) / (2.0 * m_B());
    }

    double
    BToKDilepton::AmplitudeGenerator::xi_pseudo(const double & q2) const
    {
        // cf. [BF:2001A], Eq. (22)
        return form_factors->f_p(q2);
    }

    double
    BToKDilepton::AmplitudeGenerator::normalisation(const double & q2) const
    {
        // cf. [BHP:2007A], Eq. (4.2) - (4.4)
        double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));

        return power_of<2>(g_fermi * alpha_e() * lambda_t) * sqrt(lambda(q2)) * beta_l(q2) * xi_pseudo(q2) * xi_pseudo(q2) /
                       (512.0 * power_of<5>(M_PI) * power_of<3>(m_B()));
    }
}
