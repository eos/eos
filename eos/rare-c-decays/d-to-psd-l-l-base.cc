/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026    Carolina Bolognani
 * Copyright (c) 2026    Dominik Suelmann
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
#include <eos/rare-c-decays/d-to-psd-l-l-base.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>

namespace eos
{

    DToPseudoscalarLeptonLepton::AmplitudeGenerator::AmplitudeGenerator(const Parameters & p, const Options & o) :
        model(Model::make(o.get("model"_ok, "SM"_ov), p, o)),
        form_factors(FormFactorFactory<PToP>::create(_process() + "::" + o.get("form-factors"_ok, "SSE"_ov).str(), p)),
        opt_l(o, options, "l"_ok),
        opt_q(o, options, "q"_ok),
        opt_P(o, options, "P"_ok),
        mu(p["uc::mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["WET::G_Fermi"], *this),
        m_D(p["mass::" + _D()], *this),
        m_P(p["mass::" + _P()], *this),
        m_l(p["mass::" + opt_l.str()], *this),
        opt_cp_conjugate(o, options, "cp-conjugate"_ok),
        cp_conjugate(opt_cp_conjugate.value()),
        isospin_factor(_isospin_factor()),
        q(opt_q.value()),
        lepton_flavor(opt_l.value())
    {
        Context ctx("When constructing D(s)->Pll amplitudes");

        switch (q)
        {
            case QuarkFlavor::strange:
                m_Q_msbar = [this](const double & mu) -> double { return model->m_s_msbar(mu); };
                break;
            case QuarkFlavor::down:
                m_Q_msbar = [this](const double & mu) -> double { return model->m_d_msbar(mu); };
                break;
            case QuarkFlavor::up:
                m_Q_msbar = [this](const double & mu) -> double { return model->m_u_msbar(mu); };
                break;
            default: throw InternalError("Invalid quark flavor: " + stringify(q));
        }

        if (0.0 == m_l())
        {
            throw InternalError("Zero lepton mass leads to NaNs in timelike amplitudes. Use tiny lepton mass > 0!");
        }

        this->uses(*form_factors);
        this->uses(*model);
    }

    DToPseudoscalarLeptonLepton::AmplitudeGenerator::~AmplitudeGenerator() {}

    const std::map<std::tuple<QuarkFlavor, std::string>, std::tuple<std::string, std::string, std::string, double>>
            DToPseudoscalarLeptonLepton::AmplitudeGenerator::process_map{
                {             { QuarkFlavor::up, "pi" }, { "D->pi", "D_u", "pi^0", 1.0 / std::sqrt(2.0) } },
                {           { QuarkFlavor::down, "pi" }, { "D->pi", "D_d", "pi^+" , 1.0 } },
                {         { QuarkFlavor::strange, "K" }, { "D_s->K", "D_s", "K_u" , 1.0 } },
    };

    const std::vector<OptionSpecification> DToPseudoscalarLeptonLepton::AmplitudeGenerator::options{
        Model::option_specification(),
        FormFactorFactory<PToP>::option_specification(),
        { "cp-conjugate"_ok,     { "true"_ov, "false"_ov }, "false"_ov },
        {            "l"_ok, { "e"_ov, "mu"_ov, "tau"_ov },    "mu"_ov },
    };

    double
    DToPseudoscalarLeptonLepton::AmplitudeGenerator::mqatmu() const
    {
        return model->m_s_msbar(mu);
    }

    double
    DToPseudoscalarLeptonLepton::AmplitudeGenerator::beta_l(const double & q2) const
    {
        return std::sqrt(1.0 - 4.0 * m_l() * m_l() / q2);
    }

    double
    DToPseudoscalarLeptonLepton::AmplitudeGenerator::lambda(const double & q2) const
    {
        return eos::lambda(m_D() * m_D(), m_P() * m_P(), q2);
    }

    double
    DToPseudoscalarLeptonLepton::AmplitudeGenerator::energy(const double & q2) const
    {
        return (m_D() * m_D() + m_P() * m_P() - q2) / (2.0 * m_D());
    }

    double
    DToPseudoscalarLeptonLepton::AmplitudeGenerator::xi_pseudo(const double & q2) const
    {
        // cf. [BF:2001A], Eq. (22)
        return form_factors->f_p(q2);
    }

    double
    DToPseudoscalarLeptonLepton::AmplitudeGenerator::normalisation(const double & q2) const
    {
        // cf. [BHP:2007A], Eq. (4.2) - (4.4)
        double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));

        return power_of<2>(this->isospin_factor) * power_of<2>(g_fermi * alpha_e() * lambda_t) * sqrt(lambda(q2)) * beta_l(q2) * xi_pseudo(q2) * xi_pseudo(q2) / (512.0 * power_of<5>(M_PI) * power_of<3>(m_D()));
    }
} // namespace eos
