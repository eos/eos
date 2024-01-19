/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Méril Reboud
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

#include <eos/rare-b-decays/lambda-b-to-lambda1520-gamma-base.hh>
#include <eos/utils/destringify.hh>

namespace eos
{
    LambdaBToLambda1520Gamma::AmplitudeGenerator::AmplitudeGenerator(const Parameters & p, const Options & o) :
        model(Model::make(o.get("model", "SM"), p, o)),
        form_factors(FormFactorFactory<OneHalfPlusToThreeHalfMinus>::create("Lambda_b->Lambda(1520)::" + o.get("form-factors", "ABR2022"), p)),
        hbar(p["QM::hbar"], *this),
        mu(p["sb::mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["WET::G_Fermi"], *this),
        m_Lb(p["mass::Lambda_b"], *this),
        m_Lstar(p["mass::Lambda(1520)"], *this),
        opt_cp_conjugate(o, options, "cp-conjugate"),
        cp_conjugate(opt_cp_conjugate.value())
    {
        Context ctx("When constructing Lb->L(1520)gamma amplitudes");

        this->uses(*form_factors);
        this->uses(*model);
    }

    LambdaBToLambda1520Gamma::AmplitudeGenerator::~AmplitudeGenerator() = default;

    const std::vector<OptionSpecification>
    LambdaBToLambda1520Gamma::AmplitudeGenerator::options
    {
        Model::option_specification(),
        FormFactorFactory<OneHalfPlusToThreeHalfMinus>::option_specification(),
        { "cp-conjugate", { "true", "false" },  "false" },
    };
}
