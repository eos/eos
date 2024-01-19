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

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA1520_GAMMA_BASE_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA1520_GAMMA_BASE_HH 1

#include <eos/form-factors/baryonic.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda1520-gamma.hh>
#include <eos/models/model.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    class LambdaBToLambda1520Gamma::AmplitudeGenerator :
        public ParameterUser
    {
        public:
            std::shared_ptr<Model> model;
            std::shared_ptr<FormFactors<OneHalfPlusToThreeHalfMinus>> form_factors;

            UsedParameter hbar;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;

            UsedParameter m_Lb;
            UsedParameter m_Lstar;

            BooleanOption opt_cp_conjugate;
            bool cp_conjugate;

            static const std::vector<OptionSpecification> options;

            AmplitudeGenerator(const Parameters &, const Options &);

            virtual ~AmplitudeGenerator();
            virtual LambdaBToLambda1520Gamma::Amplitudes amplitudes() const = 0;
    };

    template <typename Tag_> class LambdaBToLambda1520GammaAmplitudes;

    namespace tag
    {
        struct Naive;
    }
}

#endif
