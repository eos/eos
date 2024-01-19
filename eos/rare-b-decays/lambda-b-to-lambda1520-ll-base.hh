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

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA1520_LL_BASE_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA1520_LL_BASE_HH 1

#include <eos/form-factors/baryonic.hh>
#include <eos/models/model.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda1520-ll.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    class LambdaBToLambda1520Dilepton::AmplitudeGenerator :
        public ParameterUser
    {
        public:
            std::shared_ptr<Model> model;
            std::shared_ptr<FormFactors<OneHalfPlusToThreeHalfMinus>> form_factors;
            LeptonFlavorOption opt_l;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;
            UsedParameter hbar;

            UsedParameter m_l;
            UsedParameter m_Lb;
            UsedParameter m_Lstar;

            bool cp_conjugate;
            LeptonFlavor lepton_flavor;

            static const std::vector<OptionSpecification> options;

            AmplitudeGenerator(const Parameters &, const Options &);

            double lambda(const double & q2) const;
            double beta_l(const double & q2) const;

            virtual ~AmplitudeGenerator();
            virtual LambdaBToLambda1520Dilepton::Amplitudes amplitudes(const double & q2) const = 0;
    };

    template <typename Tag_> class LambdaBToLambda1520DileptonAmplitudes;

    namespace tag
    {
        /*
         * Naive approache with neglected lepton masses and only LO factorizable corrections.
         */
        struct Naive;
    }
}

#endif
