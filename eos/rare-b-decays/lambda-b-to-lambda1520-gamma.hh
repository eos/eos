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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA1520_GAMMA_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA1520_GAMMA_HH 1

#include <eos/maths/complex.hh>
#include <eos/rare-b-decays/decays.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /*!
     * Calculates observables in Lambda_b -> Lambda1520 gamma decays
     */
    class LambdaBToLambda1520Gamma :
        public ParameterUser,
        public PrivateImplementationPattern<LambdaBToLambda1520Gamma>
    {
        public:
            LambdaBToLambda1520Gamma(const Parameters & parameters, const Options & options);
            ~LambdaBToLambda1520Gamma();

            struct Amplitudes;
            class AmplitudeGenerator;

            /*!
             * @name Simple observables
             */
            /// @{
            /// Decay Rate
            double decay_rate() const;

            /// Branching Ratio
            double branching_ratio() const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };

    /*!
     * Amplitudes for the decay Lambda_ -> Lambda1520 gamma.
     */
    struct LambdaBToLambda1520Gamma::Amplitudes
    {
        complex<double> a_perp12;
        complex<double> a_para12;
        complex<double> a_perp32;
        complex<double> a_para32;
    };
}

#endif
