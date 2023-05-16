/*
 * Copyright (c) 2023, Méril Reboud
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA_NU_NU_IMPL_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA_NU_NU_IMPL_HH 1

#include <eos/observable.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-nu-nu.hh>

namespace eos
{
    struct LambdaBToLambdaDineutrino::AngularCoefficients
    {
        double K1ss, K1cc;

        AngularCoefficients()
        {
        }

        AngularCoefficients(const std::array<double, 2> & a) :
            K1ss(a[0]),
            K1cc(a[1])
        {
        }
    };

    class LambdaBToLambdaDineutrino::IntermediateResult :
        public CacheableObservable::IntermediateResult
    {
        public:
            LambdaBToLambdaDineutrino::AngularCoefficients ac;

            IntermediateResult()
            {
            }

            ~IntermediateResult() = default;
    };
}

#endif
