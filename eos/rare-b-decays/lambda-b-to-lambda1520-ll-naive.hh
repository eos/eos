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

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA1520_LL_NAIVE_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA1520_LL_NAIVE_HH 1

#include <eos/rare-b-decays/lambda-b-to-lambda1520-ll-base.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>

namespace eos
{
    template <>
    class LambdaBToLambda1520DileptonAmplitudes<tag::Naive> :
        public LambdaBToLambda1520Dilepton::AmplitudeGenerator
    {
        public:
            UsedParameter m_c;

            LambdaBToLambda1520DileptonAmplitudes(const Parameters & p, const Options & o);
            ~LambdaBToLambda1520DileptonAmplitudes();

            virtual LambdaBToLambda1520Dilepton::Amplitudes amplitudes(const double & q2) const;

            double m_b_PS() const;
            double mu_f() const;
            double norm(const double & q2) const;
    };
}

#endif
