/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2017 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_ZERO_RECOIL_SUM_RULE_HH
#define EOS_GUARD_EOS_FORM_FACTORS_ZERO_RECOIL_SUM_RULE_HH 1

#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    template <typename Tag_> class ZeroRecoilSumRule;

    class LambdaBToC {};

    template <>
    class ZeroRecoilSumRule<LambdaBToC> :
        public ParameterUser,
        public PrivateImplementationPattern<ZeroRecoilSumRule<LambdaBToC>>
    {
        public:
            ZeroRecoilSumRule(const Parameters &, const Options &);

            ~ZeroRecoilSumRule();

            // inclusive bounds
            double vector_current() const;
            double axialvector_current() const;

            // inelastic contributions from orbitally excited doublet Lambda_c^*
            double vector_current_inel() const;
            double axialvector_current_inel() const;

            /* Diagnostics for unit tests */
            Diagnostics diagnostics() const;
    };
}

#endif
