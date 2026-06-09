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

#include <eos/form-factors/parametric-se-impl-onehalfplus-to-onehalfplus.hh>
#include <eos/form-factors/parametric-se-impl-onehalfplus-to-threehalfminus.hh>
#include <eos/form-factors/parametric-se-impl-p-to-p.hh>
#include <eos/form-factors/parametric-se-impl-p-to-v.hh>

namespace eos
{
    // 1/2+ -> 1/2+
    template class SEFormFactors<LambdaBToLambda,  OneHalfPlusToOneHalfPlus>;
    template class SEFormFactors<LambdaCToLambda,  OneHalfPlusToOneHalfPlus>;
    template class SEFormFactors<LambdaCToNeutron, OneHalfPlusToOneHalfPlus>;

    // 1/2+ -> 3/2-
    template class SEFormFactors<LambdaBToLambda1520, OneHalfPlusToThreeHalfMinus>;

    // b -> s
    template class SEFormFactors<BToK, PToP>;
    template class SEFormFactors<BToKstar, PToV>;
    template class SEFormFactors<BsToPhi, PToV>;
    template class SEFormFactors<BsToEta, PToP>;
    template class SEFormFactors<BsToEtaPrime, PToP>;

    // b -> u
    template class SEFormFactors<BToEta, PToP>;
    template class SEFormFactors<BToEtaPrime, PToP>;
    template class SEFormFactors<BsToK, PToP>;

    // c -> s
    template class SEFormFactors<DToK, PToP>;
    template class SEFormFactors<DsToEta, PToP>;
    template class SEFormFactors<DsToEtaPrime, PToP>;
}
