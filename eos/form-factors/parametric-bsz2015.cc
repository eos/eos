/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015 Frederik Beaujean
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

#include <eos/form-factors/parametric-bsz2015-impl.hh>

namespace eos
{
    template class BSZ2015FormFactors<BToDstar, PToV>;
    template class BSZ2015FormFactors<BToKstar, PToV>;
    template class BSZ2015FormFactors<BToOmega, PToV>;
    template class BSZ2015FormFactors<BToRho, PToV>;
    template class BSZ2015FormFactors<BsToDsstar, PToV>;
    template class BSZ2015FormFactors<BsToKstar, PToV>;
    template class BSZ2015FormFactors<BsToPhi, PToV>;

    template class BSZ2015FormFactors<BToD, PToP>;
    template class BSZ2015FormFactors<BToEta, PToP>;
    template class BSZ2015FormFactors<BToEtaPrime, PToP>;
    template class BSZ2015FormFactors<BToK, PToP>;
    template class BSZ2015FormFactors<BToPi, PToP>;
    template class BSZ2015FormFactors<BsToDs, PToP>;
    template class BSZ2015FormFactors<BsToEta, PToP>;
    template class BSZ2015FormFactors<BsToEtaPrime, PToP>;
    template class BSZ2015FormFactors<BsToK, PToP>;
    template class BSZ2015FormFactors<DToEta, PToP>;
    template class BSZ2015FormFactors<DToEtaPrime, PToP>;
    template class BSZ2015FormFactors<DToK, PToP>;
    template class BSZ2015FormFactors<DToPi, PToP>;
    template class BSZ2015FormFactors<DsToEta, PToP>;
    template class BSZ2015FormFactors<DsToEtaPrime, PToP>;
    template class BSZ2015FormFactors<DsToK, PToP>;
}
