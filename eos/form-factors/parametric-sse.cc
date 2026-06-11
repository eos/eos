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

#include <eos/form-factors/parametric-sse-impl-p-to-p.hh>
#include <eos/form-factors/parametric-sse-impl-p-to-v.hh>

namespace eos
{
    template class SSEFormFactors<BToDstar, PToV>;
    template class SSEFormFactors<BToKstar, PToV>;
    template class SSEFormFactors<BToOmega, PToV>;
    template class SSEFormFactors<BToRho, PToV>;
    template class SSEFormFactors<BcToJpsi, PToV>;
    template class SSEFormFactors<BsToDsstar, PToV>;
    template class SSEFormFactors<BsToKstar, PToV>;
    template class SSEFormFactors<BsToPhi, PToV>;

    template class SSEFormFactors<BToD, PToP>;
    template class SSEFormFactors<BToEta, PToP>;
    template class SSEFormFactors<BToEtaPrime, PToP>;
    template class SSEFormFactors<BToK, PToP>;
    template class SSEFormFactors<BToPi, PToP>;
    template class SSEFormFactors<BsToDs, PToP>;
    template class SSEFormFactors<BsToEta, PToP>;
    template class SSEFormFactors<BsToEtaPrime, PToP>;
    template class SSEFormFactors<BsToK, PToP>;
    template class SSEFormFactors<DToEta, PToP>;
    template class SSEFormFactors<DToEtaPrime, PToP>;
    template class SSEFormFactors<DToK, PToP>;
    template class SSEFormFactors<DToPi, PToP>;
    template class SSEFormFactors<DsToEta, PToP>;
    template class SSEFormFactors<DsToEtaPrime, PToP>;
    template class SSEFormFactors<DsToK, PToP>;
}
