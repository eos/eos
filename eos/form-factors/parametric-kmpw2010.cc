/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013-2016, 2018 Danny van Dyk
 * Copyright (c) 2015 Christoph Bobeth
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

#include <eos/form-factors/parametric-kmpw2010-impl.hh>

namespace eos
{
    // mass B_d, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToV>::_m_B       = 5.2795;
    // mass K^*0, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToV>::_m_Kstar   = 0.89594;
    // mass B_s (0-), cf. [KMPW:2010A]
    const double KMPW2010FormFactors<PToV>::_m_Bs2_0m  = 5.366 * 5.366;
    // mass B_s (1-), cf. [KMPW:2010A]
    const double KMPW2010FormFactors<PToV>::_m_Bs2_1m  = 5.412 * 5.412;
    // mass B_s (1+), cf. [KMPW:2010A]
    const double KMPW2010FormFactors<PToV>::_m_Bs2_1p  = 5.829 * 5.829;
    const double KMPW2010FormFactors<PToV>::_tau_p     = (_m_B + _m_Kstar) * (_m_B + _m_Kstar);
    const double KMPW2010FormFactors<PToV>::_tau_m     = (_m_B - _m_Kstar) * (_m_B - _m_Kstar);
    const double KMPW2010FormFactors<PToV>::_tau_0     = _tau_p - std::sqrt(_tau_p * _tau_p - _tau_m * _tau_p);
    template class KMPW2010FormFactors<PToV>;

    // mass B_u, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToP>::_m_B      = 5.27917;
    // mass K^+, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToP>::_m_K      = 0.493677;
    // mass B_s^* (1-), cf. [KMPW:2010A]
    const double KMPW2010FormFactors<PToP>::_m_Bs2    = 5.412 * 5.412;
    const double KMPW2010FormFactors<PToP>::_tau_p    = (_m_B + _m_K) * (_m_B + _m_K);
    const double KMPW2010FormFactors<PToP>::_tau_m    = (_m_B - _m_K) * (_m_B - _m_K);
    const double KMPW2010FormFactors<PToP>::_tau_0    = _tau_p - std::sqrt(_tau_p * _tau_p - _tau_m * _tau_p);
    template class KMPW2010FormFactors<PToP>;
}
