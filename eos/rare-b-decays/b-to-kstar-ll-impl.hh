/*
 * Copyright (c) 2021, Méril Reboud
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_IMPL_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_IMPL_HH 1

#include <eos/observable.hh>
#include <eos/rare-b-decays/b-to-kstar-ll.hh>

namespace eos
{
    /*!
     * Amplitudes for the decay B -> K^* l lbar.
     */
    struct BToKstarDilepton::Amplitudes
    {
        complex<double> a_long_right, a_long_left;
        complex<double> a_perp_right, a_perp_left;
        complex<double> a_para_right, a_para_left;
        complex<double> a_time, a_scal;
        complex<double> a_para_perp, a_time_long;
        complex<double> a_time_perp, a_long_perp;
        complex<double> a_time_para, a_long_para;
    };

    struct BToKstarDilepton::AngularCoefficients
    {
        double j1s, j1c;
        double j2s, j2c;
        double j3;
        double j4;
        double j5;
        double j6s, j6c;
        double j7;
        double j8;
        double j9;

        AngularCoefficients()
        {
        }

        AngularCoefficients(const std::array<double, 12> & a) :
            j1s(a[0]),
            j1c(a[1]),
            j2s(a[2]),
            j2c(a[3]),
            j3(a[4]),
            j4(a[5]),
            j5(a[6]),
            j6s(a[7]),
            j6c(a[8]),
            j7(a[9]),
            j8(a[10]),
            j9(a[11])
        {
        }
    };

    class BToKstarDilepton::IntermediateResult :
        public CacheableObservable::IntermediateResult
    {
        public:
            BToKstarDilepton::AngularCoefficients ac;

            IntermediateResult()
            {
            }

            ~IntermediateResult() = default;
    };
}

#endif
