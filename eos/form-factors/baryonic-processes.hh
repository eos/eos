/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_BARYONIC_PROCESSES_HH
#define EOS_GUARD_EOS_FORM_FACTORS_BARYONIC_PROCESSES_HH 1

namespace eos
{
    /*
     * J=1/2^+ -> J=1/2^+ transitions
     */

    struct LambdaBToLambda {
        static constexpr const char * label = "Lambda_b->Lambda";
        // initial state mass
        static constexpr double m1 = 5.61951;
        // final state mass
        static constexpr double m2 = 1.115683;
        // semileptonic kinematic endpoint
        static constexpr double tm = (m1 - m2) * (m1 - m2);
        // pair production threshold: B + K
        static constexpr double tp = (5.279 + 0.494) * (5.279 + 0.494);
        // first resonances sorted by spin/parity
        static constexpr double mR2_0m = 5.367 * 5.367;
        static constexpr double mR2_0p = 5.711 * 5.711;
        static constexpr double mR2_1m = 5.416 * 5.416;
        static constexpr double mR2_1p = 5.750 * 5.750;
    };

    struct LambdaBToLambdaC {
        static constexpr const char * label = "Lambda_b->Lambda_c";
        // initial state mass
        static constexpr double m1 = 5.61951;
        // final state mass
        static constexpr double m2 = 2.2865;
        // semileptonic kinematic endpoint
        static constexpr double tm = (m1 - m2) * (m1 - m2);
        // first resonances sorted by spin/parity (from DLM2015 table VII)
        static constexpr double mBc = 6.276;
        static constexpr double mR2_0m = mBc * mBc;
        static constexpr double mR2_0p = (mBc + 0.449) * (mBc + 0.449);
        static constexpr double mR2_1m = (mBc + 0.056) * (mBc + 0.056);
        static constexpr double mR2_1p = (mBc + 0.492) * (mBc + 0.492);
        // see cf. DKMR2017 for the def. of t_+'s; FF specific
        static constexpr double tp_0m  = mR2_0m;
        static constexpr double tp_0p  = mR2_0p;
        static constexpr double tp_1m  = mR2_1m;
        static constexpr double tp_1p  = mR2_1p;
    };
}

#endif
