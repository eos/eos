/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <src/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <src/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil.hh>

#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace eos;

int
main(int, char **)
{
    try
    {
        Options options;
        options.set("form-factors", "BZ2004");
        BToKstarDilepton<LowRecoil> decay(Parameters::Defaults(), options);

        std::cout << "#s(GeV^2) A_0^L A_0^R A_pp^L A_pp^R A_pa^L A_pa^R" << std::endl;

        const unsigned N = 20;
        for (unsigned i = 0 ; i <= N ; ++i)
        {
            const double s_low = 14.0;
            const double s_high = 19.21;
            double s = s_low + i * (s_high - s_low) / N;

            std::cout
                << std::setprecision(5)
                << std::scientific
                << s << "\t"
                << real(decay.a_long(left_handed, s)) << "\t"
                << imag(decay.a_long(left_handed, s)) << "\t"
                << real(decay.a_long(right_handed,s)) << "\t"
                << imag(decay.a_long(right_handed,s)) << "\t"
                << real(decay.a_perp(left_handed, s)) << "\t"
                << imag(decay.a_perp(left_handed, s)) << "\t"
                << real(decay.a_perp(right_handed,s)) << "\t"
                << imag(decay.a_perp(right_handed,s)) << "\t"
                << real(decay.a_par(left_handed,  s)) << "\t"
                << imag(decay.a_par(left_handed,  s)) << "\t"
                << real(decay.a_par(right_handed, s)) << "\t"
                << imag(decay.a_par(right_handed, s)) << "\t"
                << std::endl;
        }
    }
    catch(...)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
