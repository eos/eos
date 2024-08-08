/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Viktor Kuschke
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

#include <eos/maths/power-of.hh>
#include <eos/maths/polylog.hh>
#include <eos/maths/multiplepolylog-li22.hh>

#include <eos/nonlocal-form-factors/charm-loops-impl.hh>

#include <eos/utils/exception.hh>
#include <eos/utils/log.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <complex>

#include <boost/predef.h>

#if BOOST_COMP_GNUC
#  pragma GCC optimize("no-var-tracking")
#endif

namespace eos
{
    using std::complex;
    using std::log;
    using std::sqrt;
    using std::real;

    namespace agv_2019a
    {
        complex<double> f19d(const CharmLoopsParameters & clp)
        {
            // f19d = - f29d / (2 * N_c)
            return - f29d(clp) / (6.0);
        }

        complex<double> f29d(const CharmLoopsParameters & clp)
        {
            const double lnmuhat = log(clp.muhat);
            const complex<double> xd = clp.xd;
            const complex<double> yd = clp.yd;

            const complex<double> xdinv = 1.0 / xd;
            const complex<double> ydinv = 1.0 / yd;

            const double imydinv = imag(ydinv);

            // Powers of xd and yd

            const complex<double> xd2 = power_of<2>(xd);
            const complex<double> xd3 = power_of<3>(xd);
            const complex<double> xd4 = power_of<4>(xd);
            const complex<double> xd5 = power_of<5>(xd);
            const complex<double> xd6 = power_of<6>(xd);
            const complex<double> xd7 = power_of<7>(xd);
            const complex<double> xd8 = power_of<8>(xd);
            const complex<double> xd9 = power_of<9>(xd);
            const complex<double> xd10 = power_of<10>(xd);
            const complex<double> xd11 = power_of<11>(xd);
            const complex<double> xd12 = power_of<12>(xd);
            const complex<double> xd13 = power_of<13>(xd);
            const complex<double> xd14 = power_of<14>(xd);
            const complex<double> xd15 = power_of<15>(xd);
            const complex<double> xd16 = power_of<16>(xd);
            const complex<double> xd17 = power_of<17>(xd);
            const complex<double> xd18 = power_of<18>(xd);
            const complex<double> xd19 = power_of<19>(xd);
            const complex<double> xd20 = power_of<20>(xd);
            const complex<double> xd21 = power_of<21>(xd);
            const complex<double> xd22 = power_of<22>(xd);
            const complex<double> xd24 = power_of<24>(xd);

            const complex<double> yd2 = power_of<2>(yd);
            const complex<double> yd3 = power_of<3>(yd);
            const complex<double> yd4 = power_of<4>(yd);
            const complex<double> yd5 = power_of<5>(yd);
            const complex<double> yd6 = power_of<6>(yd);
            const complex<double> yd7 = power_of<7>(yd);
            const complex<double> yd8 = power_of<8>(yd);
            const complex<double> yd9 = power_of<9>(yd);
            const complex<double> yd10 = power_of<10>(yd);
            const complex<double> yd11 = power_of<11>(yd);

            // weights appearing in the GPLs [AGV:2019A] p. 34

            const complex<double> w4 = (2.0 * xd) / power_of<2>(1.0 - xd);
            const complex<double> w5 = (2.0 * xd) / power_of<2>(1.0 + xd);
            const complex<double> w7 = (8.0 * xd2) / (1.0 - 6.0 * xd2 + xd4);

            const complex<double> w4inv = 1.0 / w4;
            const complex<double> w5inv = 1.0 / w5;
            const complex<double> w7inv = 1.0 / w7;

            // Logs of xd and yd

            const complex<double> ln1pyd = log(1.0 + yd);
            const complex<double> ln1myd = log(1.0 - yd);

            const complex<double> num1 = 32.0 * (yd8 * (-1.0 + 3.0 * yd2) + xd20 * yd8 * (-1.0 + 3.0 * yd2) - 6.0 * xd2 * yd6 * (2.0 - 9.0 * yd2 + 9.0 * yd4)
                 - 6.0 * xd18 * yd6 * (2.0 - 9.0 * yd2 + 9.0 * yd4) + xd4 * yd4 * (-40.0 + 344.0 * yd2 - 805.0 * yd4 + 335.0 * yd6)
                 + xd16 * yd4 * (-40.0 + 344.0 * yd2 - 805.0 * yd4 + 335.0 * yd6) - 8.0 * xd6 * yd2 * (8.0 - 186.0 * yd2 + 632.0 * yd4 - 479.0 * yd6
                   + 107.0 * yd8) - 8.0 * xd14 * yd2 * (8.0 - 186.0 * yd2 + 632.0 * yd4 - 479.0 * yd6 + 107.0 * yd8)
                 + 2.0 * xd8 * yd2 * (1408.0 - 5548.0 * yd2 + 5140.0 * yd4 - 2093.0 * yd6 + 407.0 * yd8)
                 + 2.0 * xd12 * yd2 * (1408.0 - 5548.0 * yd2 + 5140.0 * yd4 - 2093.0 * yd6 + 407.0 * yd8)
                 - 4.0 * xd10 * (-512.0 + 2144.0 * yd2 - 2904.0 * yd4 + 3418.0 * yd6 - 2217.0 * yd8 + 505.0 * yd10));
            const complex<double> num2 = 32.0 * power_of<2>(xd + xd3) * (-1.0 + yd) * power_of<2>(1.0 + yd) * (1.0 - yd + 9.0 * yd2 + 7.0 * yd3 + xd4 * (1.0 - yd + 9.0 * yd2 + 7.0 * yd3)
                 + xd2 * (26.0 - 26.0 * yd - 6.0 * yd2 + 38.0 * yd3));
            const complex<double> num3 = 2.0 * (yd6 * (-1.0 + yd2) + xd24 * yd6 * (-1.0 + yd2) + 32.0 * 1.0i * xd3 * yd4 * power_of<2>(-1.0 + yd2) - 32.0 * 1.0i * xd21 * yd4 * power_of<2>(-1.0 + yd2)
                 + 1920.0 * 1.0i * xd7 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2) - 1920.0 * 1.0i * xd17 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2)
                 - 32.0 * 1.0i * xd5 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) + 32.0 * 1.0i * xd19 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2)
                 - 4.0 * xd2 * yd4 * (4.0 - 5.0 * yd2 + 17.0 * yd4) - 4.0 * xd22 * yd4 * (4.0 - 5.0 * yd2 + 17.0 * yd4)
                 + 128.0 * 1.0i * xd9 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4) - 128.0 * 1.0i * xd15 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4)
                 - 64.0 * 1.0i * xd11 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4) + 64.0 * 1.0i * xd13 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4)
                 + 2.0 * xd4 * yd2 * (-40.0 + 232.0 * yd2 - 633.0 * yd4 + 633.0 * yd6)
                 + 2.0 * xd20 * yd2 * (-40.0 + 232.0 * yd2 - 633.0 * yd4 + 633.0 * yd6) + 8.0 * xd10 * (-3312.0 + 26640.0 * yd2 - 20660.0 * yd4 + 1685.0 * yd6
                   + 1135.0 * yd8) + 8.0 * xd14 * (-3312.0 + 26640.0 * yd2 - 20660.0 * yd4 + 1685.0 * yd6 + 1135.0 * yd8)
                 - 4.0 * xd6 * (32.0 - 736.0 * yd2 + 1940.0 * yd4 - 4817.0 * yd6 + 2253.0 * yd8)
                 - 4.0 * xd18 * (32.0 - 736.0 * yd2 + 1940.0 * yd4 - 4817.0 * yd6 + 2253.0 * yd8)
                 - 4.0 * xd12 * (16000.0 - 118920.0 * yd2 + 142408.0 * yd4 - 64785.0 * yd6 + 11409.0 * yd8)
                 + xd8 * (5376.0 - 21696.0 * yd2 + 111296.0 * yd4 - 95535.0 * yd6 + 21551.0 * yd8)
                 + xd16 * (5376.0 - 21696.0 * yd2 + 111296.0 * yd4 - 95535.0 * yd6 + 21551.0 * yd8));
            const complex<double> num4 = 2.0 * (yd6 * (-1.0 + yd2) + xd24 * yd6 * (-1.0 + yd2) - 32.0 * 1.0i * xd3 * yd4 * power_of<2>(-1.0 + yd2) + 32.0 * 1.0i * xd21 * yd4 * power_of<2>(-1.0 + yd2)
                 - 1920.0 * 1.0i * xd7 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2) + 1920.0 * 1.0i * xd17 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2)
                 + 32.0 * 1.0i * xd5 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) - 32.0 * 1.0i * xd19 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2)
                 - 4.0 * xd2 * yd4 * (4.0 - 5.0 * yd2 + 17.0 * yd4) - 4.0 * xd22 * yd4 * (4.0 - 5.0 * yd2 + 17.0 * yd4)
                 - 128.0 * 1.0i * xd9 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4) + 128.0 * 1.0i * xd15 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4)
                 + 64.0 * 1.0i * xd11 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4) - 64.0 * 1.0i * xd13 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4)
                 + 2.0 * xd4 * yd2 * (-40.0 + 232.0 * yd2 - 633.0 * yd4 + 633.0 * yd6)
                 + 2.0 * xd20 * yd2 * (-40.0 + 232.0 * yd2 - 633.0 * yd4 + 633.0 * yd6) + 8.0 * xd10 * (-3312.0 + 26640.0 * yd2 - 20660.0 * yd4 + 1685.0 * yd6
                   + 1135.0 * yd8) + 8.0 * xd14 * (-3312.0 + 26640.0 * yd2 - 20660.0 * yd4 + 1685.0 * yd6 + 1135.0 * yd8)
                 - 4.0 * xd6 * (32.0 - 736.0 * yd2 + 1940.0 * yd4 - 4817.0 * yd6 + 2253.0 * yd8)
                 - 4.0 * xd18 * (32.0 - 736.0 * yd2 + 1940.0 * yd4 - 4817.0 * yd6 + 2253.0 * yd8)
                 - 4.0 * xd12 * (16000.0 - 118920.0 * yd2 + 142408.0 * yd4 - 64785.0 * yd6 + 11409.0 * yd8)
                 + xd8 * (5376.0 - 21696.0 * yd2 + 111296.0 * yd4 - 95535.0 * yd6 + 21551.0 * yd8)
                 + xd16 * (5376.0 - 21696.0 * yd2 + 111296.0 * yd4 - 95535.0 * yd6 + 21551.0 * yd8));
            const complex<double> num5 = 32.0 * xd2 * power_of<2>(-1.0i + xd) * power_of<2>(1.0i + xd) * power_of<2>(-1.0 + yd) * (1.0 + yd) * (-1.0 - yd - 9.0 * yd2 + 7.0 * yd3
                 + xd4 * (-1.0 - yd - 9.0 * yd2 + 7.0 * yd3) + xd2 * (-26.0 - 26.0 * yd + 6.0 * yd2 + 38.0 * yd3));
            const complex<double> num6 = 2.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * power_of<2>(1.0 + yd) * ((-1.0 + yd) * yd2 + xd8 * (-1.0 + yd) * yd2
                 + xd4 * (-232.0 + 232.0 * yd + 2.0 * yd2 - 386.0 * yd3) - 4.0 * xd2 * (5.0 - 5.0 * yd + 24.0 * yd2 + 24.0 * yd3)
                 - 4.0 * xd6 * (5.0 - 5.0 * yd + 24.0 * yd2 + 24.0 * yd3));
            const complex<double> num7 = 256.0 * xd2 * power_of<4>(1.0 + xd2) * (-1.0 + yd) * yd3 * (1.0 + yd);
            const complex<double> num8 = 16.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 * (3.0 + yd2) + xd8 * yd2 * (3.0 + yd2)
                 + xd4 * (-408.0 + 706.0 * yd2 - 402.0 * yd4) + 4.0 * xd2 * (5.0 - 23.0 * yd2 + 6.0 * yd4) + 4.0 * xd6 * (5.0 - 23.0 * yd2 + 6.0 * yd4));
            const complex<double> num9 = 16.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 * (-1.0 + 5.0 * yd2) + xd8 * yd2 * (-1.0 + 5.0 * yd2)
                 - 12.0 * xd2 * (1.0 + yd2 + 2.0 * yd4) - 12.0 * xd6 * (1.0 + yd2 + 2.0 * yd4) - 2.0 * xd4 * (44.0 - 53.0 * yd2 + 61.0 * yd4));
            const complex<double> num10 = 2.0 * power_of<2>(1.0 + xd2) * power_of<2>(-1.0 + yd) * (1.0 + yd) * (yd2 * (1.0 + yd) + xd8 * yd2 * (1.0 + yd)
                 + xd2 * (20.0 + 20.0 * yd + 96.0 * yd2 - 96.0 * yd3) + xd6 * (20.0 + 20.0 * yd + 96.0 * yd2 - 96.0 * yd3)
                 - 2.0 * xd4 * (-116.0 - 116.0 * yd + yd2 + 193.0 * yd3));
            const complex<double> num11 = (1.0 + 4.0 * xd + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd4 - 4.0 * xd * yd4 - 4.0 * xd17 * yd4 + xd18 * yd4
                 + 16.0 * xd3 * yd2 * (-3.0 + 5.0 * yd2) + 16.0 * xd15 * yd2 * (-3.0 + 5.0 * yd2) + xd2 * yd2 * (12.0 + 13.0 * yd2)
                 + xd16 * yd2 * (12.0 + 13.0 * yd2) + xd8 * (-6848.0 + 2980.0 * yd2 - 566.0 * yd4) + xd10 * (-6848.0 + 2980.0 * yd2 - 566.0 * yd4)
                 - 16.0 * xd5 * (8.0 - 46.0 * yd2 + 21.0 * yd4) - 16.0 * xd13 * (8.0 - 46.0 * yd2 + 21.0 * yd4) - 4.0 * xd4 * (-8.0 + 43.0 * yd2 + 56.0 * yd4)
                 - 4.0 * xd14 * (-8.0 + 43.0 * yd2 + 56.0 * yd4) - 16.0 * xd7 * (-224.0 + 45.0 * yd2 + 133.0 * yd4) - 16.0 * xd11 * (-224.0 + 45.0 * yd2 + 133.0 * yd4)
                 + 4.0 * xd6 * (-216.0 - 513.0 * yd2 + 194.0 * yd4) + 4.0 * xd12 * (-216.0 - 513.0 * yd2 + 194.0 * yd4)
                 + 8.0 * xd9 * (2976.0 - 4472.0 * yd2 + 1621.0 * yd4));
            const complex<double> num12 = (1.0 - 4.0 * xd + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd4 + 4.0 * xd * yd4 + 4.0 * xd17 * yd4 + xd18 * yd4
                 + xd2 * yd2 * (12.0 + 13.0 * yd2) + xd16 * yd2 * (12.0 + 13.0 * yd2) + xd8 * (-6848.0 + 2980.0 * yd2 - 566.0 * yd4)
                 + xd10 * (-6848.0 + 2980.0 * yd2 - 566.0 * yd4) + xd3 * (48.0 * yd2 - 80.0 * yd4) + xd15 * (48.0 * yd2 - 80.0 * yd4)
                 + 16.0 * xd5 * (8.0 - 46.0 * yd2 + 21.0 * yd4) + 16.0 * xd13 * (8.0 - 46.0 * yd2 + 21.0 * yd4) - 4.0 * xd4 * (-8.0 + 43.0 * yd2 + 56.0 * yd4)
                 - 4.0 * xd14 * (-8.0 + 43.0 * yd2 + 56.0 * yd4) + 16.0 * xd7 * (-224.0 + 45.0 * yd2 + 133.0 * yd4) + 16.0 * xd11 * (-224.0 + 45.0 * yd2 + 133.0 * yd4)
                 + 4.0 * xd6 * (-216.0 - 513.0 * yd2 + 194.0 * yd4) + 4.0 * xd12 * (-216.0 - 513.0 * yd2 + 194.0 * yd4)
                 - 8.0 * xd9 * (2976.0 - 4472.0 * yd2 + 1621.0 * yd4));
            const complex<double> num13 = 4.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * power_of<2>(1.0 + yd) * ((-1.0 + yd) * yd2 + xd8 * (-1.0 + yd) * yd2
                 + 4.0 * xd2 * (-1.0 + yd + 12.0 * yd2 + 4.0 * yd3) + 4.0 * xd6 * (-1.0 + yd + 12.0 * yd2 + 4.0 * yd3)
                 + 2.0 * xd4 * (92.0 - 92.0 * yd - 47.0 * yd2 + 111.0 * yd3));
            const complex<double> num14 = 4.0 * power_of<2>(1.0 + xd2) * power_of<2>(-1.0 + yd) * (1.0 + yd) * (yd2 * (1.0 + yd) + xd8 * yd2 * (1.0 + yd)
                 + 4.0 * xd2 * (1.0 + yd - 12.0 * yd2 + 4.0 * yd3) + 4.0 * xd6 * (1.0 + yd - 12.0 * yd2 + 4.0 * yd3)
                 + 2.0 * xd4 * (-92.0 - 92.0 * yd + 47.0 * yd2 + 111.0 * yd3));
            const complex<double> num15 = 8.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 + 3.0 * yd4 + xd2 * (4.0 - 52.0 * yd2) + xd6 * (4.0 - 52.0 * yd2)
                 + xd4 * (-248.0 + 406.0 * yd2 - 262.0 * yd4) + xd8 * (yd2 + 3.0 * yd4));
            const complex<double> num16 = 32.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (-2.0 * 1.0i * xd * (-1.0 + yd) + yd + xd2 * yd) * (2.0 * 1.0i * xd * (-1.0 + yd) + yd + xd2 * yd)
                 * (yd + xd2 * yd - 2.0 * xd * (1.0 + yd)) * (yd + xd2 * yd + 2.0 * xd * (1.0 + yd));
            const complex<double> num17 = 32.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (-2.0 * xd * (-1.0 + yd) + yd + xd2 * yd) * (2.0 * xd * (-1.0 + yd) + yd + xd2 * yd)
                 * (yd + xd2 * yd - 2.0 * 1.0i * xd * (1.0 + yd)) * (yd + xd2 * yd + 2.0 * 1.0i * xd * (1.0 + yd));
            const complex<double> num18 = 2.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 + 7.0 * yd4 + xd4 * (-312.0 + 534.0 * yd2 + 128.0 * yd3 - 302.0 * yd4)
                 + 4.0 * xd2 * (1.0 - 13.0 * yd2 + 16.0 * yd3 + 4.0 * yd4) + 4.0 * xd6 * (1.0 - 13.0 * yd2 + 16.0 * yd3 + 4.0 * yd4) + xd8 * (yd2 + 7.0 * yd4));
            const complex<double> num19 = 2.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 + 7.0 * yd4 + 4.0 * xd2 * (1.0 - 13.0 * yd2 - 16.0 * yd3 + 4.0 * yd4)
                 + 4.0 * xd6 * (1.0 - 13.0 * yd2 - 16.0 * yd3 + 4.0 * yd4) + xd8 * (yd2 + 7.0 * yd4) - 2.0 * xd4 * (156.0 - 267.0 * yd2 + 64.0 * yd3 + 151.0 * yd4));
            const complex<double> num20 = 128.0 * xd2 * power_of<4>(1.0 + xd2) * power_of<2>(-1.0 + yd) * yd2 * power_of<2>(1.0 + yd);

            const complex<double> denom1 = 9.0 * (xd2 * (8.0 - 6.0 * yd) + yd + xd4 * yd) * (yd + xd4 * yd - 2.0 * xd2 * (4.0 + 3.0 * yd)) * power_of<3>(yd3 + xd4 * yd3 - 2.0 * xd2 * yd * (-2.0 + yd2));
            const complex<double> denom2 = 9.0 * power_of<3>(yd2 + xd4 * yd2 - 2.0 * xd2 * (-2.0 + yd2));
            const complex<double> denom3 = 9.0 * xd2 * (xd2 * (8.0 - 6.0 * yd) + yd + xd4 * yd) * (yd + xd4 * yd - 2.0 * xd2 * (4.0 + 3.0 * yd)) * power_of<3>(yd2 + xd4 * yd2 - 2.0 * xd2 * (-2.0 + yd2));
            const complex<double> denom4 = 9.0 * xd2 * (xd2 * (8.0 - 6.0 * yd) + yd + xd4 * yd) * (yd + xd4 * yd - 2.0 * xd2 * (4.0 + 3.0 * yd)) * power_of<2>(yd2 + xd4 * yd2 - 2.0 * xd2 * (-2.0 + yd2));
            const complex<double> denom5 = 9.0 * power_of<2>(yd3 + xd4 * yd3 - 2.0 * xd2 * yd * (-2.0 + yd2));
            const complex<double> denom6 = 9.0 * xd2 * yd2 * (xd2 * (8.0 - 6.0 * yd) + yd + xd4 * yd) * (yd + xd4 * yd - 2.0 * xd2 * (4.0 + 3.0 * yd)) * power_of<3>(yd2 + xd4 * yd2 - 2.0 * xd2 * (-2.0 + yd2));
            const complex<double> denom7 = 9.0 * (xd2 * (8.0 - 6.0 * yd) + yd + xd4 * yd) * power_of<3>(yd3 + xd4 * yd3 - 2.0 * xd2 * yd * (-2.0 + yd2));
            const complex<double> denom8 = 9.0 * (yd + xd4 * yd - 2.0 * xd2 * (4.0 + 3.0 * yd)) * power_of<3>(yd3 + xd4 * yd3 - 2.0 * xd2 * yd * (-2.0 + yd2));

            const complex<double> term1 = (2.0 * (yd4 * (-6.0 * (23.0 + 24.0 * ln2) + yd2 * (147.0 + 2.0 * pisqu + 384.0 * ln2 + 384.0 * ln2squ))
                    + xd8 * yd4 * (-6.0 * (23.0 + 24.0 * ln2) + yd2 * (147.0 + 2.0 * pisqu + 384.0 * ln2 + 384.0 * ln2squ))
                    - 4.0 * xd2 * yd2 * (6.0 * (47.0 + 48.0 * ln2) - 4.0 * yd2 * (114.0 + pisqu + 264.0 * ln2 + 192.0 * ln2squ)
                        + yd4 * (165.0 + 2.0 * pisqu + 528.0 * ln2 + 384.0 * ln2squ)) - 4.0 * xd6 * yd2 * (6.0 * (47.0 + 48.0 * ln2)
                        - 4.0 * yd2 * (114.0 + pisqu + 264.0 * ln2 + 192.0 * ln2squ) + yd4 * (165.0 + 2.0 * pisqu + 528.0 * ln2 + 384.0 * ln2squ))
                    + 2.0 * xd4 * (16.0 * yd2 * (156.0 + pisqu + 216.0 * ln2 + 192.0 * ln2squ) + 3.0 * yd6 * (171.0 + 2.0 * pisqu + 192.0 * ln2 + 384.0 * ln2squ)
                        - 2.0 * yd4 * (915.0 + 8.0 * pisqu + 1272.0 * ln2 + 1536.0 * ln2squ) - 384.0 * (3.0 + ln4)))) / (27.0 * power_of<2>(yd3 + xd4 * yd3 - 2.0 * xd2 * yd * (-2.0 + yd2)));

            const complex<double> logs1 = 32.0 * power_of<2>(xd + xd3) * (-1.0 + yd) * power_of<2>(1.0 + yd) * (-1.0 + yd + yd3 * (-7.0 + 32.0 * ln2) - yd2 * (9.0 + 32.0 * ln2) + xd4 * (-1.0 + yd + yd3 * (-7.0 + 32.0 * ln2) - yd2 * (9.0 + 32.0 * ln2))
                + xd2 * (-26.0 + 26.0 * yd + yd2 * (6.0 - 64.0 * ln2) + yd3 * (-38.0 + 64.0 * ln2)));
            const complex<double> logs2 = 32.0 * xd2 * power_of<2>(-1.0i + xd) * power_of<2>(1.0i + xd) * power_of<2>(-1.0 + yd) * (1.0 + yd) * (1.0 + yd + yd3 * (-7.0 + 32.0 * ln2) + yd2 * (9.0 + 32.0 * ln2) + xd4 * (1.0 + yd + yd3 * (-7.0 + 32.0 * ln2) + yd2 * (9.0 + 32.0 * ln2))
                + xd2 * (26.0 + 26.0 * yd + yd3 * (-38.0 + 64.0 * ln2) + yd2 * (-6.0 + 64.0 * ln2)));
            const complex<double> logs3 = 16.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 * (-1.0 + 5.0 * yd2) + xd8 * yd2 * (-1.0 + 5.0 * yd2) + 4.0 * xd2 * (-3.0 + yd4 * (-6.0 + 32.0 * ln2) - yd2 * (3.0 + 32.0 * ln2))
                + 4.0 * xd6 * (-3.0 + yd4 * (-6.0 + 32.0 * ln2) - yd2 * (3.0 + 32.0 * ln2)) + 2.0 * xd4 * (-44.0 + yd2 * (53.0 - 128.0 * ln2) + yd4 * (-61.0 + 128.0 * ln2)));
            const complex<double> logs4 = 4.0 * power_of<4>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 * (-1.0 + yd2) + xd4 * yd2 * (-1.0 + yd2) + 2.0 * xd2 * (-6.0 - 64.0 * yd3 + 3.0 * yd4 * (-7.0 + 64.0 * ln2) - yd2 * (5.0 + 192.0 * ln2)));
            const complex<double> logs5 = 4.0 * power_of<4>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 * (-1.0 + yd2) + xd4 * yd2 * (-1.0 + yd2) + 2.0 * xd2 * (-6.0 + 64.0 * yd3 + 3.0 * yd4 * (-7.0 + 64.0 * ln2) - yd2 * (5.0 + 192.0 * ln2)));
            const complex<double> logs6 = 16.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 * (3.0 + yd2) + xd8 * yd2 * (3.0 + yd2) + 4.0 * xd2 * (5.0 + yd4 * (6.0 + 32.0 * ln2) - yd2 * (23.0 + 32.0 * ln2)) + 4.0 * xd6 * (5.0 + yd4 * (6.0 + 32.0 * ln2) - yd2 * (23.0 + 32.0 * ln2))
                + xd4 * (-408.0 + yd2 * (706.0 - 256.0 * ln2) + yd4 * (-402.0 + 256.0 * ln2)));
            const complex<double> logs7 = 16.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 + 3.0 * yd4 + xd8 * (yd2 + 3.0 * yd4) + 4.0 * xd2 * (1.0 + 32.0 * yd4 * ln2 - yd2 * (13.0 + 32.0 * ln2)) + 4.0 * xd6 * (1.0 + 32.0 * yd4 * ln2 - yd2 * (13.0 + 32.0 * ln2))
                + xd4 * (-248.0 + yd2 * (406.0 - 256.0 * ln2) + yd4 * (-262.0 + 256.0 * ln2)));
            const complex<double> logs8 = 2.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 * (-1.0 + 17.0 * yd2) + xd8 * yd2 * (-1.0 + 17.0 * yd2) + 4.0 * xd2 * (-1.0 - 80.0 * yd3 + yd2 * (13.0 - 128.0 * ln2) + 4.0 * yd4 * (5.0 + 32.0 * ln2))
                + 4.0 * xd6 * (-1.0 - 80.0 * yd3 + yd2 * (13.0 - 128.0 * ln2) + 4.0 * yd4 * (5.0 + 32.0 * ln2)) + 2.0 * xd4 * (-36.0 - 320.0 * yd3 + yd2 * (117.0 - 512.0 * ln2) + yd4 * (31.0 + 512.0 * ln2)));
            const complex<double> logs9 = 2.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 * (-1.0 + 17.0 * yd2) + xd8 * yd2 * (-1.0 + 17.0 * yd2) + 4.0 * xd2 * (-1.0 + 80.0 * yd3 + yd2 * (13.0 - 128.0 * ln2) + 4.0 * yd4 * (5.0 + 32.0 * ln2))
                + 4.0 * xd6 * (-1.0 + 80.0 * yd3 + yd2 * (13.0 - 128.0 * ln2) + 4.0 * yd4 * (5.0 + 32.0 * ln2)) + 2.0 * xd4 * (-36.0 + 320.0 * yd3 + yd2 * (117.0 - 512.0 * ln2) + yd4 * (31.0 + 512.0 * ln2)));
            const complex<double> logs10 = 8.0 * (-1.0 + yd) * (1.0 + yd) * (yd8 + xd20 * yd8 + xd2 * (8.0 * yd6 - 2.0 * yd8 * (31.0 + 192.0 * ln2)) + xd18 * (8.0 * yd6 - 2.0 * yd8 * (31.0 + 192.0 * ln2))
                + 8.0 * xd6 * yd4 * (-76.0 + yd4 * (-217.0 + 192.0 * ln2) + 12.0 * yd2 * (55.0 + 256.0 * ln2)) + 8.0 * xd14 * yd4 * (-76.0 + yd4 * (-217.0 + 192.0 * ln2) + 12.0 * yd2 * (55.0 + 256.0 * ln2))
                + xd8 * (2048.0 * yd2 + 7952.0 * yd4 + yd8 * (930.0 - 27648.0 * ln2) + 96.0 * yd6 * (-29.0 + 1024.0 * ln2)) + xd12 * (2048.0 * yd2 + 7952.0 * yd4 + yd8 * (930.0 - 27648.0 * ln2) + 96.0 * yd6 * (-29.0 + 1024.0 * ln2))
                + xd4 * yd4 * (-16.0 - 288.0 * yd2 + yd4 * (605.0 + 3072.0 * ln2)) + xd16 * yd4 * (-16.0 - 288.0 * yd2 + yd4 * (605.0 + 3072.0 * ln2))
                + xd10 * (4096.0 + 4096.0 * yd2 - 7488.0 * yd4 + 16.0 * yd6 * (1003.0 + 9216.0 * ln2) - 4.0 * yd8 * (1405.0 + 12864.0 * ln2)));
            const complex<double> logs11 = 512.0 * xd2 * power_of<4>(1.0 + xd2) * (-1.0 + yd) * yd2 * (1.0 + yd) * (-3.0 * yd - ln4 + yd2 * ln4);
            const complex<double> logs12 = 512.0 * xd2 * power_of<4>(1.0 + xd2) * (-1.0 + yd) * yd2 * (1.0 + yd) * (3.0 * yd - ln4 + yd2 * ln4);
            const complex<double> logs13 = 16.0 * (-8.0 * xd6 * yd2 * (16.0 - 372.0 * yd2 + yd8 * (214.0 - 32.0 * ln2) + 32.0 * yd * ln2 - 670.0 * yd3 * ln2 + 939.0 * yd5 * ln2 - 540.0 * yd7 * ln2 + 239.0 * yd9 * ln2 + 16.0 * yd4 * (79.0 + 32.0 * ln2) - 2.0 * yd6 * (479.0 + 240.0 * ln2))
                - 8.0 * xd14 * yd2 * (16.0 - 372.0 * yd2 + yd8 * (214.0 - 32.0 * ln2) + 32.0 * yd * ln2 - 670.0 * yd3 * ln2 + 939.0 * yd5 * ln2 - 540.0 * yd7 * ln2 + 239.0 * yd9 * ln2 + 16.0 * yd4 * (79.0 + 32.0 * ln2) - 2.0 * yd6 * (479.0 + 240.0 * ln2))
                + xd4 * yd4 * (-80.0 + 688.0 * yd2 - 288.0 * yd * ln2 + 1165.0 * yd3 * ln2 - 954.0 * yd5 * ln2 + 77.0 * yd7 * ln2 - 2.0 * yd4 * (805.0 + 256.0 * ln2) + yd6 * (670.0 + 512.0 * ln2))
                + xd16 * yd4 * (-80.0 + 688.0 * yd2 - 288.0 * yd * ln2 + 1165.0 * yd3 * ln2 - 954.0 * yd5 * ln2 + 77.0 * yd7 * ln2 - 2.0 * yd4 * (805.0 + 256.0 * ln2) + yd6 * (670.0 + 512.0 * ln2))
                + 2.0 * xd8 * yd2 * (2816.0 - 11096.0 * yd2 + yd4 * (10280.0 - 8192.0 * ln2) + yd8 * (814.0 - 2304.0 * ln2) + 5632.0 * yd * ln2 - 12528.0 * yd3 * ln2 + 16057.0 * yd5 * ln2 - 11170.0 * yd7 * ln2 + 2009.0 * yd9 * ln2
                + 2.0 * yd6 * (-2093.0 + 5248.0 * ln2)) + 2.0 * xd12 * yd2 * (2816.0 - 11096.0 * yd2 + yd4 * (10280.0 - 8192.0 * ln2) + yd8 * (814.0 - 2304.0 * ln2) + 5632.0 * yd * ln2 - 12528.0 * yd3 * ln2 + 16057.0 * yd5 * ln2 - 11170.0 * yd7 * ln2
                + 2009.0 * yd9 * ln2 + 2.0 * yd6 * (-2093.0 + 5248.0 * ln2)) + 4.0 * xd10 * (1024.0 - 4288.0 * yd2 + 5808.0 * yd4 + 5760.0 * yd3 * ln2 - 15354.0 * yd5 * ln2 + 20429.0 * yd7 * ln2 - 13836.0 * yd9 * ln2 + 3001.0 * yd11 * ln2
                - 2.0 * yd10 * (505.0 + 1072.0 * ln2) - 4.0 * yd6 * (1709.0 + 1536.0 * ln2) + yd8 * (4434.0 + 8288.0 * ln2)) + yd7 * (-2.0 * yd + 6.0 * yd3 + ln2 + yd4 * ln2 - yd2 * ln4)
                + xd20 * yd7 * (-2.0 * yd + 6.0 * yd3 + ln2 + yd4 * ln2 - yd2 * ln4) + 2.0 * xd2 * yd5 * (-12.0 * yd - 33.0 * yd2 * ln2 + 28.0 * yd4 * ln2 - 2.0 * yd5 * (27.0 + 16.0 * ln2) + yd3 * (54.0 + 32.0 * ln2) + ln4 + yd6 * log(8.0))
                + 2.0 * xd18 * yd5 * (-12.0 * yd - 33.0 * yd2 * ln2 + 28.0 * yd4 * ln2 - 2.0 * yd5 * (27.0 + 16.0 * ln2) + yd3 * (54.0 + 32.0 * ln2) + ln4 + yd6 * log(8.0)));
            const complex<double> logs14 = 16.0 * (xd4 * yd4 * (80.0 - 688.0 * yd2 - 288.0 * yd * ln2 + 1165.0 * yd3 * ln2 - 954.0 * yd5 * ln2 + 77.0 * yd7 * ln2 - 2.0 * yd6 * (335.0 + 256.0 * ln2) + 2.0 * yd4 * (805.0 + 256.0 * ln2))
                + xd16 * yd4 * (80.0 - 688.0 * yd2 - 288.0 * yd * ln2 + 1165.0 * yd3 * ln2 - 954.0 * yd5 * ln2 + 77.0 * yd7 * ln2 - 2.0 * yd6 * (335.0 + 256.0 * ln2) + 2.0 * yd4 * (805.0 + 256.0 * ln2))
                - 8.0 * xd6 * yd2 * (-16.0 + 372.0 * yd2 + 32.0 * yd * ln2 - 670.0 * yd3 * ln2 + 939.0 * yd5 * ln2 - 540.0 * yd7 * ln2 + 239.0 * yd9 * ln2 + yd8 * (-214.0 + 32.0 * ln2) - 16.0 * yd4 * (79.0 + 32.0 * ln2) + yd6 * (958.0 + 480.0 * ln2))
                - 8.0 * xd14 * yd2 * (-16.0 + 372.0 * yd2 + 32.0 * yd * ln2 - 670.0 * yd3 * ln2 + 939.0 * yd5 * ln2 - 540.0 * yd7 * ln2 + 239.0 * yd9 * ln2 + yd8 * (-214.0 + 32.0 * ln2) - 16.0 * yd4 * (79.0 + 32.0 * ln2) + yd6 * (958.0 + 480.0 * ln2))
                + 2.0 * xd8 * yd2 * (-2816.0 + 11096.0 * yd2 + yd6 * (4186.0 - 10496.0 * ln2) + 5632.0 * yd * ln2 - 12528.0 * yd3 * ln2 + 16057.0 * yd5 * ln2 - 11170.0 * yd7 * ln2 + 2009.0 * yd9 * ln2 + 8.0 * yd4 * (-1285.0 + 1024.0 * ln2)
                + yd8 * (-814.0 + 2304.0 * ln2)) + 2.0 * xd12 * yd2 * (-2816.0 + 11096.0 * yd2 + yd6 * (4186.0 - 10496.0 * ln2) + 5632.0 * yd * ln2 - 12528.0 * yd3 * ln2 + 16057.0 * yd5 * ln2 - 11170.0 * yd7 * ln2 + 2009.0 * yd9 * ln2
                + 8.0 * yd4 * (-1285.0 + 1024.0 * ln2) + yd8 * (-814.0 + 2304.0 * ln2)) + 4.0 * xd10 * (-1024.0 + 4288.0 * yd2 - 5808.0 * yd4 + 5760.0 * yd3 * ln2 - 15354.0 * yd5 * ln2 + 20429.0 * yd7 * ln2 - 13836.0 * yd9 * ln2 + 3001.0 * yd11 * ln2
                + 2.0 * yd10 * (505.0 + 1072.0 * ln2) - 2.0 * yd8 * (2217.0 + 4144.0 * ln2) + yd6 * (6836.0 + 6144.0 * ln2)) + yd7 * (2.0 * yd - 6.0 * yd3 + ln2 + yd4 * ln2 - yd2 * ln4)
                + xd20 * yd7 * (2.0 * yd - 6.0 * yd3 + ln2 + yd4 * ln2 - yd2 * ln4) + 2.0 * xd2 * yd5 * (12.0 * yd - 33.0 * yd2 * ln2 + 28.0 * yd4 * ln2 - 2.0 * yd3 * (27.0 + 16.0 * ln2) + yd5 * (54.0 + 32.0 * ln2) + ln4 + yd6 * log(8.0))
                + 2.0 * xd18 * yd5 * (12.0 * yd - 33.0 * yd2 * ln2 + 28.0 * yd4 * ln2 - 2.0 * yd3 * (27.0 + 16.0 * ln2) + yd5 * (54.0 + 32.0 * ln2) + ln4 + yd6 * log(8.0)));
            const complex<double> logs15 = 16.0 * (-4.0 * xd2 * yd2 * (6.0 - 2.0 * yd2 * (11.0 + 16.0 * ln2) + yd4 * (11.0 + 16.0 * ln2)) - 4.0 * xd6 * yd2 * (6.0 - 2.0 * yd2 * (11.0 + 16.0 * ln2) + yd4 * (11.0 + 16.0 * ln2)) + yd4 * (-3.0 + 8.0 * yd2 * (1.0 + ln4))
                + xd8 * yd4 * (-3.0 + 8.0 * yd2 * (1.0 + ln4)) + 2.0 * xd4 * (-16.0 + 8.0 * yd2 * (9.0 + 16.0 * ln2) - yd4 * (53.0 + 128.0 * ln2) + 12.0 * yd6 * (1.0 + log(16.0))));
            const complex<double> logs16 = 4.0 * (yd8 * (-1.0 + yd2) * ln2 + xd24 * yd8 * (-1.0 + yd2) * ln2 - 32.0 * 1.0i * xd3 * yd6 * power_of<2>(-1.0 + yd2) * ln2 + 32.0 * 1.0i * xd21 * yd6 * power_of<2>(-1.0 + yd2) * ln2 - 1920.0 * 1.0i * xd7 * yd4 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2) * ln2
                + 1920.0 * 1.0i * xd17 * yd4 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2) * ln2 + 32.0 * 1.0i * xd5 * yd4 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) * ln2 - 32.0 * 1.0i * xd19 * yd4 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) * ln2
                - 128.0 * 1.0i * xd9 * yd2 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4) * ln2 + 128.0 * 1.0i * xd15 * yd2 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4) * ln2 + 64.0 * 1.0i * xd11 * yd2 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4) * ln2
                - 64.0 * 1.0i * xd13 * yd2 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4) * ln2 + 2.0 * xd4 * yd4 * (-40.0 * ln2 + 8.0 * yd2 * (9.0 + 29.0 * ln2) - 3.0 * yd4 * (108.0 + 211.0 * ln2) + yd6 * (312.0 + 633.0 * ln2))
                + 2.0 * xd20 * yd4 * (-40.0 * ln2 + 8.0 * yd2 * (9.0 + 29.0 * ln2) - 3.0 * yd4 * (108.0 + 211.0 * ln2) + yd6 * (312.0 + 633.0 * ln2))
                - 4.0 * xd6 * yd2 * (32.0 * ln2 - 32.0 * yd2 * (4.0 + 23.0 * ln2) + 4.0 * yd4 * (292.0 + 485.0 * ln2) + yd8 * (1144.0 + 2253.0 * ln2) - yd6 * (2599.0 + 4817.0 * ln2))
                - 4.0 * xd18 * yd2 * (32.0 * ln2 - 32.0 * yd2 * (4.0 + 23.0 * ln2) + 4.0 * yd4 * (292.0 + 485.0 * ln2) + yd8 * (1144.0 + 2253.0 * ln2) - yd6 * (2599.0 + 4817.0 * ln2))
                + 8.0 * xd10 * yd2 * (-16.0 * (304.0 + 207.0 * ln2) + 5.0 * yd8 * (-704.0 + 227.0 * ln2) + 16.0 * yd2 * (1484.0 + 1665.0 * ln2) + yd6 * (16491.0 + 1685.0 * ln2) - 4.0 * yd4 * (7534.0 + 5165.0 * ln2))
                + 8.0 * xd14 * yd2 * (-16.0 * (304.0 + 207.0 * ln2) + 5.0 * yd8 * (-704.0 + 227.0 * ln2) + 16.0 * yd2 * (1484.0 + 1665.0 * ln2) + yd6 * (16491.0 + 1685.0 * ln2) - 4.0 * yd4 * (7534.0 + 5165.0 * ln2))
                - 4.0 * xd12 * (8192.0 + yd8 * (38412.0 - 64785.0 * ln2) + 128.0 * yd2 * (-358.0 + 125.0 * ln2) - 24.0 * yd4 * (-3168.0 + 4955.0 * ln2) + yd10 * (-8168.0 + 11409.0 * ln2) + 8.0 * yd6 * (-9123.0 + 17801.0 * ln2))
                + xd8 * yd2 * (256.0 * (2.0 + 21.0 * ln2) - 64.0 * yd2 * (280.0 + 339.0 * ln2) + 64.0 * yd4 * (1045.0 + 1739.0 * ln2) - 5.0 * yd6 * (11744.0 + 19107.0 * ln2) + yd8 * (15808.0 + 21551.0 * ln2))
                + xd16 * yd2 * (256.0 * (2.0 + 21.0 * ln2) - 64.0 * yd2 * (280.0 + 339.0 * ln2) + 64.0 * yd4 * (1045.0 + 1739.0 * ln2) - 5.0 * yd6 * (11744.0 + 19107.0 * ln2) + yd8 * (15808.0 + 21551.0 * ln2))
                - 4.0 * xd2 * yd6 * (yd4 * (8.0 + 17.0 * ln2) + log(16.0) - yd2 * (3.0 + log(32.0))) - 4.0 * xd22 * yd6 * (yd4 * (8.0 + 17.0 * ln2) + log(16.0) - yd2 * (3.0 + log(32.0))));
            const complex<double> logs17 = 4.0 * (yd8 * (-1.0 + yd2) * ln2 + xd24 * yd8 * (-1.0 + yd2) * ln2 + 32.0 * 1.0i * xd3 * yd6 * power_of<2>(-1.0 + yd2) * ln2 - 32.0 * 1.0i * xd21 * yd6 * power_of<2>(-1.0 + yd2) * ln2 + 1920.0 * 1.0i * xd7 * yd4 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2) * ln2
                - 1920.0 * 1.0i * xd17 * yd4 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2) * ln2 - 32.0 * 1.0i * xd5 * yd4 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) * ln2 + 32.0 * 1.0i * xd19 * yd4 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) * ln2
                + 128.0 * 1.0i * xd9 * yd2 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4) * ln2 - 128.0 * 1.0i * xd15 * yd2 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4) * ln2 - 64.0 * 1.0i * xd11 * yd2 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4) * ln2
                + 64.0 * 1.0i * xd13 * yd2 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4) * ln2 + 2.0 * xd4 * yd4 * (-40.0 * ln2 + 8.0 * yd2 * (9.0 + 29.0 * ln2) - 3.0 * yd4 * (108.0 + 211.0 * ln2) + yd6 * (312.0 + 633.0 * ln2))
                + 2.0 * xd20 * yd4 * (-40.0 * ln2 + 8.0 * yd2 * (9.0 + 29.0 * ln2) - 3.0 * yd4 * (108.0 + 211.0 * ln2) + yd6 * (312.0 + 633.0 * ln2))
                - 4.0 * xd6 * yd2 * (32.0 * ln2 - 32.0 * yd2 * (4.0 + 23.0 * ln2) + 4.0 * yd4 * (292.0 + 485.0 * ln2) + yd8 * (1144.0 + 2253.0 * ln2) - yd6 * (2599.0 + 4817.0 * ln2))
                - 4.0 * xd18 * yd2 * (32.0 * ln2 - 32.0 * yd2 * (4.0 + 23.0 * ln2) + 4.0 * yd4 * (292.0 + 485.0 * ln2) + yd8 * (1144.0 + 2253.0 * ln2) - yd6 * (2599.0 + 4817.0 * ln2))
                + 8.0 * xd10 * yd2 * (-16.0 * (304.0 + 207.0 * ln2) + 5.0 * yd8 * (-704.0 + 227.0 * ln2) + 16.0 * yd2 * (1484.0 + 1665.0 * ln2) + yd6 * (16491.0 + 1685.0 * ln2) - 4.0 * yd4 * (7534.0 + 5165.0 * ln2))
                + 8.0 * xd14 * yd2 * (-16.0 * (304.0 + 207.0 * ln2) + 5.0 * yd8 * (-704.0 + 227.0 * ln2) + 16.0 * yd2 * (1484.0 + 1665.0 * ln2) + yd6 * (16491.0 + 1685.0 * ln2) - 4.0 * yd4 * (7534.0 + 5165.0 * ln2))
                - 4.0 * xd12 * (8192.0 + yd8 * (38412.0 - 64785.0 * ln2) + 128.0 * yd2 * (-358.0 + 125.0 * ln2) - 24.0 * yd4 * (-3168.0 + 4955.0 * ln2) + yd10 * (-8168.0 + 11409.0 * ln2) + 8.0 * yd6 * (-9123.0 + 17801.0 * ln2))
                + xd8 * yd2 * (256.0 * (2.0 + 21.0 * ln2) - 64.0 * yd2 * (280.0 + 339.0 * ln2) + 64.0 * yd4 * (1045.0 + 1739.0 * ln2) - 5.0 * yd6 * (11744.0 + 19107.0 * ln2) + yd8 * (15808.0 + 21551.0 * ln2))
                + xd16 * yd2 * (256.0 * (2.0 + 21.0 * ln2) - 64.0 * yd2 * (280.0 + 339.0 * ln2) + 64.0 * yd4 * (1045.0 + 1739.0 * ln2) - 5.0 * yd6 * (11744.0 + 19107.0 * ln2) + yd8 * (15808.0 + 21551.0 * ln2))
                - 4.0 * xd2 * yd6 * (yd4 * (8.0 + 17.0 * ln2) + log(16.0) - yd2 * (3.0 + log(32.0))) - 4.0 * xd22 * yd6 * (yd4 * (8.0 + 17.0 * ln2) + log(16.0) - yd2 * (3.0 + log(32.0))));
            const complex<double> logs18 = 32.0 * ((-yd2) * ln2 + xd4 * (yd2 * (192.0 - 543.0 * ln2) + yd6 * (60.0 - 259.0 * ln2) + 240.0 * ln2 + 2.0 * yd4 * (-96.0 + 281.0 * ln2))
                + xd8 * (yd2 * (192.0 - 543.0 * ln2) + yd6 * (60.0 - 259.0 * ln2) + 240.0 * ln2 + 2.0 * yd4 * (-96.0 + 281.0 * ln2)) - 4.0 * xd6 * (-2.0 * (32.0 + 61.0 * ln2) - 4.0 * yd4 * (18.0 + 77.0 * ln2) + yd6 * (20.0 + 131.0 * ln2) + yd2 * (96.0 + 299.0 * ln2))
                - yd4 * ln4 + yd6 * (4.0 + log(8.0)) + xd12 * ((-yd2) * ln2 - yd4 * ln4 + yd6 * (4.0 + log(8.0))) + xd2 * (-4.0 * ln2 + 54.0 * yd2 * ln2 + yd6 * (-24.0 + ln64) - 8.0 * yd4 * (-6.0 + log(128.0)))
                + xd10 * (-4.0 * ln2 + 54.0 * yd2 * ln2 + yd6 * (-24.0 + ln64) - 8.0 * yd4 * (-6.0 + log(128.0))));
            const complex<double> logs19 = 2.0 * (yd6 * (-1.0 + yd2) + xd24 * yd6 * (-1.0 + yd2) - 32.0 * 1.0i * xd3 * yd4 * power_of<2>(-1.0 + yd2) + 32.0 * 1.0i * xd21 * yd4 * power_of<2>(-1.0 + yd2) - 1920.0 * 1.0i * xd7 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2) + 1920.0 * 1.0i * xd17 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2)
                + 32.0 * 1.0i * xd5 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) - 32.0 * 1.0i * xd19 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) - 128.0 * 1.0i * xd9 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4) + 128.0 * 1.0i * xd15 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4)
                + 64.0 * 1.0i * xd11 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4) - 64.0 * 1.0i * xd13 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4)
                - 2.0 * xd4 * yd2 * (40.0 + yd4 * (633.0 - 1088.0 * ln2) - 160.0 * ln2 + yd6 * (-633.0 + 112.0 * ln2) + 8.0 * yd2 * (-29.0 + 142.0 * ln2))
                - 2.0 * xd20 * yd2 * (40.0 + yd4 * (633.0 - 1088.0 * ln2) - 160.0 * ln2 + yd6 * (-633.0 + 112.0 * ln2) + 8.0 * yd2 * (-29.0 + 142.0 * ln2))
                + 4.0 * xd6 * (-32.0 + yd6 * (4817.0 - 8936.0 * ln2) + yd2 * (736.0 - 3200.0 * ln2) + 3.0 * yd8 * (-751.0 + 836.0 * ln2) + 4.0 * yd4 * (-485.0 + 2407.0 * ln2))
                + 4.0 * xd18 * (-32.0 + yd6 * (4817.0 - 8936.0 * ln2) + yd2 * (736.0 - 3200.0 * ln2) + 3.0 * yd8 * (-751.0 + 836.0 * ln2) + 4.0 * yd4 * (-485.0 + 2407.0 * ln2))
                + 8.0 * xd10 * (yd6 * (1685.0 - 68744.0 * ln2) + 368.0 * (-9.0 + 128.0 * ln2) - 16.0 * yd2 * (-1665.0 + 7804.0 * ln2) + yd8 * (1135.0 + 10012.0 * ln2) + 4.0 * yd4 * (-5165.0 + 34123.0 * ln2))
                + 8.0 * xd14 * (yd6 * (1685.0 - 68744.0 * ln2) + 368.0 * (-9.0 + 128.0 * ln2) - 16.0 * yd2 * (-1665.0 + 7804.0 * ln2) + yd8 * (1135.0 + 10012.0 * ln2) + 4.0 * yd4 * (-5165.0 + 34123.0 * ln2))
                + xd8 * (yd8 * (21551.0 - 71296.0 * ln2) + 256.0 * (21.0 - 80.0 * ln2) + 192.0 * yd2 * (-113.0 + 932.0 * ln2) - 64.0 * yd4 * (-1739.0 + 4986.0 * ln2) + 3.0 * yd6 * (-31845.0 + 77312.0 * ln2))
                + xd16 * (yd8 * (21551.0 - 71296.0 * ln2) + 256.0 * (21.0 - 80.0 * ln2) + 192.0 * yd2 * (-113.0 + 932.0 * ln2) - 64.0 * yd4 * (-1739.0 + 4986.0 * ln2) + 3.0 * yd6 * (-31845.0 + 77312.0 * ln2))
                + 4.0 * xd12 * (128.0 * (-125.0 + 1552.0 * ln2) - 105.0 * yd6 * (-617.0 + 3904.0 * ln2) - 24.0 * yd2 * (-4955.0 + 24812.0 * ln2) + yd8 * (-11409.0 + 80816.0 * ln2) + 8.0 * yd4 * (-17801.0 + 90742.0 * ln2))
                - 4.0 * xd2 * yd4 * (4.0 - 12.0 * ln2 + yd4 * (17.0 + log(16.0)) + yd2 * (-5.0 + ln256)) - 4.0 * xd22 * yd4 * (4.0 - 12.0 * ln2 + yd4 * (17.0 + log(16.0)) + yd2 * (-5.0 + ln256)));
            const complex<double> logs20 = 2.0 * (yd6 * (-1.0 + yd2) + xd24 * yd6 * (-1.0 + yd2) + 32.0 * 1.0i * xd3 * yd4 * power_of<2>(-1.0 + yd2) - 32.0 * 1.0i * xd21 * yd4 * power_of<2>(-1.0 + yd2) + 1920.0 * 1.0i * xd7 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2) - 1920.0 * 1.0i * xd17 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2)
                - 32.0 * 1.0i * xd5 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) + 32.0 * 1.0i * xd19 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) + 128.0 * 1.0i * xd9 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4) - 128.0 * 1.0i * xd15 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4)
                - 64.0 * 1.0i * xd11 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4) + 64.0 * 1.0i * xd13 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4)
                - 2.0 * xd4 * yd2 * (40.0 + yd4 * (633.0 - 1088.0 * ln2) - 160.0 * ln2 + yd6 * (-633.0 + 112.0 * ln2) + 8.0 * yd2 * (-29.0 + 142.0 * ln2))
                - 2.0 * xd20 * yd2 * (40.0 + yd4 * (633.0 - 1088.0 * ln2) - 160.0 * ln2 + yd6 * (-633.0 + 112.0 * ln2) + 8.0 * yd2 * (-29.0 + 142.0 * ln2))
                + 4.0 * xd6 * (-32.0 + yd6 * (4817.0 - 8936.0 * ln2) + yd2 * (736.0 - 3200.0 * ln2) + 3.0 * yd8 * (-751.0 + 836.0 * ln2) + 4.0 * yd4 * (-485.0 + 2407.0 * ln2))
                + 4.0 * xd18 * (-32.0 + yd6 * (4817.0 - 8936.0 * ln2) + yd2 * (736.0 - 3200.0 * ln2) + 3.0 * yd8 * (-751.0 + 836.0 * ln2) + 4.0 * yd4 * (-485.0 + 2407.0 * ln2))
                + 8.0 * xd10 * (yd6 * (1685.0 - 68744.0 * ln2) + 368.0 * (-9.0 + 128.0 * ln2) - 16.0 * yd2 * (-1665.0 + 7804.0 * ln2) + yd8 * (1135.0 + 10012.0 * ln2) + 4.0 * yd4 * (-5165.0 + 34123.0 * ln2))
                + 8.0 * xd14 * (yd6 * (1685.0 - 68744.0 * ln2) + 368.0 * (-9.0 + 128.0 * ln2) - 16.0 * yd2 * (-1665.0 + 7804.0 * ln2) + yd8 * (1135.0 + 10012.0 * ln2) + 4.0 * yd4 * (-5165.0 + 34123.0 * ln2))
                + xd8 * (yd8 * (21551.0 - 71296.0 * ln2) + 256.0 * (21.0 - 80.0 * ln2) + 192.0 * yd2 * (-113.0 + 932.0 * ln2) - 64.0 * yd4 * (-1739.0 + 4986.0 * ln2) + 3.0 * yd6 * (-31845.0 + 77312.0 * ln2))
                + xd16 * (yd8 * (21551.0 - 71296.0 * ln2) + 256.0 * (21.0 - 80.0 * ln2) + 192.0 * yd2 * (-113.0 + 932.0 * ln2) - 64.0 * yd4 * (-1739.0 + 4986.0 * ln2) + 3.0 * yd6 * (-31845.0 + 77312.0 * ln2))
                + 4.0 * xd12 * (128.0 * (-125.0 + 1552.0 * ln2) - 105.0 * yd6 * (-617.0 + 3904.0 * ln2) - 24.0 * yd2 * (-4955.0 + 24812.0 * ln2) + yd8 * (-11409.0 + 80816.0 * ln2) + 8.0 * yd4 * (-17801.0 + 90742.0 * ln2))
                - 4.0 * xd2 * yd4 * (4.0 - 12.0 * ln2 + yd4 * (17.0 + log(16.0)) + yd2 * (-5.0 + ln256)) - 4.0 * xd22 * yd4 * (4.0 - 12.0 * ln2 + yd4 * (17.0 + log(16.0)) + yd2 * (-5.0 + ln256)));
            const complex<double> logs21 = 2.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 * (-1.0 + 9.0 * yd2) + xd8 * yd2 * (-1.0 + 9.0 * yd2) - 2.0 * xd4 * (180.0 - 704.0 * yd3 + yd4 * (233.0 + 512.0 * ln2) - yd2 * (245.0 + 512.0 * ln2))
                - 4.0 * xd2 * (5.0 - 176.0 * yd3 + yd2 * (19.0 - 128.0 * ln2) + 16.0 * yd4 * (1.0 + ln256)) - 4.0 * xd6 * (5.0 - 176.0 * yd3 + yd2 * (19.0 - 128.0 * ln2) + 16.0 * yd4 * (1.0 + ln256)));
            const complex<double> logs22 = 2.0 * power_of<2>(1.0 + xd2) * (-1.0 + yd) * (1.0 + yd) * (yd2 * (-1.0 + 9.0 * yd2) + xd8 * yd2 * (-1.0 + 9.0 * yd2) - 2.0 * xd4 * (180.0 + 704.0 * yd3 + yd4 * (233.0 + 512.0 * ln2) - yd2 * (245.0 + 512.0 * ln2))
                - 4.0 * xd2 * (5.0 + 176.0 * yd3 + yd2 * (19.0 - 128.0 * ln2) + 16.0 * yd4 * (1.0 + ln256)) - 4.0 * xd6 * (5.0 + 176.0 * yd3 + yd2 * (19.0 - 128.0 * ln2) + 16.0 * yd4 * (1.0 + ln256)));
            const complex<double> logs23 = 4.0 * (yd8 * (-11.0 - 16.0 * ln2 + 3.0 * yd2 * (7.0 + 16.0 * ln2)) + xd20 * yd8 * (-11.0 - 16.0 * ln2 + 3.0 * yd2 * (7.0 + 16.0 * ln2)) - 2.0 * xd2 * yd6 * (68.0 + 96.0 * ln2 + yd4 * (203.0 + 432.0 * ln2) - yd2 * (241.0 + 432.0 * ln2))
                - 2.0 * xd18 * yd6 * (68.0 + 96.0 * ln2 + yd4 * (203.0 + 432.0 * ln2) - yd2 * (241.0 + 432.0 * ln2)) + 2.0 * xd8 * yd2 * (256.0 * (81.0 + 88.0 * ln2) + 40.0 * yd4 * (2251.0 + 2056.0 * ln2) - 13.0 * yd6 * (3639.0 + 2576.0 * ln2)
                + yd8 * (9765.0 + 6512.0 * ln2) - 8.0 * yd2 * (9583.0 + 11096.0 * ln2)) + 2.0 * xd12 * yd2 * (256.0 * (81.0 + 88.0 * ln2) + 40.0 * yd4 * (2251.0 + 2056.0 * ln2) - 13.0 * yd6 * (3639.0 + 2576.0 * ln2) + yd8 * (9765.0 + 6512.0 * ln2)
                - 8.0 * yd2 * (9583.0 + 11096.0 * ln2)) + xd4 * yd4 * (-16.0 * (33.0 + 40.0 * ln2) + 16.0 * yd2 * (249.0 + 344.0 * ln2) + yd6 * (2977.0 + 5360.0 * ln2) - yd4 * (7263.0 + 12880.0 * ln2))
                + xd16 * yd4 * (-16.0 * (33.0 + 40.0 * ln2) + 16.0 * yd2 * (249.0 + 344.0 * ln2) + yd6 * (2977.0 + 5360.0 * ln2) - yd4 * (7263.0 + 12880.0 * ln2))
                - 4.0 * xd10 * (64.0 * yd2 * (663.0 + 536.0 * ln2) - 16.0 * yd4 * (4123.0 + 2904.0 * ln2) + yd10 * (5825.0 + 8080.0 * ln2) - 3.0 * yd8 * (9649.0 + 11824.0 * ln2) + yd6 * (59068.0 + 54688.0 * ln2) - 2048.0 * (5.0 + log(16.0)))
                - 8.0 * xd6 * yd2 * (32.0 * yd4 * (195.0 + 316.0 * ln2) - 4.0 * yd2 * (541.0 + 744.0 * ln2) + yd8 * (1309.0 + 1712.0 * ln2) - yd6 * (5055.0 + 7664.0 * ln2) + 16.0 * (5.0 + ln256))
                - 8.0 * xd14 * yd2 * (32.0 * yd4 * (195.0 + 316.0 * ln2) - 4.0 * yd2 * (541.0 + 744.0 * ln2) + yd8 * (1309.0 + 1712.0 * ln2) - yd6 * (5055.0 + 7664.0 * ln2) + 16.0 * (5.0 + ln256)));
            const complex<double> logs24 = 2.0 * (yd6 * (1.0 - 6.0 * yd + 8.0 * yd2 - 2.0 * yd3 + 15.0 * yd4) + xd16 * yd6 * (1.0 - 6.0 * yd + 8.0 * yd2 - 2.0 * yd3 + 15.0 * yd4) - 2.0 * xd8 * (512.0 + 1408.0 * yd - 6224.0 * yd2 + yd4 * (7328.0 - 4288.0 * ln2) + 8.0 * yd8 * (429.0 + 152.0 * ln2)
                + 32.0 * yd3 * (31.0 + 216.0 * ln2) - 16.0 * yd5 * (199.0 + 816.0 * ln2) - yd6 * (3059.0 + 3904.0 * ln2) + yd10 * (-269.0 + 6976.0 * ln2) - 2.0 * yd9 * (437.0 + 9344.0 * ln2) + yd7 * (-622.0 + 24832.0 * ln2))
                + 8.0 * xd2 * yd4 * (-7.0 * yd - 4.0 * ln2 + 28.0 * yd6 * ln2 - 2.0 * yd2 * (3.0 + 14.0 * ln2) - yd5 * (17.0 + 64.0 * ln2) + yd4 * (-2.0 + log(16.0)) + 16.0 * yd3 * (3.0 + log(16.0)))
                + 8.0 * xd14 * yd4 * (-7.0 * yd - 4.0 * ln2 + 28.0 * yd6 * ln2 - 2.0 * yd2 * (3.0 + 14.0 * ln2) - yd5 * (17.0 + 64.0 * ln2) + yd4 * (-2.0 + log(16.0)) + 16.0 * yd3 * (3.0 + log(16.0)))
                - 8.0 * xd6 * yd * (48.0 + 184.0 * yd + yd8 * (207.0 - 2752.0 * ln2) + yd4 * (921.0 - 1408.0 * ln2) + 48.0 * yd6 * (-31.0 + 68.0 * ln2) - 4.0 * yd3 * (356.0 + 113.0 * ln2) + 4.0 * yd9 * (-8.0 + 215.0 * ln2) + 2.0 * yd7 * (55.0 + 482.0 * ln2)
                - 2.0 * yd5 * (-737.0 + 686.0 * ln2) + 112.0 * yd2 * (-1.0 + ln256)) - 8.0 * xd10 * yd * (48.0 + 184.0 * yd + yd8 * (207.0 - 2752.0 * ln2) + yd4 * (921.0 - 1408.0 * ln2) + 48.0 * yd6 * (-31.0 + 68.0 * ln2) - 4.0 * yd3 * (356.0 + 113.0 * ln2)
                + 4.0 * yd9 * (-8.0 + 215.0 * ln2) + 2.0 * yd7 * (55.0 + 482.0 * ln2) - 2.0 * yd5 * (-737.0 + 686.0 * ln2) + 112.0 * yd2 * (-1.0 + ln256))
                + 4.0 * xd4 * yd2 * (-4.0 + yd3 * (420.0 - 448.0 * ln2) + 5.0 * yd8 * (-27.0 + 16.0 * ln2) - 8.0 * yd2 * (23.0 + 22.0 * ln2) - 6.0 * yd5 * (111.0 + 32.0 * ln2) - 8.0 * yd6 * (65.0 + 202.0 * ln2) + yd7 * (614.0 + 704.0 * ln2) + yd4 * (627.0 + 1712.0 * ln2)
                - 8.0 * yd * (5.0 + ln256)) + 4.0 * xd12 * yd2 * (-4.0 + yd3 * (420.0 - 448.0 * ln2) + 5.0 * yd8 * (-27.0 + 16.0 * ln2) - 8.0 * yd2 * (23.0 + 22.0 * ln2) - 6.0 * yd5 * (111.0 + 32.0 * ln2) - 8.0 * yd6 * (65.0 + 202.0 * ln2) + yd7 * (614.0 + 704.0 * ln2)
                + yd4 * (627.0 + 1712.0 * ln2) - 8.0 * yd * (5.0 + ln256)));
            const complex<double> logs25 = 2.0 * (yd6 * (1.0 + 6.0 * yd + 8.0 * yd2 + 2.0 * yd3 + 15.0 * yd4) + xd16 * yd6 * (1.0 + 6.0 * yd + 8.0 * yd2 + 2.0 * yd3 + 15.0 * yd4) - 2.0 * xd8 * (512.0 - 1408.0 * yd - 6224.0 * yd2 + yd7 * (622.0 - 24832.0 * ln2) + yd4 * (7328.0 - 4288.0 * ln2)
                + 8.0 * yd8 * (429.0 + 152.0 * ln2) - 32.0 * yd3 * (31.0 + 216.0 * ln2) + 16.0 * yd5 * (199.0 + 816.0 * ln2) - yd6 * (3059.0 + 3904.0 * ln2) + yd10 * (-269.0 + 6976.0 * ln2) + 2.0 * yd9 * (437.0 + 9344.0 * ln2))
                + 8.0 * xd2 * yd4 * (7.0 * yd - 4.0 * ln2 + 28.0 * yd6 * ln2 - 2.0 * yd2 * (3.0 + 14.0 * ln2) + yd5 * (17.0 + 64.0 * ln2) + yd4 * (-2.0 + log(16.0)) - 16.0 * yd3 * (3.0 + log(16.0)))
                + 8.0 * xd14 * yd4 * (7.0 * yd - 4.0 * ln2 + 28.0 * yd6 * ln2 - 2.0 * yd2 * (3.0 + 14.0 * ln2) + yd5 * (17.0 + 64.0 * ln2) + yd4 * (-2.0 + log(16.0)) - 16.0 * yd3 * (3.0 + log(16.0)))
                - 8.0 * xd6 * yd * (-48.0 + 184.0 * yd - 48.0 * yd6 * (-31.0 + 68.0 * ln2) - 4.0 * yd3 * (356.0 + 113.0 * ln2) + 4.0 * yd9 * (-8.0 + 215.0 * ln2) + 2.0 * yd7 * (55.0 + 482.0 * ln2) - 2.0 * yd5 * (-737.0 + 686.0 * ln2) + yd4 * (-921.0 + 1408.0 * ln2)
                + yd8 * (-207.0 + 2752.0 * ln2) - 112.0 * yd2 * (-1.0 + ln256)) - 8.0 * xd10 * yd * (-48.0 + 184.0 * yd - 48.0 * yd6 * (-31.0 + 68.0 * ln2) - 4.0 * yd3 * (356.0 + 113.0 * ln2) + 4.0 * yd9 * (-8.0 + 215.0 * ln2) + 2.0 * yd7 * (55.0 + 482.0 * ln2)
                - 2.0 * yd5 * (-737.0 + 686.0 * ln2) + yd4 * (-921.0 + 1408.0 * ln2) + yd8 * (-207.0 + 2752.0 * ln2) - 112.0 * yd2 * (-1.0 + ln256))
                + 4.0 * xd4 * yd2 * (-4.0 + 5.0 * yd8 * (-27.0 + 16.0 * ln2) + 28.0 * yd3 * (-15.0 + 16.0 * ln2) - 8.0 * yd2 * (23.0 + 22.0 * ln2) + 6.0 * yd5 * (111.0 + 32.0 * ln2) - 8.0 * yd6 * (65.0 + 202.0 * ln2) - 2.0 * yd7 * (307.0 + 352.0 * ln2)
                + yd4 * (627.0 + 1712.0 * ln2) + 8.0 * yd * (5.0 + ln256)) + 4.0 * xd12 * yd2 * (-4.0 + 5.0 * yd8 * (-27.0 + 16.0 * ln2) + 28.0 * yd3 * (-15.0 + 16.0 * ln2) - 8.0 * yd2 * (23.0 + 22.0 * ln2) + 6.0 * yd5 * (111.0 + 32.0 * ln2) - 8.0 * yd6 * (65.0 + 202.0 * ln2)
                - 2.0 * yd7 * (307.0 + 352.0 * ln2) + yd4 * (627.0 + 1712.0 * ln2) + 8.0 * yd * (5.0 + ln256)));
            const complex<double> logs26 = 2.0 * (yd6 * (-1.0 + yd2) + xd24 * yd6 * (-1.0 + yd2) - 32.0 * 1.0i * xd3 * yd4 * power_of<2>(-1.0 + yd2) + 32.0 * 1.0i * xd21 * yd4 * power_of<2>(-1.0 + yd2) - 1920.0 * 1.0i * xd7 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2) + 1920.0 * 1.0i * xd17 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2)
                + 32.0 * 1.0i * xd5 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) - 32.0 * 1.0i * xd19 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) - 128.0 * 1.0i * xd9 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4) + 128.0 * 1.0i * xd15 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4)
                + 64.0 * 1.0i * xd11 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4) - 64.0 * 1.0i * xd13 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4)
                + 8.0 * xd10 * (yd6 * (1685.0 - 28776.0 * ln2) + yd2 * (26640.0 - 28352.0 * ln2) + 16.0 * (-207.0 + 896.0 * ln2) + yd8 * (1135.0 + 5516.0 * ln2) + 4.0 * yd4 * (-5165.0 + 9319.0 * ln2))
                + 8.0 * xd14 * (yd6 * (1685.0 - 28776.0 * ln2) + yd2 * (26640.0 - 28352.0 * ln2) + 16.0 * (-207.0 + 896.0 * ln2) + yd8 * (1135.0 + 5516.0 * ln2) + 4.0 * yd4 * (-5165.0 + 9319.0 * ln2))
                + 4.0 * xd12 * (3200.0 * (-5.0 + 16.0 * ln2) + 3.0 * yd8 * (-3803.0 + 9552.0 * ln2) - 8.0 * yd2 * (-14865.0 + 14884.0 * ln2) - 7.0 * yd6 * (-9255.0 + 17728.0 * ln2) + 8.0 * yd4 * (-17801.0 + 20414.0 * ln2))
                + xd8 * (yd8 * (21551.0 - 17536.0 * ln2) + 768.0 * (7.0 + 16.0 * ln2) + 192.0 * yd2 * (-113.0 + 68.0 * ln2) - 64.0 * yd4 * (-1739.0 + 498.0 * ln2) + yd6 * (-95535.0 + 24064.0 * ln2))
                + xd16 * (yd8 * (21551.0 - 17536.0 * ln2) + 768.0 * (7.0 + 16.0 * ln2) + 192.0 * yd2 * (-113.0 + 68.0 * ln2) - 64.0 * yd4 * (-1739.0 + 498.0 * ln2) + yd6 * (-95535.0 + 24064.0 * ln2))
                - 4.0 * xd2 * yd4 * (4.0 + yd4 * (17.0 + 20.0 * ln2) - yd2 * (5.0 + 24.0 * ln2) + log(16.0)) - 4.0 * xd22 * yd4 * (4.0 + yd4 * (17.0 + 20.0 * ln2) - yd2 * (5.0 + 24.0 * ln2) + log(16.0))
                + 2.0 * xd4 * yd2 * (-3.0 * yd4 * (211.0 + 192.0 * ln2) + yd6 * (633.0 + 592.0 * ln2) + 8.0 * yd2 * (29.0 + log(1024.0)) - 8.0 * (5.0 + log(4096.0)))
                + 2.0 * xd20 * yd2 * (-3.0 * yd4 * (211.0 + 192.0 * ln2) + yd6 * (633.0 + 592.0 * ln2) + 8.0 * yd2 * (29.0 + log(1024.0)) - 8.0 * (5.0 + log(4096.0)))
                - 4.0 * xd6 * (32.0 + 20.0 * yd4 * (97.0 + 41.0 * ln2) + yd8 * (2253.0 + 772.0 * ln2) - yd6 * (4817.0 + 1208.0 * ln2) - 32.0 * yd2 * (23.0 + log(4096.0)))
                - 4.0 * xd18 * (32.0 + 20.0 * yd4 * (97.0 + 41.0 * ln2) + yd8 * (2253.0 + 772.0 * ln2) - yd6 * (4817.0 + 1208.0 * ln2) - 32.0 * yd2 * (23.0 + log(4096.0))));
            const complex<double> logs27 = 2.0 * (yd6 * (-1.0 + yd2) + xd24 * yd6 * (-1.0 + yd2) + 32.0 * 1.0i * xd3 * yd4 * power_of<2>(-1.0 + yd2) - 32.0 * 1.0i * xd21 * yd4 * power_of<2>(-1.0 + yd2) + 1920.0 * 1.0i * xd7 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2) - 1920.0 * 1.0i * xd17 * yd2 * (-2.0 + yd2) * power_of<2>(-1.0 + yd2)
                - 32.0 * 1.0i * xd5 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) + 32.0 * 1.0i * xd19 * yd2 * power_of<2>(-1.0 + yd2) * (-8.0 + 17.0 * yd2) + 128.0 * 1.0i * xd9 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4) - 128.0 * 1.0i * xd15 * power_of<2>(-1.0 + yd2) * (-128.0 + 26.0 * yd2 + 73.0 * yd4)
                - 64.0 * 1.0i * xd11 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4) + 64.0 * 1.0i * xd13 * power_of<2>(-1.0 + yd2) * (1280.0 - 2164.0 * yd2 + 917.0 * yd4)
                + 8.0 * xd10 * (yd6 * (1685.0 - 28776.0 * ln2) + yd2 * (26640.0 - 28352.0 * ln2) + 16.0 * (-207.0 + 896.0 * ln2) + yd8 * (1135.0 + 5516.0 * ln2) + 4.0 * yd4 * (-5165.0 + 9319.0 * ln2))
                + 8.0 * xd14 * (yd6 * (1685.0 - 28776.0 * ln2) + yd2 * (26640.0 - 28352.0 * ln2) + 16.0 * (-207.0 + 896.0 * ln2) + yd8 * (1135.0 + 5516.0 * ln2) + 4.0 * yd4 * (-5165.0 + 9319.0 * ln2))
                + 4.0 * xd12 * (3200.0 * (-5.0 + 16.0 * ln2) + 3.0 * yd8 * (-3803.0 + 9552.0 * ln2) - 8.0 * yd2 * (-14865.0 + 14884.0 * ln2) - 7.0 * yd6 * (-9255.0 + 17728.0 * ln2) + 8.0 * yd4 * (-17801.0 + 20414.0 * ln2))
                + xd8 * (yd8 * (21551.0 - 17536.0 * ln2) + 768.0 * (7.0 + 16.0 * ln2) + 192.0 * yd2 * (-113.0 + 68.0 * ln2) - 64.0 * yd4 * (-1739.0 + 498.0 * ln2) + yd6 * (-95535.0 + 24064.0 * ln2))
                + xd16 * (yd8 * (21551.0 - 17536.0 * ln2) + 768.0 * (7.0 + 16.0 * ln2) + 192.0 * yd2 * (-113.0 + 68.0 * ln2) - 64.0 * yd4 * (-1739.0 + 498.0 * ln2) + yd6 * (-95535.0 + 24064.0 * ln2))
                - 4.0 * xd2 * yd4 * (4.0 + yd4 * (17.0 + 20.0 * ln2) - yd2 * (5.0 + 24.0 * ln2) + log(16.0)) - 4.0 * xd22 * yd4 * (4.0 + yd4 * (17.0 + 20.0 * ln2) - yd2 * (5.0 + 24.0 * ln2) + log(16.0))
                + 2.0 * xd4 * yd2 * (-3.0 * yd4 * (211.0 + 192.0 * ln2) + yd6 * (633.0 + 592.0 * ln2) + 8.0 * yd2 * (29.0 + log(1024.0)) - 8.0 * (5.0 + log(4096.0)))
                + 2.0 * xd20 * yd2 * (-3.0 * yd4 * (211.0 + 192.0 * ln2) + yd6 * (633.0 + 592.0 * ln2) + 8.0 * yd2 * (29.0 + log(1024.0)) - 8.0 * (5.0 + log(4096.0)))
                - 4.0 * xd6 * (32.0 + 20.0 * yd4 * (97.0 + 41.0 * ln2) + yd8 * (2253.0 + 772.0 * ln2) - yd6 * (4817.0 + 1208.0 * ln2) - 32.0 * yd2 * (23.0 + log(4096.0)))
                - 4.0 * xd18 * (32.0 + 20.0 * yd4 * (97.0 + 41.0 * ln2) + yd8 * (2253.0 + 772.0 * ln2) - yd6 * (4817.0 + 1208.0 * ln2) - 32.0 * yd2 * (23.0 + log(4096.0))));

            const complex<double> factor1 = 16.0 / (9.0 * yd3);

            const complex<double> f29dPart1 = (64.0 * power_of<2>(log(xd))) / 9.0 + (num4 * dilog((-1.0i) * xd) + num3 * dilog(1.0i * xd) + num3 * log(1.0 - 1.0i * xd) * log(xd)
                    + num4 * log(1.0 + 1.0i * xd) * log(xd)) / denom3 + (num1 * log(xd) * ln1myd - num1 * log(xd) * ln1pyd) / denom1
                + ((-num11 - num12) * dilog(-1.0i / (-1.0i + wx3)) + (-num11 - num12) * dilog(1.0i / (1.0i + wx3))
                + (-num11 - num12) * dilog(-1.0i / (-1.0i + wx4)) + (-num11 - num12) * dilog(1.0i / (1.0i + wx4))
                + num11 * dilog((1.0i - xd) / (1.0i + wx3)) + num11 * dilog((1.0i - xd) / (1.0i + wx4)) - num11 * dilog(-(xd / wx3))
                - num12 * dilog(xd / wx3) - num11 * dilog(-(xd / wx4)) - num12 * dilog(xd / wx4)
                + num12 * dilog((-1.0i + xd) / (-1.0i + wx3)) + num12 * dilog((-1.0i + xd) / (-1.0i + wx4)) + num11 * dilog((1.0i + xd) / (1.0i - wx3))
                + num12 * dilog((1.0i + xd) / (1.0i + wx3)) + num11 * dilog((1.0i + xd) / (1.0i - wx4)) + num12 * dilog((1.0i + xd) / (1.0i + wx4))
                + num12 * log((wx3 - xd) / (1.0i + wx3)) * log(1.0 - 1.0i * xd) + num12 * log((wx3 - xd) / (-1.0i + wx3)) * log(1.0 + 1.0i * xd)
                + (-2.0 * num11 * ln2 - num11 * log(xd)) * log((wx3 + xd) / wx3) - 2.0 * num11 * ln2 * log((wx4 + xd) / wx4)
                + log(1.0 - 1.0i * xd) * (num12 * log((wx4 - xd) / (1.0i + wx4)) + num11 * log((wx3 + xd) / (-1.0i + wx3))
                    + num11 * log((wx4 + xd) / (-1.0i + wx4))) + log(1.0 + 1.0i * xd) * (num12 * log((wx4 - xd) / (-1.0i + wx4))
                    + num11 * log((wx3 + xd) / (1.0i + wx3)) + num11 * log((wx4 + xd) / (1.0i + wx4))) - 2.0 * num12 * ln2 * log(1.0 - xd / wx3)
                - 2.0 * num12 * ln2 * log(1.0 - xd / wx4) + log(xd) * (-(num11 * log((wx4 + xd) / wx4)) - num12 * log(1.0 - xd / wx3)
                    - num12 * log(1.0 - xd / wx4)) + (2.0 * 1.0i) * num12 * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3)
                + (2.0 * 1.0i) * num11 * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                + (2.0 * 1.0i) * num12 * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4)
                + (2.0 * 1.0i) * num11 * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                + (2.0 * 1.0i) * num11 * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3)
                + (2.0 * 1.0i) * num12 * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                + (2.0 * 1.0i) * num11 * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4)
                + (2.0 * 1.0i) * num12 * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4)) / denom4;

            const complex<double> f29dPart2 = ((-1.0 / 12.0) * (num10 * pisqu * ln2) - (num6 * pisqu * ln2) / 12.0 + num8 * trilog(1.0 / 2.0 - (1.0i / 2.0) * xd) + (-num8 - num9) * trilog(1.0 - 1.0i * xd) + num8 * trilog((1.0 + 1.0i * xd) / 2.0) + (-num8 - num9) * trilog(1.0 + 1.0i * xd) - 2.0 * num8 * trilog((-1.0i) * xd) - 2.0 * num8 * trilog(1.0i * xd) + (num8 * trilog(-xd2)) / 2.0 - num8 * trilog((2.0 * xd) / (-1.0i + xd))
                - num8 * trilog((2.0 * xd) / (1.0i + xd)) - num10 * trilog((1.0 - yd) / 2.0) - 7.0 * num7 * trilog(-yd) + 7.0 * num7 * trilog(yd) - num6 * trilog((1.0 + yd) / 2.0) + (num10 * power_of<3>(ln2)) / 6.0 + (num6 * power_of<3>(ln2)) / 6.0 - (num8 * power_of<3>(ln2)) / 3.0 + (num8 * pisqu * ln4) / 12.0 + num9 * dilog(1.0 - 1.0i * xd) * log(1.0 - 1.0i * xd) + num9 * dilog(1.0i * xd) * log(1.0 - 1.0i * xd)
                + (num8 * power_of<3>(log(1.0 - 1.0i * xd))) / 6.0 + num9 * dilog(1.0 + 1.0i * xd) * log(1.0 + 1.0i * xd) + num9 * dilog((-1.0i) * xd) * log(1.0 + 1.0i * xd) + (num8 * power_of<3>(log(1.0 + 1.0i * xd))) / 6.0 + ((-1.0 / 12.0) * (num8 * pisqu) + (num8 * ln2squ) / 2.0) * log((-1.0i) * xd) + ((-1.0 / 12.0) * (num8 * pisqu) + (num8 * ln2squ) / 2.0) * log(1.0i * xd) + 6.0 * num7 * dilog(-yd) * log(xd)
                - 6.0 * num7 * dilog(yd) * log(xd) + dilog(1.0 / 2.0 - (1.0i / 2.0) * xd) * (num8 * log((-1.0i) * xd) + num8 * log(1.0i * xd) - num8 * log(xd)) + dilog((1.0 + 1.0i * xd) / 2.0) * (num8 * log((-1.0i) * xd) + num8 * log(1.0i * xd) - num8 * log(xd)) + li2half * (-(num8 * log((-1.0i) * xd)) - num8 * log(1.0i * xd) + (-num2 - num5 + 2.0 * num8) * log(xd))
                + power_of<2>(log(1.0 + 1.0i * xd)) * ((num9 * log((-1.0i) * xd)) / 2.0 - (num8 * log((-2.0 * 1.0i) * xd)) / 2.0 + (num9 * log(xd)) / 2.0) + power_of<2>(log(1.0 - 1.0i * xd)) * ((num9 * log(1.0i * xd)) / 2.0 - (num8 * log((2.0 * 1.0i) * xd)) / 2.0 + (num9 * log(xd)) / 2.0)
                + log(1.0 + 1.0i * xd) * ((num8 * pisqu) / 12.0 + (num8 * ln2squ) / 2.0 - num8 * ln2 * log((-1.0i) * xd) - num8 * ln2 * log(1.0i * xd) + num8 * ln2 * log(xd)) + log(1.0 - 1.0i * xd) * ((num8 * pisqu) / 12.0 + (num8 * ln2squ) / 2.0 - num8 * ln2 * log((-1.0i) * xd) - num8 * ln2 * log(1.0i * xd)
                + log(1.0 + 1.0i * xd) * (num8 * log((-1.0i) * xd) + num8 * log(1.0i * xd)) + num8 * ln2 * log(xd)) + ((-1.0 / 12.0) * (num10 * ln64) - (num5 * log(1.0 - 1.0i * xd)) / 2.0 - (num5 * log(1.0 + 1.0i * xd)) / 2.0 + (num5 * log(xd)) / 2.0) * power_of<2>(ln1myd) - (num10 * power_of<3>(ln1myd)) / 6.0
                + ((-1.0 / 12.0) * (num6 * pisqu) + (num6 * ln2squ) / 2.0 - num2 * ln2 * log(xd)) * ln1pyd + ((-1.0 / 12.0) * (num6 * ln64) - (num2 * log(1.0 - 1.0i * xd)) / 2.0 - (num2 * log(1.0 + 1.0i * xd)) / 2.0 + (num2 * log(xd)) / 2.0) * power_of<2>(ln1pyd) - (num6 * power_of<3>(ln1pyd)) / 6.0
                + dilog((1.0 + yd) / 2.0) * (num2 * log(xd) + num10 * ln1myd + num6 * ln1pyd) + dilog((1.0 - yd) / 2.0) * (num5 * log(xd) + num10 * ln1myd + num6 * ln1pyd)
                + ln1myd * ((-1.0 / 12.0) * (num10 * pisqu) + (num10 * ln2squ) / 2.0 - num5 * ln2 * log(xd) + (num10 * log((1.0 - yd) / 2.0) + num6 * log((1.0 + yd) / 2.0)) * ln1pyd) + (7.0 * num10 * zeta3) / 8.0 + (7.0 * num6 * zeta3) / 8.0 + (num8 * zeta3) / 4.0 + 2.0 * num9 * zeta3) / denom2;

            const complex<double> f29dPart3 = (dilog(1.0 / 2.0 - (1.0i / 2.0) * xd) * (-2.0 * num14 * ln1myd - 2.0 * num13 * ln1pyd) + dilog((1.0 + 1.0i * xd) / 2.0) * (-2.0 * num14 * ln1myd - 2.0 * num13 * ln1pyd)
                + dilog((-1.0i) * xd) * (-2.0 * num14 * ln1myd - 2.0 * num13 * ln1pyd) + dilog(1.0i * xd) * (-2.0 * num14 * ln1myd - 2.0 * num13 * ln1pyd)
                + dilog((1.0i - xd) / (1.0i + wx3)) * (-(num14 * ln1myd) - num13 * ln1pyd) + dilog((1.0i - xd) / (1.0i + wx4)) * (-(num14 * ln1myd) - num13 * ln1pyd)
                + dilog((-1.0i + xd) / (-1.0i + wx3)) * (-(num14 * ln1myd) - num13 * ln1pyd) + dilog((-1.0i + xd) / (-1.0i + wx4)) * (-(num14 * ln1myd) - num13 * ln1pyd)
                + dilog((1.0i + xd) / (1.0i - wx3)) * (-(num14 * ln1myd) - num13 * ln1pyd) + dilog((1.0i + xd) / (1.0i + wx3)) * (-(num14 * ln1myd) - num13 * ln1pyd)
                + dilog((1.0i + xd) / (1.0i - wx4)) * (-(num14 * ln1myd) - num13 * ln1pyd) + dilog((1.0i + xd) / (1.0i + wx4)) * (-(num14 * ln1myd) - num13 * ln1pyd)
                + dilog(-(xd / wx3)) * (num14 * ln1myd + num13 * ln1pyd) + dilog(xd / wx3) * (num14 * ln1myd + num13 * ln1pyd) + dilog(-(xd / wx4)) * (num14 * ln1myd + num13 * ln1pyd)
                + dilog(xd / wx4) * (num14 * ln1myd + num13 * ln1pyd) + dilog(-1.0i / (-1.0i + wx3)) * (2.0 * num14 * ln1myd + 2.0 * num13 * ln1pyd) + dilog(1.0i / (1.0i + wx3)) * (2.0 * num14 * ln1myd + 2.0 * num13 * ln1pyd)
                + dilog(-1.0i / (-1.0i + wx4)) * (2.0 * num14 * ln1myd + 2.0 * num13 * ln1pyd) + dilog(1.0i / (1.0i + wx4)) * (2.0 * num14 * ln1myd + 2.0 * num13 * ln1pyd) + li2half * (4.0 * num14 * ln1myd + 4.0 * num13 * ln1pyd)
                + ln1pyd * (num13 * power_of<2>(log(1.0 - 1.0i * xd)) + num13 * power_of<2>(log(1.0 + 1.0i * xd)) + 2.0 * num13 * ln2 * log((wx3 + xd) / wx3) + 2.0 * num13 * ln2 * log((wx4 + xd) / wx4)
                + log(1.0 - 1.0i * xd) * (2.0 * num13 * ln2 - num13 * log((wx3 - xd) / (1.0i + wx3)) - num13 * log((wx4 - xd) / (1.0i + wx4)) - 2.0 * num13 * log(xd) - num13 * log((wx3 + xd) / (-1.0i + wx3)) - num13 * log((wx4 + xd) / (-1.0i + wx4)))
                + log(1.0 + 1.0i * xd) * (2.0 * num13 * ln2 - num13 * log((wx3 - xd) / (-1.0i + wx3)) - num13 * log((wx4 - xd) / (-1.0i + wx4)) - 2.0 * num13 * log(xd) - num13 * log((wx3 + xd) / (1.0i + wx3)) - num13 * log((wx4 + xd) / (1.0i + wx4)))
                + 2.0 * num13 * ln2 * log(1.0 - xd / wx3) + 2.0 * num13 * ln2 * log(1.0 - xd / wx4) + log(xd) * (num13 * log((wx3 + xd) / wx3) + num13 * log((wx4 + xd) / wx4) + num13 * log(1.0 - xd / wx3) + num13 * log(1.0 - xd / wx4))
                - (2.0 * 1.0i) * num13 * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) - (2.0 * 1.0i) * num13 * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                - (2.0 * 1.0i) * num13 * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) - (2.0 * 1.0i) * num13 * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                - (2.0 * 1.0i) * num13 * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) - (2.0 * 1.0i) * num13 * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                - (2.0 * 1.0i) * num13 * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) - (2.0 * 1.0i) * num13 * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))
                + ln1myd * (num14 * power_of<2>(log(1.0 - 1.0i * xd)) + num14 * power_of<2>(log(1.0 + 1.0i * xd)) + 2.0 * num14 * ln2 * log((wx3 + xd) / wx3) + 2.0 * num14 * ln2 * log((wx4 + xd) / wx4)
                + log(1.0 - 1.0i * xd) * (2.0 * num14 * ln2 - num14 * log((wx3 - xd) / (1.0i + wx3)) - num14 * log((wx4 - xd) / (1.0i + wx4)) - 2.0 * num14 * log(xd) - num14 * log((wx3 + xd) / (-1.0i + wx3)) - num14 * log((wx4 + xd) / (-1.0i + wx4)))
                + log(1.0 + 1.0i * xd) * (2.0 * num14 * ln2 - num14 * log((wx3 - xd) / (-1.0i + wx3)) - num14 * log((wx4 - xd) / (-1.0i + wx4)) - 2.0 * num14 * log(xd) - num14 * log((wx3 + xd) / (1.0i + wx3)) - num14 * log((wx4 + xd) / (1.0i + wx4)))
                + 2.0 * num14 * ln2 * log(1.0 - xd / wx3) + 2.0 * num14 * ln2 * log(1.0 - xd / wx4) + log(xd) * (num14 * log((wx3 + xd) / wx3) + num14 * log((wx4 + xd) / wx4) + num14 * log(1.0 - xd / wx3) + num14 * log(1.0 - xd / wx4))
                - (2.0 * 1.0i) * num14 * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) - (2.0 * 1.0i) * num14 * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                - (2.0 * 1.0i) * num14 * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) - (2.0 * 1.0i) * num14 * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                - (2.0 * 1.0i) * num14 * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) - (2.0 * 1.0i) * num14 * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                - (2.0 * 1.0i) * num14 * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) - (2.0 * 1.0i) * num14 * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))) / denom2;

            const complex<double> f29dPart41 = (num15 / denom2) * (16.0 * li3half - 6.0 * trilog(-1.0i / wx3) - 6.0 * trilog(1.0i / wx3) - 10.0 * trilog(-1.0i / (-1.0i + wx3)) - 2.0 * trilog(wx3 / (-1.0i + wx3)) + 4.0 * trilog((-1.0i + wx3) / (2.0 * wx3)) - 4.0 * trilog((-1.0i + wx3) / wx3)
                - 10.0 * trilog(1.0i / (1.0i + wx3)) + 4.0 * trilog((1.0i - wx3) / (1.0i + wx3)) - 2.0 * trilog(wx3 / (1.0i + wx3)) + 4.0 * trilog((1.0i + wx3) / (1.0i - wx3)) + 4.0 * trilog((1.0i + wx3) / (2.0 * wx3)) - 4.0 * trilog((1.0i + wx3) / wx3) - 6.0 * trilog(-1.0i / wx4)
                - 6.0 * trilog(1.0i / wx4) - 10.0 * trilog(-1.0i / (-1.0i + wx4)) - 2.0 * trilog(wx4 / (-1.0i + wx4)) + 4.0 * trilog((-1.0i + wx4) / (2.0 * wx4)) - 4.0 * trilog((-1.0i + wx4) / wx4) - 10.0 * trilog(1.0i / (1.0i + wx4)) + 4.0 * trilog((1.0i - wx4) / (1.0i + wx4))
                - 2.0 * trilog(wx4 / (1.0i + wx4)) + 4.0 * trilog((1.0i + wx4) / (1.0i - wx4)) + 4.0 * trilog((1.0i + wx4) / (2.0 * wx4)) - 4.0 * trilog((1.0i + wx4) / wx4) + 5.0 * trilog((1.0i - xd) / (1.0i + wx3)) + 5.0 * trilog((1.0i - xd) / (1.0i + wx4))
                + 2.0 * trilog((-1.0i + wx3) / (wx3 - xd)) + 2.0 * trilog((1.0i + wx3) / (wx3 - xd)) + trilog((wx3 - xd) / (-1.0i + wx3)) + trilog((wx3 - xd) / (1.0i + wx3)) + 2.0 * trilog((-1.0i + wx4) / (wx4 - xd)) + 2.0 * trilog((1.0i + wx4) / (wx4 - xd))
                + trilog((wx4 - xd) / (-1.0i + wx4)) + trilog((wx4 - xd) / (1.0i + wx4)) - 8.0 * trilog(1.0 / 2.0 - (1.0i / 2.0) * xd) - 2.0 * trilog(((-1.0i + wx3) * (1.0 - 1.0i * xd)) / (2.0 * (wx3 - xd))) - 2.0 * trilog(((-1.0i + wx4) * (1.0 - 1.0i * xd)) / (2.0 * (wx4 - xd)))
                - 8.0 * trilog((1.0 + 1.0i * xd) / 2.0) - 2.0 * trilog(((1.0i + wx3) * (1.0 + 1.0i * xd)) / (2.0 * (wx3 - xd))) - 2.0 * trilog(((1.0i + wx4) * (1.0 + 1.0i * xd)) / (2.0 * (wx4 - xd))) - trilog((1.0i * (wx3 - xd)) / ((-1.0i + wx3) * xd))
                - trilog(((-1.0i) * (wx3 - xd)) / ((1.0i + wx3) * xd)) - trilog((1.0i * (wx4 - xd)) / ((-1.0i + wx4) * xd)) - trilog(((-1.0i) * (wx4 - xd)) / ((1.0i + wx4) * xd)) - 2.0 * trilog(-(xd / wx3)) - 2.0 * trilog(xd / wx3) - 2.0 * trilog(-(xd / wx4))
                - 2.0 * trilog(xd / wx4) + 4.0 * trilog(-1.0i / (-1.0i + xd)) - trilog(((-1.0i) * (wx3 - xd)) / (wx3 * (-1.0i + xd))) - trilog(((-1.0i) * (wx4 - xd)) / (wx4 * (-1.0i + xd))) + 5.0 * trilog((-1.0i + xd) / (-1.0i + wx3)) + 5.0 * trilog((-1.0i + xd) / (-1.0i + wx4))
                - 4.0 * trilog((-1.0i + xd) * xdinv) + trilog((wx3 * (-1.0i + xd)) / ((-1.0i + wx3) * xd)) + trilog((wx3 * (-1.0i + xd)) / ((1.0i + wx3) * xd)) + trilog((wx4 * (-1.0i + xd)) / ((-1.0i + wx4) * xd)) + trilog((wx4 * (-1.0i + xd)) / ((1.0i + wx4) * xd))
                + 4.0 * trilog(1.0i / (1.0i + xd)) - trilog((1.0i * (wx3 - xd)) / (wx3 * (1.0i + xd))) - trilog((1.0i * (wx4 - xd)) / (wx4 * (1.0i + xd))) - 2.0 * trilog(((-1.0i + wx3) * (-1.0i + xd)) / ((1.0i + wx3) * (1.0i + xd)))
                - 2.0 * trilog(((1.0i + wx3) * (-1.0i + xd)) / ((-1.0i + wx3) * (1.0i + xd))) - 2.0 * trilog(((-1.0i + wx4) * (-1.0i + xd)) / ((1.0i + wx4) * (1.0i + xd))) - 2.0 * trilog(((1.0i + wx4) * (-1.0i + xd)) / ((-1.0i + wx4) * (1.0i + xd))) + 5.0 * trilog((1.0i + xd) / (1.0i - wx3))
                + 5.0 * trilog((1.0i + xd) / (1.0i + wx3)) + 5.0 * trilog((1.0i + xd) / (1.0i - wx4)) + 5.0 * trilog((1.0i + xd) / (1.0i + wx4)) + 2.0 * trilog(-((1.0i + xd) / (wx3 - xd))) + 2.0 * trilog(-((1.0i + xd) / (wx4 - xd))) - 4.0 * trilog((1.0i + xd) * xdinv)
                + trilog((wx3 * (1.0i + xd)) / ((-1.0i + wx3) * xd)) + trilog((wx3 * (1.0i + xd)) / ((1.0i + wx3) * xd)) + trilog((wx4 * (1.0i + xd)) / ((-1.0i + wx4) * xd)) + trilog((wx4 * (1.0i + xd)) / ((1.0i + wx4) * xd))
                - 2.0 * trilog(((-1.0i + wx3) * (1.0i + xd)) / ((1.0i + wx3) * (-1.0i + xd))) - 2.0 * trilog(((1.0i + wx3) * (1.0i + xd)) / ((-1.0i + wx3) * (-1.0i + xd))) - 2.0 * trilog(((-1.0i + wx4) * (1.0i + xd)) / ((1.0i + wx4) * (-1.0i + xd)))
                - 2.0 * trilog(((1.0i + wx4) * (1.0i + xd)) / ((-1.0i + wx4) * (-1.0i + xd))) + 2.0 * trilog((-1.0i + xd) / (-wx3 + xd)) + trilog((-wx3 + xd) / (-1.0i + xd)) + trilog((-wx3 + xd) / (1.0i + xd)) + 2.0 * trilog((-1.0i + wx3) / (wx3 + xd))
                + 2.0 * trilog((1.0i + wx3) / (wx3 + xd)) - 2.0 * trilog(((1.0i + wx3) * (1.0 - 1.0i * xd)) / (2.0 * (wx3 + xd))) - 2.0 * trilog(((-1.0i + wx3) * (1.0 + 1.0i * xd)) / (2.0 * (wx3 + xd))) + 2.0 * trilog((-1.0i + xd) / (wx3 + xd)) + 2.0 * trilog((1.0i + xd) / (wx3 + xd))
                + trilog((wx3 + xd) / (-1.0i + wx3)) + trilog((wx3 + xd) / (1.0i + wx3)) - trilog(((-1.0i) * (wx3 + xd)) / ((-1.0i + wx3) * xd)) - trilog((1.0i * (wx3 + xd)) / ((1.0i + wx3) * xd)) + trilog((wx3 + xd) / (-1.0i + xd))
                - trilog(((-1.0i) * (wx3 + xd)) / (wx3 * (-1.0i + xd))) + trilog((wx3 + xd) / (1.0i + xd)) - trilog((1.0i * (wx3 + xd)) / (wx3 * (1.0i + xd))) + 2.0 * trilog((-1.0i + xd) / (-wx4 + xd)) + trilog((-wx4 + xd) / (-1.0i + xd))
                + trilog((-wx4 + xd) / (1.0i + xd)) + 2.0 * trilog((-1.0i + wx4) / (wx4 + xd)) + 2.0 * trilog((1.0i + wx4) / (wx4 + xd)) - 2.0 * trilog(((1.0i + wx4) * (1.0 - 1.0i * xd)) / (2.0 * (wx4 + xd))) - 2.0 * trilog(((-1.0i + wx4) * (1.0 + 1.0i * xd)) / (2.0 * (wx4 + xd)))
                + 2.0 * trilog((-1.0i + xd) / (wx4 + xd)) + 2.0 * trilog((1.0i + xd) / (wx4 + xd)) + trilog((wx4 + xd) / (-1.0i + wx4)) + trilog((wx4 + xd) / (1.0i + wx4)) - trilog(((-1.0i) * (wx4 + xd)) / ((-1.0i + wx4) * xd))
                - trilog((1.0i * (wx4 + xd)) / ((1.0i + wx4) * xd)) + trilog((wx4 + xd) / (-1.0i + xd)) - trilog(((-1.0i) * (wx4 + xd)) / (wx4 * (-1.0i + xd))) + trilog((wx4 + xd) / (1.0i + xd)) - trilog((1.0i * (wx4 + xd)) / (wx4 * (1.0i + xd)))
                - 4.0 * dilog((-1.0i + wx3) / (2.0 * wx3)) * lnhalf - 4.0 * dilog((1.0i - wx3) / (1.0i + wx3)) * lnhalf - 4.0 * dilog((1.0i + wx3) / (1.0i - wx3)) * lnhalf - 4.0 * dilog((1.0i + wx3) / (2.0 * wx3)) * lnhalf
                - 4.0 * dilog((-1.0i + wx4) / (2.0 * wx4)) * lnhalf - 4.0 * dilog((1.0i - wx4) / (1.0i + wx4)) * lnhalf - 4.0 * dilog((1.0i + wx4) / (1.0i - wx4)) * lnhalf - 4.0 * dilog((1.0i + wx4) / (2.0 * wx4)) * lnhalf - (2.0 * pisqu * log((-2.0 * 1.0i) / wx3)) / 3.0
                - (2.0 * power_of<3>(log((-2.0 * 1.0i) / wx3))) / 3.0 - (2.0 * pisqu * log((2.0 * 1.0i) / wx3)) / 3.0 - (2.0 * power_of<3>(log((2.0 * 1.0i) / wx3))) / 3.0 + 2.0 * power_of<2>(lnhalf) * log(wx3 / (-1.0i + wx3)) - 2.0 * power_of<2>(lnhalf) * log((2.0 * wx3) / (-1.0i + wx3))
                - 2.0 * power_of<2>(lnhalf) * log((-1.0i + wx3) / (2.0 * wx3)) + 2.0 * power_of<2>(lnhalf) * log((-1.0i + wx3) / wx3) + 2.0 * power_of<2>(lnhalf) * log(wx3 / (1.0i + wx3)) - 2.0 * power_of<2>(lnhalf) * log((2.0 * wx3) / (1.0i + wx3)) + (2.0 * pisqu * log((2.0 * (-1.0i + wx3)) / (1.0i + wx3))) / 3.0
                + (2.0 * power_of<3>(log((2.0 * (-1.0i + wx3)) / (1.0i + wx3)))) / 3.0 - 2.0 * power_of<2>(lnhalf) * log((1.0i + wx3) / (2.0 * wx3)) + 2.0 * power_of<2>(lnhalf) * log((1.0i + wx3) / wx3) + (2.0 * pisqu * log((2.0 * (1.0i + wx3)) / (-1.0i + wx3))) / 3.0
                + (2.0 * power_of<3>(log((2.0 * (1.0i + wx3)) / (-1.0i + wx3)))) / 3.0 - (2.0 * pisqu * log((-2.0 * 1.0i) / wx4)) / 3.0 - (2.0 * power_of<3>(log((-2.0 * 1.0i) / wx4))) / 3.0 - (2.0 * pisqu * log((2.0 * 1.0i) / wx4)) / 3.0 - (2.0 * power_of<3>(log((2.0 * 1.0i) / wx4))) / 3.0 + 2.0 * power_of<2>(lnhalf) * log(wx4 / (-1.0i + wx4))
                - 2.0 * power_of<2>(lnhalf) * log((2.0 * wx4) / (-1.0i + wx4)) - 2.0 * power_of<2>(lnhalf) * log((-1.0i + wx4) / (2.0 * wx4)) + 2.0 * power_of<2>(lnhalf) * log((-1.0i + wx4) / wx4) + 2.0 * power_of<2>(lnhalf) * log(wx4 / (1.0i + wx4)) - 2.0 * power_of<2>(lnhalf) * log((2.0 * wx4) / (1.0i + wx4))
                + (2.0 * pisqu * log((2.0 * (-1.0i + wx4)) / (1.0i + wx4))) / 3.0 + (2.0 * power_of<3>(log((2.0 * (-1.0i + wx4)) / (1.0i + wx4)))) / 3.0 - 2.0 * power_of<2>(lnhalf) * log((1.0i + wx4) / (2.0 * wx4)) + 2.0 * power_of<2>(lnhalf) * log((1.0i + wx4) / wx4)
                + (2.0 * pisqu * log((2.0 * (1.0i + wx4)) / (-1.0i + wx4))) / 3.0 + (2.0 * power_of<3>(log((2.0 * (1.0i + wx4)) / (-1.0i + wx4)))) / 3.0 + (pisqu * log(wx3 / (1.0i - xd))) / 6.0 + power_of<3>(log(wx3 / (1.0i - xd))) / 6.0 + (pisqu * log(wx4 / (1.0i - xd))) / 6.0 + power_of<3>(log(wx4 / (1.0i - xd))) / 6.0
                + (pisqu * log((-2.0 * 1.0i) / (wx3 - xd))) / 3.0 + power_of<3>(log((-2.0 * 1.0i) / (wx3 - xd))) / 3.0 + (pisqu * log((2.0 * 1.0i) / (wx3 - xd))) / 3.0 + power_of<3>(log((2.0 * 1.0i) / (wx3 - xd))) / 3.0 + (pisqu * log((-2.0 * 1.0i) / (wx4 - xd))) / 3.0 + power_of<3>(log((-2.0 * 1.0i) / (wx4 - xd))) / 3.0
                + (pisqu * log((2.0 * 1.0i) / (wx4 - xd))) / 3.0 + power_of<3>(log((2.0 * 1.0i) / (wx4 - xd))) / 3.0 + 2.0 * dilog(((-1.0i + wx3) * (1.0 - 1.0i * xd)) / (2.0 * (wx3 - xd))) * log(1.0 / 2.0 - (1.0i / 2.0) * xd)
                + 2.0 * dilog(((-1.0i + wx4) * (1.0 - 1.0i * xd)) / (2.0 * (wx4 - xd))) * log(1.0 / 2.0 - (1.0i / 2.0) * xd) - 2.0 * dilog(-((1.0i + xd) / (wx3 - xd))) * log(1.0 / 2.0 - (1.0i / 2.0) * xd) - 2.0 * dilog(-((1.0i + xd) / (wx4 - xd))) * log(1.0 / 2.0 - (1.0i / 2.0) * xd)
                + 2.0 * dilog(((-1.0i + wx3) * (1.0i + xd)) / ((1.0i + wx3) * (-1.0i + xd))) * log(1.0 / 2.0 - (1.0i / 2.0) * xd) + 2.0 * dilog(((1.0i + wx3) * (1.0i + xd)) / ((-1.0i + wx3) * (-1.0i + xd))) * log(1.0 / 2.0 - (1.0i / 2.0) * xd)
                + 2.0 * dilog(((-1.0i + wx4) * (1.0i + xd)) / ((1.0i + wx4) * (-1.0i + xd))) * log(1.0 / 2.0 - (1.0i / 2.0) * xd) + 2.0 * dilog(((1.0i + wx4) * (1.0i + xd)) / ((-1.0i + wx4) * (-1.0i + xd))) * log(1.0 / 2.0 - (1.0i / 2.0) * xd)
                + 2.0 * dilog(((1.0i + wx3) * (1.0 - 1.0i * xd)) / (2.0 * (wx3 + xd))) * log(1.0 / 2.0 - (1.0i / 2.0) * xd) - 2.0 * dilog((1.0i + xd) / (wx3 + xd)) * log(1.0 / 2.0 - (1.0i / 2.0) * xd) + 2.0 * dilog(((1.0i + wx4) * (1.0 - 1.0i * xd)) / (2.0 * (wx4 + xd))) * log(1.0 / 2.0 - (1.0i / 2.0) * xd)
                - 2.0 * dilog((1.0i + xd) / (wx4 + xd)) * log(1.0 / 2.0 - (1.0i / 2.0) * xd) - log((1.0i + wx3) / (wx3 - xd)) * power_of<2>(log(1.0 / 2.0 - (1.0i / 2.0) * xd)) + 2.0 * dilog(((1.0i + wx3) * (1.0 + 1.0i * xd)) / (2.0 * (wx3 - xd))) * log(1.0 / 2.0 + (1.0i / 2.0) * xd)
                + 2.0 * dilog(((1.0i + wx4) * (1.0 + 1.0i * xd)) / (2.0 * (wx4 - xd))) * log(1.0 / 2.0 + (1.0i / 2.0) * xd) + 2.0 * dilog(((-1.0i + wx3) * (-1.0i + xd)) / ((1.0i + wx3) * (1.0i + xd))) * log(1.0 / 2.0 + (1.0i / 2.0) * xd)
                + 2.0 * dilog(((1.0i + wx3) * (-1.0i + xd)) / ((-1.0i + wx3) * (1.0i + xd))) * log(1.0 / 2.0 + (1.0i / 2.0) * xd) + 2.0 * dilog(((-1.0i + wx4) * (-1.0i + xd)) / ((1.0i + wx4) * (1.0i + xd))) * log(1.0 / 2.0 + (1.0i / 2.0) * xd)
                + 2.0 * dilog(((1.0i + wx4) * (-1.0i + xd)) / ((-1.0i + wx4) * (1.0i + xd))) * log(1.0 / 2.0 + (1.0i / 2.0) * xd) - 2.0 * dilog((-1.0i + xd) / (-wx3 + xd)) * log(1.0 / 2.0 + (1.0i / 2.0) * xd)
                + 2.0 * dilog(((-1.0i + wx3) * (1.0 + 1.0i * xd)) / (2.0 * (wx3 + xd))) * log(1.0 / 2.0 + (1.0i / 2.0) * xd) - 2.0 * dilog((-1.0i + xd) / (wx3 + xd)) * log(1.0 / 2.0 + (1.0i / 2.0) * xd) - 2.0 * dilog((-1.0i + xd) / (-wx4 + xd)) * log(1.0 / 2.0 + (1.0i / 2.0) * xd)
                + 2.0 * dilog(((-1.0i + wx4) * (1.0 + 1.0i * xd)) / (2.0 * (wx4 + xd))) * log(1.0 / 2.0 + (1.0i / 2.0) * xd) - 2.0 * dilog((-1.0i + xd) / (wx4 + xd)) * log(1.0 / 2.0 + (1.0i / 2.0) * xd) - log((-1.0i + wx3) / (wx3 - xd)) * power_of<2>(log(1.0 / 2.0 + (1.0i / 2.0) * xd))
                + dilog((1.0i - xd) / (1.0i + wx3)) * (-2.0 * log(1.0 / 2.0 + (1.0i / 2.0) * xd) - log(1.0 - 1.0i * xd)) + dilog((1.0i - xd) / (1.0i + wx4)) * (-2.0 * log(1.0 / 2.0 + (1.0i / 2.0) * xd) - log(1.0 - 1.0i * xd))
                + dilog((-1.0i + xd) / (-1.0i + wx3)) * (-2.0 * log(1.0 / 2.0 + (1.0i / 2.0) * xd) - log(1.0 - 1.0i * xd)) + dilog((-1.0i + xd) / (-1.0i + wx4)) * (-2.0 * log(1.0 / 2.0 + (1.0i / 2.0) * xd) - log(1.0 - 1.0i * xd)) + 4.0 * dilog((1.0i + xd) * xdinv) * log(1.0 - 1.0i * xd)
                - dilog((wx3 * (1.0i + xd)) / ((-1.0i + wx3) * xd)) * log(1.0 - 1.0i * xd) - dilog((wx3 * (1.0i + xd)) / ((1.0i + wx3) * xd)) * log(1.0 - 1.0i * xd) - dilog((wx4 * (1.0i + xd)) / ((-1.0i + wx4) * xd)) * log(1.0 - 1.0i * xd)
                - dilog((wx4 * (1.0i + xd)) / ((1.0i + wx4) * xd)) * log(1.0 - 1.0i * xd) + dilog(-1.0i / (-1.0i + wx3)) * (4.0 * lnhalf - 2.0 * log(1.0 - 1.0i * xd) - 2.0 * log(1.0 + 1.0i * xd))
                + dilog(1.0i / (1.0i + wx3)) * (4.0 * lnhalf - 2.0 * log(1.0 - 1.0i * xd) - 2.0 * log(1.0 + 1.0i * xd)) + dilog(-1.0i / (-1.0i + wx4)) * (4.0 * lnhalf - 2.0 * log(1.0 - 1.0i * xd) - 2.0 * log(1.0 + 1.0i * xd))
                + dilog(1.0i / (1.0i + wx4)) * (4.0 * lnhalf - 2.0 * log(1.0 - 1.0i * xd) - 2.0 * log(1.0 + 1.0i * xd)) + dilog((1.0i + xd) / (1.0i - wx3)) * (-2.0 * log(1.0 / 2.0 - (1.0i / 2.0) * xd) - log(1.0 + 1.0i * xd))
                + dilog((1.0i + xd) / (1.0i + wx3)) * (-2.0 * log(1.0 / 2.0 - (1.0i / 2.0) * xd) - log(1.0 + 1.0i * xd)) + dilog((1.0i + xd) / (1.0i - wx4)) * (-2.0 * log(1.0 / 2.0 - (1.0i / 2.0) * xd) - log(1.0 + 1.0i * xd))
                + dilog((1.0i + xd) / (1.0i + wx4)) * (-2.0 * log(1.0 / 2.0 - (1.0i / 2.0) * xd) - log(1.0 + 1.0i * xd)) + dilog(-1.0i / wx3) * (4.0 * lnhalf - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd))
                + dilog(1.0i / wx3) * (4.0 * lnhalf - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd)) + dilog(-1.0i / wx4) * (4.0 * lnhalf - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd)) + dilog(1.0i / wx4) * (4.0 * lnhalf - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd))
                + 4.0 * dilog((-1.0i + xd) * xdinv) * log(1.0 + 1.0i * xd) - dilog((wx3 * (-1.0i + xd)) / ((-1.0i + wx3) * xd)) * log(1.0 + 1.0i * xd) - dilog((wx3 * (-1.0i + xd)) / ((1.0i + wx3) * xd)) * log(1.0 + 1.0i * xd)
                - dilog((wx4 * (-1.0i + xd)) / ((-1.0i + wx4) * xd)) * log(1.0 + 1.0i * xd) - dilog((wx4 * (-1.0i + xd)) / ((1.0i + wx4) * xd)) * log(1.0 + 1.0i * xd) + dilog((-2.0 * 1.0i) / (-1.0i + wx3)) * (2.0 * log(1.0 - 1.0i * xd) + 2.0 * log(1.0 + 1.0i * xd))
                + dilog((2.0 * 1.0i) / (1.0i + wx3)) * (2.0 * log(1.0 - 1.0i * xd) + 2.0 * log(1.0 + 1.0i * xd)) + dilog((-2.0 * 1.0i) / (-1.0i + wx4)) * (2.0 * log(1.0 - 1.0i * xd) + 2.0 * log(1.0 + 1.0i * xd)) + dilog((2.0 * 1.0i) / (1.0i + wx4)) * (2.0 * log(1.0 - 1.0i * xd) + 2.0 * log(1.0 + 1.0i * xd))
                + log((wx3 - xd) / (1.0i + wx3)) * (power_of<2>(log(1.0 - 1.0i * xd)) / 2.0 + log(1.0 - 1.0i * xd) * log(1.0 + 1.0i * xd)) + log((wx4 - xd) / (1.0i + wx4)) * (power_of<2>(log(1.0 - 1.0i * xd)) / 2.0 + log(1.0 - 1.0i * xd) * log(1.0 + 1.0i * xd))
                + log((wx3 - xd) / (-1.0i + wx3)) * (log(1.0 - 1.0i * xd) * log(1.0 + 1.0i * xd) + power_of<2>(log(1.0 + 1.0i * xd)) / 2.0) + log((wx4 - xd) / (-1.0i + wx4)) * (log(1.0 - 1.0i * xd) * log(1.0 + 1.0i * xd) + power_of<2>(log(1.0 + 1.0i * xd)) / 2.0) - (2.0 * pisqu * log(-1.0i / xd)) / 3.0
                - (2.0 * power_of<3>(log(-1.0i / xd))) / 3.0 - (2.0 * pisqu * log(1.0i / xd)) / 3.0 - (2.0 * power_of<3>(log(1.0i / xd))) / 3.0 + dilog(-(xd / wx3)) * (-2.0 * ln2 - log(xd)) + dilog(xd / wx3) * (-2.0 * ln2 - log(xd))
                + dilog(-(xd / wx4)) * (-2.0 * ln2 - log(xd)) + dilog(xd / wx4) * (-2.0 * ln2 - log(xd)) + 2.0 * dilog((-1.0i) * xd) * log(xd) + 2.0 * dilog(1.0i * xd) * log(xd) + dilog(wx3 / (-1.0i + wx3)) * (4.0 * ln2 + 2.0 * log(xd))
                + dilog(wx3 / (1.0i + wx3)) * (4.0 * ln2 + 2.0 * log(xd)) + dilog(wx4 / (-1.0i + wx4)) * (4.0 * ln2 + 2.0 * log(xd)) + dilog(wx4 / (1.0i + wx4)) * (4.0 * ln2 + 2.0 * log(xd)) + (pisqu * log(wx3 / (-1.0i + xd))) / 6.0
                + power_of<3>(log(wx3 / (-1.0i + xd))) / 6.0 - (pisqu * log((2.0 - (2.0 * 1.0i) * wx3) / ((-1.0i + wx3) * (-1.0i + xd)))) / 3.0 - power_of<3>(log((2.0 - (2.0 * 1.0i) * wx3) / ((-1.0i + wx3) * (-1.0i + xd)))) / 3.0 - (pisqu * log((-2.0 - (2.0 * 1.0i) * wx3) / ((1.0i + wx3) * (-1.0i + xd)))) / 3.0
                - power_of<3>(log((-2.0 - (2.0 * 1.0i) * wx3) / ((1.0i + wx3) * (-1.0i + xd)))) / 3.0 + (pisqu * log(wx4 / (-1.0i + xd))) / 6.0 + power_of<3>(log(wx4 / (-1.0i + xd))) / 6.0 - (pisqu * log((2.0 - (2.0 * 1.0i) * wx4) / ((-1.0i + wx4) * (-1.0i + xd)))) / 3.0
                - power_of<3>(log((2.0 - (2.0 * 1.0i) * wx4) / ((-1.0i + wx4) * (-1.0i + xd)))) / 3.0 - (pisqu * log((-2.0 - (2.0 * 1.0i) * wx4) / ((1.0i + wx4) * (-1.0i + xd)))) / 3.0 - power_of<3>(log((-2.0 - (2.0 * 1.0i) * wx4) / ((1.0i + wx4) * (-1.0i + xd)))) / 3.0 + (pisqu * log(-(wx3 / (1.0i + xd)))) / 6.0
                + power_of<3>(log(-(wx3 / (1.0i + xd)))) / 6.0 + (pisqu * log(wx3 / (1.0i + xd))) / 6.0 + power_of<3>(log(wx3 / (1.0i + xd))) / 6.0 - (pisqu * log((2.0 + (2.0 * 1.0i) * wx3) / ((1.0i + wx3) * (1.0i + xd)))) / 3.0 - power_of<3>(log((2.0 + (2.0 * 1.0i) * wx3) / ((1.0i + wx3) * (1.0i + xd)))) / 3.0
                - (pisqu * log(((2.0 * 1.0i) * (1.0i + wx3)) / ((-1.0i + wx3) * (1.0i + xd)))) / 3.0 - power_of<3>(log(((2.0 * 1.0i) * (1.0i + wx3)) / ((-1.0i + wx3) * (1.0i + xd)))) / 3.0 + (pisqu * log(-(wx4 / (1.0i + xd)))) / 6.0 + power_of<3>(log(-(wx4 / (1.0i + xd)))) / 6.0 + (pisqu * log(wx4 / (1.0i + xd))) / 6.0
                + power_of<3>(log(wx4 / (1.0i + xd))) / 6.0 - (pisqu * log((2.0 + (2.0 * 1.0i) * wx4) / ((1.0i + wx4) * (1.0i + xd)))) / 3.0 - power_of<3>(log((2.0 + (2.0 * 1.0i) * wx4) / ((1.0i + wx4) * (1.0i + xd)))) / 3.0 - (pisqu * log(((2.0 * 1.0i) * (1.0i + wx4)) / ((-1.0i + wx4) * (1.0i + xd)))) / 3.0
                - power_of<3>(log(((2.0 * 1.0i) * (1.0i + wx4)) / ((-1.0i + wx4) * (1.0i + xd)))) / 3.0 + (pisqu * log((-2.0 * 1.0i) / (wx3 + xd))) / 3.0 + power_of<3>(log((-2.0 * 1.0i) / (wx3 + xd))) / 3.0 + (pisqu * log((2.0 * 1.0i) / (wx3 + xd))) / 3.0 + power_of<3>(log((2.0 * 1.0i) / (wx3 + xd))) / 3.0
                + dilog((wx3 + xd) / (-1.0i + wx3)) * (-2.0 * ln2 - log(xd) - log((wx3 + xd) / wx3)) + dilog((wx3 + xd) / (1.0i + wx3)) * (-2.0 * ln2 - log(xd) - log((wx3 + xd) / wx3))
                + dilog(((-1.0i) * (wx3 + xd)) / ((-1.0i + wx3) * xd)) * log((wx3 + xd) / wx3) + dilog((1.0i * (wx3 + xd)) / ((1.0i + wx3) * xd)) * log((wx3 + xd) / wx3) - dilog((wx3 + xd) / (-1.0i + xd)) * log((wx3 + xd) / wx3)
                + dilog(((-1.0i) * (wx3 + xd)) / (wx3 * (-1.0i + xd))) * log((wx3 + xd) / wx3) - dilog((wx3 + xd) / (1.0i + xd)) * log((wx3 + xd) / wx3) + dilog((1.0i * (wx3 + xd)) / (wx3 * (1.0i + xd))) * log((wx3 + xd) / wx3)
                + (-2.0 * ln2 * log((1.0i + xd) / (1.0i - wx3)) - log(xd) * log((1.0i + xd) / (1.0i - wx3))) * log((wx3 + xd) / wx3) + ((-1.0 / 2.0) * log((1.0i + wx3) / (1.0i - xd)) + log(((1.0i + wx3) * xd) / (wx3 * (-1.0i + xd))) / 2.0
                    + log((wx3 * (-1.0i + xd)) / ((1.0i + wx3) * xd)) / 2.0 - log((1.0i - wx3) / (1.0i + xd)) / 2.0 + log(((-1.0i + wx3) * xd) / (wx3 * (1.0i + xd))) / 2.0 - log((1.0i + xd) / (1.0i - wx3)) / 2.0 + log((wx3 * (1.0i + xd)) / ((-1.0i + wx3) * xd)) / 2.0) *
                power_of<2>(log((wx3 + xd) / wx3)) + log((1.0i - xd) / (1.0i + wx3)) * ((-2.0 * ln2 - log(xd)) * log((wx3 + xd) / wx3) - power_of<2>(log((wx3 + xd) / wx3)) / 2.0) + (pisqu * log((-2.0 * 1.0i) / (wx4 + xd))) / 3.0 + power_of<3>(log((-2.0 * 1.0i) / (wx4 + xd))) / 3.0
                + (pisqu * log((2.0 * 1.0i) / (wx4 + xd))) / 3.0 + power_of<3>(log((2.0 * 1.0i) / (wx4 + xd))) / 3.0 + dilog((wx4 + xd) / (-1.0i + wx4)) * (-2.0 * ln2 - log(xd) - log((wx4 + xd) / wx4))
                + dilog((wx4 + xd) / (1.0i + wx4)) * (-2.0 * ln2 - log(xd) - log((wx4 + xd) / wx4)) + dilog(((-1.0i) * (wx4 + xd)) / ((-1.0i + wx4) * xd)) * log((wx4 + xd) / wx4)
                + dilog((1.0i * (wx4 + xd)) / ((1.0i + wx4) * xd)) * log((wx4 + xd) / wx4) - dilog((wx4 + xd) / (-1.0i + xd)) * log((wx4 + xd) / wx4) + dilog(((-1.0i) * (wx4 + xd)) / (wx4 * (-1.0i + xd))) * log((wx4 + xd) / wx4)
                - dilog((wx4 + xd) / (1.0i + xd)) * log((wx4 + xd) / wx4) + dilog((1.0i * (wx4 + xd)) / (wx4 * (1.0i + xd))) * log((wx4 + xd) / wx4) - 2.0 * ln2 * log((1.0i + xd) / (1.0i - wx4)) * log((wx4 + xd) / wx4)
                + ((-1.0 / 2.0) * log((1.0i + wx4) / (1.0i - xd)) + log(((1.0i + wx4) * xd) / (wx4 * (-1.0i + xd))) / 2.0 + log((wx4 * (-1.0i + xd)) / ((1.0i + wx4) * xd)) / 2.0 - log((1.0i - wx4) / (1.0i + xd)) / 2.0 + log(((-1.0i + wx4) * xd) / (wx4 * (1.0i + xd))) / 2.0
                    - log((1.0i + xd) / (1.0i - wx4)) / 2.0 + log((wx4 * (1.0i + xd)) / ((-1.0i + wx4) * xd)) / 2.0) * power_of<2>(log((wx4 + xd) / wx4)) + log((1.0i - xd) / (1.0i + wx4)) * (-2.0 * ln2 * log((wx4 + xd) / wx4) - power_of<2>(log((wx4 + xd) / wx4)) / 2.0)
                + log(1.0 / 2.0 - (1.0i / 2.0) * xd) * (-2.0 * log((wx3 - xd) / (1.0i + wx3)) * log(1.0 + 1.0i * xd) - 2.0 * log((wx4 - xd) / (1.0i + wx4)) * log(1.0 + 1.0i * xd) + log(1.0 + 1.0i * xd) * (-2.0 * log((wx3 + xd) / (-1.0i + wx3)) - 2.0 * log((wx4 + xd) / (-1.0i + wx4))))
                + log(1.0 / 2.0 + (1.0i / 2.0) * xd) * (-2.0 * log((wx3 - xd) / (-1.0i + wx3)) * log(1.0 - 1.0i * xd) - 2.0 * log((wx4 - xd) / (-1.0i + wx4)) * log(1.0 - 1.0i * xd) + log(1.0 - 1.0i * xd) * (-2.0 * log((wx3 + xd) / (1.0i + wx3)) - 2.0 * log((wx4 + xd) / (1.0i + wx4))))
                + power_of<2>(log(1.0 / 2.0 - (1.0i / 2.0) * xd)) * (-log((wx3 - xd) / (1.0i + wx3)) - log((1.0i + wx4) / (wx4 - xd)) - log((wx4 - xd) / (1.0i + wx4)) + log(((1.0i + wx3) * (1.0 + 1.0i * xd)) / (2.0 * (wx3 - xd)))
                    + log(((1.0i + wx4) * (1.0 + 1.0i * xd)) / (2.0 * (wx4 - xd))) + log(((-2.0 * 1.0i) * (wx3 - xd)) / ((1.0i + wx3) * (-1.0i + xd))) + log(((-2.0 * 1.0i) * (wx4 - xd)) / ((1.0i + wx4) * (-1.0i + xd))) - log((-1.0i + wx3) / (wx3 + xd))
                    + log(((-1.0i + wx3) * (1.0 + 1.0i * xd)) / (2.0 * (wx3 + xd))) - log((wx3 + xd) / (-1.0i + wx3)) + log(((-2.0 * 1.0i) * (wx3 + xd)) / ((-1.0i + wx3) * (-1.0i + xd))) - log((-1.0i + wx4) / (wx4 + xd))
                    + log(((-1.0i + wx4) * (1.0 + 1.0i * xd)) / (2.0 * (wx4 + xd))) - log((wx4 + xd) / (-1.0i + wx4)) + log(((-2.0 * 1.0i) * (wx4 + xd)) / ((-1.0i + wx4) * (-1.0i + xd))))
                + power_of<2>(log(1.0 / 2.0 + (1.0i / 2.0) * xd)) * (-log((wx3 - xd) / (-1.0i + wx3)) - log((-1.0i + wx4) / (wx4 - xd)) - log((wx4 - xd) / (-1.0i + wx4)) + log(((-1.0i + wx3) * (1.0 - 1.0i * xd)) / (2.0 * (wx3 - xd)))
                    + log(((-1.0i + wx4) * (1.0 - 1.0i * xd)) / (2.0 * (wx4 - xd))) + log(((2.0 * 1.0i) * (wx3 - xd)) / ((-1.0i + wx3) * (1.0i + xd))) + log(((2.0 * 1.0i) * (wx4 - xd)) / ((-1.0i + wx4) * (1.0i + xd))) - log((1.0i + wx3) / (wx3 + xd))
                    + log(((1.0i + wx3) * (1.0 - 1.0i * xd)) / (2.0 * (wx3 + xd))) - log((wx3 + xd) / (1.0i + wx3)) + log(((2.0 * 1.0i) * (wx3 + xd)) / ((1.0i + wx3) * (1.0i + xd))) - log((1.0i + wx4) / (wx4 + xd)) + log(((1.0i + wx4) * (1.0 - 1.0i * xd)) / (2.0 * (wx4 + xd)))
                    - log((wx4 + xd) / (1.0i + wx4)) + log(((2.0 * 1.0i) * (wx4 + xd)) / ((1.0i + wx4) * (1.0i + xd)))) + dilog((wx3 - xd) / (-1.0i + wx3)) * (-2.0 * ln2 - log(xd) - log(1.0 - xd / wx3))
                + dilog((wx3 - xd) / (1.0i + wx3)) * (-2.0 * ln2 - log(xd) - log(1.0 - xd / wx3)) + dilog((1.0i * (wx3 - xd)) / ((-1.0i + wx3) * xd)) * log(1.0 - xd / wx3) + dilog(((-1.0i) * (wx3 - xd)) / ((1.0i + wx3) * xd)) * log(1.0 - xd / wx3)
                + dilog(((-1.0i) * (wx3 - xd)) / (wx3 * (-1.0i + xd))) * log(1.0 - xd / wx3) + dilog((1.0i * (wx3 - xd)) / (wx3 * (1.0i + xd))) * log(1.0 - xd / wx3) - dilog((-wx3 + xd) / (-1.0i + xd)) * log(1.0 - xd / wx3)
                - dilog((-wx3 + xd) / (1.0i + xd)) * log(1.0 - xd / wx3) - 2.0 * ln2 * log((1.0i + xd) / (1.0i + wx3)) * log(1.0 - xd / wx3) + ((-1.0 / 2.0) * log((-1.0i + wx3) / (-1.0i + xd)) + log(((-1.0i + wx3) * xd) / (wx3 * (-1.0i + xd))) / 2.0
                    + log((wx3 * (-1.0i + xd)) / ((-1.0i + wx3) * xd)) / 2.0 - log((1.0i + wx3) / (1.0i + xd)) / 2.0 + log(((1.0i + wx3) * xd) / (wx3 * (1.0i + xd))) / 2.0 - log((1.0i + xd) / (1.0i + wx3)) / 2.0 + log((wx3 * (1.0i + xd)) / ((1.0i + wx3) * xd)) / 2.0) *
                power_of<2>(log(1.0 - xd / wx3)) + log((-1.0i + xd) / (-1.0i + wx3)) * (-2.0 * ln2 * log(1.0 - xd / wx3) - power_of<2>(log(1.0 - xd / wx3)) / 2.0) + dilog((wx4 - xd) / (-1.0i + wx4)) * (-2.0 * ln2 - log(xd) - log(1.0 - xd / wx4))
                + dilog((wx4 - xd) / (1.0i + wx4)) * (-2.0 * ln2 - log(xd) - log(1.0 - xd / wx4)) + dilog((1.0i * (wx4 - xd)) / ((-1.0i + wx4) * xd)) * log(1.0 - xd / wx4) + dilog(((-1.0i) * (wx4 - xd)) / ((1.0i + wx4) * xd)) * log(1.0 - xd / wx4)
                + dilog(((-1.0i) * (wx4 - xd)) / (wx4 * (-1.0i + xd))) * log(1.0 - xd / wx4) + dilog((1.0i * (wx4 - xd)) / (wx4 * (1.0i + xd))) * log(1.0 - xd / wx4) - dilog((-wx4 + xd) / (-1.0i + xd)) * log(1.0 - xd / wx4)
                - dilog((-wx4 + xd) / (1.0i + xd)) * log(1.0 - xd / wx4) - 2.0 * ln2 * log((1.0i + xd) / (1.0i + wx4)) * log(1.0 - xd / wx4) + ((-1.0 / 2.0) * log((-1.0i + wx4) / (-1.0i + xd)) + log(((-1.0i + wx4) * xd) / (wx4 * (-1.0i + xd))) / 2.0
                    + log((wx4 * (-1.0i + xd)) / ((-1.0i + wx4) * xd)) / 2.0 - log((1.0i + wx4) / (1.0i + xd)) / 2.0 + log(((1.0i + wx4) * xd) / (wx4 * (1.0i + xd))) / 2.0 - log((1.0i + xd) / (1.0i + wx4)) / 2.0 + log((wx4 * (1.0i + xd)) / ((1.0i + wx4) * xd)) / 2.0) *
                power_of<2>(log(1.0 - xd / wx4)) + log((-1.0i + xd) / (-1.0i + wx4)) * (-2.0 * ln2 * log(1.0 - xd / wx4) - power_of<2>(log(1.0 - xd / wx4)) / 2.0) - 1.0i * M_PI * H1(1.0i * wx3, -(wx3 * xdinv)) * power_of<2>(log((-1.0i + xd) / wx3)) * my_sign(-imag(wx3 * xdinv))
                - 1.0i * M_PI * H1((-1.0i) * wx3, -(wx3 * xdinv)) * power_of<2>(log((1.0i + xd) / wx3)) * my_sign(-imag(wx3 * xdinv)) - 1.0i * M_PI * H1((-1.0i) * wx3, wx3 * xdinv) * power_of<2>(log((1.0i - xd) / wx3)) * my_sign(imag(wx3 * xdinv))
                - 1.0i * M_PI * H1(1.0i * wx3, wx3 * xdinv) * power_of<2>(log(-((1.0i + xd) / wx3))) * my_sign(imag(wx3 * xdinv)) - 1.0i * M_PI * H1(1.0i * wx4, -(wx4 * xdinv)) * power_of<2>(log((-1.0i + xd) / wx4)) * my_sign(-imag(wx4 * xdinv))
                - 1.0i * M_PI * H1((-1.0i) * wx4, -(wx4 * xdinv)) * power_of<2>(log((1.0i + xd) / wx4)) * my_sign(-imag(wx4 * xdinv)) - 1.0i * M_PI * H1((-1.0i) * wx4, wx4 * xdinv) * power_of<2>(log((1.0i - xd) / wx4)) * my_sign(imag(wx4 * xdinv))
                - 1.0i * M_PI * H1(1.0i * wx4, wx4 * xdinv) * power_of<2>(log(-((1.0i + xd) / wx4))) * my_sign(imag(wx4 * xdinv)) - (2.0 * 1.0i) * M_PI * H1((-2.0 * 1.0i) / (-1.0i + wx3), (-2.0 * 1.0i) / (-1.0i + xd)) * power_of<2>(log((1.0i / 2.0) * (wx3 - xd))) * my_sign(2.0 * real(1.0 / (1.0i - xd)))
                - (2.0 * 1.0i) * M_PI * H1((-2.0 * 1.0i) / (-1.0i + wx4), (-2.0 * 1.0i) / (-1.0i + xd)) * power_of<2>(log((1.0i / 2.0) * (wx4 - xd))) * my_sign(2.0 * real(1.0 / (1.0i - xd))) + (2.0 * 1.0i) * M_PI * H1((-1.0i + wx3) / (-1.0i + xd), (-2.0 * 1.0i) / (-1.0i + xd)) * power_of<2>(log(((1.0i + wx3) * (1.0 + 1.0i * xd)) / (2.0 * 1.0i - 2.0 * wx3))) *
                my_sign(2.0 * real(1.0 / (1.0i - xd))) + (2.0 * 1.0i) * M_PI * H1((-1.0i + wx4) / (-1.0i + xd), (-2.0 * 1.0i) / (-1.0i + xd)) * power_of<2>(log(((1.0i + wx4) * (1.0 + 1.0i * xd)) / (2.0 * 1.0i - 2.0 * wx4))) * my_sign(2.0 * real(1.0 / (1.0i - xd)))
                + (2.0 * 1.0i) * M_PI * H1((1.0i + wx3) / (1.0i - xd), (-2.0 * 1.0i) / (-1.0i + xd)) * power_of<2>(log((((-1.0 / 2.0) * 1.0i) * (-1.0i + wx3) * (-1.0i + xd)) / (1.0i + wx3))) * my_sign(2.0 * real(1.0 / (1.0i - xd))) + (2.0 * 1.0i) * M_PI * H1((1.0i + wx4) / (1.0i - xd), (-2.0 * 1.0i) / (-1.0i + xd)) *
                power_of<2>(log((((-1.0 / 2.0) * 1.0i) * (-1.0i + wx4) * (-1.0i + xd)) / (1.0i + wx4))) * my_sign(2.0 * real(1.0 / (1.0i - xd))) - (2.0 * 1.0i) * M_PI * H1((2.0 * 1.0i) / (1.0i + wx3), (-2.0 * 1.0i) / (-1.0i + xd)) * power_of<2>(log(((-1.0 / 2.0) * 1.0i) * (wx3 + xd))) * my_sign(2.0 * real(1.0 / (1.0i - xd)))
                - (2.0 * 1.0i) * M_PI * H1((2.0 * 1.0i) / (1.0i + wx4), (-2.0 * 1.0i) / (-1.0i + xd)) * power_of<2>(log(((-1.0 / 2.0) * 1.0i) * (wx4 + xd))) * my_sign(2.0 * real(1.0 / (1.0i - xd)))
                + power_of<2>(log((-1.0i - wx3inv) * xd)) * (1.0i * M_PI * H1(1.0i / xd, wx3 * xdinv) * my_sign(imag(wx3 * xdinv)) - 1.0i * M_PI * H1(-(wx3 * xdinv), -1.0i / xd) * my_sign(-real(xdinv)))
                + power_of<2>(log((-1.0i + wx3inv) * xd)) * (1.0i * M_PI * H1(1.0i / xd, -(wx3 * xdinv)) * my_sign(-imag(wx3 * xdinv)) - 1.0i * M_PI * H1(wx3 * xdinv, -1.0i / xd) * my_sign(-real(xdinv)))
                + power_of<2>(log((-1.0i - wx4inv) * xd)) * (1.0i * M_PI * H1(1.0i / xd, wx4 * xdinv) * my_sign(imag(wx4 * xdinv)) - 1.0i * M_PI * H1(-(wx4 * xdinv), -1.0i / xd) * my_sign(-real(xdinv)))
                + power_of<2>(log((-1.0i + wx4inv) * xd)) * (1.0i * M_PI * H1(1.0i / xd, -(wx4 * xdinv)) * my_sign(-imag(wx4 * xdinv)) - 1.0i * M_PI * H1(wx4 * xdinv, -1.0i / xd) * my_sign(-real(xdinv)))
                + power_of<2>(log((1.0i - wx3inv) * xd)) * (1.0i * M_PI * H1(-1.0i / xd, wx3 * xdinv) * my_sign(imag(wx3 * xdinv)) - 1.0i * M_PI * H1(-(wx3 * xdinv), 1.0i / xd) * my_sign(real(xdinv)))
                + power_of<2>(log((1.0i + wx3inv) * xd)) * (1.0i * M_PI * H1(-1.0i / xd, -(wx3 * xdinv)) * my_sign(-imag(wx3 * xdinv)) - 1.0i * M_PI * H1(wx3 * xdinv, 1.0i / xd) * my_sign(real(xdinv)))
                + power_of<2>(log((1.0i - wx4inv) * xd)) * (1.0i * M_PI * H1(-1.0i / xd, wx4 * xdinv) * my_sign(imag(wx4 * xdinv)) - 1.0i * M_PI * H1(-(wx4 * xdinv), 1.0i / xd) * my_sign(real(xdinv)))
                + power_of<2>(log((1.0i + wx4inv) * xd)) * (1.0i * M_PI * H1(-1.0i / xd, -(wx4 * xdinv)) * my_sign(-imag(wx4 * xdinv)) - 1.0i * M_PI * H1(wx4 * xdinv, 1.0i / xd) * my_sign(real(xdinv))) - (2.0 * 1.0i) * M_PI * H1((2.0 * 1.0i) / (1.0i + wx3), (2.0 * 1.0i) / (1.0i + xd)) * power_of<2>(log(((-1.0 / 2.0) * 1.0i) * (wx3 - xd))) *
                my_sign(2.0 * real(1.0 / (1.0i + xd))) - (2.0 * 1.0i) * M_PI * H1((2.0 * 1.0i) / (1.0i + wx4), (2.0 * 1.0i) / (1.0i + xd)) * power_of<2>(log(((-1.0 / 2.0) * 1.0i) * (wx4 - xd))) * my_sign(2.0 * real(1.0 / (1.0i + xd)))
                + (2.0 * 1.0i) * M_PI * H1((1.0i + wx3) / (1.0i + xd), (2.0 * 1.0i) / (1.0i + xd)) * power_of<2>(log(((1.0 + 1.0i * wx3) * (1.0i + xd)) / (2.0 * (1.0i + wx3)))) * my_sign(2.0 * real(1.0 / (1.0i + xd))) + (2.0 * 1.0i) * M_PI * H1((1.0i + wx4) / (1.0i + xd), (2.0 * 1.0i) / (1.0i + xd)) *
                power_of<2>(log(((1.0 + 1.0i * wx4) * (1.0i + xd)) / (2.0 * (1.0i + wx4)))) * my_sign(2.0 * real(1.0 / (1.0i + xd))) - (2.0 * 1.0i) * M_PI * H1((-2.0 * 1.0i) / (-1.0i + wx3), (2.0 * 1.0i) / (1.0i + xd)) * power_of<2>(log((1.0i / 2.0) * (wx3 + xd))) * my_sign(2.0 * real(1.0 / (1.0i + xd)))
                - (2.0 * 1.0i) * M_PI * H1((-2.0 * 1.0i) / (-1.0i + wx4), (2.0 * 1.0i) / (1.0i + xd)) * power_of<2>(log((1.0i / 2.0) * (wx4 + xd))) * my_sign(2.0 * real(1.0 / (1.0i + xd))) + (2.0 * 1.0i) * M_PI * H1((1.0i - wx3) / (1.0i + xd), (2.0 * 1.0i) / (1.0i + xd)) * power_of<2>(log((1.0i + wx3 + xd - 1.0i * wx3 * xd) / (2.0 * 1.0i - 2.0 * wx3))) *
                my_sign(2.0 * real(1.0 / (1.0i + xd))) + (2.0 * 1.0i) * M_PI * H1((1.0i - wx4) / (1.0i + xd), (2.0 * 1.0i) / (1.0i + xd)) * power_of<2>(log((1.0i + wx4 + xd - 1.0i * wx4 * xd) / (2.0 * 1.0i - 2.0 * wx4))) * my_sign(2.0 * real(1.0 / (1.0i + xd)))
                + (2.0 * 1.0i) * M_PI * power_of<2>(log(1.0 / 2.0 - (1.0i / 2.0) * wx3)) * my_sign(real(wx3) / 2.0) * T(1.0, 1.0 / 2.0 - (1.0i / 2.0) * xd, 1.0 / 2.0 - (1.0i / 2.0) * wx3) + (2.0 * 1.0i) * M_PI * power_of<2>(log(1.0 / 2.0 + (1.0i / 2.0) * wx3)) * my_sign((-1.0 / 2.0) * real(wx3)) * T(1.0, 1.0 / 2.0 - (1.0i / 2.0) * xd, 1.0 / 2.0 + (1.0i / 2.0) * wx3)
                + (2.0 * 1.0i) * M_PI * power_of<2>(log(((-1.0 / 2.0) * 1.0i) * (1.0i + wx4))) * my_sign(real(wx4) / 2.0) * T(1.0, 1.0 / 2.0 - (1.0i / 2.0) * xd, 1.0 / 2.0 - (1.0i / 2.0) * wx4) + (2.0 * 1.0i) * M_PI * power_of<2>(log(1.0 / 2.0 + (1.0i / 2.0) * wx4)) * my_sign((-1.0 / 2.0) * real(wx4)) * T(1.0, 1.0 / 2.0 - (1.0i / 2.0) * xd, 1.0 / 2.0 + (1.0i / 2.0) * wx4)
                - (2.0 * 1.0i) * M_PI * power_of<2>(log((wx3 + xd) / (1.0i + wx3))) * my_sign(imag((1.0i - xd) / (1.0i + wx3))) * T(1.0, 1.0 / 2.0 - (1.0i / 2.0) * xd, (wx3 + xd) / (1.0i + wx3)) - (2.0 * 1.0i) * M_PI * power_of<2>(log((wx4 + xd) / (1.0i + wx4))) * my_sign(imag((1.0i - xd) / (1.0i + wx4))) *
                T(1.0, 1.0 / 2.0 - (1.0i / 2.0) * xd, (wx4 + xd) / (1.0i + wx4)) + (2.0 * 1.0i) * M_PI * power_of<2>(log(1.0 / 2.0 - (1.0i / 2.0) * wx3)) * my_sign(real(wx3) / 2.0) * T(1.0, 1.0 / 2.0 + (1.0i / 2.0) * xd, 1.0 / 2.0 - (1.0i / 2.0) * wx3)
                + (2.0 * 1.0i) * M_PI * power_of<2>(log(1.0 / 2.0 + (1.0i / 2.0) * wx3)) * my_sign((-1.0 / 2.0) * real(wx3)) * T(1.0, 1.0 / 2.0 + (1.0i / 2.0) * xd, 1.0 / 2.0 + (1.0i / 2.0) * wx3) + (2.0 * 1.0i) * M_PI * power_of<2>(log(((-1.0 / 2.0) * 1.0i) * (1.0i + wx4))) * my_sign(real(wx4) / 2.0) * T(1.0, 1.0 / 2.0 + (1.0i / 2.0) * xd, 1.0 / 2.0 - (1.0i / 2.0) * wx4)
                + (2.0 * 1.0i) * M_PI * power_of<2>(log(1.0 / 2.0 + (1.0i / 2.0) * wx4)) * my_sign((-1.0 / 2.0) * real(wx4)) * T(1.0, 1.0 / 2.0 + (1.0i / 2.0) * xd, 1.0 / 2.0 + (1.0i / 2.0) * wx4) - (2.0 * 1.0i) * M_PI * power_of<2>(log((wx3 - xd) / (1.0i + wx3))) * my_sign(imag((1.0i + xd) / (1.0i + wx3))) *
                T(1.0, 1.0 / 2.0 + (1.0i / 2.0) * xd, (wx3 - xd) / (1.0i + wx3)) - (2.0 * 1.0i) * M_PI * power_of<2>(log((wx4 - xd) / (1.0i + wx4))) * my_sign(imag((1.0i + xd) / (1.0i + wx4))) * T(1.0, 1.0 / 2.0 + (1.0i / 2.0) * xd, (wx4 - xd) / (1.0i + wx4))
                + (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 + (1.0i / 2.0) * wx3) * log(1.0 - 1.0i * wx3) * my_sign(real(wx3)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) + 1.0i * M_PI * power_of<2>(log(1.0 - 1.0i * wx3)) * my_sign(real(wx3)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3)
                + (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 - (1.0i / 2.0) * wx3) * log(1.0 + 1.0i * wx3) * my_sign(-real(wx3)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) + 1.0i * M_PI * power_of<2>(log(1.0 + 1.0i * wx3)) * my_sign(-real(wx3)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                + (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 + (1.0i / 2.0) * wx4) * log(1.0 - 1.0i * wx4) * my_sign(real(wx4)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) + 1.0i * M_PI * power_of<2>(log(1.0 - 1.0i * wx4)) * my_sign(real(wx4)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4)
                + (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 - (1.0i / 2.0) * wx4) * log(1.0 + 1.0i * wx4) * my_sign(-real(wx4)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4) + 1.0i * M_PI * power_of<2>(log(1.0 + 1.0i * wx4)) * my_sign(-real(wx4)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                - (2.0 * 1.0i) * M_PI * power_of<2>(log((wx3 + xd) / (-1.0i + wx3))) * my_sign(imag((1.0i + xd) / (1.0i - wx3))) * T(1.0, (1.0 + 1.0i * xd) / 2.0, (wx3 + xd) / (-1.0i + wx3)) - (2.0 * 1.0i) * M_PI * power_of<2>(log((wx4 + xd) / (-1.0i + wx4))) * my_sign(imag((1.0i + xd) / (1.0i - wx4))) *
                T(1.0, (1.0 + 1.0i * xd) / 2.0, (wx4 + xd) / (-1.0i + wx4)) + (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 + (1.0i / 2.0) * wx3) * log(1.0 - 1.0i * wx3) * my_sign(real(wx3)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) + 1.0i * M_PI * power_of<2>(log(1.0 - 1.0i * wx3)) * my_sign(real(wx3)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3)
                + dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + wx3)) * (4.0 * log(wx3 / (1.0i + wx3)) - 2.0 * log((wx3 - xd) / (1.0i + wx3)) - 2.0 * log((wx3 + xd) / (1.0i + wx3)) + (4.0 * 1.0i) * M_PI * my_sign(real(wx3)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3)
                    + (4.0 * 1.0i) * M_PI * my_sign(real(wx3)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3)) + (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 - (1.0i / 2.0) * wx3) * log(1.0 + 1.0i * wx3) * my_sign(-real(wx3)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                + 1.0i * M_PI * power_of<2>(log(1.0 + 1.0i * wx3)) * my_sign(-real(wx3)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) + dilog(1.0 / 2.0 + (1.0i / 2.0) * wx3) * (4.0 * log(wx3 / (-1.0i + wx3)) - 2.0 * log((wx3 - xd) / (-1.0i + wx3)) - 2.0 * log((wx3 + xd) / (-1.0i + wx3))
                    + (4.0 * 1.0i) * M_PI * my_sign(-real(wx3)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) + (4.0 * 1.0i) * M_PI * my_sign(-real(wx3)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)) + (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 + (1.0i / 2.0) * wx4) * log(1.0 - 1.0i * wx4) * my_sign(real(wx4)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4)
                + 1.0i * M_PI * power_of<2>(log(1.0 - 1.0i * wx4)) * my_sign(real(wx4)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) + dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + wx4)) * (4.0 * log(wx4 / (1.0i + wx4)) - 2.0 * log((wx4 - xd) / (1.0i + wx4)) - 2.0 * log((wx4 + xd) / (1.0i + wx4))
                    + (4.0 * 1.0i) * M_PI * my_sign(real(wx4)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) + (4.0 * 1.0i) * M_PI * my_sign(real(wx4)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4)) + (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 - (1.0i / 2.0) * wx4) * log(1.0 + 1.0i * wx4) * my_sign(-real(wx4)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4)
                + 1.0i * M_PI * power_of<2>(log(1.0 + 1.0i * wx4)) * my_sign(-real(wx4)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4) + log(1.0 - 1.0i * xd) * (log(1.0 + 1.0i * xd) * (log((wx3 + xd) / (-1.0i + wx3)) + log((wx3 + xd) / (1.0i + wx3)) + log((wx4 + xd) / (-1.0i + wx4))
                    + log((wx4 + xd) / (1.0i + wx4))) - (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 - (1.0i / 2.0) * wx3) * my_sign(2.0 * real(1.0 / (1.0i + xd))) * T(1.0, 1.0 / 2.0 + (1.0i / 2.0) * xd, 1.0 / 2.0 - (1.0i / 2.0) * wx3) - (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 + (1.0i / 2.0) * wx3) * my_sign(2.0 * real(1.0 / (1.0i + xd))) *
                    T(1.0, 1.0 / 2.0 + (1.0i / 2.0) * xd, 1.0 / 2.0 + (1.0i / 2.0) * wx3) - (4.0 * 1.0i) * M_PI * log(((-1.0 / 2.0) * 1.0i) * (1.0i + wx4)) * my_sign(2.0 * real(1.0 / (1.0i + xd))) * T(1.0, 1.0 / 2.0 + (1.0i / 2.0) * xd, 1.0 / 2.0 - (1.0i / 2.0) * wx4)
                    - (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 + (1.0i / 2.0) * wx4) * my_sign(2.0 * real(1.0 / (1.0i + xd))) * T(1.0, 1.0 / 2.0 + (1.0i / 2.0) * xd, 1.0 / 2.0 + (1.0i / 2.0) * wx4) + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(xdinv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3)
                    + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(-real(xdinv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(xdinv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4)
                    + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(-real(xdinv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4) + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3)
                    + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4)
                    + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4)) + li2half * (-4.0 * log(wx3 / (-1.0i + wx3)) - 4.0 * log(wx3 / (1.0i + wx3)) - 4.0 * log(wx4 / (-1.0i + wx4)) - 4.0 * log(wx4 / (1.0i + wx4))
                    + 2.0 * log((wx3 - xd) / (-1.0i + wx3)) + 2.0 * log((wx3 - xd) / (1.0i + wx3)) + 2.0 * log((wx4 - xd) / (-1.0i + wx4)) + 2.0 * log((wx4 - xd) / (1.0i + wx4)) + 2.0 * log((wx3 + xd) / (-1.0i + wx3)) + 2.0 * log((wx3 + xd) / (1.0i + wx3))
                    + 2.0 * log((wx4 + xd) / (-1.0i + wx4)) + 2.0 * log((wx4 + xd) / (1.0i + wx4)) - (4.0 * 1.0i) * M_PI * my_sign(real(wx3)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) - (4.0 * 1.0i) * M_PI * my_sign(-real(wx3)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                    - (4.0 * 1.0i) * M_PI * my_sign(real(wx4)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) - (4.0 * 1.0i) * M_PI * my_sign(-real(wx4)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4) - (4.0 * 1.0i) * M_PI * my_sign(real(wx3)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3)
                    - (4.0 * 1.0i) * M_PI * my_sign(-real(wx3)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) - (4.0 * 1.0i) * M_PI * my_sign(real(wx4)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) - (4.0 * 1.0i) * M_PI * my_sign(-real(wx4)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))
                + dilog(1.0 / 2.0 + (1.0i / 2.0) * wx4) * (4.0 * log(wx4 / (-1.0i + wx4)) - 2.0 * log((wx4 - xd) / (-1.0i + wx4)) - 2.0 * log((wx4 + xd) / (-1.0i + wx4)) + (4.0 * 1.0i) * M_PI * my_sign(-real(wx4)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    + (4.0 * 1.0i) * M_PI * my_sign(-real(wx4)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4)) + log(1.0 + 1.0i * xd) * ((-4.0 * 1.0i) * M_PI * log(1.0 / 2.0 - (1.0i / 2.0) * wx3) * my_sign(2.0 * real(1.0 / (1.0i - xd))) * T(1.0, 1.0 / 2.0 - (1.0i / 2.0) * xd, 1.0 / 2.0 - (1.0i / 2.0) * wx3)
                    - (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 + (1.0i / 2.0) * wx3) * my_sign(2.0 * real(1.0 / (1.0i - xd))) * T(1.0, 1.0 / 2.0 - (1.0i / 2.0) * xd, 1.0 / 2.0 + (1.0i / 2.0) * wx3) - (4.0 * 1.0i) * M_PI * log(((-1.0 / 2.0) * 1.0i) * (1.0i + wx4)) * my_sign(2.0 * real(1.0 / (1.0i - xd))) *
                    T(1.0, 1.0 / 2.0 - (1.0i / 2.0) * xd, 1.0 / 2.0 - (1.0i / 2.0) * wx4) - (4.0 * 1.0i) * M_PI * log(1.0 / 2.0 + (1.0i / 2.0) * wx4) * my_sign(2.0 * real(1.0 / (1.0i - xd))) * T(1.0, 1.0 / 2.0 - (1.0i / 2.0) * xd, 1.0 / 2.0 + (1.0i / 2.0) * wx4)
                    + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                    + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(real(xdinv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(xdinv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                    + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(real(xdinv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(xdinv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))
                - (2.0 * 1.0i) * M_PI * power_of<2>(log((wx3 - xd) / (-1.0i + wx3))) * my_sign(imag((-1.0i + xd) / (-1.0i + wx3))) * T(1.0, ((-1.0 / 2.0) * 1.0i) * (1.0i + xd), (wx3 - xd) / (-1.0i + wx3)) - (2.0 * 1.0i) * M_PI * power_of<2>(log((wx4 - xd) / (-1.0i + wx4))) * my_sign(imag((-1.0i + xd) / (-1.0i + wx4))) *
                T(1.0, ((-1.0 / 2.0) * 1.0i) * (1.0i + xd), (wx4 - xd) / (-1.0i + wx4)) + 1.0i * M_PI * power_of<2>(log((-1.0i + wx3) / wx3)) * my_sign(real(wx3inv)) * T(1.0, (wx3 + xd) / wx3, (-1.0i + wx3) / wx3)
                - (4.0 * 1.0i) * M_PI * ln2 * log((-1.0i + wx3) / wx3) * my_sign(-real(wx3)) * T(1.0, (wx3 + xd) / wx3, (-1.0i + wx3) / wx3) + 1.0i * M_PI * power_of<2>(log((1.0i + wx3) / wx3)) * my_sign(-real(wx3inv)) * T(1.0, (wx3 + xd) / wx3, (1.0i + wx3) / wx3)
                - (4.0 * 1.0i) * M_PI * ln2 * log((1.0i + wx3) / wx3) * my_sign(real(wx3)) * T(1.0, (wx3 + xd) / wx3, (1.0i + wx3) / wx3) + 1.0i * M_PI * power_of<2>(log((-1.0i + wx4) / wx4)) * my_sign(real(wx4inv)) * T(1.0, (wx4 + xd) / wx4, (-1.0i + wx4) / wx4)
                - (4.0 * 1.0i) * M_PI * ln2 * log((-1.0i + wx4) / wx4) * my_sign(-real(wx4)) * T(1.0, (wx4 + xd) / wx4, (-1.0i + wx4) / wx4) + 1.0i * M_PI * power_of<2>(log((1.0i + wx4) / wx4)) * my_sign(-real(wx4inv)) * T(1.0, (wx4 + xd) / wx4, (1.0i + wx4) / wx4)
                - (4.0 * 1.0i) * M_PI * ln2 * log((1.0i + wx4) / wx4) * my_sign(real(wx4)) * T(1.0, (wx4 + xd) / wx4, (1.0i + wx4) / wx4) + 1.0i * M_PI * power_of<2>(log((-1.0i + wx3) / wx3)) * my_sign(real(wx3inv)) * T(1.0, 1.0 - xd / wx3, (-1.0i + wx3) / wx3)
                - (4.0 * 1.0i) * M_PI * ln2 * log((-1.0i + wx3) / wx3) * my_sign(-real(wx3)) * T(1.0, 1.0 - xd / wx3, (-1.0i + wx3) / wx3) + 1.0i * M_PI * power_of<2>(log((1.0i + wx3) / wx3)) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - xd / wx3, (1.0i + wx3) / wx3)
                - (4.0 * 1.0i) * M_PI * ln2 * log((1.0i + wx3) / wx3) * my_sign(real(wx3)) * T(1.0, 1.0 - xd / wx3, (1.0i + wx3) / wx3) + 1.0i * M_PI * power_of<2>(log((-1.0i + wx4) / wx4)) * my_sign(real(wx4inv)) * T(1.0, 1.0 - xd / wx4, (-1.0i + wx4) / wx4)
                - (4.0 * 1.0i) * M_PI * ln2 * log((-1.0i + wx4) / wx4) * my_sign(-real(wx4)) * T(1.0, 1.0 - xd / wx4, (-1.0i + wx4) / wx4) + 1.0i * M_PI * power_of<2>(log((1.0i + wx4) / wx4)) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - xd / wx4, (1.0i + wx4) / wx4)
                - (4.0 * 1.0i) * M_PI * ln2 * log((1.0i + wx4) / wx4) * my_sign(real(wx4)) * T(1.0, 1.0 - xd / wx4, (1.0i + wx4) / wx4) + log(xd) * (-(log((1.0i - xd) / (1.0i + wx4)) * log((wx4 + xd) / wx4)) - log((1.0i + xd) / (1.0i - wx4)) * log((wx4 + xd) / wx4)
                    - log((-1.0i + xd) / (-1.0i + wx3)) * log(1.0 - xd / wx3) - log((1.0i + xd) / (1.0i + wx3)) * log(1.0 - xd / wx3) - log((-1.0i + xd) / (-1.0i + wx4)) * log(1.0 - xd / wx4) - log((1.0i + xd) / (1.0i + wx4)) * log(1.0 - xd / wx4)
                    - (2.0 * 1.0i) * M_PI * log((-1.0i + wx3) / wx3) * my_sign(-imag(wx3 * xdinv)) * T(1.0, (wx3 + xd) / wx3, (-1.0i + wx3) / wx3) - (2.0 * 1.0i) * M_PI * log((1.0i + wx3) / wx3) * my_sign(-imag(wx3 * xdinv)) * T(1.0, (wx3 + xd) / wx3, (1.0i + wx3) / wx3)
                    - (2.0 * 1.0i) * M_PI * log((-1.0i + wx4) / wx4) * my_sign(-imag(wx4 * xdinv)) * T(1.0, (wx4 + xd) / wx4, (-1.0i + wx4) / wx4) - (2.0 * 1.0i) * M_PI * log((1.0i + wx4) / wx4) * my_sign(-imag(wx4 * xdinv)) * T(1.0, (wx4 + xd) / wx4, (1.0i + wx4) / wx4)
                    - (2.0 * 1.0i) * M_PI * log((-1.0i + wx3) / wx3) * my_sign(imag(wx3 * xdinv)) * T(1.0, 1.0 - xd / wx3, (-1.0i + wx3) / wx3) - (2.0 * 1.0i) * M_PI * log((1.0i + wx3) / wx3) * my_sign(imag(wx3 * xdinv)) * T(1.0, 1.0 - xd / wx3, (1.0i + wx3) / wx3)
                    - (2.0 * 1.0i) * M_PI * log((-1.0i + wx4) / wx4) * my_sign(imag(wx4 * xdinv)) * T(1.0, 1.0 - xd / wx4, (-1.0i + wx4) / wx4) - (2.0 * 1.0i) * M_PI * log((1.0i + wx4) / wx4) * my_sign(imag(wx4 * xdinv)) * T(1.0, 1.0 - xd / wx4, (1.0i + wx4) / wx4))
                + power_of<2>(log(1.0 - 1.0i * xd)) * (2.0 * log(-1.0i / xd) - log(((-1.0i) * (wx3 - xd)) / ((1.0i + wx3) * xd)) / 2.0 - log(((-1.0i) * (wx4 - xd)) / ((1.0i + wx4) * xd)) / 2.0 + log((wx3 + xd) / (-1.0i + wx3)) / 2.0 - log(((-1.0i) * (wx3 + xd)) / ((-1.0i + wx3) * xd)) / 2.0
                    + log((wx4 + xd) / (-1.0i + wx4)) / 2.0 - log(((-1.0i) * (wx4 + xd)) / ((-1.0i + wx4) * xd)) / 2.0 - 1.0i * M_PI * my_sign(real(xd)) * T(1.0, (wx3 + xd) / wx3, 1.0 - 1.0i * xd) - 1.0i * M_PI * my_sign(real(xd)) * T(1.0, (wx4 + xd) / wx4, 1.0 - 1.0i * xd)
                    - 1.0i * M_PI * my_sign(real(xd)) * T(1.0, 1.0 - xd / wx3, 1.0 - 1.0i * xd) - 1.0i * M_PI * my_sign(real(xd)) * T(1.0, 1.0 - xd / wx4, 1.0 - 1.0i * xd)) + power_of<2>(log(1.0 + 1.0i * xd)) * (2.0 * log(1.0i / xd) - log((1.0i * (wx3 - xd)) / ((-1.0i + wx3) * xd)) / 2.0
                    - log((1.0i * (wx4 - xd)) / ((-1.0i + wx4) * xd)) / 2.0 + log((wx3 + xd) / (1.0i + wx3)) / 2.0 - log((1.0i * (wx3 + xd)) / ((1.0i + wx3) * xd)) / 2.0 + log((wx4 + xd) / (1.0i + wx4)) / 2.0 - log((1.0i * (wx4 + xd)) / ((1.0i + wx4) * xd)) / 2.0
                    - 1.0i * M_PI * my_sign(-real(xd)) * T(1.0, (wx3 + xd) / wx3, 1.0 + 1.0i * xd) - 1.0i * M_PI * my_sign(-real(xd)) * T(1.0, (wx4 + xd) / wx4, 1.0 + 1.0i * xd) - 1.0i * M_PI * my_sign(-real(xd)) * T(1.0, 1.0 - xd / wx3, 1.0 + 1.0i * xd) - 1.0i * M_PI * my_sign(-real(xd)) * T(1.0, 1.0 - xd / wx4, 1.0 + 1.0i * xd)));

            const complex<double> f29dPart42 = num16 * ((-2.0 * li3half + trilog(1.0 / (1.0 - w4)) + trilog(-w4inv) - trilog((-1.0 + w4) / (2.0 * w4)) + trilog((-1.0 + w4) / w4) - trilog(-w4) + trilog(w4) - trilog((1.0 - w4) / (1.0 + w4))
                + trilog(1.0 + w5inv) + trilog(w5inv) + trilog(-w5) - trilog(w5) + trilog(1.0 / (1.0 + w5)) - trilog((1.0 + w5) / (1.0 - w5)) - trilog((1.0 + w5) / (2.0 * w5)) - 2.0 * trilog(1.0 + ydinv) + 2.0 * trilog((1.0 - yd) / 2.0)
                - trilog((1.0 - yd) / (1.0 + w5)) - trilog((-1.0 + w4) / (w4 - yd)) - trilog((-1.0 + yd) / (-1.0 + w4)) + 2.0 * trilog((-1.0 + yd) / yd) - 2.0 * trilog(-yd) + 2.0 * trilog(yd) + trilog(((-1.0 + w4) * (1.0 + yd)) / (2.0 * (w4 - yd)))
                - trilog((w4 * (1.0 + yd)) / (w4 - yd)) + trilog(((-1.0 + w4) * (1.0 + yd)) / ((1.0 + w4) * (-1.0 + yd))) + trilog(((1.0 + w5) * (1.0 + yd)) / ((-1.0 + w5) * (-1.0 + yd))) - trilog((-1.0 + yd) / (-w4 + yd)) - trilog((1.0 + w5) / (w5 + yd))
                - trilog((-1.0 + yd) / (w5 + yd)) - trilog((w5 * (1.0 + yd)) / (w5 + yd)) + trilog(((1.0 + w5) * (1.0 + yd)) / (2.0 * (w5 + yd))) + trilog((w4 - w4 * yd) / (w4 - yd)) + trilog((w5 - w5 * yd) / (w5 + yd))
                - dilog(-w4inv) * lnhalf + dilog((-1.0 + w4) / (2.0 * w4)) * lnhalf + dilog((1.0 - w4) / (1.0 + w4)) * lnhalf - dilog(w5inv) * lnhalf + dilog((1.0 + w5) / (1.0 - w5)) * lnhalf
                + dilog((1.0 + w5) / (2.0 * w5)) * lnhalf - (power_of<2>(lnhalf) * log(1.0 + w4inv)) / 2.0 + (pisqu * log(2.0 / w4)) / 6.0 + power_of<3>(log(2.0 / w4)) / 6.0 - (pisqu * log((2.0 * (-1.0 + w4)) / (1.0 + w4))) / 6.0 - power_of<3>(log((2.0 * (-1.0 + w4)) / (1.0 + w4))) / 6.0
                - (power_of<2>(lnhalf) * log(w4 / (1.0 + w4))) / 2.0 + (power_of<2>(lnhalf) * log((2.0 * w4) / (1.0 + w4))) / 2.0 + (power_of<2>(lnhalf) * log((1.0 + w4) / (2.0 * w4))) / 2.0 + (pisqu * log(-2.0 / w5)) / 6.0 + power_of<3>(log(-2.0 / w5)) / 6.0
                + (power_of<2>(lnhalf) * log((-1.0 + w5) / (2.0 * w5))) / 2.0 - (power_of<2>(lnhalf) * log((-1.0 + w5) / w5)) / 2.0 - (power_of<2>(lnhalf) * log(w5 / (-1.0 + w5))) / 2.0 + (power_of<2>(lnhalf) * log((2.0 * w5) / (-1.0 + w5))) / 2.0 - (pisqu * log((2.0 * (1.0 + w5)) / (-1.0 + w5))) / 6.0
                - power_of<3>(log((2.0 * (1.0 + w5)) / (-1.0 + w5))) / 6.0 + dilog(1.0 / (1.0 - w4)) * (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd))
                + dilog(1.0 / (1.0 + w5)) * (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd)) + dilog((1.0 - yd) / (1.0 + w5)) * (2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd))
                + dilog((-1.0 + yd) / (-1.0 + w4)) * (2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd)) + dilog(1.0 / (1.0 + w4)) * (-lnhalf + 2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd))
                + dilog(1.0 / (1.0 - w5)) * (-lnhalf + 2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd)) - dilog(2.0 / (1.0 + w4)) * ln1myd - dilog(-2.0 / (-1.0 + w5)) * ln1myd
                - 2.0 * dilog((-1.0 + yd) / yd) * ln1myd + dilog((-1.0 + yd) / (-w4 + yd)) * ln1myd + dilog((-1.0 + yd) / (w5 + yd)) * ln1myd - dilog((w4 - w4 * yd) / (w4 - yd)) * ln1myd
                - dilog((w5 - w5 * yd) / (w5 + yd)) * ln1myd + (pisqu * log(1.0 / (w4 - yd))) / 6.0 + power_of<3>(log(1.0 / (w4 - yd))) / 6.0 - (pisqu * log(2.0 / (w4 - yd))) / 6.0 - power_of<3>(log(2.0 / (w4 - yd))) / 6.0
                + (pisqu * log((-2.0 * (-1.0 + w4)) / ((1.0 + w4) * (-1.0 + yd)))) / 6.0 + power_of<3>(log((-2.0 * (-1.0 + w4)) / ((1.0 + w4) * (-1.0 + yd)))) / 6.0 + (pisqu * log((-2.0 * (1.0 + w5)) / ((-1.0 + w5) * (-1.0 + yd)))) / 6.0 + power_of<3>(log((-2.0 * (1.0 + w5)) / ((-1.0 + w5) * (-1.0 + yd)))) / 6.0
                - power_of<3>(log(-ydinv)) / 3.0 + (pisqu * log(ydinv)) / 3.0 + power_of<3>(log(ydinv)) / 3.0 - dilog(((-1.0 + w4) * (1.0 + yd)) / (2.0 * (w4 - yd))) * log((1.0 + yd) / 2.0) - dilog(((-1.0 + w4) * (1.0 + yd)) / ((1.0 + w4) * (-1.0 + yd))) * log((1.0 + yd) / 2.0)
                - dilog(((1.0 + w5) * (1.0 + yd)) / ((-1.0 + w5) * (-1.0 + yd))) * log((1.0 + yd) / 2.0) - dilog(((1.0 + w5) * (1.0 + yd)) / (2.0 * (w5 + yd))) * log((1.0 + yd) / 2.0) + dilog(-((1.0 + yd) / (w4 - yd))) * (log((1.0 + yd) / 2.0) - ln1pyd)
                + dilog((1.0 + yd) / (w5 + yd)) * (log((1.0 + yd) / 2.0) - ln1pyd) + dilog((1.0 + yd) / (1.0 + w4)) * (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd) + ln1myd + log((1.0 + yd) / 2.0)
                    - ln1pyd) + dilog((1.0 + yd) / (1.0 - w5)) * (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd) + ln1myd + log((1.0 + yd) / 2.0) - ln1pyd) + 2.0 * dilog(1.0 + ydinv) * ln1pyd
                + dilog((w4 * (1.0 + yd)) / (w4 - yd)) * ln1pyd + dilog((w5 * (1.0 + yd)) / (w5 + yd)) * ln1pyd - (log((1.0 + w4) / (w4 - yd)) * power_of<2>(ln1pyd)) / 2.0
                + log((w4 - yd) / (1.0 + w4)) * (power_of<2>(log((1.0 + yd) / 2.0)) / 2.0 + (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd)) * ln1pyd - power_of<2>(ln1pyd) / 2.0) + log(-ydinv) * ((-1.0 / 3.0) * pisqu + power_of<2>(ln1pyd))
                - (pisqu * log(1.0 / (-w4 + yd))) / 6.0 - power_of<3>(log(1.0 / (-w4 + yd))) / 6.0 - (pisqu * log(-2.0 / (w5 + yd))) / 6.0 - power_of<3>(log(-2.0 / (w5 + yd))) / 6.0 + (pisqu * log(-1.0 / (w5 + yd))) / 6.0 + power_of<3>(log(-1.0 / (w5 + yd))) / 6.0
                - (pisqu * log(1.0 / (w5 + yd))) / 6.0 - power_of<3>(log(1.0 / (w5 + yd))) / 6.0 + power_of<2>(ln1myd) * (log((-1.0 + w4) / (w4 - yd)) / 2.0 - log(ydinv) - log(((-1.0 + w4) * yd) / (w4 - yd)) / 2.0 + log((1.0 + w5) / (w5 + yd)) / 2.0
                    - log(((1.0 + w5) * yd) / (w5 + yd)) / 2.0) + (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd)) * ln1pyd * log((w5 + yd) / (-1.0 + w5)) + dilog(-w5) * log((w5 + yd) / w5)
                - dilog(w5) * log((w5 + yd) / w5) + dilog((1.0i - xd) / (1.0i + wx3)) * log((w5 + yd) / w5) + dilog((1.0i - xd) / (1.0i + wx4)) * log((w5 + yd) / w5) - dilog(-(xd / wx3)) * log((w5 + yd) / w5)
                - dilog(-(xd / wx4)) * log((w5 + yd) / w5) + dilog((1.0i + xd) / (1.0i - wx3)) * log((w5 + yd) / w5) + dilog((1.0i + xd) / (1.0i - wx4)) * log((w5 + yd) / w5)
                + power_of<2>(log((1.0 + yd) / 2.0)) * (log((1.0 + w4) / (w4 - yd)) / 2.0 - log((-2.0 * (w4 - yd)) / ((1.0 + w4) * (-1.0 + yd))) / 2.0 - log((-1.0 / 2.0) * (((1.0 + w4) * (-1.0 + yd)) / (w4 - yd))) / 2.0 + log((-1.0 + w5) / (w5 + yd)) / 2.0
                    - log((-1.0 / 2.0) * (((-1.0 + w5) * (-1.0 + yd)) / (w5 + yd))) / 2.0 + log((w5 + yd) / (-1.0 + w5)) / 2.0 - log((-2.0 * (w5 + yd)) / ((-1.0 + w5) * (-1.0 + yd))) / 2.0) + dilog(-1.0i / (-1.0i + wx3)) * (-log((w5 + yd) / w5) - log(1.0 - yd / w4))
                + dilog(1.0i / (1.0i + wx3)) * (-log((w5 + yd) / w5) - log(1.0 - yd / w4)) + dilog(-1.0i / (-1.0i + wx4)) * (-log((w5 + yd) / w5) - log(1.0 - yd / w4)) + dilog(1.0i / (1.0i + wx4)) * (-log((w5 + yd) / w5) - log(1.0 - yd / w4))
                - dilog(-w4) * log(1.0 - yd / w4) + dilog(w4) * log(1.0 - yd / w4) - dilog(xd / wx3) * log(1.0 - yd / w4) - dilog(xd / wx4) * log(1.0 - yd / w4) + dilog((-1.0i + xd) / (-1.0i + wx3)) * log(1.0 - yd / w4)
                + dilog((-1.0i + xd) / (-1.0i + wx4)) * log(1.0 - yd / w4) + dilog((1.0i + xd) / (1.0i + wx3)) * log(1.0 - yd / w4) + dilog((1.0i + xd) / (1.0i + wx4)) * log(1.0 - yd / w4)
                + dilog(1.0 / 2.0 + (1.0i / 2.0) * xd) * (log((w5 + yd) / w5) + log(1.0 - yd / w4)) + dilog((-1.0i) * xd) * (log((w5 + yd) / w5) + log(1.0 - yd / w4)) + dilog(1.0i * xd) * (log((w5 + yd) / w5) + log(1.0 - yd / w4))
                + dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + xd)) * (log((w5 + yd) / w5) + log(1.0 - yd / w4)) + power_of<2>(ln1pyd) * (log(((1.0 + w4) * yd) / (-w4 + yd)) / 2.0 - log((-1.0 + w5) / (w5 + yd)) / 2.0 - log((w5 + yd) / (-1.0 + w5)) / 2.0
                    + log((yd - w5 * yd) / (w5 + yd)) / 2.0) - 1.0i * M_PI * H1((-1.0 + w4) / (-1.0 + yd), -2.0 / (-1.0 + yd)) * power_of<2>(log(((1.0 + w4) * (-1.0 + yd)) / (2.0 * (-1.0 + w4)))) * my_sign(2.0 * imag(1.0 / (1.0 - yd)))
                - 1.0i * M_PI * H1((1.0 + w5) / (1.0 - yd), -2.0 / (-1.0 + yd)) * power_of<2>(log(((-1.0 + w5) * (-1.0 + yd)) / (2.0 * (1.0 + w5)))) * my_sign(2.0 * imag(1.0 / (1.0 - yd))) + 1.0i * M_PI * H1(-2.0 / (-1.0 + w4), -2.0 / (-1.0 + yd)) * power_of<2>(log((-w4 + yd) / 2.0)) * my_sign(2.0 * imag(1.0 / (1.0 - yd)))
                + 1.0i * M_PI * H1(2.0 / (1.0 + w5), -2.0 / (-1.0 + yd)) * power_of<2>(log((w5 + yd) / 2.0)) * my_sign(2.0 * imag(1.0 / (1.0 - yd))) - 1.0i * M_PI * H1(-w4inv, -ydinv) * power_of<2>(log(-w4 + yd)) * my_sign(-imydinv)
                - 1.0i * M_PI * H1(w5inv, -ydinv) * power_of<2>(log(w5 + yd)) * my_sign(-imydinv) + 1.0i * M_PI * H1(w4inv, ydinv) * power_of<2>(log(w4 - yd)) * my_sign(imydinv) + 1.0i * M_PI * H1(-w5inv, ydinv) * power_of<2>(log(-w5 - yd)) * my_sign(imydinv)
                + log((w5 + yd) / w5) * ((-1.0 / 2.0) * power_of<2>(log(1.0 - 1.0i * xd)) - power_of<2>(log(1.0 + 1.0i * xd)) / 2.0 - 2.0 * ln2 * log((wx3 + xd) / wx3) + log(xd) * (-log((wx3 + xd) / wx3) - log((wx4 + xd) / wx4)) - 2.0 * ln2 * log((wx4 + xd) / wx4)
                    + log(1.0 - 1.0i * xd) * (ln2 + log(xd) + log((wx3 + xd) / (-1.0i + wx3)) + log((wx4 + xd) / (-1.0i + wx4))) + log(1.0 + 1.0i * xd) * (ln2 + log(xd) + log((wx3 + xd) / (1.0i + wx3)) + log((wx4 + xd) / (1.0i + wx4)))
                    + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4))
                + log(1.0 - yd / w4) * ((-1.0 / 2.0) * power_of<2>(log(1.0 - 1.0i * xd)) - power_of<2>(log(1.0 + 1.0i * xd)) / 2.0 + log(1.0 + 1.0i * xd) * (ln2 + log((wx3 - xd) / (-1.0i + wx3)) + log((wx4 - xd) / (-1.0i + wx4)) + log(xd))
                    + log(1.0 - 1.0i * xd) * (ln2 + log((wx3 - xd) / (1.0i + wx3)) + log((wx4 - xd) / (1.0i + wx4)) + log(xd)) - 2.0 * ln2 * log(1.0 - xd / wx3) + log(xd) * (-log(1.0 - xd / wx3) - log(1.0 - xd / wx4))
                    - 2.0 * ln2 * log(1.0 - xd / wx4) + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4)
                    + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))
                + (4.0 * 1.0i) * M_PI * ln2 * log(1.0 - w4) * my_sign(imag(w4inv)) * T(1.0, 1.0 - yd, 1.0 - w4) - (2.0 * 1.0i) * M_PI * log(1.0 - w4) * log((1.0 + w4) / 2.0) * my_sign(imag(w4)) * T(1.0, 1.0 - yd, 1.0 - w4)
                + dilog((1.0 - w4) / 2.0) * (-log(w4 / (-1.0 + w4)) + log((w4 - yd) / (-1.0 + w4)) - (2.0 * 1.0i) * M_PI * my_sign(imag(w4)) * T(1.0, 1.0 - yd, 1.0 - w4)) + (4.0 * 1.0i) * M_PI * ln2 * log(1.0 + w5) * my_sign(-imag(w5inv)) * T(1.0, 1.0 - yd, 1.0 + w5)
                - (2.0 * 1.0i) * M_PI * log((1.0 - w5) / 2.0) * log(1.0 + w5) * my_sign(-imag(w5)) * T(1.0, 1.0 - yd, 1.0 + w5) + dilog((1.0 + w5) / 2.0) * (-log(w5 / (1.0 + w5)) + log((w5 + yd) / (1.0 + w5)) - (2.0 * 1.0i) * M_PI * my_sign(-imag(w5)) * T(1.0, 1.0 - yd, 1.0 + w5))
                + li2half * (log(w4 / (-1.0 + w4)) + log(w5 / (1.0 + w5)) - log((w4 - yd) / (-1.0 + w4)) - 2.0 * log((w5 + yd) / w5) - log((w5 + yd) / (1.0 + w5)) - 2.0 * log(1.0 - yd / w4) + (2.0 * 1.0i) * M_PI * my_sign(imag(w4)) * T(1.0, 1.0 - yd, 1.0 - w4)
                    + (2.0 * 1.0i) * M_PI * my_sign(-imag(w5)) * T(1.0, 1.0 - yd, 1.0 + w5)) - 1.0i * M_PI * power_of<2>(log((1.0 + w4) / 2.0)) * my_sign((-1.0 / 2.0) * imag(w4)) * T(1.0, (1.0 + yd) / 2.0, (1.0 + w4) / 2.0) - 1.0i * M_PI * power_of<2>(log((1.0 - w5) / 2.0)) * my_sign(imag(w5) / 2.0) * T(1.0, (1.0 + yd) / 2.0, (1.0 - w5) / 2.0)
                + ln1myd * ((2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd)) * log((w4 - yd) / (-1.0 + w4)) + log((w4 - yd) / (1.0 + w4)) * log((1.0 + yd) / 2.0) + log((1.0 + yd) / 2.0) * log((w5 + yd) / (-1.0 + w5))
                    + (2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd)) * log((w5 + yd) / (1.0 + w5)) + (2.0 * 1.0i) * M_PI * log((1.0 + w4) / 2.0) * my_sign(2.0 * imag(1.0 / (1.0 - yd))) * T(1.0, (1.0 + yd) / 2.0, (1.0 + w4) / 2.0)
                    + (2.0 * 1.0i) * M_PI * log((1.0 - w5) / 2.0) * my_sign(2.0 * imag(1.0 / (1.0 - yd))) * T(1.0, (1.0 + yd) / 2.0, (1.0 - w5) / 2.0)) + 1.0i * M_PI * power_of<2>(log((w4 - yd) / (-1.0 + w4))) * my_sign(imag((-1.0 + yd) / (-1.0 + w4))) * T(1.0, (1.0 + yd) / 2.0, (w4 - yd) / (-1.0 + w4))
                + 1.0i * M_PI * power_of<2>(log((w5 + yd) / (1.0 + w5))) * my_sign(imag((1.0 - yd) / (1.0 + w5))) * T(1.0, (1.0 + yd) / 2.0, (w5 + yd) / (1.0 + w5)) - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 + w4) * my_sign(-imag(w4inv)) * T(1.0, 1.0 + yd, 1.0 + w4)
                + 1.0i * M_PI * power_of<2>(log(1.0 + w4)) * my_sign(-imag(w4)) * T(1.0, 1.0 + yd, 1.0 + w4) - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 - w5) * my_sign(imag(w5inv)) * T(1.0, 1.0 + yd, 1.0 - w5) + 1.0i * M_PI * power_of<2>(log(1.0 - w5)) * my_sign(imag(w5)) * T(1.0, 1.0 + yd, 1.0 - w5)
                + log(xd) * ((2.0 * 1.0i) * M_PI * log(1.0 - w4) * my_sign(imag(w4inv)) * T(1.0, 1.0 - yd, 1.0 - w4) + (2.0 * 1.0i) * M_PI * log(1.0 + w5) * my_sign(-imag(w5inv)) * T(1.0, 1.0 - yd, 1.0 + w5)
                    - (2.0 * 1.0i) * M_PI * log(1.0 + w4) * my_sign(-imag(w4inv)) * T(1.0, 1.0 + yd, 1.0 + w4) - (2.0 * 1.0i) * M_PI * log(1.0 - w5) * my_sign(imag(w5inv)) * T(1.0, 1.0 + yd, 1.0 - w5))
                + log(1.0 - 1.0i * xd) * ((-2.0 * 1.0i) * M_PI * log(1.0 - w4) * my_sign(imag(w4inv)) * T(1.0, 1.0 - yd, 1.0 - w4) - (2.0 * 1.0i) * M_PI * log(1.0 + w5) * my_sign(-imag(w5inv)) * T(1.0, 1.0 - yd, 1.0 + w5)
                    + (2.0 * 1.0i) * M_PI * log(1.0 + w4) * my_sign(-imag(w4inv)) * T(1.0, 1.0 + yd, 1.0 + w4) + (2.0 * 1.0i) * M_PI * log(1.0 - w5) * my_sign(imag(w5inv)) * T(1.0, 1.0 + yd, 1.0 - w5))
                + log(1.0 + 1.0i * xd) * ((-2.0 * 1.0i) * M_PI * log(1.0 - w4) * my_sign(imag(w4inv)) * T(1.0, 1.0 - yd, 1.0 - w4) - (2.0 * 1.0i) * M_PI * log(1.0 + w5) * my_sign(-imag(w5inv)) * T(1.0, 1.0 - yd, 1.0 + w5)
                    + (2.0 * 1.0i) * M_PI * log(1.0 + w4) * my_sign(-imag(w4inv)) * T(1.0, 1.0 + yd, 1.0 + w4) + (2.0 * 1.0i) * M_PI * log(1.0 - w5) * my_sign(imag(w5inv)) * T(1.0, 1.0 + yd, 1.0 - w5))
                + power_of<2>(log((w5 + yd) / w5)) * (1.0i * M_PI * my_sign(-imag(yd / w5)) * T(1.0, 1.0 - yd, (w5 + yd) / w5) - 1.0i * M_PI * my_sign(-imag(yd / w5)) * T(1.0, 1.0 + yd, (w5 + yd) / w5))
                + power_of<2>(log(1.0 - yd / w4)) * (1.0i * M_PI * my_sign(imag(yd / w4)) * T(1.0, 1.0 - yd, 1.0 - yd / w4) - 1.0i * M_PI * my_sign(imag(yd / w4)) * T(1.0, 1.0 + yd, 1.0 - yd / w4))) / denom2);

            const complex<double> f29dPart43 = num17 * ((-2.0 * li3half + trilog(1.0 / (1.0 - w4)) + trilog(-w4inv) - trilog((-1.0 + w4) / (2.0 * w4)) + trilog((-1.0 + w4) / w4) - trilog(-w4) + trilog(w4) - trilog((1.0 - w4) / (1.0 + w4))
                + trilog(1.0 + w5inv) + trilog(w5inv) + trilog(-w5) - trilog(w5) + trilog(1.0 / (1.0 + w5)) - trilog((1.0 + w5) / (1.0 - w5)) - trilog((1.0 + w5) / (2.0 * w5)) + 2.0 * trilog(1.0 + ydinv) - trilog((1.0 + w5) / (w5 - yd))
                + trilog((-1.0 / 2.0) * (((1.0 + w5) * (-1.0 + yd)) / (w5 - yd))) - 2.0 * trilog((-1.0 + yd) / yd) + 2.0 * trilog(-yd) - 2.0 * trilog(yd) + trilog(((-1.0 + w4) * (-1.0 + yd)) / ((1.0 + w4) * (1.0 + yd)))
                + trilog(((1.0 + w5) * (-1.0 + yd)) / ((-1.0 + w5) * (1.0 + yd))) + 2.0 * trilog((1.0 + yd) / 2.0) - trilog((1.0 + yd) / (1.0 - w4)) - trilog((1.0 + yd) / (1.0 + w5)) - trilog(-((1.0 + yd) / (w5 - yd))) + trilog((w5 * (1.0 + yd)) / (w5 - yd))
                - trilog((-1.0 + w4) / (w4 + yd)) + trilog((-1.0 / 2.0) * (((-1.0 + w4) * (-1.0 + yd)) / (w4 + yd))) - trilog((1.0 + yd) / (w4 + yd)) + trilog((w4 * (1.0 + yd)) / (w4 + yd)) - trilog((w4 - w4 * yd) / (w4 + yd))
                - trilog((w5 - w5 * yd) / (w5 - yd)) - dilog(-w4inv) * lnhalf + dilog((-1.0 + w4) / (2.0 * w4)) * lnhalf + dilog((1.0 - w4) / (1.0 + w4)) * lnhalf - dilog(w5inv) * lnhalf
                + dilog((1.0 + w5) / (1.0 - w5)) * lnhalf + dilog((1.0 + w5) / (2.0 * w5)) * lnhalf - (power_of<2>(lnhalf) * log(1.0 + w4inv)) / 2.0 + (pisqu * log(2.0 / w4)) / 6.0 + power_of<3>(log(2.0 / w4)) / 6.0 - (pisqu * log((2.0 * (-1.0 + w4)) / (1.0 + w4))) / 6.0
                - power_of<3>(log((2.0 * (-1.0 + w4)) / (1.0 + w4))) / 6.0 - (power_of<2>(lnhalf) * log(w4 / (1.0 + w4))) / 2.0 + (power_of<2>(lnhalf) * log((2.0 * w4) / (1.0 + w4))) / 2.0 + (power_of<2>(lnhalf) * log((1.0 + w4) / (2.0 * w4))) / 2.0 + (pisqu * log(-2.0 / w5)) / 6.0 + power_of<3>(log(-2.0 / w5)) / 6.0
                + (power_of<2>(lnhalf) * log((-1.0 + w5) / (2.0 * w5))) / 2.0 - (power_of<2>(lnhalf) * log((-1.0 + w5) / w5)) / 2.0 - (power_of<2>(lnhalf) * log(w5 / (-1.0 + w5))) / 2.0 + (power_of<2>(lnhalf) * log((2.0 * w5) / (-1.0 + w5))) / 2.0 - (pisqu * log((2.0 * (1.0 + w5)) / (-1.0 + w5))) / 6.0
                - power_of<3>(log((2.0 * (1.0 + w5)) / (-1.0 + w5))) / 6.0 + dilog(1.0 / (1.0 - w4)) * (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd))
                + dilog(1.0 / (1.0 + w5)) * (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd)) + dilog((1.0 + yd) / (1.0 - w4)) * (2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd))
                + dilog((1.0 + yd) / (1.0 + w5)) * (2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd)) + dilog(1.0 / (1.0 + w4)) * (-lnhalf + 2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd))
                + dilog(1.0 / (1.0 - w5)) * (-lnhalf + 2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd)) - dilog((-1.0 / 2.0) * (((1.0 + w5) * (-1.0 + yd)) / (w5 - yd))) * log((1.0 - yd) / 2.0)
                - dilog(((-1.0 + w4) * (-1.0 + yd)) / ((1.0 + w4) * (1.0 + yd))) * log((1.0 - yd) / 2.0) - dilog(((1.0 + w5) * (-1.0 + yd)) / ((-1.0 + w5) * (1.0 + yd))) * log((1.0 - yd) / 2.0)
                - dilog((-1.0 / 2.0) * (((-1.0 + w4) * (-1.0 + yd)) / (w4 + yd))) * log((1.0 - yd) / 2.0) + dilog((-1.0 + yd) / (w4 + yd)) * (log((1.0 - yd) / 2.0) - ln1myd)
                + dilog((-1.0 + yd) / (-w5 + yd)) * (log((1.0 - yd) / 2.0) - ln1myd) + 2.0 * dilog((-1.0 + yd) / yd) * ln1myd + dilog((w4 - w4 * yd) / (w4 + yd)) * ln1myd + dilog((w5 - w5 * yd) / (w5 - yd)) * ln1myd
                - (pisqu * log(-2.0 / (w5 - yd))) / 6.0 - power_of<3>(log(-2.0 / (w5 - yd))) / 6.0 - (pisqu * log(1.0 / (w5 - yd))) / 6.0 - power_of<3>(log(1.0 / (w5 - yd))) / 6.0 - (power_of<2>(ln1myd) * log((-1.0 + w5) / (w5 - yd))) / 2.0 + power_of<3>(log(-ydinv)) / 3.0
                - (pisqu * log(ydinv)) / 3.0 - power_of<3>(log(ydinv)) / 3.0 + (pisqu * log((2.0 * (-1.0 + w4)) / ((1.0 + w4) * (1.0 + yd)))) / 6.0 + power_of<3>(log((2.0 * (-1.0 + w4)) / ((1.0 + w4) * (1.0 + yd)))) / 6.0 + (pisqu * log((2.0 * (1.0 + w5)) / ((-1.0 + w5) * (1.0 + yd)))) / 6.0
                + power_of<3>(log((2.0 * (1.0 + w5)) / ((-1.0 + w5) * (1.0 + yd)))) / 6.0 - dilog(2.0 / (1.0 + w4)) * ln1pyd - dilog(-2.0 / (-1.0 + w5)) * ln1pyd - 2.0 * dilog(1.0 + ydinv) * ln1pyd + dilog(-((1.0 + yd) / (w5 - yd))) * ln1pyd
                - dilog((w5 * (1.0 + yd)) / (w5 - yd)) * ln1pyd + dilog((1.0 + yd) / (w4 + yd)) * ln1pyd - dilog((w4 * (1.0 + yd)) / (w4 + yd)) * ln1pyd
                + dilog((1.0 - yd) / (1.0 + w4)) * (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd) + log((1.0 - yd) / 2.0) - ln1myd + ln1pyd)
                + dilog((-1.0 + yd) / (-1.0 + w5)) * (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd) + log((1.0 - yd) / 2.0) - ln1myd + ln1pyd) + log(-ydinv) * (pisqu / 3.0 - power_of<2>(ln1pyd))
                - (pisqu * log(-1.0 / (w4 + yd))) / 6.0 - power_of<3>(log(-1.0 / (w4 + yd))) / 6.0 + (pisqu * log(1.0 / (w4 + yd))) / 6.0 + power_of<3>(log(1.0 / (w4 + yd))) / 6.0 - (pisqu * log(2.0 / (w4 + yd))) / 6.0 - power_of<3>(log(2.0 / (w4 + yd))) / 6.0
                - dilog(-w4) * log((w4 + yd) / w4) + dilog(w4) * log((w4 + yd) / w4) - dilog(xd / wx3) * log((w4 + yd) / w4) - dilog(xd / wx4) * log((w4 + yd) / w4) + dilog((-1.0i + xd) / (-1.0i + wx3)) * log((w4 + yd) / w4)
                + dilog((-1.0i + xd) / (-1.0i + wx4)) * log((w4 + yd) / w4) + dilog((1.0i + xd) / (1.0i + wx3)) * log((w4 + yd) / w4) + dilog((1.0i + xd) / (1.0i + wx4)) * log((w4 + yd) / w4)
                + power_of<2>(ln1myd) * ((-1.0 / 2.0) * log((w5 - yd) / (-1.0 + w5)) + log(ydinv) + log(((-1.0 + w5) * yd) / (w5 - yd)) / 2.0 - log((1.0 + w4) / (w4 + yd)) / 2.0 + log(((1.0 + w4) * yd) / (w4 + yd)) / 2.0 - log((w4 + yd) / (1.0 + w4)) / 2.0)
                + ln1myd * ((-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd)) * log((w5 - yd) / (-1.0 + w5)) + (-2.0 * ln2 + log(1.0 - 1.0i * xd) + log(1.0 + 1.0i * xd) - log(xd)) * log((w4 + yd) / (1.0 + w4)))
                + log((1.0 - yd) / 2.0) * (log((w5 - yd) / (-1.0 + w5)) * ln1pyd + ln1pyd * log((w4 + yd) / (1.0 + w4))) + power_of<2>(log((1.0 - yd) / 2.0)) * (log((-1.0 + w5) / (w5 - yd)) / 2.0 + log((w5 - yd) / (-1.0 + w5)) / 2.0
                    - log((2.0 * (w5 - yd)) / ((-1.0 + w5) * (1.0 + yd))) / 2.0 - log(((-1.0 + w5) * (1.0 + yd)) / (2.0 * (w5 - yd))) / 2.0 + log((1.0 + w4) / (w4 + yd)) / 2.0 - log(((1.0 + w4) * (1.0 + yd)) / (2.0 * (w4 + yd))) / 2.0 + log((w4 + yd) / (1.0 + w4)) / 2.0
                    - log((2.0 * (w4 + yd)) / ((1.0 + w4) * (1.0 + yd))) / 2.0) + (pisqu * log(1.0 / (-w5 + yd))) / 6.0 + power_of<3>(log(1.0 / (-w5 + yd))) / 6.0 + power_of<2>(ln1pyd) * (log((1.0 + w5) / (w5 - yd)) / 2.0 + log((-1.0 + w4) / (w4 + yd)) / 2.0
                    - log(((1.0 + w5) * yd) / (-w5 + yd)) / 2.0 - log((yd - w4 * yd) / (w4 + yd)) / 2.0) + dilog(-1.0i / (-1.0i + wx3)) * (-log((w4 + yd) / w4) - log(1.0 - yd / w5)) + dilog(1.0i / (1.0i + wx3)) * (-log((w4 + yd) / w4) - log(1.0 - yd / w5))
                + dilog(-1.0i / (-1.0i + wx4)) * (-log((w4 + yd) / w4) - log(1.0 - yd / w5)) + dilog(1.0i / (1.0i + wx4)) * (-log((w4 + yd) / w4) - log(1.0 - yd / w5)) + dilog(-w5) * log(1.0 - yd / w5) - dilog(w5) * log(1.0 - yd / w5)
                + dilog((1.0i - xd) / (1.0i + wx3)) * log(1.0 - yd / w5) + dilog((1.0i - xd) / (1.0i + wx4)) * log(1.0 - yd / w5) - dilog(-(xd / wx3)) * log(1.0 - yd / w5) - dilog(-(xd / wx4)) * log(1.0 - yd / w5)
                + dilog((1.0i + xd) / (1.0i - wx3)) * log(1.0 - yd / w5) + dilog((1.0i + xd) / (1.0i - wx4)) * log(1.0 - yd / w5) + dilog(1.0 / 2.0 + (1.0i / 2.0) * xd) * (log((w4 + yd) / w4) + log(1.0 - yd / w5))
                + dilog((-1.0i) * xd) * (log((w4 + yd) / w4) + log(1.0 - yd / w5)) + dilog(1.0i * xd) * (log((w4 + yd) / w4) + log(1.0 - yd / w5)) + dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + xd)) * (log((w4 + yd) / w4) + log(1.0 - yd / w5))
                + 1.0i * M_PI * H1(w4inv, -ydinv) * power_of<2>(log(w4 + yd)) * my_sign(-imydinv) + 1.0i * M_PI * H1(-w5inv, -ydinv) * power_of<2>(log(-w5 + yd)) * my_sign(-imydinv) - 1.0i * M_PI * H1(-w4inv, ydinv) * power_of<2>(log(-w4 - yd)) * my_sign(imydinv)
                - 1.0i * M_PI * H1(w5inv, ydinv) * power_of<2>(log(w5 - yd)) * my_sign(imydinv) + 1.0i * M_PI * H1(-2.0 / (-1.0 + w4), 2.0 / (1.0 + yd)) * power_of<2>(log((-w4 - yd) / 2.0)) * my_sign(2.0 * imag(1.0 / (1.0 + yd)))
                + 1.0i * M_PI * H1(2.0 / (1.0 + w5), 2.0 / (1.0 + yd)) * power_of<2>(log((w5 - yd) / 2.0)) * my_sign(2.0 * imag(1.0 / (1.0 + yd))) - 1.0i * M_PI * H1((1.0 - w4) / (1.0 + yd), 2.0 / (1.0 + yd)) * power_of<2>(log(((1.0 + w4) * (1.0 + yd)) / (2.0 - 2.0 * w4))) * my_sign(2.0 * imag(1.0 / (1.0 + yd)))
                - 1.0i * M_PI * H1((1.0 + w5) / (1.0 + yd), 2.0 / (1.0 + yd)) * power_of<2>(log((-1.0 / 2.0) * (((-1.0 + w5) * (1.0 + yd)) / (1.0 + w5)))) * my_sign(2.0 * imag(1.0 / (1.0 + yd)))
                + log(1.0 - yd / w5) * ((-1.0 / 2.0) * power_of<2>(log(1.0 - 1.0i * xd)) - power_of<2>(log(1.0 + 1.0i * xd)) / 2.0 - 2.0 * ln2 * log((wx3 + xd) / wx3) + log(xd) * (-log((wx3 + xd) / wx3) - log((wx4 + xd) / wx4)) - 2.0 * ln2 * log((wx4 + xd) / wx4)
                    + log(1.0 - 1.0i * xd) * (ln2 + log(xd) + log((wx3 + xd) / (-1.0i + wx3)) + log((wx4 + xd) / (-1.0i + wx4))) + log(1.0 + 1.0i * xd) * (ln2 + log(xd) + log((wx3 + xd) / (1.0i + wx3)) + log((wx4 + xd) / (1.0i + wx4)))
                    + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4))
                + log((w4 + yd) / w4) * ((-1.0 / 2.0) * power_of<2>(log(1.0 - 1.0i * xd)) - power_of<2>(log(1.0 + 1.0i * xd)) / 2.0 + log(1.0 + 1.0i * xd) * (ln2 + log((wx3 - xd) / (-1.0i + wx3)) + log((wx4 - xd) / (-1.0i + wx4)) + log(xd))
                    + log(1.0 - 1.0i * xd) * (ln2 + log((wx3 - xd) / (1.0i + wx3)) + log((wx4 - xd) / (1.0i + wx4)) + log(xd)) - 2.0 * ln2 * log(1.0 - xd / wx3) + log(xd) * (-log(1.0 - xd / wx3) - log(1.0 - xd / wx4))
                    - 2.0 * ln2 * log(1.0 - xd / wx4) + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4)
                    + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) + (2.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))
                - 1.0i * M_PI * power_of<2>(log((1.0 + w4) / 2.0)) * my_sign((-1.0 / 2.0) * imag(w4)) * T(1.0, (1.0 - yd) / 2.0, (1.0 + w4) / 2.0) - 1.0i * M_PI * power_of<2>(log((1.0 - w5) / 2.0)) * my_sign(imag(w5) / 2.0) * T(1.0, (1.0 - yd) / 2.0, (1.0 - w5) / 2.0)
                + ln1pyd * ((2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd)) * log((w5 - yd) / (1.0 + w5)) + (2.0 * ln2 - log(1.0 - 1.0i * xd) - log(1.0 + 1.0i * xd) + log(xd)) * log((w4 + yd) / (-1.0 + w4))
                    + (2.0 * 1.0i) * M_PI * log((1.0 + w4) / 2.0) * my_sign(2.0 * imag(1.0 / (1.0 + yd))) * T(1.0, (1.0 - yd) / 2.0, (1.0 + w4) / 2.0) + (2.0 * 1.0i) * M_PI * log((1.0 - w5) / 2.0) * my_sign(2.0 * imag(1.0 / (1.0 + yd))) * T(1.0, (1.0 - yd) / 2.0, (1.0 - w5) / 2.0))
                + 1.0i * M_PI * power_of<2>(log((w5 - yd) / (1.0 + w5))) * my_sign(imag((1.0 + yd) / (1.0 + w5))) * T(1.0, (1.0 - yd) / 2.0, (w5 - yd) / (1.0 + w5)) + 1.0i * M_PI * power_of<2>(log((w4 + yd) / (-1.0 + w4))) * my_sign(imag((1.0 + yd) / (1.0 - w4))) * T(1.0, (1.0 - yd) / 2.0, (w4 + yd) / (-1.0 + w4))
                - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 + w4) * my_sign(-imag(w4inv)) * T(1.0, 1.0 - yd, 1.0 + w4) + 1.0i * M_PI * power_of<2>(log(1.0 + w4)) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd, 1.0 + w4) - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 - w5) * my_sign(imag(w5inv)) *
                T(1.0, 1.0 - yd, 1.0 - w5) + 1.0i * M_PI * power_of<2>(log(1.0 - w5)) * my_sign(imag(w5)) * T(1.0, 1.0 - yd, 1.0 - w5) + (4.0 * 1.0i) * M_PI * ln2 * log(1.0 - w4) * my_sign(imag(w4inv)) * T(1.0, 1.0 + yd, 1.0 - w4)
                - (2.0 * 1.0i) * M_PI * log(1.0 - w4) * log((1.0 + w4) / 2.0) * my_sign(imag(w4)) * T(1.0, 1.0 + yd, 1.0 - w4) + dilog((1.0 - w4) / 2.0) * (-log(w4 / (-1.0 + w4)) + log((w4 + yd) / (-1.0 + w4)) - (2.0 * 1.0i) * M_PI * my_sign(imag(w4)) * T(1.0, 1.0 + yd, 1.0 - w4))
                + (4.0 * 1.0i) * M_PI * ln2 * log(1.0 + w5) * my_sign(-imag(w5inv)) * T(1.0, 1.0 + yd, 1.0 + w5) - (2.0 * 1.0i) * M_PI * log((1.0 - w5) / 2.0) * log(1.0 + w5) * my_sign(-imag(w5)) * T(1.0, 1.0 + yd, 1.0 + w5)
                + log(1.0 - 1.0i * xd) * ((2.0 * 1.0i) * M_PI * log(1.0 + w4) * my_sign(-imag(w4inv)) * T(1.0, 1.0 - yd, 1.0 + w4) + (2.0 * 1.0i) * M_PI * log(1.0 - w5) * my_sign(imag(w5inv)) * T(1.0, 1.0 - yd, 1.0 - w5)
                    - (2.0 * 1.0i) * M_PI * log(1.0 - w4) * my_sign(imag(w4inv)) * T(1.0, 1.0 + yd, 1.0 - w4) - (2.0 * 1.0i) * M_PI * log(1.0 + w5) * my_sign(-imag(w5inv)) * T(1.0, 1.0 + yd, 1.0 + w5))
                + log(1.0 + 1.0i * xd) * ((2.0 * 1.0i) * M_PI * log(1.0 + w4) * my_sign(-imag(w4inv)) * T(1.0, 1.0 - yd, 1.0 + w4) + (2.0 * 1.0i) * M_PI * log(1.0 - w5) * my_sign(imag(w5inv)) * T(1.0, 1.0 - yd, 1.0 - w5)
                    - (2.0 * 1.0i) * M_PI * log(1.0 - w4) * my_sign(imag(w4inv)) * T(1.0, 1.0 + yd, 1.0 - w4) - (2.0 * 1.0i) * M_PI * log(1.0 + w5) * my_sign(-imag(w5inv)) * T(1.0, 1.0 + yd, 1.0 + w5))
                + log(xd) * ((-2.0 * 1.0i) * M_PI * log(1.0 + w4) * my_sign(-imag(w4inv)) * T(1.0, 1.0 - yd, 1.0 + w4) - (2.0 * 1.0i) * M_PI * log(1.0 - w5) * my_sign(imag(w5inv)) * T(1.0, 1.0 - yd, 1.0 - w5)
                    + (2.0 * 1.0i) * M_PI * log(1.0 - w4) * my_sign(imag(w4inv)) * T(1.0, 1.0 + yd, 1.0 - w4) + (2.0 * 1.0i) * M_PI * log(1.0 + w5) * my_sign(-imag(w5inv)) * T(1.0, 1.0 + yd, 1.0 + w5))
                + dilog((1.0 + w5) / 2.0) * (-log(w5 / (1.0 + w5)) + log((w5 - yd) / (1.0 + w5)) - (2.0 * 1.0i) * M_PI * my_sign(-imag(w5)) * T(1.0, 1.0 + yd, 1.0 + w5))
                + li2half * (log(w4 / (-1.0 + w4)) + log(w5 / (1.0 + w5)) - log((w5 - yd) / (1.0 + w5)) - log((w4 + yd) / (-1.0 + w4)) - 2.0 * log((w4 + yd) / w4) - 2.0 * log(1.0 - yd / w5) + (2.0 * 1.0i) * M_PI * my_sign(imag(w4)) * T(1.0, 1.0 + yd, 1.0 - w4)
                    + (2.0 * 1.0i) * M_PI * my_sign(-imag(w5)) * T(1.0, 1.0 + yd, 1.0 + w5)) + power_of<2>(log((w4 + yd) / w4)) * ((-1.0i) * M_PI * my_sign(-imag(yd / w4)) * T(1.0, 1.0 - yd, (w4 + yd) / w4) + 1.0i * M_PI * my_sign(-imag(yd / w4)) * T(1.0, 1.0 + yd, (w4 + yd) / w4))
                + power_of<2>(log(1.0 - yd / w5)) * ((-1.0i) * M_PI * my_sign(imag(yd / w5)) * T(1.0, 1.0 - yd, 1.0 - yd / w5) + 1.0i * M_PI * my_sign(imag(yd / w5)) * T(1.0, 1.0 + yd, 1.0 - yd / w5))) / denom2);

            const complex<double> f29dPart44 = num18 * ((li3half - trilog(1.0 / (1.0 - w7)) - trilog(-w7inv) + trilog((-1.0 + w7) / (2.0 * w7)) - trilog((-1.0 + w7) / w7) + trilog(-w7) - trilog(w7) + trilog((1.0 - w7) / (1.0 + w7))
                - trilog(1.0 + ydinv) + trilog((-1.0 + yd) / yd) - trilog(-yd) + trilog(yd) - trilog(((-1.0 + w7) * (-1.0 + yd)) / ((1.0 + w7) * (1.0 + yd))) - trilog((1.0 + yd) / 2.0) + trilog((1.0 + yd) / (1.0 - w7))
                + trilog((-1.0 + w7) / (w7 + yd)) - trilog((-1.0 / 2.0) * (((-1.0 + w7) * (-1.0 + yd)) / (w7 + yd))) + trilog((1.0 + yd) / (w7 + yd)) - trilog((w7 * (1.0 + yd)) / (w7 + yd)) + trilog((w7 - w7 * yd) / (w7 + yd))
                + dilog(-w7inv) * lnhalf - dilog((-1.0 + w7) / (2.0 * w7)) * lnhalf - dilog((1.0 - w7) / (1.0 + w7)) * lnhalf + (power_of<2>(lnhalf) * log(1.0 + w7inv)) / 2.0 - (pisqu * log(2.0 / w7)) / 6.0 - power_of<3>(log(2.0 / w7)) / 6.0
                + (pisqu * log((2.0 * (-1.0 + w7)) / (1.0 + w7))) / 6.0 + power_of<3>(log((2.0 * (-1.0 + w7)) / (1.0 + w7))) / 6.0 + (power_of<2>(lnhalf) * log(w7 / (1.0 + w7))) / 2.0 - (power_of<2>(lnhalf) * log((2.0 * w7) / (1.0 + w7))) / 2.0 - (power_of<2>(lnhalf) * log((1.0 + w7) / (2.0 * w7))) / 2.0
                + dilog((1.0 + yd) / (1.0 - w7)) * (-4.0 * ln2 + 2.0 * log(1.0 - 1.0i * xd) + 2.0 * log(1.0 + 1.0i * xd) - 2.0 * log(xd)) + dilog(1.0 / (1.0 + w7)) * (lnhalf - 4.0 * ln2 + 2.0 * log(1.0 - 1.0i * xd) + 2.0 * log(1.0 + 1.0i * xd) - 2.0 * log(xd))
                + dilog(1.0 / (1.0 - w7)) * (4.0 * ln2 - 2.0 * log(1.0 - 1.0i * xd) - 2.0 * log(1.0 + 1.0i * xd) + 2.0 * log(xd)) + dilog(((-1.0 + w7) * (-1.0 + yd)) / ((1.0 + w7) * (1.0 + yd))) * log((1.0 - yd) / 2.0)
                + dilog((-1.0 / 2.0) * (((-1.0 + w7) * (-1.0 + yd)) / (w7 + yd))) * log((1.0 - yd) / 2.0) - dilog((-1.0 + yd) / yd) * ln1myd - dilog((w7 - w7 * yd) / (w7 + yd)) * ln1myd
                + dilog((-1.0 + yd) / (w7 + yd)) * (-log((1.0 - yd) / 2.0) + ln1myd) - power_of<3>(log(-ydinv)) / 6.0 + (pisqu / 6.0 - power_of<2>(ln1myd) / 2.0) * log(ydinv) + power_of<3>(log(ydinv)) / 6.0
                - (pisqu * log((2.0 * (-1.0 + w7)) / ((1.0 + w7) * (1.0 + yd)))) / 6.0 - power_of<3>(log((2.0 * (-1.0 + w7)) / ((1.0 + w7) * (1.0 + yd)))) / 6.0 + dilog((1.0 - yd) / (1.0 + w7)) * (4.0 * ln2 - 2.0 * log(1.0 - 1.0i * xd) - 2.0 * log(1.0 + 1.0i * xd) + 2.0 * log(xd)
                    - log((1.0 - yd) / 2.0) + ln1myd - ln1pyd) + dilog(2.0 / (1.0 + w7)) * ln1pyd + dilog(1.0 + ydinv) * ln1pyd - dilog((1.0 + yd) / (w7 + yd)) * ln1pyd
                + dilog((w7 * (1.0 + yd)) / (w7 + yd)) * ln1pyd + log(-ydinv) * ((-1.0 / 6.0) * pisqu + power_of<2>(ln1pyd) / 2.0) + (pisqu * log(-1.0 / (w7 + yd))) / 6.0 + power_of<3>(log(-1.0 / (w7 + yd))) / 6.0 - (pisqu * log(1.0 / (w7 + yd))) / 6.0
                - power_of<3>(log(1.0 / (w7 + yd))) / 6.0 + (pisqu * log(2.0 / (w7 + yd))) / 6.0 + power_of<3>(log(2.0 / (w7 + yd))) / 6.0 + dilog(-w7) * log((w7 + yd) / w7) - dilog(w7) * log((w7 + yd) / w7) + 4.0 * dilog(-1.0i / (-1.0i + wx3)) * log((w7 + yd) / w7)
                + 4.0 * dilog(1.0i / (1.0i + wx3)) * log((w7 + yd) / w7) + 4.0 * dilog(-1.0i / (-1.0i + wx4)) * log((w7 + yd) / w7) + 4.0 * dilog(1.0i / (1.0i + wx4)) * log((w7 + yd) / w7) - 2.0 * dilog((1.0i - xd) / (1.0i + wx3)) * log((w7 + yd) / w7)
                - 2.0 * dilog((1.0i - xd) / (1.0i + wx4)) * log((w7 + yd) / w7) - 4.0 * dilog(1.0 / 2.0 + (1.0i / 2.0) * xd) * log((w7 + yd) / w7) - 4.0 * dilog((-1.0i) * xd) * log((w7 + yd) / w7) - 4.0 * dilog(1.0i * xd) * log((w7 + yd) / w7)
                + 2.0 * dilog(-(xd / wx3)) * log((w7 + yd) / w7) + 2.0 * dilog(xd / wx3) * log((w7 + yd) / w7) + 2.0 * dilog(-(xd / wx4)) * log((w7 + yd) / w7) + 2.0 * dilog(xd / wx4) * log((w7 + yd) / w7)
                - 2.0 * dilog((-1.0i + xd) / (-1.0i + wx3)) * log((w7 + yd) / w7) - 2.0 * dilog((-1.0i + xd) / (-1.0i + wx4)) * log((w7 + yd) / w7) - 4.0 * dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + xd)) * log((w7 + yd) / w7)
                - 2.0 * dilog((1.0i + xd) / (1.0i - wx3)) * log((w7 + yd) / w7) - 2.0 * dilog((1.0i + xd) / (1.0i + wx3)) * log((w7 + yd) / w7) - 2.0 * dilog((1.0i + xd) / (1.0i - wx4)) * log((w7 + yd) / w7)
                - 2.0 * dilog((1.0i + xd) / (1.0i + wx4)) * log((w7 + yd) / w7) + power_of<2>(ln1myd) * (log((1.0 + w7) / (w7 + yd)) / 2.0 - log(((1.0 + w7) * yd) / (w7 + yd)) / 2.0 + log((w7 + yd) / (1.0 + w7)) / 2.0)
                + (4.0 * ln2 - 2.0 * log(1.0 - 1.0i * xd) - 2.0 * log(1.0 + 1.0i * xd) + 2.0 * log(xd)) * ln1myd * log((w7 + yd) / (1.0 + w7))
                + power_of<2>(log((1.0 - yd) / 2.0)) * ((-1.0 / 2.0) * log((1.0 + w7) / (w7 + yd)) + log(((1.0 + w7) * (1.0 + yd)) / (2.0 * (w7 + yd))) / 2.0 - log((w7 + yd) / (1.0 + w7)) / 2.0 + log((2.0 * (w7 + yd)) / ((1.0 + w7) * (1.0 + yd))) / 2.0)
                + power_of<2>(ln1pyd) * ((-1.0 / 2.0) * log((-1.0 + w7) / (w7 + yd)) + log((yd - w7 * yd) / (w7 + yd)) / 2.0) - 1.0i * M_PI * H1(w7inv, -ydinv) * power_of<2>(log(w7 + yd)) * my_sign(-imydinv)
                + 1.0i * M_PI * H1(-w7inv, ydinv) * power_of<2>(log(-w7 - yd)) * my_sign(imydinv) - 1.0i * M_PI * H1(-2.0 / (-1.0 + w7), 2.0 / (1.0 + yd)) * power_of<2>(log((-w7 - yd) / 2.0)) * my_sign(2.0 * imag(1.0 / (1.0 + yd)))
                + 1.0i * M_PI * H1((1.0 - w7) / (1.0 + yd), 2.0 / (1.0 + yd)) * power_of<2>(log(((1.0 + w7) * (1.0 + yd)) / (2.0 - 2.0 * w7))) * my_sign(2.0 * imag(1.0 / (1.0 + yd)))
                + log((w7 + yd) / w7) * (2.0 * power_of<2>(log(1.0 - 1.0i * xd)) + 2.0 * power_of<2>(log(1.0 + 1.0i * xd)) + 4.0 * ln2 * log((wx3 + xd) / wx3) + 4.0 * ln2 * log((wx4 + xd) / wx4)
                    + log(1.0 - 1.0i * xd) * (-4.0 * ln2 - 2.0 * log((wx3 - xd) / (1.0i + wx3)) - 2.0 * log((wx4 - xd) / (1.0i + wx4)) - 4.0 * log(xd) - 2.0 * log((wx3 + xd) / (-1.0i + wx3)) - 2.0 * log((wx4 + xd) / (-1.0i + wx4)))
                    + log(1.0 + 1.0i * xd) * (-4.0 * ln2 - 2.0 * log((wx3 - xd) / (-1.0i + wx3)) - 2.0 * log((wx4 - xd) / (-1.0i + wx4)) - 4.0 * log(xd) - 2.0 * log((wx3 + xd) / (1.0i + wx3)) - 2.0 * log((wx4 + xd) / (1.0i + wx4)))
                    + 4.0 * ln2 * log(1.0 - xd / wx3) + 4.0 * ln2 * log(1.0 - xd / wx4) + log(xd) * (2.0 * log((wx3 + xd) / wx3) + 2.0 * log((wx4 + xd) / wx4) + 2.0 * log(1.0 - xd / wx3) + 2.0 * log(1.0 - xd / wx4))
                    - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                    - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                    - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))
                + 1.0i * M_PI * power_of<2>(log((1.0 + w7) / 2.0)) * my_sign((-1.0 / 2.0) * imag(w7)) * T(1.0, (1.0 - yd) / 2.0, (1.0 + w7) / 2.0) + ln1pyd * ((-4.0 * ln2 + 2.0 * log(1.0 - 1.0i * xd) + 2.0 * log(1.0 + 1.0i * xd) - 2.0 * log(xd)) * log((w7 + yd) / (-1.0 + w7))
                    - log((1.0 - yd) / 2.0) * log((w7 + yd) / (1.0 + w7)) - (2.0 * 1.0i) * M_PI * log((1.0 + w7) / 2.0) * my_sign(2.0 * imag(1.0 / (1.0 + yd))) * T(1.0, (1.0 - yd) / 2.0, (1.0 + w7) / 2.0)) - 1.0i * M_PI * power_of<2>(log((w7 + yd) / (-1.0 + w7))) * my_sign(imag((1.0 + yd) / (1.0 - w7))) *
                T(1.0, (1.0 - yd) / 2.0, (w7 + yd) / (-1.0 + w7)) + (8.0 * 1.0i) * M_PI * ln2 * log(1.0 + w7) * my_sign(-imag(w7inv)) * T(1.0, 1.0 - yd, 1.0 + w7) - 1.0i * M_PI * power_of<2>(log(1.0 + w7)) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd, 1.0 + w7)
                - (8.0 * 1.0i) * M_PI * ln2 * log(1.0 - w7) * my_sign(imag(w7inv)) * T(1.0, 1.0 + yd, 1.0 - w7) + (2.0 * 1.0i) * M_PI * log(1.0 - w7) * log((1.0 + w7) / 2.0) * my_sign(imag(w7)) * T(1.0, 1.0 + yd, 1.0 - w7)
                + log(xd) * ((4.0 * 1.0i) * M_PI * log(1.0 + w7) * my_sign(-imag(w7inv)) * T(1.0, 1.0 - yd, 1.0 + w7) - (4.0 * 1.0i) * M_PI * log(1.0 - w7) * my_sign(imag(w7inv)) * T(1.0, 1.0 + yd, 1.0 - w7))
                + log(1.0 - 1.0i * xd) * ((-4.0 * 1.0i) * M_PI * log(1.0 + w7) * my_sign(-imag(w7inv)) * T(1.0, 1.0 - yd, 1.0 + w7) + (4.0 * 1.0i) * M_PI * log(1.0 - w7) * my_sign(imag(w7inv)) * T(1.0, 1.0 + yd, 1.0 - w7))
                + log(1.0 + 1.0i * xd) * ((-4.0 * 1.0i) * M_PI * log(1.0 + w7) * my_sign(-imag(w7inv)) * T(1.0, 1.0 - yd, 1.0 + w7) + (4.0 * 1.0i) * M_PI * log(1.0 - w7) * my_sign(imag(w7inv)) * T(1.0, 1.0 + yd, 1.0 - w7))
                + li2half * (-log(w7 / (-1.0 + w7)) + log((w7 + yd) / (-1.0 + w7)) + 8.0 * log((w7 + yd) / w7) - (2.0 * 1.0i) * M_PI * my_sign(imag(w7)) * T(1.0, 1.0 + yd, 1.0 - w7))
                + dilog((1.0 - w7) / 2.0) * (log(w7 / (-1.0 + w7)) - log((w7 + yd) / (-1.0 + w7)) + (2.0 * 1.0i) * M_PI * my_sign(imag(w7)) * T(1.0, 1.0 + yd, 1.0 - w7))
                + power_of<2>(log((w7 + yd) / w7)) * (1.0i * M_PI * my_sign(-imag(yd / w7)) * T(1.0, 1.0 - yd, (w7 + yd) / w7) - 1.0i * M_PI * my_sign(-imag(yd / w7)) * T(1.0, 1.0 + yd, (w7 + yd) / w7))) / denom2);

            const complex<double> f29dPart45 = num19 * ((li3half - trilog(1.0 / (1.0 - w7)) - trilog(-w7inv) + trilog((-1.0 + w7) / (2.0 * w7)) - trilog((-1.0 + w7) / w7) + trilog(-w7) - trilog(w7) + trilog((1.0 - w7) / (1.0 + w7))
                + trilog(1.0 + ydinv) - trilog((1.0 - yd) / 2.0) + trilog((-1.0 + w7) / (w7 - yd)) + trilog((-1.0 + yd) / (-1.0 + w7)) - trilog((-1.0 + yd) / yd) + trilog(-yd) - trilog(yd) - trilog(((-1.0 + w7) * (1.0 + yd)) / (2.0 * (w7 - yd)))
                + trilog((w7 * (1.0 + yd)) / (w7 - yd)) - trilog(((-1.0 + w7) * (1.0 + yd)) / ((1.0 + w7) * (-1.0 + yd))) + trilog((-1.0 + yd) / (-w7 + yd)) - trilog((w7 - w7 * yd) / (w7 - yd)) + dilog(-w7inv) * lnhalf
                - dilog((-1.0 + w7) / (2.0 * w7)) * lnhalf - dilog((1.0 - w7) / (1.0 + w7)) * lnhalf + (power_of<2>(lnhalf) * log(1.0 + w7inv)) / 2.0 - (pisqu * log(2.0 / w7)) / 6.0 - power_of<3>(log(2.0 / w7)) / 6.0 + (pisqu * log((2.0 * (-1.0 + w7)) / (1.0 + w7))) / 6.0
                + power_of<3>(log((2.0 * (-1.0 + w7)) / (1.0 + w7))) / 6.0 + (power_of<2>(lnhalf) * log(w7 / (1.0 + w7))) / 2.0 - (power_of<2>(lnhalf) * log((2.0 * w7) / (1.0 + w7))) / 2.0 - (power_of<2>(lnhalf) * log((1.0 + w7) / (2.0 * w7))) / 2.0
                + dilog((-1.0 + yd) / (-1.0 + w7)) * (-4.0 * ln2 + 2.0 * log(1.0 - 1.0i * xd) + 2.0 * log(1.0 + 1.0i * xd) - 2.0 * log(xd)) + dilog(1.0 / (1.0 + w7)) * (lnhalf - 4.0 * ln2 + 2.0 * log(1.0 - 1.0i * xd) + 2.0 * log(1.0 + 1.0i * xd) - 2.0 * log(xd))
                + dilog(1.0 / (1.0 - w7)) * (4.0 * ln2 - 2.0 * log(1.0 - 1.0i * xd) - 2.0 * log(1.0 + 1.0i * xd) + 2.0 * log(xd)) + dilog(2.0 / (1.0 + w7)) * ln1myd + dilog((-1.0 + yd) / yd) * ln1myd
                - dilog((-1.0 + yd) / (-w7 + yd)) * ln1myd + dilog((w7 - w7 * yd) / (w7 - yd)) * ln1myd - (pisqu * log(1.0 / (w7 - yd))) / 6.0 - power_of<3>(log(1.0 / (w7 - yd))) / 6.0 + (pisqu * log(2.0 / (w7 - yd))) / 6.0
                + power_of<3>(log(2.0 / (w7 - yd))) / 6.0 - (pisqu * log((-2.0 * (-1.0 + w7)) / ((1.0 + w7) * (-1.0 + yd)))) / 6.0 - power_of<3>(log((-2.0 * (-1.0 + w7)) / ((1.0 + w7) * (-1.0 + yd)))) / 6.0 + power_of<3>(log(-ydinv)) / 6.0 - (pisqu * log(ydinv)) / 6.0 - power_of<3>(log(ydinv)) / 6.0
                + power_of<2>(ln1myd) * ((-1.0 / 2.0) * log((-1.0 + w7) / (w7 - yd)) + log(ydinv) / 2.0 + log(((-1.0 + w7) * yd) / (w7 - yd)) / 2.0) + dilog(((-1.0 + w7) * (1.0 + yd)) / (2.0 * (w7 - yd))) * log((1.0 + yd) / 2.0)
                + dilog(((-1.0 + w7) * (1.0 + yd)) / ((1.0 + w7) * (-1.0 + yd))) * log((1.0 + yd) / 2.0) + ((-1.0 / 2.0) * log((1.0 + w7) / (w7 - yd)) + log((-2.0 * (w7 - yd)) / ((1.0 + w7) * (-1.0 + yd))) / 2.0
                    + log((-1.0 / 2.0) * (((1.0 + w7) * (-1.0 + yd)) / (w7 - yd))) / 2.0) * power_of<2>(log((1.0 + yd) / 2.0)) - dilog(1.0 + ydinv) * ln1pyd - dilog((w7 * (1.0 + yd)) / (w7 - yd)) * ln1pyd
                + (log((1.0 + w7) / (w7 - yd)) * power_of<2>(ln1pyd)) / 2.0 + dilog(-((1.0 + yd) / (w7 - yd))) * (-log((1.0 + yd) / 2.0) + ln1pyd)
                + dilog((1.0 + yd) / (1.0 + w7)) * (4.0 * ln2 - 2.0 * log(1.0 - 1.0i * xd) - 2.0 * log(1.0 + 1.0i * xd) + 2.0 * log(xd) - ln1myd - log((1.0 + yd) / 2.0) + ln1pyd) + log(-ydinv) * (pisqu / 6.0 - power_of<2>(ln1pyd) / 2.0)
                + log((w7 - yd) / (1.0 + w7)) * ((-1.0 / 2.0) * power_of<2>(log((1.0 + yd) / 2.0)) + (4.0 * ln2 - 2.0 * log(1.0 - 1.0i * xd) - 2.0 * log(1.0 + 1.0i * xd) + 2.0 * log(xd)) * ln1pyd + power_of<2>(ln1pyd) / 2.0) + (pisqu * log(1.0 / (-w7 + yd))) / 6.0
                + power_of<3>(log(1.0 / (-w7 + yd))) / 6.0 - (power_of<2>(ln1pyd) * log(((1.0 + w7) * yd) / (-w7 + yd))) / 2.0 + dilog(-w7) * log(1.0 - yd / w7) - dilog(w7) * log(1.0 - yd / w7) + 4.0 * dilog(-1.0i / (-1.0i + wx3)) * log(1.0 - yd / w7)
                + 4.0 * dilog(1.0i / (1.0i + wx3)) * log(1.0 - yd / w7) + 4.0 * dilog(-1.0i / (-1.0i + wx4)) * log(1.0 - yd / w7) + 4.0 * dilog(1.0i / (1.0i + wx4)) * log(1.0 - yd / w7) - 2.0 * dilog((1.0i - xd) / (1.0i + wx3)) * log(1.0 - yd / w7)
                - 2.0 * dilog((1.0i - xd) / (1.0i + wx4)) * log(1.0 - yd / w7) - 4.0 * dilog(1.0 / 2.0 + (1.0i / 2.0) * xd) * log(1.0 - yd / w7) - 4.0 * dilog((-1.0i) * xd) * log(1.0 - yd / w7) - 4.0 * dilog(1.0i * xd) * log(1.0 - yd / w7) + 2.0 * dilog(-(xd / wx3)) * log(1.0 - yd / w7)
                + 2.0 * dilog(xd / wx3) * log(1.0 - yd / w7) + 2.0 * dilog(-(xd / wx4)) * log(1.0 - yd / w7) + 2.0 * dilog(xd / wx4) * log(1.0 - yd / w7) - 2.0 * dilog((-1.0i + xd) / (-1.0i + wx3)) * log(1.0 - yd / w7)
                - 2.0 * dilog((-1.0i + xd) / (-1.0i + wx4)) * log(1.0 - yd / w7) - 4.0 * dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + xd)) * log(1.0 - yd / w7) - 2.0 * dilog((1.0i + xd) / (1.0i - wx3)) * log(1.0 - yd / w7) - 2.0 * dilog((1.0i + xd) / (1.0i + wx3)) * log(1.0 - yd / w7)
                - 2.0 * dilog((1.0i + xd) / (1.0i - wx4)) * log(1.0 - yd / w7) - 2.0 * dilog((1.0i + xd) / (1.0i + wx4)) * log(1.0 - yd / w7) + 1.0i * M_PI * H1((-1.0 + w7) / (-1.0 + yd), -2.0 / (-1.0 + yd)) * power_of<2>(log(((1.0 + w7) * (-1.0 + yd)) / (2.0 * (-1.0 + w7)))) *
                my_sign(2.0 * imag(1.0 / (1.0 - yd))) - 1.0i * M_PI * H1(-2.0 / (-1.0 + w7), -2.0 / (-1.0 + yd)) * power_of<2>(log((-w7 + yd) / 2.0)) * my_sign(2.0 * imag(1.0 / (1.0 - yd))) + 1.0i * M_PI * H1(-w7inv, -ydinv) * power_of<2>(log(-w7 + yd)) * my_sign(-imydinv)
                - 1.0i * M_PI * H1(w7inv, ydinv) * power_of<2>(log(w7 - yd)) * my_sign(imydinv) + log(1.0 - yd / w7) * (2.0 * power_of<2>(log(1.0 - 1.0i * xd)) + 2.0 * power_of<2>(log(1.0 + 1.0i * xd)) + 4.0 * ln2 * log((wx3 + xd) / wx3) + 4.0 * ln2 * log((wx4 + xd) / wx4)
                    + log(1.0 - 1.0i * xd) * (-4.0 * ln2 - 2.0 * log((wx3 - xd) / (1.0i + wx3)) - 2.0 * log((wx4 - xd) / (1.0i + wx4)) - 4.0 * log(xd) - 2.0 * log((wx3 + xd) / (-1.0i + wx3)) - 2.0 * log((wx4 + xd) / (-1.0i + wx4)))
                    + log(1.0 + 1.0i * xd) * (-4.0 * ln2 - 2.0 * log((wx3 - xd) / (-1.0i + wx3)) - 2.0 * log((wx4 - xd) / (-1.0i + wx4)) - 4.0 * log(xd) - 2.0 * log((wx3 + xd) / (1.0i + wx3)) - 2.0 * log((wx4 + xd) / (1.0i + wx4)))
                    + 4.0 * ln2 * log(1.0 - xd / wx3) + 4.0 * ln2 * log(1.0 - xd / wx4) + log(xd) * (2.0 * log((wx3 + xd) / wx3) + 2.0 * log((wx4 + xd) / wx4) + 2.0 * log(1.0 - xd / wx3) + 2.0 * log(1.0 - xd / wx4))
                    - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                    - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                    - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))
                - (8.0 * 1.0i) * M_PI * ln2 * log(1.0 - w7) * my_sign(imag(w7inv)) * T(1.0, 1.0 - yd, 1.0 - w7) + (2.0 * 1.0i) * M_PI * log(1.0 - w7) * log((1.0 + w7) / 2.0) * my_sign(imag(w7)) * T(1.0, 1.0 - yd, 1.0 - w7)
                + li2half * (-log(w7 / (-1.0 + w7)) + log((w7 - yd) / (-1.0 + w7)) + 8.0 * log(1.0 - yd / w7) - (2.0 * 1.0i) * M_PI * my_sign(imag(w7)) * T(1.0, 1.0 - yd, 1.0 - w7))
                + dilog((1.0 - w7) / 2.0) * (log(w7 / (-1.0 + w7)) - log((w7 - yd) / (-1.0 + w7)) + (2.0 * 1.0i) * M_PI * my_sign(imag(w7)) * T(1.0, 1.0 - yd, 1.0 - w7)) + 1.0i * M_PI * power_of<2>(log((1.0 + w7) / 2.0)) * my_sign((-1.0 / 2.0) * imag(w7)) * T(1.0, (1.0 + yd) / 2.0, (1.0 + w7) / 2.0)
                + ln1myd * ((-4.0 * ln2 + 2.0 * log(1.0 - 1.0i * xd) + 2.0 * log(1.0 + 1.0i * xd) - 2.0 * log(xd)) * log((w7 - yd) / (-1.0 + w7)) - log((w7 - yd) / (1.0 + w7)) * log((1.0 + yd) / 2.0)
                    - (2.0 * 1.0i) * M_PI * log((1.0 + w7) / 2.0) * my_sign(2.0 * imag(1.0 / (1.0 - yd))) * T(1.0, (1.0 + yd) / 2.0, (1.0 + w7) / 2.0)) - 1.0i * M_PI * power_of<2>(log((w7 - yd) / (-1.0 + w7))) * my_sign(imag((-1.0 + yd) / (-1.0 + w7))) * T(1.0, (1.0 + yd) / 2.0, (w7 - yd) / (-1.0 + w7))
                + (8.0 * 1.0i) * M_PI * ln2 * log(1.0 + w7) * my_sign(-imag(w7inv)) * T(1.0, 1.0 + yd, 1.0 + w7) - 1.0i * M_PI * power_of<2>(log(1.0 + w7)) * my_sign(-imag(w7)) * T(1.0, 1.0 + yd, 1.0 + w7)
                + log(1.0 - 1.0i * xd) * ((4.0 * 1.0i) * M_PI * log(1.0 - w7) * my_sign(imag(w7inv)) * T(1.0, 1.0 - yd, 1.0 - w7) - (4.0 * 1.0i) * M_PI * log(1.0 + w7) * my_sign(-imag(w7inv)) * T(1.0, 1.0 + yd, 1.0 + w7))
                + log(1.0 + 1.0i * xd) * ((4.0 * 1.0i) * M_PI * log(1.0 - w7) * my_sign(imag(w7inv)) * T(1.0, 1.0 - yd, 1.0 - w7) - (4.0 * 1.0i) * M_PI * log(1.0 + w7) * my_sign(-imag(w7inv)) * T(1.0, 1.0 + yd, 1.0 + w7))
                + log(xd) * ((-4.0 * 1.0i) * M_PI * log(1.0 - w7) * my_sign(imag(w7inv)) * T(1.0, 1.0 - yd, 1.0 - w7) + (4.0 * 1.0i) * M_PI * log(1.0 + w7) * my_sign(-imag(w7inv)) * T(1.0, 1.0 + yd, 1.0 + w7))
                + power_of<2>(log(1.0 - yd / w7)) * ((-1.0i) * M_PI * my_sign(imag(yd / w7)) * T(1.0, 1.0 - yd, 1.0 - yd / w7) + 1.0i * M_PI * my_sign(imag(yd / w7)) * T(1.0, 1.0 + yd, 1.0 - yd / w7))) / denom2);

            const complex<double> f29dPart5 = (num20 * (16.0 * power_of<2>(li2half) + (dilog((w5 - yd) / (-1.0 + w5)) + dilog(-(yd / w5)) + dilog(yd / w5) + dilog((w5 + yd) / (-1.0 + w5)) + log((1.0 + yd) / (1.0 - w5)) * log((w5 + yd) / w5)
                    + log((-1.0 + yd) / (-1.0 + w5)) * log(1.0 - yd / w5)) * (-2.0 * power_of<2>(log(1.0 - 1.0i * xd)) - 2.0 * power_of<2>(log(1.0 + 1.0i * xd)) + log(xd) * (-4.0 * log((wx3 + xd) / wx3) - 4.0 * log((wx4 + xd) / wx4))
                    + log(1.0 - 1.0i * xd) * (-4.0 * ln2 + 4.0 * log(xd) + 4.0 * log((wx3 + xd) / (-1.0i + wx3)) + 4.0 * log((wx4 + xd) / (-1.0i + wx4))) + log(1.0 + 1.0i * xd) * (-4.0 * ln2 + 4.0 * log(xd) + 4.0 * log((wx3 + xd) / (1.0i + wx3))
                    + 4.0 * log((wx4 + xd) / (1.0i + wx4))) + (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) + (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    + (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) + (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4))
                + dilog(w5 / (-1.0 + w5)) * (8.0 * dilog(-1.0i / (-1.0i + wx3)) + 8.0 * dilog(1.0i / (1.0i + wx3)) + 8.0 * dilog(-1.0i / (-1.0i + wx4)) + 8.0 * dilog(1.0i / (1.0i + wx4)) - 8.0 * dilog((1.0i - xd) / (1.0i + wx3)) - 8.0 * dilog((1.0i - xd) / (1.0i + wx4))
                    - 8.0 * dilog(1.0 / 2.0 + (1.0i / 2.0) * xd) - 8.0 * dilog((-1.0i) * xd) - 8.0 * dilog(1.0i * xd) + 8.0 * dilog(-(xd / wx3)) + 8.0 * dilog(-(xd / wx4)) - 8.0 * dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + xd)) - 8.0 * dilog((1.0i + xd) / (1.0i - wx3))
                    - 8.0 * dilog((1.0i + xd) / (1.0i - wx4)) + 4.0 * power_of<2>(log(1.0 - 1.0i * xd)) + 4.0 * power_of<2>(log(1.0 + 1.0i * xd)) + log(xd) * (8.0 * log((wx3 + xd) / wx3) + 8.0 * log((wx4 + xd) / wx4))
                    + log(1.0 - 1.0i * xd) * (8.0 * ln2 - 8.0 * log(xd) - 8.0 * log((wx3 + xd) / (-1.0i + wx3)) - 8.0 * log((wx4 + xd) / (-1.0i + wx4))) + log(1.0 + 1.0i * xd) * (8.0 * ln2 - 8.0 * log(xd) - 8.0 * log((wx3 + xd) / (1.0i + wx3))
                    - 8.0 * log((wx4 + xd) / (1.0i + wx4))) - (16.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) - (16.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    - (16.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) - (16.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4))
                + (dilog((w7 - yd) / (1.0 + w7)) + dilog(-yd) + dilog(yd) + dilog(-(yd / w7)) + dilog(yd / w7) + dilog((w7 + yd) / (1.0 + w7)) + log((1.0 - yd) / (1.0 + w7)) * log((w7 + yd) / w7)
                    + log((1.0 + yd) / (1.0 + w7)) * log(1.0 - yd / w7) - dilog((1.0 - yd) / 2.0) - dilog((1.0 + yd) / 2.0)) * (2.0 * power_of<2>(log(1.0 - 1.0i * xd)) + 2.0 * power_of<2>(log(1.0 + 1.0i * xd))
                    + log(1.0 - 1.0i * xd) * (4.0 * ln2 - 2.0 * log((wx3 - xd) / (1.0i + wx3)) - 2.0 * log((wx4 - xd) / (1.0i + wx4)) - 4.0 * log(xd) - 2.0 * log((wx3 + xd) / (-1.0i + wx3)) - 2.0 * log((wx4 + xd) / (-1.0i + wx4)))
                    + log(1.0 + 1.0i * xd) * (4.0 * ln2 - 2.0 * log((wx3 - xd) / (-1.0i + wx3)) - 2.0 * log((wx4 - xd) / (-1.0i + wx4)) - 4.0 * log(xd) - 2.0 * log((wx3 + xd) / (1.0i + wx3)) - 2.0 * log((wx4 + xd) / (1.0i + wx4)))
                    + log(xd) * (2.0 * log((wx3 + xd) / wx3) + 2.0 * log((wx4 + xd) / wx4) + 2.0 * log(1.0 - xd / wx3) + 2.0 * log(1.0 - xd / wx4)) - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3)
                    - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4)
                    - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4) - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3)
                    - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) - (4.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4)
                    - (4.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4)) + (dilog((w4 - yd) / (1.0 + w4)) + dilog(-(yd / w4)) + dilog(yd / w4) + dilog((w4 + yd) / (1.0 + w4))
                    + log((1.0 - yd) / (1.0 + w4)) * log((w4 + yd) / w4) + log((1.0 + yd) / (1.0 + w4)) * log(1.0 - yd / w4)) * (-2.0 * power_of<2>(log(1.0 - 1.0i * xd)) - 2.0 * power_of<2>(log(1.0 + 1.0i * xd))
                    + log(1.0 + 1.0i * xd) * (-4.0 * ln2 + 4.0 * log((wx3 - xd) / (-1.0i + wx3)) + 4.0 * log((wx4 - xd) / (-1.0i + wx4)) + 4.0 * log(xd)) + log(1.0 - 1.0i * xd) * (-4.0 * ln2 + 4.0 * log((wx3 - xd) / (1.0i + wx3))
                    + 4.0 * log((wx4 - xd) / (1.0i + wx4)) + 4.0 * log(xd)) + log(xd) * (-4.0 * log(1.0 - xd / wx3) - 4.0 * log(1.0 - xd / wx4)) + (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3)
                    + (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) + (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                    + (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4)) + dilog(w7 / (1.0 + w7)) * (-8.0 * dilog(-1.0i / (-1.0i + wx3)) - 8.0 * dilog(1.0i / (1.0i + wx3)) - 8.0 * dilog(-1.0i / (-1.0i + wx4)) - 8.0 * dilog(1.0i / (1.0i + wx4))
                    + 4.0 * dilog((1.0i - xd) / (1.0i + wx3)) + 4.0 * dilog((1.0i - xd) / (1.0i + wx4)) + 8.0 * dilog(1.0 / 2.0 + (1.0i / 2.0) * xd) + 8.0 * dilog((-1.0i) * xd) + 8.0 * dilog(1.0i * xd) - 4.0 * dilog(-(xd / wx3)) - 4.0 * dilog(xd / wx3) - 4.0 * dilog(-(xd / wx4))
                    - 4.0 * dilog(xd / wx4) + 4.0 * dilog((-1.0i + xd) / (-1.0i + wx3)) + 4.0 * dilog((-1.0i + xd) / (-1.0i + wx4)) + 8.0 * dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + xd)) + 4.0 * dilog((1.0i + xd) / (1.0i - wx3)) + 4.0 * dilog((1.0i + xd) / (1.0i + wx3))
                    + 4.0 * dilog((1.0i + xd) / (1.0i - wx4)) + 4.0 * dilog((1.0i + xd) / (1.0i + wx4)) - 4.0 * power_of<2>(log(1.0 - 1.0i * xd)) - 4.0 * power_of<2>(log(1.0 + 1.0i * xd)) + log(1.0 - 1.0i * xd) * (-8.0 * ln2 + 4.0 * log((wx3 - xd) / (1.0i + wx3)) + 4.0 * log((wx4 - xd) / (1.0i + wx4))
                    + 8.0 * log(xd) + 4.0 * log((wx3 + xd) / (-1.0i + wx3)) + 4.0 * log((wx4 + xd) / (-1.0i + wx4))) + log(1.0 + 1.0i * xd) * (-8.0 * ln2 + 4.0 * log((wx3 - xd) / (-1.0i + wx3)) + 4.0 * log((wx4 - xd) / (-1.0i + wx4)) + 8.0 * log(xd)
                    + 4.0 * log((wx3 + xd) / (1.0i + wx3)) + 4.0 * log((wx4 + xd) / (1.0i + wx4))) + log(xd) * (-4.0 * log((wx3 + xd) / wx3) - 4.0 * log((wx4 + xd) / wx4) - 4.0 * log(1.0 - xd / wx3) - 4.0 * log(1.0 - xd / wx4))
                    + (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) + (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                    + (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) + (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    + (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) + (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                    + (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) + (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))
                + dilog(w4 / (1.0 + w4)) * (8.0 * dilog(-1.0i / (-1.0i + wx3)) + 8.0 * dilog(1.0i / (1.0i + wx3)) + 8.0 * dilog(-1.0i / (-1.0i + wx4)) + 8.0 * dilog(1.0i / (1.0i + wx4)) - 8.0 * dilog(1.0 / 2.0 + (1.0i / 2.0) * xd) - 8.0 * dilog((-1.0i) * xd) - 8.0 * dilog(1.0i * xd)
                    + 8.0 * dilog(xd / wx3) + 8.0 * dilog(xd / wx4) - 8.0 * dilog((-1.0i + xd) / (-1.0i + wx3)) - 8.0 * dilog((-1.0i + xd) / (-1.0i + wx4)) - 8.0 * dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + xd)) - 8.0 * dilog((1.0i + xd) / (1.0i + wx3))
                    - 8.0 * dilog((1.0i + xd) / (1.0i + wx4)) + 4.0 * power_of<2>(log(1.0 - 1.0i * xd)) + 4.0 * power_of<2>(log(1.0 + 1.0i * xd)) + log(1.0 + 1.0i * xd) * (8.0 * ln2 - 8.0 * log((wx3 - xd) / (-1.0i + wx3)) - 8.0 * log((wx4 - xd) / (-1.0i + wx4)) - 8.0 * log(xd))
                    + log(1.0 - 1.0i * xd) * (8.0 * ln2 - 8.0 * log((wx3 - xd) / (1.0i + wx3)) - 8.0 * log((wx4 - xd) / (1.0i + wx4)) - 8.0 * log(xd)) + log(xd) * (8.0 * log(1.0 - xd / wx3) + 8.0 * log(1.0 - xd / wx4))
                    - (16.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) - (16.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4)
                    - (16.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) - (16.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))
                + (ln1pyd + ln1myd) * (2.0 * ln2 * power_of<2>(log(1.0 - 1.0i * xd)) + 2.0 * ln2 * power_of<2>(log(1.0 + 1.0i * xd)) + log(1.0 - 1.0i * xd) * (4.0 * ln2squ - 2.0 * ln2 * log((wx3 - xd) / (1.0i + wx3))
                    - 2.0 * ln2 * log((wx4 - xd) / (1.0i + wx4)) - 4.0 * ln2 * log(xd) - 2.0 * ln2 * log((wx3 + xd) / (-1.0i + wx3)) - 2.0 * ln2 * log((wx4 + xd) / (-1.0i + wx4)))
                    + log(1.0 + 1.0i * xd) * (4.0 * ln2squ - 2.0 * ln2 * log((wx3 - xd) / (-1.0i + wx3)) - 2.0 * ln2 * log((wx4 - xd) / (-1.0i + wx4)) - 4.0 * ln2 * log(xd) - 2.0 * ln2 * log((wx3 + xd) / (1.0i + wx3))
                    - 2.0 * ln2 * log((wx4 + xd) / (1.0i + wx4))) + log(xd) * (2.0 * ln2 * log((wx3 + xd) / wx3) + 2.0 * ln2 * log((wx4 + xd) / wx4) + 2.0 * ln2 * log(1.0 - xd / wx3) + 2.0 * ln2 * log(1.0 - xd / wx4))
                    - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                    - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                    - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) - (4.0 * 1.0i) * M_PI * ln2 * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4))
                - 16.0 * pisqu * log(1.0 + w4inv) * log(1.0 - 1.0i * wx3) * my_sign(-imag(w4)) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv)
                - 16.0 * pisqu * log(1.0 + w4inv) * log(1.0 - 1.0i * wx4) * my_sign(-imag(w4)) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv)
                - 16.0 * pisqu * log(1.0 + w4inv) * log(1.0 + 1.0i * wx3) * my_sign(-imag(w4)) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv)
                - 16.0 * pisqu * log(1.0 + w4inv) * log(1.0 + 1.0i * wx4) * my_sign(-imag(w4)) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv)
                - 16.0 * pisqu * log((-1.0 + w5) / w5) * log(1.0 + 1.0i * wx3) * my_sign(imag(w5)) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                - 16.0 * pisqu * log((-1.0 + w5) / w5) * log(1.0 + 1.0i * wx4) * my_sign(imag(w5)) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                - 16.0 * pisqu * log((-1.0 + w5) / w5) * log(1.0 - 1.0i * wx3) * my_sign(imag(w5)) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                - 16.0 * pisqu * log((-1.0 + w5) / w5) * log(1.0 - 1.0i * wx4) * my_sign(imag(w5)) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 - 1.0i * wx3) * my_sign(-imag(w7)) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 + 1.0i * wx3) * my_sign(-imag(w7)) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 - 1.0i * wx4) * my_sign(-imag(w7)) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 + 1.0i * wx4) * my_sign(-imag(w7)) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 - 1.0i * wx3) * my_sign(-imag(w7)) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 + 1.0i * wx3) * my_sign(-imag(w7)) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 - 1.0i * wx4) * my_sign(-imag(w7)) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 + 1.0i * wx4) * my_sign(-imag(w7)) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                - 16.0 * pisqu * log(1.0 + w4inv) * log(1.0 - 1.0i * wx3) * my_sign(-imag(w4)) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv) - 16.0 * pisqu * log(1.0 + w4inv) * log(1.0 - 1.0i * wx4) * my_sign(-imag(w4)) *
                my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv) - 16.0 * pisqu * log(1.0 + w4inv) * log(1.0 + 1.0i * wx3) * my_sign(-imag(w4)) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) *
                T(1.0, 1.0 - yd / w4, 1.0 + w4inv) - 16.0 * pisqu * log(1.0 + w4inv) * log(1.0 + 1.0i * wx4) * my_sign(-imag(w4)) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv)
                - 16.0 * pisqu * log((-1.0 + w5) / w5) * log(1.0 + 1.0i * wx3) * my_sign(imag(w5)) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) - 16.0 * pisqu * log((-1.0 + w5) / w5) * log(1.0 + 1.0i * wx4) * my_sign(imag(w5)) *
                my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) - 16.0 * pisqu * log((-1.0 + w5) / w5) * log(1.0 - 1.0i * wx3) * my_sign(imag(w5)) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) *
                T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) - 16.0 * pisqu * log((-1.0 + w5) / w5) * log(1.0 - 1.0i * wx4) * my_sign(imag(w5)) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5)
                + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 - 1.0i * wx3) * my_sign(-imag(w7)) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv) + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 + 1.0i * wx3) * my_sign(-imag(w7)) *
                my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv) + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 - 1.0i * wx4) * my_sign(-imag(w7)) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) *
                T(1.0, 1.0 - yd / w7, 1.0 + w7inv) + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 + 1.0i * wx4) * my_sign(-imag(w7)) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv)
                + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 - 1.0i * wx3) * my_sign(-imag(w7)) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv) + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 + 1.0i * wx3) * my_sign(-imag(w7)) *
                my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv) + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 - 1.0i * wx4) * my_sign(-imag(w7)) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) *
                T(1.0, 1.0 - yd / w7, 1.0 + w7inv) + 8.0 * pisqu * log(1.0 + w7inv) * log(1.0 + 1.0i * wx4) * my_sign(-imag(w7)) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv)
                + (dilog((-1.0i + xd) / (-1.0i + wx3)) + dilog((-1.0i + xd) / (-1.0i + wx4)) + dilog((1.0i + xd) / (1.0i + wx3)) + dilog((1.0i + xd) / (1.0i + wx4))) * (2.0 * dilog((1.0 - yd) / 2.0) + 4.0 * dilog((w4 - yd) / (1.0 + w4))
                    - 2.0 * dilog((w7 - yd) / (1.0 + w7)) - 2.0 * dilog(-yd) - 2.0 * dilog(yd) + 4.0 * dilog(-(yd / w4)) + 4.0 * dilog(yd / w4) - 2.0 * dilog(-(yd / w7)) - 2.0 * dilog(yd / w7) + 2.0 * dilog((1.0 + yd) / 2.0) + 4.0 * dilog((w4 + yd) / (1.0 + w4))
                    - 2.0 * dilog((w7 + yd) / (1.0 + w7)) - 2.0 * ln2 * ln1myd - 2.0 * ln2 * ln1pyd + 4.0 * log((1.0 - yd) / (1.0 + w4)) * log((w4 + yd) / w4) - 2.0 * log((1.0 - yd) / (1.0 + w7)) * log((w7 + yd) / w7)
                    + 4.0 * log((1.0 + yd) / (1.0 + w4)) * log(1.0 - yd / w4) - 2.0 * log((1.0 + yd) / (1.0 + w7)) * log(1.0 - yd / w7) + (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv)
                    - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv) + (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv)
                    - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv)) + (dilog((1.0i - xd) / (1.0i + wx3)) + dilog((1.0i - xd) / (1.0i + wx4)) + dilog((1.0i + xd) / (1.0i - wx3)) + dilog((1.0i + xd) / (1.0i - wx4))) *
                (2.0 * dilog((1.0 - yd) / 2.0) + 4.0 * dilog((w5 - yd) / (-1.0 + w5)) - 2.0 * dilog((w7 - yd) / (1.0 + w7)) - 2.0 * dilog(-yd) - 2.0 * dilog(yd) + 4.0 * dilog(-(yd / w5)) + 4.0 * dilog(yd / w5) - 2.0 * dilog(-(yd / w7)) - 2.0 * dilog(yd / w7)
                    + 2.0 * dilog((1.0 + yd) / 2.0) + 4.0 * dilog((w5 + yd) / (-1.0 + w5)) - 2.0 * dilog((w7 + yd) / (1.0 + w7)) - 2.0 * ln2 * ln1myd - 2.0 * ln2 * ln1pyd + 4.0 * log((1.0 + yd) / (1.0 - w5)) * log((w5 + yd) / w5)
                    - 2.0 * log((1.0 - yd) / (1.0 + w7)) * log((w7 + yd) / w7) + 4.0 * log((-1.0 + yd) / (-1.0 + w5)) * log(1.0 - yd / w5) - 2.0 * log((1.0 + yd) / (1.0 + w7)) * log(1.0 - yd / w7)
                    + (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5) - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                    + (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                + (dilog(xd / wx3) + dilog(xd / wx4)) * (-2.0 * dilog((1.0 - yd) / 2.0) - 4.0 * dilog((w4 - yd) / (1.0 + w4)) + 2.0 * dilog((w7 - yd) / (1.0 + w7)) + 2.0 * dilog(-yd) + 2.0 * dilog(yd) - 4.0 * dilog(-(yd / w4)) - 4.0 * dilog(yd / w4)
                    + 2.0 * dilog(-(yd / w7)) + 2.0 * dilog(yd / w7) - 2.0 * dilog((1.0 + yd) / 2.0) - 4.0 * dilog((w4 + yd) / (1.0 + w4)) + 2.0 * dilog((w7 + yd) / (1.0 + w7)) + 2.0 * ln2 * ln1myd + 2.0 * ln2 * ln1pyd
                    - 4.0 * log((1.0 - yd) / (1.0 + w4)) * log((w4 + yd) / w4) + 2.0 * log((1.0 - yd) / (1.0 + w7)) * log((w7 + yd) / w7) - 4.0 * log((1.0 + yd) / (1.0 + w4)) * log(1.0 - yd / w4) + 2.0 * log((1.0 + yd) / (1.0 + w7)) * log(1.0 - yd / w7)
                    - (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) + (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                    - (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv) + (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                + (power_of<2>(log(1.0 - 1.0i * xd)) + power_of<2>(log(1.0 + 1.0i * xd))) * ((-4.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) - (4.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                    + (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv) - (4.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv)
                    - (4.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) + (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                + (dilog(-(xd / wx3)) + dilog(-(xd / wx4))) * (-2.0 * dilog((1.0 - yd) / 2.0) - 4.0 * dilog((w5 - yd) / (-1.0 + w5)) + 2.0 * dilog((w7 - yd) / (1.0 + w7)) + 2.0 * dilog(-yd) + 2.0 * dilog(yd) - 4.0 * dilog(-(yd / w5)) - 4.0 * dilog(yd / w5)
                    + 2.0 * dilog(-(yd / w7)) + 2.0 * dilog(yd / w7) - 2.0 * dilog((1.0 + yd) / 2.0) - 4.0 * dilog((w5 + yd) / (-1.0 + w5)) + 2.0 * dilog((w7 + yd) / (1.0 + w7)) + 2.0 * ln2 * ln1myd + 2.0 * ln2 * ln1pyd
                    - 4.0 * log((1.0 + yd) / (1.0 - w5)) * log((w5 + yd) / w5) + 2.0 * log((1.0 - yd) / (1.0 + w7)) * log((w7 + yd) / w7) - 4.0 * log((-1.0 + yd) / (-1.0 + w5)) * log(1.0 - yd / w5) + 2.0 * log((1.0 + yd) / (1.0 + w7)) * log(1.0 - yd / w7)
                    - (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5) + (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                    - (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) + (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                + (dilog((-1.0i) * xd) + dilog(1.0 / 2.0 + (1.0i / 2.0) * xd) + dilog(1.0i * xd) + dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + xd)) - dilog(-1.0i / (-1.0i + wx3)) - dilog(1.0i / (1.0i + wx3)) - dilog(-1.0i / (-1.0i + wx4)) - dilog(1.0i / (1.0i + wx4))) *
                (4.0 * dilog((1.0 - yd) / 2.0) + 4.0 * dilog((w4 - yd) / (1.0 + w4)) + 4.0 * dilog((w5 - yd) / (-1.0 + w5)) - 4.0 * dilog((w7 - yd) / (1.0 + w7)) - 4.0 * dilog(-yd) - 4.0 * dilog(yd) + 4.0 * dilog(-(yd / w4)) + 4.0 * dilog(yd / w4)
                    + 4.0 * dilog(-(yd / w5)) + 4.0 * dilog(yd / w5) - 4.0 * dilog(-(yd / w7)) - 4.0 * dilog(yd / w7) + 4.0 * dilog((1.0 + yd) / 2.0) + 4.0 * dilog((w4 + yd) / (1.0 + w4)) + 4.0 * dilog((w5 + yd) / (-1.0 + w5)) - 4.0 * dilog((w7 + yd) / (1.0 + w7))
                    - 4.0 * ln2 * ln1myd - 4.0 * ln2 * ln1pyd + 4.0 * log((1.0 - yd) / (1.0 + w4)) * log((w4 + yd) / w4) + 4.0 * log((1.0 + yd) / (1.0 - w5)) * log((w5 + yd) / w5) - 4.0 * log((1.0 - yd) / (1.0 + w7)) * log((w7 + yd) / w7)
                    + 4.0 * log((1.0 + yd) / (1.0 + w4)) * log(1.0 - yd / w4) + 4.0 * log((-1.0 + yd) / (-1.0 + w5)) * log(1.0 - yd / w5) - 4.0 * log((1.0 + yd) / (1.0 + w7)) * log(1.0 - yd / w7)
                    + (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) + (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                    - (8.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv) + (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv)
                    + (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) - (8.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                + li2half * (16.0 * dilog(w4 / (1.0 + w4)) + 16.0 * dilog(w5 / (-1.0 + w5)) - 16.0 * dilog(w7 / (1.0 + w7)) + 8.0 * dilog(-1.0i / (-1.0i + wx3)) + 8.0 * dilog(1.0i / (1.0i + wx3)) + 8.0 * dilog(-1.0i / (-1.0i + wx4)) + 8.0 * dilog(1.0i / (1.0i + wx4))
                    - 4.0 * dilog((1.0i - xd) / (1.0i + wx3)) - 4.0 * dilog((1.0i - xd) / (1.0i + wx4)) - 8.0 * dilog(1.0 / 2.0 + (1.0i / 2.0) * xd) - 8.0 * dilog((-1.0i) * xd) - 8.0 * dilog(1.0i * xd) + 4.0 * dilog(-(xd / wx3)) + 4.0 * dilog(xd / wx3) + 4.0 * dilog(-(xd / wx4))
                    + 4.0 * dilog(xd / wx4) - 4.0 * dilog((-1.0i + xd) / (-1.0i + wx3)) - 4.0 * dilog((-1.0i + xd) / (-1.0i + wx4)) - 8.0 * dilog(((-1.0 / 2.0) * 1.0i) * (1.0i + xd)) - 4.0 * dilog((1.0i + xd) / (1.0i - wx3)) - 4.0 * dilog((1.0i + xd) / (1.0i + wx3))
                    - 4.0 * dilog((1.0i + xd) / (1.0i - wx4)) - 4.0 * dilog((1.0i + xd) / (1.0i + wx4)) - 8.0 * dilog((1.0 - yd) / 2.0) - 8.0 * dilog((w4 - yd) / (1.0 + w4)) - 8.0 * dilog((w5 - yd) / (-1.0 + w5)) + 8.0 * dilog((w7 - yd) / (1.0 + w7)) + 8.0 * dilog(-yd)
                    + 8.0 * dilog(yd) - 8.0 * dilog(-(yd / w4)) - 8.0 * dilog(yd / w4) - 8.0 * dilog(-(yd / w5)) - 8.0 * dilog(yd / w5) + 8.0 * dilog(-(yd / w7)) + 8.0 * dilog(yd / w7) - 8.0 * dilog((1.0 + yd) / 2.0) - 8.0 * dilog((w4 + yd) / (1.0 + w4))
                    - 8.0 * dilog((w5 + yd) / (-1.0 + w5)) + 8.0 * dilog((w7 + yd) / (1.0 + w7)) + 4.0 * power_of<2>(log(1.0 - 1.0i * xd)) + 4.0 * power_of<2>(log(1.0 + 1.0i * xd)) + log(1.0 - 1.0i * xd) * (8.0 * ln2 - 4.0 * log((wx3 - xd) / (1.0i + wx3)) - 4.0 * log((wx4 - xd) / (1.0i + wx4))
                    - 8.0 * log(xd) - 4.0 * log((wx3 + xd) / (-1.0i + wx3)) - 4.0 * log((wx4 + xd) / (-1.0i + wx4))) + log(1.0 + 1.0i * xd) * (8.0 * ln2 - 4.0 * log((wx3 - xd) / (-1.0i + wx3)) - 4.0 * log((wx4 - xd) / (-1.0i + wx4)) - 8.0 * log(xd)
                    - 4.0 * log((wx3 + xd) / (1.0i + wx3)) - 4.0 * log((wx4 + xd) / (1.0i + wx4))) + log(xd) * (4.0 * log((wx3 + xd) / wx3) + 4.0 * log((wx4 + xd) / wx4) + 4.0 * log(1.0 - xd / wx3) + 4.0 * log(1.0 - xd / wx4)) + 8.0 * ln2 * ln1myd
                    + 8.0 * ln2 * ln1pyd - 8.0 * log((1.0 - yd) / (1.0 + w4)) * log((w4 + yd) / w4) - 8.0 * log((1.0 + yd) / (1.0 - w5)) * log((w5 + yd) / w5) + 8.0 * log((1.0 - yd) / (1.0 + w7)) * log((w7 + yd) / w7)
                    - 8.0 * log((1.0 + yd) / (1.0 + w4)) * log(1.0 - yd / w4) - 8.0 * log((-1.0 + yd) / (-1.0 + w5)) * log(1.0 - yd / w5) + 8.0 * log((1.0 + yd) / (1.0 + w7)) * log(1.0 - yd / w7)
                    - (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx3) - (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx3)
                    - (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 - 1.0i * wx4) - (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 - 1.0i * xd, 1.0 + 1.0i * wx4)
                    - (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx3) * my_sign(-real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx3) - (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx3) * my_sign(real(wx3inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx3)
                    - (8.0 * 1.0i) * M_PI * log(1.0 - 1.0i * wx4) * my_sign(-real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 - 1.0i * wx4) - (8.0 * 1.0i) * M_PI * log(1.0 + 1.0i * wx4) * my_sign(real(wx4inv)) * T(1.0, 1.0 + 1.0i * xd, 1.0 + 1.0i * wx4)
                    - (16.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) - (16.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                    + (16.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv) - (16.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv)
                    - (16.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) + (16.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                + log(xd) * ((log(1.0 - xd / wx3) + log(1.0 - xd / wx4)) * ((-8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) + (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                    - (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv) + (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                    + (log((wx3 + xd) / wx3) + log((wx4 + xd) / wx4)) * ((-8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5) + (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv)
                    - (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) + (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv)))
                + log(1.0 - 1.0i * xd) * ((-8.0 * 1.0i) * M_PI * ln2 * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) - (8.0 * 1.0i) * M_PI * ln2 * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                    + (8.0 * 1.0i) * M_PI * ln2 * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv) - (8.0 * 1.0i) * M_PI * ln2 * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv)
                    - (8.0 * 1.0i) * M_PI * ln2 * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) + (8.0 * 1.0i) * M_PI * ln2 * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv)
                    + (log((wx3 - xd) / (1.0i + wx3)) + log((wx4 - xd) / (1.0i + wx4))) * ((8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) *
                    T(1.0, (w7 + yd) / w7, 1.0 + w7inv) + (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv) - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                    + (log((wx3 + xd) / (-1.0i + wx3)) + log((wx4 + xd) / (-1.0i + wx4))) * ((8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5) - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) *
                    T(1.0, (w7 + yd) / w7, 1.0 + w7inv) + (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                    + log(xd) * ((8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) + (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                    - (8.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv) + (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv)
                    + (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) - (8.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv)))
                + log(1.0 + 1.0i * xd) * ((-8.0 * 1.0i) * M_PI * ln2 * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) - (8.0 * 1.0i) * M_PI * ln2 * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                    + (8.0 * 1.0i) * M_PI * ln2 * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv) - (8.0 * 1.0i) * M_PI * ln2 * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv)
                    - (8.0 * 1.0i) * M_PI * ln2 * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) + (8.0 * 1.0i) * M_PI * ln2 * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv)
                    + (log((wx3 - xd) / (-1.0i + wx3)) + log((wx4 - xd) / (-1.0i + wx4))) * ((8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) *
                    T(1.0, (w7 + yd) / w7, 1.0 + w7inv) + (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv) - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                    + (log((wx3 + xd) / (1.0i + wx3)) + log((wx4 + xd) / (1.0i + wx4))) * ((8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5) - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) *
                    T(1.0, (w7 + yd) / w7, 1.0 + w7inv) + (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) - (4.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))
                    + log(xd) * ((8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, (w4 + yd) / w4, 1.0 + w4inv) + (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, (w5 + yd) / w5, (-1.0 + w5) / w5)
                    - (8.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, (w7 + yd) / w7, 1.0 + w7inv) + (8.0 * 1.0i) * M_PI * log(1.0 + w4inv) * my_sign(-imag(w4)) * T(1.0, 1.0 - yd / w4, 1.0 + w4inv)
                    + (8.0 * 1.0i) * M_PI * log((-1.0 + w5) / w5) * my_sign(imag(w5)) * T(1.0, 1.0 - yd / w5, (-1.0 + w5) / w5) - (8.0 * 1.0i) * M_PI * log(1.0 + w7inv) * my_sign(-imag(w7)) * T(1.0, 1.0 - yd / w7, 1.0 + w7inv))))) / denom2;

            const complex<double> f29dlogs1 = term1 + (64.0 * power_of<2>(lnmuhat)) / 9.0 + (logs17 * log(1.0 - 1.0i * xd) + logs16 * log(1.0 + 1.0i * xd)) / denom6
                + ((-logs19 - logs20) * li2half + logs20 * dilog(1.0 / 2.0 - (1.0i / 2.0) * xd) + logs19 * dilog((1.0 + 1.0i * xd) / 2.0) - logs20 * ln2 * log(1.0 - 1.0i * xd) - (logs27 * power_of<2>(log(1.0 - 1.0i * xd))) / 2.0
                - logs19 * ln2 * log(1.0 + 1.0i * xd) - (logs26 * power_of<2>(log(1.0 + 1.0i * xd))) / 2.0) / denom3 + (logs15 * log(xd)) / denom5
                + (-(logs24 * li2half) + logs24 * dilog((1.0 - yd) / 2.0) - logs24 * ln2 * ln1myd + (logs24 * power_of<2>(ln1myd)) / 2.0) / denom7
                + factor1 * (-4.0 * yd * lnmuhat + 9.0 * yd3 * lnmuhat + 16.0 * yd3 * ln2 * lnmuhat - 8.0 * yd3 * lnmuhat * log(1.0 - 1.0i * xd) - 8.0 * yd3 * lnmuhat * log(1.0 + 1.0i * xd)
                + 8.0 * yd3 * lnmuhat * log(xd) + (-2.0 * lnmuhat + 6.0 * yd2 * lnmuhat) * ln1myd + (2.0 * lnmuhat - 6.0 * yd2 * lnmuhat) * ln1pyd)
                + (-(logs10 * dilog(-yd)) + logs10 * dilog(yd) + (logs23 - logs13 * log(1.0 - 1.0i * xd) - logs13 * log(1.0 + 1.0i * xd)) * ln1myd
                + (-logs23 - logs14 * log(1.0 - 1.0i * xd) - logs14 * log(1.0 + 1.0i * xd)) * ln1pyd) / denom1
                + (-(logs25 * li2half) + logs25 * dilog((1.0 + yd) / 2.0) - logs25 * ln2 * ln1pyd + (logs25 * power_of<2>(ln1pyd)) / 2.0) / denom8
                + (logs18 * dilog((-1.0i) * xd) + logs18 * dilog(1.0i * xd) - logs21 * trilog(1.0 - yd) - logs21 * trilog(yd) - logs21 * trilog(yd / (1.0 + yd)) + logs21 * trilog((2.0 * yd) / (1.0 + yd))
                - logs21 * trilog((1.0 + yd) / 2.0) - (logs21 * pisqu * ln2) / 12.0 + (logs21 * power_of<3>(ln2)) / 6.0 + dilog((1.0 + yd) / 2.0) * (logs1 * log(1.0 - 1.0i * xd) + logs1 * log(1.0 + 1.0i * xd))
                + dilog(-yd) * (-(logs12 * log(1.0 - 1.0i * xd)) - logs12 * log(1.0 + 1.0i * xd)) + li2half * ((-logs1 - logs2) * log(1.0 - 1.0i * xd) + (-logs1 - logs2) * log(1.0 + 1.0i * xd))
                + dilog((1.0 - yd) / 2.0) * (logs2 * log(1.0 - 1.0i * xd) + logs2 * log(1.0 + 1.0i * xd)) + logs21 * dilog(1.0 - yd) * ln1myd
                + (-(logs2 * ln2 * log(1.0 - 1.0i * xd)) - logs2 * ln2 * log(1.0 + 1.0i * xd)) * ln1myd + (logs21 * power_of<2>(ln1myd) * log(yd)) / 2.0
                + ((logs21 * pisqu) / 12.0 - (logs21 * ln2squ) / 2.0 - logs1 * ln2 * log(1.0 - 1.0i * xd) - logs1 * ln2 * log(1.0 + 1.0i * xd)) * ln1pyd + (logs21 * ln2 * power_of<2>(ln1pyd)) / 2.0
                + dilog(yd) * (-(logs11 * log(1.0 - 1.0i * xd)) - logs11 * log(1.0 + 1.0i * xd) + logs21 * ln1pyd) + (15.0 * logs21 * zeta3) / 8.0) / denom2;

            const complex<double> f29dlogs2 = ((logs3 * pisqu * ln4) / 12.0 + (-2.0 * logs6 + logs8 + logs9) * li3half + (logs3 - logs6 + logs7) * trilog(1.0 / 2.0 - (1.0i / 2.0) * xd) - logs7 * trilog(1.0 - 1.0i * xd)
                + (logs3 - logs6 + logs7) * trilog((1.0 + 1.0i * xd) / 2.0) - logs7 * trilog(1.0 + 1.0i * xd) + logs7 * trilog((-1.0i) * xd) + logs7 * trilog(1.0i * xd) + logs7 * trilog(xd / (-1.0i + xd))
                - logs7 * trilog((2.0 * xd) / (-1.0i + xd)) + logs7 * trilog(xd / (1.0i + xd)) - logs7 * trilog((2.0 * xd) / (1.0i + xd)) + (-logs22 + logs4 + logs5 - logs8 - 2.0 * logs9) * trilog((1.0 - yd) / 2.0)
                + (-logs4 + 2.0 * logs5) * trilog(1.0 - yd) + (-logs22 + logs4) * trilog(-yd) + logs5 * trilog(yd) + (-logs22 + logs5) * trilog(yd / (-1.0 + yd))
                + (logs22 - logs4 - logs5) * trilog((2.0 * yd) / (-1.0 + yd)) + logs4 * trilog(yd / (1.0 + yd)) + (-logs4 - logs5) * trilog((2.0 * yd) / (1.0 + yd))
                + (logs4 + logs5 - 2.0 * logs8 - logs9) * trilog((1.0 + yd) / 2.0) + (-logs22 + 2.0 * logs4 - logs5) * trilog(1.0 + yd) - (logs22 * pisqu * ln2) / 12.0 - (logs6 * pisqu * ln2) / 3.0
                + (logs7 * pisqu * ln2) / 6.0 - (logs8 * pisqu * ln2) / 6.0 - (logs9 * pisqu * ln2) / 6.0 + (logs22 * power_of<3>(ln2)) / 6.0 - (logs3 * power_of<3>(ln2)) / 3.0 - (logs4 * power_of<3>(ln2)) / 3.0
                - (logs5 * power_of<3>(ln2)) / 3.0 + (2.0 * logs6 * power_of<3>(ln2)) / 3.0 - (logs7 * power_of<3>(ln2)) / 3.0 + (logs8 * power_of<3>(ln2)) / 3.0 + (logs9 * power_of<3>(ln2)) / 3.0 + (logs4 * pisqu * ln4) / 12.0
                + (logs5 * pisqu * ln4) / 12.0 + logs7 * dilog(1.0 - 1.0i * xd) * log(1.0 - 1.0i * xd) - logs7 * dilog((-1.0i) * xd) * log(1.0 - 1.0i * xd) - (logs3 * power_of<3>(log(1.0 - 1.0i * xd))) / 6.0
                + logs7 * dilog(1.0 + 1.0i * xd) * log(1.0 + 1.0i * xd) - logs7 * dilog(1.0i * xd) * log(1.0 + 1.0i * xd) + ((logs3 * pisqu) / 12.0 - (logs6 * pisqu) / 12.0 - (logs7 * pisqu) / 12.0 - (logs3 * ln2squ) / 2.0
                + (logs6 * ln2squ) / 2.0 + (logs7 * ln2squ) / 2.0) * log(1.0 + 1.0i * xd) - (logs3 * power_of<3>(log(1.0 + 1.0i * xd))) / 6.0
                + dilog(1.0 / 2.0 - (1.0i / 2.0) * xd) * ((-logs3 + logs6) * log(1.0 - 1.0i * xd) + (-logs3 + logs6) * log(1.0 + 1.0i * xd))
                + dilog((1.0 + 1.0i * xd) / 2.0) * ((-logs3 + logs6) * log(1.0 - 1.0i * xd) + (-logs3 + logs6) * log(1.0 + 1.0i * xd))
                + log(1.0 - 1.0i * xd) * ((logs3 * pisqu) / 12.0 - (logs6 * pisqu) / 12.0 - (logs7 * pisqu) / 12.0 - (logs3 * ln2squ) / 2.0 + (logs6 * ln2squ) / 2.0 + (logs7 * ln2squ) / 2.0
                + (-2.0 * logs6 * ln2 - logs3 * log(1.0 / 2.0 - (1.0i / 2.0) * xd)) * log(1.0 + 1.0i * xd) - logs3 * log((1.0 + 1.0i * xd) / 2.0) * log(1.0 + 1.0i * xd))
                + power_of<2>(log(1.0 + 1.0i * xd)) * ((-1.0 / 2.0) * (logs7 * ln2) + (logs3 * ln64) / 12.0 + (logs6 * log(1.0 / 2.0 - (1.0i / 2.0) * xd)) / 2.0 + (logs7 * log((-1.0i) * xd)) / 2.0)
                + power_of<2>(log(1.0 - 1.0i * xd)) * ((-1.0 / 2.0) * (logs7 * ln2) + (logs3 * ln64) / 12.0 + (logs6 * log((1.0 + 1.0i * xd) / 2.0)) / 2.0 + (logs7 * log(1.0i * xd)) / 2.0) - 2.0 * logs5 * dilog(1.0 - yd) * ln1myd
                + power_of<2>(ln1myd) * ((logs22 * ln2) / 2.0 - (logs5 * ln64) / 12.0 + (logs4 * log((1.0 - yd) / 8.0)) / 6.0 + ((-1.0 / 2.0) * logs4 - logs5) * log(yd) + (logs8 / 2.0 + logs9) * log((1.0 + yd) / 2.0))
                + (logs22 - 2.0 * logs4) * dilog(1.0 + yd) * ln1pyd + ((-1.0 / 12.0) * (logs4 * pisqu) + (logs5 * pisqu) / 12.0 - (logs8 * pisqu) / 12.0 + (logs4 * ln2squ) / 2.0 + (logs5 * ln2squ) / 2.0
                + (logs8 * ln2squ) / 2.0) * ln1pyd + ((-1.0 / 12.0) * (logs4 * ln64) + (logs8 + logs9 / 2.0) * log((1.0 - yd) / 2.0) + (logs22 / 2.0 - logs4 - logs5 / 2.0) * log(-yd)
                + (logs5 * log((1.0 + yd) / 8.0)) / 6.0) * power_of<2>(ln1pyd) + dilog(-yd) * ((logs22 - logs4) * ln1myd - logs4 * ln1pyd)
                + dilog(yd) * (-(logs5 * ln1myd) - logs5 * ln1pyd) + dilog((1.0 - yd) / 2.0) * ((logs8 + 2.0 * logs9) * ln1myd + logs8 * ln1pyd)
                + dilog((1.0 + yd) / 2.0) * (logs9 * ln1myd + (2.0 * logs8 + logs9) * ln1pyd) + ln1myd * ((logs22 * pisqu) / 12.0 + (logs4 * pisqu) / 12.0 - (logs5 * pisqu) / 12.0 - (logs9 * pisqu) / 12.0
                - (logs22 * ln2squ) / 2.0 + (logs4 * ln2squ) / 2.0 + (logs5 * ln2squ) / 2.0 + (logs9 * ln2squ) / 2.0 + (-(logs8 * ln2) - logs9 * ln2) * ln1pyd)
                + (15.0 * logs22 * zeta3) / 8.0 - (7.0 * logs3 * zeta3) / 4.0 - (11.0 * logs4 * zeta3) / 4.0 - (11.0 * logs5 * zeta3) / 4.0 + (7.0 * logs6 * zeta3) / 2.0 + (logs7 * zeta3) / 4.0 + (7.0 * logs8 * zeta3) / 4.0
                + (7.0 * logs9 * zeta3) / 4.0) / denom2;

            const complex<double> f29dPart6 = (num20 / denom2) * f29d_part6(clp);

            const complex<double> f29dPart78 = (num20 / denom2) * f279d_log2_terms(clp);

            const complex<double> f29dPart9 = (num20 / denom2) * GPLweight4Parts(clp);

            return f29dPart1 + f29dPart2 + f29dPart3 + f29dPart41 + f29dPart42 + f29dPart43 + f29dPart44 + f29dPart45 + f29dPart5 + f29dPart6 + f29dPart78 + f29dPart9 + f29dlogs1 + f29dlogs2;
        }
    }
}
