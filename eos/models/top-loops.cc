/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
 * Copyright (c) 2010 Christoph Bobeth
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
#include <eos/models/top-loops.hh>

#include <gsl/gsl_sf_dilog.h>

#include <cmath>

namespace eos
{
    double
    TopLoops::A0(const double & x_t)
    {
        return (-22.0 * power_of<3>(x_t) + 153.0 * power_of<2>(x_t) - 159.0 * x_t + 46.0) / (36.0 * power_of<3>(x_t - 1.0))
               + power_of<2>(x_t) * (2.0 - 3.0 * x_t) * log(x_t) / (2.0 * power_of<4>(x_t - 1.0));
    }

    double
    TopLoops::B0(const double & x_t)
    {
        return -1.0 / (4.0 * (x_t - 1.0)) + x_t / (4.0 * power_of<2>(x_t - 1.0)) * log(x_t);
    }

    double
    TopLoops::C0(const double & x_t)
    {
        return x_t * (x_t - 6) / (8 * (x_t - 1.0)) + x_t * (3.0 * x_t + 2.0) / (8.0 * power_of<2>(x_t - 1.0)) * log(x_t);
    }

    double
    TopLoops::D0(const double & x_t)
    {
        return (47.0 * power_of<3>(x_t) - 237.0 * power_of<2>(x_t) + 312.0 * x_t - 104.0) / (108.0 * power_of<3>(x_t - 1.0))
               + (-3.0 * power_of<4>(x_t) + 30.0 * power_of<3>(x_t) - 54.0 * power_of<2>(x_t) + 32.0 * x_t - 8.0) / (18.0 * power_of<4>(x_t - 1.0)) * log(x_t);
    }

    double
    TopLoops::E0(const double & x_t)
    {
        return (7.0 * power_of<3>(x_t) + 21.0 * power_of<2>(x_t) - 42.0 * x_t - 4.0) / (36.0 * power_of<3>(x_t - 1.0))
               + (-4.0 + 16.0 * x_t - 9.0 * power_of<2>(x_t)) * log(x_t) / (6.0 * power_of<4>(x_t - 1.0));
    }

    double
    TopLoops::F0(const double & x_t)
    {
        return (8.0 - 30.0 * x_t + 9.0 * power_of<2>(x_t) - 5.0 * power_of<3>(x_t)) / (12.0 * power_of<3>(x_t - 1.0))
               + 3.0 * power_of<2>(x_t) / (2.0 * power_of<4>(x_t - 1.0)) * log(x_t);
    }

    double
    TopLoops::A1(const double & x_t, const double & log_t)
    {
        return (32 * power_of<4>(x_t) + 244 * power_of<3>(x_t) - 160 * x_t * x_t + 16 * x_t) / (9 * power_of<4>(x_t - 1.0)) * gsl_sf_dilog(1.0 - 1.0 / x_t)
               + (774 * power_of<4>(x_t) + 2826 * power_of<3>(x_t) - 1994 * x_t * x_t + 130 * x_t - 8) / (81 * power_of<5>(x_t - 1.0)) * log(x_t)
               + (-94 * power_of<4>(x_t) - 18665 * power_of<3>(x_t) + 20682 * x_t * x_t - 9113 * x_t + 2006) / (243 * power_of<4>(x_t - 1.0))
               + 2.0 * log_t
                         * ((12 * power_of<4>(x_t) + 92 * power_of<3>(x_t) - 56 * x_t * x_t) / (3 * power_of<5>(x_t - 1.0)) * log(x_t)
                            + (-68 * power_of<4>(x_t) - 202 * power_of<3>(x_t) - 804 * x_t * x_t + 794 * x_t - 152) / (27 * power_of<4>(x_t - 1.0)));
    }

    double
    TopLoops::B1(const double & x_t, const double & log_t)
    {
        return (-2.0 * x_t) / power_of<2>(x_t - 1.0) * gsl_sf_dilog(1.0 - 1.0 / x_t) + (x_t * x_t - 17 * x_t) / (3 * power_of<3>(x_t - 1.0)) * log(x_t)
               + (13 * x_t + 3) / (3 * power_of<2>(x_t - 1.0)) + 2.0 * log_t * ((-2 * x_t * x_t - 2 * x_t) / power_of<3>(x_t - 1.0) * log(x_t) + 4 * x_t / power_of<2>(x_t - 1.0));
    }

    double
    TopLoops::C1(const double & x_t, const double & log_t)
    {
        return (-power_of<3>(x_t) - 4 * x_t) / power_of<2>(x_t - 1.0) * gsl_sf_dilog(1.0 - 1.0 / x_t)
               + (-3 * power_of<3>(x_t) - 14 * x_t * x_t - 23 * x_t) / (3 * power_of<3>(x_t - 1.0)) * log(x_t)
               + (4 * power_of<3>(x_t) + 7 * x_t * x_t + 29 * x_t) / (3 * power_of<2>(x_t - 1.0))
               + 2.0 * log_t * ((-8 * x_t * x_t - 2 * x_t) / power_of<3>(x_t - 1.0) * log(x_t) + (power_of<3>(x_t) + x_t * x_t + 8 * x_t) / power_of<2>(x_t - 1.0));
    }

    double
    TopLoops::D1(const double & x_t, const double & log_t)
    {
        return (380 * power_of<4>(x_t) - 1352 * power_of<3>(x_t) + 1656 * x_t * x_t - 784 * x_t + 256) / (81 * power_of<4>(x_t - 1.0)) * gsl_sf_dilog(1.0 - 1.0 / x_t)
               + (-304 * power_of<4>(x_t) - 1716 * power_of<3>(x_t) + 4644 * x_t * x_t - 2768 * x_t + 720) / (81 * power_of<5>(x_t - 1.0)) * log(x_t)
               + (-6175 * power_of<4>(x_t) + 41608 * power_of<3>(x_t) - 66723 * x_t * x_t + 33106 * x_t - 7000) / (729 * power_of<4>(x_t - 1.0))
               + 2.0 * log_t
                         * ((-648 * power_of<4>(x_t) + 720 * power_of<3>(x_t) + 232 * x_t * x_t + 160 * x_t - 32) / (81 * power_of<5>(x_t - 1.0)) * log(x_t)
                            + (-352 * power_of<4>(x_t) + 4912 * power_of<3>(x_t) - 8280 * x_t * x_t + 3304 * x_t - 880) / (243 * power_of<4>(x_t - 1.0)));
    }

    double
    TopLoops::E1(const double & x_t, const double & log_t)
    {
        return (515 * power_of<4>(x_t) - 614 * power_of<3>(x_t) - 81 * x_t * x_t - 190 * x_t + 40) / (54 * power_of<4>(x_t - 1.0)) * gsl_sf_dilog(1.0 - 1.0 / x_t)
               + (1030 * power_of<4>(x_t) - 435 * power_of<3>(x_t) - 1373 * x_t * x_t - 1950 * x_t + 424) / (108 * power_of<5>(x_t - 1.0)) * log(x_t)
               + (-29467 * power_of<4>(x_t) + 45604 * power_of<3>(x_t) - 30237 * x_t * x_t + 66532 * x_t - 10960) / (1944 * power_of<4>(x_t - 1.0))
               + 2 * log_t
                         * ((1125 * power_of<3>(x_t) - 1685 * x_t * x_t - 380 * x_t + 76) / (54 * power_of<5>(x_t - 1.0)) * log(x_t)
                            + (133 * power_of<4>(x_t) - 2758 * power_of<3>(x_t) - 2061 * x_t * x_t + 11522 * x_t - 1652) / (324 * power_of<4>(x_t - 1.0)));
    }

    double
    TopLoops::F1(const double & x_t, const double & log_t)
    {
        return (4.0 * power_of<4>(x_t) - 40 * power_of<3>(x_t) - 41 * x_t * x_t - x_t) / (3 * power_of<4>(x_t - 1.0)) * gsl_sf_dilog(1.0 - 1.0 / x_t)
               + (144.0 * power_of<4>(x_t) - 3177 * power_of<3>(x_t) - 3661 * x_t * x_t - 250 * x_t + 32) / (108 * power_of<5>(x_t - 1.0)) * log(x_t)
               + (-247.0 * power_of<4>(x_t) + 11890 * power_of<3>(x_t) + 31779 * x_t * x_t - 2966 * x_t + 1016) / (648 * power_of<4>(x_t - 1.0))
               + 2.0 * log_t
                         * ((-17 * power_of<3>(x_t) - 31 * x_t * x_t) / power_of<5>(x_t - 1.0) * log(x_t)
                            + (-35.0 * power_of<4>(x_t) + 170 * power_of<3>(x_t) + 447 * x_t * x_t + 338 * x_t - 56) / (18 * power_of<4>(x_t - 1.0)));
    }

    double
    TopLoops::G1(const double & x_t, const double & log_t)
    {
        return (-42.0 - 161 * x_t + 293.0 * x_t * x_t + 6.0 * power_of<3>(x_t)) / (81.0 * power_of<3>(x_t - 1.0))
               + (68.0 - 332 * x_t - 42.0 * x_t * x_t + 30.0 * power_of<3>(x_t)) * log(x_t) / (81.0 * power_of<4>(x_t - 1.0))
               + 2.0 * log_t
                         * ((90.0 * x_t * x_t - 160.0 * x_t + 40) * log(x_t) / (27.0 * power_of<4>(x_t - 1.0))
                            + (-35.0 * power_of<3>(x_t) - 105.0 * x_t * x_t + 210.0 * x_t + 20.0) / (81.0 * power_of<3>(x_t - 1.0)))
               + (10.0 * (-4.0 + 16.0 * x_t + 3.0 * x_t * x_t - 10.0 * power_of<3>(x_t) + power_of<4>(x_t)) * gsl_sf_dilog(1.0 - 1.0 / x_t)) / (27.0 * power_of<4>(x_t - 1.0));
    }
} // namespace eos
