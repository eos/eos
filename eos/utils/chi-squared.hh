/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_CHI_SQUARED_HH
#define EOS_GUARD_SRC_UTILS_CHI_SQUARED_HH 1

#include <eos/utils/exception.hh>

#include <functional>

namespace eos
{
    struct ChiSquared
    {
        /*
         * Signature for Chi-Squared functions.
         */

        typedef std::function<double (double, double, double, double, double, double)> Function;

        /*
         * Chi-Squared function with theory offset.
         *
         * Calculation according to [BHvD2010], Eq. (4.20)
         *
         * @theory_min:         Minimal value of theory prediction.
         * @theory_central:     Central value of theory prediction.
         * @theory_max:         Maximal value of theory prediction.
         * @experiment_min:     Minimal value of experimental measurement.
         * @experiment_central: Central value of experimental measurement.
         * @experiment_max:     Maximal value of experimental measurement.
         */
        static double with_theory_offset(const double & theory_min, const double & theory_central, const double & theory_max,
                const double & experiment_min, const double & experiment_central, const double & experiment_max);

        /*
         * Chi-Squared function with combined uncertainty from theory and experiment.
         *
         *   chi = (th_cen - exp_cen) / sqrt(sigma^2 + delta_\pm^2)
         *
         * @theory_min:         Minimal value of theory prediction.
         * @theory_central:     Central value of theory prediction.
         * @theory_max:         Maximal value of theory prediction.
         * @experiment_min:     Minimal value of experimental measurement.
         * @experiment_central: Central value of experimental measurement.
         * @experiment_max:     Maximal value of experimental measurement.
         */
        static double with_combined_uncertainties(const double & theory_min, const double & theory_central, const double & theory_max,
                const double & experiment_min, const double & experiment_central, const double & experiment_max);

        /*
         * Factory method to create a ChiSquared::Function from a given name
         */
        static Function make(const std::string & name);
    };

    struct NoSuchChiSquaredError :
        public Exception
    {
        NoSuchChiSquaredError(const std::string & name);
    };
}

#endif
