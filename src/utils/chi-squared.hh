/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_CHI_SQUARED_HH
#define EOS_GUARD_SRC_UTILS_CHI_SQUARED_HH 1

namespace eos
{
    struct ChiSquared
    {
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
    };
}

#endif
