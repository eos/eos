/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_TOP_LOOPS_HH
#define WFITTER_GUARD_SRC_UTILS_TOP_LOOPS_HH 1

namespace wf
{
    struct TopLoops
    {
        /* One-Loop Functions
         *
         * MSbar scheme
         *
         * x_t: m_t(mu_t)^2 / m_W^2
         */
        static double A0(const double & x_t);

        static double B0(const double & x_t);

        static double C0(const double & x_t);

        static double D0(const double & x_t);

        static double E0(const double & x_t);

        static double F0(const double & x_t);

        /*
         * Two-Loop Functions
         *
         * MSbar scheme
         *
         * x_t:   m_t(mu_t)^2 / m_W^2
         * log_t: ln(mu_t / m_t(mu_t))
         */
        static double A1(const double & x_t, const double & log_t);

        static double B1(const double & x_t, const double & log_t);

        static double C1(const double & x_t, const double & log_t);

        static double D1(const double & x_t, const double & log_t);

        static double E1(const double & x_t, const double & log_t);

        static double F1(const double & x_t, const double & log_t);

        static double G1(const double & x_t, const double & log_t);
    };
}

#endif
