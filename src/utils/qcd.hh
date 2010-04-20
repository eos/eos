/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_QCD_HH
#define WFITTER_GUARD_SRC_UTILS_QCD_HH 1

namespace wf
{
    class QCD
    {
        public:
            // \alpha_s for n_f = 5
            static double alpha_s(const double & mu);

            // The quadratic casimir operator
            static const double casimir_f = 4.0 / 3.0;

            // Potential-subtracted b mass m_b^PS
            static double mb_PS(const double & mb_MSbar, const double & mu, const double & mu_PS);

            // Pole mass of the b
            static double mb_pole(const double & mb_MSbar);

            // Running MSbar mass of the b
            static double mb_MSbar(const double & mb_MSbar0, const double & mu);

            // Running MSbar mass of the c
            static double mc_MSbar(const double & mb_MSbar0, const double & mu);
    };
}

#endif
