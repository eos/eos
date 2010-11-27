/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_QCD_HH
#define EOS_GUARD_SRC_UTILS_QCD_HH 1

#include <array>

namespace eos
{
    class QCD
    {
        public:
            typedef std::array<double, 4> AnomalousMassDimension;
            typedef std::array<double, 4> BetaFunction;
            struct Parameters;

            /*
             * Running coupling alpha_s down to mu from mu_0, based on QCD parameters in params.
             *
             * Calculation according to [CKS2000].
             *
             * @mu        : scale at which alpha_s shall be evaluated
             * @alpha_s_0 : alpha_s at the initial scale mu_0
             * @mu_0      : initial scale mu_0
             * @beta      : parameters of QCD beta function that control the running
             */
            static double alpha_s(const double & mu, const double & alpha_s_0, const double & mu_0, const BetaFunction & beta);

            /*
             * Running quark mass in the MSbar scheme.
             *
             * Calculation according to [CKS2000].
             *
             * @m_q_0       : quark mass in the MSbar scheme at the scale m_q
             * @alpha_s_0   : alpha_s at the scale mu_0 = m_q
             * @alpha_s_mu  : alpha_s at the scale mu
             * @beta        : parameters of QCD beta function that control the running
             * @gamma_m     : parameters of QCD anomalous mass dimension that control the running
             */
            static double m_q_msbar(const double & m_q, const double & alpha_s_mu0, const double & alpha_s_mu, const BetaFunction & beta, const AnomalousMassDimension & gamma_m);

            /*
             * Pole mass shift from the MSbar scheme
             *
             * Calculation according to [CERN2002-002].
             *
             * @m_q     : Quark mass in the MSbar scheme at the scale m_q
             * @alpha_s : alpha_s at the scale m_q
             * @params  : QCD parameters that control the calculation
             */
            static double m_q_pole(const double & m_q, const double & alpha_s, const Parameters & params);

            /*
             * Potential-Subtracted mass (PS mass) shift from the MSbar scheme
             *
             * Calculation according to [B1998].
             *
             * @m_q     : Quark mass in the MSbar scheme at the scale m_q
             * @alpha_s : alpha_s at the scalem_q
             * @mu_f    : Factorization scale
             * @params  : QCD parameters that control the calculation
             */
            static double m_q_ps(const double & m_q, const double & alpha_s, const double & mu_f, const Parameters & params);

            /*
             * The quadratic casimir operator for the fundamental representation
             */
            static const double casimir_f;

            /*
             * The quadratic casimir operator for the ajoint representation
             */
            static const double casimir_a;
    };

    struct QCD::Parameters
    {
        /* Number of active flavors */
        double _nf;

        QCD::Parameters
        nf(const double & v)
        {
            _nf = v;
            return *this;
        }

        /* Coefficients of the QCD potential */
        double _a1;
        double _a2;

        QCD::Parameters
        a1(const double & v)
        {
            _a1 = v;
            return *this;
        }

        QCD::Parameters
        a2(const double & v)
        {
            _a2 = v;
            return *this;
        }

        /* Coefficients of the QCD beta function */
        double _beta0;
        double _beta1;
        double _beta2;
        double _beta3;

        QCD::Parameters
        beta0(const double & v)
        {
            _beta0 = v;
            return *this;
        }

        QCD::Parameters
        beta1(const double & v)
        {
            _beta1 = v;
            return *this;
        }

        QCD::Parameters
        beta2(const double & v)
        {
            _beta2 = v;
            return *this;
        }

        QCD::Parameters
        beta3(const double & v)
        {
            _beta3 = v;
            return *this;
        }

        /* Coefficients of the anomalous dimension of quark masses (MSbar scheme) */
        double _gamma0_m;
        double _gamma1_m;
        double _gamma2_m;
        double _gamma3_m;

        QCD::Parameters
        gamma0_m(const double & v)
        {
            _gamma0_m = v;
            return *this;
        }

        QCD::Parameters
        gamma1_m(const double & v)
        {
            _gamma1_m = v;
            return *this;
        }

        QCD::Parameters
        gamma2_m(const double & v)
        {
            _gamma2_m = v;
            return *this;
        }

        QCD::Parameters
        gamma3_m(const double & v)
        {
            _gamma3_m = v;
            return *this;
        }
    };
}

#endif
