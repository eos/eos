/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/form_factors.hh>

#include <cmath>

namespace wf
{
    /* Form Factors according to [ABHH1999] */
    class ABHH1999FormFactors :
        public FormFactors<BToKstar>
    {
        private:
            // cf. [ABHH1999], p. 8, Table 3
            static const double _v_f0[3];
            static const double _v_c1[3];
            static const double _v_c2[3];

            static const double _a_0_f0[3];
            static const double _a_0_c1[3];
            static const double _a_0_c2[3];

            static const double _a_1_f0[3];
            static const double _a_1_c1[3];
            static const double _a_1_c2[3];

            static const double _a_2_f0[3];
            static const double _a_2_c1[3];
            static const double _a_2_c2[3];

            unsigned _set;

            // cf. [ABHH1999], Eq. (3.7), p. 10 and following text
            double _calc(const double & f0, const double & c1, const double & c2, const double & s_hat)
            {
                return f0 * std::exp(c1 * s_hat + c2 * s_hat * s_hat);
            }

        public:
            ABHH1999FormFactors(unsigned set) :
                _set(set)
            {
            }

            ~ABHH1999FormFactors()
            {
            }

            virtual double v(const double & s_hat)
            {
                return _calc(_v_f0[_set], _v_c1[_set], _v_c2[_set], s_hat);
            }

            virtual double a_0(const double & s_hat)
            {
                return _calc(_a_0_f0[_set], _a_0_c1[_set], _a_0_c2[_set], s_hat);
            }

            virtual double a_1(const double & s_hat)
            {
                return _calc(_a_1_f0[_set], _a_1_c1[_set], _a_1_c2[_set], s_hat);
            }

            virtual double a_2(const double & s_hat)
            {
                return _calc(_a_2_f0[_set], _a_2_c1[_set], _a_2_c2[_set], s_hat);
            }
    };

    const double ABHH1999FormFactors::_v_f0[3] = { 0.399, 0.457, 0.548 };
    const double ABHH1999FormFactors::_v_c1[3] = { 1.537, 1.482, 1.462 };
    const double ABHH1999FormFactors::_v_c2[3] = { 1.123, 1.015, 0.953 };

    const double ABHH1999FormFactors::_a_0_f0[3] = { 0.412, 0.471, 0.698 };
    const double ABHH1999FormFactors::_a_0_c1[3] = { 1.543, 1.505, 1.945 };
    const double ABHH1999FormFactors::_a_0_c2[3] = { 0.954, 0.710, 0.314 };

    const double ABHH1999FormFactors::_a_1_f0[3] = { 0.294, 0.337, 0.385 };
    const double ABHH1999FormFactors::_a_1_c1[3] = { 0.656, 0.602, 0.557 };
    const double ABHH1999FormFactors::_a_1_c2[3] = { 0.456, 0.258, 0.068 };

    const double ABHH1999FormFactors::_a_2_f0[3] = { 0.246, 0.282, 0.320 };
    const double ABHH1999FormFactors::_a_2_c1[3] = { 1.237, 1.172, 1.083 };
    const double ABHH1999FormFactors::_a_2_c2[3] = { 0.822, 0.567, 0.393 };

    /* Form Factors according to [BZ2004] */
    class BZ2004FormFactors :
        public FormFactors<BToKstar>
    {
        private:
            // cf. [BZ2004], Eqs. (59)-(61)
            // Here ratio_{r,fit} is m_B^2 / m_{R,fit}^2. For Eq. (61), we use Eq. (59) with r1 = 0.
            double _calc_eq59(const double & r1, const double & r2, const double & ratio_r, const double & ratio_fit, const double & s_hat)
            {
                double x_r = s_hat * ratio_r;
                double x_fit = s_hat * ratio_fit;

                return r1 / (1.0 - x_r) + r2 / (1.0 - x_fit);
            }

            double _calc_eq60(const double & r1, const double & r2, const double & ratio_fit, const double & s_hat)
            {
                double x_fit = s_hat * ratio_fit;
                double denom = 1.0 - x_fit;

                return r1 / denom + r2 / denom / denom;
            }

        public:
            BZ2004FormFactors()
            {
            }

            ~BZ2004FormFactors()
            {
            }

            // cf. [BZ2004], Table 8, p. 28
            virtual double v(const double & s_hat)
            {
                return _calc_eq59(0.923, -0.511, 0.99248, 0.56434, s_hat);
            }

            virtual double a_0(const double & s_hat)
            {
                return _calc_eq59(1.364, -0.990, 1.0, 0.75798, s_hat);
            }

            virtual double a_1(const double & s_hat)
            {
                return _calc_eq59(0.0, 0.290, 1.0, 0.69040, s_hat);
            }

            virtual double a_2(const double & s_hat)
            {
                return _calc_eq60(-0.084, 0.342, 0.53612, s_hat);
            }
    };

    /* Form Factors according to [IKKR2007] */
    class IKKR2007FormFactors :
        public FormFactors<BToKstar>
    {
        private:
            // cf. [IKKR2007], Eq. (48)
            // Here s_hat = q^2 / m_B^2
            double _calc(const double & f0, const double & a, const double & b, const double & s_hat)
            {
                return f0 / (1.0 - s_hat * a + s_hat * s_hat * b);
            }

            // cf. [IKKR2007], Table III, p. 9
            // For rewriting the form factor to the [BZ2004] Basis, cf. [IKKR2007], Eq. (6), p. 2
            // The ^c basis corresponds to the [BZ2004] Basis.
            virtual double v(const double & s_hat)
            {
                return _calc(0.37, 2.05, 1.13, s_hat);
            }

            virtual double a_0(const double & s_hat)
            {
                static const double alpha = 2.4596;
                static const double beta = 2.5139;

                // A_0 <- alpha (A_0 - A_+) - s_hat beta A_-
                return alpha * (_calc(0.40, 0.98, 0.034, s_hat) - _calc(0.30, 1.92, 0.97, s_hat)) - s_hat * beta * _calc(-0.38, 2.10, 1.19, s_hat);
            }

            virtual double a_1(const double & s_hat)
            {
                return _calc(0.28438, 0.97, 0.034, s_hat);
            }

            virtual double a_2(const double & s_hat)
            {
                return _calc(0.30, 1.92, 0.97, s_hat);
            }
    };

    FormFactors<BToKstar>::~FormFactors()
    {
    }

    std::tr1::shared_ptr<FormFactors<BToKstar>>
    FormFactorFactory<BToKstar>::create(const std::string & label)
    {
        std::tr1::shared_ptr<FormFactors<BToKstar>> result;
        std::string name, set;

        std::string::size_type sep(label.find("::"));
        if (std::string::npos == sep)
        {
            name = label;
        }
        else
        {
            name = label.substr(0, sep);
            set = label.substr(sep + 2);
        }

        if ("ABHH1999" == name)
        {
            if (set.empty())
                set = "1";

            unsigned index(set[0] - '1');

            if (index < 3)
                result = std::tr1::shared_ptr<FormFactors<BToKstar>>(new ABHH1999FormFactors(index));
        }
        else if ("BZ2004" == name)
        {
            result = std::tr1::shared_ptr<FormFactors<BToKstar>>(new BZ2004FormFactors);
        }
        else if ("IKKR2007" == name)
        {
            result = std::tr1::shared_ptr<FormFactors<BToKstar>>(new IKKR2007FormFactors);
        }

        return result;
    }
}
