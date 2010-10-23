/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/form-factors.hh>
#include <src/utils/destringify.hh>

#include <map>

namespace eos
{
    /* P -> V Processes */

    struct BToKstar { };
    struct BsToPhi { };

    /* Form Factors according to [BZ2004] */
    template <typename Process_> class BZ2004FormFactors :
        public FormFactors<PToV>
    {
        private:
            Parameter _v_factor, _a0_factor, _a1_factor, _a2_factor;

            // fit parametrisation for P -> V according to [BZ2004]
            static const double _v_r1, _v_r2, _v_m2r, _v_m2fit;
            static const double _a0_r1, _a0_r2, _a0_m2r, _a0_m2fit;
            static const double _a1_r2, _a1_m2fit;
            static const double _a2_r1, _a2_r2, _a2_m2fit;

            BZ2004FormFactors(const Parameters & p) :
                _v_factor(p["formfactors::v_uncertainty"]),
                _a0_factor(p["formfactors::a0_uncertainty"]),
                _a1_factor(p["formfactors::a1_uncertainty"]),
                _a2_factor(p["formfactors::a2_uncertainty"])
            {
            }

            // cf. [BZ2004], Eq. 59, p. 27
            static inline double _calc_eq59(const double & s, const double & r_1, const double & r_2, const double & m2r, const double & m2fit)
            {
                return r_1 / (1.0 - s / m2r) + r_2 / (1.0 - s / m2fit);
            }

            // cf. [BZ2004], Eq. 60, p. 29
            static inline double _calc_eq60(const double & s, const double & r_1, const double & r_2, const double & m2fit)
            {
                double denom = 1.0 - s / m2fit;

                return r_1 / denom + r_2 / denom / denom;
            }

            // cf. [BZ2004], Eq. 61, p. 29
            static inline double _calc_eq61(const double & s, const double & r_2, const double & m2fit)
            {
                return r_2 / (1.0 - s / m2fit);
            }

        public:
            static FormFactors<PToV> * make(const Parameters & parameters, unsigned)
            {
                return new BZ2004FormFactors(parameters);
            }

            virtual double v(const double & s) const
            {
                return _calc_eq59(s, _v_r1, _v_r2, _v_m2r, _v_m2fit);
            }

            virtual double a_0(const double & s) const
            {
                return _calc_eq59(s, _a0_r1, _a0_r2, _a0_m2r, _a0_m2fit);
            }

            virtual double a_1(const double & s) const
            {
                return _calc_eq61(s, _a1_r2, _a1_m2fit);
            }

            virtual double a_2(const double & s) const
            {
                return _calc_eq60(s, _a2_r1, _a2_r2, _a2_m2fit);
            }
    };

    /* For the values below, cf. [BZ2004], Table 8, p. 28 */

    /* B_{u,d} -> K^* */
    template class BZ2004FormFactors<BToKstar>;
    template <> const double BZ2004FormFactors<BToKstar>::_v_r1 = +0.923;
    template <> const double BZ2004FormFactors<BToKstar>::_v_r2 = -0.511;
    template <> const double BZ2004FormFactors<BToKstar>::_v_m2r = 5.32 * 5.32;
    template <> const double BZ2004FormFactors<BToKstar>::_v_m2fit = 49.40;
    template <> const double BZ2004FormFactors<BToKstar>::_a0_r1 = +1.364;
    template <> const double BZ2004FormFactors<BToKstar>::_a0_r2 = -0.990;
    template <> const double BZ2004FormFactors<BToKstar>::_a0_m2r = 5.28 * 5.28;
    template <> const double BZ2004FormFactors<BToKstar>::_a0_m2fit = 36.78;
    template <> const double BZ2004FormFactors<BToKstar>::_a1_r2 = +0.290;
    template <> const double BZ2004FormFactors<BToKstar>::_a1_m2fit = 40.38;
    template <> const double BZ2004FormFactors<BToKstar>::_a2_r1 = -0.084;
    template <> const double BZ2004FormFactors<BToKstar>::_a2_r2 = +0.342;
    template <> const double BZ2004FormFactors<BToKstar>::_a2_m2fit = 52.00;

    /* B_s -> phi */
    template class BZ2004FormFactors<BsToPhi>;
    template <> const double BZ2004FormFactors<BsToPhi>::_v_r1 = +1.484;
    template <> const double BZ2004FormFactors<BsToPhi>::_v_r2 = -1.049;
    template <> const double BZ2004FormFactors<BsToPhi>::_v_m2r = 5.42 * 5.42;
    template <> const double BZ2004FormFactors<BsToPhi>::_v_m2fit = 39.52;
    template <> const double BZ2004FormFactors<BsToPhi>::_a0_r1 = +3.310;
    template <> const double BZ2004FormFactors<BsToPhi>::_a0_r2 = -2.835;
    template <> const double BZ2004FormFactors<BsToPhi>::_a0_m2r = 5.37 * 5.37;
    template <> const double BZ2004FormFactors<BsToPhi>::_a0_m2fit = 31.57;
    template <> const double BZ2004FormFactors<BsToPhi>::_a1_r2 = +0.308;
    template <> const double BZ2004FormFactors<BsToPhi>::_a1_m2fit = 36.54;
    template <> const double BZ2004FormFactors<BsToPhi>::_a2_r1 = -0.054;
    template <> const double BZ2004FormFactors<BsToPhi>::_a2_r2 = +0.288;
    template <> const double BZ2004FormFactors<BsToPhi>::_a2_m2fit = 48.94;

    FormFactors<PToV>::~FormFactors()
    {
    }

    std::shared_ptr<FormFactors<PToV>>
    FormFactorFactory<PToV>::create(const std::string & label, const Parameters & parameters)
    {
        std::shared_ptr<FormFactors<PToV>> result;

        typedef std::tuple<std::string, std::string> KeyType;
        typedef std::function<FormFactors<PToV> * (const Parameters &, unsigned)> ValueType;
        static const std::map<KeyType, ValueType> form_factors
        {
            { KeyType("B->K^*",     "BZ2004"), &BZ2004FormFactors<BToKstar>::make   },
            { KeyType("Bs->phi",    "BZ2004"), &BZ2004FormFactors<BsToPhi>::make    },
        };

        /*
         * Labels have the form
         *
         *   PROCESS@NAME[:SET]
         *
         * The brackets indicate the latter part to be optional.
         */

        std::string process, name, input(label);
        unsigned set(0);

        std::string::size_type sep_at(input.find('@')), sep_colon(input.find(':'));
        if (std::string::npos == sep_at)
            return result;

        if (std::string::npos != sep_colon)
        {
            set = destringify<unsigned>(input.substr(sep_colon + 1));
            input.erase(sep_colon + 1);
        }

        name = input.substr(sep_at + 1);
        process = input.substr(0, sep_at);

        auto i = form_factors.find(KeyType(process, name));
        if (form_factors.cend() == i)
            return result;

        result = std::shared_ptr<FormFactors<PToV>>(i->second(parameters, set));

        return result;
    }
}
