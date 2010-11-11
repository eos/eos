/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/form-factors.hh>
#include <src/utils/destringify.hh>
#include <src/utils/power_of.hh>

#include <map>
#include <cmath>

namespace eos
{
    /* P -> V Processes */

    struct BToKstar { };
    struct BsToPhi { };

    template <typename Process_, typename Transition_> class BZ2004FormFactors;

    /* Form Factors according to [BZ2004] */
    template <typename Process_> class BZ2004FormFactors<Process_, PToV> :
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
    template class BZ2004FormFactors<BToKstar, PToV>;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_v_r1 = +0.923;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_v_r2 = -0.511;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_v_m2r = 5.32 * 5.32;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_v_m2fit = 49.40;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a0_r1 = +1.364;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a0_r2 = -0.990;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a0_m2r = 5.28 * 5.28;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a0_m2fit = 36.78;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a1_r2 = +0.290;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a1_m2fit = 40.38;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a2_r1 = -0.084;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a2_r2 = +0.342;
    template <> const double BZ2004FormFactors<BToKstar, PToV>::_a2_m2fit = 52.00;

    /* B_s -> phi */
    template class BZ2004FormFactors<BsToPhi, PToV>;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_v_r1 = +1.484;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_v_r2 = -1.049;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_v_m2r = 5.42 * 5.42;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_v_m2fit = 39.52;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a0_r1 = +3.310;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a0_r2 = -2.835;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a0_m2r = 5.37 * 5.37;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a0_m2fit = 31.57;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a1_r2 = +0.308;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a1_m2fit = 36.54;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a2_r1 = -0.054;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a2_r2 = +0.288;
    template <> const double BZ2004FormFactors<BsToPhi, PToV>::_a2_m2fit = 48.94;

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
            { KeyType("B->K^*",     "BZ2004"), &BZ2004FormFactors<BToKstar, PToV>::make   },
            { KeyType("Bs->phi",    "BZ2004"), &BZ2004FormFactors<BsToPhi, PToV>::make    },
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

    /* P -> P Processes */

    struct BToK { };

    /* Form Factors according to [BZ2004v2] */
    template <typename Process_> class BZ2004FormFactors<Process_, PToP> :
        public FormFactors<PToP>
    {
        private:
            Parameter _f_p_factor, _f_0_factor, _f_t_factor;

            // fit parametrisation for P -> P according to [BZ2004v2]
            static const double _r1_p, _r2_p, _r1_t, _r2_t, _r2_0;
            static const double _mfit2, _m12;

            BZ2004FormFactors(const Parameters & p) :
                _f_p_factor(p["formfactors::fp_uncertainty"]),
                _f_0_factor(p["formfactors::f0_uncertainty"]),
                _f_t_factor(p["formfactors::ft_uncertainty"])
            {
            }


        public:
            static FormFactors<PToP> * make(const Parameters & parameters, unsigned)
            {
                return new BZ2004FormFactors(parameters);
            }

            virtual double f_p(const double & s) const
            {
                // [BZ2004v2] eq. (11)
                return _f_p_factor * (_r1_p / (1 - s / _m12) + _r2_p / power_of<2>(1 - s / _m12));
            }

            virtual double f_0(const double & s) const
            {
                // [BZ2004v2] eq. (12)
                return _f_0_factor * (_r2_0 / (1 - s / _mfit2));
            }

            virtual double f_t(const double & s) const
            {
                // [BZ2004v2] eq. (11)
                return _f_t_factor * (_r1_t / (1 - s / _m12) + _r2_t / power_of<2>(1 - s / _m12));
            }
    };


    /* B_{u,d} -> K */
    template class BZ2004FormFactors<BToK, PToP>;

    /* For the values below, cf. [BZ2004v2], Table 1, p. 8 */
    template <> const double BZ2004FormFactors<BToK, PToP>::_r1_p     = 0.162;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r2_p     = 0.173;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r1_t     = 0.161;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r2_t     = 0.198;
    template <> const double BZ2004FormFactors<BToK, PToP>::_r2_0     = 0.330;
    template <> const double BZ2004FormFactors<BToK, PToP>::_mfit2    = 37.46;
    template <> const double BZ2004FormFactors<BToK, PToP>::_m12      = 5.41 * 5.41;


    /* Form Factors according to [KMPW2010] */
    template <typename Process_> class KMPW2010FormFactors :
        public FormFactors<PToP>
    {
        private:
            Parameter _f_p_factor, _f_0_factor, _f_t_factor;

            // fit parametrisation for P -> P according to [KMPW2010]
            static const double _b1_p, _b1_0, _b1_t;
            static const double _f0_p, _f0_0, _f0_t;
            static const double _tau_p, _tau_m, _tau_0;
            static const double _m_B, _m_K, _m_Bs2;

            KMPW2010FormFactors(const Parameters & p) :
                _f_p_factor(p["formfactors::fp_uncertainty"]),
                _f_0_factor(p["formfactors::f0_uncertainty"]),
                _f_t_factor(p["formfactors::ft_uncertainty"])
            {
            }

            static double _calc_z(const double & s)
            {
                return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
            }


        public:
            static FormFactors<PToP> * make(const Parameters & parameters, unsigned)
            {
                return new KMPW2010FormFactors(parameters);
            }

            virtual double f_p(const double & s) const
            {
                // [KMPW2010] eq. (8.8)
                return _f_p_factor * _f0_p / (1 - s / _m_Bs2) * (1 + _b1_p * (_calc_z(s) - _calc_z(0) + 0.5 * (power_of<2>(_calc_z(s)) - power_of<2>(_calc_z(0)))));
            }

            virtual double f_0(const double & s) const
            {
                // [KMPW2010] eq. (8.8)
                return _f_0_factor * _f0_0 * (1 + _b1_0 * (_calc_z(s) - _calc_z(0) + 0.5 * (power_of<2>(_calc_z(s)) - power_of<2>(_calc_z(0)))));
            }

            virtual double f_t(const double & s) const
            {
                // [KMPW2010] eq. (8.8)
                return _f_t_factor * _f0_t / (1 - s / _m_Bs2) * (1 + _b1_t * (_calc_z(s) - _calc_z(0) + 0.5 * (power_of<2>(_calc_z(s)) - power_of<2>(_calc_z(0)))));
            }
    };

    /* B_{u,d} -> K */
    template class KMPW2010FormFactors<BToK>;

    /* For the values below, cf. [KMPW2010], Table 4, p. 31 */
    template <> const double KMPW2010FormFactors<BToK>::_f0_p     =  0.34;
    template <> const double KMPW2010FormFactors<BToK>::_f0_0     =  0.34;
    template <> const double KMPW2010FormFactors<BToK>::_f0_t     =  0.39;
    template <> const double KMPW2010FormFactors<BToK>::_b1_p     = -2.1;
    template <> const double KMPW2010FormFactors<BToK>::_b1_0     = -4.3;
    template <> const double KMPW2010FormFactors<BToK>::_b1_t     = -2.2;

    template <> const double KMPW2010FormFactors<BToK>::_tau_p    = (_m_B + _m_K) * (_m_B + _m_K);
    template <> const double KMPW2010FormFactors<BToK>::_tau_m    = (_m_B - _m_K) * (_m_B - _m_K);
    template <> const double KMPW2010FormFactors<BToK>::_tau_0    = _tau_p - std::sqrt(_tau_p * _tau_p - _tau_m * _tau_p);

    // Masses, cf. PDG 2008
    template <> const double KMPW2010FormFactors<BToK>::_m_B      = 5.280;
    template <> const double KMPW2010FormFactors<BToK>::_m_K      = 0.498;
    template <> const double KMPW2010FormFactors<BToK>::_m_Bs2    = 5.325 * 5.325;

    FormFactors<PToP>::~FormFactors()
    {
    }

    std::shared_ptr<FormFactors<PToP>>
    FormFactorFactory<PToP>::create(const std::string & label, const Parameters & parameters)
    {
        std::shared_ptr<FormFactors<PToP>> result;

        typedef std::tuple<std::string, std::string> KeyType;
        typedef std::function<FormFactors<PToP> * (const Parameters &, unsigned)> ValueType;
        static const std::map<KeyType, ValueType> form_factors
        {
            { KeyType("B->K",     "BZ2004v2"), &BZ2004FormFactors<BToK, PToP>::make   },
            { KeyType("B->K",     "KMPW2010"), &KMPW2010FormFactors<BToK>::make   }
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

        result = std::shared_ptr<FormFactors<PToP>>(i->second(parameters, set));

        return result;
    }
}
