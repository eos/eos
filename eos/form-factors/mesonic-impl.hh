/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013, 2014 Danny van Dyk
 * Copyright (c) 2010, 2011 Christian Wacker
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_FORM_FACTORS_IMPL_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_FORM_FACTORS_IMPL_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/analytic-b-to-pi.hh>
#include <eos/utils/derivative.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/power_of.hh>

namespace eos
{
    /* Form Factors according to [BZ2004] */
    template <typename Process_, typename Transition_> class BZ2004FormFactors;

    /* Form Factors according to [KMPW2010] */
    template <typename Transition_> class KMPW2010FormFactors;

    /* Form Factors according to [BFW2010] */
    template <typename Transition_> class BFW2010FormFactors;

    /* Form Factors according to [BCL2008] */
    template <typename Process_> class BCL2008FormFactors;

    /* P -> V Processes */

    struct BToKstar {
        static constexpr double mB = 5.279;
        static constexpr double mV = 0.896;
    };

    struct BsToPhi {
        static constexpr double mB = 5.366;
        static constexpr double mV = 1.020;
    };

    template <typename Process_> class BZ2004FormFactors<Process_, PToV> :
        public FormFactors<PToV>
    {
        private:
            UsedParameter _v_factor, _a0_factor, _a1_factor, _a2_factor;

            // fit parametrisation for P -> V according to [BZ2004]
            static const double _v_r1, _v_r2, _v_m2r, _v_m2fit;
            static const double _a0_r1, _a0_r2, _a0_m2r, _a0_m2fit;
            static const double _a1_r2, _a1_m2fit;
            static const double _a2_r1, _a2_r2, _a2_m2fit;

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
            BZ2004FormFactors(const Parameters & p, const Options &) :
                _v_factor(p["B->K^*::v_uncertainty@BZ2004"], *this),
                _a0_factor(p["B->K^*::a0_uncertainty@BZ2004"], *this),
                _a1_factor(p["B->K^*::a1_uncertainty@BZ2004"], *this),
                _a2_factor(p["B->K^*::a2_uncertainty@BZ2004"], *this)
            {
            }

            static FormFactors<PToV> * make(const Parameters & parameters, unsigned)
            {
                return new BZ2004FormFactors(parameters, Options());
            }

            virtual double v(const double & s) const
            {
                return _v_factor * _calc_eq59(s, _v_r1, _v_r2, _v_m2r, _v_m2fit);
            }

            virtual double a_0(const double & s) const
            {
                return _a0_factor * _calc_eq59(s, _a0_r1, _a0_r2, _a0_m2r, _a0_m2fit);
            }

            virtual double a_1(const double & s) const
            {
                return _a1_factor * _calc_eq61(s, _a1_r2, _a1_m2fit);
            }

            virtual double a_2(const double & s) const
            {
                return _a2_factor * _calc_eq60(s, _a2_r1, _a2_r2, _a2_m2fit);
            }

            virtual double a_12(const double & s) const
            {
                const double mB = Process_::mB, mB2 = mB * mB;
                const double mV = Process_::mV, mV2 = mV * mV;
                const double lambda = eos::lambda(mB2, mV2, s);

                return ((mB + mV) * (mB + mV) * (mB2 - mV2 - s) * this->a_1(s)
                    - lambda * this->a_2(s)) / (16.0 * mB * mV2 * (mB + mV));
            }
    };

    template <> class KMPW2010FormFactors<PToV> :
        public FormFactors<PToV>
    {
        private:
            // fit parametrisation for P -> V according to [KMPW2010]
            UsedParameter _f0_V, _b1_V;
            UsedParameter _f0_A0, _b1_A0;
            UsedParameter _f0_A1, _b1_A1;
            UsedParameter _f0_A2, _b1_A2;
            static const double _tau_p, _tau_m, _tau_0;
            static const double _m_B, _m_Kstar, _m_Bs2_0m, _m_Bs2_1m, _m_Bs2_1p;

            static double _calc_z(const double & s)
            {
                return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
            }

        public:
            KMPW2010FormFactors(const Parameters & p, const Options &) :
                _f0_V(p["B->K^*::F^V(0)@KMPW2010"],   *this),
                _b1_V(p["B->K^*::b^V_1@KMPW2010"],    *this),
                _f0_A0(p["B->K^*::F^A0(0)@KMPW2010"], *this),
                _b1_A0(p["B->K^*::b^A0_1@KMPW2010"],  *this),
                _f0_A1(p["B->K^*::F^A1(0)@KMPW2010"], *this),
                _b1_A1(p["B->K^*::b^A1_1@KMPW2010"],  *this),
                _f0_A2(p["B->K^*::F^A2(0)@KMPW2010"], *this),
                _b1_A2(p["B->K^*::b^A2_1@KMPW2010"],  *this)
            {
            }

            static FormFactors<PToV> * make(const Parameters & parameters, unsigned)
            {
                return new KMPW2010FormFactors(parameters, Options());
            }

            virtual double v(const double & s) const
            {
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                // cf. [KMPW2010], Eq. (8.8), p. 30
                return _f0_V() / (1.0 - s / _m_Bs2_1m) * (1.0 + _b1_V() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

            virtual double a_0(const double & s) const
            {
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                // cf. [KMPW2010], Eq. (8.8), p. 30
                return _f0_A0() / (1.0 - s / _m_Bs2_0m) * (1.0 + _b1_A0() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

            virtual double a_1(const double & s) const
            {
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                // cf. [KMPW2010], Eq. (8.8), p. 30
                return _f0_A1() / (1.0 - s / _m_Bs2_1p) * (1.0 + _b1_A1() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

            virtual double a_2(const double & s) const
            {
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                // cf. [KMPW2010], Eq. (8.8), p. 30
                return _f0_A2() / (1.0 - s / _m_Bs2_1p) * (1.0 + _b1_A2() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }
            virtual double a_12(const double & s) const
            {
                const double mB = BToKstar::mB, mB2 = mB * mB;
                const double mV = BToKstar::mV, mV2 = mV * mV;
                const double lambda = eos::lambda(mB2, mV2, s);

                return ((mB + mV) * (mB + mV) * (mB2 - mV2 - s) * this->a_1(s)
                    - lambda * this->a_2(s)) / (16.0 * mB * mV2 * (mB + mV));
            }
    };

    /* P -> P Processes */

    double FormFactors<PToP>::f_p_d1(const double & s) const
    {
        using namespace std::placeholders;

        auto f = std::function<double (const double &)>(std::bind(&FormFactors<PToP>::f_p, this, _1));

        return derivative<1u, deriv::TwoSided>(f, s);
    }

    double FormFactors<PToP>::f_p_d2(const double & s) const
    {
        using namespace std::placeholders;

        auto f = std::function<double (const double &)>(std::bind(&FormFactors<PToP>::f_p, this, _1));

        return derivative<2u, deriv::TwoSided>(f, s);
    }

    struct BToK {
        typedef PToP Transition;
        static constexpr const char * label = "B->K";
        static constexpr const double m_B = 5.279;
        static constexpr const double m_P = 0.492;
        static constexpr const double m2_Br0m = 5.366 * 5.366; // B_s
        static constexpr const double m2_Br1m = 5.415 * 5.415; // B_s^*
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
    };

    struct BToPi {
        typedef PToP Transition;
        static constexpr const char * label = "B->pi";
        static constexpr const double m_B = 5.279;
        static constexpr const double m_P = 0.135;
        static constexpr const double m2_Br0m = 5.279 * 5.279; // B_{u,d}
        static constexpr const double m2_Br1m = 5.325 * 5.325; // B_{u,d}^*
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
    };


    /*
     * P -> P Form Factors in the simplified series expansion according to
     * [BCL2008].
     */
    template <typename Process_> class BCL2008FormFactors :
        public FormFactors<typename Process_::Transition>
    {
        private:
            /*
             * Fit parametrisation for P -> P according to [BCL2008], eq. (11)
             * with K = 3. Note that we factor out the form factor at q^2 = 0
             * by setting t_0 = 0.0, thus b_k -> b_k / b_0. Note that the last
             * coefficient b_K is fixed by eq. (14).
             */
            UsedParameter _f_plus_0, _b_plus_1, _b_plus_2;
            UsedParameter            _b_zero_1, _b_zero_2;
            UsedParameter _f_t_0,    _b_t_1,    _b_t_2;

            static double _z(const double & s)
            {
                static const double m_B = Process_::m_B;
                static const double m_P = Process_::m_P;
                static const double tau_p = Process_::tau_p;
                static const double tau_0 = (m_B + m_P) * (std::sqrt(m_B) - std::sqrt(m_P)) * (std::sqrt(m_B) - std::sqrt(m_P));

                return (std::sqrt(tau_p - s) - std::sqrt(tau_p - tau_0))
                    / (std::sqrt(tau_p - s) + std::sqrt(tau_p - tau_0));
            }

        public:
            BCL2008FormFactors(const Parameters & p, const Options &) :
                _f_plus_0(p[std::string(Process_::label) + "::f_+(0)@BCL2008"], *this),
                _b_plus_1(p[std::string(Process_::label) + "::b_+^1@BCL2008"],  *this),
                _b_plus_2(p[std::string(Process_::label) + "::b_+^2@BCL2008"],  *this),
                _b_zero_1(p[std::string(Process_::label) + "::b_0^1@BCL2008"],  *this),
                _b_zero_2(p[std::string(Process_::label) + "::b_0^2@BCL2008"],  *this),
                _f_t_0(p[std::string(Process_::label)    + "::f_T(0)@BCL2008"], *this),
                _b_t_1(p[std::string(Process_::label)    + "::b_T^1@BCL2008"],  *this),
                _b_t_2(p[std::string(Process_::label)    + "::b_T^2@BCL2008"],  *this)
            {
            }

            static FormFactors<PToP> * make(const Parameters & parameters, unsigned)
            {
                return new BCL2008FormFactors(parameters, Options());
            }

            virtual double f_p(const double & s) const
            {
                const double z = _z(s), z2 = z * z, z3 = z * z2;
                const double z0 = _z(0), z02 = z0 * z0, z03 = z0 * z02;
                const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03;

                return _f_plus_0 / (1.0 - s / Process_::m2_Br1m) * (1.0 + _b_plus_1 * (zbar - z3bar / 3.0) + _b_plus_2 * (z2bar + 2.0 * z3bar / 3.0));
            }

            virtual double f_0(const double & s) const
            {
                const double z = _z(s), z2 = z * z, z3 = z * z2;
                const double z0 = _z(0), z02 = z0 * z0, z03 = z0 * z02;
                const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03;

                // note that f_0(0) = f_+(0)!
                return _f_plus_0 / (1.0 - s / Process_::m2_Br1m) * (1.0 + _b_zero_1 * (zbar - z3bar / 3.0) + _b_zero_2 * (z2bar + 2.0 * z3bar / 3.0));
            }

            virtual double f_t(const double & s) const
            {
                const double z = _z(s), z2 = z * z, z3 = z * z2;
                const double z0 = _z(0), z02 = z0 * z0, z03 = z0 * z02;
                const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03;

                return _f_t_0 / (1.0 - s / Process_::m2_Br1m) * (1.0 + _b_t_1 * (zbar - z3bar / 3.0) + _b_t_2 * (z2bar + 2.0 * z3bar / 3.0));
            }
    };

    /* Form Factors according to [BZ2004v2] */
    template <typename Process_> class BZ2004FormFactors<Process_, PToP> :
        public FormFactors<PToP>
    {
        private:
            UsedParameter _f_p_factor, _f_0_factor, _f_t_factor;

            // fit parametrisation for P -> P according to [BZ2004v2]
            static const double _r1_p, _r2_p, _r1_t, _r2_t, _r2_0;
            static const double _mfit2, _m12;

        public:
            BZ2004FormFactors(const Parameters & p, const Options &) :
                _f_p_factor(p["B->K::fp_uncertainty@BZ2004v2"], *this),
                _f_0_factor(p["B->K::f0_uncertainty@BZ2004v2"], *this),
                _f_t_factor(p["B->K::ft_uncertainty@BZ2004v2"], *this)
            {
            }

            static FormFactors<PToP> * make(const Parameters & parameters, unsigned)
            {
                return new BZ2004FormFactors(parameters, Options());
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


    // Form Factors according to [BZ2004v3]
    template <typename Process_> class BZ2004FormFactorsSplit :
        public FormFactors<PToP>
    {
        private:
            UsedParameter _f_p_factor, _f_0_factor, _f_t_factor;

            static const double _r1_p_asymptotic, _r2_p_asymptotic;
            static const double _r2_0_asymptotic;
            static const double _r1_t_asymptotic, _r2_t_asymptotic;
            static const double _mfit2_0_asymptotic;

            static const double _m12_asymptotic;

            static const double _f_p_a_1, _f_p_b_1, _f_p_c_1, _f_p_d_1;
            static const double _f_p_a_2, _f_p_b_2, _f_p_c_2, _f_p_d_2;
            static const double _f_p_a_4, _f_p_b_4, _f_p_c_4, _f_p_d_4;

            static const double _f_0_a_1, _f_0_b_1, _f_0_c_1, _f_0_d_1;
            static const double _f_0_a_2, _f_0_b_2, _f_0_c_2, _f_0_d_2;
            static const double _f_0_a_4, _f_0_b_4, _f_0_c_4, _f_0_d_4;

            static const double _f_t_a_1, _f_t_b_1, _f_t_c_1, _f_t_d_1;
            static const double _f_t_a_2, _f_t_b_2, _f_t_c_2, _f_t_d_2;
            static const double _f_t_a_4, _f_t_b_4, _f_t_c_4, _f_t_d_4;

            // Gegenbauer moments
            UsedParameter _a_1, _a_2, _a_4;

            // Polynomial of degree 3, cf. [BZ2004v3], eq. (A.6), p. 28
            double poly3(const double & s, const double & a, const double & b, const double & c, const double & d) const
            {
                return a + b * s + c * s * s + d * s * s * s;
            }

            double f_p_asymptotic(const double & s) const
            {
                // cf. [BZ2004v3], eq. (A.2), p. 26
                return _r1_p_asymptotic / (1.0 - s / _m12_asymptotic) + _r2_p_asymptotic / power_of<2>(1.0 - s / _m12_asymptotic);
            }

            double f_0_asymptotic(const double & s) const
            {
                // cf. [BZ2004v3], eq. (A.3), p. 26
                return _r2_0_asymptotic / (1.0 - s / _mfit2_0_asymptotic);
            }

            double f_t_asymptotic(const double & s) const
            {
                // cf. [BZ2004v3], eq. (A.2), p. 26
                return _r1_t_asymptotic / (1.0 - s / _m12_asymptotic) + _r2_t_asymptotic / power_of<2>(1.0 - s / _m12_asymptotic);
            }

        public:
            BZ2004FormFactorsSplit(const Parameters & p, const Options &) :
                _f_p_factor(p["B->K::fp_uncertainty@BZ2004v2"], *this),
                _f_0_factor(p["B->K::f0_uncertainty@BZ2004v2"], *this),
                _f_t_factor(p["B->K::ft_uncertainty@BZ2004v2"], *this),
                _a_1(p["B->K::a_1@2.2GeV"], *this),
                _a_2(p["B->K::a_2@2.2GeV"], *this),
                _a_4(p["B->K::a_4@2.2GeV"], *this)
            {
            }

            static FormFactors<PToP> * make(const Parameters & parameters, unsigned)
            {
                return new BZ2004FormFactorsSplit(parameters, Options());
            }

            virtual double f_p(const double & s) const
            {
                // cf. [BZ2004v3], eq. (A.5), p. 27
                return _f_p_factor * (f_p_asymptotic(s) +
                                      _a_1 * poly3(s, _f_p_a_1, _f_p_b_1, _f_p_c_1, _f_p_d_1) +
                                      _a_2 * poly3(s, _f_p_a_2, _f_p_b_2, _f_p_c_2, _f_p_d_2) +
                                      _a_4 * poly3(s, _f_p_a_4, _f_p_b_4, _f_p_c_4, _f_p_d_4));
            }

            virtual double f_0(const double & s) const
            {
                // cf. [BZ2004v3], eq. (A.5), p. 27
                return _f_0_factor * (f_0_asymptotic(s) +
                                      _a_1 * poly3(s, _f_0_a_1, _f_0_b_1, _f_0_c_1, _f_0_d_1) +
                                      _a_2 * poly3(s, _f_0_a_2, _f_0_b_2, _f_0_c_2, _f_0_d_2) +
                                      _a_4 * poly3(s, _f_0_a_4, _f_0_b_4, _f_0_c_4, _f_0_d_4));
            }

            virtual double f_t(const double & s) const
            {
                // cf. [BZ2004v3], eq. (A.5), p. 27
                return _f_t_factor * (f_t_asymptotic(s) +
                                      _a_1 * poly3(s, _f_t_a_1, _f_t_b_1, _f_t_c_1, _f_t_d_1) +
                                      _a_2 * poly3(s, _f_t_a_2, _f_t_b_2, _f_t_c_2, _f_t_d_2) +
                                      _a_4 * poly3(s, _f_t_a_4, _f_t_b_4, _f_t_c_4, _f_t_d_4));
            }
    };

    template <> class KMPW2010FormFactors<PToP> :
        public FormFactors<PToP>
    {
        private:
            // fit parametrisation for P -> P according to [KMPW2010]
            UsedParameter _b1_p, _b1_0, _b1_t;
            UsedParameter _f0_p, _f0_t;
            static const double _tau_p, _tau_m, _tau_0;
            static const double _m_B, _m_K, _m_Bs2;

            static double _calc_z(const double & s)
            {
                return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
            }

        public:
            KMPW2010FormFactors(const Parameters & p, const Options &) :
                _b1_p(p["B->K::b^p_1@KMPW2010"], *this),
                _b1_0(p["B->K::b^0_1@KMPW2010"], *this),
                _b1_t(p["B->K::b^t_1@KMPW2010"], *this),
                _f0_p(p["B->K::F^p(0)@KMPW2010"], *this),
                _f0_t(p["B->K::F^t(0)@KMPW2010"], *this)
            {
            }

            static FormFactors<PToP> * make(const Parameters & parameters, unsigned)
            {
                return new KMPW2010FormFactors(parameters, Options());
            }

            virtual double f_p(const double & s) const
            {
                // cf. [KMPW2010], Eq. (8.8), p. 30
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                return _f0_p() / (1 - s / _m_Bs2) * (1 + _b1_p() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

            virtual double f_0(const double & s) const
            {
                // cf. [KMPW2010], Eq. (8.8), p. 30
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                return _f0_p() * (1 + _b1_0() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

            virtual double f_t(const double & s) const
            {
                // cf. [KMPW2010], Eq. (8.8), p. 30
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                return _f0_t() / (1 - s / _m_Bs2) * (1 + _b1_t() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }
    };

    template <> class BFW2010FormFactors<PToP> :
          public FormFactors<PToP>
    {
          private:
              // fit parametrisation for P -> P according to [BFW2011]. We use the simple series expansion and
              // the results form LCSR only.
              UsedParameter _alpha_V0_0, _alpha_V0_1;
              UsedParameter _alpha_Vt_0np,_alpha_Vt_1np;
              UsedParameter _alpha_T0_0, _alpha_T0_1;
              static const double _tau_p, _tau_m, _tau_0;
              static const double _m_B, _m_K, _m_Bs2;

              static double _calc_z(const double & s)
              {
                  return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
              }

          public:
              BFW2010FormFactors(const Parameters & p, const Options &) :
                  _alpha_V0_0(p["B->K::alpha^V0_0@BFW2010"], *this),
                  _alpha_V0_1(p["B->K::alpha^V0_1@BFW2010"], *this),
                  _alpha_Vt_0np(p["B->K::alpha^Vt_0np@BFW2010"], *this),
                  _alpha_Vt_1np(p["B->K::alpha^Vt_1np@BFW2010"], *this),
                  _alpha_T0_0(p["B->K::alpha^T0_0@BFW2010"], *this),
                  _alpha_T0_1(p["B->K::alpha^T0_1@BFW2010"], *this)
              {
              }

              static FormFactors<PToP> * make(const Parameters & parameters, unsigned)
              {
                  return new BFW2010FormFactors(parameters, Options());
              }

              virtual double f_p(const double & s) const
              {
                  // cf. [BFW2010], Eq. (43), p. 13, replacements Eq. (45), p. 16 and Eq. (7), p. 4
                  return 1.0 / (1.0 - s / _m_Bs2) * (_alpha_V0_0 + _alpha_V0_1 * _calc_z(s));
              }

              virtual double f_0(const double & s) const
              {
                  // cf. [BFW2010], Eq. (43), p. 13, replacements Eq. (45), p. 16 and Eq. (7), p. 4
                  static const double prefactor = _m_B * _m_B / (_m_B * _m_B -_m_K * _m_K);
                  return prefactor * (_alpha_Vt_0np + _alpha_Vt_1np * _calc_z(s));
              }

              virtual double f_t(const double & s) const
              {
                  // cf. [BFW2010], Eq. (43), p. 13, replacements Eq. (45), p. 16 and Eq. (9), p. 4
                  static const double prefactor = (_m_B + _m_K) / _m_B;
                  return prefactor / (1.0 - s / _m_Bs2) * (_alpha_T0_0 + _alpha_T0_1 * _calc_z(s));
              }
    };
}

#endif
