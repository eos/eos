/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013-2016, 2018 Danny van Dyk
 * Copyright (c) 2010, 2011 Christian Wacker
 * Copyright (c) 2015 Frederik Beaujean
 * Copyright (c) 2015 Christoph Bobeth
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2018 Keri Vos
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
#include <array>

namespace eos
{
    // Process    = B -> K^*, B -> D^*, B -> rho, B_s -> phi etc.
    // Transition = B -> V or B -> P

    /* Form Factors according to [BZ2004] */
    template <typename Process_, typename Transition_> class BZ2004FormFactors;

    /* Form Factors according to [KMPW2010] */
    template <typename Transition_> class KMPW2010FormFactors;

    /* Form Factors according to [BFW2010] */
    template <typename Process_, typename Transition_> class BFW2010FormFactors;

    /* Form Factors according to [BCL2008] */
    template <typename Process_> class BCL2008FormFactors;

    /* Form Factors according to [BSZ2015] */
    template <typename Process_, typename Transition_> class BSZ2015FormFactors;

    /* P -> V Processes */

    struct BToDstar {
        static constexpr const char * label = "B->D^*";
        static constexpr double mB = 5.279;
        static constexpr double mV = 2.0103;
        static constexpr double mBc = 6.2751;
        static constexpr double mR2_0m = (mBc + 0.000) * (mBc + 0.000);
        static constexpr double mR2_1m = (mBc + 0.056) * (mBc + 0.056);
        static constexpr double mR2_1p = (mBc + 0.492) * (mBc + 0.492);
    };

    struct BToKstar {
        static constexpr const char * label = "B->K^*";
        static constexpr double mB = 5.279;
        static constexpr double mV = 0.896;
        static constexpr double mR2_0m = 5.336 * 5.336;
        static constexpr double mR2_1m = 5.412 * 5.412;
        static constexpr double mR2_1p = 5.829 * 5.829;
    };

    struct BToRho {
        static constexpr const char * label = "B->rho";
        static constexpr double mB = 5.279;
        static constexpr double mV = 0.7751;
        static constexpr double mR2_0m = 5.279 * 5.279;
        static constexpr double mR2_1m = 5.325 * 5.325;
        static constexpr double mR2_1p = 5.724 * 5.724;
    };

    struct BsToPhi {
        static constexpr const char * label = "B_s->phi";
        static constexpr double mB = 5.336;
        static constexpr double mV = 1.020;
        static constexpr double mR2_0m = 5.336 * 5.336;
        static constexpr double mR2_1m = 5.412 * 5.412;
        static constexpr double mR2_1p = 5.829 * 5.829;
    };

    struct BsToKstar {
        static constexpr const char * label = "B_s->K^*";
        static constexpr double mB = 5.336;
        static constexpr double mV = 0.896;
        static constexpr double mR2_0m = 5.279 * 5.279;
        static constexpr double mR2_1m = 5.325 * 5.325;
        static constexpr double mR2_1p = 5.723 * 5.723;
    };

    template <typename Process_> class BZ2004FormFactors<Process_, PToV> :
        public FormFactors<PToV>
    {
        private:
            UsedParameter _v_factor,
                          _a0_factor, _a1_factor, _a2_factor,
                          _t1_factor, _t2_factor, _t3_factor;

            // fit parametrisation for P -> V according to [BZ2004]
            static const double
               _v_r1, _v_r2, _v_m2r, _v_m2fit,
               _a0_r1, _a0_r2, _a0_m2r, _a0_m2fit,
               _a1_r2, _a1_m2fit,
               _a2_r1, _a2_r2, _a2_m2fit,
               _t1_r1, _t1_r2, _t1_m2r, _t1_m2fit,
               _t2_r2, _t2_m2fit,
               _t3t_r1, _t3t_r2, _t3t_m2fit;

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
                _a2_factor(p["B->K^*::a2_uncertainty@BZ2004"], *this),
                _t1_factor(p["B->K^*::t1_uncertainty@BZ2004"], *this),
                _t2_factor(p["B->K^*::t2_uncertainty@BZ2004"], *this),
                _t3_factor(p["B->K^*::t3_uncertainty@BZ2004"], *this)
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

            virtual double t_1(const double & s) const
            {
                return _t1_factor * _calc_eq59(s, _t1_r1, _t1_r2, _t1_m2r, _t1_m2fit);
            }

            virtual double t_2(const double & s) const
            {
                return _t2_factor * _calc_eq61(s, _t2_r2, _t2_m2fit);
            }

            virtual double t_3(const double & s) const
            {
                const double mB = Process_::mB, mB2 = mB * mB;
                const double mV = Process_::mV, mV2 = mV * mV;

                // cf. [BZ2004], Eq. (8), p.4
                return (mB2 - mV2) / s
                   * (_t3_factor * _calc_eq60(s, _t3t_r1, _t3t_r2, _t3t_m2fit) - this->t_2(s));
            }

            virtual double t_23(const double & s) const
            {
                const double mB = Process_::mB, mB2 = mB * mB;
                const double mV = Process_::mV, mV2 = mV * mV;
                const double lambda = eos::lambda(mB2, mV2, s);

                return ((mB2 - mV2) * (mB2 + 3.0 * mV2 - s) * this->t_2(s)
                        - lambda * this->t_3(s)) / (8.0 * mB * mV2 * (mB - mV));
            }
    };

    template <> class KMPW2010FormFactors<PToV> :
        public FormFactors<PToV>
    {
        private:
            // fit parametrisation for P -> V according to [KMPW2010]
            UsedParameter
               _f0_V, _b1_V,
               _f0_A0, _b1_A0, _f0_A1, _b1_A1, _f0_A2, _b1_A2,
               _f0_T1, _b1_T1, _f0_T2, _b1_T2, _f0_T3, _b1_T3;

            static const double _tau_p, _tau_m, _tau_0;
            static const double _m_B, _m_Kstar, _m_Bs2_0m, _m_Bs2_1m, _m_Bs2_1p;

            static double _calc_z(const double & s)
            {
                return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
            }

            static double ff_KMPW(const double & s, const double & f0, const double & b1, const double & m2)
            {
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                // cf. [KMPW2010], Eq. (8.8), p. 30
                return f0 / (1.0 - s / m2) * (1.0 + b1 * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

        // Remark: this is hardcoded B->K* FF's nuisance parameters (OK, because in [KMPW2010] only B->K^* calculated)
        public:
            KMPW2010FormFactors(const Parameters & p, const Options &) :
                _f0_V(p["B->K^*::F^V(0)@KMPW2010"],   *this),   _b1_V(p["B->K^*::b^V_1@KMPW2010"],    *this),
                _f0_A0(p["B->K^*::F^A0(0)@KMPW2010"], *this),   _b1_A0(p["B->K^*::b^A0_1@KMPW2010"],  *this),
                _f0_A1(p["B->K^*::F^A1(0)@KMPW2010"], *this),   _b1_A1(p["B->K^*::b^A1_1@KMPW2010"],  *this),
                _f0_A2(p["B->K^*::F^A2(0)@KMPW2010"], *this),   _b1_A2(p["B->K^*::b^A2_1@KMPW2010"],  *this),
                _f0_T1(p["B->K^*::F^T1(0)@KMPW2010"], *this),   _b1_T1(p["B->K^*::b^T1_1@KMPW2010"],  *this),
                _f0_T2(p["B->K^*::F^T2(0)@KMPW2010"], *this),   _b1_T2(p["B->K^*::b^T2_1@KMPW2010"],  *this),
                _f0_T3(p["B->K^*::F^T3(0)@KMPW2010"], *this),   _b1_T3(p["B->K^*::b^T3_1@KMPW2010"],  *this)
            {
            }

            static FormFactors<PToV> * make(const Parameters & parameters, unsigned)
            {
                return new KMPW2010FormFactors(parameters, Options());
            }

            virtual double v(const double & s) const
            {
                return ff_KMPW(s, _f0_V(), _b1_V(), _m_Bs2_1m);
            }

            virtual double a_0(const double & s) const
            {
                return ff_KMPW(s, _f0_A0(), _b1_A0(), _m_Bs2_0m);
            }

            virtual double a_1(const double & s) const
            {
                return ff_KMPW(s, _f0_A1(), _b1_A1(), _m_Bs2_1p);
            }

            virtual double a_2(const double & s) const
            {
                return ff_KMPW(s, _f0_A2(), _b1_A2(), _m_Bs2_1p);
            }

            virtual double a_12(const double & s) const
            {
                const double mB = BToKstar::mB, mB2 = mB * mB;
                const double mV = BToKstar::mV, mV2 = mV * mV;
                const double lambda = eos::lambda(mB2, mV2, s);

                return ((mB + mV) * (mB + mV) * (mB2 - mV2 - s) * this->a_1(s)
                    - lambda * this->a_2(s)) / (16.0 * mB * mV2 * (mB + mV));
            }

            virtual double t_1(const double & s) const
            {
                return ff_KMPW(s, _f0_T1(), _b1_T1(), _m_Bs2_1m);
            }

            virtual double t_2(const double & s) const
            {
                return ff_KMPW(s, _f0_T2(), _b1_T2(), _m_Bs2_1p);
            }

            virtual double t_3(const double & s) const
            {
                return ff_KMPW(s, _f0_T3(), _b1_T3(), _m_Bs2_1p);
            }

            virtual double t_23(const double & s) const
            {
                const double mB = BToKstar::mB, mB2 = mB * mB;
                const double mV = BToKstar::mV, mV2 = mV * mV;
                const double lambda = eos::lambda(mB2, mV2, s);

                return ((mB2 - mV2) * (mB2 + 3.0 * mV2 - s) * this->t_2(s)
                        - lambda * this->t_3(s)) / (8.0 * mB * mV2 * (mB - mV));
            }
    };

    template <typename Tag_> class BFW2010FormFactors<Tag_, PToV> :
        public FormFactors<PToV>
    {
        private:
            // fit parametrisation for B_q -> Kstar according to [BFW2011]. We use the simple series expansion and
            // the results form LCSR only.
            UsedParameter _beta_V0_0, _beta_V0_1;
            UsedParameter _beta_V1_0, _beta_V1_1;
            UsedParameter _beta_V2_0, _beta_V2_1;
            UsedParameter _beta_Vt_0,_beta_Vt_1;
            static constexpr const double _m_B = Tag_::mB, _m_V = Tag_::mV;
            static constexpr const double _tau_p = (_m_B + _m_V) * (_m_B + _m_V);
            static constexpr const double _tau_m = (_m_B - _m_V) * (_m_B - _m_V);
            static constexpr const double _m_R2_0m = Tag_::mR2_0m;
            static constexpr const double _m_R2_1m = Tag_::mR2_1m;
            static constexpr const double _m_R2_1p = Tag_::mR2_1p;

            // std::sqrt is not declared as constexpr, which prohibits static
            // initialization with more strict compiler than g++ (i.e., clang).
            static double _tau_0()
            {
                return _tau_p - std::sqrt(_tau_p * _tau_p - _tau_m * _tau_p);
            }

            static double _calc_z(const double & s)
            {
                return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0())) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0()));
            }

        public:
            BFW2010FormFactors(const Parameters & p, const Options &) :
                _beta_V0_0(p[std::string(Tag_::label) + "::beta^V0_0@BFW2010"], *this),
                _beta_V0_1(p[std::string(Tag_::label) + "::beta^V0_1@BFW2010"], *this),
                _beta_V1_0(p[std::string(Tag_::label) + "::beta^V1_0@BFW2010"], *this),
                _beta_V1_1(p[std::string(Tag_::label) + "::beta^V1_1@BFW2010"], *this),
                _beta_V2_0(p[std::string(Tag_::label) + "::beta^V2_0@BFW2010"], *this),
                _beta_V2_1(p[std::string(Tag_::label) + "::beta^V2_1@BFW2010"], *this),
                _beta_Vt_0(p[std::string(Tag_::label) + "::beta^Vt_0@BFW2010"], *this),
                _beta_Vt_1(p[std::string(Tag_::label) + "::beta^Vt_1@BFW2010"], *this)
            {
            }

            static FormFactors<PToV> * make(const Parameters & parameters, unsigned)
            {
                return new BFW2010FormFactors(parameters, Options());
            }

            virtual double a_2(const double & s) const
            {
                // cf. [BFW2010], Eq. (44), p. 16, replacements Eq. (45), p. 16 and Eq. (11), p. 5 
                return (_m_B * (_m_B + _m_V)) / ((_tau_m-s) * (_tau_p-s)) * 1.0 / (1.0 - s / _m_R2_1p) * ((_m_B * _m_B -_m_V * _m_V - s) / sqrt(2)
                        * (_beta_V2_0 + _beta_V2_1 * _calc_z(s)) - (2 * _m_B * _m_V) 
                        * (_beta_V0_0 + _beta_V0_1 * _calc_z(s)));
            }

            virtual double v(const double & s) const
            {
                // cf. [BFW2010], Eq. (44), p. 16, replacements Eq. (45), p. 16 and Eq. (11), p. 5
                static const double prefactor = (_m_B + _m_V)/(_m_B * sqrt(2));
                return prefactor * 1.0 / (1.0 - s / _m_R2_1m) * (_beta_V1_0 + _beta_V1_1 * _calc_z(s));
            }

            virtual double a_1(const double & s) const
            {
                // cf. [BFW2010], Eq. (44), p. 16, replacements Eq. (45), p. 16 and Eq. (11), p. 5
                static const double prefactor = _m_B / (sqrt(2) * (_m_B + _m_V));
                return prefactor * 1.0 / (1.0 - s / _m_R2_1p) * (_beta_V2_0 + _beta_V2_1 * _calc_z(s));
            }

            virtual double a_0(const double & s) const
            {
                // cf. [BFW2010], Eq. (44), p. 16, replacements Eq. (45), p. 16 and Eq. (11), p. 5
                return 1.0 / (1.0 - s / _m_R2_0m) * (_beta_Vt_0 + _beta_Vt_1 * _calc_z(s));
            }

            virtual double a_12(const double & s) const
            {
                const double mB = BToKstar::mB, mB2 = mB * mB;
                const double mV = BToKstar::mV, mV2 = mV * mV;
                const double lambda = eos::lambda(mB2, mV2, s);

                return ((mB + mV) * (mB + mV) * (mB2 - mV2 - s) * this->a_1(s)
                    - lambda * this->a_2(s)) / (16.0 * mB * mV2 * (mB + mV));
            }

            virtual double t_1(const double &) const
            {
                throw InternalError("BFW2010FormFactors<>::t_1: Tensor form factors not yet implemented");
            }

            virtual double t_2(const double &) const
            {
                throw InternalError("BFW2010FormFactors<>::t_2: Tensor form factors not yet implemented");
            }

            virtual double t_3(const double &) const
            {
                throw InternalError("BFW2010FormFactors<>::t_3: Tensor form factors not yet implemented");
            }

            virtual double t_23(const double &) const
            {
                throw InternalError("BFW2010FormFactors<>::t_23: Tensor form factors not yet implemented");
            }
    };

    template <typename Process_> class FMvD2015FormFactors :
        public FormFactors<PToV>
    {
        private:
            UsedParameter _Fpara_0, _Fpara_beta1;
            UsedParameter _Fperp_0, _Fperp_beta1;
            UsedParameter _Flong_0;
            UsedParameter _Ftime_0, _Ftime_beta1;

            double _lambda(const double & s) const
            {
                static const double mB2 = Process_::mB * Process_::mB, mB4 = mB2 * mB2;
                static const double mV2 = Process_::mV * Process_::mV, mV4 = mV2 * mV2;

                return s * s + mB4 + mV4 - 2.0 * (s * (mB2 + mV2) + mB2 * mV2);
            }

            double _t_plus() const
            {
                return power_of<2>(Process_::mB + Process_::mV);
            }

            double _t_minus() const
            {
                return power_of<2>(Process_::mB - Process_::mV);
            }

            double _t_zero() const
            {
                return _t_plus() - std::sqrt((_t_plus() - _t_minus()) * _t_plus());
            }

            double _z(const double & t) const
            {
                return (std::sqrt(_t_plus() - t) - std::sqrt(_t_plus() - _t_zero())) / (std::sqrt(_t_plus() - t) + std::sqrt(_t_plus() - _t_zero()));
            }

            double _Flong_beta1() const
            {
                return (1.0 - _Fpara_0 / _Flong_0 * std::sqrt(_t_minus() / (2.0 * Process_::mB * Process_::mB)) * (1.0 + _Fpara_beta1 * (_z(_t_minus()) - _z(0.0)))) / (_z(0) - _z(_t_minus()));
            }

        public:
            FMvD2015FormFactors(const Parameters & p, const Options &) :
                _Fpara_0(p[std::string(Process_::label) + "::F_para(0)@FMvD2015"], *this),
                _Fpara_beta1(p[std::string(Process_::label) + "::beta_para^1@FMvD2015"], *this),
                _Fperp_0(p[std::string(Process_::label) + "::F_perp(0)@FMvD2015"], *this),
                _Fperp_beta1(p[std::string(Process_::label) + "::beta_perp^1@FMvD2015"], *this),
                _Flong_0(p[std::string(Process_::label) + "::F_long(0)@FMvD2015"], *this),
                _Ftime_0(p[std::string(Process_::label) + "::F_time(0)@FMvD2015"], *this),
                _Ftime_beta1(p[std::string(Process_::label) + "::beta_time^1@FMvD2015"], *this)
            {
            }

            ~FMvD2015FormFactors()
            {
            }

            static FormFactors<PToV> * make(const Parameters & parameters, unsigned)
            {
                return new FMvD2015FormFactors(parameters, Options());
            }

            virtual double f_long(const double & s) const
            {
                double zs = _z(s), z0 = _z(0.0);

                return _Flong_0 / (1.0 - s / Process_::mR2_1p) * (1.0 + _Flong_beta1() * (zs - z0));
            }

            virtual double f_perp(const double & s) const
            {
                static const double mB = Process_::mB, mB2 = mB * mB;
                static const double mV = Process_::mV, mV2 = mV * mV;

                double zs = _z(s), z0 = _z(0.0);
                double kin = std::sqrt(_lambda(s)) / (mB2 - mV2);

                return _Fperp_0 * kin / (1.0 - s / Process_::mR2_1m) * (1.0 + _Fperp_beta1 * (zs - z0));
            }

            virtual double f_para(const double & s) const
            {
                double zs = _z(s), z0 = _z(0.0);

                return _Fpara_0 / (1.0 - s / Process_::mR2_1p) * (1.0 + _Fpara_beta1 * (zs - z0));
            }

            virtual double f_time(const double & s) const
            {
                static const double mB = Process_::mB, mB2 = mB * mB;
                static const double mV = Process_::mV, mV2 = mV * mV;

                double zs = _z(s), z0 = _z(0.0);
                double kin = std::sqrt(_lambda(s)) / (mB2 - mV2);

                return _Ftime_0 * kin / (1.0 - s / Process_::mR2_0m) * (1.0 + _Ftime_beta1 * (zs - z0));
            }

            virtual double v(const double & s) const
            {
                static const double mB = Process_::mB;
                static const double mV = Process_::mV;

                return this->f_perp(s) * mB * (mB + mV) / std::sqrt(2.0 * _lambda(s));
            }

            virtual double a_0(const double & s) const
            {
                static const double mB = Process_::mB, mB2 = mB * mB;

                return mB2 / std::sqrt(_lambda(s)) * f_time(s);
            }

            virtual double a_1(const double & s) const
            {
                return this->f_para(s) * Process_::mB / std::sqrt(2.0) / (Process_::mB + Process_::mV);
            }

            virtual double a_2(const double & s) const
            {
                static const double mB = Process_::mB, mB2 = mB * mB;
                static const double mV = Process_::mV, mV2 = mV * mV;

                return ((mB + mV) * mB / std::sqrt(2.0) * (mB2 - mV2 - s) * this->f_para(s) - 2.0 * mV * mB2 * (mB + mV) * this->f_long(s)) / _lambda(s);
            }

            virtual double a_12(const double & s) const
            {
                return this->f_long(s) * Process_::mB / (8.0 * Process_::mV);
            }

            virtual double t_1(const double &) const
            {
                throw InternalError("FMvD2015FormFactors<>::t_1: Tensor form factors not yet implemented");
            }

            virtual double t_2(const double &) const
            {
                throw InternalError("FMvD2015FormFactors<>::t_2: Tensor form factors not yet implemented");
            }

            virtual double t_3(const double &) const
            {
                throw InternalError("FMvD2015FormFactors<>::t_3: Tensor form factors not yet implemented");
            }

            virtual double t_23(const double &) const
            {
                throw InternalError("FMvD2015FormFactors<>::t_23: Tensor form factors not yet implemented");
            }
    };

    template <typename Process_> class BSZ2015FormFactors<Process_, PToV> :
        public FormFactors<PToV>
    {
        private:
            // fit parametrization for P -> V according to [BSZ2015]
            std::array<UsedParameter, 3> _a_A0, _a_A1, _a_V, _a_T1, _a_T23;
            // use constraint (B.6) in [BSZ2015] to remove A_12(0)
            std::array<UsedParameter, 2> _a_A12, _a_T2;

            const double _mB, _mB2, _mV, _mV2, _kin_factor;
            const double _tau_p, _tau_0;
            const double _z_0;

            static double _calc_tau_0(const double & m_B, const double & m_V)
            {
                const double tau_p = power_of<2>(m_B + m_V);
                const double tau_m = power_of<2>(m_B - m_V);
                return tau_p * (1.0 - std::sqrt(1.0 - tau_m / tau_p));
            }

            double _calc_z(const double & s) const
            {
                return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
            }

            template <typename Parameter_>
            double _calc_ff(const double & s, const double & m2_R, const std::array<Parameter_, 3> & a) const
            {
                const double diff_z = _calc_z(s) - _z_0;
                return 1.0 / (1.0 - s / m2_R) *
                       (a[0] + a[1] * diff_z + a[2] * power_of<2>(diff_z));
            }

            static std::string _par_name(const std::string & ff_name)
            {
                return std::string(Process_::label) + std::string("::alpha^") + ff_name + std::string("@BSZ2015");
            }

        public:
            BSZ2015FormFactors(const Parameters & p, const Options &) :
                _a_A0{{  UsedParameter(p[_par_name("A0_0")],  *this),
                         UsedParameter(p[_par_name("A0_1")],  *this),
                         UsedParameter(p[_par_name("A0_2")],  *this) }},
                _a_A1{{  UsedParameter(p[_par_name("A1_0")],  *this),
                         UsedParameter(p[_par_name("A1_1")],  *this),
                         UsedParameter(p[_par_name("A1_2")],  *this) }},
                _a_V{{   UsedParameter(p[_par_name("V_0")],   *this),
                         UsedParameter(p[_par_name("V_1")],   *this),
                         UsedParameter(p[_par_name("V_2")],   *this) }},
                _a_T1{{  UsedParameter(p[_par_name("T1_0")],  *this),
                         UsedParameter(p[_par_name("T1_1")],  *this),
                         UsedParameter(p[_par_name("T1_2")],  *this) }},
                _a_T23{{ UsedParameter(p[_par_name("T23_0")], *this),
                         UsedParameter(p[_par_name("T23_1")], *this),
                         UsedParameter(p[_par_name("T23_2")], *this) }},
                _a_A12{{ UsedParameter(p[_par_name("A12_1")], *this),
                         UsedParameter(p[_par_name("A12_2")], *this) }},
                _a_T2{{  UsedParameter(p[_par_name("T2_1")],  *this),
                         UsedParameter(p[_par_name("T2_2")],  *this) }},
                _mB(Process_::mB),
                _mB2(power_of<2>(_mB)),
                _mV(Process_::mV),
                _mV2(power_of<2>(_mV)),
                _kin_factor((_mB2 - _mV2) / (8.0 * _mB * _mV)),
                _tau_p(power_of<2>(_mB + _mV)),
                _tau_0(_calc_tau_0(_mB, _mV)),
                _z_0(_calc_z(0))
            {
            }

            ~BSZ2015FormFactors()
            {
            }

            static FormFactors<PToV> * make(const Parameters & parameters, unsigned)
            {
                return new BSZ2015FormFactors(parameters, Options());
            }

            virtual double v(const double & s) const
            {
                return _calc_ff(s, Process_::mR2_1m, _a_V);
            }

            virtual double a_0(const double & s) const
            {
                return _calc_ff(s, Process_::mR2_0m, _a_A0);
            }

            virtual double a_1(const double & s) const
            {
                return _calc_ff(s, Process_::mR2_1p, _a_A1);
            }

            virtual double a_2(const double & s) const
            {
                const double lambda = eos::lambda(_mB2, _mV2, s);

                return (power_of<2>(_mB + _mV) * (_mB2 - _mV2 - s) * a_1(s)
                        - 16.0 * _mB * _mV2 * (_mB + _mV) * a_12(s)) / lambda;
            }

            virtual double a_12(const double & s) const
            {
                // use constraint (B.6) in [BSZ2015] to remove A_12(0)
                std::array<double, 3> values
                {{
                    _kin_factor * _a_A0[0],
                    _a_A12[1 - 1],
                    _a_A12[2 - 1],
                }};

                return _calc_ff(s, Process_::mR2_1p, values);
            }

            virtual double t_1(const double & s) const
            {
                return _calc_ff(s, Process_::mR2_1m, _a_T1);
            }

            virtual double t_2(const double & s) const
            {
                // use constraint T_1(0) = T_2(0) to replace T_2(0)
                std::array<double, 3> values
                {{
                    _a_T1[0],
                    _a_T2[1 - 1],
                    _a_T2[2 - 1],
                }};
                return _calc_ff(s, Process_::mR2_1p, values);
            }

            virtual double t_3(const double & s) const
            {
                const double lambda = eos::lambda(_mB2, _mV2, s);

                return ((_mB2 - _mV2) * (_mB2 + 3.0 * _mV2 - s) * t_2(s)
                        - 8.0 * _mB * _mV2 * (_mB - _mV) * t_23(s)) / lambda;
            }

            virtual double t_23(const double & s) const
            {
                return _calc_ff(s, Process_::mR2_1p, _a_T23);
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
        static constexpr const double m2_Br1m = 5.415 * 5.415; // B_s^*
        static constexpr const double m2_Br0p = 5.630 * 5.630; // B_s scalar
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct BToPi {
        typedef PToP Transition;
        static constexpr const char * label = "B->pi";
        static constexpr const double m_B = 5.279;
        static constexpr const double m_P = 0.135;
        static constexpr const double m2_Br1m = 5.325 * 5.325; // B_{u,d}^*
        static constexpr const double m2_Br0p = 5.540 * 5.540; // B_{u,d} scalar: M(B_s scalar) - M(B_s^*) + M(B_{u,d}^*)
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = true;
    };

    struct BToD {
        typedef PToP Transition;
        static constexpr const char * label = "B->D";
        static constexpr const double m_B = 5.279;
        static constexpr const double m_P = 1.870;
        // resonance masses from [HPQCD2015A]
        static constexpr const double m2_Br1m = 6.330 * 6.330; // B_c^*
        static constexpr const double m2_Br0p = 6.420 * 6.420; // B_c scalar
        static constexpr const double tau_p = (m_B + m_P) * (m_B + m_P);
        static constexpr const double tau_m = (m_B - m_P) * (m_B - m_P);
        static constexpr const bool uses_tensor_form_factors = false;
    };

    /*
     * P -> P Form Factors in the simplified series expansion according to
     * [BCL2008].
     */
    template <typename Process_, bool with_tensor_> class BCL2008FormFactorBase;

    template <typename Process_> class BCL2008FormFactorBase<Process_, false> :
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

        protected:
            double _z(const double & s) const
            {
                static const double m_B = Process_::m_B;
                static const double m_P = Process_::m_P;
                static const double tau_p = Process_::tau_p;
                static const double tau_0 = (m_B + m_P) * (std::sqrt(m_B) - std::sqrt(m_P)) * (std::sqrt(m_B) - std::sqrt(m_P));

                return (std::sqrt(tau_p - s) - std::sqrt(tau_p - tau_0))
                    / (std::sqrt(tau_p - s) + std::sqrt(tau_p - tau_0));
            }

        public:
            BCL2008FormFactorBase(const Parameters & p, const Options &) :
                _f_plus_0(p[std::string(Process_::label) + "::f_+(0)@BCL2008"], *this),
                _b_plus_1(p[std::string(Process_::label) + "::b_+^1@BCL2008"],  *this),
                _b_plus_2(p[std::string(Process_::label) + "::b_+^2@BCL2008"],  *this),
                _b_zero_1(p[std::string(Process_::label) + "::b_0^1@BCL2008"],  *this),
                _b_zero_2(p[std::string(Process_::label) + "::b_0^2@BCL2008"],  *this)
            {
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
                const double z = _z(s), z2 = z * z;
                const double z0 = _z(0), z02 = z0 * z0;
                const double zbar = z - z0, z2bar = z2 - z02;

                // note that f_0(0) = f_+(0)!
                // for f_0(s) we do not have an equation of motion to express _b_zero_K in terms of the
                // other coefficients!
                return _f_plus_0 / (1.0 - s / Process_::m2_Br0p) * (1.0 + _b_zero_1 * zbar + _b_zero_2 * z2bar);
            }

            virtual double f_t(const double &) const
            {
                throw InternalError("This form factor parametrization has no inputs for tensor form factors.");

                return 0.0;
            }
    };

    template <typename Process_> class BCL2008FormFactorBase<Process_, true> :
        public BCL2008FormFactorBase<Process_, false>
    {
        private:
            /*
             * Fit parametrisation for P -> P according to [BCL2008], eq. (11)
             * with K = 3. Note that we factor out the form factor at q^2 = 0
             * by setting t_0 = 0.0, thus b_k -> b_k / b_0. Note that the last
             * coefficient b_K is fixed by eq. (14).
             *
             * Tensor form factors only. Vector and scalar form factor in BCL2008FormFactorBase<>.
             */
            UsedParameter _f_t_0,    _b_t_1,    _b_t_2;

        public:
            BCL2008FormFactorBase(const Parameters & p, const Options & o) :
                BCL2008FormFactorBase<Process_, false>(p, o),
                _f_t_0(p[std::string(Process_::label)    + "::f_T(0)@BCL2008"], *this),
                _b_t_1(p[std::string(Process_::label)    + "::b_T^1@BCL2008"],  *this),
                _b_t_2(p[std::string(Process_::label)    + "::b_T^2@BCL2008"],  *this)
            {
            }

            virtual double f_t(const double & s) const
            {
                const double z = this->_z(s), z2 = z * z, z3 = z * z2;
                const double z0 = this->_z(0), z02 = z0 * z0, z03 = z0 * z02;
                const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03;

                return _f_t_0 / (1.0 - s / Process_::m2_Br1m) * (1.0 + _b_t_1 * (zbar - z3bar / 3.0) + _b_t_2 * (z2bar + 2.0 * z3bar / 3.0));
            }
    };

    template <typename Process_> class BCL2008FormFactors :
        public BCL2008FormFactorBase<Process_, Process_::uses_tensor_form_factors>
    {
        public:
            BCL2008FormFactors(const Parameters & p, const Options & o) :
                BCL2008FormFactorBase<Process_, Process_::uses_tensor_form_factors>(p, o)
            {
            }

            static FormFactors<PToP> * make(const Parameters & parameters, unsigned)
            {
                return new BCL2008FormFactors(parameters, Options());
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

    template <> class BFW2010FormFactors<BToK, PToP> :
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

    /* P -> PP Processes */

    struct BToPiPi {
        typedef PToPP Transition;
        static constexpr const char * label = "B->pipi";
        static constexpr double mB  = 5.2795;
        static constexpr double mP1 = 0.13957;
        static constexpr double mP2 = 0.13957;

        // for pole and t_0 calculation in zhat
        static constexpr double mBst = 5.32465;

        // for pole calculation in z, depending on the current at hand
        static constexpr double mR2_1m = 5.32465;
        static constexpr double mR2_1p = 5.72590;
        static constexpr double mR2_0m = 5.27932;
    };

    template <typename Process_> class FvDV2018FormFactors :
        public FormFactors<PToPP>
    {
        private:
            UsedParameter _a_Fperp_0_0, _a_Fperp_0_1, _a_Fperp_0_2, _a_Fperp_0_3, _a_Fperp_1_0, _a_Fperp_1_1, _a_Fperp_1_2;
            UsedParameter _b_Fperp_0_0, _b_Fperp_0_1, _b_Fperp_0_2, _b_Fperp_0_3, _b_Fperp_1_0, _b_Fperp_1_1, _b_Fperp_1_2;
            UsedParameter _c_Fperp_0_0, _c_Fperp_0_1, _c_Fperp_0_2, _c_Fperp_0_3, _c_Fperp_1_0, _c_Fperp_1_1, _c_Fperp_1_2;

            UsedParameter _a_Fpara_0_0, _a_Fpara_0_1, _a_Fpara_0_2, _a_Fpara_0_3, _a_Fpara_1_0, _a_Fpara_1_1, _a_Fpara_1_2;
            UsedParameter _b_Fpara_0_0, _b_Fpara_0_1, _b_Fpara_0_2, _b_Fpara_0_3, _b_Fpara_1_0, _b_Fpara_1_1, _b_Fpara_1_2;
            UsedParameter _c_Fpara_0_0, _c_Fpara_0_1, _c_Fpara_0_2, _c_Fpara_0_3, _c_Fpara_1_0, _c_Fpara_1_1, _c_Fpara_1_2;

            UsedParameter _a_Flong_0_0, _a_Flong_0_1, _a_Flong_0_2, _a_Flong_0_3, _a_Flong_1_0, _a_Flong_1_1, _a_Flong_1_2;
            UsedParameter _b_Flong_0_0, _b_Flong_0_1, _b_Flong_0_2, _b_Flong_0_3, _b_Flong_1_0, _b_Flong_1_1, _b_Flong_1_2;
            UsedParameter _c_Flong_0_0, _c_Flong_0_1, _c_Flong_0_2, _c_Flong_0_3, _c_Flong_1_0, _c_Flong_1_1, _c_Flong_1_2;

            UsedParameter _a_Ftime_0_0, _a_Ftime_0_1, _a_Ftime_0_2, _a_Ftime_0_3, _a_Ftime_1_0, _a_Ftime_1_1, _a_Ftime_1_2;
            UsedParameter _b_Ftime_0_0, _b_Ftime_0_1, _b_Ftime_0_2, _b_Ftime_0_3, _b_Ftime_1_0, _b_Ftime_1_1, _b_Ftime_1_2;
            UsedParameter _c_Ftime_0_0, _c_Ftime_0_1, _c_Ftime_0_2, _c_Ftime_0_3, _c_Ftime_1_0, _c_Ftime_1_1, _c_Ftime_1_2;

            static double _calc_z(const double & t, const double & t_p, const double & t_0)
            {
                return (std::sqrt(t_p - t) - std::sqrt(t_p - t_0)) / (std::sqrt(t_p - t) + std::sqrt(t_p - t_0));
            }

            double _z(const double & t) const
            {
                static constexpr double mB  = Process_::mB;
                static constexpr double mP1 = Process_::mP1;
                static constexpr double mP2 = Process_::mP2;

                static constexpr double t_p = power_of<2>(mB + mP1 + mP2);
                static constexpr double t_0 = 0.0;

                return _calc_z(t, t_p, t_0);
            }

            double _zhat(const double & that) const
            {
                static constexpr double mB    = Process_::mB;
                static constexpr double mP2   = Process_::mP2;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                static constexpr double that_p = power_of<2>(mB + mP2);
                static const     double that_0 = that_p - std::sqrt(that_p * (that_p - mBst2));

                return _calc_z(that, that_p, that_0);
            }

            double _blaschke(const double & z, const double & zh) const
            {
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                const double zBst2  = _z(mBst2);
                const double zhBst2 = _zhat(mBst2);

                double result  = 1.0;
                result *= (1.0 - z  *  zBst2) / (z  -  zBst2);
                result *= (1.0 - zh * zhBst2) / (zh - zhBst2);

                return result;
            }

            double _blaschke_res_qhat2(const double & z) const
            {
                static constexpr double mB    = Process_::mB;
                static constexpr double mP2   = Process_::mP2;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                static constexpr double that_p = power_of<2>(mB + mP2);

                const double zBst2  = _z(mBst2);

                double result  = 4.0 * (mBst2 - that_p);
                result *= (1.0 - z  *  zBst2) / (z  -  zBst2);

                return result;
            }

          public:
            FvDV2018FormFactors(const Parameters & p, const Options &) :
                // perp
                _a_Fperp_0_0(p["B->pipi::a^Fperp_0_0@FvDV2018"], *this),
                _a_Fperp_0_1(p["B->pipi::a^Fperp_0_1@FvDV2018"], *this),
                _a_Fperp_0_2(p["B->pipi::a^Fperp_0_2@FvDV2018"], *this),
                _a_Fperp_0_3(p["B->pipi::a^Fperp_0_3@FvDV2018"], *this),
                _a_Fperp_1_0(p["B->pipi::a^Fperp_1_0@FvDV2018"], *this),
                _a_Fperp_1_1(p["B->pipi::a^Fperp_1_1@FvDV2018"], *this),
                _a_Fperp_1_2(p["B->pipi::a^Fperp_1_2@FvDV2018"], *this),
                _b_Fperp_0_0(p["B->pipi::b^Fperp_0_0@FvDV2018"], *this),
                _b_Fperp_0_1(p["B->pipi::b^Fperp_0_1@FvDV2018"], *this),
                _b_Fperp_0_2(p["B->pipi::b^Fperp_0_2@FvDV2018"], *this),
                _b_Fperp_0_3(p["B->pipi::b^Fperp_0_3@FvDV2018"], *this),
                _b_Fperp_1_0(p["B->pipi::b^Fperp_1_0@FvDV2018"], *this),
                _b_Fperp_1_1(p["B->pipi::b^Fperp_1_1@FvDV2018"], *this),
                _b_Fperp_1_2(p["B->pipi::b^Fperp_1_2@FvDV2018"], *this),
                _c_Fperp_0_0(p["B->pipi::c^Fperp_0_0@FvDV2018"], *this),
                _c_Fperp_0_1(p["B->pipi::c^Fperp_0_1@FvDV2018"], *this),
                _c_Fperp_0_2(p["B->pipi::c^Fperp_0_2@FvDV2018"], *this),
                _c_Fperp_0_3(p["B->pipi::c^Fperp_0_3@FvDV2018"], *this),
                _c_Fperp_1_0(p["B->pipi::c^Fperp_1_0@FvDV2018"], *this),
                _c_Fperp_1_1(p["B->pipi::c^Fperp_1_1@FvDV2018"], *this),
                _c_Fperp_1_2(p["B->pipi::c^Fperp_1_2@FvDV2018"], *this),
                // para
                _a_Fpara_0_0(p["B->pipi::a^Fpara_0_0@FvDV2018"], *this),
                _a_Fpara_0_1(p["B->pipi::a^Fpara_0_1@FvDV2018"], *this),
                _a_Fpara_0_2(p["B->pipi::a^Fpara_0_2@FvDV2018"], *this),
                _a_Fpara_0_3(p["B->pipi::a^Fpara_0_3@FvDV2018"], *this),
                _a_Fpara_1_0(p["B->pipi::a^Fpara_1_0@FvDV2018"], *this),
                _a_Fpara_1_1(p["B->pipi::a^Fpara_1_1@FvDV2018"], *this),
                _a_Fpara_1_2(p["B->pipi::a^Fpara_1_2@FvDV2018"], *this),
                _b_Fpara_0_0(p["B->pipi::b^Fpara_0_0@FvDV2018"], *this),
                _b_Fpara_0_1(p["B->pipi::b^Fpara_0_1@FvDV2018"], *this),
                _b_Fpara_0_2(p["B->pipi::b^Fpara_0_2@FvDV2018"], *this),
                _b_Fpara_0_3(p["B->pipi::b^Fpara_0_3@FvDV2018"], *this),
                _b_Fpara_1_0(p["B->pipi::b^Fpara_1_0@FvDV2018"], *this),
                _b_Fpara_1_1(p["B->pipi::b^Fpara_1_1@FvDV2018"], *this),
                _b_Fpara_1_2(p["B->pipi::b^Fpara_1_2@FvDV2018"], *this),
                _c_Fpara_0_0(p["B->pipi::c^Fpara_0_0@FvDV2018"], *this),
                _c_Fpara_0_1(p["B->pipi::c^Fpara_0_1@FvDV2018"], *this),
                _c_Fpara_0_2(p["B->pipi::c^Fpara_0_2@FvDV2018"], *this),
                _c_Fpara_0_3(p["B->pipi::c^Fpara_0_3@FvDV2018"], *this),
                _c_Fpara_1_0(p["B->pipi::c^Fpara_1_0@FvDV2018"], *this),
                _c_Fpara_1_1(p["B->pipi::c^Fpara_1_1@FvDV2018"], *this),
                _c_Fpara_1_2(p["B->pipi::c^Fpara_1_2@FvDV2018"], *this),
                // long
                _a_Flong_0_0(p["B->pipi::a^Flong_0_0@FvDV2018"], *this),
                _a_Flong_0_1(p["B->pipi::a^Flong_0_1@FvDV2018"], *this),
                _a_Flong_0_2(p["B->pipi::a^Flong_0_2@FvDV2018"], *this),
                _a_Flong_0_3(p["B->pipi::a^Flong_0_3@FvDV2018"], *this),
                _a_Flong_1_0(p["B->pipi::a^Flong_1_0@FvDV2018"], *this),
                _a_Flong_1_1(p["B->pipi::a^Flong_1_1@FvDV2018"], *this),
                _a_Flong_1_2(p["B->pipi::a^Flong_1_2@FvDV2018"], *this),
                _b_Flong_0_0(p["B->pipi::b^Flong_0_0@FvDV2018"], *this),
                _b_Flong_0_1(p["B->pipi::b^Flong_0_1@FvDV2018"], *this),
                _b_Flong_0_2(p["B->pipi::b^Flong_0_2@FvDV2018"], *this),
                _b_Flong_0_3(p["B->pipi::b^Flong_0_3@FvDV2018"], *this),
                _b_Flong_1_0(p["B->pipi::b^Flong_1_0@FvDV2018"], *this),
                _b_Flong_1_1(p["B->pipi::b^Flong_1_1@FvDV2018"], *this),
                _b_Flong_1_2(p["B->pipi::b^Flong_1_2@FvDV2018"], *this),
                _c_Flong_0_0(p["B->pipi::c^Flong_0_0@FvDV2018"], *this),
                _c_Flong_0_1(p["B->pipi::c^Flong_0_1@FvDV2018"], *this),
                _c_Flong_0_2(p["B->pipi::c^Flong_0_2@FvDV2018"], *this),
                _c_Flong_0_3(p["B->pipi::c^Flong_0_3@FvDV2018"], *this),
                _c_Flong_1_0(p["B->pipi::c^Flong_1_0@FvDV2018"], *this),
                _c_Flong_1_1(p["B->pipi::c^Flong_1_1@FvDV2018"], *this),
                _c_Flong_1_2(p["B->pipi::c^Flong_1_2@FvDV2018"], *this),
                // time
                _a_Ftime_0_0(p["B->pipi::a^Ftime_0_0@FvDV2018"], *this),
                _a_Ftime_0_1(p["B->pipi::a^Ftime_0_1@FvDV2018"], *this),
                _a_Ftime_0_2(p["B->pipi::a^Ftime_0_2@FvDV2018"], *this),
                _a_Ftime_0_3(p["B->pipi::a^Ftime_0_3@FvDV2018"], *this),
                _a_Ftime_1_0(p["B->pipi::a^Ftime_1_0@FvDV2018"], *this),
                _a_Ftime_1_1(p["B->pipi::a^Ftime_1_1@FvDV2018"], *this),
                _a_Ftime_1_2(p["B->pipi::a^Ftime_1_2@FvDV2018"], *this),
                _b_Ftime_0_0(p["B->pipi::b^Ftime_0_0@FvDV2018"], *this),
                _b_Ftime_0_1(p["B->pipi::b^Ftime_0_1@FvDV2018"], *this),
                _b_Ftime_0_2(p["B->pipi::b^Ftime_0_2@FvDV2018"], *this),
                _b_Ftime_0_3(p["B->pipi::b^Ftime_0_3@FvDV2018"], *this),
                _b_Ftime_1_0(p["B->pipi::b^Ftime_1_0@FvDV2018"], *this),
                _b_Ftime_1_1(p["B->pipi::b^Ftime_1_1@FvDV2018"], *this),
                _b_Ftime_1_2(p["B->pipi::b^Ftime_1_2@FvDV2018"], *this),
                _c_Ftime_0_0(p["B->pipi::c^Ftime_0_0@FvDV2018"], *this),
                _c_Ftime_0_1(p["B->pipi::c^Ftime_0_1@FvDV2018"], *this),
                _c_Ftime_0_2(p["B->pipi::c^Ftime_0_2@FvDV2018"], *this),
                _c_Ftime_0_3(p["B->pipi::c^Ftime_0_3@FvDV2018"], *this),
                _c_Ftime_1_0(p["B->pipi::c^Ftime_1_0@FvDV2018"], *this),
                _c_Ftime_1_1(p["B->pipi::c^Ftime_1_1@FvDV2018"], *this),
                _c_Ftime_1_2(p["B->pipi::c^Ftime_1_2@FvDV2018"], *this)
            {
            }

            static FormFactors<PToPP> * make(const Parameters & parameters, const Options & options)
            {
                return new FvDV2018FormFactors(parameters, options);
            }

            virtual complex<double> f_perp(const double & q2, const double & k2, const double & ctheta) const
            {
                static constexpr double mB  = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mP2 = Process_::mP2, mP22 = mP2 * mP2;

                const double lambda = eos::lambda(q2, k2, mB2);
                const double E2     = (mB2 + k2 - q2 - ctheta * std::sqrt(lambda)) / (4.0 * mB);
                const double qhat2  = mB2 + mP22 - 2.0 * mB * E2;

                const double z  = this->_z(q2);
                const double zh = this->_zhat(qhat2); 

                const double a = _a_Fperp_0_0 + _a_Fperp_1_0 * z + _a_Fperp_0_1 * zh + _a_Fperp_1_1 * z * zh + _a_Fperp_1_2 * z * zh * zh + _a_Fperp_0_2 * zh * zh + _a_Fperp_0_3 * zh * zh * zh;
                const double b = _b_Fperp_0_0 + _b_Fperp_1_0 * z + _b_Fperp_0_1 * zh + _b_Fperp_1_1 * z * zh + _b_Fperp_1_2 * z * zh * zh + _b_Fperp_0_2 * zh * zh + _b_Fperp_0_3 * zh * zh * zh;
                const double c = _c_Fperp_0_0 + _c_Fperp_1_0 * z + _c_Fperp_0_1 * zh + _c_Fperp_1_1 * z * zh + _c_Fperp_1_2 * z * zh * zh + _c_Fperp_0_2 * zh * zh + _c_Fperp_0_3 * zh * zh * zh;

                const double blaschke = this->_blaschke(z, zh);

                complex<double> result{ 0.0, blaschke * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * std::sqrt(lambda) / (mB * std::sqrt(k2)) };

                return result;
            }

            virtual double f_perp_im_res_qhat2(const double & q2, const double & k2) const
            {
                static constexpr double mB    = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                const double lambda = eos::lambda(q2, k2, mB2);
                const double z  = this->_z(q2);
                const double zh = this->_z(mBst2);

                const double a = _a_Fperp_0_0 + _a_Fperp_1_0 * z + _a_Fperp_0_1 * zh + _a_Fperp_1_1 * z * zh + _a_Fperp_1_2 * z * zh * zh + _a_Fperp_0_2 * zh * zh + _a_Fperp_0_3 * zh * zh * zh;
                const double b = _b_Fperp_0_0 + _b_Fperp_1_0 * z + _b_Fperp_0_1 * zh + _b_Fperp_1_1 * z * zh + _b_Fperp_1_2 * z * zh * zh + _b_Fperp_0_2 * zh * zh + _b_Fperp_0_3 * zh * zh * zh;
                const double c = _c_Fperp_0_0 + _c_Fperp_1_0 * z + _c_Fperp_0_1 * zh + _c_Fperp_1_1 * z * zh + _c_Fperp_1_2 * z * zh * zh + _c_Fperp_0_2 * zh * zh + _c_Fperp_0_3 * zh * zh * zh;

                const double blaschke_res_qhat2 = this->_blaschke_res_qhat2(z);

                double result = blaschke_res_qhat2 * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * std::sqrt(lambda) / (mB * std::sqrt(k2));

                return result;
            }

            virtual complex<double> f_para(const double & q2, const double & k2, const double & ctheta) const
            {
                static constexpr double mB  = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mP2 = Process_::mP2, mP22 = mP2 * mP2;

                const double lambda = eos::lambda(q2, k2, mB2);
                const double E2     = (mB2 + k2 - q2 - ctheta * std::sqrt(lambda)) / (4.0 * mB);
                const double qhat2  = mB2 + mP22 - 2.0 * mB * E2;

                const double z  = this->_z(q2);
                const double zh = this->_zhat(qhat2); 

                const double a = _a_Fpara_0_0 + _a_Fpara_1_0 * z + _a_Fpara_0_1 * zh + _a_Fpara_1_1 * z * zh + _a_Fpara_1_2 * z * zh * zh + _a_Fpara_0_2 * zh * zh + _a_Fpara_0_3 * zh * zh * zh;
                const double b = _b_Fpara_0_0 + _b_Fpara_1_0 * z + _b_Fpara_0_1 * zh + _b_Fpara_1_1 * z * zh + _b_Fpara_1_2 * z * zh * zh + _b_Fpara_0_2 * zh * zh + _b_Fpara_0_3 * zh * zh * zh;
                const double c = _c_Fpara_0_0 + _c_Fpara_1_0 * z + _c_Fpara_0_1 * zh + _c_Fpara_1_1 * z * zh + _c_Fpara_1_2 * z * zh * zh + _c_Fpara_0_2 * zh * zh + _c_Fpara_0_3 * zh * zh * zh;

                const double blaschke = this->_blaschke(z, zh);

                complex<double> result{ 0.0, blaschke * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB / std::sqrt(k2) };

                return result;
            }

            virtual double f_para_im_res_qhat2(const double & q2, const double & k2) const
            {
                static constexpr double mB    = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                const double z  = this->_z(q2);
                const double zh = this->_z(mBst2);

                const double a = _a_Fpara_0_0 + _a_Fpara_1_0 * z + _a_Fpara_0_1 * zh + _a_Fpara_1_1 * z * zh + _a_Fpara_1_2 * z * zh * zh + _a_Fpara_0_2 * zh * zh + _a_Fpara_0_3 * zh * zh * zh;
                const double b = _b_Fpara_0_0 + _b_Fpara_1_0 * z + _b_Fpara_0_1 * zh + _b_Fpara_1_1 * z * zh + _b_Fpara_1_2 * z * zh * zh + _b_Fpara_0_2 * zh * zh + _b_Fpara_0_3 * zh * zh * zh;
                const double c = _c_Fpara_0_0 + _c_Fpara_1_0 * z + _c_Fpara_0_1 * zh + _c_Fpara_1_1 * z * zh + _c_Fpara_1_2 * z * zh * zh + _c_Fpara_0_2 * zh * zh + _c_Fpara_0_3 * zh * zh * zh;

                const double blaschke_res_qhat2 = this->_blaschke_res_qhat2(z);

                double result = blaschke_res_qhat2 * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB / std::sqrt(k2);

                return result;
            }

            virtual complex<double> f_long(const double & q2, const double & k2, const double & ctheta) const
            {
                static constexpr double mB  = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mP2 = Process_::mP2, mP22 = mP2 * mP2;

                const double lambda = eos::lambda(q2, k2, mB2);
                const double E2     = (mB2 + k2 - q2 - ctheta * std::sqrt(lambda)) / (4.0 * mB);
                const double qhat2  = mB2 + mP22 - 2.0 * mB * E2;

                const double z  = this->_z(q2);
                const double zh = this->_zhat(qhat2); 

                const double a = _a_Flong_0_0 + _a_Flong_1_0 * z + _a_Flong_0_1 * zh + _a_Flong_1_1 * z * zh + _a_Flong_1_2 * z * zh * zh + _a_Flong_0_2 * zh * zh + _a_Flong_0_3 * zh * zh * zh;
                const double b = _b_Flong_0_0 + _b_Flong_1_0 * z + _b_Flong_0_1 * zh + _b_Flong_1_1 * z * zh + _b_Flong_1_2 * z * zh * zh + _b_Flong_0_2 * zh * zh + _b_Flong_0_3 * zh * zh * zh;
                const double c = _c_Flong_0_0 + _c_Flong_1_0 * z + _c_Flong_0_1 * zh + _c_Flong_1_1 * z * zh + _c_Flong_1_2 * z * zh * zh + _c_Flong_0_2 * zh * zh + _c_Flong_0_3 * zh * zh * zh;

                const double blaschke = this->_blaschke(z, zh);

                complex<double> result{ 0.0, blaschke * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB / std::sqrt(q2) * mB2 / std::sqrt(lambda) * mB2 / k2 };

                return result;
            }

            virtual double f_long_im_res_qhat2(const double & q2, const double & k2) const
            {
                static constexpr double mB    = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                const double lambda = eos::lambda(q2, k2, mB2);
                const double z  = this->_z(q2);
                const double zh = this->_z(mBst2);

                const double a = _a_Flong_0_0 + _a_Flong_1_0 * z + _a_Flong_0_1 * zh + _a_Flong_1_1 * z * zh + _a_Flong_1_2 * z * zh * zh + _a_Flong_0_2 * zh * zh + _a_Flong_0_3 * zh * zh * zh;
                const double b = _b_Flong_0_0 + _b_Flong_1_0 * z + _b_Flong_0_1 * zh + _b_Flong_1_1 * z * zh + _b_Flong_1_2 * z * zh * zh + _b_Flong_0_2 * zh * zh + _b_Flong_0_3 * zh * zh * zh;
                const double c = _c_Flong_0_0 + _c_Flong_1_0 * z + _c_Flong_0_1 * zh + _c_Flong_1_1 * z * zh + _c_Flong_1_2 * z * zh * zh + _c_Flong_0_2 * zh * zh + _c_Flong_0_3 * zh * zh * zh;

                const double blaschke_res_qhat2 = this->_blaschke_res_qhat2(z);

                double result = blaschke_res_qhat2 * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB / std::sqrt(q2) * mB2 / std::sqrt(lambda) * mB2 / k2;

                return result;
            }

            virtual complex<double> f_time(const double & q2, const double & k2, const double & ctheta) const
            {
                static constexpr double mB  = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mP2 = Process_::mP2, mP22 = mP2 * mP2;

                const double lambda = eos::lambda(q2, k2, mB2);
                const double E2     = (mB2 + k2 - q2 - ctheta * std::sqrt(lambda)) / (4.0 * mB);
                const double qhat2  = mB2 + mP22 - 2.0 * mB * E2;

                const double z  = this->_z(q2);
                const double zh = this->_zhat(qhat2);

                const double a = _a_Ftime_0_0 + _a_Ftime_1_0 * z + _a_Ftime_0_1 * zh + _a_Ftime_1_1 * z * zh + _a_Ftime_1_2 * z * zh * zh + _a_Ftime_0_2 * zh * zh + _a_Ftime_0_3 * zh * zh * zh;
                const double b = _b_Ftime_0_0 + _b_Ftime_1_0 * z + _b_Ftime_0_1 * zh + _b_Ftime_1_1 * z * zh + _b_Ftime_1_2 * z * zh * zh + _b_Ftime_0_2 * zh * zh + _b_Ftime_0_3 * zh * zh * zh;
                const double c = _c_Ftime_0_0 + _c_Ftime_1_0 * z + _c_Ftime_0_1 * zh + _c_Ftime_1_1 * z * zh + _c_Ftime_1_2 * z * zh * zh + _c_Ftime_0_2 * zh * zh + _c_Ftime_0_3 * zh * zh * zh;

                const double blaschke = this->_blaschke(z, zh);

                complex<double> result{ 0.0, blaschke * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB * mB2 / std::sqrt(q2) / k2 };

                return result;
            }

            virtual double f_time_im_res_qhat2(const double & q2, const double & k2) const
            {
                static constexpr double mB    = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                const double z  = this->_z(q2);
                const double zh = this->_z(mBst2);

                const double a = _a_Ftime_0_0 + _a_Ftime_1_0 * z + _a_Ftime_0_1 * zh + _a_Ftime_1_1 * z * zh + _a_Ftime_1_2 * z * zh * zh + _a_Ftime_0_2 * zh * zh + _a_Ftime_0_3 * zh * zh * zh;
                const double b = _b_Ftime_0_0 + _b_Ftime_1_0 * z + _b_Ftime_0_1 * zh + _b_Ftime_1_1 * z * zh + _b_Ftime_1_2 * z * zh * zh + _b_Ftime_0_2 * zh * zh + _b_Ftime_0_3 * zh * zh * zh;
                const double c = _c_Ftime_0_0 + _c_Ftime_1_0 * z + _c_Ftime_0_1 * zh + _c_Ftime_1_1 * z * zh + _c_Ftime_1_2 * z * zh * zh + _c_Ftime_0_2 * zh * zh + _c_Ftime_0_3 * zh * zh * zh;

                const double blaschke_res_qhat2 = this->_blaschke_res_qhat2(z);

                double result = blaschke_res_qhat2 * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB * mB2 / std::sqrt(q2) / k2;

                return result;
            }
    };
}

#endif
