/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <src/rare-b-decays/form-factors.hh>
#include <src/utils/options.hh>
#include <src/utils/power_of.hh>

namespace eos
{
    /* Form Factors according to [BZ2004] */
    template <typename Process_, typename Transition_> class BZ2004FormFactors;

    /* Form Factors according to [KMPW2010] */
    template <typename Transition_> class KMPW2010FormFactors;


    /* P -> V Processes */

    struct BToKstar { };
    struct BsToPhi { };

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
                _v_factor(p["formfactors::v_uncertainty"], *this),
                _a0_factor(p["formfactors::a0_uncertainty"], *this),
                _a1_factor(p["formfactors::a1_uncertainty"], *this),
                _a2_factor(p["formfactors::a2_uncertainty"], *this)
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
    };

    /* P -> P Processes */

    struct BToK { };

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
                _f_p_factor(p["formfactors::fp_uncertainty"], *this),
                _f_0_factor(p["formfactors::f0_uncertainty"], *this),
                _f_t_factor(p["formfactors::ft_uncertainty"], *this)
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
                _f_p_factor(p["formfactors::fp_uncertainty"], *this),
                _f_0_factor(p["formfactors::f0_uncertainty"], *this),
                _f_t_factor(p["formfactors::ft_uncertainty"], *this),
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
            UsedParameter _f0_p, _f0_0, _f0_t;
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
                _f0_0(p["B->K::F^0(0)@KMPW2010"], *this),
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

                return _f0_0() * (1 + _b1_0() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

            virtual double f_t(const double & s) const
            {
                // cf. [KMPW2010], Eq. (8.8), p. 30
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                return _f0_t() / (1 - s / _m_Bs2) * (1 + _b1_t() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }
    };
}

#endif
