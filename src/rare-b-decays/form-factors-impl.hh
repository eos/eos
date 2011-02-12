/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
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
                _v_factor(p["formfactors::v_uncertainty"]),
                _a0_factor(p["formfactors::a0_uncertainty"]),
                _a1_factor(p["formfactors::a1_uncertainty"]),
                _a2_factor(p["formfactors::a2_uncertainty"])
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
            Parameter _f_p_factor, _f_0_factor, _f_t_factor;

            // fit parametrisation for P -> P according to [BZ2004v2]
            static const double _r1_p, _r2_p, _r1_t, _r2_t, _r2_0;
            static const double _mfit2, _m12;

        public:
            BZ2004FormFactors(const Parameters & p, const Options &) :
                _f_p_factor(p["formfactors::fp_uncertainty"]),
                _f_0_factor(p["formfactors::f0_uncertainty"]),
                _f_t_factor(p["formfactors::ft_uncertainty"])
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

            static double _calc_z(const double & s)
            {
                return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
            }

        public:
            KMPW2010FormFactors(const Parameters & p, const Options &) :
                _f_p_factor(p["formfactors::fp_uncertainty"]),
                _f_0_factor(p["formfactors::f0_uncertainty"]),
                _f_t_factor(p["formfactors::ft_uncertainty"])
            {
            }

            static FormFactors<PToP> * make(const Parameters & parameters, unsigned)
            {
                return new KMPW2010FormFactors(parameters, Options());
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
}

#endif
