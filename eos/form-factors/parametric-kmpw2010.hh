/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013-2016, 2018 Danny van Dyk
 * Copyright (c) 2015 Christoph Bobeth
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KMPW2010_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KMPW2010_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>

namespace eos
{
    /* Form Factors according to [KMPW:2010A] */
    template <typename Transition_> class KMPW2010FormFactors;

    template <> class KMPW2010FormFactors<PToV> :
        public FormFactors<PToV>
    {
        private:
            // fit parametrisation for P -> V according to [KMPW:2010]
            UsedParameter
               _f0_V, _b1_V,
               _f0_A0, _b1_A0, _f0_A1, _b1_A1, _f0_A2, _b1_A2,
               _f0_T1, _b1_T1, _f0_T2, _b1_T2, _f0_T3, _b1_T3;

            static const double _tau_p, _tau_m, _tau_0;
            static const double _m_B, _m_Kstar, _m_Bs2_0m, _m_Bs2_1m, _m_Bs2_1p;

            static double _calc_z(const double & s);

            static double ff_KMPW(const double & s, const double & f0, const double & b1, const double & m2);

        public:
            KMPW2010FormFactors(const Parameters & p, const Options &);

            static FormFactors<PToV> * make(const Parameters & parameters, const Options & options);

            virtual double v(const double & s) const;

            virtual double a_0(const double & s) const;

            virtual double a_1(const double & s) const;

            virtual double a_2(const double & s) const;

            virtual double a_12(const double & s) const;

            virtual double t_1(const double & s) const;

            virtual double t_2(const double & s) const;

            virtual double t_3(const double & s) const;

            virtual double t_23(const double & s) const;

            virtual double f_perp(const double &) const;

            virtual double f_para(const double &) const;

            virtual double f_long(const double &) const;

            virtual double f_perp_T(const double &) const;

            virtual double f_para_T(const double &) const;

            virtual double f_long_T(const double &) const;
    };

    extern template class KMPW2010FormFactors<PToV>;

    template <> class KMPW2010FormFactors<PToP> :
        public FormFactors<PToP>
    {
        private:
            // fit parametrisation for P -> P according to [KMPW:2010]
            UsedParameter _b1_p, _b1_0, _b1_t;
            UsedParameter _f0_p, _f0_t;
            static const double _tau_p, _tau_m, _tau_0;
            static const double _m_B, _m_K, _m_Bs2;

            static double _calc_z(const double & s);

        public:
            KMPW2010FormFactors(const Parameters & p, const Options &);

            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options);

            virtual double f_p(const double & s) const;

            virtual double f_0(const double & s) const;

            virtual double f_t(const double & s) const;

            virtual double f_plus_T(const double &) const;
    };

    extern template class KMPW2010FormFactors<PToP>;
}

#endif
