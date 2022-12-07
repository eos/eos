/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014, 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2017 Elena Graverini
 * Copyright (c) 2017, 2018 Marzia Bordone
 * Copyright (c) 2018 Ahmet Kokulu
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BBGORVD2018_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BBGORVD2018_HH 1

#include <eos/form-factors/baryonic.hh>
#include <eos/form-factors/baryonic-processes.hh>
#include <eos/form-factors/hqet-b-to-c.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/kinematic.hh>
#include <eos/models/model.hh>
#include <eos/utils/options.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    template <typename Transition_, typename Process_> class HQETFormFactors;

    /*
     * J=1/2^+ -> J=1/2^- transitions
     */

    template <typename Process_>
    class HQETFormFactors<OneHalfPlusToOneHalfMinus, Process_> :
        public FormFactors<OneHalfPlusToOneHalfMinus>
    {
        private:
            HQETBToC _b_to_c;

            UsedParameter _zeta_max, _rho, _delta_3b, _rho_3b;

            static constexpr double mLb = Process_::m1;
            static constexpr double mLcs = Process_::m2;
            static constexpr double mLb2 = mLb * mLb;
            static constexpr double mLcs2 = mLcs * mLcs;

            static constexpr double m_b_pole = 4.8;
            static constexpr double m_c_pole = 1.4;

            static constexpr double lambdabar = mLb - m_b_pole;
            static constexpr double lambdabarprime = mLcs - m_c_pole;

            static constexpr double s_max = (mLb - mLcs) * (mLb - mLcs);

            // auxiliary kinematics functions
            static constexpr double _s_plus(const double & s)
            {
                return power_of<2>(mLb + mLcs) - s;
            }
            static constexpr double _s_minus(const double & s)
            {
                return power_of<2>(mLb - mLcs) - s;
            }

            // parametrization of the Isgur-Wise functions
            double _z(const double & s) const
            {
                return _zeta_max * (1.0 + _rho * (s / s_max - 1.0));
            }

            double _z3b(const double & s) const
            {
                return _zeta_max * (_delta_3b + _rho_3b * (s / s_max - 1.0));
            }

            inline double omega(const double & s) const
            {
                return (mLb2 + mLcs2 - s) / (2.0 * mLb * mLcs);
            }

            inline double omegabar(const double & s) const
            {
                return omega(s) * (1.0 + lambdabar / m_b_pole + lambdabarprime / m_c_pole)
                    - (lambdabar / m_c_pole + lambdabarprime / m_b_pole);
            }

        public:
            HQETFormFactors(const Parameters & p);
            ~HQETFormFactors();

            static FormFactors<OneHalfPlusToOneHalfMinus> * make(const Parameters &, const Options &);

            // vector current
            virtual double f_time_v(const double & s) const;
            virtual double f_long_v(const double & s) const;
            virtual double f_perp_v(const double & s) const;

            // axial vector current
            virtual double f_time_a(const double & s) const;
            virtual double f_long_a(const double & s) const;
            virtual double f_perp_a(const double & s) const;

            Diagnostics diagnostics() const;
    };

    extern template class HQETFormFactors<OneHalfPlusToOneHalfMinus, LambdaBToLambdaC2595>;

    /*
     * J=1/2^+ -> J=3/2^- transitions
     */

    template <typename Process_>
    class HQETFormFactors<OneHalfPlusToThreeHalfMinus, Process_> :
        public FormFactors<OneHalfPlusToThreeHalfMinus>
    {
        private:
            HQETBToC _b_to_c;

            UsedParameter _zeta_max, _rho, _delta_3b, _rho_3b;

            static constexpr double mLb = Process_::m1;
            static constexpr double mLcs = Process_::m2;
            static constexpr double mLb2 = mLb * mLb;
            static constexpr double mLcs2 = mLcs * mLcs;

            static constexpr double m_b_pole = 4.8;
            static constexpr double m_c_pole = 1.4;

            static constexpr double lambdabar = mLb - m_b_pole;
            static constexpr double lambdabarprime = mLcs - m_c_pole;

            static constexpr double s_max = (mLb - mLcs) * (mLb - mLcs);

            // auxiliary kinematics functions
            static constexpr double _s_plus(const double & s)
            {
                return power_of<2>((mLb + mLcs)) - s;
            }
            static constexpr double _s_minus(const double & s)
            {
                return power_of<2>((mLb - mLcs)) - s;
            }

            // parametrization of the Isgur-Wise functions
            double _z(const double & s) const
            {
                return _zeta_max * (1.0 + _rho * (s / s_max - 1.0));
            }

            double _z3b(const double & s) const
            {
                return _zeta_max * (_delta_3b + _rho_3b * (s / s_max - 1.0));
            }

            inline double omega(const double & s) const
            {
                return (mLb2 + mLcs2 - s) / (2.0 * mLb * mLcs);
            }

            inline double omegabar(const double & s) const
            {
                return omega(s) * (1.0 + lambdabar / m_b_pole + lambdabarprime / m_c_pole)
                    - (lambdabar / m_c_pole + lambdabarprime / m_b_pole);
            }

        public:
            HQETFormFactors(const Parameters & p);
            ~HQETFormFactors();

            static FormFactors<OneHalfPlusToThreeHalfMinus> * make(const Parameters &, const Options &);

            // vector current
            virtual double f_time12_v(const double & s) const;
            virtual double f_long12_v(const double & s) const;
            virtual double f_perp12_v(const double & s) const;
            virtual double f_perp32_v(const double & s) const;

            // axial vector current
            virtual double f_time12_a(const double & s) const;
            virtual double f_long12_a(const double & s) const;
            virtual double f_perp12_a(const double & s) const;
            virtual double f_perp32_a(const double & s) const;

            // tensor current
            virtual double f_long12_t(const double &) const { throw InternalError("HQETFormFactors::f_long12_t(): not implemented"); }
            virtual double f_perp12_t(const double &) const { throw InternalError("HQETFormFactors::f_perp12_t(): not implemented"); }
            virtual double f_perp32_t(const double &) const { throw InternalError("HQETFormFactors::f_perp32_t(): not implemented"); }
            virtual double f_long12_t5(const double &) const { throw InternalError("HQETFormFactors::f_long12_t5(): not implemented"); }
            virtual double f_perp12_t5(const double &) const { throw InternalError("HQETFormFactors::f_perp12_t5(): not implemented"); }
            virtual double f_perp32_t5(const double &) const { throw InternalError("HQETFormFactors::f_perp32_t5(): not implemented"); }

            virtual Diagnostics diagnostics() const;
    };

    extern template class HQETFormFactors<OneHalfPlusToThreeHalfMinus, LambdaBToLambdaC2625>;
}

#endif
