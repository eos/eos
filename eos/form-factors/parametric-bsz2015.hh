/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015 Frederik Beaujean
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BSZ2015_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BSZ2015_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>

#include <array>

namespace eos
{
    /* Form Factors according to [BSZ2015] */
    template <typename Process_, typename Transition_> class BSZ2015FormFactors;

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

            static double _calc_tau_0(const double & m_B, const double & m_V);

            complex<double> _calc_z(const complex<double> & s) const;

            double _calc_z(const double & s) const;

            template <typename Parameter_>
            complex<double> _calc_ff(const complex<double> & s, const double & m2_R, const std::array<Parameter_, 3> & a) const;

            static std::string _par_name(const std::string & ff_name);

        public:
            BSZ2015FormFactors(const Parameters & p, const Options &);

            ~BSZ2015FormFactors();

            static FormFactors<PToV> * make(const Parameters & parameters, const Options & options);

            virtual complex<double> v(const complex<double> & s) const;

            virtual complex<double> a_0(const complex<double> & s) const;

            virtual complex<double> a_1(const complex<double> & s) const;

            virtual complex<double> a_12(const complex<double> & s) const;

            virtual complex<double> a_2(const complex<double> & s) const;

            virtual complex<double> t_1(const complex<double> & s) const;

            virtual complex<double> t_2(const complex<double> & s) const;

            virtual complex<double> t_23(const complex<double> & s) const;

            virtual double v(const double & s) const;

            virtual double a_0(const double & s) const;

            virtual double a_1(const double & s) const;

            virtual double a_12(const double & s) const;

            virtual double a_2(const double & s) const;

            virtual double t_1(const double & s) const;

            virtual double t_2(const double & s) const;

            virtual double t_23(const double & s) const;

            virtual double t_3(const double & s) const;

            virtual double f_perp(const double & s) const;

            virtual double f_para(const double & s) const;

            virtual double f_long(const double & s) const;

            virtual double f_perp_T(const double & s) const;

            virtual double f_para_T(const double & s) const;

            virtual double f_long_T(const double & s) const;

            virtual double f_long_T_Normalized(const double & s) const;
    };

    extern template class BSZ2015FormFactors<BToRho, PToV>;
    extern template class BSZ2015FormFactors<BToOmega, PToV>;
    extern template class BSZ2015FormFactors<BToKstar, PToV>;
    extern template class BSZ2015FormFactors<BToDstar, PToV>;
    extern template class BSZ2015FormFactors<BsToKstar, PToV>;
    extern template class BSZ2015FormFactors<BsToPhi, PToV>;
    extern template class BSZ2015FormFactors<BsToDsstar, PToV>;

    template <typename Process_> class BSZ2015FormFactors<Process_, PToP> :
        public FormFactors<PToP>
    {
        private:
            // fit parametrization for P -> P inspired by [BSZ2015]
            std::array<UsedParameter, 3> _a_fp, _a_ft;
            // use equation of motion to remove f_0(0) as a free parameter
            std::array<UsedParameter, 2> _a_fz;

            const double _mB, _mB2, _mP, _mP2;
            const double _tau_p, _tau_0;
            const double _z_0;

            static double _calc_tau_0(const double & m_B, const double & m_P);

            complex<double> _calc_z(const complex<double> & s) const;

            double _calc_z(const double & s) const;

            template <typename Parameter_>
            complex<double> _calc_ff(const complex<double> & s, const double & m2_R, const std::array<Parameter_, 3> & a) const;

            static std::string _par_name(const std::string & ff_name);

        public:
            BSZ2015FormFactors(const Parameters & p, const Options &);

            ~BSZ2015FormFactors();

            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options);

            virtual complex<double> f_p(const complex<double> & s) const;

            virtual complex<double> f_t(const complex<double> & s) const;

            virtual complex<double> f_0(const complex<double> & s) const;

            virtual complex<double> f_plus_T(const complex<double> & s) const;

            virtual double f_p(const double & s) const;

            virtual double f_t(const double & s) const;

            virtual double f_0(const double & s) const;

            virtual double f_plus_T(const double & s) const;
    };

    extern template class BSZ2015FormFactors<BToPi, PToP>;
    extern template class BSZ2015FormFactors<BToK, PToP>;
    extern template class BSZ2015FormFactors<BToD, PToP>;
    extern template class BSZ2015FormFactors<BsToK, PToP>;
    extern template class BSZ2015FormFactors<BsToDs, PToP>;
    extern template class BSZ2015FormFactors<DToPi, PToP>;
    extern template class BSZ2015FormFactors<DsToK, PToP>;
    extern template class BSZ2015FormFactors<DToK, PToP>;
}

#endif
