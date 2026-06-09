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
#include <eos/maths/power-of.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>

#include <array>

namespace eos
{
    /* Form Factors according to [BSZ:2015A] */
    template <typename Process_, typename Transition_> class BSZ2015FormFactors;

    template <typename Process_, typename Transition_> class BSZ2015FormFactorTraits;


    // P -> V
    template <typename Process_>
    class BSZ2015FormFactorTraits<Process_, PToV> :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and should match the
            // the ones used for the extraction of the coefficients of the z-expension
            UsedParameter m_B, m_V;
            UsedParameter m_R_0m, m_R_1m, m_R_1p;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1p_names;

            BSZ2015FormFactorTraits(const Parameters & p) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@BSZ2015"], *this)),
                m_V(UsedParameter(p[std::string(Process_::name_V) + "@BSZ2015"], *this)),
                m_R_0m(UsedParameter(p[resonance_0m_names.at(Process_::partonic_transition)], *this)),
                m_R_1m(UsedParameter(p[resonance_1m_names.at(Process_::partonic_transition)], *this)),
                m_R_1p(UsedParameter(p[resonance_1p_names.at(Process_::partonic_transition)], *this))
            {
            }

            double tp() const
            {
                return power_of<2>(m_B + m_V);
            }

            double tm() const
            {
                return power_of<2>(m_B - m_V);
            }

            double t0() const
            {
                return tp() * (1.0 - std::sqrt(1.0 - tm() / tp()));
            }

            complex<double> calc_z(const complex<double> & s) const
            {
                return (std::sqrt(tp() - s) - std::sqrt(tp() - t0())) / (std::sqrt(tp() - s) + std::sqrt(tp() - t0()));
            }

            double calc_z(const double & s) const
            {
                return real(calc_z(complex<double>(s, 0.0)));
            }
    };

    template <typename Process_> class BSZ2015FormFactors<Process_, PToV> :
        public FormFactors<PToV>
    {
        private:
            // fit parametrization for P -> V according to [BSZ:2015A]
            std::array<UsedParameter, 3> _a_A0, _a_A1, _a_V, _a_T1, _a_T23;
            // use constraint (B.6) in [BSZ:2015A] to remove A_12(0)
            std::array<UsedParameter, 2> _a_A12, _a_T2;

            const BSZ2015FormFactorTraits<Process_, PToV> _traits;

            const UsedParameter & _mB, _mV;

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
    };

    extern template class BSZ2015FormFactors<BToDstar, PToV>;
    extern template class BSZ2015FormFactors<BToKstar, PToV>;
    extern template class BSZ2015FormFactors<BToOmega, PToV>;
    extern template class BSZ2015FormFactors<BToRho, PToV>;
    extern template class BSZ2015FormFactors<BcToJpsi, PToV>;
    extern template class BSZ2015FormFactors<BsToDsstar, PToV>;
    extern template class BSZ2015FormFactors<BsToKstar, PToV>;
    extern template class BSZ2015FormFactors<BsToPhi, PToV>;


    // P -> P
    template <typename Process_>
    class BSZ2015FormFactorTraits<Process_, PToP> :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and should match the
            // the ones used for the extraction of the coefficients of the z-expension
            UsedParameter m_B, m_P;
            UsedParameter m_R_0p, m_R_1m;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0p_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1m_names;

            BSZ2015FormFactorTraits(const Parameters & p) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@BSZ2015"], *this)),
                m_P(UsedParameter(p[std::string(Process_::name_P) + "@BSZ2015"], *this)),
                m_R_0p(UsedParameter(p[resonance_0p_names.at(Process_::partonic_transition)], *this)),
                m_R_1m(UsedParameter(p[resonance_1m_names.at(Process_::partonic_transition)], *this))
            {
            }

            double tp() const
            {
                return power_of<2>(m_B + m_P);
            }

            double tm() const
            {
                return power_of<2>(m_B - m_P);
            }

            double t0() const
            {
                return tp() * (1.0 - std::sqrt(1.0 - tm() / tp()));
            }

            complex<double> calc_z(const complex<double> & s) const
            {
                return (std::sqrt(tp() - s) - std::sqrt(tp() - t0())) / (std::sqrt(tp() - s) + std::sqrt(tp() - t0()));
            }

            double calc_z(const double & s) const
            {
                return real(calc_z(complex<double>(s, 0.0)));
            }
    };

    template <typename Process_> class BSZ2015FormFactors<Process_, PToP> :
        public FormFactors<PToP>
    {
        private:
            // fit parametrization for P -> P inspired by [BSZ:2015A]
            std::array<UsedParameter, 3> _a_fp, _a_ft;
            // use equation of motion to remove f_0(0) as a free parameter
            std::array<UsedParameter, 2> _a_fz;

            const BSZ2015FormFactorTraits<Process_, PToP> _traits;

            const UsedParameter & _mB, _mP;

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

    extern template class BSZ2015FormFactors<BToD, PToP>;
    extern template class BSZ2015FormFactors<BToEta, PToP>;
    extern template class BSZ2015FormFactors<BToEtaPrime, PToP>;
    extern template class BSZ2015FormFactors<BToK, PToP>;
    extern template class BSZ2015FormFactors<BToPi, PToP>;
    extern template class BSZ2015FormFactors<BsToDs, PToP>;
    extern template class BSZ2015FormFactors<BsToEta, PToP>;
    extern template class BSZ2015FormFactors<BsToEtaPrime, PToP>;
    extern template class BSZ2015FormFactors<BsToK, PToP>;
    extern template class BSZ2015FormFactors<DToK, PToP>;
    extern template class BSZ2015FormFactors<DToEta, PToP>;
    extern template class BSZ2015FormFactors<DToEtaPrime, PToP>;
    extern template class BSZ2015FormFactors<DToPi, PToP>;
    extern template class BSZ2015FormFactors<DsToK, PToP>;
    extern template class BSZ2015FormFactors<DsToEta, PToP>;
    extern template class BSZ2015FormFactors<DsToEtaPrime, PToP>;
}

#endif
