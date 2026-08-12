/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_SSE_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_SSE_HH 1

#include <eos/form-factors/baryonic.hh>
#include <eos/form-factors/baryonic-processes.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>

#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <string>
#include <tuple>

namespace eos
{
    /* Form Factors according to [BSZ:2015A] */
    template <typename Process_, typename Transition_> class SSEFormFactors;

    template <typename Process_, typename Transition_> class SSEFormFactorTraits;


    // P -> V
    template <typename Process_>
    class SSEFormFactorTraits<Process_, PToV> :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and should match the
            // the ones used for the extraction of the coefficients of the z-expension
            UsedParameter m_B, m_V;
            UsedParameter m_R_0m, m_R_1m, m_R_1p;
            // Pair-production thresholds: tp_v for the vector/scalar current (V, T_1),
            // tp_a for the (pseudo)scalar/axial current (A_0, A_1, A_12, T_2, T_23).
            UsedParameter tp_v, tp_a;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1p_names;

            SSEFormFactorTraits(const Parameters & p) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@HME"], *this)),
                m_V(UsedParameter(p[std::string(Process_::name_V) + "@HME"], *this)),
                m_R_0m(UsedParameter(p[resonance_0m_names.at(Process_::partonic_transition)], *this)),
                m_R_1m(UsedParameter(p[resonance_1m_names.at(Process_::partonic_transition)], *this)),
                m_R_1p(UsedParameter(p[resonance_1p_names.at(Process_::partonic_transition)], *this)),
                tp_v(UsedParameter(p[std::string(Process_::label) + "::tp_v@SSE"], *this)),
                tp_a(UsedParameter(p[std::string(Process_::label) + "::tp_a@SSE"], *this))
            {
            }

            double tm() const
            {
                return power_of<2>(m_B - m_V);
            }

            // Optimized expansion point t_0 for the given pair-production threshold t_+.
            double t0(const double & tp) const
            {
                return tp * (1.0 - std::sqrt(1.0 - tm() / tp));
            }

            complex<double> calc_z(const complex<double> & s, const double & tp) const
            {
                const double t_0 = t0(tp);
                return (std::sqrt(tp - s) - std::sqrt(tp - t_0)) / (std::sqrt(tp - s) + std::sqrt(tp - t_0));
            }

            double calc_z(const double & s, const double & tp) const
            {
                return real(calc_z(complex<double>(s, 0.0), tp));
            }
    };

    template <typename Process_> class SSEFormFactors<Process_, PToV> :
        public FormFactors<PToV>
    {
        private:
            // Form factors with a free leading coefficient (a_0, a_1, a_2).
            std::array<UsedParameter, 3> _a_A0, _a_V, _a_T1;
            // Form factors whose leading coefficient is fixed by a kinematic relation
            // and therefore only carry (a_1, a_2):
            //   A_12(0)   = R * A_0(0)    [constraint (B.6) in [BSZ:2015A]]
            //   A_12(t_-) = R * A_1(t_-)  [required for A_2 to be finite at t_-]
            //   T_2(0)    = T_1(0)
            //   T_23(t_-) = K * T_2(t_-)  [required for T_3 to be finite at t_-]
            std::array<UsedParameter, 2> _a_A1, _a_A12, _a_T2, _a_T23;

            const SSEFormFactorTraits<Process_, PToV> _traits;

            const UsedParameter & _mB, _mV;

            template <typename Parameter_>
            complex<double> _calc_ff(const complex<double> & s, const double & m2_R, const double & tp, const std::array<Parameter_, 3> & a) const;

            static std::string _par_name(const std::string & ff_name);

            // Kinematic factors appearing in the relations above.
            double _R() const;
            double _K() const;

            // Leading coefficients removed by the kinematic relations.
            double _a_A12_0() const;
            double _a_A1_0() const;
            double _a_T2_0() const;
            double _a_T23_0() const;

        public:
            SSEFormFactors(const Parameters & p, const Options &);

            ~SSEFormFactors();

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

    extern template class SSEFormFactors<BToDstar, PToV>;
    extern template class SSEFormFactors<BToKstar, PToV>;
    extern template class SSEFormFactors<BToOmega, PToV>;
    extern template class SSEFormFactors<BToRho, PToV>;
    extern template class SSEFormFactors<BcToJpsi, PToV>;
    extern template class SSEFormFactors<BsToDsstar, PToV>;
    extern template class SSEFormFactors<BsToKstar, PToV>;
    extern template class SSEFormFactors<BsToPhi, PToV>;


    // P -> P
    template <typename Process_>
    class SSEFormFactorTraits<Process_, PToP> :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and should match the
            // the ones used for the extraction of the coefficients of the z-expension
            UsedParameter m_B, m_P;
            UsedParameter m_R_0p, m_R_1m;
            // Pair-production threshold: all P -> P form factors couple to the
            // vector/scalar current and share the lowest 2-particle threshold.
            UsedParameter tp;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0p_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1m_names;

            SSEFormFactorTraits(const Parameters & p) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@HME"], *this)),
                m_P(UsedParameter(p[std::string(Process_::name_P) + "@HME"], *this)),
                m_R_0p(UsedParameter(p[resonance_0p_names.at(Process_::partonic_transition)], *this)),
                m_R_1m(UsedParameter(p[resonance_1m_names.at(Process_::partonic_transition)], *this)),
                tp(UsedParameter(p[std::string(Process_::label) + "::tp@SSE"], *this))
            {
            }

            double tm() const
            {
                return power_of<2>(m_B - m_P);
            }

            // Optimized expansion point t_0 for the given pair-production threshold t_+.
            double t0(const double & tp) const
            {
                return tp * (1.0 - std::sqrt(1.0 - tm() / tp));
            }

            complex<double> calc_z(const complex<double> & s, const double & tp) const
            {
                const double t_0 = t0(tp);
                return (std::sqrt(tp - s) - std::sqrt(tp - t_0)) / (std::sqrt(tp - s) + std::sqrt(tp - t_0));
            }

            double calc_z(const double & s, const double & tp) const
            {
                return real(calc_z(complex<double>(s, 0.0), tp));
            }
    };

    template <typename Process_> class SSEFormFactors<Process_, PToP> :
        public FormFactors<PToP>
    {
        private:
            // fit parametrization for P -> P inspired by [BSZ:2015A]
            std::array<UsedParameter, 3> _a_fp, _a_ft;
            // use equation of motion to remove f_0(0) as a free parameter
            std::array<UsedParameter, 2> _a_fz;

            const SSEFormFactorTraits<Process_, PToP> _traits;

            const UsedParameter & _mB, _mP;

            template <typename Parameter_>
            complex<double> _calc_ff(const complex<double> & s, const double & m2_R, const double & tp, const std::array<Parameter_, 3> & a) const;

            static std::string _par_name(const std::string & ff_name);

        public:
            SSEFormFactors(const Parameters & p, const Options &);

            ~SSEFormFactors();

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

    extern template class SSEFormFactors<BToD, PToP>;
    extern template class SSEFormFactors<BToEta, PToP>;
    extern template class SSEFormFactors<BToEtaPrime, PToP>;
    extern template class SSEFormFactors<BToK, PToP>;
    extern template class SSEFormFactors<BToPi, PToP>;
    extern template class SSEFormFactors<BsToDs, PToP>;
    extern template class SSEFormFactors<BsToEta, PToP>;
    extern template class SSEFormFactors<BsToEtaPrime, PToP>;
    extern template class SSEFormFactors<BsToK, PToP>;
    extern template class SSEFormFactors<DToK, PToP>;
    extern template class SSEFormFactors<DToEta, PToP>;
    extern template class SSEFormFactors<DToEtaPrime, PToP>;
    extern template class SSEFormFactors<DToPi, PToP>;
    extern template class SSEFormFactors<DsToK, PToP>;
    extern template class SSEFormFactors<DsToEta, PToP>;
    extern template class SSEFormFactors<DsToEtaPrime, PToP>;


    // 1/2^+ -> 1/2^+
    template <typename Process_>
    class SSEFormFactorTraits<Process_, OneHalfPlusToOneHalfPlus> :
        public virtual ParameterUser
    {
        public:
            // m_1 is the mass of the heavier baryon, m_2 the mass of the lighter one.
            // The baryon masses are shared with the [SE] parametrization; the resonance
            // masses use the same @HME inputs as the mesonic SSE parametrization.
            UsedParameter m_1, m_2;
            UsedParameter m_R_0m, m_R_0p, m_R_1m, m_R_1p;
            // Pair-production thresholds: tp_v for the vector/scalar and tensor (T)
            // currents, tp_a for the (pseudo)scalar/axial and pseudo-tensor (T5) currents.
            UsedParameter tp_v, tp_a;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0p_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1p_names;

            SSEFormFactorTraits(const Parameters & p) :
                m_1(UsedParameter(p[std::string(Process_::name_1) + "@HME"], *this)),
                m_2(UsedParameter(p[std::string(Process_::name_2) + "@HME"], *this)),
                m_R_0m(UsedParameter(p[resonance_0m_names.at(Process_::partonic_transition)], *this)),
                m_R_0p(UsedParameter(p[resonance_0p_names.at(Process_::partonic_transition)], *this)),
                m_R_1m(UsedParameter(p[resonance_1m_names.at(Process_::partonic_transition)], *this)),
                m_R_1p(UsedParameter(p[resonance_1p_names.at(Process_::partonic_transition)], *this)),
                tp_v(UsedParameter(p[std::string(Process_::label) + "::tp_v@SSE"], *this)),
                tp_a(UsedParameter(p[std::string(Process_::label) + "::tp_a@SSE"], *this))
            {
            }

            double tm() const
            {
                return power_of<2>(m_1 - m_2);
            }

            // Optimized expansion point t_0 for the given pair-production threshold t_+.
            double t0(const double & tp) const
            {
                return tp * (1.0 - std::sqrt(1.0 - tm() / tp));
            }

            complex<double> calc_z(const complex<double> & s, const double & tp) const
            {
                const double t_0 = t0(tp);
                return (std::sqrt(tp - s) - std::sqrt(tp - t_0)) / (std::sqrt(tp - s) + std::sqrt(tp - t_0));
            }

            double calc_z(const double & s, const double & tp) const
            {
                return real(calc_z(complex<double>(s, 0.0), tp));
            }
    };

    template <typename Process_> class SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus> :
        public FormFactors<OneHalfPlusToOneHalfPlus>
    {
        private:
            // Form factors with a free leading coefficient (a_0, a_1, a_2).
            std::array<UsedParameter, 3> _a_long_v, _a_perp_v, _a_long_a, _a_long_t, _a_perp_t5;
            // Form factors whose leading coefficient is fixed by an equation of motion
            // and therefore only carry (a_1, a_2).
            std::array<UsedParameter, 2> _a_time_v, _a_time_a, _a_perp_a, _a_perp_t, _a_long_t5;

            const SSEFormFactorTraits<Process_, OneHalfPlusToOneHalfPlus> _traits;

            const UsedParameter & _m_1, _m_2;

            template <typename Parameter_>
            double _calc_ff(const double & s, const double & m_R, const double & tp, const std::array<Parameter_, 3> & a) const;

            static std::string _par_name(const std::string & pol, const std::string & current, unsigned idx);

            // Leading coefficients removed by the equations of motion.
            double _a_time_v_0() const;
            double _a_time_a_0() const;
            double _a_perp_a_0() const;
            double _a_perp_t_0() const;
            double _a_long_t5_0() const;

        public:
            SSEFormFactors(const Parameters & p, const Options &);

            ~SSEFormFactors();

            static FormFactors<OneHalfPlusToOneHalfPlus> * make(const Parameters & parameters, const Options & options);

            virtual double f_time_v(const double & s) const;
            virtual double f_long_v(const double & s) const;
            virtual double f_perp_v(const double & s) const;

            virtual double f_time_a(const double & s) const;
            virtual double f_long_a(const double & s) const;
            virtual double f_perp_a(const double & s) const;

            virtual double f_long_t(const double & s) const;
            virtual double f_perp_t(const double & s) const;

            virtual double f_long_t5(const double & s) const;
            virtual double f_perp_t5(const double & s) const;
    };

    extern template class SSEFormFactors<LambdaBToLambda,  OneHalfPlusToOneHalfPlus>;
    extern template class SSEFormFactors<LambdaCToLambda,  OneHalfPlusToOneHalfPlus>;
    extern template class SSEFormFactors<LambdaCToNeutron, OneHalfPlusToOneHalfPlus>;
    extern template class SSEFormFactors<LambdaCToProton,  OneHalfPlusToOneHalfPlus>;
}

#endif
