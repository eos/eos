/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Méril Reboud
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BFW2010_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BFW2010_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/power-of.hh>
#include <eos/maths/szego-polynomial.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>
#include <eos/utils/reference-name.hh>

#include <array>

namespace eos
{
    /* Form Factors according to [BFW:2010A] */
    template <typename Process_, typename Transition_> class BFW2010FormFactors;

    template <typename Process_, typename Transition_> class BFW2010FormFactorTraits;


    // P -> V
    template <typename Process_>
    class BFW2010FormFactorTraits<Process_, PToV> :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and should match the
            // the ones used for the extraction of the coefficients of the z-expension
            UsedParameter m_B, m_V;
            UsedParameter m_R_0m, m_R_1m, m_R_1p;
            UsedParameter tp_a, tp_v, t0;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1p_names;

            BFW2010FormFactorTraits(const Parameters & p) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@BSZ2015"], *this)),
                m_V(UsedParameter(p[std::string(Process_::name_V) + "@BSZ2015"], *this)),
                m_R_0m(UsedParameter(p[resonance_0m_names.at(Process_::partonic_transition)], *this)),
                m_R_1m(UsedParameter(p[resonance_1m_names.at(Process_::partonic_transition)], *this)),
                m_R_1p(UsedParameter(p[resonance_1p_names.at(Process_::partonic_transition)], *this)),
                tp_a(UsedParameter(p[std::string(Process_::label) + "::tp_a@BFW2010"], *this)),
                tp_v(UsedParameter(p[std::string(Process_::label) + "::tp_v@BFW2010"], *this)),
                t0(UsedParameter(p[std::string(Process_::label) + "::t0@BFW2010"], *this))
            {
            }

            double tm() const
            {
                return power_of<2>(m_B - m_V);
            }

            complex<double> calc_z(const complex<double> & s, const complex<double> & sp, const complex<double> & s0) const
            {
                return (std::sqrt(sp - s) - std::sqrt(sp - s0)) / (std::sqrt(sp - s) + std::sqrt(sp - s0));
            }

            double calc_z(const double & s, const double & sp, const double & s0) const
            {
                if (s > sp)
                    throw InternalError("The real conformal mapping is used above threshold: " + stringify(s) + " > " + stringify(sp));

                return real(calc_z(complex<double>(s, 0.0), complex<double>(sp, 0.0), complex<double>(s0, 0.0)));
            }

            std::array<double, 6> orthonormal_polynomials_v(const double & z) const
            {
                const double measure = 2 * std::arg(calc_z(complex<double>(power_of<2>(m_B + m_V)), complex<double>(tp_v), complex<double>(t0)));
                const SzegoPolynomial<5> polynomials_set(SzegoPolynomial<5>::FlatMeasure(measure));

                return polynomials_set(z);
            }

            std::array<double, 6> orthonormal_polynomials_a(const double & z) const
            {
                const double measure = 2 * std::arg(calc_z(complex<double>(power_of<2>(m_B + m_V)), complex<double>(tp_a), complex<double>(t0)));
                const SzegoPolynomial<5> polynomials_set(SzegoPolynomial<5>::FlatMeasure(measure));

                return polynomials_set(z);
            }
    };

    template <typename Process_> class BFW2010FormFactors<Process_, PToV> :
        public FormFactors<PToV>
    {
        private:
            // fit parametrization for P -> V according to [BFW:2010A]
            std::array<UsedParameter, 5> _a_A0, _a_V, _a_T1;
            // use endpoint relations (see eq. (3.2) in [HLMW:2015A]) to remove parameters
            std::array<UsedParameter, 4> _a_A12, _a_T2, _a_A1, _a_T23;

            const BFW2010FormFactorTraits<Process_, PToV> _traits;

            const UsedParameter & _mB, _mV;

            QualifiedName _par_name(const std::string & ff_name, unsigned idx) const;
            double _phi(const double & t, const double & t_p, const double & chi,
                        const int & A, const unsigned B, const unsigned C, const unsigned k,
                        const unsigned p, const unsigned n, const unsigned m) const;

            inline double _phi_v(const double & q2) const;
            inline double _phi_a_0(const double & q2) const;
            inline double _phi_a_1(const double & q2) const;
            inline double _phi_a_12(const double & q2) const;
            inline double _phi_t_1(const double & q2) const;
            inline double _phi_t_2(const double & q2) const;
            inline double _phi_t_23(const double & q2) const;

            // End-point relations
            double _a_A12_0() const;
            double _a_T2_0() const;
            double _a_A1_0() const;
            double _a_T23_0() const;

        public:
            BFW2010FormFactors(const Parameters & p, const Options &);

            ~BFW2010FormFactors();

            static FormFactors<PToV> * make(const Parameters & parameters, const Options & options);

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

            // Saturations of the dispersive bounds
            // J = 0
            double saturation_0p_v() const;
            double saturation_0m_a() const;
            // J = 1
            double saturation_1m_v_0() const;
            double saturation_1m_v_perp() const;
            double saturation_1m_v_para() const;
            double saturation_1m_v() const;
            double saturation_1p_a_0() const;
            double saturation_1p_a_perp() const;
            double saturation_1p_a_para() const;
            double saturation_1p_a() const;
            double saturation_1m_t_0() const;
            double saturation_1m_t_perp() const;
            double saturation_1m_t_para() const;
            double saturation_1m_t() const;
            double saturation_1p_t5_0() const;
            double saturation_1p_t5_perp() const;
            double saturation_1p_t5_para() const;
            double saturation_1p_t5() const;

            Diagnostics diagnostics() const;


            // Auxilliary functions: series and derivative of the series
            double v_series(const double & s) const;
            double a_0_series(const double & s) const;
            double a_1_series(const double & s) const;
            double a_12_series(const double & s) const;
            double t_1_series(const double & s) const;
            double t_2_series(const double & s) const;
            double t_23_series(const double & s) const;

            double v_series_prime(const double & s) const;
            double a_0_series_prime(const double & s) const;
            double a_1_series_prime(const double & s) const;
            double a_12_series_prime(const double & s) const;
            double t_1_series_prime(const double & s) const;
            double t_2_series_prime(const double & s) const;
            double t_23_series_prime(const double & s) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> options;
    };

    extern template class BFW2010FormFactors<BToKstar, PToV>;
    extern template class BFW2010FormFactors<BsToPhi, PToV>;


    // P -> P
    template <typename Process_>
    class BFW2010FormFactorTraits<Process_, PToP> :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and should match the
            // the ones used for the extraction of the coefficients of the z-expension
            UsedParameter m_B, m_P;
            UsedParameter m_R_0p, m_R_1m;
            UsedParameter tp, t0;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0p_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1m_names;

            BFW2010FormFactorTraits(const Parameters & p) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@BSZ2015"], *this)),
                m_P(UsedParameter(p[std::string(Process_::name_P) + "@BSZ2015"], *this)),
                m_R_0p(UsedParameter(p[resonance_0p_names.at(Process_::partonic_transition)], *this)),
                m_R_1m(UsedParameter(p[resonance_1m_names.at(Process_::partonic_transition)], *this)),
                tp(UsedParameter(p[std::string(Process_::label) + "::tp@BFW2010"], *this)),
                t0(UsedParameter(p[std::string(Process_::label) + "::t0@BFW2010"], *this))
            {
            }

            double tm() const
            {
                return power_of<2>(m_B - m_P);
            }

            complex<double> calc_z(const complex<double> & s, const complex<double> & sp, const complex<double> & s0) const
            {
                return (std::sqrt(sp - s) - std::sqrt(sp - s0)) / (std::sqrt(sp - s) + std::sqrt(sp - s0));
            }

            double calc_z(const double & s, const double & sp, const double & s0) const
            {
                if (s > sp)
                    throw InternalError("The real conformal mapping is used above threshold: " + stringify(s) + " > " + stringify(sp));

                return real(calc_z(complex<double>(s, 0.0), complex<double>(sp, 0.0), complex<double>(s0, 0.0)));
            }

            std::array<double, 6> orthonormal_polynomials(const double & z) const
            {
                // Since BK is the lowest b\bar{s} pair production threshold, those are just monomials
                const SzegoPolynomial<5> polynomials_set(SzegoPolynomial<5>::FlatMeasure(2 * M_PI));

                return polynomials_set(z);
            }

            std::array<complex<double>, 6> orthonormal_polynomials(const complex<double> & z) const
            {
                const double measure = 2 * std::arg(calc_z(complex<double>(power_of<2>(m_B + m_P)), complex<double>(tp), complex<double>(t0)));
                const SzegoPolynomial<5> polynomials_set(SzegoPolynomial<5>::FlatMeasure(measure));

                return polynomials_set(z);
            }

            std::array<complex<double>, 6> orthonormal_polynomials_derivatives(const complex<double> & z) const
            {
                const double measure = 2 * std::arg(calc_z(complex<double>(power_of<2>(m_B + m_P)), complex<double>(tp), complex<double>(t0)));
                const SzegoPolynomial<5> polynomials_set(SzegoPolynomial<5>::FlatMeasure(measure));

                return polynomials_set.derivatives(z);
            }
    };

    template <typename Process_> class BFW2010FormFactors<Process_, PToP> :
        public FormFactors<PToP>
    {
        private:
            // fit parametrization for P -> P inspired by [BFW:2010A]
            std::array<UsedParameter, 5> _a_fp, _a_ft;
            // use equation of motion to remove f_0(0) as a free parameter
            std::array<UsedParameter, 4> _a_f0;

            const BFW2010FormFactorTraits<Process_, PToP> _traits;

            const UsedParameter & _mB, _mP;

            QualifiedName _par_name(const std::string & ff_name, unsigned idx) const;

            double _phi(const double & s, const double & t_p, const double & chi,
                        const int & A, const unsigned B, const unsigned C, const unsigned k,
                        const unsigned p, const unsigned n, const unsigned m) const;

            inline double _phi_f_p(const double & q2) const;
            inline double _phi_f_0(const double & q2) const;
            inline double _phi_f_t(const double & q2) const;

            // End-point relations
            double _a_f0_0() const;

        public:
            BFW2010FormFactors(const Parameters & p, const Options &);

            ~BFW2010FormFactors();

            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options);

            virtual double f_p(const double & s) const;
            virtual double f_0(const double & s) const;
            virtual double f_t(const double & s) const;

            virtual double f_plus_T(const double & s) const;

            // Saturations of the dispersive bounds
            // J = 0
            double saturation_0p_v() const;
            double saturation_0m_a() const;
            // J = 1
            double saturation_1m_v_0() const;
            double saturation_1m_v_perp() const;
            double saturation_1m_v_para() const;
            double saturation_1m_v() const;
            double saturation_1p_a_0() const;
            double saturation_1p_a_perp() const;
            double saturation_1p_a_para() const;
            double saturation_1p_a() const;
            double saturation_1m_t_0() const;
            double saturation_1m_t_perp() const;
            double saturation_1m_t_para() const;
            double saturation_1m_t() const;
            double saturation_1p_t5_0() const;
            double saturation_1p_t5_perp() const;
            double saturation_1p_t5_para() const;
            double saturation_1p_t5() const;

            Diagnostics diagnostics() const;

            // Auxilliary functions: series and derivative of the series
            double f_p_series(const double & s) const;
            double f_0_series(const double & s) const;
            double f_t_series(const double & s) const;

            double f_p_series_prime(const double & s) const;
            double f_0_series_prime(const double & s) const;
            double f_t_series_prime(const double & s) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> options;
    };

    extern template class BFW2010FormFactors<BToK,  PToP>;
    extern template class BFW2010FormFactors<BsToK, PToP>;
}

#endif
