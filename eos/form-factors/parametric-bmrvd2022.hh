/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

/*
 * Copyright (c) 2021-2022 Danny van Dyk
 * Copyright (c) 2021-2022 Muslem Rahimi
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BMRvD2022_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BMRvD2022_HH 1

#include <eos/form-factors/form-factors-fwd.hh>
#include <eos/form-factors/baryonic.hh>
#include <eos/form-factors/baryonic-processes.hh>
#include <eos/models/model.hh>
#include <eos/utils/reference-name.hh>

#include <set>
#include <vector>

namespace eos
{
    template <typename Process_>
    class BMRvD2022FormFactorTraits :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and should match the
            // the ones used for the extraction of the coefficients of the z-expension
            UsedParameter m_1, m_2; // m_1 is the mass of the heavier particle, m_2 the mass of the lighter particle
            UsedParameter m_R_0m, m_R_0p, m_R_1m, m_R_1p;
            UsedParameter t0; // z(t_0) = 0, t_m is the endpoint of the semileptonic process
            UsedParameter tp_a, tp_v; // pair production thresholds

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0p_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1p_names;

            BMRvD2022FormFactorTraits(const Parameters & p) :
                m_1(UsedParameter(p[std::string(Process_::name_1) + "@BMRvD2022"], *this)),
                m_2(UsedParameter(p[std::string(Process_::name_2) + "@BMRvD2022"], *this)),
                m_R_0m(UsedParameter(p[resonance_0m_names.at(Process_::partonic_transition)], *this)),
                m_R_0p(UsedParameter(p[resonance_0p_names.at(Process_::partonic_transition)], *this)),
                m_R_1m(UsedParameter(p[resonance_1m_names.at(Process_::partonic_transition)], *this)),
                m_R_1p(UsedParameter(p[resonance_1p_names.at(Process_::partonic_transition)], *this)),
                t0(UsedParameter(p[std::string(Process_::label) + "::t0@BMRvD2022"], *this)),
                tp_a(UsedParameter(p[std::string(Process_::label) + "::tp_a@BMRvD2022"], *this)),
                tp_v(UsedParameter(p[std::string(Process_::label) + "::tp_v@BMRvD2022"], *this))
            {
            }

            double tm() const
            {
                return power_of<2>(m_1 - m_2);
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
                const double measure = 2 * std::arg(calc_z(complex<double>(power_of<2>(m_1 + m_2)), complex<double>(tp_v), complex<double>(t0)));
                const SzegoPolynomial<5> polynomials_set(SzegoPolynomial<5>::FlatMeasure(measure));

                return polynomials_set(z);
            }

            std::array<double, 6> orthonormal_polynomials_a(const double & z) const
            {
                const double measure = 2 * std::arg(calc_z(complex<double>(power_of<2>(m_1 + m_2)), complex<double>(tp_a), complex<double>(t0)));
                const SzegoPolynomial<5> polynomials_set(SzegoPolynomial<5>::FlatMeasure(measure));

                return polynomials_set(z);
            }
    };

    template <typename Process_>
    class BMRvD2022FormFactors :
        public FormFactors<OneHalfPlusToOneHalfPlus>
    {
        private:
            const std::array<UsedParameter, 4> _a_time_v;  // a_0^(time,V)  is obtained from the EoM f_t^V(q2 = 0) = f_0^V(q2 = 0)
            const std::array<UsedParameter, 5> _a_long_v;
            const std::array<UsedParameter, 5> _a_perp_v;
            const std::array<UsedParameter, 4> _a_time_a;  // a_0^(time,A)  is obtained from the EoM f_t^A(q2 = 0) = f_0^A(q2 = 0)
            const std::array<UsedParameter, 5> _a_long_a;
            const std::array<UsedParameter, 4> _a_perp_a;  // a_0^(perp,A)  is obtained from the EoM f_perp^A(q2 = t_-) = f_0^A(q2 = t_-)
            const std::array<UsedParameter, 5> _a_long_t;
            const std::array<UsedParameter, 4> _a_perp_t;  // a_0^(perp,T)  is obtained from the EoM f_perp^T(q2 = 0) = f_perp^T5(q2 = 0)

            const std::array<UsedParameter, 4> _a_long_t5; // a_0^(long,T5) is obtained from the EoM f_long^T5(q2 = t_-) = f_perp^T5(q2 = t_-)
            const std::array<UsedParameter, 5> _a_perp_t5;

            const BMRvD2022FormFactorTraits<Process_> _traits;

            const UsedParameter & _m_1, _m_2; // m_1 is the mass of the heavier particle, m_2 the mass of the lighter particle

            QualifiedName _par_name(const std::string & pol, const std::string & current, unsigned idx) const;
            double _phi(const double & s, const double & chi, const double & s_p, const double & a, const double & b, const double & c,
                    const double & d, const double & e, const double & f, const double & g) const;

            inline double _phi_time_v(const double & q2) const;
            inline double _phi_long_v(const double & q2) const;
            inline double _phi_perp_v(const double & q2) const;
            inline double _phi_time_a(const double & q2) const;
            inline double _phi_long_a(const double & q2) const;
            inline double _phi_perp_a(const double & q2) const;
            inline double _phi_long_t(const double & q2) const;
            inline double _phi_perp_t(const double & q2) const;
            inline double _phi_long_t5(const double & q2) const;
            inline double _phi_perp_t5(const double & q2) const;

            double _a_time_v_0() const;
            double _a_time_a_0() const;
            double _a_perp_a_0() const;
            double _a_perp_t_0() const;
            double _a_long_t5_0() const;

        public:
            BMRvD2022FormFactors(const Parameters & parameters, const Options & options);
            virtual ~BMRvD2022FormFactors() = default;

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

    extern template class BMRvD2022FormFactors<LambdaBToLambda>;
    extern template class BMRvD2022FormFactors<LambdaCToLambda>;
}

#endif
