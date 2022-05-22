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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_ABR2022_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_ABR2022_HH 1

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
    class ABR2022FormFactors :
        public FormFactors<OneHalfPlusToThreeHalfMinus>
    {
        private:
            const double _m_1, _m_2; // m_1 is the mass of the heavier particle, m_2 the mass of the lighter particle
            const double _t_0, _t_m, _t_p; // z(t_0) = 0, t_m is the endpoint of the semileptonic process, and t_p is the pair production threshold,

            const std::array<UsedParameter, 4> _a_time12_v;  // a_0^(time12,V) is obtained from the EoM f_time12^V(q2 = 0) \propto f_long12^V(q2 = 0)
            const std::array<UsedParameter, 4> _a_long12_v;  // a_0^(long12,V) is obtained from f_long12^V(q2 = q2max) \propto f_perp32^V(q2 = q2max)
            const std::array<UsedParameter, 4> _a_perp12_v;  // a_0^(perp12,V) is obtained from f_perp12^V(q2 = q2max) = - f_perp32^V(q2 = q2max)
            const std::array<UsedParameter, 5> _a_perp32_v;
            const std::array<UsedParameter, 4> _a_time12_a;  // a_0^(time12,A) is obtained from f_time12^A(q2 = q2max) = 0
            const std::array<UsedParameter, 4> _a_long12_a;  // a_0^(long12,A) is obtained from the EoM f_time12^A(q2 = 0) \propto f_long12^A(q2 = 0)
            const std::array<UsedParameter, 4> _a_perp12_a;  // a_0^(perp12,A) is obtained from f_perp12^A(q2 = q2max) = f_long12^A(q2 = q2max) + f_perp32^A(q2 = q2max)
            const std::array<UsedParameter, 5> _a_perp32_a;
            const std::array<UsedParameter, 4> _a_long12_t;  // a_0^(long12,T) is obtained from f_long12^T(q2 = q2max) \propto f_perp32^T(q2 = q2max)
            const std::array<UsedParameter, 4> _a_perp12_t;  // a_0^(perp12,T) is obtained from f_perp12^T(q2 = q2max) = - f_perp32^T(q2 = q2max)
            const std::array<UsedParameter, 5> _a_perp32_t;
            const std::array<UsedParameter, 4> _a_long12_t5; // a_0^(long12,T5) is obtained from f_long12^T5(q2 = q2max) = f_perp12^T5(q2 = q2max) + f_perp32^T5(q2 = q2max)
            const std::array<UsedParameter, 4> _a_perp12_t5; // a_0^(perp12,T5) is obtained from the EoM f_perp12^T5(q2 = 0) \propto f_perp12^T(q2 = 0)
            const std::array<UsedParameter, 4> _a_perp32_t5; // a_0^(perp32,T5) is obtained from the EoM f_perp32^T5(q2 = 0) \propto f_perp32^T(q2 = 0)

            QualifiedName _par_name(const std::string & pol, const std::string & current, unsigned idx) const;
            double _z(const double & t, const double & t_0) const;
            double _phi(const double & s, const double & chi, const double & A, const double & B, const double & d, const double & e,
                        const double & f, const double & g, const double & n) const;

            inline double _phi_time12_v(const double & q2) const;
            inline double _phi_long12_v(const double & q2) const;
            inline double _phi_perp12_v(const double & q2) const;
            inline double _phi_perp32_v(const double & q2) const;
            inline double _phi_time12_a(const double & q2) const;
            inline double _phi_long12_a(const double & q2) const;
            inline double _phi_perp12_a(const double & q2) const;
            inline double _phi_perp32_a(const double & q2) const;
            inline double _phi_long12_t(const double & q2) const;
            inline double _phi_perp12_t(const double & q2) const;
            inline double _phi_perp32_t(const double & q2) const;
            inline double _phi_long12_t5(const double & q2) const;
            inline double _phi_perp12_t5(const double & q2) const;
            inline double _phi_perp32_t5(const double & q2) const;

            // End-point relations
            double _a_long12_v_0() const;
            double _a_perp12_v_0() const;
            double _a_time12_a_0() const;
            double _a_long12_t_0() const;
            double _a_perp12_t_0() const;
            double _a_perp32_t5_0() const;
            double _a_time12_v_0() const;
            double _a_long12_a_0() const;
            double _a_perp12_t5_0() const;
            double _a_perp12_a_0() const;
            double _a_long12_t5_0() const;

        public:
            ABR2022FormFactors(const Parameters & parameters, const Options & options);
            virtual ~ABR2022FormFactors() = default;

            static FormFactors<OneHalfPlusToThreeHalfMinus> * make(const Parameters & parameters, const Options & options);

            virtual double f_time12_v(const double & s) const; // a.k.a zero
            virtual double f_long12_v(const double & s) const; // a.k.a plus
            virtual double f_perp12_v(const double & s) const;
            virtual double f_perp32_v(const double & s) const;

            virtual double f_time12_a(const double & s) const; // a.k.a zero
            virtual double f_long12_a(const double & s) const; // a.k.a plus
            virtual double f_perp12_a(const double & s) const;
            virtual double f_perp32_a(const double & s) const;

            virtual double f_long12_t(const double & s) const; // a.k.a plus
            virtual double f_perp12_t(const double & s) const;
            virtual double f_perp32_t(const double & s) const;

            virtual double f_long12_t5(const double & s) const; // a.k.a plus
            virtual double f_perp12_t5(const double & s) const;
            virtual double f_perp32_t5(const double & s) const;

            double saturation_0p_v() const;
            double saturation_1m_v() const;
            double saturation_0m_a() const;
            double saturation_1p_a() const;
            double saturation_1m_t() const;
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

    extern template class ABR2022FormFactors<LambdaBToLambda1520>;
}

#endif
