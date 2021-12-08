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
    class BMRvD2022FormFactors :
        public FormFactors<OneHalfPlusToOneHalfPlus>
    {
        private:
            const double _m_1, _m_2; // m_1 is the mass of the heavier particle, m_2 the mass of the lighter particle
            const double _t_0, _t_m, _t_p; // z(t_0) = 0, t_m is the endpoint of the semileptonic process, and t_p is the pair production threshold,

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

            QualifiedName _par_name(const std::string & pol, const std::string & current, unsigned idx) const;
            double _z(const double & t, const double & t_0) const;
            double _phi(const double & s, const double & chi, const double & a, const double & b, const double & c,
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

            double bound_0p() const;
            double bound_1m() const;
            double bound_0m() const;
            double bound_1p() const;
            double bound_T() const;
            double bound_T5() const;

            double bound_0p_prior() const;
            double bound_1m_prior() const;
            double bound_0m_prior() const;
            double bound_1p_prior() const;
            double bound_T_prior() const;
            double bound_T5_prior() const;

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
}

#endif
