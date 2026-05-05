/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Nico Gubernari
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_G2026_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_G2026_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>
#include <eos/utils/reference-name.hh>

#include <array>
#include <vector>

namespace eos
{
    /* Form Factors according to [G:2026A] */
    template <typename Process_, typename Transition_> class G2026FormFactors;

    template <typename Process_, typename Transition_> class G2026FormFactorTraits;


    // P -> V
    template <typename Process_>
    class G2026FormFactorTraits<Process_, PToV> :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and should match the
            // the ones used for the extraction of the coefficients of the z-expansion
            UsedParameter m_B, m_V;
            UsedParameter sV, sA, s0V, s0A, Q2, eta;
            UsedParameter tchi_A0, tchi_A1, tchi_V1, tchi_T1, tchi_AT1; //tchi_1m_v, tchi_0m_a, tchi_1p_a, tchi_1m_t, tchi_1p_t5;

            std::vector<double> m_R_A0, m_R_V1, m_R_A1;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>> pole_A0_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>> pole_V1_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>> pole_A1_names;


            G2026FormFactorTraits(const Parameters & p) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@BSZ2015"], *this)),
                m_V(UsedParameter(p[std::string(Process_::name_V) + "@BSZ2015"], *this)),
                sV(UsedParameter(p[std::string(Process_::label) + "::tp_v@G2026"], *this)),
                sA(UsedParameter(p[std::string(Process_::label) + "::tp_a@G2026"], *this)),
                s0V(UsedParameter(p[std::string(Process_::label) + "::t0_v@G2026"], *this)),
                s0A(UsedParameter(p[std::string(Process_::label) + "::t0_a@G2026"], *this)),
                Q2(UsedParameter(p[std::string(Process_::label) + "::Q2@G2026"], *this)),
                eta(UsedParameter(p[std::string(Process_::label) + "::eta@G2026"], *this)),
                tchi_A0(UsedParameter(p[std::string(Process_::label) + "::tchi_A0@G2026"], *this)),
                tchi_A1(UsedParameter(p[std::string(Process_::label) + "::tchi_A1@G2026"], *this)),
                tchi_V1(UsedParameter(p[std::string(Process_::label) + "::tchi_V1@G2026"], *this)),
                tchi_T1(UsedParameter(p[std::string(Process_::label) + "::tchi_T1@G2026"], *this)),
                tchi_AT1(UsedParameter(p[std::string(Process_::label) + "::tchi_AT1@G2026"], *this)),
                m_R_A0(_resonance_masses_or_empty(*this, p, pole_A0_names, Process_::partonic_transition)),
                m_R_V1(_resonance_masses_or_empty(*this, p, pole_V1_names, Process_::partonic_transition)),
                m_R_A1(_resonance_masses_or_empty(*this, p, pole_A1_names, Process_::partonic_transition))
            {
            }

            double sm() const
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

            std::array<double, 6> monomials_v(const double & z) const
            {
                return { 1.0, z, power_of<2>(z), power_of<3>(z), power_of<4>(z), power_of<5>(z) };
            }

            std::array<double, 6> monomials_a(const double & z) const
            {
                return { 1.0, z, power_of<2>(z), power_of<3>(z), power_of<4>(z), power_of<5>(z) };
            }

            double blaschke_product(const double & q2, const double & sG, const std::vector<double> & masses) const
            {
                double product = 1.0;
                for (const auto & m : masses)
                {
                    if (power_of<2>(m) < sG)
                        product *= calc_z(q2, sG, power_of<2>(m));
                }
                return product;
            }

        private:
            static std::vector<double> _resonance_masses_or_empty(
                ParameterUser & user,
                const Parameters & p,
                const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>> & names,
                const std::tuple<QuarkFlavor, QuarkFlavor> & key)
            {
                const auto it = names.find(key);
                if (it == names.end())
                    return {};

                std::vector<double> result;
                for (const auto & name : it->second)
                {
                    if (const_cast<Parameters &>(p).has(name))
                    {
                        const auto parameter = p[name];
                        user.uses(parameter.id());
                        result.push_back(parameter.evaluate());
                    }
                }
                return result;
            }
    };

    template <typename Process_> class G2026FormFactors<Process_, PToV> :
        public FormFactors<PToV>
    {
        private:
            // fit parametrization for P -> V according to [G:2026A]
            std::array<UsedParameter, 5> _a_A0, _a_V, _a_T1;
            // use endpoint relations (see eq. (3.2) in [HLMW:2015A]) to remove parameters
            std::array<UsedParameter, 4> _a_A12, _a_T2, _a_A1, _a_T23;

            const G2026FormFactorTraits<Process_, PToV> _traits;

            const UsedParameter & _mB, _mV;

            QualifiedName _par_name(const std::string & ff_name, unsigned idx) const;

            double _phi(
                const double & s,
                const double & sG,
                const double & s0,
                const std::vector<double> & Mres,
                const double & tchi,
                const unsigned Kn,
                const int Ksp,
                const int Ksm,
                const int Kspm,
                const unsigned a,
                const unsigned b,
                const unsigned c,
                const unsigned d
            ) const;

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
            G2026FormFactors(const Parameters & p, const Options &);

            ~G2026FormFactors();

            static FormFactors<PToV> * make(const Parameters & parameters, const Options & options);

            virtual double v(const double & q2) const;
            virtual double a_0(const double & q2) const;
            virtual double a_1(const double & q2) const;
            virtual double a_12(const double & q2) const;
            virtual double a_2(const double & q2) const;
            virtual double t_1(const double & q2) const;
            virtual double t_2(const double & q2) const;
            virtual double t_23(const double & q2) const;
            virtual double t_3(const double & q2) const;

            // Unused but needed to satisfy the FormFactors interface
            virtual double f_perp(const double & q2) const;
            virtual double f_para(const double & q2) const;
            virtual double f_long(const double & q2) const;
            virtual double f_perp_T(const double & q2) const;
            virtual double f_para_T(const double & q2) const;
            virtual double f_long_T(const double & q2) const;

            // Saturations of the dispersive bounds.
            // EOS channel labels; in [G:2026A] these are denoted
            // 0p_v = V0, 0m_a = A0, 1m_v = V1, 1p_a = A1, 1m_t = T1, 1p_t5 = AT1.
            // J = 0
            double saturation_0p_v() const;   // [G:2026A]: saturation_V0
            double saturation_0m_a() const;   // [G:2026A]: saturation_A0
            // J = 1
            double saturation_1m_v() const;   // [G:2026A]: saturation_V1
            double saturation_1p_a() const;   // [G:2026A]: saturation_A1
            double saturation_1m_t() const;   // [G:2026A]: saturation_T1
            double saturation_1p_t5() const;  // [G:2026A]: saturation_AT1

            Diagnostics diagnostics() const;


            // Auxilliary functions: series and derivative of the series
            double v_series(const double & q2) const;
            double a_0_series(const double & q2) const;
            double a_1_series(const double & q2) const;
            double a_12_series(const double & q2) const;
            double t_1_series(const double & q2) const;
            double t_2_series(const double & q2) const;
            double t_23_series(const double & q2) const;

            double v_series_prime(const double & q2) const;
            double a_0_series_prime(const double & q2) const;
            double a_1_series_prime(const double & q2) const;
            double a_12_series_prime(const double & q2) const;
            double t_1_series_prime(const double & q2) const;
            double t_2_series_prime(const double & q2) const;
            double t_23_series_prime(const double & q2) const;

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

    extern template class G2026FormFactors<BsToKstar,  PToV>;
    extern template class G2026FormFactors<BToKstar,   PToV>;
    extern template class G2026FormFactors<BsToPhi,    PToV>;
    extern template class G2026FormFactors<BToDstar,   PToV>;
    extern template class G2026FormFactors<BsToDsstar, PToV>;


    // P -> P
    template <typename Process_>
    class G2026FormFactorTraits<Process_, PToP> :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and should match the
            // the ones used for the extraction of the coefficients of the z-expansion
            UsedParameter m_B, m_P;
            UsedParameter sV, s0V, Q2, eta;
            UsedParameter tchi_V0, tchi_V1, tchi_T1; // tchi_1m_v, tchi_0p_v, tchi_1m_t

            std::vector<double> m_R_V0, m_R_V1;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>> pole_V0_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>> pole_V1_names;

            G2026FormFactorTraits(const Parameters & p) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@BSZ2015"], *this)),
                m_P(UsedParameter(p[std::string(Process_::name_P) + "@BSZ2015"], *this)),
                sV(UsedParameter(p[std::string(Process_::label) + "::tp_v@G2026"], *this)),
                s0V(UsedParameter(p[std::string(Process_::label) + "::t0_v@G2026"], *this)),
                Q2(UsedParameter(p[std::string(Process_::label) + "::Q2@G2026"], *this)),
                eta(UsedParameter(p[std::string(Process_::label) + "::eta@G2026"], *this)),
                tchi_V0(UsedParameter(p[std::string(Process_::label) + "::tchi_V0@G2026"], *this)),
                tchi_V1(UsedParameter(p[std::string(Process_::label) + "::tchi_V1@G2026"], *this)),
                tchi_T1(UsedParameter(p[std::string(Process_::label) + "::tchi_T1@G2026"], *this)),
                m_R_V0(_resonance_masses_or_empty(*this, p, pole_V0_names, Process_::partonic_transition)),
                m_R_V1(_resonance_masses_or_empty(*this, p, pole_V1_names, Process_::partonic_transition))
            {
            }

            double sm() const
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

            std::array<double, 6> monomials(const double & z) const
            {
                return { 1.0, z, power_of<2>(z), power_of<3>(z), power_of<4>(z), power_of<5>(z) };
            }

            std::array<complex<double>, 6> monomials(const complex<double> & z) const
            {
                return { 1.0, z, power_of<2>(z), power_of<3>(z), power_of<4>(z), power_of<5>(z) };
            }

            std::array<complex<double>, 6> monomials_derivatives(const complex<double> & z) const
            {
                return {
                    complex<double>(0.0, 0.0), 1.0,
                    2.0 * z,
                    3.0 * power_of<2>(z),
                    4.0 * power_of<3>(z),
                    5.0 * power_of<4>(z)
                };
            }

            double blaschke_product(const double & q2, const double & sG, const std::vector<double> & masses) const
            {
                double product = 1.0;
                for (const auto & m : masses)
                {
                    if (power_of<2>(m) < sG)
                        product *= calc_z(q2, sG, power_of<2>(m));
                }
                return product;
            }

        private:
            static std::vector<double> _resonance_masses_or_empty(
                ParameterUser & user,
                const Parameters & p,
                const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>> & names,
                const std::tuple<QuarkFlavor, QuarkFlavor> & key)
            {
                const auto it = names.find(key);
                if (it == names.end())
                    return {};

                std::vector<double> result;
                for (const auto & name : it->second)
                {
                    if (const_cast<Parameters &>(p).has(name))
                    {
                        const auto parameter = p[name];
                        user.uses(parameter.id());
                        result.push_back(parameter.evaluate());
                    }
                }
                return result;
            }

    };

    template <typename Process_> class G2026FormFactors<Process_, PToP> :
        public FormFactors<PToP>
    {
        private:
            // fit parametrization for P -> P inspired by [G:2026A]
            std::array<UsedParameter, 5> _a_fp, _a_ft;
            // use equation of motion to remove f_0(0) as a free parameter
            std::array<UsedParameter, 4> _a_f0;

            const G2026FormFactorTraits<Process_, PToP> _traits;

            const UsedParameter & _mB, _mP;

            QualifiedName _par_name(const std::string & ff_name, unsigned idx) const;

            double _phi(
                const double & s,
                const double & sG,
                const std::vector<double> & Mres,
                const double & tchi,
                const unsigned Kn,
                const int Ksp,
                const int Ksm,
                const int Kspm,
                const unsigned a,
                const unsigned b,
                const unsigned c,
                const unsigned d
            ) const;

            inline double _phi_f_p(const double & q2) const;
            inline double _phi_f_0(const double & q2) const;
            inline double _phi_f_t(const double & q2) const;

            // End-point relations
            double _a_f0_0() const;

        public:
            G2026FormFactors(const Parameters & p, const Options &);

            ~G2026FormFactors();

            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options);

            virtual double f_p(const double & q2) const;
            virtual double f_0(const double & q2) const;
            virtual double f_t(const double & q2) const;

            // Unused but needed to satisfy the FormFactors interface
            virtual double f_plus_T(const double & q2) const;

            // Saturations of the dispersive bounds.
            // EOS channel labels; in [G:2026A] these are denoted
            // 0p_v = V0, 0m_a = A0, 1m_v = V1, 1p_a = A1, 1m_t = T1, 1p_t5 = AT1.
            // J = 0
            double saturation_0p_v() const;   // [G:2026A]: saturation_V0
            double saturation_0m_a() const;   // [G:2026A]: saturation_A0
            // J = 1
            double saturation_1m_v() const;   // [G:2026A]: saturation_V1
            double saturation_1p_a() const;   // [G:2026A]: saturation_A1
            double saturation_1m_t() const;   // [G:2026A]: saturation_T1
            double saturation_1p_t5() const;  // [G:2026A]: saturation_AT1

            Diagnostics diagnostics() const;

            // Auxilliary functions: series and derivative of the series
            double f_p_series(const double & q2) const;
            double f_0_series(const double & q2) const;
            double f_t_series(const double & q2) const;

            double f_p_series_prime(const double & q2) const;
            double f_0_series_prime(const double & q2) const;
            double f_t_series_prime(const double & q2) const;

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

    extern template class G2026FormFactors<BToPi,  PToP>;
    extern template class G2026FormFactors<BsToK,  PToP>;
    extern template class G2026FormFactors<BToK,   PToP>;
    extern template class G2026FormFactors<BToD,   PToP>;
    extern template class G2026FormFactors<BsToDs, PToP>;
}

#endif
