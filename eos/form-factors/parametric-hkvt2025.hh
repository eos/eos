/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_HKVT2025_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_HKVT2025_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/power-of.hh>
#include <eos/maths/szego-polynomial.hh>
#include <eos/scattering/scattering-amplitudes.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/reference-name.hh>

#include <array>

namespace eos
{
    /* Form Factors according to [HKvT:2025A] */
    template <typename Process_, typename Transition_> class HKVT2025FormFactors;

    template <typename Process_, typename Transition_> class HKVT2025FormFactorTraits;


    // P -> PP
    template <typename Process_>
    class HKVT2025FormFactorTraits<Process_, PToPP> :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and should match the
            // the ones used for the extraction of the coefficients of the x-y-z-expansion
            UsedParameter m_B, m_P1, m_P2;
            UsedParameter m_R_0m, m_R_1m, m_R_1p;
            UsedParameter chi_0m_a, chi_1m_v, chi_1p_a;
            UsedParameter q2p_a, q2p_v, q20, k20;
            std::array<UsedParameter, 2> k2in;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1p_names;
            static const std::map<std::string, unsigned> charge_map;

            HKVT2025FormFactorTraits(const Parameters & p) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@HKvT2025"], *this)),
                m_P1(UsedParameter(p[std::string(Process_::name_P1) + "@HKvT2025"], *this)),
                m_P2(UsedParameter(p[std::string(Process_::name_P2) + "@HKvT2025"], *this)),
                m_R_0m(UsedParameter(p[resonance_0m_names.at(Process_::partonic_transition)], *this)),
                m_R_1m(UsedParameter(p[resonance_1m_names.at(Process_::partonic_transition)], *this)),
                m_R_1p(UsedParameter(p[resonance_1p_names.at(Process_::partonic_transition)], *this)),
                chi_0m_a(UsedParameter(p[std::string(Process_::label) + "::chi_0m_a@HKvT2025"], *this)),
                chi_1m_v(UsedParameter(p[std::string(Process_::label) + "::chi_1m_v@HKvT2025"], *this)),
                chi_1p_a(UsedParameter(p[std::string(Process_::label) + "::chi_1p_a@HKvT2025"], *this)),
                q2p_a(UsedParameter(p[std::string(Process_::label) + "::tp_a@HKvT2025"], *this)),
                q2p_v(UsedParameter(p[std::string(Process_::label) + "::tp_v@HKvT2025"], *this)),
                q20(UsedParameter(p[std::string(Process_::label) + "::t0@HKvT2025"], *this)),
                k20(UsedParameter(p[std::string(Process_::label) + "::s0@HKvT2025"], *this)),
                k2in{ UsedParameter(p[std::string(Process_::label) + "::sin0@HKvT2025"], *this), UsedParameter(p[std::string(Process_::label) + "::sin1@HKvT2025"], *this) }
            {
            }

            double tm(const double & k2) const
            {
                return power_of<2>(m_B - std::sqrt(k2));
            }

            double lam_b(const double & k2, const double & q2) const
            {
                return (q2 - power_of<2>(m_B + std::sqrt(k2))) * (q2 - power_of<2>(m_B - std::sqrt(k2)));
            }

            double kappa(const double & k2, const double & q2) const
            {
                const double lamq3 = lam_b(k2, q2);
                const double lams12 = (k2 - power_of<2>(m_P1 + m_P2)) * (k2 - power_of<2>(m_P1 - m_P2));
                if (lamq3 < 0.0 || lams12 < 0.0)
                    return 0.0;
                return std::sqrt(lams12 * lamq3) / k2;
            }

            complex<double> calc_y(const complex<double> & k2, const complex<double> & k2in, const complex<double> & k20) const
            {
                return (std::sqrt(k2in - k2) - std::sqrt(k2in - k20)) / (std::sqrt(k2in - k2) + std::sqrt(k2in - k20));
            }

            complex<double> calc_z(const complex<double> & q2, const complex<double> & q2p, const complex<double> & q20) const
            {
                return (std::sqrt(q2p - q2) - std::sqrt(q2p - q20)) / (std::sqrt(q2p - q2) + std::sqrt(q2p - q20));
            }

            double calc_y(const double & k2, const double & k2in, const double & k20) const
            {
                if (k2 > k2in)
                    throw InternalError("The real conformal mapping is used above threshold: " + stringify(k2) + " > " + stringify(k2in));

                return real(calc_y(complex<double>(k2, 0.0), complex<double>(k2in, 0.0), complex<double>(k20, 0.0)));
            }

            double calc_z(const double & q2, const double & q2p, const double & q20) const
            {
                if (q2 > q2p)
                    throw InternalError("The real conformal mapping is used above threshold: " + stringify(q2) + " > " + stringify(q2p));

                return real(calc_z(complex<double>(q2, 0.0), complex<double>(q2p, 0.0), complex<double>(q20, 0.0)));
            }

            std::array<double, 3> orthonormal_polynomials_v(const double & z, const double & k2) const
            {
                const double measure = 2 * std::arg(calc_z(complex<double>(power_of<2>(m_B + std::sqrt(k2))), complex<double>(q2p_v), complex<double>(q20)));
                const SzegoPolynomial<2> polynomials_set(SzegoPolynomial<2>::FlatMeasure(measure));

                return polynomials_set(z);
            }

            std::array<double, 3> orthonormal_polynomials_a(const double & z, const double & k2) const
            {
                const double measure = 2 * std::arg(calc_z(complex<double>(power_of<2>(m_B + std::sqrt(k2))), complex<double>(q2p_a), complex<double>(q20)));
                const SzegoPolynomial<2> polynomials_set(SzegoPolynomial<2>::FlatMeasure(measure));

                return polynomials_set(z);
            }

            std::array<complex<double>, 3> threshold_improved_polynomials(const complex<double> & y, const unsigned & l) const
            {
                if (l == 1)
                {
                    return std::array<complex<double>, 3>{ 1.0, y - power_of<3>(y) / 3.0, power_of<2>(y) + 2.0 * power_of<3>(y) / 3.0 };
                }
                else if (l == 2)
                {
                    return std::array<complex<double>, 3>{ 1.0, y - (15.0 * power_of<3>(y) / 3.0 + 2.0 * power_of<4>(y)) / 7.0, power_of<2>(y) + 2.0 * (4.0 * power_of<3>(y) + 5.0 * power_of<4>(y) / 4.0) / 7.0 };
                }
                return std::array<complex<double>, 3>{ 1.0, y, power_of<2>(y) };
            }

            std::array<double, 3> threshold_improved_polynomials(const double & y, const unsigned & l) const
            {
                if (l == 1)
                {
                    return std::array<double, 3>{ 1.0, y - power_of<3>(y) / 3.0, power_of<2>(y) + 2.0 * power_of<3>(y) / 3.0 };
                }
                else if (l == 2)
                {
                    return std::array<double, 3>{ 1.0, y - (15.0 * power_of<3>(y) / 3.0 + 2.0 * power_of<4>(y)) / 7.0, power_of<2>(y) + 2.0 * (4.0 * power_of<3>(y) + 5.0 * power_of<4>(y) / 4.0) / 7.0 };
                }
                return std::array<double, 3>{ 1.0, y, power_of<2>(y) };
            }
    };

    template <typename Process_> class HKVT2025FormFactors<Process_, PToPP> :
        public FormFactors<PToPP>, public virtual ParameterUser
    {
        private:
            // Switches for enabled partial waves and isospin configurations
            std::array<double, 4> switch_L;
            std::array<double, 2> switch_I;

            // For most processes, we only have one possible isospin configuration in the final state, not so for Pi Pi, where we have isoscalar and isovector
            const std::array<std::array<std::array< std::array<UsedParameter, 3>, 3>, 3>, 2> _a_g, _a_f, _a_F1, _a_F2;

            const HKVT2025FormFactorTraits<Process_, PToPP> traits;

            const UsedParameter & mB, mP1, mP2;

            IsospinOption opt_I;
            PartialWaveOption opt_L;
            SwitchOption opt_C;
            IntegerOption opt_int_points;

            std::shared_ptr<ScatteringAmplitudes<PPToPP>> scattering_amplitudes;

            const unsigned charge;

            QualifiedName _exp_par_name(const std::string & ff_name, unsigned iso, unsigned wave, unsigned z_order, unsigned y_order) const;
            double _phi(const int & l, const double & k2, const double & q2, const double & q2p, const double & chi,
                        const unsigned a, const unsigned b, const double N) const;

            inline double _phi_g(const double & q2, const double & k2, const unsigned & l) const;
            inline double _phi_f(const double & q2, const double & k2, const unsigned & l) const;
            inline double _phi_F1(const double & q2, const double & k2, const unsigned & l) const;
            inline double _phi_F2(const double & q2, const double & k2, const unsigned & l) const;

            double _unitarity_integrand_0m(const double & s, const unsigned & l, const unsigned & iso) const;
            double _unitarity_integrand_1p(const double & s, const unsigned & l, const unsigned & iso) const;
            double _unitarity_integrand_1m(const double & s, const unsigned & l, const unsigned & iso) const;

        public:
            HKVT2025FormFactors(const Parameters & p, const Options & o);

            ~HKVT2025FormFactors();

            static FormFactors<PToPP> * make(const Parameters & parameters, const Options & options);

            complex<double> v_perp(const double & q2, const double & k2, const unsigned & l, const unsigned & iso) const;
            complex<double> a_t(const double & q2, const double & k2, const unsigned & l, const unsigned & iso) const;
            complex<double> a_0(const double & q2, const double & k2, const unsigned & l, const unsigned & iso) const;
            complex<double> a_par(const double & q2, const double & k2, const unsigned & l, const unsigned & iso) const;

            complex<double> g_tilde(const double & q2, const double & k2, const unsigned & l, const unsigned & iso) const;
            complex<double> F2_tilde(const double & q2, const double & k2, const unsigned & l, const unsigned & iso) const;
            complex<double> F1_tilde(const double & q2, const double & k2, const unsigned & l, const unsigned & iso) const;
            complex<double> f_tilde(const double & q2, const double & k2, const unsigned & l, const unsigned & iso) const;

            // Saturations of the dispersive bounds
            // J = 0
            double saturation_0p_v() const;
            double saturation_0m_a() const;
            // J = 1
            double saturation_1m_v() const;
            double saturation_1p_a() const;

            virtual std::array<complex<double>, 4> f_perp(const double & q2, const double & k2) const override;
            virtual std::array<complex<double>, 4> f_para(const double & q2, const double & k2) const override;
            virtual std::array<complex<double>, 4> f_long(const double & q2, const double & k2) const override;
            virtual std::array<complex<double>, 4> f_time(const double & q2, const double & k2) const override;

            virtual complex<double> f_perp(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_para(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_long(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_time(const double & q2, const double & k2, const double & z) const override;

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

    extern template class HKVT2025FormFactors<BToPiPi, PToPP>;

}

#endif
