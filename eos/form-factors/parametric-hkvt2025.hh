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
#include <eos/utils/kinematic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>
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
            UsedParameter m_P1, m_P2, m_P3;
            UsedParameter m_R_0m, m_R_1m, m_R_1p;
            UsedParameter q2p_a, q2p_v, q20, sin0, sin1, s0;

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_0m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1m_names;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string> resonance_1p_names;

            HKVT2025FormFactorTraits(const Parameters & p) :
                m_P1(UsedParameter(p[std::string(Process_::name_P1) + "@HKvT2025"], *this)),
                m_P2(UsedParameter(p[std::string(Process_::name_P2) + "@HKvT2025"], *this)),
                m_P3(UsedParameter(p[std::string(Process_::name_P3) + "@HKvT2025"], *this)),
                m_R_0m(UsedParameter(p[resonance_0m_names.at(Process_::partonic_transition)], *this)),
                m_R_1m(UsedParameter(p[resonance_1m_names.at(Process_::partonic_transition)], *this)),
                m_R_1p(UsedParameter(p[resonance_1p_names.at(Process_::partonic_transition)], *this)),
                q2p_a(UsedParameter(p[std::string(Process_::label) + "::tp_a@HKvT2025"], *this)),
                q2p_v(UsedParameter(p[std::string(Process_::label) + "::tp_v@HKvT2025"], *this)),
                q20(UsedParameter(p[std::string(Process_::label) + "::t0@HKvT2025"], *this)),
                sin0(UsedParameter(p[std::string(Process_::label) + "::sin0@HKvT2025"], *this)),
                sin1(UsedParameter(p[std::string(Process_::label) + "::sin1@HKvT2025"], *this)),
                s0(UsedParameter(p[std::string(Process_::label) + "::s0@HKvT2025"], *this))
            {
            }

            double tm(const double & s) const
            {
                return power_of<2>(m_P1 - std::sqrt(s));
            }

            double kappa(const double & s, const double & q2) const
            {
                const double lamq3 = (q2 - power_of<2>(m_P1 + std::sqrt(s))) * (q2 - power_of<2>(m_P1 - std::sqrt(s)));
                const double lams12 = (s - power_of<2>(m_P2 + m_P3)) * (s - power_of<2>(m_P2 - m_P3));
                if (lamq3 < 0.0 || lams12 < 0.0)
                    return 0.0;
                return std::sqrt(lams12 * lamq3) / s;
            }

            complex<double> calc_y(const complex<double> & s, const complex<double> & sin, const complex<double> & s0) const
            {
                return (std::sqrt(sin - s) - std::sqrt(sin - s0)) / (std::sqrt(sin - s) + std::sqrt(sin - s0));
            }

            complex<double> calc_z(const complex<double> & q2, const complex<double> & q2p, const complex<double> & q20) const
            {
                return (std::sqrt(q2p - q2) - std::sqrt(q2p - q20)) / (std::sqrt(q2p - q2) + std::sqrt(q2p - q20));
            }

            double calc_y(const double & s, const double & sin, const double & s0) const
            {
                if (s > sin)
                    throw InternalError("The real conformal mapping is used above threshold: " + stringify(s) + " > " + stringify(sin));

                return real(calc_y(complex<double>(s, 0.0), complex<double>(sin, 0.0), complex<double>(s0, 0.0)));
            }

            double calc_z(const double & q2, const double & q2p, const double & q20) const
            {
                if (q2 > q2p)
                    throw InternalError("The real conformal mapping is used above threshold: " + stringify(q2) + " > " + stringify(q2p));

                return real(calc_z(complex<double>(q2, 0.0), complex<double>(q2p, 0.0), complex<double>(q20, 0.0)));
            }

            std::array<double, 3> orthonormal_polynomials_v(const double & z, const double & s) const
            {
                const double measure = 2 * std::arg(calc_z(complex<double>(power_of<2>(m_P1 + std::sqrt(s))), complex<double>(q2p_v), complex<double>(q20)));
                const SzegoPolynomial<2> polynomials_set(SzegoPolynomial<2>::FlatMeasure(measure));

                return polynomials_set(z);
            }

            std::array<double, 3> orthonormal_polynomials_a(const double & z, const double & s) const
            {
                const double measure = 2 * std::arg(calc_z(complex<double>(power_of<2>(m_P1 + std::sqrt(s))), complex<double>(q2p_a), complex<double>(q20)));
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
        public FormFactors<PToPP>
    {
        private:
            // We want more expansion terms in z than in y for now
            // For most processes, we only have one possible isospin configuration in the final state, not so for Pi Pi, where we have isoscalar and isovector
            std::vector<std::array< std::array<UsedParameter, 3>, 3>> _a_g_1, _a_f_1, _a_F1_1, _a_F2_1;
            std::vector<std::array< std::array<UsedParameter, 3>, 3>> _a_g_2, _a_f_2, _a_F1_2, _a_F2_2;

            const HKVT2025FormFactorTraits<Process_, PToPP> _traits;

            const UsedParameter & _mP1, _mP2, _mP3;
            const unsigned _numWaves;

            QualifiedName _par_name(const std::string & ff_name, unsigned wave, unsigned idx) const;
            double _phi(const int & l, const double & s, const double & q2, const double & q2p, const double & chi,
                        const unsigned a, const unsigned b, const double N) const;

            inline double _phi_g(const double & q2, const double & s, const unsigned & l) const;
            inline double _phi_f(const double & q2, const double & s, const unsigned & l) const;
            inline double _phi_F1(const double & q2, const double & s, const unsigned & l) const;
            inline double _phi_F2(const double & q2, const double & s, const unsigned & l) const;

        public:
            HKVT2025FormFactors(const Parameters & p, const Options &);

            ~HKVT2025FormFactors();

            static FormFactors<PToPP> * make(const Parameters & parameters, const Options & options);

            virtual complex<double> v_perp(const double & q2, const double & s, const unsigned & l, const bool & iso) const override;
            virtual complex<double> a_t(const double & q2, const double & s, const unsigned & l, const bool & iso) const override;
            virtual complex<double> a_0(const double & q2, const double & s, const unsigned & l, const bool & iso) const override;
            virtual complex<double> a_par(const double & q2, const double & s, const unsigned & l, const bool & iso) const override;

            virtual unsigned numWaves() const override { return _numWaves; }

            virtual double unitarity_integrand_0m(const double & s, const unsigned & l, const bool & iso) const override;
            virtual double unitarity_integrand_1p(const double & s, const unsigned & l, const bool & iso) const override;
            virtual double unitarity_integrand_1m(const double & s, const unsigned & l, const bool & iso) const override;

            /* Placeholders */
            virtual complex<double> f_perp(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_para(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_long(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_time(const double & q2, const double & k2, const double & z) const override;

            virtual double f_perp_im_res_qhat2(const double & q2, const double & k2) const override;
            virtual double f_para_im_res_qhat2(const double & q2, const double & k2) const override;
            virtual double f_long_im_res_qhat2(const double & q2, const double & k2) const override;
            virtual double f_time_im_res_qhat2(const double & q2, const double & k2) const override;

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
