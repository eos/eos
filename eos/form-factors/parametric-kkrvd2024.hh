/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

/*
 * Copyright (c) 2020-2024 Danny van Dyk
 * Copyright (c) 2024 Matthew J. Kirk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KKRvD2024_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KKRvD2024_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/reference-name.hh>

#include <array>

namespace eos
{
    template <typename Process_> class KKRvD2024FormFactors;

    template <> class KKRvD2024FormFactors<VacuumToPiPi> :
        public FormFactors<VacuumToPP>
    {
        private:
            // parameters for form factor f_+ (I=1 projection)
            std::array<UsedParameter, 8u> _b_fp_I1;
            UsedParameter _M_fp_I1;
            UsedParameter _G_fp_I1;

            // hadron masses
            UsedParameter _m_pi;

            // parameter for zero point of z
            UsedParameter _t_0;

            UsedParameter _hbar;

            inline std::string _par_name(const std::string & ff, const std::string & isospin, const std::string & index) const
            {
                return "0->pipi::b_(" + ff + "," + isospin + ")^" + index + "@KKRvD2024";
            }

            inline double _hbarc() const
            {
                return _hbar() * 299792458 * 1e15; // GeV fm
            }

            inline double _t_p() const
            {
                return 4.0 * _m_pi() * _m_pi();
            }

            inline complex<double> _z(const complex<double> & q2, const double & t_0) const
            {
                const auto t_p = _t_p();

                return (sqrt(t_p - q2) - sqrt(t_p - t_0)) / (sqrt(t_p - q2) + sqrt(t_p - t_0));
            }

            inline complex<double> _zr(const double & M, const double & Gamma) const
            {
                return 1.0 / _z(power_of<2>(complex<double>(M, -Gamma/2)), _t_0());
            }

            double _b0_fp_I1(const double & chi, const complex<double> & zr) const;
            double _b1_fp_I1(const double & chi, const complex<double> & zr) const;

        public:
            KKRvD2024FormFactors(const Parameters & p, const Options & o);
            ~KKRvD2024FormFactors();

            static FormFactors<VacuumToPP> * make(const Parameters & p, const Options & o);

            /* auxiliary functions */
            complex<double> z(const complex<double> & q2) const;
            complex<double> dzdq2(const complex<double> & q2) const;
            complex<double> dzdq2_II(const complex<double> & q2) const;
            complex<double> w(const complex<double> & z) const;
            complex<double> phitilde_p(const complex<double> & z, const double & chi) const;
            complex<double> phitildeprime_p(const complex<double> & z, const double & chi) const;
            complex<double> series_m(const complex<double> & z, const std::array<double, 10u> & c) const;
            double dFdq2_q2eq0() const;

            /* form factors on the real axis */
            virtual complex<double> f_p(const double & q2) const override;
            virtual complex<double> f_0(const double & q2) const override;
            virtual complex<double> f_t(const double & q2) const override;

            /* form factor in the complex q2 plane */
            virtual complex<double> f_p(const complex<double> & q2) const override;
            virtual complex<double> f_0(const complex<double> & q2) const override;
            virtual complex<double> f_t(const complex<double> & q2) const override;

            /* auxiliary observables */
            double r_pi_squared() const; // squared charge radius of the pion

            /* auxiliary pseudo observables */
            double b_0() const; // value of the series coefficient b_0
            double b_1() const; // value of the series coefficient b_0
            complex<double> residue_rho() const; // complex-valued residue of the form factor on the rho pole
            double re_residue_rho() const;
            double im_residue_rho() const;
            complex<double> residue_rho_q2() const; // complex-valued residue in q2 of the form factor on the rho pole
            double re_residue_rho_q2() const;
            double im_residue_rho_q2() const;

            /* saturation of the dispersive bound */
            double dispersive_integrand(const double & alpha) const;
            double saturation() const;

            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> option_specifications;

            static const std::set<ReferenceName> references;
    };

    extern template class KKRvD2024FormFactors<VacuumToPiPi>;
}

#endif
