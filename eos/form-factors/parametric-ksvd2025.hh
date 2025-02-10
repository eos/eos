/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

/*
 * Copyright (c) 2025 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KSVD2025_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KSVD2025_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/reference-name.hh>

#include <array>

namespace eos
{
    template <typename Process_> class KSvD2025FormFactors;

    template <> class KSvD2025FormFactors<VacuumToKPi> :
        public FormFactors<VacuumToPP>
    {
        private:
            // parameters for form factor f_+
            std::array<UsedParameter, 8u> _b_fp;
            std::array<UsedParameter, 2u> _M_fp;
            std::array<UsedParameter, 2u> _G_fp;

            // hadron masses
            UsedParameter _m_K;
            UsedParameter _m_pi;

            // parameter for zero point of z
            UsedParameter _t_0;

            UsedParameter _hbar;

            inline std::string _coeff_name(const std::string & ff, const std::string & index) const
            {
                return "0->Kpi::b_" + ff + "^" + index + "@KSvD2025";
            }

            inline double _hbarc() const
            {
                return _hbar() * 299792458 * 1e15; // GeV fm
            }

            inline double _t_p() const
            {
                return power_of<2>(_m_pi() + _m_K());
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

            double _b0_fp(const double & chi_p) const;
            double _b0_fz(const double & chi_p, const double & chi_z) const;

        public:
            KSvD2025FormFactors(const Parameters & p, const Options & o);
            ~KSvD2025FormFactors();

            static FormFactors<VacuumToPP> * make(const Parameters & p, const Options & o);

            /* auxiliary functions */
            complex<double> z(const complex<double> & q2) const;
            complex<double> dzdq2(const complex<double> & q2) const;
            complex<double> dzdq2_II(const complex<double> & q2) const;
            complex<double> series_m(const complex<double> & z, const std::array<double, 10u> & c) const;

            /* functions pertaining to f_p */
            complex<double> w_p(const complex<double> & z) const;
            complex<double> phitilde_p(const complex<double> & z, const double & chi) const;
            complex<double> phitildeprime_p(const complex<double> & z, const double & chi) const;

            /* functions pertaining to f_z */
            complex<double> w_z(const complex<double> & z) const;
            complex<double> phitilde_z(const complex<double> & z, const double & chi) const;
            complex<double> phitildeprime_z(const complex<double> & z, const double & chi) const;

            /* form factors on the real axis */
            virtual complex<double> f_p(const double & q2) const override;
            virtual complex<double> f_0(const double & q2) const override;
            virtual complex<double> f_t(const double & q2) const override;

            /* form factor in the complex q2 plane */
            virtual complex<double> f_p(const complex<double> & q2) const override;
            virtual complex<double> f_0(const complex<double> & q2) const override;
            virtual complex<double> f_t(const complex<double> & q2) const override;

            /* saturation of the dispersive bound */
            double dispersive_integrand(const double & alpha) const;
            double saturation() const;

            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> option_specifications;

            static const std::set<ReferenceName> references;
    };

    extern template class KSvD2025FormFactors<VacuumToKPi>;
}

#endif
