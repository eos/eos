/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Stephan Kürten
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KKVDZ2022_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KKVDZ2022_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/reference-name.hh>

#include <memory>

namespace eos
{
    class KKvDZ2022FormFactors :
        public FormFactors<PToGammaOffShell>
    {
        private:
            SwitchOption opt_subtracted;
            double switch_subtracted;
            UsedParameter _s_0;
            static std::string _par_name(const std::string & ff_name);
            double _width_omega(const double & q2) const;
            double _width_rho(const double & q2) const;
            double _z_omega(const double & k2) const;
            double _z_rho(const double & k2) const;
            double _subtraction_polynomial(const double & k2, const std::array<UsedParameter, 3> & c) const;
            complex<double> _calc_ff_contribution_omega(const double & q2, const double & k2,
                const double & m_r_sq, const std::array<UsedParameter, 3> & a,
                const UsedParameter & s_0) const;
            complex<double> _calc_ff_contribution_rho(const double & q2, const double & k2,
                const double & m_r_sq, const std::array<UsedParameter, 3> & a,
                const UsedParameter & s_0) const;
            std::array<std::array<UsedParameter, 3>, 4> _a_omega;
            std::array<std::array<UsedParameter, 3>, 4> _a_rho;
            std::array<std::array<UsedParameter, 3>, 4> _c_subtraction;

        public:
            KKvDZ2022FormFactors(const Parameters &, const Options &);

            static FormFactors<PToGammaOffShell> * make(const Parameters &, const Options &);

            virtual complex<double> F_1(const double & q2, const double & k2) const override;
            virtual complex<double> F_2(const double & q2, const double & k2) const override;
            virtual complex<double> F_3(const double & q2, const double & k2) const override;
            virtual complex<double> F_4(const double & q2, const double & k2) const override;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> options;
    };
}

#endif
