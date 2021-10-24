/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014, 2015, 2020 Danny van Dyk
 * Copyright (c) 2019, 2020 Domagoj Leljak
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_PI_HH
#define EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_PI_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    class AnalyticFormFactorBToPiDKMMO2008 :
        public FormFactors<PToP>,
        PrivateImplementationPattern<AnalyticFormFactorBToPiDKMMO2008>
    {
        public:
            AnalyticFormFactorBToPiDKMMO2008(const Parameters &, const Options &);

            ~AnalyticFormFactorBToPiDKMMO2008();

            static FormFactors<PToP> * make(const Parameters &, const Options &);

            /* Leading-order terms */
            double F_lo_tw2(const double & q2) const;
            double F_lo_tw3(const double & q2) const;
            double F_lo_tw4(const double & q2) const;
            double Ftil_lo_tw3(const double & q2) const;
            double Ftil_lo_tw4(const double & q2) const;
            double FT_lo_tw2(const double & q2) const;
            double FT_lo_tw3(const double & q2) const;
            double FT_lo_tw4(const double & q2) const;

            /* Next-to-leading-order terms */
            double F_nlo_tw2(const double & q2) const;
            double F_nlo_tw3(const double & q2) const;
            double Ftil_nlo_tw2(const double & q2) const;
            double Ftil_nlo_tw3(const double & q2) const;
            double FT_nlo_tw2(const double & q2) const;
            double FT_nlo_tw3(const double & q2) const;

            /* Form factors */
            virtual double f_p(const double & q2) const;
            virtual double f_0(const double & q2) const;
            virtual double f_t(const double & q2) const;

            virtual double f_plus_T(const double & q2) const;

            /* B mass from the LCSR and the SVZ sum rule, respectively */
            double MBp_lcsr(const double & q2) const;
            double MB0_lcsr(const double & q2) const;
            double MBT_lcsr(const double & q2) const;
            double MB_svz() const;

            /* B meson decay constant at NLO */
            double decay_constant() const;

            /* Diagnostics for unit tests */
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
    };
}

#endif
