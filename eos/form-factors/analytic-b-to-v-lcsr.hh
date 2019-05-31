/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2018 Nico Gubernari
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_V_LCSR_HH
#define EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_V_LCSR_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    namespace lcsr
    {
        struct BToRho;
        struct BToKstar;
        struct BToDstar;
        struct BsToKstar;
        struct BsToPhi;
        struct BsToDsstar;
    }

    template <typename Process_>
    class AnalyticFormFactorBToVLCSR :
        public FormFactors<PToV>,
        PrivateImplementationPattern<AnalyticFormFactorBToVLCSR<Process_>>
    {
        public:
            AnalyticFormFactorBToVLCSR(const Parameters &, const Options &);

            ~AnalyticFormFactorBToVLCSR();

            static FormFactors<PToV> * make(const Parameters &, const Options &);

            /* Form factors */
            virtual double v(const double & q2) const;
            virtual double a_0(const double & q2) const;
            virtual double a_1(const double & q2) const;
            virtual double a_2(const double & q2) const;
            virtual double a_12(const double & q2) const;

            virtual double t_1(const double & q2) const;
            virtual double t_2(const double & q2) const;
            virtual double t_3(const double & q2) const;
            virtual double t_23(const double & q2) const;

            /* First moments of the sum rules */
            double normalized_moment_1_a_1(const double & q2) const;
            double normalized_moment_1_a_2(const double & q2) const;
            double normalized_moment_1_a_30(const double & q2) const;
            double normalized_moment_1_v(const double & q2) const;
            double normalized_moment_1_t_1(const double & q2) const;
            double normalized_moment_1_t_23A(const double & q2) const;
            double normalized_moment_1_t_23B(const double & q2) const;

            /* Diagnostics for unit tests */
            Diagnostics diagnostics() const;
    };

    extern template class AnalyticFormFactorBToVLCSR<lcsr::BToRho>;
    extern template class AnalyticFormFactorBToVLCSR<lcsr::BToKstar>;
    extern template class AnalyticFormFactorBToVLCSR<lcsr::BToDstar>;
    extern template class AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>;
    extern template class AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>;
    extern template class AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>;
}
#endif
