/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_MPV2025_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_MPV2025_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    class MPV2025FormFactors :
        public FormFactors<PToPP>,
        PrivateImplementationPattern<MPV2025FormFactors>
    {
        public:
            MPV2025FormFactors(const Parameters &, const Options &);

            ~MPV2025FormFactors();

            static FormFactors<PToPP> * make(const Parameters &, const Options &);

            /* Form factors */
            virtual complex<double> f_perp(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_para(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_long(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_time(const double & q2, const double & k2, const double & z) const override;

            /* Diagnostics for unit tests */
            Diagnostics diagnostics() const;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };
}

#endif
