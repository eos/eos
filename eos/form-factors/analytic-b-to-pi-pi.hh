/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_PI_PI_HH
#define EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_PI_PI_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    class AnalyticFormFactorBToPiPiBFvD2016 :
        public FormFactors<PToPP>,
        PrivateImplementationPattern<AnalyticFormFactorBToPiPiBFvD2016>
    {
        public:
            AnalyticFormFactorBToPiPiBFvD2016(const Parameters &, const Options &);

            ~AnalyticFormFactorBToPiPiBFvD2016();

            static FormFactors<PToPP> * make(const Parameters &, const Options &);

            /* Leading-order terms */
            complex<double> F_perp_lo_tw2(const double & q2, const double & k2, const double & z) const;
            complex<double> F_perp_lo_tw3(const double & q2, const double & k2, const double & z) const;

            // Partial waves
            virtual std::array<complex<double>, 4> f_perp(const double & q2, const double & k2) const override;
            virtual std::array<complex<double>, 4> f_para(const double & q2, const double & k2) const override;
            virtual std::array<complex<double>, 4> f_long(const double & q2, const double & k2) const override;
            virtual std::array<complex<double>, 4> f_time(const double & q2, const double & k2) const override;

            /* Form factors */
            virtual complex<double> f_perp(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_para(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_long(const double & q2, const double & k2, const double & z) const override;
            virtual complex<double> f_time(const double & q2, const double & k2, const double & z) const override;

            double re_f_perp(const double & q2, const double & k2, const double & z) const;
            double im_f_perp(const double & q2, const double & k2, const double & z) const;
            double re_f_para(const double & q2, const double & k2, const double & z) const;
            double im_f_para(const double & q2, const double & k2, const double & z) const;
            double re_f_long(const double & q2, const double & k2, const double & z) const;
            double im_f_long(const double & q2, const double & k2, const double & z) const;
            double re_f_time(const double & q2, const double & k2, const double & z) const;
            double im_f_time(const double & q2, const double & k2, const double & z) const;

            /* Diagnostics for unit tests */
            Diagnostics diagnostics() const;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };

    class AnalyticFormFactorBToPiPiFvDV2018 :
        public FormFactors<PToPP>,
        PrivateImplementationPattern<AnalyticFormFactorBToPiPiFvDV2018>
    {
        public:
            AnalyticFormFactorBToPiPiFvDV2018(const Parameters &, const Options &);

            ~AnalyticFormFactorBToPiPiFvDV2018();

            static FormFactors<PToPP> * make(const Parameters &, const Options &);

            // Partial waves
            virtual std::array<complex<double>, 4> f_perp(const double & q2, const double & k2) const;
            virtual std::array<complex<double>, 4> f_para(const double & q2, const double & k2) const;
            virtual std::array<complex<double>, 4> f_long(const double & q2, const double & k2) const;
            virtual std::array<complex<double>, 4> f_time(const double & q2, const double & k2) const;

            /* Form factors */
            virtual complex<double> f_perp(const double & q2, const double & k2, const double & z) const;
            virtual complex<double> f_para(const double & q2, const double & k2, const double & z) const;
            virtual complex<double> f_long(const double & q2, const double & k2, const double & z) const;
            virtual complex<double> f_time(const double & q2, const double & k2, const double & z) const;

            /* Form factor residues */
            double f_perp_im_res_qhat2(const double & q2, const double & k2) const;
            double f_para_im_res_qhat2(const double & q2, const double & k2) const;
            double f_long_im_res_qhat2(const double & q2, const double & k2) const;
            double f_time_im_res_qhat2(const double & q2, const double & k2) const;

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
            static const std::vector<OptionSpecification> options;
    };
}

#endif
