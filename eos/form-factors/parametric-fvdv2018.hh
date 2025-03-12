/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2022 Danny van Dyk
 * Copyright (c) 2018 Keri Vos
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_FVDV2018_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_FVDV2018_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/utils/options.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    template <typename Process_> class FvDV2018FormFactors :
        public FormFactors<PToPP>
    {
        private:
            UsedParameter _a_Fperp_0_0, _a_Fperp_0_1, _a_Fperp_0_2, _a_Fperp_0_3, _a_Fperp_1_0, _a_Fperp_1_1, _a_Fperp_1_2;
            UsedParameter _b_Fperp_0_0, _b_Fperp_0_1, _b_Fperp_0_2, _b_Fperp_0_3, _b_Fperp_1_0, _b_Fperp_1_1, _b_Fperp_1_2;
            UsedParameter _c_Fperp_0_0, _c_Fperp_0_1, _c_Fperp_0_2, _c_Fperp_0_3, _c_Fperp_1_0, _c_Fperp_1_1, _c_Fperp_1_2;

            UsedParameter _a_Fpara_0_0, _a_Fpara_0_1, _a_Fpara_0_2, _a_Fpara_0_3, _a_Fpara_1_0, _a_Fpara_1_1, _a_Fpara_1_2;
            UsedParameter _b_Fpara_0_0, _b_Fpara_0_1, _b_Fpara_0_2, _b_Fpara_0_3, _b_Fpara_1_0, _b_Fpara_1_1, _b_Fpara_1_2;
            UsedParameter _c_Fpara_0_0, _c_Fpara_0_1, _c_Fpara_0_2, _c_Fpara_0_3, _c_Fpara_1_0, _c_Fpara_1_1, _c_Fpara_1_2;

            UsedParameter _a_Flong_0_0, _a_Flong_0_1, _a_Flong_0_2, _a_Flong_0_3, _a_Flong_1_0, _a_Flong_1_1, _a_Flong_1_2;
            UsedParameter _b_Flong_0_0, _b_Flong_0_1, _b_Flong_0_2, _b_Flong_0_3, _b_Flong_1_0, _b_Flong_1_1, _b_Flong_1_2;
            UsedParameter _c_Flong_0_0, _c_Flong_0_1, _c_Flong_0_2, _c_Flong_0_3, _c_Flong_1_0, _c_Flong_1_1, _c_Flong_1_2;

            UsedParameter _a_Ftime_0_0, _a_Ftime_0_1, _a_Ftime_0_2, _a_Ftime_0_3, _a_Ftime_1_0, _a_Ftime_1_1, _a_Ftime_1_2;
            UsedParameter _b_Ftime_0_0, _b_Ftime_0_1, _b_Ftime_0_2, _b_Ftime_0_3, _b_Ftime_1_0, _b_Ftime_1_1, _b_Ftime_1_2;
            UsedParameter _c_Ftime_0_0, _c_Ftime_0_1, _c_Ftime_0_2, _c_Ftime_0_3, _c_Ftime_1_0, _c_Ftime_1_1, _c_Ftime_1_2;

            static double _calc_z(const double & t, const double & t_p, const double & t_0);
            inline double _z(const double & t) const;
            inline double _zhat(const double & that) const;
            inline double _blaschke(const double & z, const double & zh) const;
            inline double _blaschke_res_qhat2(const double & z) const;

          public:
            FvDV2018FormFactors(const Parameters & p, const Options &);
            ~FvDV2018FormFactors();

            static FormFactors<PToPP> * make(const Parameters & parameters, const Options & options);

            virtual complex<double> f_perp(const double & q2, const double & k2, const double & ctheta) const override;
            virtual complex<double> f_para(const double & q2, const double & k2, const double & ctheta) const override;
            virtual complex<double> f_long(const double & q2, const double & k2, const double & ctheta) const override;
            virtual complex<double> f_time(const double & q2, const double & k2, const double & ctheta) const override;

            double f_perp_im_res_qhat2(const double & q2, const double & k2) const;
            double f_para_im_res_qhat2(const double & q2, const double & k2) const;
            double f_long_im_res_qhat2(const double & q2, const double & k2) const;
            double f_time_im_res_qhat2(const double & q2, const double & k2) const;

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

    extern template class FvDV2018FormFactors<BToPiPi>;
}

#endif
