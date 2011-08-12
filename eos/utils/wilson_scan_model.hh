/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_WILSON_SCAN_MODEL_HH
#define EOS_GUARD_SRC_UTILS_WILSON_SCAN_MODEL_HH 1

#include <eos/utils/model.hh>
#include <eos/utils/standard-model.hh>

namespace eos
{
    class WilsonScanComponent :
        public virtual ModelComponent<components::DeltaB1>
    {
        private:
            /* QCD parameters */
            UsedParameter _alpha_s_Z__deltab1;
            UsedParameter _mu_b__deltab1;

            /* Masses */
            UsedParameter _m_Z__deltab1;

            /* Renormalization scale */
            UsedParameter _mu__deltab1;

            /* b->s Wilson coefficients */
            UsedParameter _c1;
            UsedParameter _c2;
            UsedParameter _c3;
            UsedParameter _c4;
            UsedParameter _c5;
            UsedParameter _c6;
            Parameter _abs_c7, _arg_c7;
            Parameter _re_c7, _im_c7;
            std::function<complex<double> ()> _c7;
            UsedParameter _c8;
            Parameter _abs_c9, _arg_c9;
            Parameter _re_c9, _im_c9;
            std::function<complex<double> ()> _c9;
            Parameter _abs_c10, _arg_c10;
            Parameter _re_c10, _im_c10;
            std::function<complex<double> ()> _c10;
            Parameter _abs_c7prime, _arg_c7prime;
            Parameter _re_c7prime, _im_c7prime;
            std::function<complex<double> ()> _c7prime;
            Parameter _abs_c9prime, _arg_c9prime;
            Parameter _re_c9prime, _im_c9prime;
            std::function<complex<double> ()> _c9prime;
            Parameter _abs_c10prime, _arg_c10prime;
            Parameter _re_c10prime, _im_c10prime;
            std::function<complex<double> ()> _c10prime;

        public:
            WilsonScanComponent(const Parameters &, const Options &, ParameterUser &);

            /* b->s Wilson coefficients */
            virtual WilsonCoefficients<BToS> wilson_coefficients_b_to_s(const bool & cp_conjugate) const;
    };

    class WilsonScanModel :
        public Model,
        public SMComponent<components::CKM>,
        public SMComponent<components::QCD>,
        public WilsonScanComponent
    {
        public:
            WilsonScanModel(const Parameters &, const Options &);
            virtual ~WilsonScanModel();

            static std::shared_ptr<Model> make(const Parameters &, const Options &);
    };
}

#endif
