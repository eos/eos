/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2015 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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
    template <typename Tag_> class WilsonScanComponent;

    template <>
    class WilsonScanComponent<components::DeltaBS1> :
        public virtual ModelComponent<components::DeltaBS1>
    {
        protected:
            /* QCD parameters */
            UsedParameter _alpha_s_Z__deltabs1;
            UsedParameter _mu_b__deltabs1;

            /* Masses */
            UsedParameter _m_Z__deltabs1;

            /* Renormalization scale */
            UsedParameter _mu__deltabs1;

            /* b->s Wilson coefficients */
            UsedParameter _c1;
            UsedParameter _c2;
            UsedParameter _c3;
            UsedParameter _c4;
            UsedParameter _c5;
            UsedParameter _c6;
            UsedParameter _re_c7, _im_c7;
            UsedParameter _c8;
            UsedParameter _re_c9, _im_c9;
            UsedParameter _re_c10, _im_c10;
            UsedParameter _re_c7prime, _im_c7prime;
            UsedParameter _c8prime;
            UsedParameter _re_c9prime, _im_c9prime;
            UsedParameter _re_c10prime, _im_c10prime;
            UsedParameter _re_cS, _im_cS;
            UsedParameter _re_cSprime, _im_cSprime;
            UsedParameter _re_cP, _im_cP;
            UsedParameter _re_cPprime, _im_cPprime;
            UsedParameter _re_cT, _im_cT;
            UsedParameter _re_cT5, _im_cT5;

            std::function<complex<double> ()> _c7;
            std::function<complex<double> ()> _c9;
            std::function<complex<double> ()> _c10;
            std::function<complex<double> ()> _c7prime;
            std::function<complex<double> ()> _c9prime;
            std::function<complex<double> ()> _c10prime;
            std::function<complex<double> ()> _cS;
            std::function<complex<double> ()> _cSprime;
            std::function<complex<double> ()> _cP;
            std::function<complex<double> ()> _cPprime;
            std::function<complex<double> ()> _cT;
            std::function<complex<double> ()> _cT5;

        public:
            WilsonScanComponent(const Parameters &, const Options &, ParameterUser &);

            /*! b->s Wilson coefficients */
            virtual WilsonCoefficients<BToS> wilson_coefficients_b_to_s(const bool & cp_conjugate) const;
    };

    template <>
    class WilsonScanComponent<components::DeltaBU1> :
        public virtual ModelComponent<components::DeltaBU1>
    {
        private:
            /* b->u Wilson coefficients */
            UsedParameter _re_csl;
            UsedParameter _im_csl;
            UsedParameter _re_csr;
            UsedParameter _im_csr;
            UsedParameter _re_cvl;
            UsedParameter _im_cvl;
            UsedParameter _re_cvr;
            UsedParameter _im_cvr;
            UsedParameter _re_ct;
            UsedParameter _im_ct;
            std::function<complex<double> ()> _csl;
            std::function<complex<double> ()> _csr;
            std::function<complex<double> ()> _cvl;
            std::function<complex<double> ()> _cvr;
            std::function<complex<double> ()> _ct;

        public:
            WilsonScanComponent(const Parameters &, const Options &, ParameterUser &);

            /* b->u Wilson coefficients */
            virtual WilsonCoefficients<BToU> wilson_coefficients_b_to_u(const bool & cp_conjugate) const;
    };

    /*!
     * A model with all possible operators; their Wilson coefficients
     * are allowed to have arbitrary values.
     */
    class WilsonScanModel :
        public Model,
        public SMComponent<components::CKM>,
        public SMComponent<components::QCD>,
        public WilsonScanComponent<components::DeltaBS1>,
        public WilsonScanComponent<components::DeltaBU1>
    {
        public:
            WilsonScanModel(const Parameters &, const Options &);
            virtual ~WilsonScanModel();

            static std::shared_ptr<Model> make(const Parameters &, const Options &);
    };

    class ConstrainedWilsonScanComponent :
        public WilsonScanComponent<components::DeltaBS1>
    {
        public:
            ConstrainedWilsonScanComponent(const Parameters &, const Options &, ParameterUser &);
    };

    /*!
     * Special case of @see WilsonScanModel with C_S = - C_P, C'_S = C'_P, and C_T = C_T5 = 0.
     *
     * As shown in arXiv:1407.7044 eq. (8), the Wilson coefficients are not
     * independent if new physics is well above the electro-weak scale,
     * respects the SM gauge symmetry, and only dim. 6 operators contribute.
     */
    class ConstrainedWilsonScanModel :
        public Model,
        public SMComponent<components::CKM>,
        public SMComponent<components::QCD>,
        public ConstrainedWilsonScanComponent,
        public WilsonScanComponent<components::DeltaBU1>
    {
        public:
            ConstrainedWilsonScanModel(const Parameters &, const Options &);
            virtual ~ConstrainedWilsonScanModel();

            static std::shared_ptr<Model> make(const Parameters &, const Options &);
    };
}

#endif
