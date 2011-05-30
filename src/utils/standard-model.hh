/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_STANDARD_MODEL_HH
#define EOS_GUARD_SRC_UTILS_STANDARD_MODEL_HH 1

#include <src/utils/model.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace eos
{
    template <typename Tag> class SMComponent;

    template <> class SMComponent<components::CKM> :
        public virtual ModelComponent<components::CKM>
    {
        private:
            /* CKM Wolfenstein parameters */
            UsedParameter _A__ckm;
            UsedParameter _lambda__ckm;
            UsedParameter _rhobar__ckm;
            UsedParameter _etabar__ckm;

        public:
            SMComponent(const Parameters &, ParameterUser &);

            /* CKM matrix elements */
            virtual complex<double> ckm_cd() const;
            virtual complex<double> ckm_cs() const;
            virtual complex<double> ckm_cb() const;
            virtual complex<double> ckm_ud() const;
            virtual complex<double> ckm_us() const;
            virtual complex<double> ckm_ub() const;
            virtual complex<double> ckm_td() const;
            virtual complex<double> ckm_ts() const;
            virtual complex<double> ckm_tb() const;
    };

    template <> class SMComponent<components::QCD> :
        public virtual ModelComponent<components::QCD>
    {
        private:
            /* QCD parameters */
            UsedParameter _alpha_s_Z__qcd;
            UsedParameter _mu_t__qcd;
            UsedParameter _mu_b__qcd;
            UsedParameter _mu_c__qcd;
            UsedParameter _lambda_qcd__qcd;

            /* Masses */
            UsedParameter _m_t_pole__qcd;
            UsedParameter _m_b_MSbar__qcd;
            UsedParameter _m_c_MSbar__qcd;
            UsedParameter _m_Z__qcd;

        public:
            SMComponent(const Parameters &, ParameterUser &);

            /* QCD */
            virtual double alpha_s(const double &) const;
            virtual double m_t_msbar(const double & mu) const;
            virtual double m_t_pole() const;
            virtual double m_b_msbar(const double & mu) const;
            virtual double m_b_pole() const;
            virtual double m_b_ps(const double & mu_f) const;
            virtual double m_c_msbar(const double & mu) const;
            virtual double m_c_pole() const;
    };

    template <> class SMComponent<components::DeltaB1> :
        public virtual ModelComponent<components::DeltaB1>
    {
        private:
            /* QCD parameters */
            UsedParameter _alpha_s_Z__deltab1;
            UsedParameter _mu_t__deltab1;
            UsedParameter _mu_b__deltab1;
            UsedParameter _mu_c__deltab1;

            /* GSW parameters */
            UsedParameter _sw2__deltab1;

            /* Masses */
            UsedParameter _m_t_pole__deltab1;
            UsedParameter _m_W__deltab1;
            UsedParameter _m_Z__deltab1;

            /* Matching scales */
            UsedParameter _mu_0c__deltab1;
            UsedParameter _mu_0t__deltab1;

            /* Renormalization scale */
            UsedParameter _mu__deltab1;

        public:
            SMComponent(const Parameters &, ParameterUser &);

            /* b->s Wilson coefficients */
            virtual WilsonCoefficients<BToS> wilson_coefficients_b_to_s(const bool & cp_conjugate) const;
    };

    class StandardModel :
        public Model,
        public SMComponent<components::CKM>,
        public SMComponent<components::QCD>,
        public SMComponent<components::DeltaB1>
    {
        public:
            StandardModel(const Parameters &);
            virtual ~StandardModel();

            static std::shared_ptr<Model> make(const Parameters &, const Options &);
    };
}

#endif
