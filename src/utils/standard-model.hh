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
    class StandardModel :
        public Model
    {
        private:
            /* CKM Wolfenstein parameters */
            Parameter _A;
            Parameter _lambda;
            Parameter _rhobar;
            Parameter _etabar;

            /* QCD parameters */
            Parameter _alpha_s_Z;
            Parameter _mu_t;
            Parameter _mu_b;
            Parameter _mu_c;
            Parameter _lambda_qcd;

            /* Masses */
            Parameter _m_t_pole;
            Parameter _m_b_MSbar;
            Parameter _m_c_MSbar;
            Parameter _m_W;
            Parameter _m_Z;

            /* GSW parameters */
            Parameter _sw2;

            /* Renormalization scales */
            Parameter _mu;
            Parameter _mu_0c;
            Parameter _mu_0t;

        public:
            StandardModel(const Parameters &);
            virtual ~StandardModel();

            static std::shared_ptr<Model> make(const Parameters &);

            /* QCD */
            virtual double alpha_s(const double &) const;
            virtual double m_t_msbar(const double & mu) const;
            virtual double m_t_pole() const;
            virtual double m_b_msbar(const double & mu) const;
            virtual double m_b_pole() const;
            virtual double m_b_ps(const double & mu_f) const;
            virtual double m_c_msbar(const double & mu) const;
            virtual double m_c_pole() const;

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

            /* b->s Wilson coefficients */
            virtual WilsonCoefficients<BToS> wilson_coefficients_b_to_s() const;
    };
}

#endif
