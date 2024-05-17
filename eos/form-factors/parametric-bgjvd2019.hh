/* vim: set sw=4 sts=4 tw=120 et foldmethod=syntax : */

/*
 * Copyright (c) 2018-2023 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGJvD2019_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGJvD2019_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/derivative.hh>
#include <eos/models/model.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/reference-name.hh>
#include <eos/maths/polylog.hh>
#include <eos/maths/power-of.hh>

#include <cmath>
#include <limits>

namespace eos
{
    /* HQET Form Factors, based on [BLPR2017] and [JS2018] */
    template <typename Process_, typename Transition_> class HQETFormFactors;

    class HQETFormFactorBase :
        public virtual ParameterUser
    {
        protected:
            std::shared_ptr<Model> _model;

            // spin avaraged mB mass
            UsedParameter _mBar;

            // parameter for modifying the z function
            UsedParameter _a;

            // option to determine the model for the leading-power IW function
            SwitchOption _opt_lp_model;
            std::function<double (const double &)> _xi;

            // option to determine if we use z^3 terms in the leading-power IW function
            SwitchOption _opt_lp_zorder;
            double _enable_lp_z3;
            double _enable_lp_z4;
            double _enable_lp_z5;

            // option to determine if we use z^2 terms in the subleading-power IW function
            SwitchOption _opt_slp_zorder;
            double _enable_slp_z2;

            // option to determine if we use z^2 terms in the subsubleading-power IW function
            SwitchOption _opt_sslp_zorder;
            double _enable_sslp_z1;
            double _enable_sslp_z2;

            // option to determine if we use the SU3_F-symmetry limit for the subsubleading-power IW functions
            SwitchOption _opt_sslp_limit;

            // parameters for the leading Isgur-Wise function xi
            UsedParameter _xipone, _xippone, _xipppone, _xippppone, _xipppppone;

            // parameters for the subleading Isgur-Wise function chi_2
            UsedParameter _chi2one, _chi2pone, _chi2ppone;

            // parameters for the subleading Isgur-Wise function chi_3
            UsedParameter _chi3pone, _chi3ppone;

            // parameters for the subleading Isgur-Wise function eta
            UsedParameter _etaone, _etapone, _etappone;

            // parameters for subsubleading 1/m_c corrections in h_+ (B->D), equal to delta_{h_+}
            UsedParameter _l1one, _l1pone, _l1ppone;

            // parameters for subsubleading 1/m_c corrections in h_A1 (B->D^*), equal to delta_{A_1}
            UsedParameter _l2one, _l2pone, _l2ppone;

            // parameters for subsubleading 1/m_c corrections
            UsedParameter _l3one, _l3pone, _l3ppone;
            UsedParameter _l4one, _l4pone, _l4ppone;
            UsedParameter _l5one, _l5pone, _l5ppone;
            UsedParameter _l6one, _l6pone, _l6ppone;

            std::string _sslp_prefix(const std::string & prefix);

        public:
            HQETFormFactorBase(const Parameters & p, const Options & o, const std::string & prefix);
            ~HQETFormFactorBase();

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options specifications for the HQETFormFactorBase class.
             */
            static const std::vector<OptionSpecification> option_specifications;

        protected:
            /*
             * HQET parameters following [BLPR2017]
             */
            inline double _mu() const { return 2.31; } // mu^2 = m_b * m_c
            inline double _alpha_s() const { return 0.26; }
            inline double _m_b_1S() const { return 4.71; }
            inline double _m_b_pole() const { return _m_b_1S() * (1 + 2.0 / 9.0 * power_of<2>(_alpha_s())); }
            inline double _m_c_pole() const { return _m_b_pole() - 3.40; }
            inline double _lambda_1() const { return -0.30; }
            inline double _LambdaBar() const { return _mBar - _m_b_pole() + _lambda_1() / (2.0 * _m_b_1S()); }

            /*
             * Interface to Process_-specific kinematics.
             */
            virtual double _w(const double & q2) const = 0;
            virtual double _q2(const double & w) const = 0;

            inline double _zw(const double & w) const
            {
                return (std::sqrt(w + 1.0) - std::sqrt(2.0) * _a()) / (std::sqrt(w + 1.0) + std::sqrt(2.0) * _a());
            }

            inline double _z(const double & q2) const
            {
                const double w = _w(q2);

                return _zw(w);
            }

            /*
             * Isgur-Wise functions
             */
            // uses a power series ansatz
            double _xi_power_series(const double & q2) const;

            // uses an exponential ansatz and expands in (w-1) first, then in z*
            double _xi_exponential(const double & q2) const;

            double _chi2(const double & q2) const;

            double _chi3(const double & q2) const;

            double _eta(const double & q2) const;

            /*
             * Auxilliary functions for the HQET Wilson coefficients
             *
             * We use a fixed scale mu = sqrt(m_b * m_c), with m_b = 4.2 and m_c = 1.27,
             * which yields mu = 2.31 GeV.
             */

            inline double _wz(const double & z) const
            {
                return 0.5 * (z + 1.0 / z);
            }

            inline double _wp(const double & w) const { return w + std::sqrt(w * w - 1.0); }
            inline double _wm(const double & w) const { return w - std::sqrt(w * w - 1.0); }

            inline double _r(const double & w) const
            {
                if (w < 1.0)
                    return std::numeric_limits<double>::quiet_NaN();

                if (w - 1.0 < 1.0e-5)
                    return 1.0 - (w - 1.0) / 3.0;

                return std::log(_wp(w)) / std::sqrt(w * w - 1.0);
            }

            inline double _Omega(const double & w, const double & z) const
            {
                if (w < 1.0)
                    return std::numeric_limits<double>::quiet_NaN();

                const double lnz = std::log(z);

                if (w - 1.0 < 1.0e-5)
                    return -1.0 - (1.0 + z) / (1.0 - z) * lnz;

                const double wm = _wm(w);
                const double wp = _wp(w);

                const complex<double> li2wmz = dilog(1.0 - wm * z);
                const complex<double> li2wpz = dilog(1.0 - wp * z);
                const complex<double> li2wm2 = dilog(1.0 - wm * wm);
                const complex<double> li2wp2 = dilog(1.0 - wp * wp);

                return w * real(2.0 * (li2wmz - li2wpz) + li2wp2 - li2wm2) / (2.0 * std::sqrt(w * w - 1.0))
                    - w * _r(w) * lnz + 1.0;
            }

            /* Power corrections */
            double _l1(const double & w) const;
            double _l2(const double & w) const;
            double _l3(const double & w) const;
            double _l4(const double & w) const;
            double _l5(const double & w) const;
            double _l6(const double & w) const;

            /* Wilson Coefficients */
            double _CS(const double & w, const double & z) const;
            double _CP(const double & w, const double & z) const;
            double _CV1(const double & w, const double & z) const;
            double _CV2(const double & w, const double & z) const;
            double _CV3(const double & w, const double & z) const;
            double _CA1(const double & w, const double & z) const;
            double _CA2(const double & w, const double & z) const;
            double _CA3(const double & w, const double & z) const;
            double _CT1(const double & w, const double & z) const;
            double _CT2(const double & w, const double & z) const;
            double _CT3(const double & w, const double & z) const;
    };

    template <typename Process_> class HQETFormFactors<Process_, PToP> :
        public virtual HQETFormFactorBase,
        public virtual FormFactors<PToP>
    {
        private:
            UsedParameter _m_B;
            UsedParameter _m_P;

            /*
             * Kinematics
             */
            virtual double _w(const double & q2) const override;

            virtual double _q2(const double & w) const override;

            /* HQET form factors h_i */
            double _h_p(const double & q2) const;
            double _h_m(const double & q2) const;
            double _h_S(const double & q2) const;
            double _h_T(const double & q2) const;

        public:
            HQETFormFactors(const Parameters & p, const Options & o);
            ~HQETFormFactors();

            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options);

            virtual double f_p(const double & q2) const override;
            virtual double f_0(const double & q2) const override;
            virtual double f_t(const double & q2) const override;
            virtual double f_plus_T(const double & q2) const override;
            double f_m(const double & q2) const;

            /* HQET form factors h_i */
            inline double h_p(const double & q2) const { return _h_p(q2); };
            inline double h_m(const double & q2) const { return _h_m(q2); };
            inline double h_S(const double & q2) const { return _h_S(q2); };
            inline double h_T(const double & q2) const { return _h_T(q2); };

            Diagnostics diagnostics() const;

            using HQETFormFactorBase::references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };

    extern template class HQETFormFactors<BToD,   PToP>;
    extern template class HQETFormFactors<BsToDs, PToP>;

    template <typename Process_> class HQETFormFactors<Process_, PToV> :
        public virtual HQETFormFactorBase,
        public virtual FormFactors<PToV>
    {
        private:
            UsedParameter _m_B;
            UsedParameter _m_V;

            /*
             * Kinematics
             */
            virtual double _w(const double & q2) const override;

            virtual double _q2(const double & w) const override;

            /* HQET form factors h_i */

            double _h_a1(const double & q2) const;
            double _h_a2(const double & q2) const;
            double _h_a3(const double & q2) const;
            double _h_v(const double & q2) const;
            double _h_t1(const double & q2) const;
            double _h_t2(const double & q2) const;
            double _h_t3(const double & q2) const;

        public:
            HQETFormFactors(const Parameters & p, const Options & o);
            ~HQETFormFactors();

            static FormFactors<PToV> * make(const Parameters & parameters, const Options & options);

            virtual double v(const double & q2) const override;
            virtual double a_0(const double & q2) const override;
            virtual double a_1(const double & q2) const override;
            virtual double a_2(const double & q2) const override;
            virtual double a_12(const double & q2) const override;
            virtual double t_1(const double & q2) const override;
            virtual double t_2(const double & q2) const override;
            virtual double t_3(const double & q2) const override;
            virtual double t_23(const double & q2) const override;

            double a_3(const double & q2) const;

            /* HQET form factors h_i */
            inline double h_a1(const double & q2) const { return _h_a1(q2); };
            inline double h_a2(const double & q2) const { return _h_a2(q2); };
            inline double h_a3(const double & q2) const { return _h_a3(q2); };
            inline double h_v(const double & q2)  const { return _h_v(q2);  };
            inline double h_t1(const double & q2) const { return _h_t1(q2); };
            inline double h_t2(const double & q2) const { return _h_t2(q2); };
            inline double h_t3(const double & q2) const { return _h_t3(q2); };


            virtual double f_perp(const double &) const override;
            virtual double f_para(const double &) const override;
            virtual double f_long(const double &) const override;
            virtual double f_perp_T(const double &) const override;
            virtual double f_para_T(const double &) const override;
            virtual double f_long_T(const double &) const override;

            Diagnostics diagnostics() const;

            using HQETFormFactorBase::references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };

    template <typename Process_> class HQETFormFactors<Process_, VToP> :
        public HQETFormFactorBase,
        public FormFactors<VToP>
    {
        private:
            UsedParameter _m_Bst;
            UsedParameter _m_P;

            /*
             * Kinematics
             */
            virtual double _w(const double & q2) const override;
            virtual double _q2(const double & w) const override;

        public:
            HQETFormFactors(const Parameters & p, const Options & o);
            ~HQETFormFactors();

            static FormFactors<VToP> * make(const Parameters & parameters, const Options & options);

            /* HQET form factors h_i */
            virtual double h_abar_1(const double & q2) const override;
            virtual double h_abar_2(const double & q2) const override;
            virtual double h_abar_3(const double & q2) const override;
            virtual double h_vbar(const double & q2) const override;
            virtual double h_tbar_1(const double & q2) const override;
            virtual double h_tbar_2(const double & q2) const override;
            virtual double h_tbar_3(const double & q2) const override;

            Diagnostics diagnostics() const;

            using HQETFormFactorBase::references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };

    template <typename Process_> class HQETFormFactors<Process_, VToV> :
        public HQETFormFactorBase,
        public FormFactors<VToV>
    {
        private:
            /*
             * Kinematics
             */
            virtual double _w(const double & q2) const override;
            virtual double _q2(const double & w) const override;

        public:
            HQETFormFactors(const Parameters & p, const Options & o);
            ~HQETFormFactors();

            static FormFactors<VToV> * make(const Parameters & parameters, const Options & options);

            /* HQET form factors h_i */

            // vector current
            virtual double h_1(const double & q2) const override;
            virtual double h_2(const double & q2) const override;
            virtual double h_3(const double & q2) const override;
            virtual double h_4(const double & q2) const override;
            virtual double h_5(const double & q2) const override;
            virtual double h_6(const double & q2) const override;

            // axial current
            virtual double h_7(const double & q2) const override;
            virtual double h_8(const double & q2) const override;
            virtual double h_9(const double & q2) const override;
            virtual double h_10(const double & q2) const override;

            Diagnostics diagnostics() const;

            using HQETFormFactorBase::references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };
}

#endif
