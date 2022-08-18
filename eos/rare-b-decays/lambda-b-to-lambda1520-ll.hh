/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Méril Reboud
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA1520_LL_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA1520_LL_HH 1

#include <eos/observable.hh>

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

#include <array>

namespace eos
{
    /*
     * Decay: Lb -> L(1520) l lbar.
     */
    class LambdaBToLambda1520Dilepton :
        public ParameterUser,
        public PrivateImplementationPattern<LambdaBToLambda1520Dilepton>
    {
        public:
            LambdaBToLambda1520Dilepton(const Parameters & parameters, const Options & options);
            ~LambdaBToLambda1520Dilepton();

            struct AngularCoefficients;
            struct Amplitudes;
            class AmplitudeGenerator;
            struct DipoleFormFactors;
            struct FormFactorCorrections;

            /*!
             * @name Signal PDFs
             */
            // @{
            double decay_width(const double & s, const double & c_theta_l, const double & c_theta_Lstar, const double & phi) const;
            // @}

            /*!
             * @name Simple observables (@f$q^2@f$-differential)
             *
             * These observables are differential in q^2, the dilepton
             * invariant mass.
             */
            // @{
            double differential_decay_width(const double & q2) const;
            double differential_branching_ratio(const double & q2) const;
            double differential_forward_backward_asymmetry(const double & q2) const;
            double differential_longitudinal_polarisation(const double & q2) const;
            double differential_transversal_polarisation(const double & q2) const;
            // @}

            /*!
             * @name Angular observables (@f$q^2@f$-differential)
             */
            // @{
            double differential_L_1c(const double & s) const;
            double differential_L_1cc(const double & s) const;
            double differential_L_1ss(const double & s) const;
            double differential_L_2c(const double & s) const;
            double differential_L_2cc(const double & s) const;
            double differential_L_2ss(const double & s) const;
            double differential_L_3ss(const double & s) const;
            double differential_L_4ss(const double & s) const;
            double differential_L_5s(const double & s) const;
            double differential_L_5sc(const double & s) const;
            double differential_L_6s(const double & s) const;
            double differential_L_6sc(const double & s) const;
            // @}

            /*!
             * @name Simple observables (@f$q^2@f$-integrated)
             *
             * These observables are integrated over q^2, the dilepton
             * invariant mass.
             */
            // @{
            double integrated_decay_width(const double & s_min, const double & s_max) const;
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const;
            double integrated_transversal_polarisation(const double & s_min, const double & s_max) const;
            // @}

            /*!
             * @name Angular observables (@f$q^2@f$-integrated)
             */
            // @{
            double integrated_L_1c(const double & s_min, const double & s_max) const;
            double integrated_L_1cc(const double & s_min, const double & s_max) const;
            double integrated_L_1ss(const double & s_min, const double & s_max) const;
            double integrated_L_2c(const double & s_min, const double & s_max) const;
            double integrated_L_2cc(const double & s_min, const double & s_max) const;
            double integrated_L_2ss(const double & s_min, const double & s_max) const;
            double integrated_L_3ss(const double & s_min, const double & s_max) const;
            double integrated_L_4ss(const double & s_min, const double & s_max) const;
            double integrated_L_5s(const double & s_min, const double & s_max) const;
            double integrated_L_5sc(const double & s_min, const double & s_max) const;
            double integrated_L_6s(const double & s_min, const double & s_max) const;
            double integrated_L_6sc(const double & s_min, const double & s_max) const;
            // @}


            /*!
             * Descriptions of the process and its kinematics.
             */
            // @{
            static const std::string description;
            static const std::string kinematics_description_s;
            static const std::string kinematics_description_c_theta_l;
            static const std::string kinematics_description_c_theta_Lstar;
            static const std::string kinematics_description_phi;
            // @}

            /*!
             * Auxiliary methods for unit tests and diagnostic purposes.
             */
            Amplitudes amplitudes(const double & q2) const;

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

    /*!
     * Amplitudes for the decay Lb -> L(1520) l lbar.
     */
    struct LambdaBToLambda1520Dilepton::Amplitudes
    {
        complex<double> b_perp1_right, b_perp1_left;
        complex<double> b_para1_right, b_para1_left;
        complex<double> a_perp1_right, a_perp1_left;
        complex<double> a_para1_right, a_para1_left;
        complex<double> a_perp0_right, a_perp0_left;
        complex<double> a_para0_right, a_para0_left;
    };

    struct LambdaBToLambda1520Dilepton::AngularCoefficients
    {
        double L1c, L1cc, L1ss,
            L2c, L2cc, L2ss,
            L3ss, L4ss,
            L5s, L5sc,
            L6s, L6sc;

        AngularCoefficients()
        {
        }

        AngularCoefficients(const std::array<double, 12> & a) :
            L1c(a[0]),
            L1cc(a[1]),
            L1ss(a[2]),
            L2c(a[3]),
            L2cc(a[4]),
            L2ss(a[5]),
            L3ss(a[6]),
            L4ss(a[7]),
            L5s(a[8]),
            L5sc(a[9]),
            L6s(a[10]),
            L6sc(a[11])
        {
        }
    };
}

#endif
