/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022-2024 Danny van Dyk
 * Copyright (c) 2022-2024 Philip Lüghausen
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

#ifndef EOS_GUARD_SRC_FORM_FACTORS_ANALYTIC_B_TO_GAMMA_QCDF_HH
#define EOS_GUARD_SRC_FORM_FACTORS_ANALYTIC_B_TO_GAMMA_QCDF_HH 1

#include <eos/form-factors/heavy-meson-lcdas.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/models/model.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>

#include <array>
#include <string>
#include <memory>

namespace eos
{
    template <typename Process_>
    struct AnalyticFormFactorPToGammaQCDFTraits;

    template <>
    struct AnalyticFormFactorPToGammaQCDFTraits<BToGamma>
    {
        std::shared_ptr<HeavyMesonLCDAs> blcdas;
        std::shared_ptr<Model> model;

        static const qnp::Prefix prefix;
        static const qnp::Prefix hadronic_prefix;
        static const qnp::Prefix process;
        static const QualifiedName decay_constant;
        static const QualifiedName mass;

        static const constexpr double e_spectator =  2.0 / 3.0;
        static const constexpr double e_heavy     = -1.0 / 3.0;

        AnalyticFormFactorPToGammaQCDFTraits(const Parameters & p, const Options & o);

        double m_heavy_pole(unsigned int loop_order) const;

    };

    /*!
     * P->gamma form factors, for P = B^-, D^+, D_s^+ a heavy-light pseudoscalar meson
     *
     * We use the results obtained in QCD factorization with subleading
     * power corrections according to Ref. [BBJW:2018A].
     *
     * We further parametrize the leading LCDA phi_+ as described in
     * Ref. [FLvD:2022A] and presently omit higher-twist contributions.
     */
    template <typename Process_>
    class AnalyticFormFactorPToGammaQCDF:
        public FormFactors<PToGamma>
    {
        private:
            using Traits = AnalyticFormFactorPToGammaQCDFTraits<Process_>;

            /*
             * The form factors receive contributions in terms of integral
             * convolutions of the B-meson LCDAs.
             * Using the parametrization for the leading LCDA,
             *
             *      phi_+(w) = sum[k] a_k f_k(w) ,
             *
             * any integral involving phi_+ is expressed as a weighted sum of
             * LCDA coefficients a_k, or, in vector notation, as the inner product
             * of two vectors weights w and coefficients a:
             *
             *      int[w, ...] phi(w) kernel(w) = sum[k] a_k w_k = a.w .
             *
             * We implement the weights as fixed-size arrays.
             */
            static const unsigned int number_of_parameters = 9u;
            using Weights = std::array<double, number_of_parameters>;

            Traits traits;
            std::shared_ptr<Model> model;

            UsedParameter mu;
            UsedParameter omega_0;
            UsedParameter f_B;
            UsedParameter m_B;
            UsedParameter m_rho;
            UsedParameter lambda_bar;
            UsedParameter lambda_E2;
            UsedParameter lambda_H2;
            UsedParameter M2;
            UsedParameter s_0;
            UsedParameter mu_h1;
            UsedParameter mu_h2;

            SwitchOption opt_contributions;
            double switch_ht;
            double switch_soft;
            double switch_soft_tw_3_4;

            static const constexpr double e_spectator = AnalyticFormFactorPToGammaQCDFTraits<Process_>::e_spectator;
            static const constexpr double e_heavy     = AnalyticFormFactorPToGammaQCDFTraits<Process_>::e_heavy;
            static const constexpr double C_F         =  4.0 / 3.0;
            static const constexpr double n_l         =  4.0;

            std::string par_qname(const std::string & _name) const;

            /*!
             * Decomposition of the form factors into three terms, see Ref. [BBJW:2018A]:
             *  1. leading-power contribution from HQET with radiative corrections
             *  and two terms that are power-suppressed in 1/Egamma, 1/m_b:
             *  2. a symmetry-preserving term xi(Egamma) and
             *  3. a symmetry-breaking term delta_xi(Egamma)
             */
            ///@{
            std::tuple<double, double, double> C_K_inv_U(const double & Egamma) const;
            double F_leading_power(const double & Egamma) const;
            double xi(const double & Egamma) const;
            double delta_xi(const double & Egamma) const;
            ///@}

            /*!
             * Leading-order ingredients
             *
             * Functionals of phi_+
             */
            ///@{

            /*!
             * The inverse moment
             * int[w,0,inf] 1/w phi_+(w)
             */
            double L0() const;

            /*!
             * The incomplete inverse moment
             * int[w,0,omega_cut] 1/w phi_+(w)
             */
            double L0_incomplete(const double & omega_cut) const;

            /*!
             * The incomplete normalization
             * int[w,0,omega_cut] phi_+(w)
             */
            double norm_incomplete(const double & omega_cut) const;

            /*!
             * The incomplete Laplace transform
             * int[w,0,omega_cut] exp(-sigma w) phi_+(w)
             */
            double lapltr_incomplete(const double & omega_cut, const double & sigma) const;

            /*!
             * The derivative (-sigma) d/(d sigma) of the incomplete Laplace transform
             * int[w,0,omega_cut] exp(-sigma w) (-sigma w) phi_+(w)
             */
            double lapltr_incomplete_dsigma(const double & omega_cut, const double & sigma) const;
            ///@}


            /*!
             * Next-to-leading order radiative contributions
             *
             * Functionals of Delta phi_+^eff
             * Note: this does *not* include the leading order term
             * and this *omits* the factor alpha_s * C_F / (4.0 * pi)
             */
            ///@{

            /*!
             * The effective inverse moment
             * int[w,0,inf] 1/w Delta phi_+^eff
             */
            double L0_effective(const double & Egamma) const;

            /*!
             * The effective incomplete inverse moment
             * int[w,0,omega_cut] 1/w Delta phi_+^eff(w)
             */
            double L0_incomplete_effective(const double & Egamma, const double & omega_cut) const;

            /*!
             * The effective incomplete Laplace transform
             * int[w,0,omega_cut] exp(-sigma w) Delta phi_+^eff(w)
             */
            double lapltr_incomplete_effective(const double & Egamma, const double & omega_cut, const double & sigma, bool use_approxmiation = true) const;
            ///@}

        public:
            AnalyticFormFactorPToGammaQCDF(const Parameters &, const Options &);

            static FormFactors<PToGamma> * make(const Parameters &, const Options &);

            virtual double F_A(const double & Egamma) const;
            virtual double F_V(const double & Egamma) const;

            /* Diagnostics for unit tests */
            Diagnostics diagnostics() const;
    };

    extern template class AnalyticFormFactorPToGammaQCDF<BToGamma>;
}

#endif
