/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2023 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_MODELS_WILSON_COEFFICIENTS_HH
#define EOS_GUARD_EOS_MODELS_WILSON_COEFFICIENTS_HH 1

#include <eos/maths/complex.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qcd.hh>

#include <array>
#include <cmath>

namespace eos
{
    template <typename Tag_> struct WilsonCoefficients;

    namespace bern
    {
        struct ClassII {};  // |Delta B| = 1 semileptonic operators, cf. [AFGV:2017A], eq. (2.5), p. 6.
    }

    using ChargedCurrent = bern::ClassII;

    template <> struct WilsonCoefficients<ChargedCurrent>
    {
        /*
         * We follow the definition of [FMvD2013], eqs. (1) and (2), p. 2.
         * This coincides with the Bern basis of class II operators in Ref.
         * [AFGV:2017A], eq. (2.5), p. 6, up to a factor of V_qb.
         */

        /* Order: C_V,L, C_V,R, C_S,L, C_S,R, C_T, or equivalently: 1, 1', 5, 5', 7'. */
        std::array<complex<double>, 5> _coefficients;

        // cf. [FMvD2015], eqs. (1) and (2)
        inline complex<double> cvl() const { return _coefficients[0]; }
        inline complex<double> cvr() const { return _coefficients[1]; }
        inline complex<double> csl() const { return _coefficients[2]; }
        inline complex<double> csr() const { return _coefficients[3]; }
        inline complex<double> ct()  const { return _coefficients[4]; }
    };

    struct BToS {};

    template <> struct WilsonCoefficients<BToS>
    {
        /* Order: c1..c6, cq3..cq6, c2b, c7..c10 */
        std::array<complex<double>, 15> _sm_like_coefficients;

        /* Same order as above, with helicity flip */
        std::array<complex<double>, 15> _primed_coefficients;

        /* Scalar, pseudoscalar, and tensor coefficients */
        std::array<complex<double>, 6> _scalar_tensor_coefficients;

        double _alpha_s;

        /*! Default ctor */
        WilsonCoefficients();

        // Misiak basis, cf. [BMU1999], Eq. (2), p. 3
        inline complex<double> c1() const { return _sm_like_coefficients[0]; }
        inline complex<double> c2() const { return _sm_like_coefficients[1]; }
        inline complex<double> c3() const { return _sm_like_coefficients[2]; }
        inline complex<double> c4() const { return _sm_like_coefficients[3]; }
        inline complex<double> c5() const { return _sm_like_coefficients[4]; }
        inline complex<double> c6() const { return _sm_like_coefficients[5]; }

        inline complex<double> cq3() const { return _sm_like_coefficients[6]; }
        inline complex<double> cq4() const { return _sm_like_coefficients[7]; }
        inline complex<double> cq5() const { return _sm_like_coefficients[8]; }
        inline complex<double> cq6() const { return _sm_like_coefficients[9]; }

        inline complex<double> c2b() const { return _sm_like_coefficients[10]; }

        inline complex<double> c7() const { return 4.0 * M_PI / _alpha_s * _sm_like_coefficients[11]; }
        inline complex<double> c8() const { return 4.0 * M_PI / _alpha_s * _sm_like_coefficients[12]; }
        inline complex<double> c9() const { return 4.0 * M_PI / _alpha_s * _sm_like_coefficients[13]; }
        inline complex<double> c10() const { return 4.0 * M_PI / _alpha_s * _sm_like_coefficients[14]; }

        inline complex<double> c7prime() const { return 4.0 * M_PI / _alpha_s * _primed_coefficients[11]; }
        inline complex<double> c8prime() const { return 4.0 * M_PI / _alpha_s * _primed_coefficients[12]; }
        inline complex<double> c9prime() const { return 4.0 * M_PI / _alpha_s * _primed_coefficients[13]; }
        inline complex<double> c10prime() const { return 4.0 * M_PI / _alpha_s * _primed_coefficients[14]; }

        inline complex<double> cS() const { return  _scalar_tensor_coefficients[0]; }
        inline complex<double> cSprime() const { return  _scalar_tensor_coefficients[1]; }
        inline complex<double> cP() const { return  _scalar_tensor_coefficients[2]; }
        inline complex<double> cPprime() const { return  _scalar_tensor_coefficients[3]; }
        inline complex<double> cT() const { return  _scalar_tensor_coefficients[4]; }
        inline complex<double> cT5() const { return  _scalar_tensor_coefficients[5]; }
    };

    /*!
     * Evolution of b -> s Wilson coefficients
     *
     * Calculation according to [BMU1999], Eq. (25). The helicity-flipped, scalar,
     * pseudoscalar, and tensor components are implicitly kept at zero.
     *
     * @param wc_qcd_0  The initial scale Wilson coefficients at O(alpha_s^0)
     * @param wc_qcd_1  The initial scale Wilson coefficients at O(alpha_s^1)
     * @param wc_qcd_2  The initial scale Wilson coefficients at O(alpha_s^2)
     * @param alpha_s_0 The strong coupling constant at the initial scale
     * @param alpha_s   The strong coupling constant at the low scale
     * @param nf        The number of active flavors
     * @param beta      Coefficients of the beta function of QCD for nf active flavors.
     */
    WilsonCoefficients<BToS> evolve(const std::array<complex<double>, 15> & wc_qcd_0,
            const std::array<complex<double>, 15> & wc_qcd_1,
            const std::array<complex<double>, 15> & wc_qcd_2,
            const double & alpha_s_0, const double & alpha_s,
            const double & nf, const QCD::BetaFunction & beta);

    namespace wc
    {
        struct SBSB {};
        struct SBNuNu {};
    }

    /* Wilson coefficients for |Delta B| = |Delta S| = 2 operators */
    template <> struct WilsonCoefficients<wc::SBSB>
    {
        /*
         * The eight operators are defined as in [AFGV:2017A], eq. (2.4), p. 6
         *
         * C++ idx -> operator
         * 0       -> O_1
         * 1       -> O_2
         * 2       -> O_3
         * 3       -> O_4
         * 4       -> O_5
         * 5       -> O_1'
         * 6       -> O_2'
         * 7       -> O_3'
         */
        std::array<complex<double>, 8> _coefficients;

        /*! Default ctor */
        WilsonCoefficients();

        // AFGV basis
        inline complex<double> c1()   const { return _coefficients[0]; }
        inline complex<double> c2()   const { return _coefficients[1]; }
        inline complex<double> c3()   const { return _coefficients[2]; }
        inline complex<double> c4()   const { return _coefficients[3]; }
        inline complex<double> c5()   const { return _coefficients[4]; }
        inline complex<double> c1p()  const { return _coefficients[5]; }
        inline complex<double> c2p()  const { return _coefficients[6]; }
        inline complex<double> c3p()  const { return _coefficients[7]; }
    };

    template <> struct WilsonCoefficients<wc::SBNuNu>
    {
        /*
         * Assuming only left-handed neutrinos, only two vector operators exist beyond the SM.
         * These two operators are defined via
         *
         * H^eff_vec = 4 G_F / sqrt(2) V_tb V_ts^* alpha_e / (2 pi) (
         *             C_VL [sbar gamma^mu P_L b][nubar gamma_mu P_L nu]
         *             C_VR [sbar gamma^mu P_R b][nubar gamma_mu P_L nu]
         *         )
         *       = 4 G_F / sqrt(2) V_tb V_ts^* alpha_e / (2 pi) ( C_VL O_VL + C_VR O_VR )
         *
         * O_VR extends the SM to include coupling to RH fermions.
         *
         * Furthermore, we can add the scalar and tensor operators according to [FLS:2021A], eq. (2).
         * These operators go beyond the SM, by including Majorana neutrinos.
         *
         * The scalar operator contributions are
         *
         * H^eff_sca = 4 G_F / sqrt(2) V_tb V_ts^* alpha_e / (2 pi) ( C_SL O_SL + C_SR O_SR )
         *
         * where
         * O_SL = [sbar P_L b][bar(nu_L^c) P_L nu], where nu_L^c = C nubar_L^T, C = i gamma^2 gamma^0
         * O_SR = [sbar P_R b][bar(nu_L^c) P_L nu], where nu_L^c = C nubar_L^T, C = i gamma^2 gamma^0
         *
         * The tensor operator contribution is
         *
         * H^eff_ten = 4 G_F / sqrt(2) V_tb V_ts^* alpha_e / (2 pi) C_TL O_TL
         *
         * where
         * O_TL = [sbar sigma_{mu nu} P_L b][bar(nu_L^c) sigma^{mu nu} P_L nu], where nu_L^c = C nubar_L^T, C = i gamma^2 gamma^0
         *
         *
         * In total we then have
         * H^eff = H^eff_vec + H^eff_sca + H^eff_ten
         *
         *
         * C++ idx -> operator
         * 0       -> O_VL
         * 1       -> O_VR
         * 2       -> O_SL
         * 3       -> O_SR
         * 4       -> O_TL
         */
        std::array<complex<double>, 5> _coefficients;

        /*! Default ctor */
        WilsonCoefficients();

        inline complex<double> cVL()   const { return _coefficients[0]; }
        inline complex<double> cVR()   const { return _coefficients[1]; }
        inline complex<double> cSL()   const { return _coefficients[2]; }
        inline complex<double> cSR()   const { return _coefficients[3]; }
        inline complex<double> cTL()   const { return _coefficients[4]; }
    };}

#endif
