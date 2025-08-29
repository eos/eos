/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2014 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_NONLOCAL_FORM_FACTORS_CHARM_LOOPS_HH
#define EOS_GUARD_SRC_NONLOCAL_FORM_FACTORS_CHARM_LOOPS_HH 1

#include <eos/maths/complex.hh>
#include <eos/maths/interpolation.hh>
#include <eos/models/model.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/reference-name.hh>

#include <vector>

namespace eos
{
    class CharmLoopsInterpolation
    {
        private:
        const CSplineInterpolation real_part;
        const CSplineInterpolation imag_part;

        public:
        CharmLoopsInterpolation(std::vector<double> x, std::vector<double> y_real, std::vector<double> y_imag):
            real_part(x, y_real),
            imag_part(x, y_imag)
        {}
        ~CharmLoopsInterpolation() = default;

        complex<double> operator()(const double & s) const { return complex<double>(real_part(s), imag_part(s)); };
    };

    struct CharmLoops
    {
        /* One-loop functions */
        // cf. [BFS:2001], Eq. (11), p. 4 with m_q -> 0
        static complex<double> h(const double & mu, const double & s);
        // cf. [BFS:2001], Eq. (11), p. 4 with m_q -> 0
        static complex<double> h(const double & mu, const double & s, const double & m_q);

        /* Two-loop functions */
        // cf. [S:2004], Eq. (29), p. 8
        static complex<double> A(const double & mu, const double & s, const double & m_b);
        // cf. [S:2004], Eq. (30), pp. 8-9
        static complex<double> B(const double & mu, const double & s, const double & m_b);
        // cf. [S:2004], Eq. (31), p. 9
        static complex<double> C(const double & mu, const double & s);

        /* Non-factorizing two loop contributions */
        // massless case, cf. [S:2004], Eq. (16), p. 6
        static complex<double> F17_massless(const double & mu, const double & s, const double & m_b);
        static complex<double> F19_massless(const double & mu, const double & s, const double & m_b);
        static complex<double> F27_massless(const double & mu, const double & s, const double & m_b);
        static complex<double> F29_massless(const double & mu, const double & s, const double & m_b);
        // massless case, cf. [BFS:2001], Eqs. (82)-(83), p. 30
        static complex<double> F87_massless(const double & mu, const double & s, const double & m_b);
        static complex<double> F89_massless(                   const double & s, const double & m_b);

        // massive case, cf. [ABGW:2003], Eq. (7), p. 8
        static complex<double> F17_massive(const double & mu, const double & s, const double & m_b, const double & m_c);
        static complex<double> F19_massive(const double & mu, const double & s, const double & m_b, const double & m_c);
        static complex<double> F27_massive(const double & mu, const double & s, const double & m_b, const double & m_c);
        static complex<double> F29_massive(const double & mu, const double & s, const double & m_b, const double & m_c);
        static complex<double> delta_F29_massive(const double & mu, const double & s, const double & m_c);

        // helper functions for F8j, cf. [BFS:2001], Eqs. (29) and (84), pp. 8 and 30
        static complex<double> B0(const double & s, const double & m_q);
        static complex<double> C0(const double & s, const double & m_q);

        // Contributions proportional to Q_s = Q_b arising from intermediate ccbar pairs at NLO in alpha_s.
        // These contributions are interpolated from reference values from [AGV:2019]
        // In particular mu = m_b = 4.18 and alpha_s(M_Z) = 0.1185
        static complex<double> F17_massive_Qsb(const double & s);
        static complex<double> F19_massive_Qsb(const double & s);
        static complex<double> F27_massive_Qsb(const double & s);
        static complex<double> F29_massive_Qsb(const double & s);
    };

    struct ShortDistanceLowRecoil
    {
        /*!
         * Effective Wilson coefficient c7 in the region of low hadronic recoil.
         *
         * @param s             dilepton invariant mass
         * @param mu            renormalization scale
         * @param alpha_s       strong coupling evaluated at the scale mu
         * @param m_b_PS        PS mass of the bottom quark
         * @param use_nlo       true, if NLO contributions shall be used
         * @param wc            the Wilson coefficients
         *
         * For the calculation, cf. [GP:2004], Eq. (56)
         */
        static complex<double> c7eff(const double & s, const double & mu, const double & alpha_s, const double & m_b_PS, bool use_nlo,
                const WilsonCoefficients<BToS> & wc);

        /*!
         * Effective Wilson coefficient c8 at leading order.
         *
         * @param wc            the Wilson coefficients
         *
         * For the calculation, cf. [BFS:2001], below eq. (26)
         */
        static complex<double> c8eff(const WilsonCoefficients<BToS> & wc);

        /*!
         * Effective Wilson coefficient c9 in the region of low hadronic recoil.
         *
         * @param s                     dilepton invariant mass
         * @param mu                    renormalization scale
         * @param alpha_s               strong coupling evaluated at the scale mu
         * @param m_b_PS                PS mass of the bottom quark
         * @param m_c                   MSbar mass of the charm quark
         * @param use_nlo               true, if NLO contributions shall be used
         * @param ccbar_resonance       true, if phenomenological data from e^+e^- -> ccbar resonance -> hadrons shall be used
         * @param lambda_hat_u          certain combination of CKM matrix elements: V_ub V_us^* / (V_tb V_ts^*)
         * @param wc                    the Wilson coefficients
         *
         * For the calculation, cf. [GP:2004], Eq. (55), p. 10
         */
        static complex<double> c9eff(const double & s, const double & mu, const double & alpha_s, const double & m_b_PS, const double & m_c_MSbar,
                bool use_nlo, bool ccbar_resonance, const complex<double> & lambda_hat_u,
                const WilsonCoefficients<BToS> & wc);
    };

    class CharmLoopsAdapter :
        public ParameterUser,
        public PrivateImplementationPattern<CharmLoopsAdapter>
    {
        public:
            CharmLoopsAdapter(const Parameters & parameters, const Options & options);
            ~CharmLoopsAdapter();

            // Observables
            double real_F17(const double & re_q2, const double & im_q2) const;
            double imag_F17(const double & re_q2, const double & im_q2) const;
            double real_F19(const double & re_q2, const double & im_q2) const;
            double imag_F19(const double & re_q2, const double & im_q2) const;
            double real_F27(const double & re_q2, const double & im_q2) const;
            double imag_F27(const double & re_q2, const double & im_q2) const;
            double real_F29(const double & re_q2, const double & im_q2) const;
            double imag_F29(const double & re_q2, const double & im_q2) const;

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
}

#endif
