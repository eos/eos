/*
 * Copyright (c) 2023 Méril Reboud
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

#ifndef EOS_GUARD_EOS_SCATTERING_EETOCCBAR_HH
#define EOS_GUARD_EOS_SCATTERING_EETOCCBAR_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/kmatrix.hh>
#include <eos/utils/concrete-cacheable-observable.hh>

#include <complex>
#include <vector>

namespace eos

{
    /*
    Channels follow the following convention
    #   name          type         Nf      copy
    0   ee            ee (S)       3       -
    1   effJpsi       PP (P)       3       -
    2   eff(2S)       PP (P)       3       -
    3   D0   D0bar    PP (P)       3       -
    4   D+   D-       PP (P)       3       3 (isospin)
    5   eff(3770)     PP (P)       3       -
    6   eff(4040)     PP (P)       3       -
    7   Ds+  Ds-      PP (P)       3       -
    8   D*0  D*0bar   PP (P)       3       -
    9   D*+  D*-      PP (P)       3       8 (isospin)
    10  D0   D*0bar   PV (S)       3       -
    11  D+   D*-      PV (S)       3       10 (isospin)
    */

    // e^+e^- channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct EEChannel :
    public KMatrix<nchannels_, nresonances_>::Channel
    {
        EEChannel(std::string name, Parameter m1, Parameter m2, Parameter q0, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, 0, q0, g0s)
        {
            if (m1 != m2)
            {
                InternalError("K-matrix channels with different masses are not yet implemented.");
            }
        };

        inline double mp() { return this->_m1 + this->_m2; }
        using KMatrix<nchannels_, nresonances_>::Channel::_q0;

        const double pi = M_PI;
        const complex<double> i = complex<double>(0.0, 1.0);

        complex<double> rho(const complex<double> & s)
        {
            const double mp = this->mp();

            return (real(s) < mp * mp) ? 0.0 : std::sqrt((s - mp * mp) * s) / 16.0 / pi / s;
        }

        // Analytic continuation of i * rho * n * n
        complex<double> chew_mandelstam(const complex<double> & S)
        {
            const double mp = this->mp();
            // Adapt s to match Mathematica's behaviour on the branch cut
            const complex<double> s = S + complex<double>(0.0, 1e-15);

            return -1.0 / 8.0 / pi / pi * std::sqrt(mp * mp - s) *
                std::atan(s / std::sqrt(s * (mp * mp - s))) / std::sqrt(s);
        }
    };

    // Effective channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct EffChannel :
    public KMatrix<nchannels_, nresonances_>::Channel
    {

        EffChannel(std::string name, Parameter m1, Parameter m2, Parameter q0, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, 1, q0, g0s)
        {
            if (m1 != m2)
            {
                InternalError("K-matrix channels with different masses are not yet implemented.");
            }
        };

        inline double mp() { return this->_m1 + this->_m2; }
        using KMatrix<nchannels_, nresonances_>::Channel::_q0;

        const double pi = M_PI;
        const complex<double> i = complex<double>(0.0, 1.0);

        complex<double> rho(const complex<double> & s)
        {
            const double mp = this->mp();

            return (real(s) < mp * mp) ? 0.0 : std::sqrt((s - mp * mp) * s) / 16.0 / pi / s;
        }

        // Analytic continuation of i * rho * n * n
        complex<double> chew_mandelstam(const complex<double> & S)
        {
            const double mp = this->mp();
            const double q0 = this->_q0();
            // Adapt s to match Mathematica's behaviour on the branch cut
            const complex<double> s = S + complex<double>(0.0, 1e-15);
            const complex<double> delta = mp * mp - 4.0 * q0 * q0;

            // Blatt-Weisskopf factors, cf eq. (50.26)
            const complex<double> Fsq = power_of<2>(kmatrix_utils::blatt_weisskopf_factor(1, std::sqrt(s - mp * mp) / 2.0 / q0));

            complex<double> leading_term;
            // Fix the behavior near threshold by Taylor expanding to second order
            if (std::abs(s - mp * mp) < 1e-7)
            {
                leading_term = Fsq * (mp * mp - s) / 16.0 / mp / mp / pi / pi *
                    (-2.0 * (mp * mp - s) + mp * pi * std::sqrt(mp * mp - s));
            }
            else
            {
                leading_term  = Fsq * power_of<3>(std::sqrt(mp * mp - s)) *
                    std::atan(s / std::sqrt(s * (mp * mp - s))) / 8.0 / pi / pi / std::sqrt(s);
            }

            // atan(sqrt(delta) / (2 q0)) / sqrt(delta) is finite as delta -> 0 (limit 1 / (2 q0));
            // take the limit explicitly to avoid the 0/0 NaN at the pseudothreshold delta = 0
            // (reached e.g. by the bare defaults effective_mass == q0).
            const complex<double> atan_over_sqrt_delta = (std::abs(delta) < 1e-12)
                ? complex<double>(1.0 / (2.0 * q0))
                : std::atan(std::sqrt(delta) / 2.0 / q0) / std::sqrt(delta);

            const complex<double> loop_correction = -power_of<3>(q0) * (mp * mp - s) *
                    atan_over_sqrt_delta / pi / pi / (s - delta);

            return (leading_term + loop_correction) / 4.0 / q0 / q0;
        }
    };


    // V -> PP channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct PWavePPChannel :
    public KMatrix<nchannels_, nresonances_>::Channel
    {

        PWavePPChannel(std::string name, Parameter m1, Parameter m2, Parameter q0, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, 1, q0, g0s)
        {
        };

        inline double mp() { return this->_m1 + this->_m2; }
        using KMatrix<nchannels_, nresonances_>::Channel::_q0;

        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        complex<double> rho(const complex<double> & s)
        {
            const double mp = this->mp();

            return (real(s) < mp * mp) ? 0.0 : std::sqrt((s - mp * mp) * s) / 16.0 / pi / s;
        }

        // Analytic continuation of i * rho * n * n
        complex<double> chew_mandelstam(const complex<double> & S)
        {
            const double mp = this->mp();
            const double q0 = this->_q0();
            // Adapt s to match Mathematica's behaviour on the branch cut
            const complex<double> s = S + complex<double>(0.0, 1e-15);
            const complex<double> delta = mp * mp - 4.0 * q0 * q0;

            // Blatt-Weisskopf factors, cf eq. (50.26)
            const complex<double> Fsq = power_of<2>(kmatrix_utils::blatt_weisskopf_factor(1, std::sqrt(s - mp * mp) / 2.0 / q0));

            complex<double> leading_term;
            // Fix the behavior near threshold by Taylor expanding to second order
            if (std::abs(s - mp * mp) < 1e-7)
            {
                leading_term = Fsq * (mp * mp - s) / 16.0 / mp / mp / pi / pi *
                    (-2.0 * (mp * mp - s) + mp * pi * std::sqrt(mp * mp - s));
            }
            else
            {
                leading_term  = Fsq * power_of<3>(std::sqrt(mp * mp - s)) *
                    std::atan(s / std::sqrt(s * (mp * mp - s))) / 8.0 / pi / pi / std::sqrt(s);
            }

            // atan(sqrt(delta) / (2 q0)) / sqrt(delta) is finite as delta -> 0 (limit 1 / (2 q0));
            // take the limit explicitly to avoid the 0/0 NaN at the pseudothreshold delta = 0.
            const complex<double> atan_over_sqrt_delta = (std::abs(delta) < 1e-12)
                ? complex<double>(1.0 / (2.0 * q0))
                : std::atan(std::sqrt(delta) / 2.0 / q0) / std::sqrt(delta);

            const complex<double> loop_correction = -power_of<3>(q0) * (mp * mp - s) *
                    atan_over_sqrt_delta / pi / pi / (s - delta);

            return (leading_term + loop_correction) / 4.0 / q0 / q0;
        }
    };


    // V -> PV channel for two unequal masses (e.g. D Dbar^*), S-wave (l = 0).
    // Uses the full Kaellen phase space and the unequal-mass l=0 Chew-Mandelstam
    // function. The closed form was validated against the equal-mass limit
    // (reduces to EEChannel) and a numerical once-subtracted dispersion integral
    // of rho; see PHASE2-HARD-SCOPE.md and analytic/chew-mandelstam.py. l = 0 =>
    // the barrier factor n = 1, so the Chew-Mandelstam is the analytic
    // continuation of i * rho with no n factor. It is threshold-subtracted so that
    // CM((m1+m2)^2) = 0, matching the convention of EEChannel (l=0) and
    // PWavePPChannel (l=1); the subtracted constant is the real number
    // -(m1^2-m2^2) ln(m1/m2) / [16 pi^2 (m1+m2)^2] (it vanishes for equal masses).
    template <unsigned nchannels_, unsigned nresonances_>
    struct SWavePVChannel :
    public KMatrix<nchannels_, nresonances_>::Channel
    {
        SWavePVChannel(std::string name, Parameter m1, Parameter m2, Parameter q0, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, 0, q0, g0s)
        {
        };

        const double pi = M_PI;
        const complex<double> i = complex<double>(0.0, 1.0);

        // Kaellen function lambda(s, m1^2, m2^2) = (s - (m1+m2)^2)(s - (m1-m2)^2),
        // with a threshold at (m1+m2)^2 and a pseudothreshold at (m1-m2)^2.
        complex<double> kallen(const complex<double> & s)
        {
            const double m1 = this->_m1();
            const double m2 = this->_m2();

            return (s - power_of<2>(m1 + m2)) * (s - power_of<2>(m1 - m2));
        }

        complex<double> rho(const complex<double> & s)
        {
            const double mthr = this->_m1() + this->_m2();

            // Phase space vanishes below threshold (same convention as the
            // equal-mass channels).
            return (real(s) < mthr * mthr) ? complex<double>(0.0, 0.0) : std::sqrt(kallen(s)) / 16.0 / pi / s;
        }

        // Analytic continuation of i * rho for unequal masses, l = 0.
        complex<double> chew_mandelstam(const complex<double> & S)
        {
            const double m1 = this->_m1();
            const double m2 = this->_m2();
            // Adapt s to match the branch-cut prescription of the other channels.
            const complex<double> s = S + complex<double>(0.0, 1e-15);
            const complex<double> w = std::sqrt(kallen(s));

            const complex<double> cm = 1.0 / 16.0 / pi / pi * (
                    w / s * std::log((m1 * m1 + m2 * m2 - s + w) / (2.0 * m1 * m2))
                    - (m1 * m1 - m2 * m2) / s * std::log(m1 / m2)
                    );

            // Threshold subtraction: at s = s_th = (m1+m2)^2 the Kaellen factor
            // (hence w) vanishes, so only the asymmetry term survives. Subtracting
            // this real constant enforces CM(s_th) = 0, consistent with the other
            // channels (it cancels identically for equal masses).
            const double sth = power_of<2>(m1 + m2);
            const double cm_threshold = -1.0 / 16.0 / pi / pi * (m1 * m1 - m2 * m2) / sth * std::log(m1 / m2);

            return cm - cm_threshold;
        }
    };


    template <unsigned nchannels_, unsigned nresonances_>
    struct CharmoniumResonance :
    public KMatrix<nchannels_, nresonances_>::Resonance
    {
        CharmoniumResonance(std::string name, Parameter m) :
        KMatrix<nchannels_, nresonances_>::Resonance(name, m)
        {
        };
    };

    class EEToCCBar :
        public ParameterUser,
        public PrivateImplementationPattern<EEToCCBar>
    {
        public:

            const static long unsigned nchannels = 12;
            const static long unsigned nresonances = 4;

            struct IntermediateResult :
                public CacheableObservable::IntermediateResult
            {
                std::shared_ptr<KMatrix<nchannels, nresonances>> K;
                // Amplitude on the first RS
                std::array<complex<double>, nchannels> tmatrix_row_0;
                // Amplitude on the second RS
                std::array<complex<double>, nchannels> tmatrix2_row_0;

                complex<double> E;
                complex<double> s;
            };

            EEToCCBar(const Parameters & parameters, const Options & options);
            ~EEToCCBar();

            // Observables
            const IntermediateResult * prepare(const double & E) const;
            const IntermediateResult * prepare_complex(const double & reE, const double & imE) const;

            // double evaluate(const IntermediateResult *) const;

            // resonances widths
            double Jpsi_ee_width() const;
            double Jpsi_eff_width() const;
            double Jpsi_total_width() const;
            double psi2S_ee_width() const;
            double psi2S_eff_width() const;
            double psi2S_total_width() const;
            double psi3770_D0Dbar0_width() const;
            double psi3770_DpDm_width() const;
            double psi3770_eff_width() const;
            double psi3770_total_width() const;
            double psi4040_eff_width() const;
            double psi4040_total_width() const;

            // Phase space factor
            double rho_ee(const IntermediateResult *) const;
            double rho_eff(const IntermediateResult *) const;
            double rho_D0Dbar0(const IntermediateResult *) const;
            double rho_DpDm(const IntermediateResult *) const;

            // Chew-Mandelstam function on the first Riemann sheet
            double re_chew_mandelstam_ee(const IntermediateResult *) const;
            double im_chew_mandelstam_ee(const IntermediateResult *) const;
            double re_chew_mandelstam_eff(const IntermediateResult *) const;
            double im_chew_mandelstam_eff(const IntermediateResult *) const;
            double re_chew_mandelstam_D0Dbar0(const IntermediateResult *) const;
            double im_chew_mandelstam_D0Dbar0(const IntermediateResult *) const;
            double re_chew_mandelstam_DpDm(const IntermediateResult *) const;
            double im_chew_mandelstam_DpDm(const IntermediateResult *) const;

            // Chew-Mandelstam function on the second Riemann sheet
            double re_chew_mandelstam_II_ee(const IntermediateResult *) const;
            double im_chew_mandelstam_II_ee(const IntermediateResult *) const;
            double re_chew_mandelstam_II_eff(const IntermediateResult *) const;
            double im_chew_mandelstam_II_eff(const IntermediateResult *) const;
            double re_chew_mandelstam_II_D0Dbar0(const IntermediateResult *) const;
            double im_chew_mandelstam_II_D0Dbar0(const IntermediateResult *) const;
            double re_chew_mandelstam_II_DpDm(const IntermediateResult *) const;
            double im_chew_mandelstam_II_DpDm(const IntermediateResult *) const;


            // amplitudes on the first RS
            double re_T_eetoee(const IntermediateResult *) const;
            double im_T_eetoee(const IntermediateResult *) const;
            double re_T_eetoeff(const IntermediateResult *) const;
            double im_T_eetoeff(const IntermediateResult *) const;
            double re_T_eetoDpDm(const IntermediateResult *) const;
            double im_T_eetoDpDm(const IntermediateResult *) const;
            double re_T_eetoD0Dbar0(const IntermediateResult *) const;
            double im_T_eetoD0Dbar0(const IntermediateResult *) const;

            // amplitudes on the second RS
            double re_T_II_eetoee(const IntermediateResult *) const;
            double im_T_II_eetoee(const IntermediateResult *) const;
            double re_T_II_eetoeff(const IntermediateResult *) const;
            double im_T_II_eetoeff(const IntermediateResult *) const;
            double re_T_II_eetoDpDm(const IntermediateResult *) const;
            double im_T_II_eetoDpDm(const IntermediateResult *) const;
            double re_T_II_eetoD0Dbar0(const IntermediateResult *) const;
            double im_T_II_eetoD0Dbar0(const IntermediateResult *) const;

            // Spectral function
            double psi3770_spectral_function(const double & E) const;

            // sigma(ee -> channel)
            double sigma_eetoee(const IntermediateResult *) const;
            double sigma_eetoeff(const IntermediateResult *) const;
            double sigma_eetoD0Dbar0(const IntermediateResult *) const;
            double sigma_eetoDpDm(const IntermediateResult *) const;
            double sigma_eetoDsDs(const IntermediateResult *) const;
            double sigma_eetoDstar0Dstarbar0(const IntermediateResult *) const;
            double sigma_eetoDstarpDstarm(const IntermediateResult *) const;
            double sigma_eetoD0Dbarstar0(const IntermediateResult *) const;
            double sigma_eetoDpDstarm(const IntermediateResult *) const;

            // R ratios
            double R(const IntermediateResult *) const;

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
