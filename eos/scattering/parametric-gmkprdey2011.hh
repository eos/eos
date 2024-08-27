/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024-2025 Florian Herren
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

#ifndef EOS_GUARD_EOS_SCATTERING_PARAMETRIC_GMKPRDEY2011_HH
#define EOS_GUARD_EOS_SCATTERING_PARAMETRIC_GMKPRDEY2011_HH 1

#include <eos/scattering/single-channel.hh>
#include <eos/scattering/single-channel-processes.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/omnes-factor.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>
#include <eos/utils/quantum-numbers.hh>
#include <eos/utils/reference-name.hh>

#include <array>

namespace eos
{
    class GMKPRDEY2011ScatteringAmplitudes :
        public ScatteringAmplitudes<PPToPP>
    {
        private:
            // S0-wave has 9 parameters aside from masses
            std::array<UsedParameter, 9> _params_S0;
            // P-wave has 4 parameters aside from masses
            std::array<UsedParameter, 4> _params_P1;
            // D0-wave has 4 parameters aside from masses
            std::array<UsedParameter, 3> _params_D0;

            // Masses
            const UsedParameter _mPi, _mK, _mEta, _mRho, _mF2;
            const UsedParameter _mOmega, _GammaOmega, _kappa;
            // Parameters for conformal mappings
            const UsedParameter _sM_S0, _s0_P1, _s0_D0, _sh_D0;
            // Parameters controlling evolution of phases above \sqrt{s} = 1.42 GeV
            const UsedParameter _cont_pow_S0, _cont_pow_P1, _cont_pow_D0;

            // Omnes factors, we do not have the S0 wave since we need a two-channel treatment
            std::array<double, 4> _intervals_P1;
            std::array<double, 5> _intervals_D0;
            std::function<double(const double &)> _f_phase_P1, _f_phase_D0;
            OmnesFactor<30, 4> _omnes_P1;
            OmnesFactor<40, 5> _omnes_D0;

            QualifiedName _par_name(const std::string & partial_wave, const std::string & par_name, unsigned idx) const;
            QualifiedName _par_name(const std::string & partial_wave, const std::string & par_name) const;

            double _calc_w(const double & s, const double & s0) const;
            double _calc_s(const complex<double> & z, const double & sp, const double & s0) const;
            complex<double> _calc_z(const double & s, const double & sp, const double & s0) const;

            // Scattering phases
            double _phase_S0(const double & s) const;
            double _phase_P1(const double & s) const;
            double _phase_D0(const double & s) const;

        public:
            GMKPRDEY2011ScatteringAmplitudes(const Parameters & p, const Options &);

            ~GMKPRDEY2011ScatteringAmplitudes();

            static ScatteringAmplitudes<PPToPP> * make(const Parameters & parameters, const Options & options);

            virtual complex<double> scattering_amplitude(const double & s, const unsigned & l, const IsospinRepresentation & i) const;
            virtual complex<double> omnes_factor(const double & s, const unsigned & l, const IsospinRepresentation & i) const;
            virtual complex<double> isospin_breaking(const double & s, const unsigned & l, const IsospinRepresentation & i) const;
            virtual complex<double> omnes_outer_function(const double & s, const double & sp, const double & s0, const unsigned & npoints, const unsigned & l, const IsospinRepresentation & i) const;

            Diagnostics diagnostics() const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> options;
    };
}

#endif
