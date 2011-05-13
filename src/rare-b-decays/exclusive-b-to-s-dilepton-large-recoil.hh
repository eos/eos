/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_LARGE_RECOIL_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_LARGE_RECOIL_HH 1

#include <src/rare-b-decays/decays.hh>
#include <src/utils/complex.hh>
#include <src/utils/options.hh>
#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace eos
{
    struct LargeRecoil
    {
    };

    /*
     * Decay: B -> K^* l lbar at Large Recoil, cf. [BHP2008]
     */
    template <>
    class BToKstarDilepton<LargeRecoil> :
        public ParameterUser,
        public PrivateImplementationPattern<BToKstarDilepton<LargeRecoil>>
    {
        public:
            BToKstarDilepton(const Parameters & parameters, const Options & options);
            ~BToKstarDilepton();

            // [BHP2008], Appendix C
            complex<double> a_long(const Helicity & h, const double & s) const;
            complex<double> a_perp(const Helicity & h, const double & s) const;
            complex<double> a_par(const Helicity & h, const double & s) const;

            // Differential Observables
            double differential_branching_ratio(const double & s) const;
            double differential_decay_width(const double & s) const;
            double differential_forward_backward_asymmetry(const double & s) const;
            double differential_longitudinal_polarisation(const double & s) const;
            double differential_transverse_asymmetry_2(const double & s) const;
            double differential_transverse_asymmetry_3(const double & s) const;
            double differential_transverse_asymmetry_4(const double & s) const;
            double differential_transverse_asymmetry_5(const double & s) const;

            // Integrated Observables
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_unnormalized_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const;
            double integrated_unnormalized_longitudinal_polarisation(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_5(const double & s_min, const double & s_max) const;
    };

    /*
     * Decay: B -> K l lbar at Large Recoil, cf. [BFS2001], [BHP2007]
     */
    template <>
    class BToKDilepton<LargeRecoil> :
        public ParameterUser,
        public PrivateImplementationPattern<BToKDilepton<LargeRecoil>>
    {
        public:
            BToKDilepton(const Parameters & parameters, const Options & options);
            ~BToKDilepton();

            double a_l(const double & s) const;
            double c_l(const double & s) const;

            // Differential Observables
            double differential_branching_ratio(const double & s) const;
            double differential_flat_term(const double & s) const;

            // Integrated Observables
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_flat_term(const double & s_min, const double & s_max) const;
    };
}

#endif
