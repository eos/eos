/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_HH
#define WFITTER_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_HH 1

#include <src/rare-b-decays/decays.hh>
#include <src/utils/complex.hh>
#include <src/utils/observable.hh>
#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace wf
{
    /*
     * Decay: B -> K l lbar
     */


    // Large Recoil, cf. [BHP2008]
    struct LargeRecoil
    {
    };

    template <>
    class BToKstarDilepton<LargeRecoil> :
        public PrivateImplementationPattern<BToKstarDilepton<LargeRecoil>>
    {
        public:
            BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options);
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
            double integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_5(const double & s_min, const double & s_max) const;
    };

    // Low Recoil, cf. [BHvD2010]
    struct LowRecoil
    {
    };

    template <>
    class BToKstarDilepton<LowRecoil> :
        public PrivateImplementationPattern<BToKstarDilepton<LowRecoil>>
    {
        public:
            BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options);
            ~BToKstarDilepton();

            // [BHvD2010] Eqs. (??-??)
            complex<double> a_long(const Helicity & h, const double & s) const;
            complex<double> a_perp(const Helicity & h, const double & s) const;
            complex<double> a_par(const Helicity & h, const double & s) const;

            // [BHvD2010-2] Eqs. (??)
            double real_c(const double & s) const;
            double imag_c(const double & s) const;
            double real_c9eff(const double & s) const;
            double imag_c9eff(const double & s) const;

            // [BHvD2010] Eqs. (??-??)
            double rho_1(const double & s) const;
            double rho_2(const double & s) const;

            // Differential Observables
            double differential_branching_ratio(const double & s) const;
            double differential_decay_width(const double & s) const;
            double differential_forward_backward_asymmetry(const double & s) const;
            double differential_longitudinal_polarisation(const double & s) const;
            double differential_transverse_asymmetry_2(const double & s) const;
            double differential_transverse_asymmetry_3(const double & s) const;
            double differential_transverse_asymmetry_4(const double & s) const;
            double differential_h_1(const double & s) const;
            double differential_h_2(const double & s) const;
            double differential_h_3(const double & s) const;

            // Integrated Observables
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry_naive(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_longitudinal_polarisation_naive(const double & s_min, const double & s_max) const;
            double integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_2_naive(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_3_naive(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const;
            double integrated_transverse_asymmetry_4_naive(const double & s_min, const double & s_max) const;
            double integrated_h_1(const double & s_min, const double & s_max) const;
            double integrated_h_1_naive(const double & s_min, const double & s_max) const;
            double integrated_h_2(const double & s_min, const double & s_max) const;
            double integrated_h_2_naive(const double & s_min, const double & s_max) const;
            double integrated_h_3(const double & s_min, const double & s_max) const;
            double integrated_h_3_naive(const double & s_min, const double & s_max) const;
    };
}

#endif
