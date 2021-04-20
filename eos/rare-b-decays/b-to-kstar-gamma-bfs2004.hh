/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_GAMMA_BFS2004_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_GAMMA_BFS2004_HH 1

#include <eos/rare-b-decays/b-to-kstar-gamma-base.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>

namespace eos
{
    template <>
    class BToKstarGammaAmplitudes<tag::BFS2004> :
        public BToKstarGamma::AmplitudeGenerator
    {
        public:
            UsedParameter hbar;

            UsedParameter m_b_MSbar;
            UsedParameter m_c;
            UsedParameter m_s_MSbar;

            UsedParameter f_B;
            UsedParameter f_Kstar_par;
            UsedParameter f_Kstar_perp;
            UsedParameter lambda_B_p_inv;
            UsedParameter a_1_para;
            UsedParameter a_2_para;
            UsedParameter a_1_perp;
            UsedParameter a_2_perp;

            UsedParameter uncertainty_para;
            UsedParameter uncertainty_perp;

            std::shared_ptr<FormFactors<PToV>> form_factors;

            UsedParameter mu;

            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &)> qcdf_photon_massless_case;
            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &)> qcdf_photon_charm_case;
            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &)> qcdf_photon_bottom_case;

            BToKstarGammaAmplitudes(const Parameters & p, const Options & o);
            ~BToKstarGammaAmplitudes() = default;

            virtual BToKstarGamma::Amplitudes amplitudes() const;

            double xi_perp() const;
            double mu_f() const;
            double m_b_PS() const;
    };
}

#endif
