/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_GAMMA_BASE_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_GAMMA_BASE_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/rare-b-decays/b-to-kstar-gamma.hh>
#include <eos/models/model.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    class BToKstarGamma::AmplitudeGenerator :
        public ParameterUser
    {
        public:
            std::shared_ptr<Model> model;
            std::shared_ptr<FormFactors<PToV>> form_factors;

            UsedParameter hbar;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;

            QuarkFlavorOption q;
            UsedParameter tau;
            UsedParameter m_B;
            UsedParameter m_Kstar;

            LeptonFlavorOption l;
            UsedParameter m_l;

            BooleanOption opt_cp_conjugate;
            bool cp_conjugate;
            double e_q;

            static const std::vector<OptionSpecification> options;

            AmplitudeGenerator(const Parameters &, const Options &);

            virtual ~AmplitudeGenerator();
            virtual BToKstarGamma::Amplitudes amplitudes() const = 0;
    };

    template <typename Tag_> class BToKstarGammaAmplitudes;

    namespace tag
    {
        struct BFS2004;
    }
}

#endif
