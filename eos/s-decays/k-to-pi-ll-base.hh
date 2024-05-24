/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Danny van Dyk
 * Copyright (c) 2021 Méril Reboud
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

#ifndef MASTER_GUARD_EOS_S_DECAYS_K_TO_PI_LL_BASE_HH
#define MASTER_GUARD_EOS_S_DECAYS_K_TO_PI_LL_BASE_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/models/model.hh>
#include <eos/s-decays/k-to-pi-ll.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    class KToPiDilepton::AmplitudeGenerator :
        public ParameterUser
    {
        public:
            std::shared_ptr<Model> model;
            std::shared_ptr<FormFactors<PToP>> form_factors;
            LeptonFlavorOption opt_l;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;
            UsedParameter hbar;

            QuarkFlavorOption opt_q;
            UsedParameter tau;
            UsedParameter m_K;
            UsedParameter m_pi;
            UsedParameter m_l;

            BooleanOption opt_cp_conjugate;
            bool cp_conjugate;
            LeptonFlavor lepton_flavor;

            static const std::vector<OptionSpecification> options;

            AmplitudeGenerator(const Parameters &, const Options &);

            double beta_l(const double & q2) const;
            double energy(const double & q2) const;
            double lambda(const double & q2) const;
            double normalisation(const double & q2) const;

            virtual ~AmplitudeGenerator();
            virtual KToPiDilepton::Amplitudes amplitudes(const double & q2) const = 0;
    };

    struct KToPiDilepton::DipoleFormFactors
    {
        complex<double> calT;
    };

    template <typename Tag_> class KToPiDileptonAmplitudes;

    namespace tag
    {
        /*
         * Approaches that ...
         */
        struct ISU2004;
    }
}

#endif
