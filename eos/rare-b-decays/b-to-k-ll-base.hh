/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2015, 2016 Danny van Dyk
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

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_K_LL_BASE_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_K_LL_BASE_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/models/model.hh>
#include <eos/rare-b-decays/b-to-k-ll.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    class BToKDilepton::AmplitudeGenerator :
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
            UsedParameter tau;

            UsedParameter m_B;
            UsedParameter m_K;
            UsedParameter m_l;

            bool cp_conjugate;
            LeptonFlavor lepton_flavor;

            static const std::vector<OptionSpecification> options;

            AmplitudeGenerator(const Parameters &, const Options &);

            double beta_l(const double & q2) const;
            double energy(const double & q2) const;
            double lambda(const double & q2) const;
            double xi_pseudo(const double & q2) const;
            double normalisation(const double & q2) const;

            virtual ~AmplitudeGenerator();
            virtual BToKDilepton::Amplitudes amplitudes(const double & q2) const = 0;
    };

    struct BToKDilepton::DipoleFormFactors
    {
        complex<double> calT;
    };

    template <typename Tag_> class BToKDileptonAmplitudes;

    namespace tag
    {
        /*
         * Approaches that work at small q^2, or large recoil.
         */
        struct BFS2004;
        struct GvDV2020;

        /*
         * Approaches that work at large q^2, or low recoil.
         */
        struct GP2004;
    }
}

#endif
