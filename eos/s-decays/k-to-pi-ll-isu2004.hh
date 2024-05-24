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

#ifndef MASTER_GUARD_EOS_S_DECAYS_K_TO_PI_LL_ISU2004_HH
#define MASTER_GUARD_EOS_S_DECAYS_K_TO_PI_LL_ISU2004_HH 1

#include <eos/s-decays/k-to-pi-ll-base.hh>
#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    template <>
    class KToPiDileptonAmplitudes<tag::ISU2004> :
        public KToPiDilepton::AmplitudeGenerator
    {
        public:
            UsedParameter f_K;

            QuarkFlavorOption q;

            static const std::vector<OptionSpecification> options;

            KToPiDileptonAmplitudes(const Parameters & p, const Options & o);
            ~KToPiDileptonAmplitudes();

            virtual KToPiDilepton::Amplitudes amplitudes(const double & q2) const;
    };
}

#endif
