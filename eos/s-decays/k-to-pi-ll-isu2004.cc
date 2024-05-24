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

#include <eos/maths/power-of.hh>
#include <eos/s-decays/k-to-pi-ll-isu2004.hh>

namespace eos
{
    using namespace std::placeholders;

    KToPiDileptonAmplitudes<tag::ISU2004>::KToPiDileptonAmplitudes(const Parameters & p,
            const Options & o) :
        AmplitudeGenerator(p, o),
        f_K(p["decay-constant::K_" + o.get("q", "d")], *this),
        q(o, options, "q")
    {
        Context ctx("When constructing K->pill ISU2004 amplitudes");
    }

    KToPiDileptonAmplitudes<tag::ISU2004>::~KToPiDileptonAmplitudes()
    {
    }

    const std::vector<OptionSpecification>
    KToPiDileptonAmplitudes<tag::ISU2004>::options
    {
        { "q", { "d", "u" },  "d" },
    };

    /* Amplitudes */
    KToPiDilepton::Amplitudes
    KToPiDileptonAmplitudes<tag::ISU2004>::amplitudes(const double & s) const
    {
        KToPiDilepton::Amplitudes result;

        //WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        result.F_A  = -0.1;
        result.F_T  = 0.0;
        result.F_T5 = 0.0;
        result.F_S  = 0.0;
        result.F_P  = 0.0;
        result.F_V  = +0.1;

        return result;
    }
}
