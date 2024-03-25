/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Marta Burgos
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

#ifndef EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_SU3_AMPLITUDES_HH
#define EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_SU3_AMPLITUDES_HH 1

#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes-fwd.hh>
#include <eos/maths/complex.hh>
#include <eos/models/model.hh>
#include <eos/utils/parameters.hh>

#include <array>
#include <map>

namespace eos
{
    template <typename Transition_> class SU3FRepresentation;

    template <>
    class SU3FRepresentation<PToPP> :
        public NonleptonicAmplitudes<PToPP>
    {
        private:
            std::shared_ptr<Model> model;
            QuarkFlavorOption opt_q;
            LightMesonOption opt_p1;
            LightMesonOption opt_p2;
            BooleanOption opt_cp_conjugate;

            UsedParameter theta_18;

            su3f::rank1 B;
            mutable su3f::rank1 H3bar, H3tilde;
            mutable su3f::rank2 P1, P2;
            mutable su3f::rank3 H6bar, H6tilde, H15bar, H15tilde;

            UsedParameter Gfermi;

            UsedParameter re_AT3, im_AT3;
            UsedParameter re_CT3, im_CT3;
            UsedParameter re_AT6, im_AT6;
            UsedParameter re_CT6, im_CT6;
            UsedParameter re_AT15, im_AT15;
            UsedParameter re_CT15, im_CT15;
            UsedParameter re_BT3, im_BT3;
            UsedParameter re_BT6, im_BT6;
            UsedParameter re_BT15, im_BT15;
            UsedParameter re_DT3, im_DT3;

            UsedParameter re_AP3, im_AP3;
            UsedParameter re_CP3, im_CP3;
            UsedParameter re_AP6, im_AP6;
            UsedParameter re_CP6, im_CP6;
            UsedParameter re_AP15, im_AP15;
            UsedParameter re_CP15, im_CP15;
            UsedParameter re_BP3, im_BP3;
            UsedParameter re_BP6, im_BP6;
            UsedParameter re_BP15, im_BP15;
            UsedParameter re_DP3, im_DP3;

            static const std::vector<OptionSpecification> options;

            complex<double> lamdu;
            complex<double> lamsu;
            complex<double> lamdt;
            complex<double> lamst;

        public:
            SU3FRepresentation(const Parameters & p, const Options & o);

            ~SU3FRepresentation() {};

            inline void update() const
            {
                const double theta_18 = this->theta_18.evaluate();
                su3f::psd_octet.find(opt_p1.value())->second(theta_18, P1);
                su3f::psd_octet.find(opt_p2.value())->second(theta_18, P2);
            }

            static NonleptonicAmplitudes<PToPP> * make(const Parameters &, const Options &);

            // Helper functions
            complex<double> tree_amplitude(su3f::rank2 & p1, su3f::rank2 & p2) const;
            complex<double> penguin_amplitude(su3f::rank2 & p1, su3f::rank2 & p2) const;

            // Diagnostic functions
            complex<double> tree_amplitude() const { update(); return tree_amplitude(P1, P2); }
            complex<double> penguin_amplitude() const { update(); return penguin_amplitude(P1, P2); };

            // Amplitude for B -> P1 P2
            complex<double> ordered_amplitude() const;
            // Amplitude for B -> P2 P1
            complex<double> inverse_amplitude() const;
    };
}
#endif
