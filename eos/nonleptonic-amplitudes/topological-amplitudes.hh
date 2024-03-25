/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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

#ifndef EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_TOPOLOGICAL_AMPLITUDES_HH
#define EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_TOPOLOGICAL_AMPLITUDES_HH 1

#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes-fwd.hh>
#include <eos/maths/complex.hh>
#include <eos/models/model.hh>
#include <eos/utils/parameters.hh>

#include <array>
#include <map>

namespace eos
{
    template <typename Transition_> class TopologicalRepresentation;

    template <>
    class TopologicalRepresentation<PToPP> :
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
            mutable su3f::rank1 H1tilde;
            mutable su3f::rank2 P1, P2;
            mutable su3f::rank3 Hbar, H3tilde;

            UsedParameter Gfermi;

            UsedParameter re_T, im_T;
            UsedParameter re_C, im_C;
            UsedParameter re_A, im_A;
            UsedParameter re_E, im_E;
            UsedParameter re_TES, im_TES;
            UsedParameter re_TAS, im_TAS;
            UsedParameter re_TS, im_TS;
            UsedParameter re_TPA, im_TPA;
            UsedParameter re_TP, im_TP;
            UsedParameter re_TSS, im_TSS;

            UsedParameter re_P, im_P;
            UsedParameter re_PT, im_PT;
            UsedParameter re_S, im_S;
            UsedParameter re_PC, im_PC;
            UsedParameter re_PTA, im_PTA;
            UsedParameter re_PA, im_PA;
            UsedParameter re_PTE, im_PTE;
            UsedParameter re_PAS, im_PAS;
            UsedParameter re_PSS, im_PSS;
            UsedParameter re_PES, im_PES;

            static const std::vector<OptionSpecification> options;

            complex<double> lamdu;
            complex<double> lamsu;
            complex<double> lamdt;
            complex<double> lamst;

        public:
            TopologicalRepresentation(const Parameters & p, const Options & o);

            ~TopologicalRepresentation() {};

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
