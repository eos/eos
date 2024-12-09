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
            BooleanOption opt_B_bar;

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
            TopologicalRepresentation(const Parameters & p, const Options & o) :
                model(Model::make(o.get("model", "SM"), p, o)),
                opt_q(o, options, "q"),
                opt_p1(o, options, "P1"),
                opt_p2(o, options, "P2"),
                opt_cp_conjugate(o, options, "cp-conjugate"),
                opt_B_bar(o, options, "B_bar"),
                theta_18(p["eta::theta_18"], *this),
                B(su3f::psd_b_triplet.find(opt_q.value())->second),
                H1tilde({}),
                P1({}),
                P2({}),
                Hbar({}),
                H3tilde({}),
                Gfermi(p["WET::G_Fermi"], *this),
                re_T(  p["nonleptonic::Re{T}@Topological"],   *this),
                im_T(  p["nonleptonic::Im{T}@Topological"],   *this),
                re_C(  p["nonleptonic::Re{C}@Topological"],   *this),
                im_C(  p["nonleptonic::Im{C}@Topological"],   *this),
                re_A(  p["nonleptonic::Re{A}@Topological"],   *this),
                im_A(  p["nonleptonic::Im{A}@Topological"],   *this),
                re_E(  p["nonleptonic::Re{E}@Topological"],   *this),
                im_E(  p["nonleptonic::Im{E}@Topological"],   *this),
                re_TES(p["nonleptonic::Re{TES}@Topological"], *this),
                im_TES(p["nonleptonic::Im{TES}@Topological"], *this),
                re_TAS(p["nonleptonic::Re{TAS}@Topological"], *this),
                im_TAS(p["nonleptonic::Im{TAS}@Topological"], *this),
                re_TS( p["nonleptonic::Re{TS}@Topological"],  *this),
                im_TS( p["nonleptonic::Im{TS}@Topological"],  *this),
                re_TPA(p["nonleptonic::Re{TPA}@Topological"], *this),
                im_TPA(p["nonleptonic::Im{TPA}@Topological"], *this),
                re_TP( p["nonleptonic::Re{TP}@Topological"],  *this),
                im_TP( p["nonleptonic::Im{TP}@Topological"],  *this),
                re_TSS(p["nonleptonic::Re{TSS}@Topological"], *this),
                im_TSS(p["nonleptonic::Im{TSS}@Topological"], *this),
                re_P(  p["nonleptonic::Re{P}@Topological"],   *this),
                im_P(  p["nonleptonic::Im{P}@Topological"],   *this),
                re_PT( p["nonleptonic::Re{PT}@Topological"],  *this),
                im_PT( p["nonleptonic::Im{PT}@Topological"],  *this),
                re_S(  p["nonleptonic::Re{S}@Topological"],   *this),
                im_S(  p["nonleptonic::Im{S}@Topological"],   *this),
                re_PC( p["nonleptonic::Re{PC}@Topological"],  *this),
                im_PC( p["nonleptonic::Im{PC}@Topological"],  *this),
                re_PTA(p["nonleptonic::Re{PTA}@Topological"], *this),
                im_PTA(p["nonleptonic::Im{PTA}@Topological"], *this),
                re_PA( p["nonleptonic::Re{PA}@Topological"],  *this),
                im_PA( p["nonleptonic::Im{PA}@Topological"],  *this),
                re_PTE(p["nonleptonic::Re{PTE}@Topological"], *this),
                im_PTE(p["nonleptonic::Im{PTE}@Topological"], *this),
                re_PAS(p["nonleptonic::Re{PAS}@Topological"], *this),
                im_PAS(p["nonleptonic::Im{PAS}@Topological"], *this),
                re_PSS(p["nonleptonic::Re{PSS}@Topological"], *this),
                im_PSS(p["nonleptonic::Im{PSS}@Topological"], *this),
                re_PES(p["nonleptonic::Re{PES}@Topological"], *this),
                im_PES(p["nonleptonic::Im{PES}@Topological"], *this)
            {
                Context ctx("When constructing B->PP topological amplitudes");

                if (opt_cp_conjugate.value() ^ opt_B_bar.value())
                {
                    lamdu = model->ckm_ub() * conj(model->ckm_ud());
                    lamsu = model->ckm_ub() * conj(model->ckm_us());
                    lamdt = model->ckm_tb() * conj(model->ckm_td());
                    lamst = model->ckm_tb() * conj(model->ckm_ts());
                }
                else
                {
                    lamdu = conj(model->ckm_ub()) * model->ckm_ud();
                    lamsu = conj(model->ckm_ub()) * model->ckm_us();
                    lamdt = conj(model->ckm_tb()) * model->ckm_td();
                    lamst = conj(model->ckm_tb()) * model->ckm_ts();
                }

                H1tilde[1] = lamdt;
                H1tilde[2] = lamst;
                Hbar[0][1][0] = lamdu;
                Hbar[0][2][0] = lamsu;
                H3tilde[0][1][0] = lamdt;
                H3tilde[0][2][0] = lamst;
            };

            void update() const
            {
                const double theta_18 = this->theta_18.evaluate();
                su3f::psd_octet.find(opt_p1.value())->second(theta_18, P1);
                su3f::psd_octet.find(opt_p2.value())->second(theta_18, P2);
            }

            ~TopologicalRepresentation() {};

            static NonleptonicAmplitudes<PToPP> * make(const Parameters &, const Options &);

            // Helper functions
            complex<double> tree_amplitude(su3f::rank2 p1, su3f::rank2 p2) const;
            complex<double> penguin_amplitude(su3f::rank2 p1, su3f::rank2 p2) const;

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
