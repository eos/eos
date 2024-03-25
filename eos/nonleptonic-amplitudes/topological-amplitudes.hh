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
    namespace su3f
    {
        // Meson isospin structures
        typedef std::array<std::array<std::array<complex<double>, 3>, 3>, 3> rank3;
        typedef std::array<std::array<complex<double>, 3>, 3> rank2;
        typedef std::array<complex<double>, 3> rank1;

        const std::map<QuarkFlavor, rank1> psd_b_triplet
        {
            { QuarkFlavor::up,      {{1.0, 0.0, 0.0}}},
            { QuarkFlavor::down,    {{0.0, 1.0, 0.0}}},
            { QuarkFlavor::strange, {{0.0, 0.0, 1.0}}},
        };

        const std::map<LightMeson, rank2> psd_octet
        {
            { LightMeson::pi0,    {{{1.0 / sqrt(2.0), 0.0, 0.0}, {0.0, -1.0 / sqrt(2.0), 0.0}, {0.0, 0.0, 0.0}}} },
            { LightMeson::piplus, {{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}} },
            //...
        };
    }

    template <typename Transition_> class TopologicalRepresentation;
    // template <typename Transition_> class SU3Representation;

    template <>
    class TopologicalRepresentation<PToPP> :
        public NonleptonicAmplitudes<PToPP>
    {
        private:
            std::shared_ptr<Model> model;
            QuarkFlavorOption opt_q;
            LightMesonOption opt_p1;
            LightMesonOption opt_p2;

            su3f::rank1 B, H1tilde;
            su3f::rank2 P1, P2;
            su3f::rank3 Hbar, H3tilde;

            UsedParameter Gfermi;

            UsedParameter T;
            UsedParameter C;
            UsedParameter A;
            UsedParameter E;
            UsedParameter TES;
            UsedParameter TAS;
            UsedParameter TS;
            UsedParameter TPA;
            UsedParameter TP;
            UsedParameter TSS;

            UsedParameter P;
            UsedParameter PT;
            UsedParameter S;
            UsedParameter PC;
            UsedParameter PTA;
            UsedParameter PA;
            UsedParameter PTE;
            UsedParameter PAS;
            UsedParameter PSS;
            UsedParameter PES;

            static const std::vector<OptionSpecification> options;

        public:
            TopologicalRepresentation(const Parameters & p, const Options & o) :
                model(Model::make(o.get("model", "SM"), p, o)),
                opt_q(o, options, "q"),
                opt_p1(o, options, "P1"),
                opt_p2(o, options, "P2"),
                B(su3f::psd_b_triplet.find(opt_q.value())->second),
                H1tilde({}),
                P1(su3f::psd_octet.find(opt_p1.value())->second),
                P2(su3f::psd_octet.find(opt_p2.value())->second),
                Hbar({}),
                H3tilde({}),
                Gfermi(p["WET::G_Fermi"], *this),
                T(p["nonleptonic::T@Topological"], *this),
                C(p["nonleptonic::C@Topological"], *this),
                A(p["nonleptonic::A@Topological"], *this),
                E(p["nonleptonic::E@Topological"], *this),
                TES(p["nonleptonic::TES@Topological"], *this),
                TAS(p["nonleptonic::TAS@Topological"], *this),
                TS(p["nonleptonic::TS@Topological"], *this),
                TPA(p["nonleptonic::TPA@Topological"], *this),
                TP(p["nonleptonic::TP@Topological"], *this),
                TSS(p["nonleptonic::TSS@Topological"], *this),
                P(p["nonleptonic::P@Topological"], *this),
                PT(p["nonleptonic::PT@Topological"], *this),
                S(p["nonleptonic::S@Topological"], *this),
                PC(p["nonleptonic::PC@Topological"], *this),
                PTA(p["nonleptonic::PTA@Topological"], *this),
                PA(p["nonleptonic::PA@Topological"], *this),
                PTE(p["nonleptonic::PTE@Topological"], *this),
                PAS(p["nonleptonic::PAS@Topological"], *this),
                PSS(p["nonleptonic::PSS@Topological"], *this),
                PES(p["nonleptonic::PES@Topological"], *this)
            {
                H1tilde[2] = model->ckm_tb() * conj(model->ckm_td());
                H1tilde[3] = model->ckm_tb() * conj(model->ckm_ts());
                Hbar[1][2][1] = model->ckm_ub() * conj(model->ckm_ud());
                Hbar[1][3][1] = model->ckm_ub() * conj(model->ckm_us());
                H3tilde[1][2][1] = model->ckm_tb() * conj(model->ckm_td());
                H3tilde[1][3][1] = model->ckm_tb() * conj(model->ckm_ts());
            };

            ~TopologicalRepresentation() {};

            static NonleptonicAmplitudes<PToPP> * make(const Parameters &, const Options &);

            complex<double> tree_amplitude() const;
            complex<double> penguin_amplitude() const;

            complex<double> amplitude() const;
    };

    // template <>
    // class SU3Representation<PToPP> :
    //     public NonleptonicAmplitudes<PToPP>
    // {
    //     public:
    //         SU3Representation(const Parameters & p, const Options &);

    //         ~SU3Representation() {};

    //         complex<double> tree_amplitude() const;
    //         complex<double> penguin_amplitude() const;

    //         complex<double> amplitude() const;
    // };
}
#endif
