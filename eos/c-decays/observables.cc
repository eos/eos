/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2023-2025 Danny van Dyk
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

#include <eos/observable-impl.hh>
#include <eos/c-decays/dq-to-l-nu.hh>
#include <eos/c-decays/dstarq-to-l-nu.hh>
#include <eos/c-decays/d-to-psd-l-nu.hh>
#include <eos/c-decays/lambdac-to-lambda-l-nu.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Leptonic D decays
    // {{{
    ObservableGroup
    make_dq_to_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $D_q^{(*)+}\to \ell^+\nu$ decays)",
            R"(The option "l" selects the charged lepton flavor.)",
            {
                make_observable("D->lnu::BR", R"(\mathcal{B}(D^+ \to \ell^+\nu))",
                        Unit::None(),
                        &DqToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q"_ok, "d" } }),

                make_observable("D^*->lnu::BR", R"(\mathcal{B}(D^{*+} \to \ell^+\nu))",
                        Unit::None(),
                        &DstarqToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q"_ok, "d" } }),

                make_observable("D_s->lnu::BR", R"(\mathcal{B}(D_s^+ \to \ell^+\nu))",
                        Unit::None(),
                        &DqToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } }),

                make_observable("D_s^*->lnu::BR", R"(\mathcal{B}(D_s^{*+} \to \ell^+\nu))",
                        Unit::None(),
                        &DstarqToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q"_ok, "s" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Semileptonic D -> P(seudoscalar) decays
    // {{{

    // D -> K l nu
    // {{{
    ObservableGroup
    make_d_to_k_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $D\to K \ell^+ \nu$ decays)",
            R"(The option "l" selects the charged lepton flavor. The option "q" selects the spectator quark flavor. )"
            R"(The option "form-factors" selects the form factor parametrization.)",
            {
                make_observable("D->Klnu::dBR/dq2", R"(d\mathcal{B}(D\to K\ell^+ \nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &DToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "Q"_ok, "s" }, { "I"_ok, "1/2" } }),

                make_observable("D->Klnu::BR", R"(\mathcal{B}(D\to K\ell^+ \nu))",
                        Unit::None(),
                        &DToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "Q"_ok, "s" }, { "I"_ok, "1/2" } }),

                make_observable("D->Klnu::width", R"(\Gamma(D\to K\ell^+ \nu))",
                        Unit::None(),
                        &DToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "Q"_ok, "s" }, { "I"_ok, "1/2" } }),

                make_observable("D->Klnu::P(q2_min,q2_max)", R"(P(D\to K\ell^+ \nu))",
                        Unit::None(),
                        &DToPseudoscalarLeptonNeutrino::integrated_pdf_q2,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "Q"_ok, "s" }, { "I"_ok, "1/2" } }),

                make_observable("D->Klnu::P(q2)", R"(dP(D\to K\ell^+ \nu)/dq^2)",
                        Unit::InverseGeV2(),
                        &DToPseudoscalarLeptonNeutrino::differential_pdf_q2,
                        std::make_tuple("q2"),
                        Options{ { "Q"_ok, "s" }, { "I"_ok, "1/2" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Lambda_c decays
    // {{{
    ObservableGroup
    make_lambdac_to_lambda_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\Lambda_c \to \Lambda \ell^+ \nu$ decays)",
            R"(The option "l" selects the charged lepton flavor.)",
            {
                make_observable("Lambda_c->Lambdalnu::BR", R"(\mathcal{B}(\Lambda_c^+ \to \Lambda \ell^+ \nu))",
                        Unit::None(),
                        &LambdaCToLambdaLeptonNeutrino::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{}),

                make_observable("Lambda_c->Lambdalnu::dBR/dq2", R"(d\mathcal{B}/dq^2(\Lambda_c^+ \to \Lambda \ell^+ \nu))",
                        Unit::InverseGeV2(),
                        &LambdaCToLambdaLeptonNeutrino::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{}),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_c_decays_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in (semi)leptonic $c$-hadron decays",
            "",
            {
                // D_q^+ -> l^+ nu
                make_dq_to_l_nu_group(),

                // D -> K l^+ nu
                make_d_to_k_l_nu_group(),

                // Lc -> L l^+ nu
                make_lambdac_to_lambda_l_nu_group()
            }
        );

        return ObservableSection(imp);
    }
}
