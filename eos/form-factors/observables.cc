/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2019-2025 Danny van Dyk
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
#include <eos/form-factors/form-factor-adapter.hh>
#include <eos/form-factors/analytic-b-to-psd-dkmmo2008.hh>
#include <eos/form-factors/analytic-b-to-pi-pi.hh>
#include <eos/form-factors/analytic-b-to-v-lcsr.hh>
#include <eos/form-factors/analytic-b-to-p-lcsr.hh>
#include <eos/form-factors/heavy-meson-lcdas.hh>
#include <eos/form-factors/heavy-meson-lcdas-flvd2022.hh>
#include <eos/form-factors/observables.hh>
#include <eos/form-factors/parametric-abr2022.hh>
#include <eos/form-factors/parametric-bgjvd2019.hh>
#include <eos/form-factors/parametric-bgl1997.hh>
#include <eos/form-factors/parametric-bfw2010.hh>
#include <eos/form-factors/parametric-bmrvd2022.hh>
#include <eos/form-factors/parametric-fvdv2018.hh>
#include <eos/form-factors/parametric-hkvt2025.hh>
#include <eos/form-factors/parametric-kkrvd2024.hh>
#include <eos/form-factors/parametric-kkvdz2022.hh>
#include <eos/form-factors/parametric-ksvd2025.hh>
#include <eos/form-factors/unitarity-bounds.hh>
#include <eos/form-factors/zero-recoil-sum-rule.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    /* form factors as observables */
    template <typename Transition_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_form_factor_adapter(const char * name,
            const char * latex,
            double (FormFactors<Transition_>::* _function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        QualifiedName qn(name);
        qnp::Prefix pp = qn.prefix_part();
        std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> function(_function);

        auto result = std::make_pair(qn, std::make_shared<FormFactorAdapterEntry<Transition_, Args_ ...>>(qn, latex, Unit::None(), pp, function, kinematics_names));

        impl::observable_entries.insert(result);

        return result;
    }

    template <typename Transition_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_form_factor_adapter(const char * name,
            double (FormFactors<Transition_>::* _function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        QualifiedName qn(name);
        qnp::Prefix pp = qn.prefix_part();
        std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> function(_function);

        auto result = std::make_pair(qn, std::make_shared<FormFactorAdapterEntry<Transition_, Args_ ...>>(qn, "", Unit::None(), pp, function, kinematics_names));

        impl::observable_entries.insert(result);

        return result;
    }

    // B -> P(seudoscalar)
    // {{{

    // B -> pi
    // {{{
    ObservableGroup
    make_b_to_pi_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \pi$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to \pi$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->pi::f_+(q2)", R"(f_+^{B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B->pi::f_+'(q2)", R"(f_+^{',B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_p_d1, std::make_tuple("q2")),

                make_form_factor_adapter("B->pi::f_+''(q2)", R"(f_+^{'',B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_p_d2, std::make_tuple("q2")),

                make_form_factor_adapter("B->pi::f_T(q2)", R"(f_T^{B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B->pi::f_0(q2)", R"(f_0^{B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->pi::f_-(q2)", R"(f_-^{B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                make_expression_observable("B->pi::f_0(q2)/f_+(q2)", R"(f_0(q^2)/f_+(q^2))",
                        Unit::None(),
                        R"( <<B->pi::f_0(q2)>> / <<B->pi::f_+(q2)>> )"),

                // auxiliary variables, e.g. for determining the B-LCSR threshold parameters
                make_observable("B->pi::f_+[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BToPi>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B->pi::f_0[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BToPi>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B->pi::f_T[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BToPi>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),

                // auxiliary variables, e.g. for determining the pi-LCSR/SVZ threshold parameters
                make_observable("B->pi::M_B(f_+,LCSR)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down>::MBp_lcsr,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "d" } }),

                make_observable("B->pi::M_B(f_0,LCSR)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down>::MB0_lcsr,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "d" } }),

                make_observable("B->pi::M_B(f_T,LCSR)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down>::MBT_lcsr,
                        std::make_tuple("q2"),
                        Options{ { "q"_ok, "d" } }),

                make_observable("B->pi::M_B(SVZ)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down>::MB_svz,
                        std::make_tuple(),
                        Options{ { "q"_ok, "d" } }),

                make_observable("B->pi::f_B@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down>::decay_constant,
                        std::make_tuple(),
                        Options{ { "q"_ok, "d" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> eta
    // {{{
    ObservableGroup
    make_b_to_eta_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \eta$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to \eta$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->eta::f_+(q2)", R"(f_+^{B\to\eta}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B->eta::f_T(q2)", R"(f_T^{B\to\eta}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B->eta::f_0(q2)", R"(f_0^{B\to\eta}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->eta_prime::f_+(q2)", R"(f_+^{B\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B->eta_prime::f_T(q2)", R"(f_T^{B\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B->eta_prime::f_0(q2)", R"(f_0^{B\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_observable("B->eta::Saturation[0^+_V]", R"(\textrm{Saturation}[0^+_V])", Unit::None(),
                        &BFW2010FormFactors<BToEta, PToP>::saturation_0p_v),

                make_observable("B->eta::Saturation[0^-_A]", R"(\textrm{Saturation}[0^-_A])", Unit::None(),
                        &BFW2010FormFactors<BToEta, PToP>::saturation_0m_a),

                make_observable("B->eta::Saturation[1^-_V]", R"(\textrm{Saturation}[1^-_V])", Unit::None(),
                        &BFW2010FormFactors<BToEta, PToP>::saturation_1m_v),

                make_observable("B->eta::Saturation[1^+_A]", R"(\textrm{Saturation}[1^+_A])", Unit::None(),
                        &BFW2010FormFactors<BToEta, PToP>::saturation_1p_a),

                make_observable("B->eta::Saturation[1^-_T]", R"(\textrm{Saturation}[1^-_T])", Unit::None(),
                        &BFW2010FormFactors<BToEta, PToP>::saturation_1m_t),

                make_observable("B->eta::Saturation[1^+_T5]", R"(\textrm{Saturation}[1^+_{T5}])", Unit::None(),
                        &BFW2010FormFactors<BToEta, PToP>::saturation_1p_t5),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> K
    // {{{
    ObservableGroup
    make_b_to_k_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to K$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to K$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->K::f_+(q2)", R"(f_+^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B->K::f_0(q2)", R"(f_0^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->K::f_T(q2)", R"(f_T^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B->K::f_-(q2)", R"(f_-^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                make_observable("B->K::f_+[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BToK>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B->K::f_0[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BToK>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B->K::f_T[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BToK>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->K::F_plus(q2)", R"(F_+^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B->K::F_plus_T(q2)", R"(F_T^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_plus_T, std::make_tuple("q2")),

                make_expression_observable("B->K::F_T(q2)/F_plus(q2)", R"(F_T(q^2)/F_+(q^2))",
                        Unit::None(),
                        R"( <<B->K::f_T(q2)>> / <<B->K::F_plus(q2)>> )"),

                make_expression_observable("B->K::F_plus_T(q2)/F_plus(q2)", R"(F_{+,T}(q^2)/F_+(q^2))",
                        Unit::None(),
                        R"( <<B->K::F_plus_T(q2)>> / <<B->K::F_plus(q2)>> )"),

                make_observable("B->K::Saturation[0^+_V]", R"(\textrm{Saturation}[0^+_V])", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::saturation_0p_v),

                make_observable("B->K::Saturation[0^-_A]", R"(\textrm{Saturation}[0^-_A])", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::saturation_0m_a),

                make_observable("B->K::Saturation[1^-_V]", R"(\textrm{Saturation}[1^-_V])", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::saturation_1m_v),

                make_observable("B->K::Saturation[1^+_A]", R"(\textrm{Saturation}[1^+_A])", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::saturation_1p_a),

                make_observable("B->K::Saturation[1^-_T]", R"(\textrm{Saturation}[1^-_T])", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::saturation_1m_t),

                make_observable("B->K::Saturation[1^+_T5]", R"(\textrm{Saturation}[1^+_{T5}])", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::saturation_1p_t5),

                // Auxiliary functions for [BFW:2010A]
                make_observable("B->K::f_+_series(q2)@BFW2010", R"(\hat{f}_+^{B\to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::f_p_series, std::make_tuple("q2")),

                make_observable("B->K::f_+_series_prime(q2)@BFW2010", R"(\hat{f}_+^{\prime B\to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::f_p_series_prime, std::make_tuple("q2")),

                make_observable("B->K::f_0_series(q2)@BFW2010", R"(\hat{f}_0^{B\to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::f_0_series, std::make_tuple("q2")),

                make_observable("B->K::f_0_series_prime(q2)@BFW2010", R"(\hat{f}_0^{\prime B\to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::f_0_series_prime, std::make_tuple("q2")),

                make_observable("B->K::f_T_series(q2)@BFW2010", R"(\hat{f}_T^{B\to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::f_t_series, std::make_tuple("q2")),

                make_observable("B->K::f_T_series_prime(q2)@BFW2010", R"(\hat{f}_T^{\prime B\to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<BToK, PToP>::f_t_series_prime, std::make_tuple("q2"))
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> D
    // {{{
    ObservableGroup
    make_b_to_d_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \bar{D}$ transitions and related pseudo-observables)",
            R"(Pseudo observables representing the full basis of $B\to\bar{D}$ form factors. )"
            R"(For most pseudo-observables, the specific parametrization can be chosen via the "form-factors" option.)",
            {
                // B -> D Form Factors
                make_form_factor_adapter("B->D::f_+(q2)", R"(f_+^{B\to \bar{D}}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B->D::f_0(q2)", R"(f_0^{B\to \bar{D}}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->D::f_T(q2)", R"(f_T^{B\to \bar{D}}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B->D::f_-(q2)", R"(f_-^{B\to \bar{D}}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                make_observable("B->D::f_+[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BToD>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B->D::f_0[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BToD>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B->D::f_T[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BToD>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),

                make_observable("B->D::a_0[S_1]@HQE", R"(a_0^{S_1})",
                        Unit::None(),
                        &BGLCoefficients::S1_a0),

                make_observable("B->D::a_1[S_1]@HQE", R"(a_1^{S_1})",
                        Unit::None(),
                        &BGLCoefficients::S1_a1),

                make_observable("B->D::a_2[S_1]@HQE", R"(a_2^{S_1})",
                        Unit::None(),
                        &BGLCoefficients::S1_a2),

                make_expression_observable("B->D::a_1/a_0[S_1]@HQE", R"(a_1^{S_1}/a_0^{S_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_1[S_1]@HQE>>
                        /
                        <<B->D::a_0[S_1]@HQE>>
                        )"),

                make_expression_observable("B->D::a_2/a_0[S_1]@HQE", R"(a_2^{S_1}/a_0^{S_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_2[S_1]@HQE>>
                        /
                        <<B->D::a_0[S_1]@HQE>>
                        )"),

                make_observable("B->D::a_0[V_1]@HQE", R"(a_0^{V_1})",
                        Unit::None(),
                        &BGLCoefficients::V1_a0),

                make_observable("B->D::a_1[V_1]@HQE", R"(a_1^{V_1})",
                        Unit::None(),
                        &BGLCoefficients::V1_a1),

                make_observable("B->D::a_2[V_1]@HQE", R"(a_2^{V_1})",
                        Unit::None(),
                        &BGLCoefficients::V1_a2),

                make_expression_observable("B->D::a_1/a_0[V_1]@HQE", R"(a_1^{V_1}/a_0^{V_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_1[V_1]@HQE>>
                        /
                        <<B->D::a_0[V_1]@HQE>>
                        )"),

                make_expression_observable("B->D::a_2/a_0[V_1]@HQE", R"(a_2^{V_1}/a_0^{V_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_2[V_1]@HQE>>
                        /
                        <<B->D::a_0[V_1]@HQE>>
                        )"),

                make_observable("B->D::a_0[f_T]@HQE", R"(a_0^{f_T})",
                        Unit::None(),
                        &BGLCoefficients::fT_a0),

                make_observable("B->D::a_1[f_T]@HQE", R"(a_1^{f_T})",
                        Unit::None(),
                        &BGLCoefficients::fT_a1),

                make_observable("B->D::a_2[f_T]@HQE", R"(a_2^{f_T})",
                        Unit::None(),
                        &BGLCoefficients::fT_a2),

                make_expression_observable("B->D::a_1/a_0[f_T]@HQE", R"(a_1^{f_T}/a_0^{f_T})",
                        Unit::None(),
                        R"(
                        <<B->D::a_1[f_T]@HQE>>
                        /
                        <<B->D::a_0[f_T]@HQE>>
                        )"),

                make_expression_observable("B->D::a_2/a_0[f_T]@HQE", R"(a_2^{f_T}/a_0^{f_T})",
                        Unit::None(),
                        R"(
                        <<B->D::a_2[f_T]@HQE>>
                        /
                        <<B->D::a_0[f_T]@HQE>>
                        )"),

                make_expression_observable("B->D::f_T(q2)/f_+(q2)", R"(f_T(q^2)/f_+(q^2))",
                        Unit::None(),
                        R"( <<B->D::f_T(q2)>> / <<B->D::f_+(q2)>> )"),

                make_observable("B->D::h_+(q2)", R"(h_+^{B\to \bar{D}}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToD, PToP>::h_p,
                        std::make_tuple("q2")),

                make_observable("B->D::h_-(q2)", R"(h_-^{B\to \bar{D}}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToD, PToP>::h_m,
                        std::make_tuple("q2")),

                make_observable("B->D::h_S(q2)", R"(h_S^{B\to \bar{D}}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToD, PToP>::h_S,
                        std::make_tuple("q2")),

                make_observable("B->D::h_T(q2)", R"(h_T^{B\to \bar{D}}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToD, PToP>::h_T,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B_s -> P(seudoscalar)
    // {{{

    // B_s -> eta
    // {{{
    ObservableGroup
    make_bs_to_eta_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \eta^{(\prime)}$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to \eta^{(\prime)}$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B_s->eta::f_+(q2)", R"(f_+^{B_s\to\eta}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->eta::f_T(q2)", R"(f_T^{B_s\to\eta}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->eta::f_0(q2)", R"(f_0^{B_s\to\eta}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->eta_prime::f_+(q2)", R"(f_+^{B_s\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->eta_prime::f_T(q2)", R"(f_T^{B_s\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->eta_prime::f_0(q2)", R"(f_0^{B_s\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> K
    // {{{
    ObservableGroup
    make_bs_to_k_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \bar{K}$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to \bar{K}$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B_s->K::f_+(q2)", R"(f_+^{B_s\to \bar{K}}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K::f_0(q2)", R"(f_0^{B_s\to \bar{K}}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K::f_T(q2)", R"(f_T^{B_s\to \bar{K}}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K::f_-(q2)", R"(f_-^{B_s\to \bar{K}}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                // auxiliary variables, e.g. for determining the B-LCSR threshold parameters
                make_observable("B_s->K::f_+[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BsToK>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B_s->K::f_0[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BsToK>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B_s->K::f_T[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BsToK>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),

                // auxiliary variables, e.g. for determining the pi-LCSR/SVZ threshold parameters
                make_observable("B_s->K::M_B(f_+,LCSR)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::strange>::MBp_lcsr,
                        std::make_tuple("q2"),
                        Options{ }),

                make_observable("B_s->K::M_B(f_0,LCSR)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::strange>::MB0_lcsr,
                        std::make_tuple("q2"),
                        Options{ }),

                make_observable("B_s->K::M_B(f_T,LCSR)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::strange>::MBT_lcsr,
                        std::make_tuple("q2"),
                        Options{ }),

                make_observable("B_s->K::M_B(SVZ)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::strange>::MB_svz,
                        std::make_tuple(),
                        Options{ }),

                make_observable("B_s->K::f_B@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::strange>::decay_constant,
                        std::make_tuple(),
                        Options{ }),

                make_observable("B_s->K::Saturation[0^+_V]", R"(\textrm{Saturation}[0^+_V])", Unit::None(),
                        &BFW2010FormFactors<BsToK, PToP>::saturation_0p_v),

                make_observable("B_s->K::Saturation[1^-_V]", R"(\textrm{Saturation}[1^-_V])", Unit::None(),
                        &BFW2010FormFactors<BsToK, PToP>::saturation_1m_v),

                make_observable("B_s->K::Saturation[1^+_A]", R"(\textrm{Saturation}[1^+_A])", Unit::None(),
                        &BFW2010FormFactors<BsToK, PToP>::saturation_1p_a),

                make_observable("B_s->K::Saturation[1^-_T]", R"(\textrm{Saturation}[1^-_T])", Unit::None(),
                        &BFW2010FormFactors<BsToK, PToP>::saturation_1m_t),

                make_observable("B_s->K::Saturation[1^+_T5]", R"(\textrm{Saturation}[1^+_{T5}])", Unit::None(),
                        &BFW2010FormFactors<BsToK, PToP>::saturation_1p_t5),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> D_s
    // {{{
    ObservableGroup
    make_bs_to_ds_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \bar{D_s}$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to\bar{D}_s$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                // B -> D Form Factors
                make_form_factor_adapter("B_s->D_s::f_+(q2)", R"(f_+^{B_s\to \bar{D}_s}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s::f_0(q2)", R"(f_0^{B_s\to \bar{D}_s}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s::f_T(q2)", R"(f_T^{B_s\to \bar{D}_s}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s::f_-(q2)", R"(f_-^{B_s\to \bar{D}_s}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                make_observable("B_s->D_s::f_+[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BsToDs>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s::f_0[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BsToDs>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s::f_T[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<BsToDs>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),

                make_expression_observable("B(_s)->D(_s)::f_0(q2_num)/f_0(q2_denom)", R"(f_0(q^2_\mathrm{num})/f_+(q^2_\mathrm{denom}))",
                        Unit::None(),
                        R"(
                        <<B_s->D_s::f_0(q2)>>[q2=>q2_num]
                        /
                        <<B->D::f_0(q2)>>[q2=>q2_denom]
                        )"),

                make_expression_observable("B_s->D_s::f_T(q2)/f_+(q2)", R"(f_T(q^2)/f_+(q^2))",
                        Unit::None(),
                        R"( <<B_s->D_s::f_T(q2)>> / <<B_s->D_s::f_+(q2)>> )"),


                make_observable("B_s->D_s::a_0[S_1]@HQE", R"(a_0^{S_1})",
                        Unit::None(),
                        &BGLCoefficients::S1s_a0),

                make_observable("B_s->D_s::a_1[S_1]@HQE", R"(a_1^{S_1})",
                        Unit::None(),
                        &BGLCoefficients::S1s_a1),

                make_observable("B_s->D_s::a_2[S_1]@HQE", R"(a_2^{S_1})",
                        Unit::None(),
                        &BGLCoefficients::S1s_a2),

                make_expression_observable("B_s->D_s::a_1/a_0[S_1]@HQE", R"(a_1^{S_1}/a_0^{S_1})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s::a_1[S_1]@HQE>>
                        /
                        <<B_s->D_s::a_0[S_1]@HQE>>
                        )"),

                make_expression_observable("B_s->D_s::a_2/a_0[S_1]@HQE", R"(a_2^{S_1}/a_0^{S_1})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s::a_2[S_1]@HQE>>
                        /
                        <<B_s->D_s::a_0[S_1]@HQE>>
                        )"),

                make_observable("B_s->D_s::a_0[V_1]@HQE", R"(a_0^{V_1})",
                        Unit::None(),
                        &BGLCoefficients::V1s_a0),

                make_observable("B_s->D_s::a_1[V_1]@HQE", R"(a_1^{V_1})",
                        Unit::None(),
                        &BGLCoefficients::V1s_a1),

                make_observable("B_s->D_s::a_2[V_1]@HQE", R"(a_2^{V_1})",
                        Unit::None(),
                        &BGLCoefficients::V1s_a2),

                make_expression_observable("B_s->D_s::a_1/a_0[V_1]@HQE", R"(a_1^{V_1}/a_0^{V_1})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s::a_1[V_1]@HQE>>
                        /
                        <<B_s->D_s::a_0[V_1]@HQE>>
                        )"),

                make_expression_observable("B_s->D_s::a_2/a_0[V_1]@HQE", R"(a_2^{V_1}/a_0^{V_1})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s::a_2[V_1]@HQE>>
                        /
                        <<B_s->D_s::a_0[V_1]@HQE>>
                        )"),

                make_observable("B_s->D_s::a_0[f_T]@HQE", R"(a_0^{f_T})",
                        Unit::None(),
                        &BGLCoefficients::fTs_a0),

                make_observable("B_s->D_s::a_1[f_T]@HQE", R"(a_1^{f_T})",
                        Unit::None(),
                        &BGLCoefficients::fTs_a1),

                make_observable("B_s->D_s::a_2[f_T]@HQE", R"(a_2^{f_T})",
                        Unit::None(),
                        &BGLCoefficients::fTs_a2),

                make_expression_observable("B_s->D_s::a_1/a_0[f_T]@HQE", R"(a_1^{f_T}/a_0^{f_T})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s::a_1[f_T]@HQE>>
                        /
                        <<B_s->D_s::a_0[f_T]@HQE>>
                        )"),

                make_expression_observable("B_s->D_s::a_2/a_0[f_T]@HQE", R"(a_2^{f_T}/a_0^{f_T})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s::a_2[f_T]@HQE>>
                        /
                        <<B_s->D_s::a_0[f_T]@HQE>>
                        )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B -> gamma^*
    // {{{
    ObservableGroup
    make_b_to_gammastar_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B \to \gamma^*$ transitions)",
            R"(Pseudo observables representing the $B \to \gamma^*$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->gamma^*::Arg{F_1}(q2,k2)", R"(\text{Arg}\,F_1^{B\to \gamma^*}(q^2,k^2))",
                        &FormFactors<PToGammaOffShell>::arg_F_1, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->gamma^*::Arg{F_2}(q2,k2)", R"(\text{Arg}\,F_2^{B\to \gamma^*}(q^2,k^2))",
                        &FormFactors<PToGammaOffShell>::arg_F_2, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->gamma^*::Arg{F_3}(q2,k2)", R"(\text{Arg}\,F_3^{B\to \gamma^*}(q^2,k^2))",
                        &FormFactors<PToGammaOffShell>::arg_F_3, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->gamma^*::Arg{F_4}(q2,k2)", R"(\text{Arg}\,F_4^{B\to \gamma^*}(q^2,k^2))",
                        &FormFactors<PToGammaOffShell>::arg_F_4, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->gamma^*::Abs{F_1}(q2,k2)", R"(\text{Abs}\,F_1^{B\to \gamma^*}(q^2,k^2))",
                        &FormFactors<PToGammaOffShell>::abs_F_1, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->gamma^*::Abs{F_2}(q2,k2)", R"(\text{Abs}\,F_2^{B\to \gamma^*}(q^2,k^2))",
                        &FormFactors<PToGammaOffShell>::abs_F_2, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->gamma^*::Abs{F_3}(q2,k2)", R"(\text{Abs}\,F_3^{B\to \gamma^*}(q^2,k^2))",
                        &FormFactors<PToGammaOffShell>::abs_F_3, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->gamma^*::Abs{F_4}(q2,k2)", R"(\text{Abs}\,F_4^{B\to \gamma^*}(q^2,k^2))",
                        &FormFactors<PToGammaOffShell>::abs_F_4, std::make_tuple("q2", "k2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> V(ector)
    // {{{

    // B -> omega
    // {{{
    ObservableGroup
    make_b_to_omega_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \omega$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to \omega$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->omega::V(q2)", R"(V^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::A_0(q2)", R"(A_0^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::A_1(q2)", R"(A_1^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::A_2(q2)", R"(A_2^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::A_12(q2)", R"(A_{12}^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::T_1(q2)", R"(T_1^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::T_2(q2)", R"(T_2^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::T_3(q2)", R"(T_3^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::T_23(q2)", R"(T_{23}^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> gamma
    // {{{
    ObservableGroup
    make_b_to_gamma_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \gamma$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to \gamma$ form factors.)",
            {
                make_form_factor_adapter("B->gamma::F_V(E_gamma)", R"(F_V^{B\to \gamma}(E_\gamma))",
                        &FormFactors<PToGamma>::F_V, std::make_tuple("E_gamma")),

                make_form_factor_adapter("B->gamma::F_A(E_gamma)", R"(F_A^{B\to \gamma}(E_\gamma))",
                        &FormFactors<PToGamma>::F_A, std::make_tuple("E_gamma")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> rho
    // {{{
    ObservableGroup
    make_b_to_rho_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \rho$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to \rho$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->rho::V(q2)", R"(V^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::A_0(q2)", R"(A_0^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::A_1(q2)", R"(A_1^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::A_2(q2)", R"(A_2^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::A_12(q2)", R"(A_{12}^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::T_1(q2)", R"(T_1^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::T_2(q2)", R"(T_2^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::T_3(q2)", R"(T_3^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::T_23(q2)", R"(T_{23}^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_expression_observable("B->rho::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->rho::V(q2)>> / <<B->rho::A_1(q2)>> )"),

                make_expression_observable("B->rho::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->rho::A_2(q2)>> / <<B->rho::A_1(q2)>> )"),

                make_expression_observable("B->rho::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->rho::A_12(q2)>> / <<B->rho::A_1(q2)>> )"),

                make_expression_observable("B->rho::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B->rho::T_23(q2)>> / <<B->rho::T_2(q2)>> )"),

                make_observable("B->rho::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToRho>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B->rho::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToRho>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B->rho::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToRho>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B->rho::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToRho>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B->rho::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToRho>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B->rho::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToRho>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B->rho::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToRho>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> K^*
    // {{{
    ObservableGroup
    make_b_to_kstar_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to K^*$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to K^*$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->K^*::V(q2)", R"(V^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_0(q2)", R"(A_0^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_1(q2)", R"(A_1^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_2(q2)", R"(A_2^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_12(q2)", R"(A_{12}^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_1(q2)", R"(T_1^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_2(q2)", R"(T_2^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_3(q2)", R"(T_3^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_23(q2)", R"(T_{23}^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_perp(q2)", R"(\mathcal{F}_\perp^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_perp, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_para(q2)",  R"(\mathcal{F}_\parallel^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_para, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_long(q2)",  R"(\mathcal{F}_0^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_long, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_perp_T(q2)", R"(\mathcal{F}_{\perp,T}^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_perp_T, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_para_T(q2)", R"(\mathcal{F}_{\parallel,T}^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_para_T, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_long_T(q2)", R"(\mathcal{F}_{0,T}^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_long_T, std::make_tuple("q2")),

                make_expression_observable("B->K^*::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::V(q2)>> / <<B->K^*::A_1(q2)>> )"),

                make_expression_observable("B->K^*::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::A_2(q2)>> / <<B->K^*::A_1(q2)>> )"),

                make_expression_observable("B->K^*::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::A_12(q2)>> / <<B->K^*::A_1(q2)>> )"),

                make_expression_observable("B->K^*::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::T_23(q2)>> / <<B->K^*::T_2(q2)>> )"),

                make_observable("B->K^*::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToKstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B->K^*::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToKstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B->K^*::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToKstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B->K^*::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToKstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B->K^*::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToKstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B->K^*::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToKstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B->K^*::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToKstar>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),

                make_expression_observable("B->K^*::F_perp_T(q2)/F_perp(q2)", R"(\mathcal{F}_{\perp,T}(q^2)/\mathcal{F}_\perp(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::F_perp_T(q2)>> / <<B->K^*::F_perp(q2)>> )"),

                make_expression_observable("B->K^*::F_para_T(q2)/F_para(q2)", R"(\mathcal{F}_{\parallel,T}(q^2)/\mathcal{F}_\parallel(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::F_para_T(q2)>> / <<B->K^*::F_para(q2)>> )"),

                make_expression_observable("B->K^*::F_long_T(q2)/F_long(q2)", R"(\mathcal{F}_{0,T}(q^2)/\mathcal{F}_0(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::F_long_T(q2)>> / <<B->K^*::F_long(q2)>> )"),

                make_observable("B->K^*::Saturation[0^+_V]", R"(\textrm{Saturation}[0^+_V])", Unit::None(),
                        &BFW2010FormFactors<BToKstar, PToV>::saturation_0p_v),

                make_observable("B->K^*::Saturation[0^-_A]", R"(\textrm{Saturation}[0^-_A])", Unit::None(),
                        &BFW2010FormFactors<BToKstar, PToV>::saturation_0m_a),

                make_observable("B->K^*::Saturation[1^-_V]", R"(\textrm{Saturation}[1^-_V])", Unit::None(),
                        &BFW2010FormFactors<BToKstar, PToV>::saturation_1m_v),

                make_observable("B->K^*::Saturation[1^+_A]", R"(\textrm{Saturation}[1^+_A])", Unit::None(),
                        &BFW2010FormFactors<BToKstar, PToV>::saturation_1p_a),

                make_observable("B->K^*::Saturation[1^-_T]", R"(\textrm{Saturation}[1^-_T])", Unit::None(),
                        &BFW2010FormFactors<BToKstar, PToV>::saturation_1m_t),

                make_observable("B->K^*::Saturation[1^+_T5]", R"(\textrm{Saturation}[1^+_{T5}])", Unit::None(),
                        &BFW2010FormFactors<BToKstar, PToV>::saturation_1p_t5)
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> D^*
    // {{{
    ObservableGroup
    make_b_to_dstar_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \bar{D}^*$ transitions and related pseudo-observables)",
            R"(Pseudo observables representing the full basis of $B\to \bar{D}^*$ form factors. )"
            R"(For most pseudo-observables, the specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->D^*::V(q2)", R"(V^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_0(q2)", R"(A_0^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_1(q2)", R"(A_1^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_2(q2)", R"(A_2^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_12(q2)", R"(A_{12}^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_1(q2)", R"(T_1^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_2(q2)", R"(T_2^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_3(q2)", R"(T_3^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_23(q2)", R"(T_{23}^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_expression_observable("B->D^*::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->D^*::V(q2)>> / <<B->D^*::A_1(q2)>> )"),

                make_expression_observable("B->D^*::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->D^*::A_2(q2)>> / <<B->D^*::A_1(q2)>> )"),

                make_expression_observable("B->D^*::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->D^*::A_12(q2)>> / <<B->D^*::A_1(q2)>> )"),

                make_expression_observable("B->D^*::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B->D^*::T_23(q2)>> / <<B->D^*::T_2(q2)>> )"),

                make_observable("B->D^*::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToDstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B->D^*::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToDstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B->D^*::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                         &AnalyticFormFactorBToVLCSR<BToDstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B->D^*::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToDstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B->D^*::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToDstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B->D^*::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToDstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B->D^*::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BToDstar>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),

                make_observable("B->D^*::a_0[A_1]@HQE", R"(a_0^{A_1})",
                        Unit::None(),
                        &BGLCoefficients::A1_a0),

                make_observable("B->D^*::a_1[A_1]@HQE", R"(a_1^{A_1})",
                        Unit::None(),
                        &BGLCoefficients::A1_a1),

                make_observable("B->D^*::a_2[A_1]@HQE", R"(a_2^{A_1})",
                        Unit::None(),
                        &BGLCoefficients::A1_a2),

                make_expression_observable("B->D^*::a_1/a_0[A_1]@HQE", R"(a_1^{A_1}/a_0^{A_1})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_1[A_1]@HQE>>
                        /
                        <<B->D^*::a_0[A_1]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[A_1]@HQE", R"(a_2^{A_1}/a_0^{A_1})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_2[A_1]@HQE>>
                        /
                        <<B->D^*::a_0[A_1]@HQE>>
                        )"),

                make_observable("B->D^*::a_0[A_5]@HQE", R"(a_0^{A_5})",
                        Unit::None(),
                        &BGLCoefficients::A5_a0),

                make_observable("B->D^*::a_1[A_5]@HQE", R"(a_1^{A_5})",
                        Unit::None(),
                        &BGLCoefficients::A5_a1),

                make_observable("B->D^*::a_2[A_5]@HQE", R"(a_2^{A_5})",
                        Unit::None(),
                        &BGLCoefficients::A5_a2),

                make_expression_observable("B->D^*::a_1/a_0[A_5]@HQE", R"(a_1^{A_5}/a_0^{A_5})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_1[A_5]@HQE>>
                        /
                        <<B->D^*::a_0[A_5]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[A_5]@HQE", R"(a_2^{A_5}/a_0^{A_5})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_2[A_5]@HQE>>
                        /
                        <<B->D^*::a_0[A_5]@HQE>>
                        )"),

                make_observable("B->D^*::a_0[P_1]@HQE", R"(a_0^{P_1})",
                        Unit::None(),
                        &BGLCoefficients::P1_a0),

                make_observable("B->D^*::a_1[P_1]@HQE", R"(a_1^{P_1})",
                        Unit::None(),
                        &BGLCoefficients::P1_a1),

                make_observable("B->D^*::a_2[P_1]@HQE", R"(a_2^{P_1})",
                        Unit::None(),
                        &BGLCoefficients::P1_a2),

                make_expression_observable("B->D^*::a_1/a_0[P_1]@HQE", R"(a_1^{P_1}/a_0^{P_1})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_1[P_1]@HQE>>
                        /
                        <<B->D^*::a_0[P_1]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[P_1]@HQE", R"(a_2^{P_1}/a_0^{P_1})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_2[P_1]@HQE>>
                        /
                        <<B->D^*::a_0[P_1]@HQE>>
                        )"),

                make_observable("B->D^*::a_0[V_4]@HQE", R"(a_0^{V_4})",
                        Unit::None(),
                        &BGLCoefficients::V4_a0),

                make_observable("B->D^*::a_1[V_4]@HQE", R"(a_1^{V_4})",
                        Unit::None(),
                        &BGLCoefficients::V4_a1),

                make_observable("B->D^*::a_2[V_4]@HQE", R"(a_2^{V_4})",
                        Unit::None(),
                        &BGLCoefficients::V4_a2),

                make_expression_observable("B->D^*::a_1/a_0[V_4]@HQE", R"(a_1^{V_4}/a_0^{V_4})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_1[V_4]@HQE>>
                        /
                        <<B->D^*::a_0[V_4]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[V_4]@HQE", R"(a_2^{V_4}/a_0^{V_4})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_2[V_4]@HQE>>
                        /
                        <<B->D^*::a_0[V_4]@HQE>>
                        )"),

                make_observable("B->D^*::a_0[T_1]@HQE", R"(a_0^{T_1})",
                        Unit::None(),
                        &BGLCoefficients::T1_a0),

                make_observable("B->D^*::a_1[T_1]@HQE", R"(a_1^{T_1})",
                        Unit::None(),
                        &BGLCoefficients::T1_a1),

                make_observable("B->D^*::a_2[T_1]@HQE", R"(a_2^{T_1})",
                        Unit::None(),
                        &BGLCoefficients::T1_a2),

                make_expression_observable("B->D^*::a_1/a_0[T_1]@HQE", R"(a_1^{T_1}/a_0^{T_1})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_1[T_1]@HQE>>
                        /
                        <<B->D^*::a_0[T_1]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[T_1]@HQE", R"(a_2^{T_1}/a_0^{T_1})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_2[T_1]@HQE>>
                        /
                        <<B->D^*::a_0[T_1]@HQE>>
                        )"),

                make_observable("B->D^*::a_0[T_2]@HQE", R"(a_0^{T_2})",
                        Unit::None(),
                        &BGLCoefficients::T2_a0),

                make_observable("B->D^*::a_1[T_2]@HQE", R"(a_1^{T_2})",
                        Unit::None(),
                        &BGLCoefficients::T2_a1),

                make_observable("B->D^*::a_2[T_2]@HQE", R"(a_2^{T_2})",
                        Unit::None(),
                        &BGLCoefficients::T2_a2),

                make_expression_observable("B->D^*::a_1/a_0[T_2]@HQE", R"(a_1^{T_2}/a_0^{T_2})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_1[T_2]@HQE>>
                        /
                        <<B->D^*::a_0[T_2]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[T_2]@HQE", R"(a_2^{T_2}/a_0^{T_2})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_2[T_2]@HQE>>
                        /
                        <<B->D^*::a_0[T_2]@HQE>>
                        )"),

                make_observable("B->D^*::a_0[T_23]@HQE", R"(a_0^{T_{23}})",
                        Unit::None(),
                        &BGLCoefficients::T23_a0),

                make_observable("B->D^*::a_1[T_23]@HQE", R"(a_1^{T_{23}})",
                        Unit::None(),
                        &BGLCoefficients::T23_a1),

                make_observable("B->D^*::a_2[T_23]@HQE", R"(a_2^{T_{23}})",
                        Unit::None(),
                        &BGLCoefficients::T23_a2),

                make_expression_observable("B->D^*::a_1/a_0[T_23]@HQE", R"(a_1^{T_{23}}/a_0^{T_{23}})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_1[T_23]@HQE>>
                        /
                        <<B->D^*::a_0[T_23]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[T_23]@HQE", R"(a_2^{T_{23}}/a_0^{T_{23}})",
                        Unit::None(),
                        R"(
                        <<B->D^*::a_2[T_23]@HQE>>
                        /
                        <<B->D^*::a_0[T_23]@HQE>>
                        )"),

                make_observable("B->D^*::a_0[F_1]@BGL", R"(a_0^{F_1})",
                        Unit::None(),
                        &BGL1997FormFactors<BToDstar, PToV>::a_F1_0),

                make_observable("B->D^*::a_0[F_2]@BGL", R"(a_0^{F_2})",
                        Unit::None(),
                        &BGL1997FormFactors<BToDstar, PToV>::a_F2_0),

                make_observable("B->D^*::a_0[T_2]@BGL", R"(a_0^{T_2})",
                        Unit::None(),
                        &BGL1997FormFactors<BToDstar, PToV>::a_T2_0),

                make_observable("B->D^*::a_0[T_23]@BGL", R"(a_0^{T_{23}})",
                        Unit::None(),
                        &BGL1997FormFactors<BToDstar, PToV>::a_T23_0),

                make_observable("B->D^*::h_A1(q2)", R"(h_{A_1}^{\bar{B}\to D^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToDstar, PToV>::h_a1,
                        std::make_tuple("q2")),

                make_observable("B->D^*::h_A2(q2)", R"(h_{A_2}^{\bar{B}\to D^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToDstar, PToV>::h_a2,
                        std::make_tuple("q2")),

                make_observable("B->D^*::h_A3(q2)", R"(h_{A_3}^{\bar{B}\to D^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToDstar, PToV>::h_a3,
                        std::make_tuple("q2")),

                make_observable("B->D^*::h_V(q2)", R"(h_V^{\bar{B}\to D^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToDstar, PToV>::h_v,
                        std::make_tuple("q2")),

                make_observable("B->D^*::h_T1(q2)", R"(h_{T_1}^{\bar{B}\to D^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToDstar, PToV>::h_t1,
                        std::make_tuple("q2")),

                make_observable("B->D^*::h_T2(q2)", R"(h_{T_2}^{\bar{B}\to D^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToDstar, PToV>::h_t2,
                        std::make_tuple("q2")),

                make_observable("B->D^*::h_T3(q2)", R"(h_{T_3}^{\bar{B}\to D^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BToDstar, PToV>::h_t3,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B_s -> V(ector)
    // {{{

    // B_s -> K^*
    // {{{
    ObservableGroup
    make_bs_to_kstar_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \bar{K}^*$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to \bar{K}^*$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B_s->K^*::V(q2)", R"(V^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::A_0(q2)", R"(A_0^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::A_1(q2)", R"(A_1^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::A_2(q2)", R"(A_2^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::A_12(q2)", R"(A_{12}^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::T_1(q2)", R"(T_1^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::T_2(q2)", R"(T_2^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::T_3(q2)", R"(T_3^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::T_23(q2)", R"(T_{23}^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_expression_observable("B_s->K^*::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->K^*::V(q2)>> / <<B_s->K^*::A_1(q2)>> )"),

                make_expression_observable("B_s->K^*::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->K^*::A_2(q2)>> / <<B_s->K^*::A_1(q2)>> )"),

                make_expression_observable("B_s->K^*::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->K^*::A_12(q2)>> / <<B_s->K^*::A_1(q2)>> )"),

                make_expression_observable("B_s->K^*::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B_s->K^*::T_23(q2)>> / <<B_s->K^*::T_2(q2)>> )"),

                make_observable("B_s->K^*::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToKstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToKstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToKstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToKstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToKstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToKstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToKstar>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> phi
    // {{{
    ObservableGroup
    make_bs_to_phi_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \phi$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to \phi$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B_s->phi::V(q2)", R"(V^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::A_0(q2)", R"(A_0^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::A_1(q2)", R"(A_1^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::A_2(q2)", R"(A_2^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::A_12(q2)", R"(A_{12}^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::T_1(q2)", R"(T_1^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::T_2(q2)", R"(T_2^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::T_3(q2)", R"(T_3^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::T_23(q2)", R"(T_{23}^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_perp(q2)", R"(\mathcal{F}_\perp^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_perp, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_para(q2)", R"(\mathcal{F}_\parallel^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_para, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_long(q2)", R"(\mathcal{F}_0^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_long, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_perp_T(q2)", R"(\mathcal{F}_{\perp,T}^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_perp_T, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_para_T(q2)", R"(\mathcal{F}_{\parallel,T}^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_para_T, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_long_T(q2)", R"(\mathcal{F}_{0,T}^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_long_T, std::make_tuple("q2")),

                make_expression_observable("B_s->phi::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::V(q2)>> / <<B_s->phi::A_1(q2)>> )"),

                make_expression_observable("B_s->phi::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::A_2(q2)>> / <<B_s->phi::A_1(q2)>> )"),

                make_expression_observable("B_s->phi::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::A_12(q2)>> / <<B_s->phi::A_1(q2)>> )"),

                make_expression_observable("B_s->phi::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::T_23(q2)>> / <<B_s->phi::T_2(q2)>> )"),

                make_observable("B_s->phi::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToPhi>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToPhi>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToPhi>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToPhi>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToPhi>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToPhi>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToPhi>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),

                make_expression_observable("B_s->phi::F_perp_T(q2)/F_perp(q2)", R"(\mathcal{F}_{\perp,T}(q^2)/\mathcal{F}_\perp(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::F_perp_T(q2)>> / <<B_s->phi::F_perp(q2)>> )"),

                make_expression_observable("B_s->phi::F_para_T(q2)/F_para(q2)", R"(\mathcal{F}_{\parallel,T}(q^2)/\mathcal{F}_\parallel(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::F_para_T(q2)>> / <<B_s->phi::F_para(q2)>> )"),

                make_expression_observable("B_s->phi::F_long_T(q2)/F_long(q2)", R"(\mathcal{F}_{0,T}(q^2)/\mathcal{F}_0(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::F_long_T(q2)>> / <<B_s->phi::F_long(q2)>> )"),


                make_observable("B_s->phi::Saturation[0^+_V]", R"(\textrm{Saturation}[0^+_V])", Unit::None(),
                        &BFW2010FormFactors<BsToPhi, PToV>::saturation_0p_v),

                make_observable("B_s->phi::Saturation[0^-_A]", R"(\textrm{Saturation}[0^-_A])", Unit::None(),
                        &BFW2010FormFactors<BsToPhi, PToV>::saturation_0m_a),

                make_observable("B_s->phi::Saturation[1^-_V]", R"(\textrm{Saturation}[1^-_V])", Unit::None(),
                        &BFW2010FormFactors<BsToPhi, PToV>::saturation_1m_v),

                make_observable("B_s->phi::Saturation[1^+_A]", R"(\textrm{Saturation}[1^+_A])", Unit::None(),
                        &BFW2010FormFactors<BsToPhi, PToV>::saturation_1p_a),

                make_observable("B_s->phi::Saturation[1^-_T]", R"(\textrm{Saturation}[1^-_T])", Unit::None(),
                        &BFW2010FormFactors<BsToPhi, PToV>::saturation_1m_t),

                make_observable("B_s->phi::Saturation[1^+_T5]", R"(\textrm{Saturation}[1^+_{T5}])", Unit::None(),
                        &BFW2010FormFactors<BsToPhi, PToV>::saturation_1p_t5)
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> D_s^*
    // {{{
    ObservableGroup
    make_bs_to_dsstar_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \bar{D}_s^*$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to \bar{D}_s^*$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B_s->D_s^*::V(q2)", R"(V^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_0(q2)", R"(A_0^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_1(q2)", R"(A_1^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_2(q2)", R"(A_2^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_12(q2)", R"(A_{12}^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_1(q2)", R"(T_1^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_2(q2)", R"(T_2^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_3(q2)", R"(T_3^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_23(q2)", R"(T_{23}^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_observable("B_s->D_s^*::h_A1(q2)", R"(h_{A_1}^{\bar{B}_s\to D_s^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BsToDsstar, PToV>::h_a1,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::h_A2(q2)", R"(h_{A_2}^{\bar{B}_s\to D_s^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BsToDsstar, PToV>::h_a2,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::h_A3(q2)", R"(h_{A_3}^{\bar{B}_s\to D_s^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BsToDsstar, PToV>::h_a3,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::h_V(q2)", R"(h_V^{\bar{B}_s\to D_s^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BsToDsstar, PToV>::h_v,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::h_T1(q2)", R"(h_{T_1}^{\bar{B}_s\to D_s^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BsToDsstar, PToV>::h_t1,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::h_T2(q2)", R"(h_{T_2}^{\bar{B}_s\to D_s^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BsToDsstar, PToV>::h_t2,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::h_T3(q2)", R"(h_{T_3}^{\bar{B}_s\to D_s^*}(q^2))",
                        Unit::None(),
                        &HQETFormFactors<BsToDsstar, PToV>::h_t3,
                        std::make_tuple("q2")),

                make_expression_observable("B_s->D_s^*::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->D_s^*::V(q2)>> / <<B_s->D_s^*::A_1(q2)>> )"),

                make_expression_observable("B_s->D_s^*::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->D_s^*::A_2(q2)>> / <<B_s->D_s^*::A_1(q2)>> )"),

                make_expression_observable("B_s->D_s^*::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->D_s^*::A_12(q2)>> / <<B_s->D_s^*::A_1(q2)>> )"),

                make_expression_observable("B_s->D_s^*::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B_s->D_s^*::T_23(q2)>> / <<B_s->D_s^*::T_2(q2)>> )"),

                make_observable("B_s->D_s^*::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToDsstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToDsstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                         &AnalyticFormFactorBToVLCSR<BsToDsstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToDsstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToDsstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToDsstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<BsToDsstar>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::a_0[A_1]@HQE", R"(a_0^{A_1})",
                        Unit::None(),
                        &BGLCoefficients::A1s_a0),

                make_observable("B_s->D_s^*::a_1[A_1]@HQE", R"(a_1^{A_1})",
                        Unit::None(),
                        &BGLCoefficients::A1s_a1),

                make_observable("B_s->D_s^*::a_2[A_1]@HQE", R"(a_2^{A_1})",
                        Unit::None(),
                        &BGLCoefficients::A1s_a2),

                make_expression_observable("B_s->D_s^*::a_1/a_0[A_1]@HQE", R"(a_1^{A_1}/a_0^{A_1})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_1[A_1]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[A_1]@HQE>>
                        )"),

                make_expression_observable("B_s->D_s^*::a_2/a_0[A_1]@HQE", R"(a_2^{A_1}/a_0^{A_1})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_2[A_1]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[A_1]@HQE>>
                        )"),

                make_observable("B_s->D_s^*::a_0[A_5]@HQE", R"(a_0^{A_5})",
                        Unit::None(),
                        &BGLCoefficients::A5s_a0),

                make_observable("B_s->D_s^*::a_1[A_5]@HQE", R"(a_1^{A_5})",
                        Unit::None(),
                        &BGLCoefficients::A5s_a1),

                make_observable("B_s->D_s^*::a_2[A_5]@HQE", R"(a_2^{A_5})",
                        Unit::None(),
                        &BGLCoefficients::A5s_a2),

                make_expression_observable("B_s->D_s^*::a_1/a_0[A_5]@HQE", R"(a_1^{A_5}/a_0^{A_5})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_1[A_5]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[A_5]@HQE>>
                        )"),

                make_expression_observable("B_s->D_s^*::a_2/a_0[A_5]@HQE", R"(a_2^{A_5}/a_0^{A_5})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_2[A_5]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[A_5]@HQE>>
                        )"),

                make_observable("B_s->D_s^*::a_0[P_1]@HQE", R"(a_0^{P_1})",
                        Unit::None(),
                        &BGLCoefficients::P1s_a0),

                make_observable("B_s->D_s^*::a_1[P_1]@HQE", R"(a_1^{P_1})",
                        Unit::None(),
                        &BGLCoefficients::P1s_a1),

                make_observable("B_s->D_s^*::a_2[P_1]@HQE", R"(a_2^{P_1})",
                        Unit::None(),
                        &BGLCoefficients::P1s_a2),

                make_expression_observable("B_s->D_s^*::a_1/a_0[P_1]@HQE", R"(a_1^{P_1}/a_0^{P_1})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_1[P_1]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[P_1]@HQE>>
                        )"),

                make_expression_observable("B_s->D_s^*::a_2/a_0[P_1]@HQE", R"(a_2^{P_1}/a_0^{P_1})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_2[P_1]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[P_1]@HQE>>
                        )"),

                make_observable("B_s->D_s^*::a_0[V_4]@HQE", R"(a_0^{V_4})",
                        Unit::None(),
                        &BGLCoefficients::V4s_a0),

                make_observable("B_s->D_s^*::a_1[V_4]@HQE", R"(a_1^{V_4})",
                        Unit::None(),
                        &BGLCoefficients::V4s_a1),

                make_observable("B_s->D_s^*::a_2[V_4]@HQE", R"(a_2^{V_4})",
                        Unit::None(),
                        &BGLCoefficients::V4s_a2),

                make_expression_observable("B_s->D_s^*::a_1/a_0[V_4]@HQE", R"(a_1^{V_4}/a_0^{V_4})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_1[V_4]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[V_4]@HQE>>
                        )"),

                make_expression_observable("B_s->D_s^*::a_2/a_0[V_4]@HQE", R"(a_2^{V_4}/a_0^{V_4})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_2[V_4]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[V_4]@HQE>>
                        )"),

                make_observable("B_s->D_s^*::a_0[T_1]@HQE", R"(a_0^{T_1})",
                        Unit::None(),
                        &BGLCoefficients::T1s_a0),

                make_observable("B_s->D_s^*::a_1[T_1]@HQE", R"(a_1^{T_1})",
                        Unit::None(),
                        &BGLCoefficients::T1s_a1),

                make_observable("B_s->D_s^*::a_2[T_1]@HQE", R"(a_2^{T_1})",
                        Unit::None(),
                        &BGLCoefficients::T1s_a2),

                make_expression_observable("B_s->D_s^*::a_1/a_0[T_1]@HQE", R"(a_1^{T_1}/a_0^{T_1})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_1[T_1]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[T_1]@HQE>>
                        )"),

                make_expression_observable("B_s->D_s^*::a_2/a_0[T_1]@HQE", R"(a_2^{T_1}/a_0^{T_1})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_2[T_1]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[T_1]@HQE>>
                        )"),

                make_observable("B_s->D_s^*::a_0[T_2]@HQE", R"(a_0^{T_2})",
                        Unit::None(),
                        &BGLCoefficients::T2s_a0),

                make_observable("B_s->D_s^*::a_1[T_2]@HQE", R"(a_1^{T_2})",
                        Unit::None(),
                        &BGLCoefficients::T2s_a1),

                make_observable("B_s->D_s^*::a_2[T_2]@HQE", R"(a_2^{T_2})",
                        Unit::None(),
                        &BGLCoefficients::T2s_a2),

                make_expression_observable("B_s->D_s^*::a_1/a_0[T_2]@HQE", R"(a_1^{T_2}/a_0^{T_2})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_1[T_2]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[T_2]@HQE>>
                        )"),

                make_expression_observable("B_s->D_s^*::a_2/a_0[T_2]@HQE", R"(a_2^{T_2}/a_0^{T_2})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_2[T_2]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[T_2]@HQE>>
                        )"),

                make_observable("B_s->D_s^*::a_0[T_23]@HQE", R"(a_0^{T_{23}})",
                        Unit::None(),
                        &BGLCoefficients::T23s_a0),

                make_observable("B_s->D_s^*::a_1[T_23]@HQE", R"(a_1^{T_{23}})",
                        Unit::None(),
                        &BGLCoefficients::T23s_a1),

                make_observable("B_s->D_s^*::a_2[T_23]@HQE", R"(a_2^{T_{23}})",
                        Unit::None(),
                        &BGLCoefficients::T23s_a2),

                make_expression_observable("B_s->D_s^*::a_1/a_0[T_23]@HQE", R"(a_1^{T_{23}}/a_0^{T_{23}})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_1[T_23]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[T_23]@HQE>>
                        )"),

                make_expression_observable("B_s->D_s^*::a_2/a_0[T_23]@HQE", R"(a_2^{T_{23}}/a_0^{T_{23}})",
                        Unit::None(),
                        R"(
                        <<B_s->D_s^*::a_2[T_23]@HQE>>
                        /
                        <<B_s->D_s^*::a_0[T_23]@HQE>>
                        )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B -> P P
    // {{{

    // B -> pi pi
    // {{{
    ObservableGroup
    make_b_to_pi_pi_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B \to \pi \pi$ transitions)",
            R"(Pseudo observables representing the $B \to \pi \pi$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->pipi::Im{F_perp}(q2,k2,z)", R"(\mathrm{Im}\,F_\perp^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_perp, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{F_para}(q2,k2,z)", R"(\mathrm{Im}\,F_\parallel^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_para, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{F_long}(q2,k2,z)", R"(\mathrm{Im}\,F_0^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_long, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{F_time}(q2,k2,z)", R"(\mathrm{Im}\,F_t^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_time, std::make_tuple("q2", "k2", "z")),

                make_observable("B->pipi::Im{Res{F_perp}}(q2,k2)@FvDV2018A", R"(\mathrm{Res}\,\mathrm{Im}\,F_\perp^{B\to \pi\pi}(q^2,k^2,z))",
                        Unit::None(),
                        &AnalyticFormFactorBToPiPiFvDV2018::f_perp_im_res_qhat2,
                        std::make_tuple("q2", "k2")
                        ),

                make_observable("B->pipi::Im{Res{F_para}}(q2,k2)@FvDV2018A", R"(\mathrm{Res}\,\mathrm{Im}\,F_\parallel^{B\to \pi\pi}(q^2,k^2,z))",
                        Unit::None(),
                        &AnalyticFormFactorBToPiPiFvDV2018::f_para_im_res_qhat2,
                        std::make_tuple("q2", "k2")
                        ),

                make_observable("B->pipi::Im{Res{F_long}}(q2,k2)@FvDV2018A", R"(\mathrm{Res}\,\mathrm{Im}\,F_0^{B\to \pi\pi}(q^2,k^2,z))",
                        Unit::None(),
                        &AnalyticFormFactorBToPiPiFvDV2018::f_long_im_res_qhat2,
                        std::make_tuple("q2", "k2")
                        ),

                make_observable("B->pipi::Im{Res{F_time}}(q2,k2)@FvDV2018A", R"(\mathrm{Res}\,\mathrm{Im}\,F_t^{B\to \pi\pi}(q^2,k^2,z))",
                        Unit::None(),
                        &AnalyticFormFactorBToPiPiFvDV2018::f_time_im_res_qhat2,
                        std::make_tuple("q2", "k2")
                        ),

                make_observable("B->pipi::Im{Res{F_perp}}(q2,k2)@FvDV2018D", R"(\mathrm{Res}\,\mathrm{Im}\,F_\perp^{B\to \pi\pi}(q^2,k^2,z))",
                        Unit::None(),
                        &FvDV2018FormFactors<BToPiPi>::f_perp_im_res_qhat2,
                        std::make_tuple("q2", "k2")
                        ),

                make_observable("B->pipi::Im{Res{F_para}}(q2,k2)@FvDV2018D", R"(\mathrm{Res}\,\mathrm{Im}\,F_\parallel^{B\to \pi\pi}(q^2,k^2,z))",
                        Unit::None(),
                        &FvDV2018FormFactors<BToPiPi>::f_para_im_res_qhat2,
                        std::make_tuple("q2", "k2")
                        ),

                make_observable("B->pipi::Im{Res{F_long}}(q2,k2)@FvDV2018D", R"(\mathrm{Res}\,\mathrm{Im}\,F_0^{B\to \pi\pi}(q^2,k^2,z))",
                        Unit::None(),
                        &FvDV2018FormFactors<BToPiPi>::f_long_im_res_qhat2,
                        std::make_tuple("q2", "k2")
                        ),

                make_observable("B->pipi::Im{Res{F_time}}(q2,k2)@FvDV2018D", R"(\mathrm{Res}\,\mathrm{Im}\,F_t^{B\to \pi\pi}(q^2,k^2,z))",
                        Unit::None(),
                        &FvDV2018FormFactors<BToPiPi>::f_time_im_res_qhat2,
                        std::make_tuple("q2", "k2")
                        ),

                make_observable("B->pipi::Saturation[1^-_V]", R"(\textrm{Saturation}[1^-_V])", Unit::None(),
                        &HKVT2025FormFactors<BToPiPi, PToPP>::saturation_1m_v),

                make_observable("B->pipi::Saturation[0^-_A]", R"(\textrm{Saturation}[0^-_A])", Unit::None(),
                        &HKVT2025FormFactors<BToPiPi, PToPP>::saturation_0m_a),

                make_observable("B->pipi::Saturation[1^+_A]", R"(\textrm{Saturation}[1^+_A])", Unit::None(),
                        &HKVT2025FormFactors<BToPiPi, PToPP>::saturation_1p_a)
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // 1/2^+ -> 1/2^+
    // {{{

    // Lambda_b -> Lambda
    // {{{
    ObservableGroup
    make_lambdab_to_lambda_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $\Lambda_b \to \Lambda$ transitions)",
            R"(Pseudo observables representing the full basis of $\Lambda_b \to \Lambda$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("Lambda_b->Lambda::f_time^V(q2)", R"(f_t^{V,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_long^V(q2)", R"(f_0^{V,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_perp^V(q2)", R"(f_\perp^{V,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_time^A(q2)", R"(f_t^{A,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_long^A(q2)", R"(f_0^{A,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_perp^A(q2)", R"(f_\perp^{A,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_long^T(q2)", R"(f_0^{T,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_perp^T(q2)", R"(f_\perp^{T,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_long^T5(q2)", R"(f_0^{T5,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t5, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_perp^T5(q2)", R"(f_\perp^{T5,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t5, std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda::Saturation[0^+_V]", R"(\textrm{Saturation}[0^+_V])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::saturation_0p_v),

                make_observable("Lambda_b->Lambda::Saturation[1^-_V]", R"(\textrm{Saturation}[1^-_V])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::saturation_1m_v),

                make_observable("Lambda_b->Lambda::Saturation[0^-_A]", R"(\textrm{Saturation}[0^-_A])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::saturation_0m_a),

                make_observable("Lambda_b->Lambda::Saturation[1^+_A]", R"(\textrm{Saturation}[1^+_A])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::saturation_1p_a),

                make_observable("Lambda_b->Lambda::Saturation[1^-_T]", R"(\textrm{Saturation}[1^-_T])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::saturation_1m_t),

                make_observable("Lambda_b->Lambda::Saturation[1^+_T5]", R"(\textrm{Saturation}[1^+_{T_5}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::saturation_1p_t5),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Lambda_c -> Lambda
    // {{{
    ObservableGroup
    make_lambdac_to_lambda_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $\Lambda_c \to \Lambda$ transitions)",
            R"(Pseudo observables representing the full basis of $\Lambda_c \to \Lambda$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("Lambda_c->Lambda::f_time^V(q2)", R"(f_t^{V,\Lambda_c\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_c->Lambda::f_long^V(q2)", R"(f_0^{V,\Lambda_c\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_c->Lambda::f_perp^V(q2)", R"(f_\perp^{V,\Lambda_c\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_c->Lambda::f_time^A(q2)", R"(f_t^{A,\Lambda_c\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_c->Lambda::f_long^A(q2)", R"(f_0^{A,\Lambda_c\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_c->Lambda::f_perp^A(q2)", R"(f_\perp^{A,\Lambda_c\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_c->Lambda::f_long^T(q2)", R"(f_0^{T,\Lambda_c\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_c->Lambda::f_perp^T(q2)", R"(f_\perp^{T,\Lambda_c\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_c->Lambda::f_long^T5(q2)", R"(f_0^{T5,\Lambda_c\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t5, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_c->Lambda::f_perp^T5(q2)", R"(f_\perp^{T5,\Lambda_c\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t5, std::make_tuple("q2")),

                make_observable("Lambda_c->Lambda::Saturation[0^+_V]", R"(\textrm{Saturation}[0^+_V])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_0p_v),

                make_observable("Lambda_c->Lambda::Saturation[0^-_A]", R"(\textrm{Saturation}[0^-_A])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_0m_a),

                make_observable("Lambda_c->Lambda::Saturation[1^-_V,0]", R"(\textrm{Saturation}[1^-_{V,0}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1m_v_0),
                make_observable("Lambda_c->Lambda::Saturation[1^-_V,perp]", R"(\textrm{Saturation}[1^-_{V,\perp}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1m_v_perp),
                make_observable("Lambda_c->Lambda::Saturation[1^-_V,para]", R"(\textrm{Saturation}[1^-_{V,\parallel}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1m_v_para),
                make_observable("Lambda_c->Lambda::Saturation[1^-_V]", R"(\textrm{Saturation}[1^-_V])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1m_v),

                make_observable("Lambda_c->Lambda::Saturation[1^+_A,0]", R"(\textrm{Saturation}[1^+_{A,0}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1p_a_0),
                make_observable("Lambda_c->Lambda::Saturation[1^+_A,perp]", R"(\textrm{Saturation}[1^+_{A,\perp}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1p_a_perp),
                make_observable("Lambda_c->Lambda::Saturation[1^+_A,para]", R"(\textrm{Saturation}[1^+_{A,\parallel}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1p_a_para),
                make_observable("Lambda_c->Lambda::Saturation[1^+_A]", R"(\textrm{Saturation}[1^+_A])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1p_a),

                make_observable("Lambda_c->Lambda::Saturation[1^-_T,0]", R"(\textrm{Saturation}[1^-_{T,0}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1m_t_0),
                make_observable("Lambda_c->Lambda::Saturation[1^-_T,perp]", R"(\textrm{Saturation}[1^-_{T,\perp}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1m_t_perp),
                make_observable("Lambda_c->Lambda::Saturation[1^-_T,para]", R"(\textrm{Saturation}[1^-_{T,\parallel}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1m_t_para),
                make_observable("Lambda_c->Lambda::Saturation[1^-_T]", R"(\textrm{Saturation}[1^-_T])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1m_t),

                make_observable("Lambda_c->Lambda::Saturation[1^+_T5,0]", R"(\textrm{Saturation}[1^+_{T5,0}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1p_t5_0),
                make_observable("Lambda_c->Lambda::Saturation[1^+_T5,perp]", R"(\textrm{Saturation}[1^+_{T5,\perp}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1p_t5_perp),
                make_observable("Lambda_c->Lambda::Saturation[1^+_T5,para]", R"(\textrm{Saturation}[1^+_{T5,\parallel}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1p_t5_para),
                make_observable("Lambda_c->Lambda::Saturation[1^+_T5]", R"(\textrm{Saturation}[1^+_{T5}])", Unit::None(),
                        &BMRvD2022FormFactors<LambdaCToLambda>::saturation_1p_t5)
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Lambda_b -> Lambda_c
    // {{{
    ObservableGroup
    make_lambdab_to_lambdac_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $\Lambda_b \to \Lambda_c$ transitions)",
            R"(Pseudo observables representing the full basis of $\Lambda_b \to \Lambda_c$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("Lambda_b->Lambda_c::f_time^V(q2)", R"(f_t^{V,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_long^V(q2)", R"(f_0^{V,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_perp^V(q2)", R"(f_\perp^{V,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_time^A(q2)", R"(f_t^{A,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_long^A(q2)", R"(f_0^{A,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_perp^A(q2)", R"(f_\perp^{A,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_long^T(q2)", R"(f_0^{T,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_perp^T(q2)", R"(f_\perp^{T,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_long^T5(q2)", R"(f_0^{T5,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t5, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_perp^T5(q2)", R"(f_\perp^{T5,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t5, std::make_tuple("q2")),

                // Zero-Recoil Sum Rule for the Lambda_b -> Lambda_c Form Factors
                make_observable("Lambda_b->Lambda_c::F(1)",
                        Unit::None(),
                        &ZeroRecoilSumRule<LambdaBToC>::vector_current),

                make_observable("Lambda_b->Lambda_c::G(1)",
                        Unit::None(),
                        &ZeroRecoilSumRule<LambdaBToC>::axialvector_current),

                make_observable("Lambda_b->Lambda_c::F_inel(1)",
                        Unit::None(),
                        &ZeroRecoilSumRule<LambdaBToC>::vector_current_inel),

                make_observable("Lambda_b->Lambda_c::G_inel(1)",
                        Unit::None(),
                        &ZeroRecoilSumRule<LambdaBToC>::axialvector_current_inel),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // 1/2^+ -> 3/2^-
    // {{{

    // Lambda_b -> Lambda(1520)
    // {{{
    ObservableGroup
    make_lambdab_to_threehalf_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $\Lambda_b \to \Lambda^*$ transitions)",
            R"(Pseudo observables representing the full basis of $\Lambda_b \to \Lambda^*$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_time12^V(q2)", R"(f_t^{V,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_time12_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_long12^V(q2)", R"(f_0^{V,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_long12_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_perp12^V(q2)", R"(f_\perp^{V,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_perp12_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_perp32^V(q2)", R"(f_{\perp'}^{V,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_perp32_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_time12^A(q2)", R"(f_t^{A,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_time12_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_long12^A(q2)", R"(f_0^{A,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_long12_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_perp12^A(q2)", R"(f_\perp^{A,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_perp12_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_perp32^A(q2)", R"(f_{\perp'}^{A,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_perp32_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_long12^T(q2)", R"(f_0^{T,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_long12_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_perp12^T(q2)", R"(f_\perp^{T,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_perp12_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_perp32^T(q2)", R"(f_{\perp'}^{T,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_perp32_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_long12^T5(q2)", R"(f_0^{T5,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_long12_t5, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_perp12^T5(q2)", R"(f_\perp^{T5,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_perp12_t5, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda(1520)::f_perp32^T5(q2)", R"(f_{\perp'}^{T5,\Lambda_b\to\Lambda(1520)}(q^2))",
                        &FormFactors<OneHalfPlusToThreeHalfMinus>::f_perp32_t5, std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda(1520)::Saturation[0^+_V]", R"(\textrm{Saturation}[0^+_V])", Unit::None(),
                        &ABR2022FormFactors<LambdaBToLambda1520>::saturation_0p_v),

                make_observable("Lambda_b->Lambda(1520)::Saturation[1^-_V]", R"(\textrm{Saturation}[1^-_V])", Unit::None(),
                        &ABR2022FormFactors<LambdaBToLambda1520>::saturation_1m_v),

                make_observable("Lambda_b->Lambda(1520)::Saturation[0^-_A]", R"(\textrm{Saturation}[0^-_A])", Unit::None(),
                        &ABR2022FormFactors<LambdaBToLambda1520>::saturation_0m_a),

                make_observable("Lambda_b->Lambda(1520)::Saturation[1^+_A]", R"(\textrm{Saturation}[1^+_A])", Unit::None(),
                        &ABR2022FormFactors<LambdaBToLambda1520>::saturation_1p_a),

                make_observable("Lambda_b->Lambda(1520)::Saturation[1^-_T]", R"(\textrm{Saturation}[1^-_T])", Unit::None(),
                        &ABR2022FormFactors<LambdaBToLambda1520>::saturation_1m_t),

                make_observable("Lambda_b->Lambda(1520)::Saturation[1^+_T5]", R"(\textrm{Saturation}[1^+_{T_5}])", Unit::None(),
                        &ABR2022FormFactors<LambdaBToLambda1520>::saturation_1p_t5),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}


    // unitarity bounds
    // {{{

    ObservableGroup
    make_unitarity_bounds_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Unitarity Bounds)",
            R"(Pseudo observables arising in the various unitarity bounds of semileptonic form factors.)",
            {
                make_observable("b->c::Bound[0^+]@CLN", R"(B^{b\to c}_{0^+})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_0p),

                make_observable("b->c::Bound[0^-]@CLN", R"(B^{b\to c}_{0^-})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_0m),

                make_observable("b->c::Bound[1^+]@CLN", R"(B^{b\to c}_{1^+})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_1p),

                make_observable("b->c::Bound[1^-]@CLN", R"(B^{b\to c}_{1^-})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_1m),

                make_observable("b->c::Bound[1^+,T]@CLN", R"(B^{b\to c}_{1^+,T})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_1p_T),

                make_observable("b->c::Bound[1^-,T]@CLN", R"(B^{b\to c}_{1^-,T})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_1m_T),

                make_observable("b->c::Bound[0^+]@OPE", R"(B^{b\to c}_{0^+})",
                        Unit::None(),
                        &OPEUnitarityBounds::bound_0p),

                make_observable("b->c::Bound[0^-]@OPE", R"(B^{b\to c}_{0^-})",
                        Unit::None(),
                        &OPEUnitarityBounds::bound_0m),

                make_observable("b->c::Bound[1^+]@OPE", R"(B^{b\to c}_{1^+})",
                        Unit::None(),
                        &OPEUnitarityBounds::bound_1p),

                make_observable("b->c::Bound[1^-]@OPE", R"(B^{b\to c}_{1^-})",
                        Unit::None(),
                        &OPEUnitarityBounds::bound_1m),

                make_observable("b->c::Bound[0^+]@BGL", R"(B^{b\to c}_{0^+})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_0p),

                make_observable("b->c::Bound[0^-]@BGL", R"(B^{b\to c}_{0^-})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_0m),

                make_observable("b->c::Bound[1^+]@BGL", R"(B^{b\to c}_{1^+})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_1p),

                make_observable("b->c::Bound[1^-]@BGL", R"(B^{b\to c}_{1^-})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_1m),

                // cf. [BMRvD:2021A] eq. (31-33)
                // sb states
                make_expression_observable("B_s0::Saturation[0^+_V]", R"(\textrm{Saturation}_{B_{s,0}}[0^+_V])",
                        Unit::None(),
                        R"(<<decay-constant::B_s,0>>^2 / <<mass::B_s,0>>^2 / <<b->s::chiOPE[0^+_V]>>)"),

                make_expression_observable("B_s::Saturation[0^-_A]", R"(\textrm{Saturation}_{B_s^0}[0^-_A])",
                        Unit::None(),
                        R"(<<decay-constant::B_s>>^2 / <<mass::B_s>>^2 / <<b->s::chiOPE[0^-_A]>>)"),

                make_expression_observable("B_s^*::Saturation[1^-_V,0]", R"(\textrm{Saturation}_{B_s^*}[1^-_{V,0}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s^*>>^2 / <<mass::B_s^*>>^4 / <<b->s::chiOPE[1^-_V]>> / 3.0)"),

                make_expression_observable("B_s^*::Saturation[1^-_V,perp]", R"(\textrm{Saturation}_{B_s^*}[1^-_{V,\perp}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s^*>>^2 / <<mass::B_s^*>>^4 / <<b->s::chiOPE[1^-_V]>> / 3.0)"),

                make_expression_observable("B_s^*::Saturation[1^-_V.para]", R"(\textrm{Saturation}_{B_s^*}[1^-_{V,\parallel}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s^*>>^2 / <<mass::B_s^*>>^4 / <<b->s::chiOPE[1^-_V]>> / 3.0)"),

                make_expression_observable("B_s^*::Saturation[1^-_V]", R"(\textrm{Saturation}_{B_s^*}[1^-_V])",
                        Unit::None(),
                        R"(<<decay-constant::B_s^*>>^2 / <<mass::B_s^*>>^4 / <<b->s::chiOPE[1^-_V]>>)"),

                make_expression_observable("B_s1::Saturation[1^+_A,0]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{A,0}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s,1>>^2 / <<mass::B_s,1>>^4 / <<b->s::chiOPE[1^+_A]>> / 3.0)"),

                make_expression_observable("B_s1::Saturation[1^+_A,perp]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{A,\perp}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s,1>>^2 / <<mass::B_s,1>>^4 / <<b->s::chiOPE[1^+_A]>> / 3.0)"),

                make_expression_observable("B_s1::Saturation[1^+_A,para]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{A,\parallel}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s,1>>^2 / <<mass::B_s,1>>^4 / <<b->s::chiOPE[1^+_A]>> / 3.0)"),

                make_expression_observable("B_s1::Saturation[1^+_A]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_A])",
                        Unit::None(),
                        R"(<<decay-constant::B_s,1>>^2 / <<mass::B_s,1>>^4 / <<b->s::chiOPE[1^+_A]>>)"),

                make_expression_observable("B_s^*::Saturation[1^-_T,0]", R"(\textrm{Saturation}_{B_s^*}[1^-_{T,0}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s^*,T>>^2 / <<mass::B_s^*>>^4 / <<b->s::chiOPE[1^-_T]>> / 3.0)"),

                make_expression_observable("B_s^*::Saturation[1^-_T,perp]", R"(\textrm{Saturation}_{B_s^*}[1^-_{T,\perp}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s^*,T>>^2 / <<mass::B_s^*>>^4 / <<b->s::chiOPE[1^-_T]>> / 3.0)"),

                make_expression_observable("B_s^*::Saturation[1^-_T,para]", R"(\textrm{Saturation}_{B_s^*}[1^-_{T,\parallel}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s^*,T>>^2 / <<mass::B_s^*>>^4 / <<b->s::chiOPE[1^-_T]>> / 3.0)"),

                make_expression_observable("B_s^*::Saturation[1^-_T]", R"(\textrm{Saturation}_{B_s^*}[1^-_T])",
                        Unit::None(),
                        R"(<<decay-constant::B_s^*,T>>^2 / <<mass::B_s^*>>^4 / <<b->s::chiOPE[1^-_T]>>)"),

                make_expression_observable("B_s1::Saturation[1^+_T5,0]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{T_5,0}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s,1^T>>^2 / <<mass::B_s,1>>^4 / <<b->s::chiOPE[1^+_T5]>> / 3.0)"),

                make_expression_observable("B_s1::Saturation[1^+_T5,perp]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{T_5,\perp}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s,1^T>>^2 / <<mass::B_s,1>>^4 / <<b->s::chiOPE[1^+_T5]>> / 3.0)"),

                make_expression_observable("B_s1::Saturation[1^+_T5,para]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{T_5,\parallel}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s,1^T>>^2 / <<mass::B_s,1>>^4 / <<b->s::chiOPE[1^+_T5]>> / 3.0)"),

                make_expression_observable("B_s1::Saturation[1^+_T5]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{T_5}])",
                        Unit::None(),
                        R"(<<decay-constant::B_s,1^T>>^2 / <<mass::B_s,1>>^4 / <<b->s::chiOPE[1^+_T5]>>)"),

                // cs states
                make_expression_observable("D_s0::Saturation[0^+_V]", R"(\textrm{Saturation}_{B_{s,0}}[0^+_V])",
                        Unit::None(),
                        R"(<<decay-constant::D_s,0>>^2 / <<mass::D_s,0>>^2 / <<c->s::chiOPE[0^+_V]>>)"),

                make_expression_observable("D_s::Saturation[0^-_A]", R"(\textrm{Saturation}_{D_s^0}[0^-_A])",
                        Unit::None(),
                        R"(<<decay-constant::D_s>>^2 / <<mass::D_s>>^2 / <<c->s::chiOPE[0^-_A]>>)"),

                make_expression_observable("D_s^*::Saturation[1^-_V,0]", R"(\textrm{Saturation}_{D_s^*}[1^-_{V,0}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s^*>>^2 / <<mass::D_s^*>>^4 / <<c->s::chiOPE[1^-_V]>> / 3.0)"),

                make_expression_observable("D_s^*::Saturation[1^-_V,perp]", R"(\textrm{Saturation}_{D_s^*}[1^-_{V,\perp}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s^*>>^2 / <<mass::D_s^*>>^4 / <<c->s::chiOPE[1^-_V]>> / 3.0)"),

                make_expression_observable("D_s^*::Saturation[1^-_V.para]", R"(\textrm{Saturation}_{D_s^*}[1^-_{V,\parallel}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s^*>>^2 / <<mass::D_s^*>>^4 / <<c->s::chiOPE[1^-_V]>> / 3.0)"),

                make_expression_observable("D_s^*::Saturation[1^-_V]", R"(\textrm{Saturation}_{D_s^*}[1^-_V])",
                        Unit::None(),
                        R"(<<decay-constant::D_s^*>>^2 / <<mass::D_s^*>>^4 / <<c->s::chiOPE[1^-_V]>>)"),

                make_expression_observable("D_s1::Saturation[1^+_A,0]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{A,0}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s,1>>^2 / <<mass::D_s,1>>^4 / <<c->s::chiOPE[1^+_A]>> / 3.0)"),

                make_expression_observable("D_s1::Saturation[1^+_A,perp]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{A,\perp}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s,1>>^2 / <<mass::D_s,1>>^4 / <<c->s::chiOPE[1^+_A]>> / 3.0)"),

                make_expression_observable("D_s1::Saturation[1^+_A,para]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{A,\parallel}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s,1>>^2 / <<mass::D_s,1>>^4 / <<c->s::chiOPE[1^+_A]>> / 3.0)"),

                make_expression_observable("D_s1::Saturation[1^+_A]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_A])",
                        Unit::None(),
                        R"(<<decay-constant::D_s,1>>^2 / <<mass::D_s,1>>^4 / <<c->s::chiOPE[1^+_A]>>)"),

                make_expression_observable("D_s^*::Saturation[1^-_T,0]", R"(\textrm{Saturation}_{D_s^*}[1^-_{T,0}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s^*,T>>^2 / <<mass::D_s^*>>^4 / <<c->s::chiOPE[1^-_T]>> / 3.0)"),

                make_expression_observable("D_s^*::Saturation[1^-_T,perp]", R"(\textrm{Saturation}_{D_s^*}[1^-_{T,\perp}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s^*,T>>^2 / <<mass::D_s^*>>^4 / <<c->s::chiOPE[1^-_T]>> / 3.0)"),

                make_expression_observable("D_s^*::Saturation[1^-_T,para]", R"(\textrm{Saturation}_{D_s^*}[1^-_{T,\parallel}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s^*,T>>^2 / <<mass::D_s^*>>^4 / <<c->s::chiOPE[1^-_T]>> / 3.0)"),

                make_expression_observable("D_s^*::Saturation[1^-_T]", R"(\textrm{Saturation}_{D_s^*}[1^-_T])",
                        Unit::None(),
                        R"(<<decay-constant::D_s^*,T>>^2 / <<mass::D_s^*>>^4 / <<c->s::chiOPE[1^-_T]>>)"),

                make_expression_observable("D_s1::Saturation[1^+_T5,0]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{T_5,0}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s,1^T>>^2 / <<mass::D_s,1>>^4 / <<c->s::chiOPE[1^+_T5]>> / 3.0)"),

                make_expression_observable("D_s1::Saturation[1^+_T5,perp]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{T_5,\perp}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s,1^T>>^2 / <<mass::D_s,1>>^4 / <<c->s::chiOPE[1^+_T5]>> / 3.0)"),

                make_expression_observable("D_s1::Saturation[1^+_T5,para]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{T_5,\parallel}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s,1^T>>^2 / <<mass::D_s,1>>^4 / <<c->s::chiOPE[1^+_T5]>> / 3.0)"),

                make_expression_observable("D_s1::Saturation[1^+_T5]", R"(\textrm{Saturation}_{B_{s,1}}[1^+_{T_5}])",
                        Unit::None(),
                        R"(<<decay-constant::D_s,1^T>>^2 / <<mass::D_s,1>>^4 / <<c->s::chiOPE[1^+_T5]>>)"),

                // ub states
                make_expression_observable("B_u^*::Saturation[1^-_V]", R"(\textrm{Saturation}_{B_u^*}[1^-_{V}])",
                        Unit::None(),
                        R"(<<decay-constant::B_u^*>>^2 / <<mass::B_u^*>>^4 / <<b->u::chiOPE[1^-_V]>>)"),

                make_expression_observable("B_u::Saturation[0^-_A]", R"(\textrm{Saturation}_{B_u}[0^-_A])",
                        Unit::None(),
                        R"(<<decay-constant::B_u>>^2 / <<mass::B_u>>^2 / <<b->u::chiOPE[0^-_A]>>)")
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B-meson LCDAs
    // {{{
    ObservableGroup
    make_b_meson_lcdas_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"($B$-meson LCDAs)",
            R"(Pseudo observables arising in the description of $B$-meson Light-Cone Distribution Amplitudes (LCDAs).)",
            {
                make_observable("B::phitilde_+(-i*tau,mu)@FLvD2022", R"(\tilde\phi_{B,+}(-i \tau, \mu))",
                        Unit::None(),
                        &heavy_meson_lcdas::FLvD2022::phitilde_plus,
                        std::make_tuple("tau", "mu"),
                        Options{ { "Q"_ok, "b" }, { "q"_ok, "u" } }),
                make_observable("B::tau*d_dtau_phitilde_+(-i*tau,mu)@FLvD2022", R"(-i \tau \, \tilde\phi^{\prime}_{B,+}(-i \tau, \mu))",
                        Unit::None(),
                        &heavy_meson_lcdas::FLvD2022::t_d_dt_phitilde_plus,
                        std::make_tuple("tau", "mu"),
                        Options{ { "Q"_ok, "b" }, { "q"_ok, "u" } }),
                make_observable("B::tau^2*d2_d2tau_phitilde_+(-i*tau,mu)@FLvD2022", R"(-\tau^2 \, \tilde\phi^{\prime\prime}_{B,+}(-i \tau, \mu))",
                        Unit::None(),
                        &heavy_meson_lcdas::FLvD2022::t2_d2_d2t_phitilde_plus,
                        std::make_tuple("tau", "mu"),
                        Options{ { "Q"_ok, "b" }, { "q"_ok, "u" } }),
                make_observable("B::L0@FLvD2022", R"(L_0)",
                        Unit::InverseGeV(),
                        &heavy_meson_lcdas::FLvD2022::inverse_moment,
                        std::make_tuple("mu"),
                        Options{ { "Q"_ok, "b" }, { "q"_ok, "u" } }),
                make_observable("B::L1@FLvD2022", R"(L_1)",
                        Unit::InverseGeV(),
                        &heavy_meson_lcdas::FLvD2022::logarithmic_moment_1,
                        std::make_tuple("mu"),
                        Options{ { "Q"_ok, "b" }, { "q"_ok, "u" } }),
                make_observable("B::L2@FLvD2022", R"(L_2)",
                        Unit::InverseGeV(),
                        &heavy_meson_lcdas::FLvD2022::logarithmic_moment_2,
                        std::make_tuple("mu"),
                        Options{ { "Q"_ok, "b" }, { "q"_ok, "u" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // }}}
    // D -> P(seudoscalar)

    // D -> K
    // {{{
    ObservableGroup
    make_d_to_k_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $D \to K$ transitions)",
            R"(Pseudo observables representing the full basis of $D \to K$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("D->K::f_+(q2)", R"(f_+^{D \to K}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("D->K::f_0(q2)", R"(f_0^{D \to K}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("D->K::f_T(q2)", R"(f_T^{D \to K}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("D->K::f_-(q2)", R"(f_-^{D \to K}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                make_form_factor_adapter("D->K::F_plus(q2)", R"(F_+^{D \to K}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("D->K::F_plus_T(q2)", R"(F_T^{D \to K}(q^2))",
                        &FormFactors<PToP>::f_plus_T, std::make_tuple("q2")),

                make_expression_observable("D->K::F_T(q2)/F_plus(q2)", R"(F_T(q^2)/F_+(q^2))",
                        Unit::None(),
                        R"( <<D->K::f_T(q2)>> / <<D->K::F_plus(q2)>> )"),

                make_expression_observable("D->K::F_plus_T(q2)/F_plus(q2)", R"(F_{+,T}(q^2)/F_+(q^2))",
                        Unit::None(),
                        R"( <<D->K::F_plus_T(q2)>> / <<D->K::F_plus(q2)>> )"),

                make_observable("D->K::Saturation[0^+_V]", R"(\textrm{Saturation}[0^+_V])", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::saturation_0p_v),

                make_observable("D->K::Saturation[0^-_A]", R"(\textrm{Saturation}[0^-_A])", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::saturation_0m_a),

                make_observable("D->K::Saturation[1^-_V]", R"(\textrm{Saturation}[1^-_V])", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::saturation_1m_v),

                make_observable("D->K::Saturation[1^+_A]", R"(\textrm{Saturation}[1^+_A])", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::saturation_1p_a),

                make_observable("D->K::Saturation[1^-_T]", R"(\textrm{Saturation}[1^-_T])", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::saturation_1m_t),

                make_observable("D->K::Saturation[1^+_T5]", R"(\textrm{Saturation}[1^+_{T5}])", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::saturation_1p_t5),

                // Auxiliary functions for [BFW:2010A]
                make_observable("D->K::f_+_series(q2)@BFW2010", R"(\hat{f}_+^{D \to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::f_p_series, std::make_tuple("q2")),

                make_observable("D->K::f_+_series_prime(q2)@BFW2010", R"(\hat{f}_+^{\prime D \to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::f_p_series_prime, std::make_tuple("q2")),

                make_observable("D->K::f_0_series(q2)@BFW2010", R"(\hat{f}_0^{D \to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::f_0_series, std::make_tuple("q2")),

                make_observable("D->K::f_0_series_prime(q2)@BFW2010", R"(\hat{f}_0^{\prime D \to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::f_0_series_prime, std::make_tuple("q2")),

                make_observable("D->K::f_T_series(q2)@BFW2010", R"(\hat{f}_T^{D \to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::f_t_series, std::make_tuple("q2")),

                make_observable("D->K::f_T_series_prime(q2)@BFW2010", R"(\hat{f}_T^{\prime D \to K}(q^2))", Unit::None(),
                        &BFW2010FormFactors<DToK, PToP>::f_t_series_prime, std::make_tuple("q2"))
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // D -> eta
    // {{{
    ObservableGroup
    make_d_to_eta_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $D\to \eta^{(\prime)}$ transitions)",
            R"(Pseudo observables representing the full basis of $D\to \eta^{(\prime)}$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("D->eta::f_+(q2)", R"(f_+^{D\to\eta}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("D->eta::f_T(q2)", R"(f_T^{D\to\eta}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("D->eta::f_0(q2)", R"(f_0^{D\to\eta}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("D->eta_prime::f_+(q2)", R"(f_+^{D\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("D->eta_prime::f_T(q2)", R"(f_T^{D\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("D->eta_prime::f_0(q2)", R"(f_0^{D\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // D_s -> eta
    // {{{
    ObservableGroup
    make_ds_to_eta_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $D_s\to \eta^{(\prime)}$ transitions)",
            R"(Pseudo observables representing the full basis of $D_s\to \eta^{(\prime)}$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("D_s->eta::f_+(q2)", R"(f_+^{D_s\to\eta}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("D_s->eta::f_T(q2)", R"(f_T^{D_s\to\eta}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("D_s->eta::f_0(q2)", R"(f_0^{D_s\to\eta}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("D_s->eta_prime::f_+(q2)", R"(f_+^{D_s\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("D_s->eta_prime::f_T(q2)", R"(f_T^{D_s\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("D_s->eta_prime::f_0(q2)", R"(f_0^{D_s\to\eta'}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // 0 -> PP
    // {{{

    // 0 -> pi pi
    // {{{
    ObservableGroup
    make_vacuum_to_pipi_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>{
            R"(Form factors for $0 \to \pi \pi$ transitions)",
            R"(Pseudo observables representing the full basis of $0 \to \pi \pi$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("0->pipi::Abs{f_+}^2(q2)", R"(|f_+^{0\to \pi\pi}(q^2)|^2)",
                        &FormFactors<VacuumToPP>::abs2_f_p, std::make_tuple("q2")),

                make_form_factor_adapter("0->pipi::Arg{f_+}(q2)", R"(\textrm{arg}(f_+^{0\to\pi\pi}(q^2)))",
                        &FormFactors<VacuumToPP>::arg_f_p, std::make_tuple("q2")),

                make_form_factor_adapter("0->pipi::Re{f_+}(Re{q2},Im{q2})", R"(\textrm{Re}(f_+^{0\to\pi\pi}(q^2)))",
                        &FormFactors<VacuumToPP>::re_f_p, std::make_tuple("Re{q2}", "Im{q2}")),

                make_form_factor_adapter("0->pipi::Im{f_+}(Re{q2},Im{q2})", R"(\textrm{Im}(f_+^{0\to\pi\pi}(q^2)))",
                        &FormFactors<VacuumToPP>::im_f_p, std::make_tuple("Re{q2}", "Im{q2}")),

                make_observable("0->pipi::b_0@KKRvD2024", R"(b_0^{0 \to \pi\pi})", Unit::None(),
                        &KKRvD2024FormFactors<VacuumToPiPi>::b_0),

                make_observable("0->pipi::b_1@KKRvD2024", R"(b_1^{0 \to \pi\pi})", Unit::None(),
                        &KKRvD2024FormFactors<VacuumToPiPi>::b_1),

                make_observable("0->pipi::Re{Res_z{f_+}}@KKRvD2024", R"(\textrm{Re}\,\textrm{Res}_{z} (M_\rho^2)\,f_+^{0\to\pi\pi})", Unit::None(),
                        &KKRvD2024FormFactors<VacuumToPiPi>::re_residue_rho),

                make_observable("0->pipi::Im{Res_z{f_+}}@KKRvD2024", R"(\textrm{Im}\,\textrm{Res}_{z} (M_\rho^2)\,f_+^{0\to\pi\pi})", Unit::None(),
                        &KKRvD2024FormFactors<VacuumToPiPi>::im_residue_rho),

                make_observable("0->pipi::Re{Res_q2{f_+}}@KKRvD2024", R"(\textrm{Re}\,\textrm{Res}_{q^2} (M_\rho^2)\,f_+^{0\to\pi\pi})", Unit::None(),
                        &KKRvD2024FormFactors<VacuumToPiPi>::re_residue_rho_q2),

                make_observable("0->pipi::Im{Res_q2{f_+}}@KKRvD2024", R"(\textrm{Im}\,\textrm{Res}_{q^2} (M_\rho^2)\,f_+^{0\to\pi\pi})", Unit::None(),
                        &KKRvD2024FormFactors<VacuumToPiPi>::im_residue_rho_q2),

                make_observable("0->pipi::r_pi^2@KKRvD2024", R"(\langle r_\pi^2 \rangle)", Unit::Femtometer2(),
                        &KKRvD2024FormFactors<VacuumToPiPi>::r_pi_squared),

                make_observable("0->pipi::Saturation@KKRvD2024", R"(\textrm{Saturation})", Unit::None(),
                        &KKRvD2024FormFactors<VacuumToPiPi>::saturation)
            }
        };

        return ObservableGroup(imp);
    }
    // }}}

    // 0 -> K pi
    // {{{
    ObservableGroup
    make_vacuum_to_Kpi_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>{
            R"(Form factors for $0 \to K \pi$ transitions)",
            R"(Pseudo observables representing the full basis of $0 \to K \pi$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("0->Kpi::Abs{f_+}^2(q2)", R"(|f_+^{0\to K\pi}(q^2)|^2)",
                        &FormFactors<VacuumToPP>::abs2_f_p, std::make_tuple("q2")),

                make_form_factor_adapter("0->Kpi::Arg{f_+}(q2)", R"(\textrm{arg}(f_+^{0\to K\pi}(q^2)))",
                        &FormFactors<VacuumToPP>::arg_f_p, std::make_tuple("q2")),

                make_form_factor_adapter("0->Kpi::Re{f_+}(Re{q2},Im{q2})", R"(\textrm{Re}(f_+^{0\to K\pi}(q^2)))",
                        &FormFactors<VacuumToPP>::re_f_p, std::make_tuple("Re{q2}", "Im{q2}")),

                make_form_factor_adapter("0->Kpi::Im{f_+}(Re{q2},Im{q2})", R"(\textrm{Im}(f_+^{0\to K\pi}(q^2)))",
                        &FormFactors<VacuumToPP>::im_f_p, std::make_tuple("Re{q2}", "Im{q2}")),

                make_observable("0->Kpi::f_+(0)", R"(f_+^{0 \to K\pi}(0))", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::fp_at_0),

                make_observable("0->Kpi::lambda_prime_+", R"(\lambda_+^{('), K \to \pi})", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::lambda_prime_plus),

                make_observable("0->Kpi::lambda_doubleprime_+", R"(\lambda^{\prime\prime, 0 \to K \pi}_+)", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::lambda_doubleprime_plus),

                make_observable("0->Kpi::b0_f+@KSvD2025", R"(b_0^{+, 0 \to K\pi})", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::b0_fp),

                make_observable("0->Kpi::Saturation_f+@KSvD2025", R"(\textrm{Saturation})", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::saturation_p),

                make_form_factor_adapter("0->Kpi::Abs{f_0}^2(q2)", R"(|f_0^{0\to K\pi}(q^2)|^2)",
                        &FormFactors<VacuumToPP>::abs2_f_0, std::make_tuple("q2")),

                make_form_factor_adapter("0->Kpi::Arg{f_0}(q2)", R"(\textrm{arg}(f_0^{0\to K\pi}(q^2)))",
                        &FormFactors<VacuumToPP>::arg_f_0, std::make_tuple("q2")),

                make_form_factor_adapter("0->Kpi::Re{f_0}(Re{q2},Im{q2})", R"(\textrm{Re}(f_0^{0\to K\pi}(q^2)))",
                        &FormFactors<VacuumToPP>::re_f_0, std::make_tuple("Re{q2}", "Im{q2}")),

                make_form_factor_adapter("0->Kpi::Im{f_0}(Re{q2},Im{q2})", R"(\textrm{Im}(f_0^{0\to K\pi}(q^2)))",
                        &FormFactors<VacuumToPP>::im_f_0, std::make_tuple("Re{q2}", "Im{q2}")),

                make_observable("0->Kpi::lambda_prime_0", R"(\lambda_0^{('), K \to \pi})", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::lambda_prime_zero),

                make_observable("0->Kpi::lambda_doubleprime_0", R"(\lambda^{\prime\prime}_0^{K \to \pi})", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::lambda_doubleprime_zero),

                make_observable("0->Kpi::b0_f0@KSvD2025", R"(b_0^{0, 0 \to K\pi})", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::b0_f0),

                make_observable("0->Kpi::Delta_CT@KSvD2025", R"(\Delta_{\textrm{CT}})", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::Delta_CT),

                make_observable("0->Kpi::Delta_CTtilde@KSvD2025", R"(\Delta_{\widetilde{\textrm{CT}}})", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::Delta_CTtilde),

                make_observable("0->Kpi::Saturation_f0@KSvD2025", R"(\textrm{Saturation})", Unit::None(),
                        &KSvD2025FormFactors<VacuumToKPi>::saturation_z)
            }
        };

        return ObservableGroup(imp);
    }
    // }}}

    // }}}



    ObservableSection
    make_form_factors_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Form factors",
            "",
            {
                // B -> P
                make_b_to_pi_form_factors_group(),
                make_b_to_eta_form_factors_group(),
                make_b_to_k_form_factors_group(),
                make_b_to_d_form_factors_group(),

                // B_s -> P
                make_bs_to_eta_form_factors_group(),
                make_bs_to_k_form_factors_group(),
                make_bs_to_ds_form_factors_group(),

                // B -> gamma
                make_b_to_gammastar_form_factors_group(),

                // B -> V
                make_b_to_gamma_form_factors_group(),
                make_b_to_omega_form_factors_group(),
                make_b_to_rho_form_factors_group(),
                make_b_to_kstar_form_factors_group(),
                make_b_to_dstar_form_factors_group(),

                // B_s -> V
                make_bs_to_kstar_form_factors_group(),
                make_bs_to_phi_form_factors_group(),
                make_bs_to_dsstar_form_factors_group(),

                // B -> P P
                make_b_to_pi_pi_form_factors_group(),

                // Lb -> 1/2^+
                make_lambdab_to_lambda_form_factors_group(),
                make_lambdab_to_lambdac_form_factors_group(),
                make_lambdac_to_lambda_form_factors_group(),

                // Lb -> 3/2^-
                make_lambdab_to_threehalf_form_factors_group(),

                // unitarity bounds
                make_unitarity_bounds_group(),

                // B-meson LCDAs
                make_b_meson_lcdas_group(),

                // D -> P
                make_d_to_k_form_factors_group(),
                make_d_to_eta_form_factors_group(),

                // D_s -> P
                make_ds_to_eta_form_factors_group(),

                // 0 -> PP
                make_vacuum_to_pipi_form_factors_group(),
                make_vacuum_to_Kpi_form_factors_group(),
            }
        );

        return ObservableSection(imp);
    }



}
