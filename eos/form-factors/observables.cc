/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2019 Danny van Dyk
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
#include <eos/form-factors/analytic-b-to-pi-pi.hh>
#include <eos/form-factors/analytic-b-to-v-lcsr.hh>
#include <eos/form-factors/analytic-b-to-p-lcsr.hh>
#include <eos/form-factors/b-lcdas.hh>
#include <eos/form-factors/baryonic-impl.hh>
#include <eos/form-factors/mesonic-impl.hh>
#include <eos/form-factors/observables.hh>
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

        return std::make_pair(qn, std::make_shared<FormFactorAdapterEntry<Transition_, Args_ ...>>(qn, latex, pp, function, kinematics_names));
    }

    template <typename Transition_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_form_factor_adapter(const char * name,
            double (FormFactors<Transition_>::* _function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        QualifiedName qn(name);
        qnp::Prefix pp = qn.prefix_part();
        std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> function(_function);

        return std::make_pair(qn, std::make_shared<FormFactorAdapterEntry<Transition_, Args_ ...>>(qn, "", pp, function, kinematics_names));
    }

    template <typename Transition_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_form_factor_adapter(const char * name,
            double (FormFactors<Transition_>::* _numerator)(const Args_ & ...) const,
            double (FormFactors<Transition_>::* _denominator)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        QualifiedName qn(name);
        std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> numerator(_numerator);
        std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> denominator(_denominator);

        return std::make_pair(qn, std::make_shared<FormFactorRatioAdapterEntry<Transition_, Args_ ...>>(qn, "", numerator, denominator, kinematics_names));
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

                // auxiliary variables, e.g. for determining the B-LCSR threshold parameters
                make_observable("B->pi::f_+[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BToPi>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B->pi::f_0[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BToPi>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B->pi::f_T[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BToPi>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),

                // auxiliary variables, e.g. for determining the pi-LCSR/SVZ threshold parameters
                make_observable("B->pi::M_B(LCSR)@DKMMO2008",
                        &AnalyticFormFactorBToPiDKMMO2008::MB_lcsr,
                        std::make_tuple("s")),

                make_observable("B->pi::M_B(SVZ)@DKMMO2008",
                        &AnalyticFormFactorBToPiDKMMO2008::MB_svz),

                make_observable("B->pi::f_B@DKMMO2008",
                        &AnalyticFormFactorBToPiDKMMO2008::decay_constant),
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
                        &AnalyticFormFactorBToPLCSR<lcsr::BToK>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B->K::f_0[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BToK>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B->K::f_T[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BToK>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),
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
            R"(Form factors for $B\to \bar{D}$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to\bar{D}$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
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
                        &AnalyticFormFactorBToPLCSR<lcsr::BToD>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B->D::f_0[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BToD>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B->D::f_T[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BToD>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B_s -> P(seudoscalar)
    // {{{

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

                make_observable("B_s->K::f_+[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToK>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B_s->K::f_0[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToK>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B_s->K::f_T[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToK>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),
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
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToDs>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s::f_0[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToDs>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s::f_T[s^1/s^0](q2)",
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToDs>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B -> V(ector)
    // {{{

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

                make_form_factor_adapter("B->rho::V(q2)/A_1(q2)",
                        &FormFactors<PToV>::v, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::A_2(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_2, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::A_12(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_12, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::T_23(q2)/T_2(q2)",
                        &FormFactors<PToV>::t_23, &FormFactors<PToV>::t_2,
                        std::make_tuple("q2")),

                make_observable("B->rho::A_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B->rho::A_2[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B->rho::A_30[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B->rho::V[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B->rho::T_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B->rho::T_23A[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B->rho::T_23B[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_t_23B,
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
                make_form_factor_adapter("B->K^*::V(q2)",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_0(q2)",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_1(q2)",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_2(q2)",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_12(q2)",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_1(q2)",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_2(q2)",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_3(q2)",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_23(q2)",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::V(q2)/A_1(q2)",
                        &FormFactors<PToV>::v, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_2(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_2, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_12(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_12, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_23(q2)/T_2(q2)",
                        &FormFactors<PToV>::t_23, &FormFactors<PToV>::t_2,
                        std::make_tuple("q2")),

                make_observable("B->K^*::A_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B->K^*::A_2[s^1/s^0](q2)",
                            &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_a_2,
                                std::make_tuple("q2")),

                make_observable("B->K^*::A_30[s^1/s^0](q2)",
                         &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_a_30,
                                std::make_tuple("q2")),

                make_observable("B->K^*::V[s^1/s^0](q2)",
                                &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_v,
                                std::make_tuple("q2")),

                make_observable("B->K^*::T_1[s^1/s^0](q2)",
                                &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_t_1,
                                std::make_tuple("q2")),

                make_observable("B->K^*::T_23A[s^1/s^0](q2)",
                                &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_t_23A,
                                std::make_tuple("q2")),

                make_observable("B->K^*::T_23B[s^1/s^0](q2)",
                                &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_t_23B,
                                std::make_tuple("q2")),
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
            R"(Form factors for $B\to \bar{D}^*$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to \bar{D}^*$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->D^*::V(q2)",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_0(q2)",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_1(q2)",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_2(q2)",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_12(q2)",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_1(q2)",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_2(q2)",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_3(q2)",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_23(q2)",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::V(q2)/A_1(q2)",
                        &FormFactors<PToV>::v, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_2(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_2, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_12(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_12, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_23(q2)/T_2(q2)",
                        &FormFactors<PToV>::t_23, &FormFactors<PToV>::t_2,
                        std::make_tuple("q2")),

                make_observable("B->D^*::A_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B->D^*::A_2[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B->D^*::A_30[s^1/s^0](q2)",
                         &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B->D^*::V[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B->D^*::T_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B->D^*::T_23A[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B->D^*::T_23B[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_t_23B,
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

                make_form_factor_adapter("B_s->K^*::V(q2)/A_1(q2)",
                        &FormFactors<PToV>::v, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::A_2(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_2, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::A_12(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_12, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::T_23(q2)/T_2(q2)",
                        &FormFactors<PToV>::t_23, &FormFactors<PToV>::t_2,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::A_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::A_2[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::A_30[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::V[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::T_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::T_23A[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::T_23B[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_t_23B,
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

                make_form_factor_adapter("B_s->phi::V(q2)/A_1(q2)",
                        &FormFactors<PToV>::v, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::A_2(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_2, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::A_12(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_12, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::T_23(q2)/T_2(q2)",
                        &FormFactors<PToV>::t_23, &FormFactors<PToV>::t_2,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::A_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::A_2[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::A_30[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::V[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::T_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::T_23A[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::T_23B[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),
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
                make_form_factor_adapter("B_s->D_s^*::V(q2)",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_0(q2)",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_1(q2)",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_2(q2)",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_12(q2)",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_1(q2)",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_2(q2)",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_3(q2)",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_23(q2)",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::V(q2)/A_1(q2)",
                        &FormFactors<PToV>::v, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_2(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_2, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_12(q2)/A_1(q2)",
                        &FormFactors<PToV>::a_12, &FormFactors<PToV>::a_1,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_23(q2)/T_2(q2)",
                        &FormFactors<PToV>::t_23, &FormFactors<PToV>::t_2,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::A_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::A_2[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::A_30[s^1/s^0](q2)",
                         &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::V[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::T_1[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::T_23A[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::T_23B[s^1/s^0](q2)",
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),
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
            R"(Form factors for $B_s\to \bar{K}^*$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to \bar{K}^*$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->pipi::Im{F_perp}(q2,k2,z)", R"(\text{Im}\,F_\perp^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_perp, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{F_para}(q2,k2,z)", R"(\text{Im}\,F_\parallel^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_para, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{F_long}(q2,k2,z)", R"(\text{Im}\,F_0^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_long, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{F_time}(q2,k2,z)", R"(\text{Im}\,F_t^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_time, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{Res{F_perp}}(q2,k2)", R"(\text{Res}\,\text{Im}\,F_\perp^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::f_perp_im_res_qhat2, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->pipi::Im{Res{F_para}}(q2,k2)", R"(\text{Res}\,\text{Im}\,F_\parallel^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::f_para_im_res_qhat2, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->pipi::Im{Res{F_long}}(q2,k2)", R"(\text{Res}\,\text{Im}\,F_0^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::f_long_im_res_qhat2, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->pipi::Im{Res{F_time}}(q2,k2)", R"(\text{Res}\,\text{Im}\,F_t^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::f_time_im_res_qhat2, std::make_tuple("q2", "k2")),
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
                        &ZeroRecoilSumRule<LambdaBToC>::vector_current),

                make_observable("Lambda_b->Lambda_c::G(1)",
                        &ZeroRecoilSumRule<LambdaBToC>::axialvector_current),

                make_observable("Lambda_b->Lambda_c::F_inel(1)",
                        &ZeroRecoilSumRule<LambdaBToC>::vector_current_inel),

                make_observable("Lambda_b->Lambda_c::G_inel(1)",
                        &ZeroRecoilSumRule<LambdaBToC>::axialvector_current_inel),
            }
        );

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
                make_b_to_k_form_factors_group(),
                make_b_to_d_form_factors_group(),

                // B -> V
                make_b_to_rho_form_factors_group(),
                make_b_to_kstar_form_factors_group(),
                make_b_to_dstar_form_factors_group(),

                // B_s -> V
                make_bs_to_kstar_form_factors_group(),

                // B -> P P
                make_b_to_pi_pi_form_factors_group(),

                // Lb -> 1/2^+
                make_lambdab_to_lambda_form_factors_group(),
                make_lambdab_to_lambdac_form_factors_group(),
            }
        );

        return ObservableSection(imp);
    }
}
