/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022-2024 Danny van Dyk
 * Copyright (c) 2022-2023 Philip Lüghausen
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

#include <eos/form-factors/heavy-meson-lcdas-flvd2022.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/matrix.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <limits>
#include <map>
#include <numeric>

#include <gsl/gsl_sf_gamma.h>

namespace eos
{
    namespace heavy_meson_lcdas
    {
        FLvD2022::FLvD2022(const Parameters & p, const Options & o) :
            model(Model::make("SM", p, o)),
            opt_Q(o, options, "Q"),
            opt_q(o, options, "q"),
            opt_gminus(o, options, "gminus"),
            switch_gminus(1.0),
            opt_alpha_s(o, options, "alpha_s"),
            mu_0(p[parameter("mu_0")], *this),
            omega_0(p[parameter("omega_0")], *this),
            a({
                UsedParameter(p[parameter("a^phi+_0")], *this),
                UsedParameter(p[parameter("a^phi+_1")], *this),
                UsedParameter(p[parameter("a^phi+_2")], *this),
                UsedParameter(p[parameter("a^phi+_3")], *this),
                UsedParameter(p[parameter("a^phi+_4")], *this),
                UsedParameter(p[parameter("a^phi+_5")], *this),
                UsedParameter(p[parameter("a^phi+_6")], *this),
                UsedParameter(p[parameter("a^phi+_7")], *this),
                UsedParameter(p[parameter("a^phi+_8")], *this)
            })
        {
            // Verify the size of Weights used internally
            Weights weights;
            if (weights.size() < a.size())
            {
                throw InternalError("The number of weights implemented is smaller than the number of coefficients of phi_+");
            }

            if (opt_gminus.value() == "zero")
            {
                switch_gminus = 0.0;
            }

            if (opt_alpha_s.value() == "full")
            {
                alpha_s = [this](const double & mu) -> double
                {
                    return model->alpha_s(mu);
                };
            }
            else
            {
                alpha_s = [](const double & mu)
                {
                    // This hardcodes the RGE for single-heavy hadrons containing a b quark
                    const double C_A = 3.0;
                    const double T_F = 1.0 / 2.0;
                    const double n_f = 5.0;
                    const double LambdaQCD = 0.213;

                    const double L = 2.0 * log(mu / LambdaQCD);
                    const double beta_0 = (11.0 / 3.0 * C_A - 4.0 / 3.0 * T_F * n_f);
                    return 4.0 * M_PI / beta_0 * 1.0 / L;
                };
            }
        }

        std::string
        FLvD2022::parameter(const char * _name) const
        {
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, qnp::Prefix> prefixes
            {
                { { QuarkFlavor::bottom, QuarkFlavor::up },      qnp::Prefix("B_u") },
                { { QuarkFlavor::bottom, QuarkFlavor::strange }, qnp::Prefix("B_s") }
            };

            auto it = prefixes.find(std::make_tuple(opt_Q.value(), opt_q.value()));
            if (it == prefixes.end())
                throw InternalError("Combination of options Q=" + opt_Q.str() + ", q=" + opt_q.str() + " is not supported");

            return QualifiedName(it->second, qnp::Name(_name), qnp::Suffix("FLvD2022")).str();
        }

        HeavyMesonLCDAs *
        FLvD2022::make(const Parameters & parameters, const Options & options)
        {
            return new FLvD2022(parameters, options);
        }

        std::tuple<HeavyMesonLCDAs::CoefficientIterator, HeavyMesonLCDAs::CoefficientIterator>
        FLvD2022::coefficient_range(const double & mu) const
        {
            // copy values to array of doubles
            static thread_local std::array<double, number_of_parameters> values;
            for (size_t i = 0; i < values.size(); i++)
            {
                values[i] = a[i]; // evaluates UsedParameter
            }

            // Perform RG evolution if mu != mu_0
            // Relation between a_k(mu) and a_k(mu_0) given in Eq. (56), [FLvD:2022A]
            if (std::abs(mu_0 - mu) > std::numeric_limits<double>::epsilon())
            {
                // This hardcodes the RGE for single-heavy hadrons containing a b quark
                const double C_A = 3.0;
                const double C_F = 4.0 / 3.0;
                const double T_F = 1.0 / 2.0;
                const double n_f = 5.0;

                auto gamma_cusp = [&](const double & alpha_s)
                {
                    const double a = alpha_s / (4.0 * M_PI);
                    return (
                            a                * 4.0 * C_F
                            + power_of<2>(a) * 4.0 * C_F
                                * ( (67.0 / 9.0 - power_of<2>(M_PI) / 3.0) * C_A - 20.0 / 9.0 * T_F * n_f )
                    );
                };

                auto beta = [&](const double & alpha_s)
                {
                    const double a = alpha_s / (4.0 * M_PI);
                    return -2.0 * alpha_s * (
                            a                * (11.0 / 3.0 * C_A - 4.0 / 3.0 * T_F * n_f)
                            + power_of<2>(a) * (34.0 / 3.0 * power_of<2>(C_A) - 20.0 / 3.0 * C_A * T_F * n_f - 4.0 * C_F * T_F * n_f)
                    );
                };

                auto gamma_plus = [&](const double alpha_s)
                {
                    return -2.0 * alpha_s * C_F / (4.0 * M_PI);
                };

                const double alpha_s_0 = this->alpha_s(mu_0);
                const double alpha_s_h = this->alpha_s(mu);

                auto g_integrand = [&] (const double alpha_s)
                {
                    return gamma_cusp(alpha_s) / beta(alpha_s);
                };
                const double g = integrate<GSL::QAGS>(g_integrand, alpha_s_0, alpha_s_h);

                auto V_integrand = [&] (const double alpha_s)
                {
                    auto inner_integrand = [&] (const double alpha_s)
                    {
                        return 1.0 / beta(alpha_s);
                    };
                    const double inner = integrate<GSL::QAGS>(inner_integrand, alpha_s_0, alpha_s);

                    return -1.0 / beta(alpha_s) * (gamma_cusp(alpha_s) * inner + gamma_plus(alpha_s));
                };
                const double V = integrate<GSL::QAGS>(V_integrand, alpha_s_0, alpha_s_h);

                const double g2  = g  * g;
                const double g3  = g2  * g;

                const std::array<std::array<double, 9>, 9> rge_matrix =
                {{
                     { 1., -0.5 * g, 0.16666666666666666 * g * (1. + g), -0.041666666666666664 * g * (1. + g) * (2. + g), 0.008333333333333333 * g * (1. + g) * (2. + g) * (3. + g), -0.001388888888888889 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g), 0.0001984126984126984 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g), -0.0000248015873015873 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g) * (6. + g), 2.7557319223985893e-6 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g) * (6. + g) * (7. + g) },
                     { -1. * g, 1. + 0.5 * (-1. + g) * g, -0.16666666666666666 * g * (4. + (-1. + g) * g), 0.041666666666666664 * g * (6. + 5. * g + g3), -0.008333333333333333 * g * (1. + g) * (2. + g) * (8. + (-1. + g) * g), 0.001388888888888889 * g * (1. + g) * (2. + g) * (3. + g) * (10. + (-1. + g) * g), -0.0001984126984126984 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (12. + (-1. + g) * g), 0.0000248015873015873 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g) * (14. + (-1. + g) * g), -2.7557319223985893e-6 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g) * (6. + g) * (16. + (-1. + g) * g) },
                     { 0.5 * g * (1. + g), -0.25 * g * (4. + (-1. + g) * g), 1. + 0.08333333333333333 * (-1. + g) * g * (10. + (-1. + g) * g), -0.020833333333333332 * g * (36. + (-1. + g) * g * (16. + (-1. + g) * g)), 0.004166666666666667 * g * (1. + g) * (4. + (-1. + g) * g) * (18. + (-1. + g) * g), -0.0006944444444444445 * g * (1. + g) * (2. + g) * (120. + (-1. + g) * g * (28. + (-1. + g) * g)), 0.0000992063492063492 * g * (1. + g) * (2. + g) * (3. + g) * (180. + (-1. + g) * g * (34. + (-1. + g) * g)), -0.00001240079365079365 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (252. + (-1. + g) * g * (40. + (-1. + g) * g)), 1.3778659611992946e-6 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g) * (336. + (-1. + g) * g * (46. + (-1. + g) * g)) },
                     { -0.16666666666666666 * g * (1. + g) * (2. + g), 0.08333333333333333 * g * (6. + 5. * g + g3), -0.027777777777777776 * g * (36. + (-1. + g) * g * (16. + (-1. + g) * g)), 0.006944444444444444 * (6. + (-1. + g) * g) * (24. + (-1. + g) * g * (22. + (-1. + g) * g)), -0.001388888888888889 * g * (576. + (-1. + g) * g * (348. + (-1. + g) * g * (40. + (-1. + g) * g))), 0.0002314814814814815 * g * (1. + g) * (1440. + (-1. + g) * g * (18. + (-1. + g) * g) * (34. + (-1. + g) * g)), -0.00003306878306878307 * g * (1. + g) * (2. + g) * (16. + (-1. + g) * g) * (180. + (-1. + g) * g * (48. + (-1. + g) * g)), 4.133597883597884e-6 * g * (1. + g) * (2. + g) * (3. + g) * (5040. + (-1. + g) * g * (1356. + (-1. + g) * g * (76. + (-1. + g) * g))), -4.592886537330982e-7 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (6. + (-1. + g) * g) * (1344. + (-1. + g) * g * (82. + (-1. + g) * g)) },
                     { 0.041666666666666664 * g * (1. + g) * (2. + g) * (3. + g), -0.020833333333333332 * g * (1. + g) * (2. + g) * (8. + (-1. + g) * g), 0.006944444444444444 * g * (1. + g) * (4. + (-1. + g) * g) * (18. + (-1. + g) * g), -0.001736111111111111 * g * (576. + (-1. + g) * g * (348. + (-1. + g) * g * (40. + (-1. + g) * g))), 1. + 0.00034722222222222224 * (-1. + g) * g * (14. + (-1. + g) * g) * (264. + (-1. + g) * g * (46. + (-1. + g) * g)), -0.00005787037037037037 * g * (14400. + (-1. + g) * g * (18. + (-1. + g) * g) * (592. + (-1. + g) * g * (62. + (-1. + g) * g))), 8.267195767195768e-6 * g * (1. + g) * (43200. + (-1. + g) * g * (18. + (-1. + g) * g) * (1272. + (-1. + g) * g * (82. + (-1. + g) * g))), -1.033399470899471e-6 * g * (1. + g) * (2. + g) * (100800. + (-1. + g) * g * (41856. + (-1. + g) * g * (4028. + (-1. + g) * g * (120. + (-1. + g) * g)))), 1.1482216343327455e-7 * g * (1. + g) * (2. + g) * (3. + g) * (201600. + (-1. + g) * g * (68976. + (-1. + g) * g * (76. + (-4. + g) * g) * (73. + g * (2. + g)))) },
                     { -0.008333333333333333 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g), 0.004166666666666667 * g * (1. + g) * (2. + g) * (3. + g) * (10. + (-1. + g) * g), -0.001388888888888889 * g * (1. + g) * (2. + g) * (120. + (-1. + g) * g * (28. + (-1. + g) * g)), 0.00034722222222222224 * g * (1. + g) * (1440. + (-1. + g) * g * (18. + (-1. + g) * g) * (34. + (-1. + g) * g)), -1. * g - 0.00006944444444444444 * (-1. + g) * (18. + (-1. + g) * g) * (592. + (-1. + g) * g * (62. + (-1. + g) * g)) * g2, 1. + 0.000011574074074074073 * (-1. + g) * g * (125280. + (-1. + g) * g * (37896. + (-1. + g) * g * (3508. + (-1. + g) * g * (110. + (-1. + g) * g)))), -1.6534391534391535e-6 * g * (518400. + (-1. + g) * g * (444960. + (-1. + g) * g * (89136. + (-1. + g) * g * (5908. + (-1. + g) * g * (140. + (-1. + g) * g))))), 2.066798941798942e-7 * g * (1. + g) * (1.8144e6 + (-1. + g) * g * (1.13184e6 + (-1. + g) * g * (171576. + (-1. + g) * g * (8908. + (-1. + g) * g * (170. + (-1. + g) * g))))), -2.296443268665491e-8 * g * (1. + g) * (2. + g) * (4.8384e6 + (-1. + g) * g * (2.38752e6 + (-1. + g) * g * (292416. + (-1. + g) * g * (12508. + (-1. + g) * g * (200. + (-1. + g) * g))))) },
                     { 0.001388888888888889 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g), -0.0006944444444444445 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (12. + (-1. + g) * g), 0.0002314814814814815 * g * (1. + g) * (2. + g) * (3. + g) * (180. + (-1. + g) * g * (34. + (-1. + g) * g)), -0.00005787037037037037 * g * (1. + g) * (2. + g) * (16. + (-1. + g) * g) * (180. + (-1. + g) * g * (48. + (-1. + g) * g)), 0.000011574074074074073 * g * (1. + g) * (43200. + (-1. + g) * g * (18. + (-1. + g) * g) * (1272. + (-1. + g) * g * (82. + (-1. + g) * g))), -1. * g - 1.9290123456790124e-6 * (-1. + g) * (444960. + (-1. + g) * g * (89136. + (-1. + g) * g * (5908. + (-1. + g) * g * (140. + (-1. + g) * g)))) * g2, 1. + 2.755731922398589e-7 * (-1. + g) * g * (5.78016e6 + (-1. + g) * g * (2.036592e6 + (-1. + g) * g * (236472. + (-1. + g) * g * (10528. + (-1. + g) * g * (182. + (-1. + g) * g))))), -0.875 * g - 3.444664902998236e-8 * (-1. + g) * (2.446848e7 + (-1. + g) * g * (5.780304e6 + (-1. + g) * g * (484608. + (-1. + g) * g * (16408. + (-1. + g) * g * (224. + (-1. + g) * g))))) * g2, 3.827405447775818e-9 * g * (1. + g) * (1.016064e8 + (-1. + g) * g * (7.200576e7 + (-1. + g) * g * (1.2986496e7 + (-1. + g) * g * (858744. + (-1. + g) * g * (23548. + (-1. + g) * g * (266. + (-1. + g) * g)))))) },
                     { -0.0001984126984126984 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g) * (6. + g), 0.0000992063492063492 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g) * (14. + (-1. + g) * g), -0.00003306878306878307 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (252. + (-1. + g) * g * (40. + (-1. + g) * g)), 8.267195767195768e-6 * g * (1. + g) * (2. + g) * (3. + g) * (5040. + (-1. + g) * g * (1356. + (-1. + g) * g * (76. + (-1. + g) * g))), -1.6534391534391535e-6 * g * (1. + g) * (2. + g) * (100800. + (-1. + g) * g * (41856. + (-1. + g) * g * (4028. + (-1. + g) * g * (120. + (-1. + g) * g)))), 2.755731922398589e-7 * g * (1. + g) * (1.8144e6 + (-1. + g) * g * (1.13184e6 + (-1. + g) * g * (171576. + (-1. + g) * g * (8908. + (-1. + g) * g * (170. + (-1. + g) * g))))), -1. * g - 3.936759889140842e-8 * (-1. + g) * (2.446848e7 + (-1. + g) * g * (5.780304e6 + (-1. + g) * g * (484608. + (-1. + g) * g * (16408. + (-1. + g) * g * (224. + (-1. + g) * g))))) * g2, 1. + 4.920949861426052e-9 * (-1. + g) * g * (3.4909056e8 + (-1. + g) * g * (1.38517632e8 + (-1. + g) * g * (1.9022736e7 + (-1. + g) * g * (1.074176e6 + (-1. + g) * g * (26600. + (-1. + g) * g * (280. + (-1. + g) * g)))))), -0.8888888888888888 * g - 5.467722068251169e-10 * (-1. + g) * (1.71932544e9 + (-1. + g) * g * (4.62214656e8 + (-1. + g) * g * (4.6160784e7 + (-1. + g) * g * (1.993024e6 + (-1. + g) * g * (39144. + (-1. + g) * g * (336. + (-1. + g) * g)))))) * g2 },
                     { 0.0000248015873015873 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g) * (6. + g) * (7. + g), -0.00001240079365079365 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g) * (6. + g) * (16. + (-1. + g) * g), 4.133597883597884e-6 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (5. + g) * (336. + (-1. + g) * g * (46. + (-1. + g) * g)), -1.033399470899471e-6 * g * (1. + g) * (2. + g) * (3. + g) * (4. + g) * (6. + (-1. + g) * g) * (1344. + (-1. + g) * g * (82. + (-1. + g) * g)), 2.066798941798942e-7 * g * (1. + g) * (2. + g) * (3. + g) * (201600. + (-1. + g) * g * (68976. + (-1. + g) * g * (76. + (-4. + g) * g) * (73. + g * (2. + g)))), -3.444664902998236e-8 * g * (1. + g) * (2. + g) * (4.8384e6 + (-1. + g) * g * (2.38752e6 + (-1. + g) * g * (292416. + (-1. + g) * g * (12508. + (-1. + g) * g * (200. + (-1. + g) * g))))), 4.920949861426052e-9 * g * (1. + g) * (1.016064e8 + (-1. + g) * g * (7.200576e7 + (-1. + g) * g * (1.2986496e7 + (-1. + g) * g * (858744. + (-1. + g) * g * (23548. + (-1. + g) * g * (266. + (-1. + g) * g)))))), -1. * g - 6.151187326782565e-10 * (-1. + g) * (1.71932544e9 + (-1. + g) * g * (4.62214656e8 + (-1. + g) * g * (4.6160784e7 + (-1. + g) * g * (1.993024e6 + (-1. + g) * g * (39144. + (-1. + g) * g * (336. + (-1. + g) * g)))))) * g2, 1. + 6.834652585313961e-11 * (-1. + g) * g * (2.676022272e10 + (-1. + g) * g * (1.1692594944e10 + (-1. + g) * g * (1.831851648e9 + (-1. + g) * g * (1.23501456e8 + (-1. + g) * g * (3.908224e6 + (-1. + g) * g * (59304. + (-1. + g) * g * (408. + (-1. + g) * g))))))) },
                }};

                const double exp_2eulergamma = 3.1722189581254505;
                values = gsl_sf_gamma(2.0 - g) * std::exp(V) * std::pow(mu_0 * exp_2eulergamma / (2.0 * omega_0), -g) * rge_matrix * values;
            }

            return {values.begin(), values.end()};
        }


        /* Leading twist two-particle LCDAs */

        double
        FLvD2022::phi_plus(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_minus(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }


        double
        FLvD2022::phi_bar(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar_d1(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phitilde_plus(const double & tau, const double & mu) const
        {
            const double x = tau * omega_0;

            const double p1 = (x - 1.0) / (x + 1.0);
            const double p2 = p1 * p1;
            const double p3 = p2 * p1;
            const double p4 = p3 * p1;
            const double p5 = p4 * p1;
            const double p6 = p5 * p1;
            const double p7 = p6 * p1;
            const double p8 = p7 * p1;

            const Weights c = {
                1.0, p1, p2, p3, p4, p5, p6, p7, p8
            };

            auto [a_begin, a_end] = this->coefficient_range(mu);
            return 1.0 / power_of<2>(1.0 + x) * std::inner_product(a_begin, a_end, c.begin(), 0.0);
        }

        double
        FLvD2022::t_d_dt_phitilde_plus(const double & tau, const double & mu) const
        {
            const double x = tau * omega_0;

            const double p1 = (x - 1.0) / (x + 1.0);
            const double p2 = p1 * p1;
            const double p3 = p2 * p1;
            const double p4 = p3 * p1;
            const double p5 = p4 * p1;
            const double p6 = p5 * p1;
            const double p7 = p6 * p1;
            const double p8 = p7 * p1;

            const Weights c = {
                1.0 - x,
                (2.0 - x) * p1,
                (3.0 - x) * p2,
                (4.0 - x) * p3,
                (5.0 - x) * p4,
                (6.0 - x) * p5,
                (7.0 - x) * p6,
                (8.0 - x) * p7,
                (9.0 - x) * p8
            };

            auto [a_begin, a_end] = this->coefficient_range(mu);
            return 2.0 * x / power_of<3>(x + 1.0) / (x - 1.0) * std::inner_product(a_begin, a_end, c.begin(), 0.0);
        }

        double
        FLvD2022::t2_d2_d2t_phitilde_plus(const double & tau, const double & mu) const
        {
            const double x = tau * omega_0;

            const double p1 = (x - 1.0) / (x + 1.0);
            const double p2 = p1 * p1;
            const double p3 = p2 * p1;
            const double p4 = p3 * p1;
            const double p5 = p4 * p1;
            const double p6 = p5 * p1;
            const double p7 = p6 * p1;
            const double p8 = p7 * p1;

            const double xminus2 = (1.0 - x) * (1.0 - x);

            const Weights c = {
                3 * xminus2,
                p1 * (6 - 6 * x + 3 * xminus2),
                p2 * (8 + 2 * (4 - 6 * x) + 3 * xminus2),
                p3 * (18 + 3 * (4 - 6 * x) + 3 * xminus2),
                p4 * (32 + 4 * (4 - 6 * x) + 3 * xminus2),
                p5 * (50 + 5 * (4 - 6 * x) + 3 * xminus2),
                p6 * (72 + 6 * (4 - 6 * x) + 3 * xminus2),
                p7 * (98 + 7 * (4 - 6 * x) + 3 * xminus2),
                p8 * (128 + 8 * (4 - 6 * x) + 3 * xminus2)
            };

            auto [a_begin, a_end] = this->coefficient_range(mu);
            return 2.0 * power_of<2>(x) / power_of<2>(x - 1.0) / power_of<4>(x + 1.0) * std::inner_product(a_begin, a_end, c.begin(), 0.0);
        }

        /* Next-to-leading twist two-particle LCDAs */

        double
        FLvD2022::g_minusWW(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_minusWW_d1(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_minusWW_d2(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_plus(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_plus_d1(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_plus_d2(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_bar(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_bar_d1(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_bar_d2(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_bar_d3(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        /* Leading twist three-particle LCDAs */

        double
        FLvD2022::phi_3(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar_3(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar2_3(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar2_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar_bar_3(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::psi_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::psi_bar_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }


        double
        FLvD2022::chi_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::chi_bar_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::inverse_lambda_plus() const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::inverse_moment(const double & mu) const
        {
            // Cp. [FLvD:2022A], Eq. (43)
            const Weights c = {
                1.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 5.0, 0.0, 1.0 / 7.0, 0.0, 1.0 / 9.0
            };

            auto [a_begin, a_end] = this->coefficient_range(mu);
            return 1.0 / omega_0 * std::inner_product(a_begin, a_end, c.begin(), 0.0);
        }

        double
        FLvD2022::logarithmic_moment_1(const double & mu) const
        {
            // Cp. [FLvD:2022A], Eq. (44)
            const Weights c = {
                0, -1.0, 0, -2.0 / 3.0, 0.0, -23.0 / 45.0, 0.0, -44.0 / 105.0, 0.0
            };

            auto [a_begin, a_end] = this->coefficient_range(mu);
            return 1.0 / omega_0 * std::inner_product(a_begin, a_end, c.begin(), 0.0);
        }

        double
        FLvD2022::logarithmic_moment_2(const double & mu) const
        {
            // Cp. [FLvD:2022A], Eq. (45)
            const Weights c = {
                0.0, 0.0, 4.0 / 3.0, 0.0, 4.0 / 3.0, 0.0, 56.0 / 45.0, 0.0, 3272.0 / 2835.0
            };

            auto [a_begin, a_end] = this->coefficient_range(mu);
            return power_of<2>(M_PI) / 6.0 * inverse_moment(mu) + 1.0 / omega_0 * std::inner_product(a_begin, a_end, c.begin(), 0.0);
        }

        double
        FLvD2022::psi_A(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::psi_V(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::X_A(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::Y_A(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::Xbar_A(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::Ybar_A(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        Diagnostics
        FLvD2022::diagnostics() const
        {
            Diagnostics results;
            // add diagnostic results here
            return results;
        }

        const std::set<ReferenceName>
        FLvD2022::references
        {
            "FLvD:2022A"_rn
        };

        const std::vector<OptionSpecification>
        FLvD2022::options
        {
            { "Q",       { "b" },                "b"        },
            { "q",       { "u", "s" },           "u"        },
            { "gminus",  { "zero", "WW-limit" }, "WW-limit" },
            { "alpha_s", { "naive", "full"  },   "full"     },
        };

        std::vector<OptionSpecification>::const_iterator
        FLvD2022::begin_options()
        {
            return options.cbegin();
        }

        std::vector<OptionSpecification>::const_iterator
        FLvD2022::end_options()
        {
            return options.cend();
        }
    }

    template <>
    struct WrappedForwardIteratorTraits<HeavyMesonLCDAs::CoefficientIteratorTag>
    {
        using UnderlyingIterator = std::array<double, 9>::const_iterator;
    };
    template class WrappedForwardIterator<HeavyMesonLCDAs::CoefficientIteratorTag, const double &>;
}
