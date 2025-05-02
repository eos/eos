/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
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

#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/power-of.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda1520-ll-base.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda1520-ll-naive.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    /*!
     * Implementation for the decay @f$\bar{\Lambda_b} \to \bar{\Lambda}(1520) \ell^+ \ell^-@f$.
     */
    template <>
    struct Implementation<LambdaBToLambda1520Dilepton>
    {
        std::shared_ptr<LambdaBToLambda1520Dilepton::AmplitudeGenerator> amplitude_generator;

        std::shared_ptr<Model> model;

        LeptonFlavorOption opt_l;

        UsedParameter hbar;
        UsedParameter m_l;
        UsedParameter tau;
        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "WET"), p, o)),
            opt_l(o, options, "l"_ok),
            hbar(p["QM::hbar"], u),
            m_l(p["mass::" + opt_l.str()], u),
            tau(p["life_time::Lambda_b"], u),
            mu(p["sb" + opt_l.str() + opt_l.str() + "::mu"], u)
        {
            Context ctx("When constructing Lb->L(1520)ll observables");

            std::string tag = o.get("tag"_ok, "");

            if ("Naive" == tag)
            {
                amplitude_generator.reset(new LambdaBToLambda1520DileptonAmplitudes<tag::Naive>(p, o));
            }
            else
            {
                throw InternalError("LambdaBToLambda1520Dilepton: Unknown tag or no valid tag specified (tag = '" + tag + "')!");
            }

            u.uses(*amplitude_generator);
        }

        ~Implementation()
        {
        }

        inline std::array<double, 12> angular_coefficients_array(const LambdaBToLambda1520Dilepton::Amplitudes & A, const double & s) const
        {
            // cf. [DD:2020A], app. G, which agrees with [DN:2019A], eq. (4.2) for massless leptons
            std::array<double, 12> result;

            double z = 4.0 * power_of<2>(m_l()) / s;
            double y = m_l / std::sqrt(s);
            double beta2 = 1.0 - z;
            double beta = std::sqrt(beta2);

            // L1c
            result[0] = - 2.0 * beta * real(
                A.a_perp1_left  * conj(A.a_para1_left) - A.a_perp1_right * conj(A.a_para1_right)
                + y * (
                      A.a_paraS_left  * conj(A.a_para0_left)  + A.a_paraS_right * conj(A.a_para0_left)
                    + A.a_perpS_left  * conj(A.a_perp0_left)  + A.a_perpS_right * conj(A.a_perp0_left)
                    + A.a_paraS_right * conj(A.a_para0_right) + A.a_paraS_left  * conj(A.a_para0_right)
                    + A.a_perpS_right * conj(A.a_perp0_right) + A.a_perpS_left  * conj(A.a_perp0_right)
                )
            );

            // L1cc
            result[1] =
                  norm(A.a_para1_left) + norm(A.a_perp1_left) + norm(A.a_para1_right) + norm(A.a_perp1_right)
                + norm(A.a_paraS_left) + norm(A.a_perpS_left) + norm(A.a_paraS_right) + norm(A.a_perpS_right)
                + 2.0 * y * real(
                    - A.a_parat_right * conj(A.a_paraS_left)  + A.a_paraS_left  * conj(A.a_parat_left)
                    - A.a_perpt_right * conj(A.a_perpS_left)  + A.a_perpS_left  * conj(A.a_perpt_left)
                    - A.a_parat_left  * conj(A.a_paraS_right) + A.a_paraS_right * conj(A.a_parat_right)
                    - A.a_perpt_left  * conj(A.a_perpS_right) + A.a_perpS_right * conj(A.a_perpt_right)
                ) +  2.0 * y * y * (
                      norm(A.a_para0_left) - norm(A.a_para1_left) - norm(A.a_paraS_left) + norm(A.a_parat_left)
                    + norm(A.a_perp0_left) - norm(A.a_perp1_left) - norm(A.a_perpS_left) + norm(A.a_perpt_left)
                    + norm(A.a_para0_right) - norm(A.a_para1_right) - norm(A.a_paraS_right) + norm(A.a_parat_right)
                    + norm(A.a_perp0_right) - norm(A.a_perp1_right) - norm(A.a_perpS_right) + norm(A.a_perpt_right)
                ) + 2.0 * y * y * real(
                    A.a_para0_right * conj(A.a_para0_left)  + A.a_para1_right * conj(A.a_para1_left)
                    - A.a_paraS_right * conj(A.a_paraS_left)  - A.a_parat_right * conj(A.a_parat_left)
                    + A.a_perp0_right * conj(A.a_perp0_left)  + A.a_perp1_right * conj(A.a_perp1_left)
                    - A.a_perpS_right * conj(A.a_perpS_left)  - A.a_perpt_right * conj(A.a_perpt_left)
                    + A.a_para0_left  * conj(A.a_para0_right) + A.a_para1_left  * conj(A.a_para1_right)
                    - A.a_paraS_left  * conj(A.a_paraS_right) - A.a_parat_left  * conj(A.a_parat_right)
                    + A.a_perp0_left  * conj(A.a_perp0_right) + A.a_perp1_left  * conj(A.a_perp1_right)
                    - A.a_perpS_left  * conj(A.a_perpS_right) - A.a_perpt_left  * conj(A.a_perpt_right)
            );

            // L1ss
            result[2] = 0.5 * (
                  2.0 * norm(A.a_para0_left)  + 2.0 * norm(A.a_perp0_left)  + norm(A.a_para1_left)  + norm(A.a_perp1_left)
                + 2.0 * norm(A.a_para0_right) + 2.0 * norm(A.a_perp0_right) + norm(A.a_para1_right) + norm(A.a_perp1_right)
                + 2.0 * norm(A.a_paraS_left)  + 2.0 * norm(A.a_perpS_left)
                + 2.0 * norm(A.a_paraS_right) + 2.0 * norm(A.a_perpS_right)
                ) + 2.0 * y * real(
                    - A.a_parat_right * conj(A.a_paraS_left)  + A.a_paraS_left  * conj(A.a_parat_left)
                    - A.a_perpt_right * conj(A.a_perpS_left)  + A.a_perpS_left  * conj(A.a_perpt_left)
                    - A.a_parat_left  * conj(A.a_paraS_right) + A.a_paraS_right * conj(A.a_parat_right)
                    - A.a_perpt_left  * conj(A.a_perpS_right) + A.a_perpS_right * conj(A.a_perpt_right)
                ) + 2.0 * y * y * (
                    - norm(A.a_para0_left)  - norm(A.a_paraS_left)  + norm(A.a_parat_left)
                    - norm(A.a_perp0_left)  - norm(A.a_perpS_left)  + norm(A.a_perpt_left)
                    - norm(A.a_para0_right) - norm(A.a_paraS_right) + norm(A.a_parat_right)
                    - norm(A.a_perp0_right) - norm(A.a_perpS_right) + norm(A.a_perpt_right)
                ) + 2.0 * y * y * real(
                      A.a_para0_right * conj(A.a_para0_left)  + A.a_para1_right * conj(A.a_para1_left)
                    - A.a_paraS_right * conj(A.a_paraS_left)  - A.a_parat_right * conj(A.a_parat_left)
                    + A.a_perp0_right * conj(A.a_perp0_left)  + A.a_perp1_right * conj(A.a_perp1_left)
                    - A.a_perpS_right * conj(A.a_perpS_left)  - A.a_perpt_right * conj(A.a_perpt_left)
                    + A.a_para0_left  * conj(A.a_para0_right) + A.a_para1_left  * conj(A.a_para1_right)
                    - A.a_paraS_left  * conj(A.a_paraS_right) - A.a_parat_left  * conj(A.a_parat_right)
                    + A.a_perp0_left  * conj(A.a_perp0_right) + A.a_perp1_left  * conj(A.a_perp1_right)
                    - A.a_perpS_left  * conj(A.a_perpS_right) - A.a_perpt_left  * conj(A.a_perpt_right)
            );

            // L2c
            result[3] = - 0.5 * beta * real(
                  A.a_perp1_left  * conj(A.a_para1_left)  + 3.0 * A.b_perp1_left  * conj(A.b_para1_left)
                - A.a_perp1_right * conj(A.a_para1_right) - 3.0 * A.b_perp1_right * conj(A.b_para1_right)
                + y * (
                      A.a_paraS_left  * conj(A.a_para0_left)  + A.a_paraS_right * conj(A.a_para0_left)
                    + A.a_perpS_left  * conj(A.a_perp0_left)  + A.a_perpS_right * conj(A.a_perp0_left)
                    + A.a_paraS_right * conj(A.a_para0_right) + A.a_paraS_left  * conj(A.a_para0_right)
                    + A.a_perpS_right * conj(A.a_perp0_right) + A.a_perpS_left  * conj(A.a_perp0_right)
                )
            );


            // L2cc
            result[4] = 0.25 * (
                  norm(A.a_para1_left)  + norm(A.a_perp1_left)  + 3.0 * norm(A.b_para1_left)  + 3.0 * norm(A.b_perp1_left)
                + norm(A.a_para1_right) + norm(A.a_perp1_right) + 3.0 * norm(A.b_para1_right) + 3.0 * norm(A.b_perp1_right)
                + norm(A.a_paraS_left)  + norm(A.a_perpS_left)  + norm(A.a_paraS_right)  + norm(A.a_perpS_right)
                ) + 0.5 * y * real(
                    - A.a_parat_right * conj(A.a_paraS_left)  + A.a_paraS_left  * conj(A.a_parat_left)
                    - A.a_perpt_right * conj(A.a_perpS_left)  + A.a_perpS_left  * conj(A.a_perpt_left)
                    - A.a_parat_left  * conj(A.a_paraS_right) + A.a_paraS_right * conj(A.a_parat_right)
                    - A.a_perpt_left  * conj(A.a_perpS_right) + A.a_perpS_right * conj(A.a_perpt_right)
                ) + 0.5 * y * y * (
                      norm(A.a_para0_left)  - norm(A.a_para1_left)  - norm(A.a_paraS_left)  + norm(A.a_parat_left)
                    + norm(A.a_perp0_left)  - norm(A.a_perp1_left)  - norm(A.a_perpS_left)  + norm(A.a_perpt_left)
                    - 3.0 * norm(A.b_para1_left)  - 3.0 * norm(A.b_perp1_left)
                    + norm(A.a_para0_right) - norm(A.a_para1_right) - norm(A.a_paraS_right) + norm(A.a_parat_right)
                    + norm(A.a_perp0_right) - norm(A.a_perp1_right) - norm(A.a_perpS_right) + norm(A.a_perpt_right)
                    - 3.0 * norm(A.b_para1_right) - 3.0 * norm(A.b_perp1_right)
                ) + 0.5 * y * y * real(
                      A.a_para0_right * conj(A.a_para0_left)  + A.a_para1_right * conj(A.a_para1_left)
                    - A.a_paraS_right * conj(A.a_paraS_left)  - A.a_parat_right * conj(A.a_parat_left)
                    + A.a_perp0_right * conj(A.a_perp0_left)  + A.a_perp1_right * conj(A.a_perp1_left)
                    - A.a_perpS_right * conj(A.a_perpS_left)  - A.a_perpt_right * conj(A.a_perpt_left)
                    + 3.0 * A.b_para1_right * conj(A.b_para1_left) + 3.0 * A.b_perp1_right * conj(A.b_perp1_left)
                    + A.a_para0_left  * conj(A.a_para0_right) + A.a_para1_left  * conj(A.a_para1_right)
                    - A.a_paraS_left  * conj(A.a_paraS_right) - A.a_parat_left  * conj(A.a_parat_right)
                    + A.a_perp0_left  * conj(A.a_perp0_right) + A.a_perp1_left  * conj(A.a_perp1_right)
                    - A.a_perpS_left  * conj(A.a_perpS_right) - A.a_perpt_left  * conj(A.a_perpt_right)
                    + 3.0 * A.b_para1_left * conj(A.b_para1_right) + 3.0 * A.b_perp1_left * conj(A.b_perp1_right)
            );

            // L2ss
            result[5] = 0.125 * (
                  2.0 * norm(A.a_para0_left) + norm(A.a_para1_left) + 2.0 * norm(A.a_perp0_left) + norm(A.a_perp1_left)
                + 2.0 * norm(A.a_paraS_left) + 2.0 * norm(A.a_perpS_left) + 3.0 * norm(A.b_para1_left) + 3.0 * norm(A.b_perp1_left)
                + 2.0 * norm(A.a_para0_right) + norm(A.a_para1_right) + 2.0 * norm(A.a_perp0_right) + norm(A.a_perp1_right)
                + 2.0 * norm(A.a_paraS_right) + 2.0 * norm(A.a_perpS_right) + 3.0 * norm(A.b_para1_right) + 3.0 * norm(A.b_perp1_right)
                - 2.0 * sqrt(3.0) * real(
                    A.b_para1_left  * conj(A.a_para1_left)  - A.b_perp1_left *  conj(A.a_perp1_left)
                  + A.b_para1_right * conj(A.a_para1_right) - A.b_perp1_right * conj(A.a_perp1_right)
                )
            ) + 0.5 * y * real(
                    - A.a_parat_right * conj(A.a_paraS_left)  + A.a_paraS_left  * conj(A.a_parat_left)
                    - A.a_perpt_right * conj(A.a_perpS_left)  + A.a_perpS_left  * conj(A.a_perpt_left)
                    - A.a_parat_left  * conj(A.a_paraS_right) + A.a_paraS_right * conj(A.a_parat_right)
                    - A.a_perpt_left  * conj(A.a_perpS_right) + A.a_perpS_right * conj(A.a_perpt_right)
            ) + 0.5 * y * y * (
                    - norm(A.a_para0_left)  - norm(A.a_paraS_left)  + norm(A.a_parat_left)
                    - norm(A.a_perp0_left)  - norm(A.a_perpS_left)  + norm(A.a_perpt_left)
                    - norm(A.a_para0_right) - norm(A.a_paraS_right) + norm(A.a_parat_right)
                    - norm(A.a_perp0_right) - norm(A.a_perpS_right) + norm(A.a_perpt_right)
            ) + 0.5 * y * y * real(
                      A.a_para0_right * conj(A.a_para0_left)  + A.a_para1_right * conj(A.a_para1_left)
                    - A.a_paraS_right * conj(A.a_paraS_left)  - A.a_parat_right * conj(A.a_parat_left)
                    + A.a_perp0_right * conj(A.a_perp0_left)  + A.a_perp1_right * conj(A.a_perp1_left)
                    - A.a_perpS_right * conj(A.a_perpS_left)  - A.a_perpt_right * conj(A.a_perpt_left)
                    + 2.0 * sqrt(3.0) * (A.b_para1_left * conj(A.a_para1_left) - A.b_perp1_left * conj(A.a_perp1_left))
                    + 3.0 * A.b_para1_right * conj(A.a_para1_left) + 3.0 * A.b_perp1_right * conj(A.a_perp1_left)
                    + A.a_para0_left  * conj(A.a_para0_right) + A.a_para1_left  * conj(A.a_para1_right)
                    - A.a_paraS_left  * conj(A.a_paraS_right) - A.a_parat_left  * conj(A.a_parat_right)
                    + A.a_perp0_left  * conj(A.a_perp0_right) + A.a_perp1_left  * conj(A.a_perp1_right)
                    - A.a_perpS_left  * conj(A.a_perpS_right) - A.a_perpt_left  * conj(A.a_perpt_right)
                    + 2.0 * sqrt(3.0) * (A.b_para1_right * conj(A.a_para1_right) - A.b_perp1_right * conj(A.a_perp1_right))
                    + 3.0 * A.b_para1_left * conj(A.a_para1_right) + 3.0 * A.b_perp1_left * conj(A.a_perp1_right)
            );

            // L3ss
            result[6] = sqrt(3.0) / 2.0 * beta2 * real(
                  A.b_para1_left  * conj(A.a_para1_left)  - A.b_perp1_left  * conj(A.a_perp1_left)
                + A.b_para1_right * conj(A.a_para1_right) - A.b_perp1_right * conj(A.a_perp1_right)
            );

            // L4ss
            result[7] = sqrt(3.0) / 2.0 * beta2 * imag(
                  A.b_perp1_left  * conj(A.a_para1_left)  - A.b_para1_left  * conj(A.a_perp1_left)
                + A.b_perp1_right * conj(A.a_para1_right) - A.b_para1_right * conj(A.a_perp1_right)
            );

            // L5s
            result[8] = sqrt(3.0 / 2.0) * beta * real(
                  A.b_perp1_left  * conj(A.a_para0_left)  - A.b_para1_left  * conj(A.a_perp0_left)
                - A.b_perp1_right * conj(A.a_para0_right) - A.b_para1_right * conj(A.a_perp0_right)
                - y * (
                      A.b_para1_right * conj(A.a_paraS_left)  - A.b_perp1_right * conj(A.a_perpS_left)
                    + A.a_paraS_left  * conj(A.b_para1_left)  - A.a_perpS_left  * conj(A.b_perp1_left)
                    + A.b_para1_left  * conj(A.a_paraS_right) - A.b_perp1_left  * conj(A.a_perpS_right)
                    + A.a_paraS_right * conj(A.b_para1_right) - A.a_perpS_right * conj(A.b_perp1_right)
                )
            );

            // L5sc
            result[9] = - sqrt(3.0 / 2.0) * beta2 * real(
                  A.b_para1_left  * conj(A.a_para0_left)  - A.b_perp1_left  * conj(A.a_perp0_left)
                + A.b_para1_right * conj(A.a_para0_right) - A.b_perp1_right * conj(A.a_perp0_right)
            );

            // L6s
            result[10] = sqrt(3.0 / 2.0) * beta * imag(
                  A.b_para1_left  * conj(A.a_para0_left)  - A.b_perp1_left  * conj(A.a_perp0_left)
                - A.b_perp1_right * conj(A.a_para0_right) + A.b_para1_right * conj(A.a_perp0_right)
                - y * (
                      A.b_perp1_right * conj(A.a_paraS_left)  - A.b_para1_right * conj(A.a_perpS_left)
                    + A.a_perpS_left  * conj(A.b_para1_left)  - A.a_paraS_left  * conj(A.b_perp1_left)
                    + A.b_perp1_left  * conj(A.a_paraS_right) - A.b_para1_left  * conj(A.a_perpS_right)
                    + A.a_perpS_right * conj(A.b_para1_right) - A.a_paraS_right * conj(A.b_perp1_right)
                )
            );

            // L6sc
            result[11] = - sqrt(3.0 / 2.0) * beta2 * imag(
                  A.b_perp1_left  * conj(A.a_para0_left)  - A.b_para1_left * conj(A.a_perp0_left)
                + A.b_perp1_right * conj(A.a_para0_right) - A.b_para1_right * conj(A.a_perp0_right)
            );

            return result;
        }

        inline std::array<double, 12> differential_angular_coefficients_array(const double & s) const
        {
            return angular_coefficients_array(amplitude_generator->amplitudes(s), s);
        }

        inline LambdaBToLambda1520Dilepton::AngularCoefficients differential_angular_coefficients(const double & s) const
        {
            return LambdaBToLambda1520Dilepton::AngularCoefficients(differential_angular_coefficients_array(s));
        }

        LambdaBToLambda1520Dilepton::AngularCoefficients integrated_angular_coefficients(const double & s_min, const double & s_max) const
        {
            std::function<std::array<double, 12> (const double &)> integrand =
                    std::bind(&Implementation<LambdaBToLambda1520Dilepton>::differential_angular_coefficients_array, this, std::placeholders::_1);
            std::array<double, 12> integrated_angular_coefficients_array = integrate1D(integrand, 64, s_min, s_max);

            return LambdaBToLambda1520Dilepton::AngularCoefficients(integrated_angular_coefficients_array);
        }

        inline double decay_width(const LambdaBToLambda1520Dilepton::AngularCoefficients & a_c)
        {
            // cf. [DN:2019A], eq. (4.4)
            return 1.0 / 3.0 * (a_c.L1cc + 2.0 * a_c.L1ss + 2.0 * a_c.L2cc + 4.0 * a_c.L2ss + 2.0 * a_c.L3ss);
        }
    };

    LambdaBToLambda1520Dilepton::LambdaBToLambda1520Dilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<LambdaBToLambda1520Dilepton>(new Implementation<LambdaBToLambda1520Dilepton>(parameters, options, *this))
    {
    }

    LambdaBToLambda1520Dilepton::~LambdaBToLambda1520Dilepton()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<LambdaBToLambda1520Dilepton>::options
    {
        Model::option_specification(),
        {"l"_ok, { "e", "mu", "tau" }, "mu"}
    };

    double
    LambdaBToLambda1520Dilepton::decay_width(const double & s, const double & c_theta_l, const double & c_theta_Lstar, const double & phi) const
    {
        // compute d^4 Gamma, cf. [DN:2019A], eq. (4.1)
        // Cosine squared of the angles
        double c_theta_Lstar_2 = c_theta_Lstar * c_theta_Lstar;
        double c_theta_l_2 = c_theta_l * c_theta_l;
        double c_phi = cos(phi);
        double c_phi_2 = c_phi * c_phi;
        // Sine squared of the angles
        double s_theta_Lstar_2 = 1.0 - c_theta_Lstar_2;
        double s_theta_l_2 = 1.0 - c_theta_l_2;
        // Sine of the angles
        double s_theta_Lstar = sqrt(s_theta_Lstar_2);
        double s_theta_l = sqrt(s_theta_l_2);
        double s_phi = sin(phi);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);

        double result = 3.0 / 8.0 / M_PI * (
                c_theta_Lstar_2 * (
                    a_c.L1c * c_theta_l + a_c.L1cc * c_theta_l_2 + a_c.L1ss * s_theta_l_2
                  )
                  + s_theta_Lstar_2 * (
                      a_c.L2c * c_theta_l + a_c.L2cc * c_theta_l_2 + a_c.L2ss * s_theta_l_2
                    + a_c.L3ss * s_theta_l_2 * c_phi_2 + a_c.L4ss * s_theta_l_2 * s_phi * c_phi
                  )
                  + s_theta_Lstar * c_theta_Lstar * c_phi * (
                      a_c.L5s * s_theta_l + a_c.L5sc * s_theta_l * c_theta_l
                  )
                  + s_theta_Lstar * c_theta_Lstar * s_phi * (
                      a_c.L6s * s_theta_l + a_c.L6sc * s_theta_l * c_theta_l
                  )
                );

        return result;
    }

    double
    LambdaBToLambda1520Dilepton::differential_decay_width(const double & s) const
    {
        return _imp->decay_width(_imp->differential_angular_coefficients(s));
    }

    double
    LambdaBToLambda1520Dilepton::differential_branching_ratio(const double & s) const
    {
        return differential_decay_width(s) * _imp->tau() / _imp->hbar();
    }

    double
    LambdaBToLambda1520Dilepton::differential_forward_backward_asymmetry(const double & s) const
    {
        // cf. [DN:2019A], eq. (4.7)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return (a_c.L1c + 2.0 * a_c.L2c) / 2.0 / _imp->decay_width(a_c);
    }

    double
    LambdaBToLambda1520Dilepton::differential_longitudinal_polarisation(const double & s) const
    {
        // cf. [DN:2019A], eq. (4.6)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 1 - 2.0 * (a_c.L1cc + 2.0 * a_c.L2cc) / 3.0 / _imp->decay_width(a_c);
    }

    double
    LambdaBToLambda1520Dilepton::differential_transversal_polarisation(const double & s) const
    {
        // cf. [DN:2019A], eq. (4.6)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 2.0 * (a_c.L1cc + 2.0 * a_c.L2cc) / 3.0 / _imp->decay_width(a_c);
    }

    // differential angular coefficients
    double
    LambdaBToLambda1520Dilepton::differential_L_1c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L1c;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_1cc(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L1cc;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_1ss(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L1ss;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_2c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L2c;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_2cc(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L2cc;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_2ss(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L2ss;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_3ss(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L3ss;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_4ss(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L4ss;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_5s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L5s;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_5sc(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L5sc;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_6s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L6s;
    }

    double
    LambdaBToLambda1520Dilepton::differential_L_6sc(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.L6sc;
    }


    double
    LambdaBToLambda1520Dilepton::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        return _imp->decay_width(_imp->integrated_angular_coefficients(s_min, s_max));
    }

    double
    LambdaBToLambda1520Dilepton::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrated_decay_width(s_min, s_max) * _imp->tau() / _imp->hbar();
    }

    double
    LambdaBToLambda1520Dilepton::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        // cf. [DN:2019A], eq. (4.7)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return (a_c.L1c + 2.0 * a_c.L2c) / 2.0 / _imp->decay_width(a_c);
    }

    double
    LambdaBToLambda1520Dilepton::integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [DN:2019A], eq. (4.6)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 1 - 2.0 * (a_c.L1cc + 2.0 * a_c.L2cc) / 3.0 / _imp->decay_width(a_c);
    }

    double
    LambdaBToLambda1520Dilepton::integrated_transversal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [DN:2019A], eq. (4.6)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 2.0 * (a_c.L1cc + 2.0 * a_c.L2cc) / 3.0 / _imp->decay_width(a_c);
    }

    // integrated angular coefficients
    double
    LambdaBToLambda1520Dilepton::integrated_L_1c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L1c;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_1cc(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L1cc;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_1ss(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L1ss;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_2c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L2c;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_2cc(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L2cc;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_2ss(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L2ss;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_3ss(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L3ss;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_4ss(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L4ss;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_5s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L5s;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_5sc(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L5sc;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_6s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L6s;
    }

    double
    LambdaBToLambda1520Dilepton::integrated_L_6sc(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.L6sc;
    }


    const std::string
    LambdaBToLambda1520Dilepton::description = "\
The decay \bar{Lambda_b}->\bar{Lambda}(1520) l^+ l^-, with l=e,mu,tau \
a charged lepton and the \bar{Lambda}(1520) further \
decaying to \bar{N} K. Various theory models can be selected using the \
'tag' option";

    const std::string
    LambdaBToLambda1520Dilepton::kinematics_description_s = "\
The invariant mass of the charged lepton pair in GeV^2.";

    const std::string
    LambdaBToLambda1520Dilepton::kinematics_description_c_theta_l = "\
The cosine of the negatively-charged lepton l^-'s helicity angle theta_l in the l^+l^- rest frame.";

    const std::string
    LambdaBToLambda1520Dilepton::kinematics_description_c_theta_Lstar = "\
The cosine of the nucleon's helicity angle theta_Lstar in the Nbar-K rest frame.";

    const std::string
    LambdaBToLambda1520Dilepton::kinematics_description_phi = "\
The azimuthal angle between the Nbar-K plane and the l^+l^- plane.";

    /*
     * For diagnostic purposes only!
     */
    LambdaBToLambda1520Dilepton::Amplitudes
    LambdaBToLambda1520Dilepton::amplitudes(const double & q2) const
    {
        return _imp->amplitude_generator->amplitudes(q2);
    }

    const std::set<ReferenceName>
    LambdaBToLambda1520Dilepton::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambda1520Dilepton::begin_options()
    {
        return Implementation<LambdaBToLambda1520Dilepton>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambda1520Dilepton::end_options()
    {
        return Implementation<LambdaBToLambda1520Dilepton>::options.cend();
    }
}
