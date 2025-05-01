/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016-2025 Danny van Dyk
 * Copyright (c) 2018      Keri Vos
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

#include <eos/form-factors/analytic-b-to-pi-pi.hh>
#include <eos/form-factors/pi-lcdas.hh>
#include <eos/maths/derivative.hh>
#include <eos/utils/exception.hh>
#include <eos/maths/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/models/model.hh>
#include <eos/maths/polylog.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

#include <functional>

namespace eos
{
    template <>
    struct Implementation<AnalyticFormFactorBToPiPiBFvD2016>
    {
        struct Traces
        {
            double s1, s2, s3, s4, s5, s6, s7, s8;
        };

        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToP>> b_to_pi_ff;

        // hadronic parameters
        UsedParameter m_B;
        UsedParameter f_pi;

        // renormalization scale
        UsedParameter _mu;

        // further hadronic inputs
        PionLCDAs pi;

        // routine to determine renormlization scale
        std::function<double (const double &)> mu;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            b_to_pi_ff(FormFactorFactory<PToP>::create("B->pi::" + o.get("soft-form-factor", "BCL2008"), p, o)),
            m_B(p["mass::B_d"], u),
            f_pi(p["decay-constant::pi"], u),
            _mu(p["B->pipi::mu@BFvD2016"], u),
            pi(p, o)
        {
            std::string scale = o.get("scale", "fixed");

            if ("fixed" == scale)
            {
                mu = std::bind(&Implementation<AnalyticFormFactorBToPiPiBFvD2016>::mu_fixed, this, std::placeholders::_1);
            }
            else if ("variable" == scale)
            {
                mu = std::bind(&Implementation<AnalyticFormFactorBToPiPiBFvD2016>::mu_variable, this, std::placeholders::_1);
            }
            else
            {
                throw InvalidOptionValueError("scale", scale, "fixed, variable");
            }
#if 0
            static const double q2 = 0.6, k2 = 18.6;
            std::cout << "q2 = " << q2 << std::endl;
            std::cout << "k2 = " << k2 << std::endl;

            auto tr = traces_long(q2, k2, 0.0);
            std::cout << "s_1 = " << tr.s1 << std::endl;
            std::cout << "s_2 = " << tr.s2 << std::endl;
            std::cout << "s_3 = " << tr.s3 << std::endl;
            std::cout << "s_4 = " << tr.s4 << std::endl;
            std::cout << "s_5 = " << tr.s5 << std::endl;
            std::cout << "s_6 = " << tr.s6 << std::endl;
            std::cout << "s_7 = " << tr.s7 << std::endl;
            std::cout << "s_8 = " << tr.s8 << std::endl;

            check_thorsten();
#endif
        }

#if 0
        static double gsl_adapter(double * parameters, size_t dim, void * _imp)
        {
            if (dim != 3u)
                throw InternalError("Implemenation<BToPiPiLeptonNeutrino::normalized_differentia_decay_width_gsl_adapter(): wrong number of parameters!");

            const double q2 = parameters[0];
            const double k2 = parameters[1];
            const double z  = parameters[2];

            auto imp = reinterpret_cast<const Implementation<AnalyticFormFactorBToPiPiBFvD2016> *>(_imp);

            auto lambda = eos::lambda(q2, k2, imp->m_B() * imp->m_B());
            auto lambda_pi = eos::lambda(q2, 0.135 * 0.135, imp->m_B() * imp->m_B());
            if (lambda <= 0)
                return 0.0;

            auto tr = imp->traces_long(q2, k2, z);

            const double E1 = imp->energy_1(q2, k2, z), E2 = imp->energy_2(q2, k2, z), m_B = imp->m_B();
            const double c12 = 2.0 * E1 * m_B / k2 - 1.0, c13 = 0.5;
            const double c21 = 1.0, c22 = 1.0, c25 = -m_B / (2.0 * E2), c27 = -0.5;

            complex<double> M_long = complex<double>(1.0, 0.0) * (
                        imp->integral_lo_tw2_f1(q2, k2, z) * (c12 * tr.s2 + c13 * tr.s3)
                    +   imp->integral_lo_tw2_f2(q2, k2, z) * (c21 * tr.s1 + c22 * tr.s2 + c25 * tr.s5 + c27 * tr.s7)
                    );
            return std::norm(M_long / k2) * q2 ;
        }

        inline double R_int_num()
        {
            // Yields a numerical error of approximately 0.2%.
            static const size_t calls = 50000;

            const double x_min[3] = { 0.02, 18.60, -1.0 };
            const double x_max[3] = { 1.00, 24.60, +1.0 };

            gsl_monte_function integrand{ &gsl_adapter, 3u, const_cast<void *>(reinterpret_cast<const void *>(this)) };

            double result, error;

            gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_monte_miser_state * state = gsl_monte_miser_alloc(3u);

            gsl_monte_miser_integrate(&integrand, x_min, x_max, 3u, calls, rng, state, &result, &error);

            gsl_monte_miser_free(state);
            gsl_rng_free(rng);

            return result;
        }

        inline double R_int_denom_integrand(const double & q2)
        {
            return eos::lambda(m_B() * m_B(), 0.135 * 0.135, q2) / power_of<2>(b_to_pi_ff->f_p(q2));
        }

        inline double R_int_denom()
        {
            std::function<double (const double &)> integrand(std::bind(&Implementation<AnalyticFormFactorBToPiPiBFvD2016>::R_int_denom_integrand, *this, std::placeholders::_1));

            return integrate<GSL::QNG>(integrand, 0.02, 0.95);
        }

        inline void check_thorsten()
        {
            std::cout << "R_int = " << R_int_num() / R_int_denom() << std::endl;
        }
#endif

        /*!
         * Returns the renormalization scale value as a parameter value.
         */
        inline double mu_fixed(const double &) const
        {
            return _mu();
        }

        /*!
         * Return the renormalization scale value as a function of k2.
         */
        inline double mu_variable(const double & k2) const
        {
            return _mu() / power_of<2>(m_B()) * k2;
        }

        inline double xi_pi(const double & E2) const
        {
            // qtilde2 is the momentum transfer (squared) in the B->pi(2) system,
            // qtilde = p - k2. Therefore qtilde2 = M_B^2 - 2 E2 M_B.
            const double m_B = this->m_B(), m_B2 = m_B * m_B;
            const double qtilde2 = m_B2 - 2.0 * E2 * m_B;

            return b_to_pi_ff->f_p(qtilde2);
        }

        inline double energy_1(const double & q2, const double & k2, const double & z) const
        {
            const double m_B = this->m_B(), m_B2 = m_B * m_B;
            const double lambda = eos::lambda(m_B2, q2, k2);
            const double sqrt_lambda = std::sqrt(lambda);

            return (m_B2 + k2 - q2 + z * sqrt_lambda) / (4.0 * m_B);
        }

        inline double energy_2(const double & q2, const double & k2, const double & z) const
        {
            const double m_B = this->m_B(), m_B2 = m_B * m_B;
            const double sqrt_lambda = std::sqrt(lambda(m_B2, q2, k2));

            return (m_B2 + k2 - q2 - z * sqrt_lambda) / (4.0 * m_B);
        }

        Traces traces_perp(const double & q2, const double & k2, const double & /*z*/) const
        {
            const double m_B2 = m_B() * m_B();
            const double sqrt_k2 = std::sqrt(k2);
            const double sqrt_lambda = std::sqrt(lambda(m_B2, q2, k2));
            const double s_perp = sqrt_k2 * sqrt_lambda / (2.0 * m_B2);

            return Traces{ 0.0, 0.0, 0.0, 0.0, -s_perp, +s_perp, 0.0, 0.0 };
        }

        Traces traces_para(const double & q2, const double & k2, const double & /*z*/) const
        {
            const double m_B2 = m_B() * m_B();
            const double sqrt_k2 = std::sqrt(k2);
            const double s_para = (m_B2 + k2 - q2) / (2.0 * m_B2) * sqrt_k2;

            return Traces{ sqrt_k2, -sqrt_k2, 2.0 * sqrt_k2, -2.0 * sqrt_k2, +s_para, -s_para, 0.0, 0.0 };
        }

        Traces traces_long(const double & q2, const double & k2, const double & z) const
        {
            const double m_B = this->m_B(), m_B2 = m_B * m_B;
            const double sqrt_lambda = sqrt(lambda(m_B2, q2, k2));
            const double sqrt_q2 = std::sqrt(q2);
            const double a = (m_B2 - k2 - q2) / (2.0 * sqrt_q2);
            const double b = sqrt_lambda / (2.0 * sqrt_q2);
            const double c = k2 * (m_B2 - k2 + q2) / (2.0 * m_B2 * sqrt_q2);
            const double d = k2 / m_B2 * b;

            return Traces{ +a * z + b, -a * z + b, (+a * z + b) * 2.0, (-a * z + b) * 2.0, +c * z + d, -c * z + d, 0.0, 0.0 };
        }

        Traces traces_time(const double & q2, const double & k2, const double & z) const
        {
            const double m_B = this->m_B();
            const double E1 = energy_1(q2, k2, z), E2 = energy_2(q2, k2, z);
            const double sqrt_q2 = std::sqrt(q2);
            const double a = (2.0 * E1 * m_B - k2) / sqrt_q2;
            const double b = (2.0 * E2 * m_B - k2) / sqrt_q2;
            const double c = k2 * (m_B - 2.0 * E2) / (sqrt_q2 * m_B);
            const double d = k2 * (m_B - 2.0 * E1) / (sqrt_q2 * m_B);

            return Traces{ a, b, +2.0 * a, 2.0 * b, c, d, 0.0, 0.0 };
        }

        /* Twist 2, leading order in alpha_s */
        inline double integral_lo_tw2_f1(const double & q2, const double & k2, const double & z) const
        {
            const double mu = this->mu(k2);
            const double r1 = 2.0 * energy_1(q2, k2, z) * m_B() / k2, r1m1 = r1 - 1.0;
            const double r2 = 2.0 * energy_2(q2, k2, z) * m_B() / k2;

            const double L = std::log((r1m1 + r2) / r2);

            const double result_a0 = (3.0 * r1m1 * (r1m1 + 2.0 * r2)
                - 6.0 * r2 * (r1m1 + r2) * L) / power_of<3>(r1m1);
            const double result_a2 = 3.0 *
                (r1m1 * (r1m1 + 2.0 * r2) * (r1m1 * r1m1 + 30.0 * r1m1 * r2 + 30.0 * r2 * r2)
                - 12.0 * r2 * (r1m1 + r2) * (r1m1 * r1m1 + 5.0 * r1m1 * r2 + 5.0 * r2 * r2) * L)
                / power_of<5>(r1m1);

            return result_a0 + pi.a2(mu) * result_a2;
        }

        inline double integral_lo_tw2_f2(const double & q2, const double & k2, const double & z) const
        {
            const double mu = this->mu(k2);
            const double r1 = 2.0 * energy_1(q2, k2, z) * m_B() / k2, r1m1 = r1 - 1.0;
            const double r2 = 2.0 * energy_2(q2, k2, z) * m_B() / k2;

            const double L = std::log((r1m1 + r2) / r2);

            const double result_a0 = 6.0 * r2 * ((r1m1 + r2) * L - r1m1) / power_of<2>(r1m1);
            const double result_a2 = -6.0 * r2 *
                (r1m1 * (16.0 * r1m1 * r1m1 + 45.0 * r1m1 * r2 + 30.0 * r2 * r2)
                 - 6.0 * (r1m1 + r2) * (r1m1 * r1m1 + 5.0 * r1m1 * r2 + 5.0 * r2 * r2) * L)
                / power_of<4>(r1m1);

            return result_a0 + pi.a2(mu) * result_a2;
        }

        /*
         * Using the traces s1 through s8, any of the form factors can be cast
         * in the form given in eq. (3.13), [BFvD2016].
         */
        complex<double> ff_lo_tw2(const Traces & tr, const double & q2, const double & k2, const double & z) const
        {
            static const double CF = 4.0 / 3.0, NC = 3.0;
            const double mu = this->mu(k2);
            const double E1 = energy_1(q2, k2, z), E2 = energy_2(q2, k2, z), m_B = this->m_B();

            const double c12 = 2.0 * E1 * m_B / k2 - 1.0, c13 = 0.5;
            const double c21 = 1.0, c22 = 1.0, c25 = -m_B / (2.0 * E2), c27 = -0.5;

            const double prefactor = 2.0 * M_PI * f_pi() / k2 * xi_pi(E2) * model->alpha_s(mu) * CF / NC;

            return complex<double>(0.0, 1.0) * prefactor * (
                        integral_lo_tw2_f1(q2, k2, z) * (c12 * tr.s2 + c13 * tr.s3)
                    +   integral_lo_tw2_f2(q2, k2, z) * (c21 * tr.s1 + c22 * tr.s2 + c25 * tr.s5 + c27 * tr.s7)
                    );
        }

        /* Twist 3, leading order in alpha_s */
        inline double integral_lo_tw3_sigma1(const double & q2, const double & k2, const double & z) const
        {
            const double mu = this->mu(k2);
            const double r1 = 2.0 * energy_1(q2, k2, z) * m_B() / k2, r1m1 = r1 - 1.0;
            const double r2 = 2.0 * energy_2(q2, k2, z) * m_B() / k2;

            const double L = std::log((r1m1 + r2) / r2);

            const double result = 2.0 * m_B() * pi.mu3(mu) / k2
                * (r1m1 - (r1m1 + r2) * L) / power_of<2>(r1m1);

            return result;
        }

        inline double integral_lo_tw3_sigma2(const double & q2, const double & k2, const double & z) const
        {
            const double mu = this->mu(k2);
            const double r1 = 2.0 * energy_1(q2, k2, z) * m_B() / k2, r1m1 = r1 - 1.0;
            const double r2 = 2.0 * energy_2(q2, k2, z) * m_B() / k2;

            const double L = std::log((r1m1 + r2) / r2);

            const double result = -2.0 * m_B() * pi.mu3(mu) / k2
                * (r1m1 - r2 * L) / power_of<2>(r1m1);

            return result;
        }

        inline double integral_lo_tw3_sigma3(const double & q2, const double & k2, const double & z) const
        {
            return -1.0 * integral_lo_tw3_sigma2(q2, k2, z);
        }

        inline double integral_lo_tw3_sigma4(const double & q2, const double & k2, const double & z) const
        {
            const double r1 = 2.0 * energy_1(q2, k2, z) * m_B() / k2;
            const double r2 = 2.0 * energy_2(q2, k2, z) * m_B() / k2;

            return r2 / (r1 - 1.0) * (integral_lo_tw3_sigma2(q2, k2, z) - integral_lo_tw3_sigma1(q2, k2, z));
        }

        inline double integral_lo_tw3_finite(const double & q2, const double & k2, const double & z) const
        {
            const double mu = this->mu(k2);
            const double r1 = 2.0 * energy_1(q2, k2, z) * m_B() / k2, r1m1 = r1 - 1.0;
            const double r2 = 2.0 * energy_2(q2, k2, z) * m_B() / k2;

            const double L = std::log((r1m1 + r2) / r2);

            const double result = +2.0 * m_B() * pi.mu3(mu) / k2
                * (r1m1 + r2) * L / r1m1;

            return result;
        }

        /*
         * Using the traces s1 through s8, any of the form factors can be cast
         * in the form given in eq. (3.20) and (3.21), [BFvD2016].
         */
        complex<double> ff_lo_tw3(const Traces & tr, const double & q2, const double & k2, const double & z) const
        {
            static const double CF = 4.0 / 3.0, NC = 3.0;
            const double mu = this->mu(k2);
            const double E1 = energy_1(q2, k2, z), E2 = energy_2(q2, k2, z), m_B = this->m_B();

            const double c12 = -1.0, c14 = -E2 / m_B, c16 = 2.0 * E2 * m_B / k2, c17 = 0.5;
            const double c22 = -1.0, c24 = -E2 / m_B, c26 = 2.0 * E2 * m_B / k2, c28 = 0.5;
            const double c33 = E2 / m_B, c35 = -1.0 + (4.0 * E1 * E2 - k2) * m_B / (2.0 * E2 * k2), c37 = 0.5;
            const double c41 = -k2 / (2.0 * E2 * m_B), c43 = -0.5 * c41, c47 = E1 / E2 * 0.5;

            const double prefactor = 2.0 * M_PI * f_pi() / k2 * xi_pi(E2) * model->alpha_s(mu) * CF / NC;

            return complex<double>(0.0, 1.0) * prefactor * (
                        integral_lo_tw3_sigma1(q2, k2, z) * (c12 * tr.s2 + c14 * tr.s4 + c16 * tr.s6 + c17 * tr.s7)
                    +   integral_lo_tw3_sigma2(q2, k2, z) * (c22 * tr.s2 + c24 * tr.s4 + c26 * tr.s6 + c28 * tr.s8)
                    +   integral_lo_tw3_sigma3(q2, k2, z) * (c33 * tr.s3 + c35 * tr.s5 + c37 * tr.s7)
                    +   integral_lo_tw3_sigma4(q2, k2, z) * (c41 * tr.s1 + c43 * tr.s3 + c47 * tr.s7)
                    +   integral_lo_tw3_finite(q2, k2, z) * tr.s5
                    );
        }

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            const double m_B2 = m_B() * m_B();

            // Integral over f_1, cf. [BFvD2016], eq. (3.11)
            {
                results.add(Diagnostics::Entry{ integral_lo_tw2_f1(0.0112245 * m_B2, 0.6666667 * m_B2, -1.0), "I_1(q2: 0.0112245, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw2_f1(0.0112245 * m_B2, 0.6666667 * m_B2,  0.0), "I_1(q2: 0.0112245, k2: 0.6666667, z:  0), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw2_f1(0.0112245 * m_B2, 0.6666667 * m_B2, +1.0), "I_1(q2: 0.0112245, k2: 0.6666667, z: +1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw2_f1(0.0224490 * m_B2, 0.6666667 * m_B2, -1.0), "I_1(q2: 0.0224490, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw2_f1(0.0224490 * m_B2, 0.6666667 * m_B2,  0.0), "I_1(q2: 0.0224490, k2: 0.6666667, z:  0), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw2_f1(0.0224490 * m_B2, 0.6666667 * m_B2, +1.0), "I_1(q2: 0.0224490, k2: 0.6666667, z: +1), [BFvD2016]" });
            }

            // Integral over f_2, cf. [BFvD2016], eq. (3.11)
            {
                results.add(Diagnostics::Entry{ integral_lo_tw2_f2(0.0112245 * m_B2, 0.6666667 * m_B2, -1.0), "I_2(q2: 0.0112245, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw2_f2(0.0112245 * m_B2, 0.6666667 * m_B2,  0.0), "I_2(q2: 0.0112245, k2: 0.6666667, z:  0), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw2_f2(0.0112245 * m_B2, 0.6666667 * m_B2, +1.0), "I_2(q2: 0.0112245, k2: 0.6666667, z: +1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw2_f2(0.0224490 * m_B2, 0.6666667 * m_B2, -1.0), "I_2(q2: 0.0224490, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw2_f2(0.0224490 * m_B2, 0.6666667 * m_B2,  0.0), "I_2(q2: 0.0224490, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw2_f2(0.0224490 * m_B2, 0.6666667 * m_B2, +1.0), "I_2(q2: 0.0224490, k2: 0.6666667, z: -1), [BFvD2016]" });
            }

            // Integral over f_{sigma,1}, cf. [BFvD2016], eq. (3.21)
            {
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma1(0.0112245 * m_B2, 0.6666667 * m_B2, -1.0), "I_{sigma_1}(q2: 0.0112245, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma1(0.0112245 * m_B2, 0.6666667 * m_B2,  0.0), "I_{sigma_1}(q2: 0.0112245, k2: 0.6666667, z:  0), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma1(0.0112245 * m_B2, 0.6666667 * m_B2, +1.0), "I_{sigma_1}(q2: 0.0112245, k2: 0.6666667, z: +1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma1(0.0224490 * m_B2, 0.6666667 * m_B2, -1.0), "I_{sigma_1}(q2: 0.0224490, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma1(0.0224490 * m_B2, 0.6666667 * m_B2,  0.0), "I_{sigma_1}(q2: 0.0224490, k2: 0.6666667, z:  0), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma1(0.0224490 * m_B2, 0.6666667 * m_B2, +1.0), "I_{sigma_1}(q2: 0.0224490, k2: 0.6666667, z: +1), [BFvD2016]" });
            }

            // Integral over f_{sigma,2}, cf. [BFvD2016], eq. (3.21)
            {
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma2(0.0112245 * m_B2, 0.6666667 * m_B2, -1.0), "I_{sigma_2}(q2: 0.0112245, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma2(0.0112245 * m_B2, 0.6666667 * m_B2,  0.0), "I_{sigma_2}(q2: 0.0112245, k2: 0.6666667, z:  0), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma2(0.0112245 * m_B2, 0.6666667 * m_B2, +1.0), "I_{sigma_2}(q2: 0.0112245, k2: 0.6666667, z: +1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma2(0.0224490 * m_B2, 0.6666667 * m_B2, -1.0), "I_{sigma_2}(q2: 0.0224490, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma2(0.0224490 * m_B2, 0.6666667 * m_B2,  0.0), "I_{sigma_2}(q2: 0.0224490, k2: 0.6666667, z:  0), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_sigma2(0.0224490 * m_B2, 0.6666667 * m_B2, +1.0), "I_{sigma_2}(q2: 0.0224490, k2: 0.6666667, z: +1), [BFvD2016]" });
            }

            // Integral over f_{sigma,finite}, cf. [BFvD2016], eq. (3.21)
            {
                results.add(Diagnostics::Entry{ integral_lo_tw3_finite(0.0112245 * m_B2, 0.6666667 * m_B2, -1.0), "I_{sigma,finite}(q2: 0.0112245, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_finite(0.0112245 * m_B2, 0.6666667 * m_B2,  0.0), "I_{sigma,finite}(q2: 0.0112245, k2: 0.6666667, z:  0), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_finite(0.0112245 * m_B2, 0.6666667 * m_B2, +1.0), "I_{sigma,finite}(q2: 0.0112245, k2: 0.6666667, z: +1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_finite(0.0224490 * m_B2, 0.6666667 * m_B2, -1.0), "I_{sigma,finite}(q2: 0.0224490, k2: 0.6666667, z: -1), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_finite(0.0224490 * m_B2, 0.6666667 * m_B2,  0.0), "I_{sigma,finite}(q2: 0.0224490, k2: 0.6666667, z:  0), [BFvD2016]" });
                results.add(Diagnostics::Entry{ integral_lo_tw3_finite(0.0224490 * m_B2, 0.6666667 * m_B2, +1.0), "I_{sigma,finite}(q2: 0.0224490, k2: 0.6666667, z: +1), [BFvD2016]" });
            }

            return results;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<AnalyticFormFactorBToPiPiBFvD2016>::options
    {
    };

    AnalyticFormFactorBToPiPiBFvD2016::AnalyticFormFactorBToPiPiBFvD2016(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<AnalyticFormFactorBToPiPiBFvD2016>(new Implementation<AnalyticFormFactorBToPiPiBFvD2016>(p, o, *this))
    {
    }

    AnalyticFormFactorBToPiPiBFvD2016::~AnalyticFormFactorBToPiPiBFvD2016()
    {
    }

    FormFactors<PToPP> *
    AnalyticFormFactorBToPiPiBFvD2016::make(const Parameters & p, const Options & o)
    {
        return new AnalyticFormFactorBToPiPiBFvD2016(p, o);
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::re_f_perp(const double & q2, const double & k2, const double & z) const
    {
        return std::real(this->f_perp(q2, k2, z));
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::im_f_perp(const double & q2, const double & k2, const double & z) const
    {
        return std::imag(this->f_perp(q2, k2, z));
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::re_f_para(const double & q2, const double & k2, const double & z) const
    {
        return std::real(this->f_para(q2, k2, z));
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::im_f_para(const double & q2, const double & k2, const double & z) const
    {
        return std::imag(this->f_para(q2, k2, z));
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::re_f_long(const double & q2, const double & k2, const double & z) const
    {
        return std::real(this->f_long(q2, k2, z));
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::im_f_long(const double & q2, const double & k2, const double & z) const
    {
        return std::imag(this->f_long(q2, k2, z));
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::re_f_time(const double & q2, const double & k2, const double & z) const
    {
        return std::real(this->f_time(q2, k2, z));
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::im_f_time(const double & q2, const double & k2, const double & z) const
    {
        return std::imag(this->f_time(q2, k2, z));
    }

    complex<double>
    AnalyticFormFactorBToPiPiBFvD2016::f_perp(const double & q2, const double & k2, const double & z) const
    {
        auto tr = _imp->traces_perp(q2, k2, z);

        return _imp->ff_lo_tw2(tr, q2, k2, z) + _imp->ff_lo_tw3(tr, q2, k2, z);
    }

    complex<double>
    AnalyticFormFactorBToPiPiBFvD2016::f_para(const double & q2, const double & k2, const double & z) const
    {
        auto tr = _imp->traces_para(q2, k2, z);

        return _imp->ff_lo_tw2(tr, q2, k2, z) + _imp->ff_lo_tw3(tr, q2, k2, z);
    }

    complex<double>
    AnalyticFormFactorBToPiPiBFvD2016::f_long(const double & q2, const double & k2, const double & z) const
    {
        auto tr = _imp->traces_long(q2, k2, z);

        return _imp->ff_lo_tw2(tr, q2, k2, z) + _imp->ff_lo_tw3(tr, q2, k2, z);
    }

    complex<double>
    AnalyticFormFactorBToPiPiBFvD2016::f_time(const double & q2, const double & k2, const double & z) const
    {
        auto tr = _imp->traces_time(q2, k2, z);

        return _imp->ff_lo_tw2(tr, q2, k2, z) + _imp->ff_lo_tw3(tr, q2, k2, z);
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::f_perp_im_res_qhat2(const double & /*q2*/, const double & /*k2*/) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::f_para_im_res_qhat2(const double & /*q2*/, const double & /*k2*/) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::f_long_im_res_qhat2(const double & /*q2*/, const double & /*k2*/) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    double
    AnalyticFormFactorBToPiPiBFvD2016::f_time_im_res_qhat2(const double & /*q2*/, const double & /*k2*/) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    Diagnostics
    AnalyticFormFactorBToPiPiBFvD2016::diagnostics() const
    {
        return _imp->diagnostics();
    }

    std::vector<OptionSpecification>::const_iterator
    AnalyticFormFactorBToPiPiBFvD2016::begin_options()
    {
        return Implementation<AnalyticFormFactorBToPiPiBFvD2016>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    AnalyticFormFactorBToPiPiBFvD2016::end_options()
    {
        return Implementation<AnalyticFormFactorBToPiPiBFvD2016>::options.cend();
    }

    template <>
    struct Implementation<AnalyticFormFactorBToPiPiFvDV2018>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToP>> b_to_pi_ff;

        // hadronic parameters
        UsedParameter m_B;
        UsedParameter m_Bst;
        UsedParameter g_BstBpi;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            b_to_pi_ff(FormFactorFactory<PToP>::create("B->pi::" + o.get("soft-form-factor", "BCL2008"), p)),
            m_B(p["mass::B_d"], u),
            m_Bst(p["mass::B_d^*"], u),
            g_BstBpi(p["decay-constant::g_{B^*Bpi}"], u)
        {
        }

        inline double xi_pi(const double & q2) const
        {
            return b_to_pi_ff->f_p(q2);
        }

        inline double lambda(const double & q2, const double & k2) const
        {
            return eos::lambda(q2, k2, m_B() * m_B());
        }

        inline double f_perp_im_res_qhat2(const double & q2, const double & k2) const
        {
            const double m_B2 = power_of<2>(this->m_B());
            const double m_Bst2 = power_of<2>(this->m_Bst());
            /* divided by i */
            const double Im_contracted_T_perp = -1.0 * sqrt(k2 * lambda(q2, k2)) * (m_B2 + m_Bst2)
                    / (4.0 * m_B() * m_Bst2);

            return xi_pi(q2) * g_BstBpi() * Im_contracted_T_perp;
        }

        inline double f_para_im_res_qhat2(const double & q2, const double & k2) const
        {
            const double m_B2 = power_of<2>(this->m_B()), m_B4 = power_of<4>(this->m_B());
            const double m_Bst2 = power_of<2>(this->m_Bst());
            /* divided by i */
            const double Im_contracted_T_para = -1.0 * sqrt(k2) * (m_B4 + m_Bst2 * (q2 - k2) + m_B2 * (q2 - 3.0 * m_Bst2 - k2))
                    / (4.0 * m_B() * m_Bst2);

            return xi_pi(q2) * g_BstBpi() * Im_contracted_T_para;
        }

        inline double f_long_im_res_qhat2(const double & q2, const double & k2) const
        {
            const double m_B2 = power_of<2>(this->m_B());
            const double m_Bst2 = power_of<2>(this->m_Bst());
            /* divided by i */
            const double Im_contracted_T_long = -1.0 * (k2 * (m_B2 + m_Bst2) - (m_B2 - q2) * (m_B2 - m_Bst2))
                    * (k2 * m_Bst2 + m_B2 * (q2 - m_Bst2))
                    / (2.0 * m_B() * m_Bst2 * std::sqrt(q2 * lambda(q2, k2)));

            return xi_pi(q2) * g_BstBpi() * Im_contracted_T_long;
        }

        inline double f_time_im_res_qhat2(const double & q2, const double & k2) const
        {
            const double m_B2 = power_of<2>(this->m_B());
            const double m_Bst2 = power_of<2>(this->m_Bst());
            /* divided by i */
            const double Im_contracted_T_time = -1.0 * (m_B2 * (m_B2 - m_Bst2) * (m_Bst2 - q2) - k2 * m_Bst2 * (m_B2 + m_Bst2))
                    / (2.0 * m_B() * m_Bst2 * std::sqrt(q2));

            return xi_pi(q2) * g_BstBpi() * Im_contracted_T_time;
        }

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            return results;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<AnalyticFormFactorBToPiPiFvDV2018>::options
    {
    };

    AnalyticFormFactorBToPiPiFvDV2018::AnalyticFormFactorBToPiPiFvDV2018(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<AnalyticFormFactorBToPiPiFvDV2018>(new Implementation<AnalyticFormFactorBToPiPiFvDV2018>(p, o, *this))
    {
    }

    AnalyticFormFactorBToPiPiFvDV2018::~AnalyticFormFactorBToPiPiFvDV2018()
    {
    }

    FormFactors<PToPP> *
    AnalyticFormFactorBToPiPiFvDV2018::make(const Parameters & p, const Options & o)
    {
        return new AnalyticFormFactorBToPiPiFvDV2018(p, o);
    }

    complex<double>
    AnalyticFormFactorBToPiPiFvDV2018::f_perp(const double & /*q2*/, const double & /*k2*/, const double & /*z*/) const
    {
        throw InternalError("Not available");
    }

    complex<double>
    AnalyticFormFactorBToPiPiFvDV2018::f_para(const double & /*q2*/, const double & /*k2*/, const double & /*z*/) const
    {
        throw InternalError("Not available");
    }

    complex<double>
    AnalyticFormFactorBToPiPiFvDV2018::f_long(const double & /*q2*/, const double & /*k2*/, const double & /*z*/) const
    {
        throw InternalError("Not available");
    }

    complex<double>
    AnalyticFormFactorBToPiPiFvDV2018::f_time(const double & /*q2*/, const double & /*k2*/, const double & /*z*/) const
    {
        throw InternalError("Not available");
    }

    double
    AnalyticFormFactorBToPiPiFvDV2018::f_perp_im_res_qhat2(const double & q2, const double & k2) const
    {
        return _imp->f_perp_im_res_qhat2(q2, k2);
    }

    double
    AnalyticFormFactorBToPiPiFvDV2018::f_para_im_res_qhat2(const double & q2, const double & k2) const
    {
        return _imp->f_para_im_res_qhat2(q2, k2);
    }

    double
    AnalyticFormFactorBToPiPiFvDV2018::f_long_im_res_qhat2(const double & q2, const double & k2) const
    {
        return _imp->f_long_im_res_qhat2(q2, k2);
    }

    double
    AnalyticFormFactorBToPiPiFvDV2018::f_time_im_res_qhat2(const double & q2, const double & k2) const
    {
        return _imp->f_time_im_res_qhat2(q2, k2);
    }

    Diagnostics
    AnalyticFormFactorBToPiPiFvDV2018::diagnostics() const
    {
        return _imp->diagnostics();
    }

    std::vector<OptionSpecification>::const_iterator
    AnalyticFormFactorBToPiPiFvDV2018::begin_options()
    {
        return Implementation<AnalyticFormFactorBToPiPiFvDV2018>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    AnalyticFormFactorBToPiPiFvDV2018::end_options()
    {
        return Implementation<AnalyticFormFactorBToPiPiFvDV2018>::options.cend();
    }
}
