/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Stefan Meiser
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

#include <eos/b-decays/bq-to-dq-vec.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/vec-lcdas.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/gegenbauer-polynomial.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/maths/polylog.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    using std::norm;
    using std::real;
    using std::imag;

    /*
     * Decay: B_q -> D_q V, cf. [BBNS:2000A] (class I only, V = rho^- or Kstar^-)
     */
    template <>
    struct Implementation<BqToDqVector>
    {
        SpecifiedOption opt_model;

        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter g_fermi;

        QuarkFlavorOption opt_q;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_D;

        UsedParameter m_V;

        UsedParameter f_V;

        UsedParameter alpha_s;

        std::function<double ()> mu_P;

        std::shared_ptr<UsedParameter> ff_f_p;

        std::shared_ptr<VectorLCDAs> lcdas;

        SpecifiedOption opt_cp_conjugate;

        bool cp_conjugate;

        UsedParameter mu;

        SpecifiedOption opt_accuracy;
        double switch_lo;
        double switch_nlo;
        double switch_nlp;

        std::function<complex<double> ()> ckm_factor;
        std::function<WilsonCoefficients<bern::ClassIII> (bool)> wc;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_model(o, options, "model"),
            model(Model::make(opt_model.value(), p, o)),
            opt_q(o, options, "q"),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_B(p["mass::B_" + opt_q.str()], u),
            m_D(p["mass::D_" + opt_q.str()], u),
            m_V(p["mass::" + stringify(opt_q.value() == QuarkFlavor::down ? "K_u^*" : "rho^+")], u),
            f_V(p["decay-constant::" + stringify(opt_q.value() == QuarkFlavor::down ? "K_u^*" : "rho")], u),
            alpha_s(p["QCD::alpha_s(MZ)"], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            opt_cp_conjugate(o, options, "cp-conjugate"),
            cp_conjugate(destringify<bool>(opt_cp_conjugate.value())),
            mu(p[stringify(opt_q.value() == QuarkFlavor::down ? "s" : "d") + "bcu::mu"], u),
            opt_accuracy(o, options, "accuracy")
        {
            Context ctx("When constructing B_q->D_q V observable");

            // handle the spectator quark flavor q
            switch (opt_q.value())
            {
                case QuarkFlavor::strange:
                    ckm_factor = [this]() { return conj(model->ckm_ud()) * model->ckm_cb(); };
                    wc         = [this](const bool & cp_conjugate) { return model->wet_dbcu(cp_conjugate); };
                    ff_f_p     = std::make_shared<UsedParameter>(p["B_s->D_srho::f_p(Mrho2)"], u);
                    lcdas      = VectorLCDAs::make("rho", p, o);
                    break;
                case QuarkFlavor::down:
                    ckm_factor = [this]() { return conj(model->ckm_us()) * model->ckm_cb(); };
                    wc         = [this](const bool & cp_conjugate) { return model->wet_sbcu(cp_conjugate); };
                    ff_f_p     = std::make_shared<UsedParameter>(p["B->DK^*::f_p(MKstar2)"], u);
                    lcdas      = VectorLCDAs::make("Kstar", p, o);
                    break;
                default:
                    throw InternalError("Invalid quark flavor: " + stringify(opt_q.value()));
            }

            if (opt_accuracy.value() == "LO")
            {
                switch_lo = 1.0;
                switch_nlo = 0.0;
                switch_nlp = 0.0;
            }
            else if (opt_accuracy.value() == "NLO")
            {
                switch_lo = 0.0;
                switch_nlo = 1.0;
                switch_nlp = 0.0;
            }
            if (opt_accuracy.value() == "LO+NLO")
            {
                switch_lo = 1.0;
                switch_nlo = 1.0;
                switch_nlp = 0.0;
            }
            else if (opt_accuracy.value() == "NLP")
            {
                switch_lo = 0.0;
                switch_nlo = 0.0;
                switch_nlp = 1.0;
            }
            else if (opt_accuracy.value() == "all")
            {
                switch_lo = 1.0;
                switch_nlo = 1.0;
                switch_nlp = 1.0;
            }

            u.uses(*model);
        }

        // auxiliary functions for hard-scattering kernels
        complex<double> AVLL(const double & u, const double & z) const
        {
            // evaluate AVLL chosing the correct branch for physical z = m_c / m_b < 1
            if((0.0 < z) && (z < 1.0))
            {
                return
                    log(u * (1.0 - z * z)) / (1.0 - u * (1.0 - z * z)) - pow(log(u * (1.0 - z * z)), 2) - dilog(1.0 - u * (1.0 - z * z));
            }
            // evaluate AVLL chosing the correct branch for unphysical z = m_b / m_c > 1
            else if (z > 1.0)
            {
                return
                    (log(-(u * (1.0 - z * z))) - 1.0i * M_PI) / (1.0 - u * (1.0 - z * z)) +
                    pow(log(1 - u * (1.0 - (z * z))), 2) / 2.0 - pow(log(-(u * (1.0 - z * z))) - 1.0i * M_PI, 2) +
                    dilog(1.0 / (1.0 - u * (1.0 - z * z))) +
                    -1.0 / 3.0 * M_PI * M_PI - 1.0i * M_PI * log(1 - u * (1.0 - z * z));
            }

            throw InternalError("Invalid value for z: " + stringify(z) + " in AVLL.");
            return 0.0;
        };

        complex<double> fVLL(const double & u, const double & z) const
        {
            // evaluate fVLL chosing the correct branch for physical z = m_c / m_b < 1
            if ((0.0 < z) && (z < 1.0))
            {
                return
                    2.0 * (AVLL(u, z) - AVLL(1.0 - u, z)) - (z / (1.0 - u * (1.0 - z * z))) -
                    (u * (1.0 - z * z) * (z + 3.0 * (1.0 - u * (1.0 - z * z))) * log(u * (1.0 - z * z))) / pow(1.0 - u * (1.0 - z * z), 2);
            }
            // evaluate fVLL chosing the correct branch for unphysical z = m_b / m_c > 1
            else if (z > 1.0)
            {
                return
                    2.0 * (AVLL(u, z) - AVLL(1.0 - u, z)) - (z / (1.0 - u * (1.0 - z * z))) -
                    (u * (1.0 - z * z) * (z + 3.0 * (1.0 - u * (1.0 - z * z))) * (log(-u * (1.0 - z * z)) - 1.0i * M_PI)) / pow(1.0 - u * (1.0 - z * z), 2);
            }

            throw InternalError("Invalid value for z: " + stringify(z) + " in fVLL.");
            return 0.0;
        };

        complex<double> ASLR(const double & u, const double & z) const
        {
            // evaluate ASLR chosing the correct branch for physical z = m_c / m_b < 1
            if ((0.0 < z) && (z < 1.0))
            {
                return
                    z * z / (pow(1.0 + z, 2) * (1.0 - u * (1.0 - z * z))) +
                    ((-2.0 + u * u * pow(-1.0 + z, 2) * (2.0 + 4.0 * z + 3.0 * z * z)) * log(u * (1.0 - z * z))) / pow(1.0 - u * (1 - z * z), 2) +
                    2.0 * ((2.0 * log(u * (1.0 - z * z))) / (1.0 - u * (1 - z * z)) - pow(log(u * (1.0 - z * z)), 2) - dilog(1.0 - u * (1.0 - z * z)));
            }
            // evaluate ASLR chosing the correct branch for unphysical z = m_b / m_c > 1
            else if (z > 1.0)
            {
                return
                    z * z / (pow(1.0 + z, 2) * (1.0 - u * (1.0 - z * z))) +
                    ((-2.0 + u * u * pow(-1.0 + z, 2) * (2.0 + 4.0 * z + 3.0 * z * z)) * (-1.0i * M_PI + log(-(u * (1 - z * z))))) / pow(1.0 - u * (1 - z * z), 2) +
                    2.0 * (-1.0 / 3.0 * M_PI * M_PI + (2.0 * (-1.0i * M_PI + log(-(u * (1.0 - z * z))))) / (1.0 - u * (1 - z * z)) -
                    pow(-1.0i * M_PI + log(-(u * (1.0 - z * z))), 2) - 1.0i * M_PI * log(1.0 - u * (1.0 - z * z)) + pow(log(1.0 - u * (1.0 - z * z)), 2) / 2.0 +
                    dilog(1.0 / (1.0 - u * (1.0 - z * z))));
            }

            throw InternalError("Invalid value for z: " + stringify(z) + " in ASLR.");
            return 0.0;
        };

        complex<double> fSLR(const double & u, const double & z) const
        {
            return
                ASLR(u, z) - ASLR(1.0 - u, z);
        };

        complex<double> ASLL(const double & u, const double & z) const
        {
            // evaluate ASLL chosing the correct branch for physical z = m_c / m_b < 1
            if ((0.0 < z) && (z < 1.0))
            {
                return
                    -2.0 * ((5.0 * u) / (1.0 + z) + ((-1.0 + u * (1.0 - z) * (u * (1.0 - z) + 2.0 * z)) * log(u * (1.0 - z * z))) / (1.0 - u * (1.0 - z * z)) +
                    pow(log(u * (1.0 - z * z)), 2) + dilog(1.0 - u * (1.0 - z * z)));
            }
            // evaluate ASLL chosing the correct branch for unphysical z = m_b / m_c > 1
            else if (z > 1.0)
            {
                return
                    -2.0 * ((5.0 * u) / (1.0 + z) + ((-1.0 + u * (1.0 - z) * (u * (1.0 - z) + 2.0 * z)) * (-1.0i * M_PI + log(-(u * (1.0 - z * z))))) / (1.0 - u * (1.0 - z * z)) +
                    pow(-1.0i * M_PI + log(-(u * (1.0 - z * z))), 2) - 1.0 / 6.0 * M_PI * M_PI - pow(-1.0i * M_PI + log(1.0 - u * (1.0 - z * z)), 2) / 2.0 -
                    dilog(1.0 / (1.0 - u * (1.0 - z * z))));
            }

            throw InternalError("Invalid value for z: " + stringify(z) + " in ASLL.");
            return 0.0;
        };

        complex<double> fSLL(const double & u, const double & z) const
        {
            return
                ASLL(u, z) - ASLL(1.0 - u, z);
        };

        complex<double> ATLL(const double & u, const double & z) const
        {
            // evaluate ATLL chosing the correct branch for physical z = m_c / m_b < 1
            if ((0.0 < z) && (z < 1.0))
            {
                return
                    ((-1.0 + u * (2.0 - u - 2.0 * z  + (-2.0 + u) * z * z)) * log(u * (1.0 - z * z))) / (1.0 - u * (1.0 - z * z)) +
                    (1.0 - 2.0 * u) * (pow(log(u * (1.0 - z * z)), 2) + dilog(1.0 - u * (1.0 - z * z)));
            }
            // evaluate ATLL chosing the correct branch for unphysical z = m_b / m_c > 1
            else if (z > 1.0)
            {
                return
                    ((-1.0 + u * (2.0 - u - 2.0 * z + (-2.0 + u) * z * z)) * (-1.0i * M_PI + log(-(u * (1.0 - z * z))))) / (1.0 - u * (1.0 - z * z)) +
                    (1.0 - 2.0 * u) * (M_PI * M_PI / 3.0 + pow(-1.0i * M_PI + log(-(u * (1.0 - z * z))), 2) + 1.0i * M_PI * log(1.0 - u * (1.0 - z * z)) -
                    pow(log(1.0 - u * (1.0 - z * z)), 2) / 2.0 - dilog(1.0 / (1.0 - u * (1.0 - z * z))));
            }

            throw InternalError("Invalid value for z: " + stringify(z) + " in ATLL.");
            return 0.0;
        };

        complex<double> fTLL(const double & u, const double & z) const
        {
            return
                -((8.0 * (4.0 * u + 3.0)) / (1.0 + z)) + (8.0 * (1.0 - z)) / (1.0 + z) * (ATLL(u, z) + ATLL(1.0 - u, z));
        };

        complex<double> phiV(const double & u, const double & mu) const
        {
            // Legendre polynomials P_n
            static const GegenbauerPolynomial lp_1(1, 1.0 / 2.0);
            static const GegenbauerPolynomial lp_2(2, 1.0 / 2.0);
            static const GegenbauerPolynomial lp_3(3, 1.0 / 2.0);
            static const GegenbauerPolynomial lp_4(4, 1.0 / 2.0);
            static const GegenbauerPolynomial lp_5(5, 1.0 / 2.0);

            const double x = 2.0 * u - 1.0;
            const double P1 = lp_1.evaluate(x);
            const double P2 = lp_2.evaluate(x);
            const double P3 = lp_3.evaluate(x);
            const double P4 = lp_4.evaluate(x);
            const double P5 = lp_5.evaluate(x);

            return 3.0 * (P1 + lcdas->a1perp(mu) * P2 + lcdas->a2perp(mu) * P3 + lcdas->a3perp(mu) * P4 + lcdas->a4perp(mu) * P5);
        };

        complex<double> a_1() const
        {
            const WilsonCoefficients<bern::ClassIII> wc = this->wc(cp_conjugate);

            // cf. [BBNS:2000A], converted to the Bern basis
            const double mb = model->m_b_msbar(mu());
            const double mc = model->m_c_msbar(mu());
            const double z = mc / mb;

            const double mu_L = m_V * lcdas->fperp(mu()) / f_V;

            const double a_s_mu = model->alpha_s(mu()) / (4.0 * M_PI);

            const complex<double> a_1_lo =
                -1.0 / 3.0 * (wc.c1() + wc.c1p())  - 4.0 / 9.0 * (wc.c2() + wc.c2p()) -
                16.0 / 3.0 * (wc.c3() + wc.c3p()) - 64.0 / 9.0 * (wc.c4() + wc.c4p()) +
                1.0 / 6.0 * (wc.c5() + wc.c5p()) + 2.0 / 9.0 * (wc.c6() + wc.c6p()) +
                8.0 / 3.0 * (wc.c9() + wc.c9p()) + 32.0 / 9.0 * (wc.c10() + wc.c10p());

            // compute convolutions
            auto integrand_VLL = [&](const double & u) -> complex<double>
            {
                return (-18.0 - 6.0 * 2.0 * log(mu() / mb) + fVLL(1.0 - u, 1.0 / z) + fVLL(u, z) + (3.0 + 2.0 * log(u / (1.0 - u))) * log(z * z)) * this->lcdas->phipara(u, mu());
            };

            auto integrand_VLR = [&](const double & u) -> complex<double>
            {
                return (6.0 + 6.0 * 2.0 * log(mu() / mb) - (3.0 + 2.0 * log((1.0 - u) / u)) * log(z * z) - fVLL((1.0 - u), z) - fVLL(u, 1.0 / z)) * this->lcdas->phipara(u, mu());
            };

            auto integrand_SLL = [&](const double & u) -> complex<double>
            {
                return (-4.0 * (2.0 * u - 1.0) * (1.0 - z) / (1.0 + z) * 2.0 * log(mu() / mb) + 2.0 * ((2.0 * u - 1.0) * (1.0 - z) / (1.0 + z) + log(u / (1.0 - u))) * log(z * z) + fSLL(u, z) + fSLL((1.0 - u), 1.0 / z)) * phiV(u, mu());
            };

            auto integrand_SLR = [&](const double & u) -> complex<double>
            {
                return (2.0 * log(u / (1.0 - u)) * log(z * z) - 6.0 + fSLR(u, z) + fSLR((1.0 - u), 1.0 / z)) * phiV(u, mu());
            };

            auto integrand_TLL = [&](const double & u) -> complex<double>
            {
                return (-48.0 * 2.0 * log(mu() / mb) + 8.0 * (3.0 + ((u - (1.0 - u)) * (1.0 - z)) / (1.0 + z) * log(u / (1.0 - u))) * log(z * z) + fTLL(u, z) + fTLL(1.0 - u, 1.0 / z)) * phiV(u, mu());
            };

            auto re_integrand_VLL = [this, & integrand_VLL](const double & u) -> const double
            {
                return real(integrand_VLL(u));
            };

            auto im_integrand_VLL = [this, & integrand_VLL](const double & u) -> const double
            {
                return imag(integrand_VLL(u));
            };

            auto re_integrand_VLR = [this, & integrand_VLR](const double & u) -> const double
            {
                return real(integrand_VLR(u));
            };

            auto im_integrand_VLR = [this, & integrand_VLR](const double & u) -> const double
            {
                return imag(integrand_VLR(u));
            };

            auto re_integrand_SLL = [this, & integrand_SLL](const double & u) -> const double
            {
                return real(integrand_SLL(u));
            };

            auto im_integrand_SLL = [this, & integrand_SLL](const double & u) -> const double
            {
                return imag(integrand_SLL(u));
            };

            auto re_integrand_SLR = [this, & integrand_SLR](const double & u) -> const double
            {
                return real(integrand_SLR(u));
            };

            auto im_integrand_SLR = [this, & integrand_SLR](const double & u) -> const double
            {
                return imag(integrand_SLR(u));
            };

            auto re_integrand_TLL = [this, & integrand_TLL](const double & u) -> const double
            {
                return real(integrand_TLL(u));
            };

            const double TVLL_re = integrate<GSL::QAGS>(re_integrand_VLL, 0.0, 1.0);
            const double TVLL_im = integrate<GSL::QAGS>(im_integrand_VLL, 0.0, 1.0);

            const double TVLR_re = integrate<GSL::QAGS>(re_integrand_VLR, 0.0, 1.0);
            const double TVLR_im = integrate<GSL::QAGS>(im_integrand_VLR, 0.0, 1.0);

            const double TSLL_re = integrate<GSL::QAGS>(re_integrand_SLL, 0.0, 1.0);
            const double TSLL_im = integrate<GSL::QAGS>(im_integrand_SLL, 0.0, 1.0);

            const double TSLR_re = integrate<GSL::QAGS>(re_integrand_SLR, 0.0, 1.0);
            const double TSLR_im = integrate<GSL::QAGS>(im_integrand_SLR, 0.0, 1.0);

            const double TTLL_re = integrate<GSL::QAGS>(re_integrand_TLL, 0.0, 1.0);
            const double TTLL_im = 0.0;

            const complex<double> TVLL = TVLL_re + TVLL_im * 1.0i;
            const complex<double> TVLR = TVLR_re + TVLR_im * 1.0i;
            const complex<double> TSLL = TSLL_re + TSLL_im * 1.0i;
            const complex<double> TSLR = TSLR_re + TSLR_im * 1.0i;
            const complex<double> TTLL = TTLL_re + TTLL_im * 1.0i;

            const complex<double> a_1_nlo =
                    4.0 / 9.0  * (wc.c1() + wc.c1p()) * (-((2.0 * mu_L * TSLR) / (mb - mc)) - TVLL + 4.0) +
                    2.0 / 27.0 * (wc.c2() + wc.c2p()) * ((2.0 * mu_L * TSLR) / (mb - mc) + TVLL + 14.0) +
                    32.0 / 9.0 * (wc.c3() + wc.c3p()) * (-((mu_L * TSLR) / (mb - mc)) - 2.0 * (TVLL + 5.0)) +
                    16.0 / 27.0 * (wc.c4() + wc.c4p()) * ((mu_L * TSLR) / (mb - mc) + 2.0 * TVLL + 19.0) +
                    1.0 / 18.0 * (wc.c5() + wc.c5p()) * ((mu_L * (TTLL - 4.0 * TSLL)) / (mb - mc) + 4.0 * (TVLR + 6.0)) +
                    1.0 / 108.0 * (wc.c6() + wc.c6p()) * ((mu_L * (4.0 * TSLL - TTLL)) / (mb - mc) - 4.0 * TVLR + 84.0) +
                    1.0 / 9.0 * (wc.c7() + wc.c7p()) * (32.0 - (2.0 * mu_L * (12.0 * TSLL + TTLL)) / (mb - mc)) +
                    1.0 / 27.0 * (wc.c8() + wc.c8p()) * ((mu_L * (12.0 * TSLL + TTLL)) / (mb - mc) + 56.0) +
                    32.0 / 9.0 * (wc.c9() + wc.c9p()) * ((2.0 * mu_L * (4.0 * TSLL + TTLL)) / (mb - mc) + TVLR - 40.0) +
                    16.0 / 27.0 * (wc.c10() + wc.c10p()) * (-((2.0 * mu_L * (4.0 * TSLL + TTLL)) / (mb - mc)) - TVLR + 76.0);
            return switch_lo * a_1_lo + switch_nlo * a_s_mu * a_1_nlo;
        };


        double decay_width() const
        {
            // cf. [BBNS:2000A], eq. (212), p. 80
            const complex<double> amplitude = g_fermi() / sqrt(2.0) * ckm_factor() * f_V() * ff_f_p->evaluate()
                * sqrt(lambda(m_B * m_B, m_D * m_D, m_V * m_V)) * this->a_1();
            // cf. [BBNS:2000A], eq. (216), p. 80
            const double breakup_momentum = sqrt(lambda(m_B * m_B, m_D * m_D, m_V * m_V)) / (2.0 * m_B);

            // cf. [BBNS:2000A], eq. (221), p. 81
            return norm(amplitude) * breakup_momentum / (8.0 * M_PI * m_B * m_B);
        }

        double branching_ratio() const
        {
            return decay_width() * tau_B / hbar;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BqToDqVector>::options
    {
        Model::option_specification(),
        { "accuracy",     { "LO", "NLO", "NLP", "LO+NLO", "all" }, "all" },
        { "cp-conjugate", { "true", "false" },  "false"                  },
        { "q",            { "s", "d" }                                   },
        { "P",            { "Kstar", "rho"}                              }
    };

    BqToDqVector::BqToDqVector(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BqToDqVector>(new Implementation<BqToDqVector>(parameters, options, *this))
    {
    }

    BqToDqVector::~BqToDqVector()
    {
    }

    double
    BqToDqVector::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    BqToDqVector::decay_width() const
    {
        return _imp->decay_width();
    }

    double
    BqToDqVector::re_a_1() const
    {
        return real(_imp->a_1());
    }

    double
    BqToDqVector::im_a_1() const
    {
        return imag(_imp->a_1());
    }

    const std::set<ReferenceName>
    BqToDqVector::references
    {
        "BBNS:2000A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BqToDqVector::begin_options()
    {
        return Implementation<BqToDqVector>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BqToDqVector::end_options()
    {
        return Implementation<BqToDqVector>::options.cend();
    }
}
