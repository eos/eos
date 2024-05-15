/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023-2024 Danny van Dyk
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

#include <eos/b-decays/bq-to-dstarq-psd.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/k-lcdas.hh>
#include <eos/form-factors/pi-lcdas.hh>
#include <eos/maths/complex.hh>
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
     * Decay: B_q -> D*_q P, cf. [BBNS:2000A] (class I only, P = pi^- or K^-)
     */
    template <>
    struct Implementation<BqToDstarqPseudoscalar>
    {
        SpecifiedOption opt_model;

        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter g_fermi;

        QuarkFlavorOption opt_q;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_Dstar;

        UsedParameter m_P;

        UsedParameter f_P;

        std::function<double ()> mu_P;

        std::shared_ptr<UsedParameter> ff_a_0;

        std::shared_ptr<PseudoscalarLCDAs> lcdas;

        BooleanOption opt_cp_conjugate;

        UsedParameter mu;

        SpecifiedOption opt_accuracy;
        double switch_lo;
        double switch_nlo;
        double switch_nlp;

        std::function<complex<double> ()> ckm_factor;
        std::function<WilsonCoefficients<bern::ClassIII> ()> wc;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_model(o, options, "model"),
            model(Model::make(opt_model.value(), p, o)),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            opt_q(o, options, "q"),
            m_B(p["mass::B_" + opt_q.str()], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            m_Dstar(p["mass::D_" + opt_q.str() + "^*"], u),
            m_P(p["mass::" + stringify(opt_q.value() == QuarkFlavor::down ? "K_u" : "pi^+")], u),
            f_P(p["decay-constant::" + stringify(opt_q.value() == QuarkFlavor::down ? "K_u" : "pi")], u),
            opt_cp_conjugate(o, options, "cp-conjugate"),
            mu(p[stringify(opt_q.value() == QuarkFlavor::down ? "s" : "d") + "bcu::mu"], u),
            opt_accuracy(o, options, "accuracy")
        {
            Context ctx("When constructing B_q->Dstar_q P observable");

            // handle the spectator quark flavor q
            switch (opt_q.value())
            {
                case QuarkFlavor::strange:
                    ckm_factor = [this]() { return conj(model->ckm_ud()) * model->ckm_cb(); };
                    wc         = [this]() { return model->wet_dbcu(opt_cp_conjugate.value()); };
                    ff_a_0     = std::make_shared<UsedParameter>(p["B_s->D_s^*pi::A_0(Mpi2)"], u);
                    lcdas      = PseudoscalarLCDAs::make("pi", p, o);
                    break;
                case QuarkFlavor::down:
                    ckm_factor = [this]() { return conj(model->ckm_us()) * model->ckm_cb(); };
                    wc         = [this]() { return model->wet_sbcu(opt_cp_conjugate.value()); };
                    ff_a_0     = std::make_shared<UsedParameter>(p["B->D^*K::A_0(MK2)"], u);
                    lcdas      = PseudoscalarLCDAs::make("Kbar", p, o);
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
            // evaluate AVLL chosing the correct branch for z = m_c / m_b < 1 and z = - m_c / m_b > -1
            if((-1.0 < z) && (z < 1.0))
            {
                return
                    log(u * (1.0 - z * z)) / (1.0 - u * (1.0 - z * z)) - pow(log(u * (1.0 - z * z)), 2) - dilog(1.0 - u * (1.0 - z * z));
            }
            // evaluate AVLL chosing the correct branch for z = m_b / m_c > 1 and z = - m_b / m_c < -1
            else if ((z > 1.0) || (z < -1.0))
            {
                return
                    (log(-(u * (1.0 - z * z))) - 1.0i * M_PI) / (1.0 - u * (1.0 - z * z)) +
                    pow(log(1.0 - u * (1.0 - (z * z))), 2) / 2.0 - pow(log(-(u * (1.0 - z * z))) - 1.0i * M_PI, 2) +
                    dilog(1.0 / (1.0 - u * (1.0 - z * z))) +
                    -1.0 / 3.0 * M_PI * M_PI - 1.0i * M_PI * log(1.0 - u * (1.0 - z * z));
            }

            throw InternalError("Invalid value for z: " + stringify(z) + " in AVLL.");
            return 0.0;
        };

        complex<double> fVLL(const double & u, const double & z) const
        {
            // evaluate fVLL chosing the correct branch for z = m_c / m_b < 1 and z = - m_c / m_b > -1
            if ((-1.0 < z) && (z < 1.0))
            {
                return
                    2.0 * (AVLL(u, z) - AVLL(1.0 - u, z)) - (z / (1.0 - u * (1.0 - z * z))) -
                    (u * (1.0 - z * z) * (z + 3.0 * (1.0 - u * (1.0 - z * z))) * log(u * (1.0 - z * z))) / pow(1.0 - u * (1.0 - z * z), 2);
            }
            // evaluate fVLL chosing the correct branch for z = m_b / m_c > 1 and z = - m_b / m_c < -1
            else if ((z > 1.0) || (z < -1.0))
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
             // evaluate ASLR chosing the correct branch for z = m_c / m_b < 1 and z = - m_c / m_b > -1
            if ((-1.0 < z) && (z < 1.0))
            {
                return
                    z * z / (pow(1.0 + z, 2) * (1.0 - u * (1.0 - z * z))) +
                    ((-2.0 + u * u * pow(-1.0 + z, 2) * (2.0 + 4.0 * z + 3.0 * z * z)) * log(u * (1.0 - z * z))) / pow(1.0 - u * (1.0 - z * z), 2) +
                    2.0 * ((2.0 * log(u * (1.0 - z * z))) / (1.0 - u * (1.0 - z * z)) - pow(log(u * (1.0 - z * z)), 2) - dilog(1.0 - u * (1.0 - z * z)));
            }
            // evaluate ASLR chosing the correct branch for z = m_b / m_c > 1 and z = - m_b / m_c < -1
            if ((z > 1.0) || (z < -1.0))
            {
                return
                    z * z / (pow(1.0 + z, 2) * (1.0 - u * (1.0 - z * z))) +
                    ((-2.0 + u * u * pow(-1.0 + z, 2) * (2.0 + 4.0 * z + 3.0 * z * z))* (-1.0i * M_PI + log(-(u * (1.0 - z * z))))) / pow(1.0 - u * (1.0 - z * z), 2) +
                    2.0 * (-1.0 / 3.0 * M_PI * M_PI + (2.0 * (-1.0i * M_PI + log(-(u * (1.0 - z * z))))) / (1.0 - u * (1.0 - z * z)) -
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

        complex<double> ATLL(const double & u, const double & z) const
        {
            // evaluate ATLL chosing the correct branch for z = m_c / m_b < 1 and z = - m_c / m_b > -1
            if ((-1.0 < z) && (z < 1.0))
            {
                return
                    ((-1.0 + u * (2.0 - u - 2.0 * z  + (-2.0 + u) * z * z)) * log(u * (1.0 - z * z))) / (1.0 - u * (1.0 - z * z)) +
                    (1.0 - 2.0 * u) * (pow(log(u * (1.0 - z * z)), 2) + dilog(1.0 - u * (1.0 - z * z)));
            }
            // evaluate ATLL chosing the correct branch for z = m_b / m_c > 1 and z = - m_b / m_c < -1
            else if ((z > 1.0) || (z < -1.0))
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

        complex<double> a_1() const
        {
            const WilsonCoefficients<bern::ClassIII> wc = this->wc();

            // cf. [BBNS:2000A], converted to the Bern basis
            const double mb = model->m_b_msbar(mu());
            const double mc = model->m_c_msbar(mu());
            const double z = mc / mb;

            const double mu_L = lcdas->mu3(mu());
            const double f_3P = lcdas->f3(mu());

            const double a_s_mu = model->alpha_s(mu()) / (4.0 * M_PI);

            const complex<double> a_1_lo =
                (wc.c1() + wc.c1p()) * (-1.0 / 3.0 + (2.0 * mu_L) / (3.0 * (mb + mc))) +
                (wc.c2() + wc.c2p()) * (-4.0 / 9.0 + (8.0 * mu_L) / (9.0 * (mb + mc))) +
                ( 8.0 * (wc.c3() + wc.c3p()) * (-2.0 + mu_L / (mb + mc))) / 3.0 +
                (32.0 * (wc.c4() + wc.c4p()) * (-2.0 + mu_L / (mb + mc))) / 9.0 -
                ((wc.c5() + wc.c5p()) * (mb + mc + mu_L)) / (6.0 * (mb + mc)) -
                ( 2.0 * (wc.c6() + wc.c6p()) * (mb + mc + mu_L)) / (9.0 * (mb + mc)) +
                (-2.0 * (wc.c7() + wc.c7p()) * mu_L) / (mb + mc) -
                ( 8.0 * (wc.c8() + wc.c8p()) * mu_L) / (3.0 * (mb + mc)) +
                ( 8.0 * (wc.c9() + wc.c9p()) * (-1.0 + (8.0 * mu_L) / (mb + mc))) / 3.0 +
                (32.0 * (wc.c10() + wc.c10p()) * (-1.0 + (8.0 * mu_L) / (mb + mc))) / 9.0;

            auto a_1_nlo_integrand = [&](const double & u) -> complex<double>
            {
                static const double eps = 1.0e-10;
                complex<double> TVLL ;
                if ((u < eps) || (u > 1.0 - eps))
                {
                    TVLL = 0.0;
                }
                else
                {
                    TVLL = (-18.0 - 6.0 * 2.0 * log(mu() / mb) + fVLL(1.0 - u, -1.0 / z) + fVLL(u, -z) + (3.0 + 2.0 * log(u / (1.0 - u))) * log(z * z)) * this->lcdas->phi(u, mu());
                };

                complex<double> TVLR ;
                if ((u < eps) || (u > 1.0 - eps))
                {
                    TVLR = 0.0;
                }
                else
                {
                    TVLR = (6.0 + 6.0 * 2.0 * log(mu() / mb) - (3.0 + 2.0 * log((1.0 - u) / u)) * log(z * z) - fVLL((1.0 - u), -z) - fVLL(u, -1.0 / z)) * this->lcdas->phi(u, mu());
                };

                // Integration of TSLR gives -6.0, since all u-dependent terms are manifestly symmetric under exchange u <-> ubar = 1 - u
                const double TSLR = -6.0;

                // TSLL vanishes after integration, since the LCDA in two-particle limit is just unity and the hard-scattering kernel is antisymetric
                const double TSLL = 0.0;

                complex<double> TTLL ;
                if ((u < eps) || (u > 1.0 - eps))
                {
                    TTLL = 0.0;
                }
                else
                {
                    TTLL = -48.0 * 2.0 * log(mu() / mb) + 8.0 * (3.0 + ((u - (1.0 - u)) * (1.0 + z)) / (1.0 - z) * log(u / (1.0 - u))) * log(z * z) + fTLL(u, -z) + fTLL(1.0 - u, -1.0 / z);
                };
                return
                    (4.0 * (wc.c1() + wc.c1p()) * (4.0 - TVLL + (2.0 * (4.0 + TSLR) * mu_L) / (mb + mc))) / 9.0 +
                    (2.0 *( wc.c2() + wc.c2p()) * (14.0 + TVLL - (2.0 * (-14.0 + TSLR) * mu_L) / (mb + mc))) / 27.0 +
                    (32.0 * (wc.c3() + wc.c3p()) * (-2.0 * (5.0 + TVLL) + ((-8.0 + TSLR) * mu_L) / (mb + mc))) / 9.0 +
                    (16.0 * (wc.c4() + wc.c4p()) * ((mb + mc) * (19.0 + 2.0 * TVLL) - (-80.0 + TSLR) * mu_L)) / (27.0 * (mb + mc)) +
                    ((wc.c5() + wc.c5p()) * (-4.0 * (6.0 + TVLR) + ((-80.0 - 4.0 * TSLL + TTLL) * mu_L) / (mb + mc))) / 18.0 +
                    ((wc.c6() + wc.c6p()) * (4.0 * (-21.0 + TVLR) + ((-172.0 + 4.0 * TSLL - TTLL) * mu_L) / (mb + mc))) / 108.0 +
                    (wc.c7() + wc.c7p()) * (-32.0 / 9.0 - (2.0 * (-112.0 + 12.0 * TSLL + TTLL) * mu_L) / (9.0 * (mb + mc))) +
                    ((wc.c8() + wc.c8p()) * (-56.0 + ((140.0 + 12.0 * TSLL + TTLL) * mu_L) / (mb + mc))) / 27.0 +
                    (32.0 * (wc.c9() + wc.c9p()) * (40.0 - TVLR + (2.0 * (-48.0 + 4.0 * TSLL + TTLL) * mu_L) / (mb + mc))) / 9.0 +
                    (16.0 * (wc.c10() + wc.c10p()) * ((mb + mc) * (-76.0 + TVLR) - 2.0 * (204.0 + 4.0 * TSLL + TTLL) * mu_L)) / (27.0 * (mb + mc));
            };

            auto a_1_nlo_integrand_re = [this, & a_1_nlo_integrand](const double & u) -> double
            {
                return real(a_1_nlo_integrand(u));
            };
            auto a_1_nlo_integrand_im = [this, & a_1_nlo_integrand](const double & u) -> double
            {
                return imag(a_1_nlo_integrand(u));
            };

            const double a_1_nlo_re = integrate<GSL::QAGS>(a_1_nlo_integrand_re, 0.0, 1.0);
            const double a_1_nlo_im = integrate<GSL::QAGS>(a_1_nlo_integrand_im, 0.0, 1.0);

            const complex<double> a_1_nlo = a_1_nlo_re + a_1_nlo_im * 1.0i;

            // convoluted 3-particle hard-scattering kernels
            const double TVLL_nlp = +(5.0 * lcdas->kappa4(mu()) * m_P * m_P) / (3.0 * (mb * mb - mc * mc));
            const double TTLL_nlp = -(3.0 - lcdas->omega3(mu())) * 2.0 / pow(1.0 - z, 2);

            // calculate contributions from three-particle light-meson states
            const complex<double> a_1_nlp =
                    -(4.0 * (wc.c1() + wc.c1p()) * TVLL_nlp) /3.0 +
                    (2.0 * (wc.c2() + wc.c2p()) * TVLL_nlp) / 9.0 -
                    (64.0 * (wc.c3() + wc.c3p()) * TVLL_nlp) / 3.0 +
                    (32.0 * (wc.c4() + wc.c4p()) * TVLL_nlp) / 9.0 +
                    ((wc.c5() + wc.c5p()) * ((f_3P * m_P * m_P * TTLL_nlp) / (f_P * mb* mb *(mb + mc)) - 4.0 * TVLL_nlp)) / 6.0 +
                    (wc.c6() + wc.c6p()) * (-1.0 / 36.0 * (f_3P * m_P * m_P * TTLL_nlp) / (f_P * mb * mb * (mb + mc)) + TVLL_nlp / 9.0) +
                    (-2.0 * (wc.c7() + wc.c7p()) * f_3P * m_P() * m_P() * TTLL_nlp) / (3.0 * f_P * mb * mb * (mb + mc)) +
                    ((wc.c8() + wc.c8p()) * f_3P * m_P() * m_P() * TTLL_nlp) / (9.0 * f_P * mb * mb * (mb + mc)) +
                    (32.0 * (wc.c9() + wc.c9p()) * ((2.0 * f_3P * m_P * m_P * TTLL_nlp) / (f_P * mb * mb * (mb + mc)) - TVLL_nlp)) / 3.0 +
                    (16.0 * (wc.c10() + wc.c10p()) * ((-2.0 * f_3P * m_P * m_P * TTLL_nlp) / (f_P * mb * mb * (mb + mc)) + TVLL_nlp)) / 9.0;

            // return sum of all contributions
            return switch_lo * a_1_lo + switch_nlo * a_s_mu * a_1_nlo + switch_nlp * a_1_nlp;
        };

        double decay_width() const
        {
            // cf. [BBNS:2000A], eq. (210), p. 80
            const complex<double> amplitude = g_fermi() / sqrt(2.0) * ckm_factor() * f_P() * ff_a_0->evaluate()
                * sqrt(lambda(m_B * m_B, m_Dstar * m_Dstar, m_P * m_P)) * this->a_1();
            // cf. [BBNS:2000A], eq. (216), p. 80
            const double breakup_momentum = sqrt(lambda(m_B * m_B, m_Dstar * m_Dstar, m_P * m_P)) / (2.0 * m_B);

            // cf. [BBNS:2000A], eq. (221), p. 81
            return norm(amplitude) * breakup_momentum / (8.0 * M_PI * m_B * m_B);
        }

        double branching_ratio() const
        {
            return decay_width() * tau_B / hbar;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BqToDstarqPseudoscalar>::options
    {
        Model::option_specification(),
        { "accuracy",     { "LO", "NLO", "NLP", "LO+NLO", "all" }, "all"   },
        { "cp-conjugate", { "true", "false" },                     "false" },
        { "q",            { "s", "d" },                            ""      }
    };

    BqToDstarqPseudoscalar::BqToDstarqPseudoscalar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BqToDstarqPseudoscalar>(new Implementation<BqToDstarqPseudoscalar>(parameters, options, *this))
    {
    }

    BqToDstarqPseudoscalar::~BqToDstarqPseudoscalar()
    {
    }

    double
    BqToDstarqPseudoscalar::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    BqToDstarqPseudoscalar::decay_width() const
    {
        return _imp->decay_width();
    }

    double
    BqToDstarqPseudoscalar::re_a_1() const
    {
        return real(_imp->a_1());
    }

    double
    BqToDstarqPseudoscalar::im_a_1() const
    {
        return imag(_imp->a_1());
    }

    const std::set<ReferenceName>
    BqToDstarqPseudoscalar::references
    {
        "BBNS:2000A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BqToDstarqPseudoscalar::begin_options()
    {
        return Implementation<BqToDstarqPseudoscalar>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BqToDstarqPseudoscalar::end_options()
    {
        return Implementation<BqToDstarqPseudoscalar>::options.cend();
    }
}
