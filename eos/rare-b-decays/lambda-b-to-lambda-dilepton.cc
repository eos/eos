/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014-2025 Danny van Dyk
 * Copyright (c) 2017      Thomas Blake
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

#include <eos/form-factors/baryonic.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-dilepton.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <functional>
#include <map>
#include <string>

namespace eos
{
    using std::abs;
    using std::norm;
    using std::sqrt;

    namespace lambdab_to_lambda_dilepton
    {
        struct Amplitudes
        {
            complex<double> a_perp_0_L, a_perp_0_R;
            complex<double> a_para_0_L, a_para_0_R;
            complex<double> a_perp_1_L, a_perp_1_R;
            complex<double> a_para_1_L, a_para_1_R;
            double alpha;
            double polarisation;
        };

        struct AngularObservables
        {
            std::array<double, 34> _k;

            AngularObservables(const Amplitudes & a)
            {
                // extend to include full basis
                _k[0] = ( norm(a.a_perp_1_R) + norm(a.a_para_1_R) +
                          norm(a.a_perp_1_L) + norm(a.a_para_1_L) +
                          2.0 * norm(a.a_perp_0_R) +
                          2.0 * norm(a.a_para_0_R) +
                          2.0 * norm(a.a_perp_0_L) +
                          2.0 * norm(a.a_para_0_L) ) / 4.0;

                _k[1] = ( norm(a.a_perp_1_R) +
                          norm(a.a_para_1_R) +
                          norm(a.a_perp_1_L) +
                          norm(a.a_para_1_L) ) / 2.0;

                _k[2] = -std::real( a.a_perp_1_R*std::conj(a.a_para_1_R) -
                                    a.a_perp_1_L*std::conj(a.a_para_1_L) );

                _k[3] = std::real( a.a_perp_1_R*std::conj(a.a_para_1_R) +
                                   a.a_perp_1_L*std::conj(a.a_para_1_L) +
                                   2.0 * a.a_perp_0_R*std::conj(a.a_para_0_R) +
                                   2.0 * a.a_perp_0_L*std::conj(a.a_para_0_L) ) * a.alpha / 2.0;

                _k[4] = std::real( a.a_perp_1_R*std::conj(a.a_para_1_R) +
                                   a.a_perp_1_L*std::conj(a.a_para_1_L)) * a.alpha;

                _k[5] = -( norm(a.a_perp_1_R) +
                           norm(a.a_para_1_R) -
                           norm(a.a_perp_1_L) -
                           norm(a.a_para_1_L) ) * a.alpha / 2.0 ;

                _k[6] = -std::real( a.a_para_1_R*std::conj(a.a_perp_0_R) -
                                    a.a_perp_1_R*std::conj(a.a_para_0_R) +
                                    a.a_para_1_L*std::conj(a.a_perp_0_L) -
                                    a.a_perp_1_L*std::conj(a.a_para_0_L) ) * a.alpha / sqrt(2.0);

                _k[7] = -std::real( a.a_para_1_R*std::conj(a.a_para_0_R) -
                                    a.a_perp_1_R*std::conj(a.a_perp_0_R) -
                                    a.a_para_1_L*std::conj(a.a_para_0_L) +
                                    a.a_perp_1_L*std::conj(a.a_perp_0_L) ) * a.alpha /  sqrt(2.0);

                _k[8] = std::imag( a.a_perp_1_R*std::conj(a.a_perp_0_R) -
                                   a.a_para_1_R*std::conj(a.a_para_0_R) +
                                   a.a_perp_1_L*std::conj(a.a_perp_0_L) -
                                   a.a_para_1_L*std::conj(a.a_para_0_L) ) * a.alpha / sqrt(2.0);

                _k[9] = std::imag( a.a_perp_1_R*std::conj(a.a_para_0_R) -
                                   a.a_para_1_R*std::conj(a.a_perp_0_R) -
                                   a.a_perp_1_L*std::conj(a.a_para_0_L) +
                                   a.a_para_1_L*std::conj(a.a_perp_0_L) ) * a.alpha / sqrt(2.0);

                _k[10] = -std::real( a.a_perp_1_R*std::conj(a.a_para_1_R) -
                                    2.0 * a.a_perp_0_R*std::conj(a.a_para_0_R) +
                                    a.a_perp_1_L*std::conj(a.a_para_1_L) -
                                    2.0 * a.a_perp_0_L*std::conj(a.a_para_0_L) ) * a.polarisation / 2.0 ;

                _k[11] = -std::real( a.a_perp_1_R*std::conj(a.a_para_1_R) +
                                    a.a_perp_1_L*std::conj(a.a_para_1_L) ) * a.polarisation;

                _k[12] = ( norm(a.a_perp_1_R) +
                           norm(a.a_para_1_R ) -
                           norm(a.a_perp_1_L) -
                           norm(a.a_para_1_L) ) * a.polarisation / 2.0;

                _k[13] = -( norm(a.a_perp_1_R) +
                            norm(a.a_para_1_R) -
                            2.0 * norm(a.a_perp_0_R) -
                            2.0 * norm(a.a_para_0_R) +
                            norm(a.a_perp_1_L) +
                            norm(a.a_para_1_L) -
                            2.0 * norm(a.a_perp_0_L) -
                            2.0 * norm(a.a_para_0_L) ) * a.alpha * a.polarisation / 4.0;

                _k[14] = -( norm(a.a_perp_1_R) +
                            norm(a.a_para_1_R) +
                            norm(a.a_perp_1_L) +
                            norm(a.a_para_1_L) ) * a.alpha * a.polarisation / 2.0;

                _k[15] = std::real( a.a_perp_1_R*std::conj(a.a_para_1_R) -
                                    a.a_perp_1_L*std::conj(a.a_para_1_L) ) * a.alpha * a.polarisation ;


                _k[16] = std::real( a.a_perp_1_R*std::conj(a.a_perp_0_R) -
                                    a.a_para_1_R*std::conj(a.a_para_0_R) +
                                    a.a_perp_1_L*std::conj(a.a_perp_0_L) -
                                    a.a_para_1_L*std::conj(a.a_para_0_L) ) * a.alpha * a.polarisation / sqrt(2.0);

                _k[17] = std::real( a.a_perp_1_R*std::conj(a.a_para_0_R) -
                                    a.a_para_1_R*std::conj(a.a_perp_0_R) -
                                    a.a_perp_1_L*std::conj(a.a_para_0_L) +
                                    a.a_para_1_L*std::conj(a.a_perp_0_L) ) * a.alpha * a.polarisation / sqrt(2.0);

                _k[18] = -std::imag( a.a_para_1_R*std::conj(a.a_perp_0_R) -
                                     a.a_perp_1_R*std::conj(a.a_para_0_R) +
                                     a.a_para_1_L*std::conj(a.a_perp_0_L) -
                                     a.a_perp_1_L*std::conj(a.a_para_0_L) ) * a.alpha * a.polarisation / sqrt(2.0);

                _k[19] = -std::imag( a.a_para_1_R*std::conj(a.a_para_0_R) -
                                    a.a_perp_1_R*std::conj(a.a_perp_0_R) -
                                    a.a_para_1_L*std::conj(a.a_para_0_L) +
                                    a.a_perp_1_L*std::conj(a.a_perp_0_L) ) * a.alpha * a.polarisation / sqrt(2.0);

                _k[20] =  std::imag( a.a_para_1_R*std::conj(a.a_para_0_R) +
                                     a.a_perp_1_R*std::conj(a.a_perp_0_R) +
                                     a.a_para_1_L*std::conj(a.a_para_0_L) +
                                     a.a_perp_1_L*std::conj(a.a_perp_0_L) ) * a.polarisation / sqrt(2.0);

                _k[21] = -std::imag( a.a_perp_1_R*std::conj(a.a_para_0_R) +
                                     a.a_para_1_R*std::conj(a.a_perp_0_R) -
                                     a.a_perp_1_L*std::conj(a.a_para_0_L) -
                                     a.a_para_1_L*std::conj(a.a_perp_0_L) ) * a.polarisation / sqrt(2.0);

                _k[22] = -std::real( a.a_perp_1_R*std::conj(a.a_para_0_R) +
                                     a.a_para_1_R*std::conj(a.a_perp_0_R) +
                                     a.a_perp_1_L*std::conj(a.a_para_0_L) +
                                     a.a_para_1_L*std::conj(a.a_perp_0_L) ) * a.polarisation / sqrt(2.0);

                _k[23] = std::real( a.a_perp_1_R*std::conj(a.a_perp_0_R) +
                                    a.a_para_1_R*std::conj(a.a_para_0_R) -
                                    a.a_perp_1_L*std::conj(a.a_perp_0_L) -
                                    a.a_para_1_L*std::conj(a.a_para_0_L) ) * a.polarisation / sqrt(2.0);

                _k[24] = std::imag( a.a_perp_1_R*std::conj(a.a_para_0_R) +
                                    a.a_para_1_R*std::conj(a.a_perp_0_R) +
                                    a.a_perp_1_L*std::conj(a.a_para_0_L) +
                                    a.a_para_1_L*std::conj(a.a_perp_0_L) ) * a.alpha * a.polarisation / sqrt(2.0);

                _k[25] = -std::imag( a.a_perp_1_R*std::conj(a.a_perp_0_R) +
                                     a.a_para_1_R*std::conj(a.a_para_0_R) -
                                     a.a_perp_1_L*std::conj(a.a_perp_0_L) -
                                     a.a_para_1_L*std::conj(a.a_para_0_L) ) * a.alpha * a.polarisation / sqrt(2.0);

                _k[26] = -std::real( a.a_perp_1_R*std::conj(a.a_perp_0_R) +
                                     a.a_para_1_R*std::conj(a.a_para_0_R) +
                                     a.a_perp_1_L*std::conj(a.a_perp_0_L) +
                                     a.a_para_1_L*std::conj(a.a_para_0_L) ) * a.alpha * a.polarisation /  sqrt(2.0);

                _k[27] = std::real( a.a_perp_1_R*std::conj(a.a_para_0_R) +
                                    a.a_para_1_R*std::conj(a.a_perp_0_R) -
                                    a.a_perp_1_L*std::conj(a.a_para_0_L) -
                                    a.a_para_1_L*std::conj(a.a_perp_0_L) ) * a.alpha * a.polarisation / sqrt(2.0);

                _k[28] = 0.0;

                _k[29] = std::imag( a.a_perp_0_R*std::conj(a.a_para_0_R) +
                                    a.a_perp_0_L*std::conj(a.a_para_0_L) ) * a.alpha * a.polarisation;

                _k[30] = 0.0;

                _k[31] = ( norm(a.a_perp_0_R) -
                           norm(a.a_para_0_R) +
                           norm(a.a_perp_0_L) -
                           norm(a.a_para_0_L) ) * a.alpha * a.polarisation / 2.0;

                _k[32] = ( norm(a.a_perp_1_R) -
                           norm(a.a_para_1_R) +
                           norm(a.a_perp_1_L) -
                           norm(a.a_para_1_L) ) * a.alpha * a.polarisation / 4.0;

                _k[33] = std::imag( a.a_perp_1_R*std::conj(a.a_para_1_R) +
                                    a.a_perp_1_L*std::conj(a.a_para_1_L) ) * a.alpha * a.polarisation / 2.0;
            }

            AngularObservables(const std::array<double, 34> & k) :
                _k(k)
            {
            }

            inline double k1ss() const { return _k[0]; }
            inline double k1cc() const { return _k[1]; }
            inline double k1c()  const { return _k[2]; }
            inline double k2ss() const { return _k[3]; }
            inline double k2cc() const { return _k[4]; }
            inline double k2c()  const { return _k[5]; }
            inline double k3sc() const { return _k[8]; }
            inline double k3s()  const { return _k[9]; }
            inline double k4sc() const { return _k[6]; }
            inline double k4s()  const { return _k[7]; }


            inline double k1()  const { return _k[0]; }
            inline double k2()  const { return _k[1]; }
            inline double k3()  const { return _k[2]; }
            inline double k4()  const { return _k[3]; }
            inline double k5()  const { return _k[4]; }
            inline double k6()  const { return _k[5]; }
            inline double k7()  const { return _k[6]; }
            inline double k8()  const { return _k[7]; }
            inline double k9()  const { return _k[8]; }
            inline double k10() const { return _k[9]; }
            inline double k11() const { return _k[10]; }
            inline double k12() const { return _k[11]; }
            inline double k13() const { return _k[12]; }
            inline double k14() const { return _k[13]; }
            inline double k15() const { return _k[14]; }
            inline double k16() const { return _k[15]; }
            inline double k17() const { return _k[16]; }
            inline double k18() const { return _k[17]; }
            inline double k19() const { return _k[18]; }
            inline double k20() const { return _k[19]; }
            inline double k21() const { return _k[20]; }
            inline double k22() const { return _k[21]; }
            inline double k23() const { return _k[22]; }
            inline double k24() const { return _k[23]; }
            inline double k25() const { return _k[24]; }
            inline double k26() const { return _k[25]; }
            inline double k27() const { return _k[26]; }
            inline double k28() const { return _k[27]; }
            inline double k29() const { return _k[28]; }
            inline double k30() const { return _k[29]; }
            inline double k31() const { return _k[30]; }
            inline double k32() const { return _k[31]; }
            inline double k33() const { return _k[32]; }
            inline double k34() const { return _k[33]; }

            inline double decay_width() const
            {
                return 2.0 * k1ss() + k1cc();
            }

            inline double a_fb_leptonic() const
            {
                return 3.0 / 2.0 * k1c() / decay_width();
            }

            inline double a_fb_hadronic() const
            {
                return 1.0 / 2.0 * (2.0 * k2ss() + k2cc()) / decay_width();
            }

            inline double a_fb_combined() const
            {
                return 3.0 / 4.0 * k2c() / decay_width();
            }

            inline double f_zero() const
            {
                return (2.0 * k1ss() - k1cc()) / decay_width();
            }
        };
    }

    /* Large Recoil */

    template <> struct Implementation<LambdaBToLambdaDilepton<LargeRecoil>>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;
        UsedParameter tau_Lambda_b;

        UsedParameter g_fermi;

        UsedParameter m_Lambda_b;
        UsedParameter m_Lambda;
        UsedParameter alpha;
        UsedParameter polarisation;

        UsedParameter alpha_e;

        LeptonFlavorOption opt_l;
        UsedParameter mu;

        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            hbar(p["QM::hbar"], u),
            tau_Lambda_b(p["life_time::Lambda_b"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_Lambda_b(p["mass::Lambda_b"], u),
            m_Lambda(p["mass::Lambda"], u),
            alpha(p["Lambda::alpha"], u),
            polarisation(p["Lambda_b::polarisation@" + o.get("production-polarisation"_ok,"unpolarised") ], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            opt_l(o, options, "l"_ok),
            mu(p["sb" + opt_l.str() + opt_l.str() + "::mu"], u),
            form_factors(FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda::" + o.get("form-factors"_ok, "BFvD2014"), p, o))
        {
            Context ctx("When constructing L_b->Lll observables");

            u.uses(*form_factors);
            u.uses(*model);
        }

        double norm(const double & s) const
        {
            return g_fermi() * alpha_e() * abs(model->ckm_tb() * conj(model->ckm_ts()))
                * sqrt(s / 3.0 / 2048 / power_of<5>(M_PI) / power_of<3>(m_Lambda_b())
                * sqrt(lambda(m_Lambda_b * m_Lambda_b, m_Lambda * m_Lambda, s))); // cf. [BFvD:2014], Eq. (?), p. ??
        }

        double kappa() const
        {
            return (1.0 - 2.0 * model->alpha_s(mu) / (3.0 * M_PI) * std::log(mu / model->m_b_msbar(mu)));
        }

        lambdab_to_lambda_dilepton::Amplitudes amplitudes(const double & s)
        {
            lambdab_to_lambda_dilepton::Amplitudes result;

            double alpha_s = model->alpha_s(mu());
            double m_b_MSbar = model->m_b_msbar(mu), m_b_PS = model->m_b_ps(2.0), m_b_PS2 = m_b_PS * m_b_PS;
            double m_c_pole = model->m_c_pole();

            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), opt_l.value());

            complex<double> lambda_hat_u = model->ckm_ub() * conj(model->ckm_us()) / std::abs(model->ckm_tb() * conj(model->ckm_ts()));
            double sqrtsminus = sqrt(power_of<2>(m_Lambda_b - m_Lambda) - s), sqrtsplus = sqrt(power_of<2>(m_Lambda_b + m_Lambda) - s), sqrts = sqrt(s);
            double N = norm(s);

            /* Y(s) for the up and the top sector */
            // cf. [BFS:2001], Eq. (10), p. 4
            complex<double> Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2() + 6.0 * wc.c3() + 60.0 * wc.c5();
            complex<double> Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
            complex<double> Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64 / 3.0 * wc.c6());
            complex<double> Y_top_ = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());
            // Use b pole mass according to [BFS:2001], Sec. 3.1, paragraph Quark Masses,
            // then replace b pole mass by the PS mass.
            complex<double> Y_top = Y_top_c * CharmLoops::h(mu, s, m_c_pole)
                 + Y_top_b * CharmLoops::h(mu, s, m_b_PS)
                 + Y_top_0 * CharmLoops::h(mu, s)
                 + Y_top_;
            // cf. [BFS:2004], Eq. (43), p. 24
            complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (CharmLoops::h(mu, s, m_c_pole) - CharmLoops::h(mu, s));

            // calculate effective wilson coefficients
            // cf. [BFS:2001], below Eq. (9), p. 4
            complex<double> c7eff = wc.c7() - 1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6();
            // cf. [BFS:2001], below Eq. (26), p. 8
            complex<double> c8eff = wc.c8() + wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

            // two loop virtual corrections, cf. [AAGW:2001]
            // charm quarks
            complex<double> F27c = CharmLoops::F27_massive(mu(), s, m_b_PS, m_c_pole);
            complex<double> F17c = -F27c / 6.0;
            complex<double> F19c = CharmLoops::F19_massive(mu(), s, m_b_PS, m_c_pole);
            complex<double> F29c = CharmLoops::F29_massive(mu(), s, m_b_PS, m_c_pole);
            // up quarks
            complex<double> F27u = CharmLoops::F27_massless(mu(), s, m_b_PS);
            complex<double> F17u = -F27u / 6.0;
            complex<double> F19u = CharmLoops::F19_massless(mu(), s, m_b_PS);
            complex<double> F29u = CharmLoops::F29_massless(mu(), s, m_b_PS);
            // gluon
            complex<double> F87  = CharmLoops::F87_massless(mu(), s, m_b_PS);
            complex<double> F89  = CharmLoops::F89_massless(s, m_b_PS);

            // integredients for form factor relations
            // cf. [FY:2011]
            double L = -1.0 * (m_b_PS2 - s) / s * std::log(1.0 - s / m_b_PS2);

            // ratio of tensor to vector form factors
            // cf. [BFvD:2014], eqs. (??)-(??)
            double R1p = 1.0 + alpha_s / (3.0 * M_PI) * (2.0 * std::log(m_b_PS / mu()) - 2.0 - L);
            double R1m = 1.0 + alpha_s / (3.0 * M_PI) * (2.0 * std::log(m_b_PS / mu()) - 2.0 - L);
            double R0p = 1.0 + alpha_s / (3.0 * M_PI) * (2.0 * std::log(m_b_PS / mu()) - 2.0 + 2.0 * L);
            double R0m = 1.0 + alpha_s / (3.0 * M_PI) * (2.0 * std::log(m_b_PS / mu()) - 2.0 + 2.0 * L);

            // functions tau
            // cf. [BFvD:2014], eqs. (??)-(??)
            complex<double> tau_1p = (m_Lambda_b + m_Lambda) / m_Lambda_b * (
                    c7eff + wc.c7prime()
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F17c + wc.c2() * F27c + c8eff * F87)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F17c - F17u) + wc.c2() * (F27c - F27u) + c8eff * F87) * lambda_hat_u
                ) * R1p
                + s / (2.0 * m_b_MSbar * m_Lambda_b) * (
                    Y_top + lambda_hat_u * Y_up
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F19c + wc.c2() * F29c + wc.c8() * F89)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F19c - F19u) + wc.c2() * (F29c - F29u) + wc.c8() * F89) * lambda_hat_u
                );
            complex<double> tau_1m = (m_Lambda_b - m_Lambda) / m_Lambda_b * (
                    c7eff - wc.c7prime()
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F17c + wc.c2() * F27c)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F17c - F17u) + wc.c2() * (F27c - F27u) + c8eff * F87) * lambda_hat_u
                ) * R1m
                + s / (2.0 * m_b_MSbar * m_Lambda_b) * (
                    Y_top + lambda_hat_u * Y_up
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F19c + wc.c2() * F29c + c8eff * F87)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F19c - F19u) + wc.c2() * (F29c - F29u) + wc.c8() * F89) * lambda_hat_u
                );
            complex<double> tau_0p = m_Lambda_b / (m_Lambda_b + m_Lambda) * (
                    c7eff + wc.c7prime()
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F17c + wc.c2() * F27c)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F17c - F17u) + wc.c2() * (F27c - F27u) + c8eff * F87) * lambda_hat_u
                ) * R0p
                + m_Lambda_b / (2.0 * m_b_MSbar) * (
                    Y_top + lambda_hat_u * Y_up
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F19c + wc.c2() * F29c)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F19c - F19u) + wc.c2() * (F29c - F29u) + wc.c8() * F89) * lambda_hat_u
                );
            complex<double> tau_0m = m_Lambda_b / (m_Lambda_b - m_Lambda) * (
                    c7eff - wc.c7prime()
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F17c + wc.c2() * F27c)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F17c - F17u) + wc.c2() * (F27c - F27u) + c8eff * F87) * lambda_hat_u
                ) * R0m
                + m_Lambda_b / (2.0 * m_b_MSbar) * (
                    Y_top + lambda_hat_u * Y_up
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * F19c + wc.c2() * F29c)
                    - alpha_s / (4.0 * M_PI) * (wc.c1() * (F19c - F19u) + wc.c2() * (F29c - F29u) + wc.c8() * F89) * lambda_hat_u
                );

            // cf. [BFvD:2014], eqs. (??)-(??)
            result.a_perp_1_R = -2.0 *       N * (wc.c9() + wc.c9prime() + (wc.c10() + wc.c10prime()) + 2.0 * m_b_MSbar * m_Lambda_b / s * tau_1p) * form_factors->f_perp_v(s) * sqrtsminus;
            result.a_perp_1_L = -2.0 *       N * (wc.c9() + wc.c9prime() - (wc.c10() + wc.c10prime()) + 2.0 * m_b_MSbar * m_Lambda_b / s * tau_1p) * form_factors->f_perp_v(s) * sqrtsminus;

            result.a_para_1_R = +2.0 *       N * (wc.c9() - wc.c9prime() + (wc.c10() - wc.c10prime()) + 2.0 * m_b_MSbar * m_Lambda_b / s * tau_1m) * form_factors->f_perp_a(s) * sqrtsplus;
            result.a_para_1_L = +2.0 *       N * (wc.c9() - wc.c9prime() - (wc.c10() - wc.c10prime()) + 2.0 * m_b_MSbar * m_Lambda_b / s * tau_1m) * form_factors->f_perp_a(s) * sqrtsplus;

            result.a_perp_0_R = +sqrt(2.0) * N * (wc.c9() + wc.c9prime() + (wc.c10() + wc.c10prime()) + 2.0 * m_b_MSbar / m_Lambda_b * tau_0p) * form_factors->f_long_v(s) * (m_Lambda_b + m_Lambda)  / sqrts * sqrtsminus;
            result.a_perp_0_L = +sqrt(2.0) * N * (wc.c9() + wc.c9prime() - (wc.c10() + wc.c10prime()) + 2.0 * m_b_MSbar / m_Lambda_b * tau_0p) * form_factors->f_long_v(s) * (m_Lambda_b + m_Lambda)  / sqrts * sqrtsminus;

            result.a_para_0_R = -sqrt(2.0) * N * (wc.c9() - wc.c9prime() + (wc.c10() - wc.c10prime()) + 2.0 * m_b_MSbar / m_Lambda_b * tau_0m) * form_factors->f_long_a(s) * (m_Lambda_b - m_Lambda)  / sqrts * sqrtsplus;
            result.a_para_0_L = -sqrt(2.0) * N * (wc.c9() - wc.c9prime() - (wc.c10() - wc.c10prime()) + 2.0 * m_b_MSbar / m_Lambda_b * tau_0m) * form_factors->f_long_a(s) * (m_Lambda_b - m_Lambda)  / sqrts * sqrtsplus;

            result.alpha = this->alpha();
            result.polarisation = this->polarisation();

            return result;
        }

        std::array<double, 34> _differential_angular_observables(const double & s)
        {
            return lambdab_to_lambda_dilepton::AngularObservables(this->amplitudes(s))._k;
        }

        std::array<double, 34> _integrated_angular_observables(const double & s_min, const double & s_max)
        {
            std::function<std::array<double, 34> (const double &)> integrand(std::bind(&Implementation::_differential_angular_observables, this, std::placeholders::_1));

            return integrate1D(integrand, 64, s_min, s_max);
        }

        inline lambdab_to_lambda_dilepton::AngularObservables differential_angular_observables(const double & s)
        {
            return lambdab_to_lambda_dilepton::AngularObservables{ _differential_angular_observables(s) };
        }

        inline lambdab_to_lambda_dilepton::AngularObservables integrated_angular_observables(const double & s_min, const double & s_max)
        {
            return lambdab_to_lambda_dilepton::AngularObservables{ _integrated_angular_observables(s_min, s_max) };
        }
    };

    LambdaBToLambdaDilepton<LargeRecoil>::LambdaBToLambdaDilepton(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<LambdaBToLambdaDilepton<LargeRecoil>>(new Implementation<LambdaBToLambdaDilepton<LargeRecoil>>(p, o, *this))
    {
    }

    LambdaBToLambdaDilepton<LargeRecoil>::~LambdaBToLambdaDilepton()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<LambdaBToLambdaDilepton<LargeRecoil>>::options
    {
        { "l"_ok, {"e", "mu", "tau"}, "mu" }
    };

    /* q^2-differential observables */
    double
    LambdaBToLambdaDilepton<LargeRecoil>::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_angular_observables(s).decay_width() * _imp->tau_Lambda_b / _imp->hbar;
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_leptonic();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_hadronic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_hadronic();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_combined(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_combined();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::differential_fzero(const double & s) const
    {
        return _imp->differential_angular_observables(s).f_zero();
    }

    /* q^2-integrated observables */
    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).decay_width() * _imp->tau_Lambda_b / _imp->hbar;
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_hadronic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_hadronic();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_combined(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_combined();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_fzero(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_zero();
    }

    /* Polarised angular observables */
    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m1(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1() / o.decay_width();
    }


    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m2(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m3(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k3() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m4(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k4() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m5(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k5() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m6(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k6() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m7(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k7() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m8(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k8() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m9(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k9() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m10(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k10() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m11(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k11() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m12(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k12() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m13(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k13() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m14(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k14() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m15(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k15() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m16(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k16() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m17(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k17() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m18(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k18() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m19(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k19() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m20(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k20() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m21(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k21() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m22(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k22() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m23(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k23() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m24(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k24() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m25(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k25() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m26(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k26() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m27(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k27() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m28(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k28() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m29(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k29() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m30(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k30() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m31(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k31() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m32(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k32() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m33(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k33() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LargeRecoil>::integrated_m34(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k34() / o.decay_width();
    }

    const std::set<ReferenceName>
    LambdaBToLambdaDilepton<LargeRecoil>::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambdaDilepton<LargeRecoil>::begin_options()
    {
        return Implementation<LambdaBToLambdaDilepton<LargeRecoil>>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambdaDilepton<LargeRecoil>::end_options()
    {
        return Implementation<LambdaBToLambdaDilepton<LargeRecoil>>::options.cend();
    }


    /* Low Recoil */

    template <> struct Implementation<LambdaBToLambdaDilepton<LowRecoil>>
    {
        std::shared_ptr<Model> model;

        LeptonFlavorOption opt_l;

        UsedParameter hbar;
        UsedParameter tau_Lambda_b;

        UsedParameter g_fermi;

        UsedParameter m_Lambda_b;
        UsedParameter m_Lambda;
        UsedParameter alpha;
        UsedParameter polarisation;

        UsedParameter alpha_e;
        UsedParameter mu;

        UsedParameter r_perp_0, r_perp_1;
        UsedParameter r_para_0, r_para_1;

        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            opt_l(o, options, "l"_ok),
            hbar(p["QM::hbar"], u),
            tau_Lambda_b(p["life_time::Lambda_b"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_Lambda_b(p["mass::Lambda_b"], u),
            m_Lambda(p["mass::Lambda"], u),
            alpha(p["Lambda::alpha"], u),
            polarisation(p["Lambda_b::polarisation@" + o.get("production-polarisation"_ok,"unpolarised") ], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            mu(p["sb" + opt_l.str() + opt_l.str() + "::mu"], u),
            r_perp_0(p["Lambda_b->Lambdall::r_perp_0@MvD2016"], u),
            r_perp_1(p["Lambda_b->Lambdall::r_perp_1@MvD2016"], u),
            r_para_0(p["Lambda_b->Lambdall::r_para_0@MvD2016"], u),
            r_para_1(p["Lambda_b->Lambdall::r_para_1@MvD2016"], u),
            form_factors(FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda::" + o.get("form-factors"_ok, "DM2016"), p, o))
        {
            u.uses(*form_factors);
            u.uses(*model);
        }

        double norm(const double & s) const
        {
            return g_fermi() * alpha_e() * abs(model->ckm_tb() * conj(model->ckm_ts()))
                * sqrt(s / 3.0 / 2048 / power_of<5>(M_PI) / power_of<3>(m_Lambda_b())
                * sqrt(lambda(m_Lambda_b * m_Lambda_b, m_Lambda * m_Lambda, s))); // cf. [BFvD:2014], Eq. (3.18), p. 6
        }

        double kappa() const
        {
            return (1.0 - 2.0 * model->alpha_s(mu) / (3.0 * M_PI) * std::log(mu / model->m_b_msbar(mu)));
        }

        lambdab_to_lambda_dilepton::Amplitudes amplitudes(const double & s)
        {
            lambdab_to_lambda_dilepton::Amplitudes result;

            double alpha_s = model->alpha_s(mu()), m_b = model->m_b_ps(2.0), m_c = model->m_c_msbar(mu());
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), opt_l.value());
            complex<double> lambda_hat_u = model->ckm_ub() * conj(model->ckm_us()) / abs(model->ckm_tb() * conj(model->ckm_ts()));
            double sqrtsminus = sqrt(power_of<2>(m_Lambda_b - m_Lambda) - s), sqrtsplus = sqrt(power_of<2>(m_Lambda_b + m_Lambda) - s), sqrts = sqrt(s);
            double N = norm(s), kappa = this->kappa();

            // calculate effective wilson coefficients
            complex<double> c7eff = ShortDistanceLowRecoil::c7eff(s, mu(), alpha_s, m_b, true, wc);
            complex<double> c9eff = ShortDistanceLowRecoil::c9eff(s, mu(), alpha_s, m_b, m_c, true, false, lambda_hat_u, wc);

            // cf. [BFvD:2014], eq.s (??), p. ??
            double zeta_perp_V = (m_Lambda_b + m_Lambda) / m_Lambda_b * form_factors->f_perp_t(s)  / form_factors->f_perp_v(s);
            double zeta_perp_A = (m_Lambda_b - m_Lambda) / m_Lambda_b * form_factors->f_perp_t5(s) / form_factors->f_perp_a(s);
            double zeta_long_V = s / ((m_Lambda_b + m_Lambda) * m_Lambda_b) * form_factors->f_long_t(s)  / form_factors->f_long_v(s);
            double zeta_long_A = s / ((m_Lambda_b - m_Lambda) * m_Lambda_b) * form_factors->f_long_t5(s) / form_factors->f_long_a(s);

            // parametrize subleading power corrections, cf. [MvD:2016], eq. (B1), p. ??
            complex<double> x_perp_0 = (4.0 / 3.0 * wc.c1() + wc.c2()) * r_perp_0();
            complex<double> x_perp_1 = (4.0 / 3.0 * wc.c1() + wc.c2()) * r_perp_1();
            complex<double> x_para_0 = (4.0 / 3.0 * wc.c1() + wc.c2()) * r_para_0();
            complex<double> x_para_1 = (4.0 / 3.0 * wc.c1() + wc.c2()) * r_para_1();

            // cf. [BFvD:2014], eqs. (4.9)-(4.10), p. 11
            result.a_perp_1_R = -2.0 *       N * (c9eff + wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff + wc.c7prime()) * zeta_perp_V + (wc.c10() + wc.c10prime()) + x_perp_1) * form_factors->f_perp_v(s) * sqrtsminus;
            result.a_perp_1_L = -2.0 *       N * (c9eff + wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff + wc.c7prime()) * zeta_perp_V - (wc.c10() + wc.c10prime()) + x_perp_1) * form_factors->f_perp_v(s) * sqrtsminus;

            result.a_para_1_R = +2.0 *       N * (c9eff - wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff - wc.c7prime()) * zeta_perp_A + (wc.c10() - wc.c10prime()) + x_para_1) * form_factors->f_perp_a(s) * sqrtsplus;
            result.a_para_1_L = +2.0 *       N * (c9eff - wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff - wc.c7prime()) * zeta_perp_A - (wc.c10() - wc.c10prime()) + x_para_1) * form_factors->f_perp_a(s) * sqrtsplus;

            result.a_perp_0_R = +sqrt(2.0) * N * (c9eff + wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff + wc.c7prime()) * zeta_long_V + (wc.c10() + wc.c10prime()) + x_perp_0) * form_factors->f_long_v(s)
                * (m_Lambda_b + m_Lambda)  / sqrts * sqrtsminus;
            result.a_perp_0_L = +sqrt(2.0) * N * (c9eff + wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff + wc.c7prime()) * zeta_long_V - (wc.c10() + wc.c10prime()) + x_perp_0) * form_factors->f_long_v(s)
                * (m_Lambda_b + m_Lambda)  / sqrts * sqrtsminus;

            result.a_para_0_R = -sqrt(2.0) * N * (c9eff - wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff - wc.c7prime()) * zeta_long_A + (wc.c10() - wc.c10prime()) + x_para_0) * form_factors->f_long_a(s)
                * (m_Lambda_b - m_Lambda)  / sqrts * sqrtsplus;
            result.a_para_0_L = -sqrt(2.0) * N * (c9eff - wc.c9prime() + (2.0 * kappa * m_b * m_Lambda_b / s) * (c7eff - wc.c7prime()) * zeta_long_A - (wc.c10() - wc.c10prime()) + x_para_0) * form_factors->f_long_a(s)
                * (m_Lambda_b - m_Lambda)  / sqrts * sqrtsplus;

            result.alpha = this->alpha();
            result.polarisation = this->polarisation();

            return result;
        }

        std::array<double, 34> _differential_angular_observables(const double & s)
        {
            return lambdab_to_lambda_dilepton::AngularObservables(this->amplitudes(s))._k;
        }

        std::array<double, 34> _integrated_angular_observables(const double & s_min, const double & s_max)
        {
            std::function<std::array<double, 34> (const double &)> integrand(std::bind(&Implementation::_differential_angular_observables, this, std::placeholders::_1));

            return integrate1D(integrand, 64, s_min, s_max);
        }

        inline lambdab_to_lambda_dilepton::AngularObservables differential_angular_observables(const double & s)
        {
            return lambdab_to_lambda_dilepton::AngularObservables{ _differential_angular_observables(s) };
        }

        inline lambdab_to_lambda_dilepton::AngularObservables integrated_angular_observables(const double & s_min, const double & s_max)
        {
            return lambdab_to_lambda_dilepton::AngularObservables{ _integrated_angular_observables(s_min, s_max) };
        }
    };

    LambdaBToLambdaDilepton<LowRecoil>::LambdaBToLambdaDilepton(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<LambdaBToLambdaDilepton<LowRecoil>>(new Implementation<LambdaBToLambdaDilepton<LowRecoil>>(p, o, *this))
    {
    }

    LambdaBToLambdaDilepton<LowRecoil>::~LambdaBToLambdaDilepton()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<LambdaBToLambdaDilepton<LowRecoil>>::options
    {
        Model::option_specification(),
        { "l"_ok, {"e", "mu", "tau"}, "mu" }
    };

    /* q^2-differential observables */
    double
    LambdaBToLambdaDilepton<LowRecoil>::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_angular_observables(s).decay_width() * _imp->tau_Lambda_b / _imp->hbar;
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_leptonic();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_hadronic(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_hadronic();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_combined(const double & s) const
    {
        return _imp->differential_angular_observables(s).a_fb_combined();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::differential_fzero(const double & s) const
    {
        return _imp->differential_angular_observables(s).f_zero();
    }

    /* q^2-integrated observables */
    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).decay_width() * _imp->tau_Lambda_b / _imp->hbar;
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_hadronic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_hadronic();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_combined(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_combined();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_fzero(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_zero();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_k1ss(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1ss() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_k1cc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1cc() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_k1c(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1c() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_k2ss(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2ss() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_k2cc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2cc() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_k2c(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2c() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_k3sc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k3sc() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_k3s(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k3s() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_k4sc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k4sc() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_k4s(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k4s() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m1(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1() / o.decay_width();
    }


    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m2(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m3(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k3() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m4(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k4() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m5(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k5() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m6(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k6() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m7(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k7() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m8(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k8() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m9(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k9() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m10(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k10() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m11(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k11() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m12(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k12() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m13(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k13() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m14(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k14() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m15(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k15() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m16(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k16() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m17(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k17() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m18(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k18() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m19(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k19() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m20(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k20() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m21(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k21() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m22(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k22() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m23(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k23() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m24(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k24() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m25(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k25() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m26(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k26() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m27(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k27() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m28(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k28() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m29(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k29() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m30(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k30() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m31(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k31() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m32(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k32() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m33(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k33() / o.decay_width();
    }

    double
    LambdaBToLambdaDilepton<LowRecoil>::integrated_m34(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k34() / o.decay_width();
    }

    const std::set<ReferenceName>
    LambdaBToLambdaDilepton<LowRecoil>::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambdaDilepton<LowRecoil>::begin_options()
    {
        return Implementation<LambdaBToLambdaDilepton<LowRecoil>>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    LambdaBToLambdaDilepton<LowRecoil>::end_options()
    {
        return Implementation<LambdaBToLambdaDilepton<LowRecoil>>::options.cend();
    }
}
