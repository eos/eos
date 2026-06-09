/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
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
#include <eos/maths/polylog.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/rare-b-decays/inclusive-b-to-s-gamma.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

#include <iostream>

namespace eos
{
    /* Minimal */

    template <>
    struct Implementation<BToXsGamma<Minimal>>
    {
        std::shared_ptr<Model> model;

        UsedParameter m_b_MSbar;

        UsedParameter alpha_e;

        UsedParameter br_bcsl;

        UsedParameter uncertainty;

        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            m_b_MSbar(p["mass::b(MSbar)"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            br_bcsl(p["exp::BR(B->X_clnu)"], u),
            uncertainty(p["B->X_sgamma::uncertainty"], u),
            mu(p["sb::mu"], u)
        {
            Context ctx("When constructing B->X_sgamma observables");

            if ("SM" != o.get("model"_ok, "SM"))
            {
                Log::instance()->message("B->X_sgamma.model", ll_error)
                    << "B->X_sgamma is not yet capable to handle models beyond SM, e.g. for helicity flipped operators; use it carefully";
            }

            u.uses(*model);
        }

        double m_c_pole() const
        {
            return 1.6;
        }

        double branching_ratio() const
        {
            static const double sm = 3.15e-4;
            static const double sm_delta = 0.23e-4;
            static const double c7sm = -0.3;

            double m_c_hat = m_c_pole() / model->m_b_pole();
            double z = power_of<2>(m_c_hat);
            double z2 = power_of<2>(z), z3 = z * z2, z4 = z3 * z, lnz = log(z);

            // cf. [BMU:1999A], Eq. (46), p. 16
            double g = 1.0 - 8.0 * z + 8.0 * z3 - z4 - 12.0 * z2 * lnz;
            double kappa = 1.0 - 2.0/3.0 * model->alpha_s(model->m_b_pole()) / M_PI * (1.5 + (M_PI * M_PI - 31.0 / 4.0) * power_of<2>(1.0 - m_c_hat));

            double ckm = norm(model->ckm_tb() * conj(model->ckm_ts()) / model->ckm_cb());
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), LeptonFlavor::muon /* fake lepton flavor */);
            complex<double> c7np = wc.c7() - c7sm;

            double result = (sm + sm_delta * uncertainty)
                + 6.0 * alpha_e / M_PI * br_bcsl * ckm / g / kappa * (norm(c7np) + 2.0 * real(c7np * c7sm));

            // Make sure the approximate BR is positive definite.
            return std::max(result, 0.0);
        }
    };

    BToXsGamma<Minimal>::BToXsGamma(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToXsGamma<Minimal>>(new Implementation<BToXsGamma<Minimal>>(parameters, options, *this))
    {
    }

    BToXsGamma<Minimal>::~BToXsGamma()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<BToXsGamma<Minimal>>::options
    {
    };

    double
    BToXsGamma<Minimal>::integrated_branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    const std::set<ReferenceName>
    BToXsGamma<Minimal>::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BToXsGamma<Minimal>::begin_options()
    {
        return Implementation<BToXsGamma<Minimal>>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToXsGamma<Minimal>::end_options()
    {
        return Implementation<BToXsGamma<Minimal>>::options.cend();
    }

    /* NLO
     *
     * Implementation according to [BCMU:2002A],
     * which bases on [CMM:1996A].
     */

    template <>
    struct Implementation<BToXsGamma<NLO>>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter mu;

        UsedParameter m_B;

        UsedParameter mu2_g;

        UsedParameter mu2_pi;

        UsedParameter m_b_MSbar;

        UsedParameter alpha_e;

        UsedParameter gfermi;

        UsedParameter tau;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            hbar(p["QM::hbar"], u),
            mu(p["sb::mu"], u),
            m_B(p["mass::B_d"], u),
            mu2_g(p["B->B::mu_G^2@1GeV"], u),
            mu2_pi(p["B->B::mu_pi^2@1GeV"], u),
            m_b_MSbar(p["mass::b(MSbar)"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            gfermi(p["WET::G_Fermi"], u),
            tau(p["life_time::B" + (destringify<bool>(o.get("admixture"_ok, "true")) ? ("@Y(4S)") : ("_" + o.get("q"_ok, "d")))], u)
        {
            u.uses(*model);
        }

        static const std::array<double, 225> f22int1coeffs, f22int2coeffs;
        static const std::array<double, 225> f27int1coeffs, f27int2coeffs;

        struct Fij
        {
            complex<double> f11, f12, f17, f18;
            complex<double> f22, f27, f28;
            complex<double> f77, f78;
            complex<double> f88;
        };

        /*
         * Compute functions f_ij, cf. [CMM:1996A], p.9, Eq. (39) and below.
         * We use the fit to a polynomial for the results of integrals
         * that arise in f_22 and f_27.
         */
        Fij f_ij(const double & z, const double & lnz, const double & delta) const
        {
            // Constants
            static const double pi = M_PI, pi2 = pi * pi;

            double ztilde = z - 0.08, deltatilde = delta - 0.1, onemdelta = 1.0 - delta;
            double lndelta = std::log(delta), ln1mdelta = std::log(onemdelta);
            double dl1mdelta = std::real(dilog(onemdelta));
            double delta2 = delta * delta, delta3 = delta2 * delta;
            std::array<double, 5> zz{{ 1.0, ztilde, power_of<2>(ztilde), power_of<3>(ztilde), power_of<4>(ztilde) }};
            std::array<double, 5> lz{{ 1.0, lnz, lnz * lnz, lnz * lnz * lnz, lnz * lnz * lnz * lnz }};
            std::array<double, 5> dd{{ 1.0, deltatilde, power_of<2>(deltatilde), power_of<3>(deltatilde), power_of<4>(deltatilde) }};
            std::array<double, 5> ld{{ 1.0, lndelta, power_of<2>(lndelta), power_of<3>(lndelta), power_of<4>(lndelta) }};

            // f_22, f_27 are computed from a fit to polynomials of the form
            // (1 - z)^i ln(z)^j delta^k ln(1 - delta)^l, with j <= i, l <= k and i,j <= 4
            // Results are valid for 0.03 <= z <= 0.12 and 1/36 <= delta <= 1/3.
            double f22int1 = 0.0, f22int2 = 0.0, f27int1 = 0.0, f27int2 = 0.0;
            {
                unsigned idx = 0;
                for (unsigned i = 0 ; i <= 4 ; ++i)
                {
                    for (unsigned j = 0 ; j <= i ; ++j)
                    {
                        for (unsigned k = 0 ; k <= 4 ; ++k)
                        {
                            for (unsigned l = 0 ; l <= k ; ++l)
                            {
                                f22int1 += f22int1coeffs[idx] * zz[i] * lz[j] * dd[k] * ld[l];
                                f22int2 += f22int2coeffs[idx] * zz[i] * lz[j] * dd[k] * ld[l];
                                f27int1 += f27int1coeffs[idx] * zz[i] * lz[j] * dd[k] * ld[l];
                                f27int2 += f27int2coeffs[idx] * zz[i] * lz[j] * dd[k] * ld[l];

                                ++idx;
                            }
                        }
                    }
                }
            }

            Fij result;

            // cf. [CMM:1996A], f_22 - f_28, Eq. (39), p. 9
            result.f22 = 16.0 * z / 27.0 * (delta * f22int1 + f22int2);
            result.f27 = -8.0 * z * z / 9.0 * (delta * f27int1 + f27int2);
            result.f28 = -result.f27;
            // cf. [CMM:1996A], f_11 - f_18, below Eq. (40), p. 9
            result.f11 = result.f22 / 36.0;
            result.f12 = -result.f22 / 3.0;
            result.f17 = -result.f27 / 6.0;
            result.f18 = -result.f28 / 6.0;
            // cf. [CMM:1996A], f_77, Eq. (37), p. 8
            result.f77 = 10.0 * delta / 3.0 + delta2 / 3.0 - 2.0 * delta3 / 9.0 + delta * (delta - 4.0) / 3.0 * lndelta;
            // cf. [CMM:1996A], f_78, Eq. (39), p. 9
            result.f78 = 8.0 / 9.0 * (dl1mdelta - pi2 / 6.0 - delta * lndelta + 9.0 * delta / 4.0 - delta2 / 4.0 + delta3 / 12.0);
            // cf. [CMM:1996A], f_88, Eq. (38), pp. 8-9
            result.f88 = (-2.0 * std::log(m_b_MSbar() / model->m_s_msbar(mu())) * (delta2 + 2.0 * delta + 4.0 * ln1mdelta)
                    + 4.0 * dl1mdelta - 2.0 * pi2 / 3.0 - delta * (2.0 + delta) * lndelta + 8.0 * ln1mdelta
                    - 2.0 / 3.0 * delta3 + 3.0 * delta2 + 7.0 * delta) / 27.0;

            return result;
        }

        // Auxiliary function to compute the functions r_i
        inline std::array<complex<double>, 8> r_i(const double & z) const
        {
            // Constants
            static const double pi = M_PI, pi2 = pi * pi;
            static const double zeta3 = 1.2020569031595943;
            static const complex<double> iu(0.0, 1.0);

            // kinematics
            double z2 = z * z, z3 = z2 * z, z4 = z2 * z2, z5 = z3 * z2, z6 = z3 * z3;
            double sqrtz = std::sqrt(z), L = std::log(z), L2 = L * L, L3 = L2 * L;

            // Auxilliary function X_b, cf. [BCMU:2002A], Eq. (3.2), p. 3
            static const double Xb = -0.1684408;
            // Auxilliary function a and b for bottom quarks (z = 1)
            // cf. [BCMU:2002A], Eqs. (3.6) and (3.7), p. 4
            static complex<double> a1(4.0859, 4.0 / 9.0 * pi);
            static complex<double> b1(0.0316, 4.0 / 81.0 * pi);
            // Auxilliary function a and b for charm quarks, z ~= 1/4,
            // as an expansion in z^n ln(z)^m
            // cf. [BCMU:2002A], Eqs. (3.8) and (3.9), pp. 4, 5
            complex<double> az = 16.0 / 9.0 * complex<double>(
                        z * (5.0 / 2.0 - pi2 / 3.0 - 3.0 * zeta3 + (5.0 / 2.0 - 3.0 * pi2 / 4.0) * L + L2 / 4.0 + L3 / 12.0)
                        + z2 * (7.0 / 4.0 + 2.0 / 3.0 * pi2 - pi2 / 2.0 * L - L2 / 4.0 + L3 / 12.0)
                        + z3 * (-7.0 / 6.0 - pi2 / 4.0 + 2.0 * L - 3.0 / 4.0 * L2)
                        + z4 * (457.0 / 216.0 - 5.0 / 18.0 * pi2 - 1.0 / 72.0 * L - 5.0 / 6.0 * L2)
                        + z5 * (35101.0 / 8640.0 - 35.0 / 72.0 * pi2 - 185.0 / 144.0 * L - 35.0 / 24.0 * L2)
                        + z6 * (67801.0 / 8000.0 - 21.0 / 20.0 * pi2 - 3303.0 / 800.0 * L - 63.0 / 20.0 * L2),
                pi * (  z * (2.0 - pi2 / 6.0 + 0.5 * L + 0.5 * L2)
                        + z2 * (0.5 - pi2 / 6.0 - L + 0.5 * L2)
                        + z3
                        + z4 * 5.0 / 9.0
                        + z5 * 49.0 / 72.0
                        + z6 * 231.0 / 200.0));
            complex<double> bz = -8.0 / 9.0 * complex<double>(
                        z * (-3.0 + pi2 / 6.0 - L)
                        + z * sqrtz * (-2.0 / 3.0 * pi2)
                        + z2 * (0.5 + pi2 - 2.0 * L - 0.5 * L2)
                        + z3 * (-25.0 / 12.0 - pi2 / 9.0 - 19.0 / 18.0 * L + 2.0 * L2)
                        + z4 * (-1376.0 / 225.0 + 137.0 / 30.0 * L + 2.0 * L2 + 2.0 / 3.0 * pi2)
                        + z5 * (-131317.0 / 11760.0 + 887.0 / 84.0 * L + 5.0 * L2 + 5.0 / 3.0 * pi2)
                        + z6 * (-2807617.0 / 97200.0 + 16597.0 / 540.0 * L + 14.0 * L2 + 14.0 / 3.0 * pi2),
                pi * (  z * -1.0
                        + z2 * (1.0 - 2.0 * L)
                        + z3 * (-10.0 / 9.0 + 4.0 / 3.0 * L)
                        + z4
                        + z5 * 2.0 / 3.0
                        + z6 * 7.0 / 9.0));

            std::array<complex<double>, 8> result
            {{
                833.0 / 729.0 - (az + bz) / 3.0 + iu * 40.0 / 243.0 * pi,
                -1666.0 / 243.0 + 2.0 * (az + bz) - iu * 80.0 / 81.0 * pi,
                2392.0 / 243.0 + 8.0 * pi / 3 / std::sqrt(3.0) + 32.0 / 9.0 * Xb - a1 + 2.0 * b1 + iu * 56.0 / 81.0 * pi,
                -761.0 / 729.0 - 4.0 * pi / 9.0 / std::sqrt(3.0) - 16.0 / 27.0 * Xb + a1 / 6.0 + 5.0 * b1 / 3.0 + 2.0 * bz - iu * 148.0 / 243.0 * pi,
                56680.0 / 243.0 + 32.0 * pi / 3.0 / std::sqrt(3.0) + 128.0 / 9.0 * Xb - 16.0 * a1 + 32.0 * b1 + iu * 896.0 / 81.0 * pi,
                5710.0 / 729.0 - 16.0 * pi / 9.0 / std::sqrt(3.0) - 64.0 / 27.0 * Xb - 10.0 * a1 / 3.0 + 44.0 * b1 / 3.0 + 12.0 * az + 20.0 * bz - iu * 2296.0 / 243.0 * pi,
                -10.0 / 3.0 - 8.0 / 9.0 * pi2,
                44.0 / 9.0 - 8.0 / 27.0 * pi2 + iu * 8.0 / 9.0 * pi
            }};

            return result;
        }

        inline double phase_space_g(const double & z) const
        {
            double z2 = z * z, z3 = z2 * z, z4 = z3 * z, lnz = std::log(z);

            return 1.0 - 8.0 * z + 8.0 * z3 - z4 - 12.0 * z2 * lnz;
        }

        inline double phase_space_h(const double & z) const
        {
            static const double pi = M_PI, pi2 = pi * pi;
            double z2 = z * z, z3 = z2 * z, z4 = z3 * z, lnz = std::log(z), sqrtz = std::sqrt(z);

            return -(1.0 - z2) * (25.0 / 4.0 - 239.0/3.0 * z + 25.0 / 4.0 * z2)
                + z * lnz * (20.0 + 90.0 * z - 4.0/3.0 * z2 + 17.0/3.0 * z3)
                + z2 * lnz * lnz * (36.0 + z2)
                + (1.0 - z2) * (17.0 / 3.0 - 64.0 / 3.0 * z + 17.0 / 3.0 * z2) * std::log(1.0 - z)
                - 4.0 * (1.0 + 30.0 * z2 + z4) * lnz * std::log(1.0 - z)
                - (1.0 + 16.0 * z2 + z4) * (6.0 * std::real(dilog(z)) - pi2)
                - 32. * z * std::sqrt(z) * (1.0 + z) * (pi * pi - 4.0 * std::real(dilog(sqrtz)) + 4.0 * std::real(dilog(-sqrtz)) - 2.0 * lnz * std::log((1.0 - sqrtz)/(1.0 + sqrt(z))));
        }

        // effective Wilson coefficients c7eff according to [CMM:1996A], p. 2, Eq. (5)
        inline complex<double> c7eff(const WilsonCoefficients<BToS> & wc) const
        {
            return wc.c7() - wc.c3() / 3.0 - 4.0 * wc.c4() / 9.0 - 20.0 * wc.c5() / 3.0 - 80.0 * wc.c6() / 9.0;
        }

        // effective Wilson coefficient c8eff according to [CMM:1996A], p. 2, Eq. (5)
        inline complex<double> c8eff(const WilsonCoefficients<BToS> & wc) const
        {
            return wc.c8() + wc.c3() - wc.c4() / 6.0 + 20.0 * wc.c5() - 10.0 * wc.c6() / 3.0;
        }

        // D: perturbative b->sgamma contribution, cf. [CMM:1996A], p. 7, Eq. (31)
        inline complex<double> perturbative_bsgamma(const double & z, const WilsonCoefficients<BToS> & wc, const double & alpha_s, const double & lnmu) const
        {
            static const double pi = M_PI;

            // Elements of the anomalous mass dimension gamma_i7, cf. [CMM:1996A], Eq. (8), p. 3
            static const std::array<double, 8> gamma7
            {{
                 -208.0 / 243.0, 416.0 / 81.0, -176.0 / 81.0, -152.0 / 243.0, -6872.0 / 81.0, 4624.0 / 243.0, 32.0 / 3.0, -32.0 / 9.0
            }};

            // Functions r_i
            std::array<complex<double>, 8> r = r_i(z);

            // Reduced strong coupling
            double a_s = alpha_s / (4.0 * pi);

            // Partonic contribution, cf. [CMM:1996A], Eq. (31), p. 7
            complex<double> D = c7eff(wc);
            D += a_s * wc.c1() * (r[1 - 1] + gamma7[1 - 1] * lnmu);
            D += a_s * wc.c2() * (r[2 - 1] + gamma7[2 - 1] * lnmu);
            D += a_s * wc.c3() * (r[3 - 1] + gamma7[3 - 1] * lnmu);
            D += a_s * wc.c4() * (r[4 - 1] + gamma7[4 - 1] * lnmu);
            D += a_s * wc.c5() * (r[5 - 1] + gamma7[5 - 1] * lnmu);
            D += a_s * wc.c6() * (r[6 - 1] + gamma7[6 - 1] * lnmu);
            D += a_s * c7eff(wc) * (r[7 - 1] + gamma7[7 - 1] * lnmu);
            D += a_s * c8eff(wc) * (r[8 - 1] + gamma7[8 - 1] * lnmu);

            return D;
        }

        // Dprime: perturbative b->sgamma contribution, chirality flipped
        // based on the results given in [CMM:1996A], p. 7, Eq. (31)
        inline complex<double> perturbative_bsgamma_prime(const double & z, const WilsonCoefficients<BToS> & wc, const double & alpha_s, const double & lnmu) const
        {
            static const double pi = M_PI;

            // Elements of the anomalous mass dimension gamma_i7, cf. [CMM:1996A], Eq. (8), p. 3
            static const std::array<double, 8> gamma7
            {{
                 -208.0 / 243.0, 416.0 / 81.0, -176.0 / 81.0, -152.0 / 243.0, -6872.0 / 81.0, 4624.0 / 243.0, 32.0 / 3.0, -32.0 / 9.0
            }};

            // Functions r_i
            std::array<complex<double>, 8> r = r_i(z);

            // Reduced strong coupling
            double a_s = alpha_s / (4.0 * pi);

            // Partonic contribution, cf. [CMM:1996A], Eq. (31), p. 7
            // use chirality flipped Wilson coefficient only and
            // keep c1 to c6 at their SM values (i.e. 0).
            complex<double> Dprime = wc.c7prime();
            Dprime += a_s * wc.c7prime() * (r[7 - 1] + gamma7[7 - 1] * lnmu);
            Dprime += a_s * wc.c8prime() * (r[8 - 1] + gamma7[8 - 1] * lnmu);

            return Dprime;
        }

        // A: perturbative b->sgammagluon contribution, cf. [CMM:1996A], p. 7, Eq. (32)
        inline double perturbative_bsgammagluon(const double & z, const double & delta, const WilsonCoefficients<BToS> & wc, const double & alpha_s) const
        {
            static const double pi = M_PI;
            double lndelta = std::log(delta);

            // Functions f_ij
            Fij fij = f_ij(z, std::log(z), delta);

            // Reduced strong coupling
            double a_s = alpha_s / (4.0 * pi);

            // b->sgammagluon Bremstrahlung corrections, cf. [CMM:1996A], Eq. (32), p. 7
            double A = (std::exp(-alpha_s * lndelta * (7.0 + 2.0 * lndelta) / (3.0 * pi)) - 1.0) * (std::norm(c7eff(wc)) + std::norm(wc.c7prime()));
            A += 4.0 * a_s * std::norm(wc.c1()) * std::real(fij.f11);
            A += 4.0 * a_s * std::real(wc.c1() * std::conj(wc.c2()) * fij.f12);
            A += 4.0 * a_s * std::real(wc.c1() * std::conj(wc.c7()) * fij.f17);
            A += 4.0 * a_s * std::real(wc.c1() * std::conj(wc.c8()) * fij.f18);
            A += 4.0 * a_s * std::norm(wc.c2()) * std::real(fij.f22);
            A += 4.0 * a_s * std::real(wc.c2() * std::conj(wc.c7()) * fij.f27);
            A += 4.0 * a_s * std::real(wc.c2() * std::conj(wc.c8()) * fij.f28);
            A += 4.0 * a_s * std::norm(wc.c7()) * std::real(fij.f77);
            A += 4.0 * a_s * std::real(wc.c7() * std::conj(wc.c8()) * fij.f78);
            A += 4.0 * a_s * std::norm(wc.c8()) * std::real(fij.f88);
            A += 4.0 * a_s * std::norm(wc.c7prime()) * std::real(fij.f77);
            A += 4.0 * a_s * std::norm(wc.c8prime()) * std::real(fij.f88);
            A += 4.0 * a_s * std::real(wc.c7prime() * std::conj(wc.c8prime()) * fij.f78);

            return A;
        }

        // ToDo: this function is only used in Diagnostics
        double ratio_quark(const double & z, const double & delta, const double & ckm, const double & alpha_s, const double & mu) const
        {
            // Constants
            static const double pi = M_PI;
            static const complex<double> iu(0.0, 1.0);

            // masses
            double m_b_pole = 4.8;                  // ToDo: check whether hardcoded numerical value could be removed
            double lnmu = std::log(m_b_pole / mu);

            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu, LeptonFlavor::muon /* fake lepton flavor */);

            // Perturbative contributions
            complex<double> D = perturbative_bsgamma(z, wc, alpha_s, lnmu);
            complex<double> Dprime = perturbative_bsgamma_prime(z, wc, alpha_s, lnmu);
            double A = perturbative_bsgammagluon(z, delta, wc, alpha_s);

            // phase space factors, cf. [CMM:1996A], between Eqs. (43) and (44), p. 10
            double g = phase_space_g(z), h = phase_space_h(z);

            // cf. [CMM:1996A], Eq. (43), p. 10
            double kappa = 1.0 - (2.0 * alpha_s * h) / (3.0 * pi * g);
            // cf. [CMM:1996A], Eq. (45), p. 10
            double F = (1.0 - (8.0 * alpha_s) / (3.0 * pi)) / kappa;
            // cf. [CMM:1996A], Eq. (30), p. 7
            double result = power_of<2>(ckm)
                * 6.0 * alpha_e / (pi * g) * F
                * (std::norm(D) + std::norm(Dprime) + A);

            return result;
        }

        double decay_rate(const double & emin) const
        {
            // Constants
            static const double pi = M_PI, pi2 = pi * pi, pi4 = pi2 * pi2;
            static const complex<double> iu(0.0, 1.0);

            // masses
            double m_b_msbar = model->m_b_msbar(mu());
            double m_b_kin = model->m_b_kin(1.0); // for the HQE correction we use the kinetic mass at mu_kin = 1GeV
            double m_b_pole = model->m_b_pole();
            double m_c = model->m_c_msbar(mu());
            double lnmu = std::log(m_b_MSbar() / mu());

            // Kinematics
            // Use m_c(mu) instead of m_c_pole, cf. [GM:2001A], pp. 2-3.
            double m_c_hat = m_c / m_b_msbar;
            double z = m_c_hat * m_c_hat;

            // cutoff parameter delta
            double delta = 1.0 - 2.0 * emin / m_b_pole, delta2 = delta * delta, delta3 = delta * delta2, delta4 = delta2 * delta2;
            double lndelta = std::log(delta), ln2delta = lndelta * lndelta;

            // We only support evaluations with 0.03 <= z <= 0.12 and 1/36 <= delta <= 1/3, see also calculation of function f_ij(z, delta).
            if (((z < 0.03) || (0.12 < z))
                || ((delta < 1.0 / 36.0) || (1.0 / 3.0 < delta)))
            {
                throw InternalError("Photon Energy cutoff E_min yields z = " + stringify(z, 4u) + ", delta = " + stringify(delta, 4u) + ", which is unsupported by the current implementation of B->X_sgamma::BR(E_min)");
            }

            // Strong coupling
            double alpha_s = model->alpha_s(mu()), a_s = alpha_s / (4.0 * pi);

            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), LeptonFlavor::muon /* fake lepton flavor */);

            // Perturbative contributions
            complex<double> D = perturbative_bsgamma(z, wc, alpha_s, lnmu);
            complex<double> Dprime = perturbative_bsgamma_prime(z, wc, alpha_s, lnmu);
            double A = perturbative_bsgammagluon(z, delta, wc, alpha_s);

            /*
             * We use lambda_1 = -mu^2_pi and lambda_2 = mu^2_G / 3.0 and neglect higher order
             * terms in the 1/m_b expansion. For the alpha_s Lambda / m_b corrections we use the
             * findings of [EGN:2010A], p. 10, eq. (3.17).
             *
             * The 1/m_b terms to leading order in alpha_s are not used, since they are
             * already included in the NLO calculation by [CMM:1996A].
             *
             * The remainder enters the rate as (|C7eff|^2 + |C7'|^2) * hqe.
             */
            double hqe = /*1.0 + a_s * (16.0 / 9.0 * (4.0 - pi2 + 3.0 * lnmu) - 8.0 / 3.0 * ln2delta
                    - 4.0 / 3.0 * (7.0 + 4.0 * delta - delta2) * lndelta - 4.0 / 9.0 * (31.0 - 30.0 * delta - 3.0 * delta2 + 2.0 * delta3))*/
                - (0.5 + a_s * (8.0 / 9.0 * (4.0 - pi2 + 3.0 * lnmu) - 4.0 / 3.0 * ln2delta
                            - 2.0 * (delta - 1.0) * (-3.0 + 14.0 * delta - 46.0 * delta2 + 4.0 * delta3) / (9.0 * delta2)
                            - 2.0 * (4.0 - 12.0 * delta + 27.0 * delta2 + 14.0 * delta3 - 3.0 * delta4) / (9.0 * delta2) * lndelta)) * (mu2_pi / (m_b_kin * m_b_kin))
                + (-9.0 / 2.0 + a_s * (-(87.0 + 32.0 * pi2 - 27.0 * lnmu) / 9.0 - 16.0 / 3.0 * ln2delta
                            + (-1.0 + 20.0 / delta - 4.0 * delta - 15.0 * delta2)
                            + (13.0 + 8.0 / delta + 92.0 / 3.0 * delta + 7.0 / 3.0 * delta2) * lndelta)) * (mu2_g / (m_b_kin * m_b_kin));

            // cf. [CMM:1996A], Eq. (30), p. 7
            double lambda_t2 = std::norm(model->ckm_tb() * conj(model->ckm_ts()));
            double result = power_of<2>(gfermi()) * alpha_e / (32.0 * pi4) * lambda_t2
                * power_of<3>(m_b_pole) * power_of<2>(model->m_b_msbar(mu()))
                * (std::norm(D) + std::norm(Dprime) + A + (std::norm(c7eff(wc)) + std::norm(wc.c7prime())) * hqe);

            return result;
        }

        double branching_ratio(const double & emin) const
        {
            return decay_rate(emin) * tau() / hbar();
        }

        double photon_energy_moment_1(const double & emin) const
        {
            // Calculation is based on results by [EGN:2010A].

            const double a = model->alpha_s(mu()) / (4.0 * M_PI);
            const double m_b = model->m_b_kin(1.0), m_b2 = m_b * m_b;
            const double z = 2.0 * emin / m_b, z2 = z * z, z3 = z2 * z, z4 = z2 * z2;
            const double ln1mz = std::log(1.0 - z), lnmu = std::log(mu() / m_b);

            double c_0  = 1.0 + a * (-46.0 / 27.0 + 8.0 / 9.0 * (8.0 - 9.0 * z + z3) * ln1mz
                    + 2.0 / 27.0 * z * (96.0 - 60.0 * z - 22.0 * z2 + 9.0 * z3));
            double c_pi = 0.5
                + a * (z * (96.0 - 120.0 * z - 52.0 * z2 + 25.0 * z3 - 9.0 * z4) / (27.0 * (1.0 - z))
                        - 23.0 / 27.0 + 8.0 * (4.0 - 7.0 * z + z3) / (9.0 * (1.0 - z)) * ln1mz);
            double c_g  = -0.5 + a * (-565.0 / 81.0 - 3.0 * lnmu
                    + 2.0 * (80.0 - 33.0 * z - 81.0 * z2 + 25.0 * z3) / 27.0 * ln1mz
                    + z * (480.0 + 42.0 * z - 425.0 * z2 + 81.0 * z3) / 81.0);

            return m_b / 2.0 * (c_0 + c_pi * mu2_pi / m_b2 + c_g * mu2_g / m_b2);
        }

        double photon_energy_moment_2(const double & emin) const
        {
            // Calculation is based on results by [EGN:2010A].

            const double a = model->alpha_s(mu()) / (4.0 * M_PI);
            const double m_b = model->m_b_kin(1.0), m_b2 = m_b * m_b;
            const double z = 2.0 * emin / m_b, z2 = z * z, z3 = z2 * z, z4 = z2 * z2, z5 = z4 * z;
            const double ln1mz = std::log(1.0 - z);

            double c_0  = a * (-2.0 / 9.0 * (1.0 - 2.0 * z + z2) * (17.0 - 2.0 * z - 3.0 * z2) * ln1mz
                        + 1.0 / 270.0 * (1.0 - 2.0 * z + z2) * (61.0 - 898.0 * z - 207.0 * z2 + 144.0 * z3));
            double c_pi = 1.0 / 3.0 + a * (-4.0 / 27.0 * (2.0 - 27.0 * z + 7.0 * z3) * ln1mz
                        - 2.0 / 405.0 * (169.0 + 60.0 * z - 780.0 * z2 - 385.0 * z3 + 225.0 * z4 - 54.0 * z5));
            double c_g  = a * ((43.0 - 84.0 * z + 258.0 * z2 - 292.0 * z3 + 75.0 * z4) * ln1mz / 54.0
                        - (707.0 - 2580.0 * z + 3750.0 * z2 - 13820.0 * z3 + 14535.0 * z4 - 2592.0 * z5) / 3240.0);

            return m_b2 / 4.0 * (c_0 + c_pi * mu2_pi / m_b2 + c_g * mu2_g / m_b2);
        }

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            // Functions f_ij, cf. [CMM:1996A], p. 9
            {
                // z = 0.22^2, delta = 1/6
                Fij f_ij = this->f_ij(0.22 * 0.22, std::log(0.22 * 0.22), 1.0 / 6.0);
                results.add(Diagnostics::Entry{ std::real(f_ij.f22), "Re(f_22(z = 0.22^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f22), "Im(f_22(z = 0.22^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::real(f_ij.f27), "Re(f_27(z = 0.22^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f27), "Im(f_27(z = 0.22^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::real(f_ij.f78), "Re(f_78(z = 0.22^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f78), "Im(f_78(z = 0.22^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                // z = 0.29^2, delta = 1/6
                f_ij = this->f_ij(0.29 * 0.29, std::log(0.29 * 0.29), 1.0 / 6.0);
                results.add(Diagnostics::Entry{ std::real(f_ij.f22), "Re(f_22(z = 0.29^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f22), "Im(f_22(z = 0.29^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::real(f_ij.f27), "Re(f_27(z = 0.29^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f27), "Im(f_27(z = 0.29^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::real(f_ij.f78), "Re(f_78(z = 0.29^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f78), "Im(f_78(z = 0.29^2, delta = 1/6)), [CMM:1996A], p. 9, Eq. (39)" });
                // z = 0.22^2, delta = 1/9
                f_ij = this->f_ij(0.22 * 0.22, std::log(0.22 * 0.22), 1.0 / 9.0);
                results.add(Diagnostics::Entry{ std::real(f_ij.f22), "Re(f_22(z = 0.22^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f22), "Im(f_22(z = 0.22^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::real(f_ij.f27), "Re(f_27(z = 0.22^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f27), "Im(f_27(z = 0.22^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::real(f_ij.f78), "Re(f_78(z = 0.22^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f78), "Im(f_78(z = 0.22^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                // z = 0.29^2, delta = 1/9
                f_ij = this->f_ij(0.29 * 0.29, std::log(0.29 * 0.29), 1.0 / 9.0);
                results.add(Diagnostics::Entry{ std::real(f_ij.f22), "Re(f_22(z = 0.29^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f22), "Im(f_22(z = 0.29^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::real(f_ij.f27), "Re(f_27(z = 0.29^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f27), "Im(f_27(z = 0.29^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::real(f_ij.f78), "Re(f_78(z = 0.29^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
                results.add(Diagnostics::Entry{ std::imag(f_ij.f78), "Im(f_78(z = 0.29^2, delta = 1/9)), [CMM:1996A], p. 9, Eq. (39)" });
            }

            // Function r_i, cf. [CMM:1996A], p. 8, Eqs. (33)-(36) and [BCMU:2002A], p. 3, Eq. (3.1).
            {
                // z = 0.22^2
                std::array<complex<double>, 8> r_i = this->r_i(0.22 * 0.22);
                results.add(Diagnostics::Entry{ std::real(r_i[0]), "Re(r_1(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[0]), "Im(r_1(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[1]), "Re(r_2(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[1]), "Im(r_2(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[2]), "Re(r_3(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[2]), "Im(r_3(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[3]), "Re(r_4(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[3]), "Im(r_4(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[4]), "Re(r_5(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[4]), "Im(r_5(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[5]), "Re(r_6(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[5]), "Im(r_6(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[6]), "Re(r_7(z = 0.22^2)), [CMM:1996A], p. 8, Eq. (36)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[6]), "Im(r_7(z = 0.22^2)), [CMM:1996A], p. 8, Eq. (36)" });
                results.add(Diagnostics::Entry{ std::real(r_i[7]), "Re(r_8(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[7]), "Im(r_8(z = 0.22^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });

                // z = 0.29^2
                r_i = this->r_i(0.29 * 0.29);
                results.add(Diagnostics::Entry{ std::real(r_i[0]), "Re(r_1(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[0]), "Im(r_1(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[1]), "Re(r_2(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[1]), "Im(r_2(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[2]), "Re(r_3(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[2]), "Im(r_3(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[3]), "Re(r_4(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[3]), "Im(r_4(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[4]), "Re(r_5(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[4]), "Im(r_5(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[5]), "Re(r_6(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[5]), "Im(r_6(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::real(r_i[6]), "Re(r_7(z = 0.29^2)), [CMM:1996A], p. 8, Eq. (36)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[6]), "Im(r_7(z = 0.29^2)), [CMM:1996A], p. 8, Eq. (36)" });
                results.add(Diagnostics::Entry{ std::real(r_i[7]), "Re(r_8(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
                results.add(Diagnostics::Entry{ std::imag(r_i[7]), "Im(r_8(z = 0.29^2)), [BCMU:2002A], p. 3, Eq. (3.1)" });
            }

            // Phase space functions g and h, cf. [CMM:1996A], p. 10, Eqs. (42),(43)f
            {
                // g
                results.add(Diagnostics::Entry{ phase_space_g(0.22 * 0.22), "h(z = 0.22^2), CMM[:1996], p. 10, below Eq. (43)" });
                results.add(Diagnostics::Entry{ phase_space_g(0.29 * 0.29), "h(z = 0.29^2), CMM[:1996], p. 10, below Eq. (43)" });

                // h
                results.add(Diagnostics::Entry{ phase_space_h(0.22 * 0.22), "h(z = 0.22^2), CMM[:1996], p. 10, below Eq. (43)" });
                results.add(Diagnostics::Entry{ phase_space_h(0.29 * 0.29), "h(z = 0.29^2), CMM[:1996], p. 10, below Eq. (43)" });
            }

            // R_quark, cf. [CMM:1996A], p. 10, Eq. (44)
            {
                results.add(Diagnostics::Entry{ ratio_quark(0.29 * 0.29, 0.29 * 0.29, 0.976, 0.212, 5.0), "R_quark(z = 0.29^, delta = z,   ckm = 0.976, alpha_s = 0.212, mu_b = 5.0), CMM[:1996], p. 10, Eq. (44)" });
                results.add(Diagnostics::Entry{ ratio_quark(0.29 * 0.29, 1.0 / 3.0,   0.976, 0.212, 5.0), "R_quark(z = 0.29^, delta = 1/3, ckm = 0.976, alpha_s = 0.212, mu_b = 5.0), CMM[:1996], p. 10, Eq. (44)" });
            }

            return results;
        }
    };

    BToXsGamma<NLO>::BToXsGamma(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToXsGamma<NLO>>(new Implementation<BToXsGamma<NLO>>(parameters, options, *this))
    {
    }

    BToXsGamma<NLO>::~BToXsGamma()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<BToXsGamma<NLO>>::options
    {
    };

    double
    BToXsGamma<NLO>::integrated_branching_ratio(const double & emin) const
    {
        return _imp->branching_ratio(emin);
    }

    double
    BToXsGamma<NLO>::photon_energy_moment_1(const double & emin) const
    {
        return _imp->photon_energy_moment_1(emin);
    }

    double
    BToXsGamma<NLO>::photon_energy_moment_2(const double & emin) const
    {
        return _imp->photon_energy_moment_2(emin);
    }

    Diagnostics
    BToXsGamma<NLO>::diagnostics() const
    {
        return _imp->diagnostics();
    }

    const std::set<ReferenceName>
    BToXsGamma<NLO>::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BToXsGamma<NLO>::begin_options()
    {
        return Implementation<BToXsGamma<NLO>>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToXsGamma<NLO>::end_options()
    {
        return Implementation<BToXsGamma<NLO>>::options.cend();
    }

    const std::array<double, 225>
    Implementation<BToXsGamma<NLO>>::f22int1coeffs
    {{
   1.4646809648996456491,-13.624222527147306666,-5.6735048900322766854,28.731416748209033946,
   134.11885755958407797,62.980440698647723107,-132.38425038219737408,166.6719677392034733,
   143.23305924837139056,-84.495279334261267381,128.78261323930827656,-467.31044471411214904,
   427.79783239077857985,-139.45029294942700437,35.246235394683755229,-487428.7525669911977,
   1.3688328663169634371e6,1.2555750144023651334e6,1.551017490629471732e7,-2.9503034942111583783e7,
   3.2403193212802561119e7,-7.0604921470993827594e7,4.4074805036709588197e7,1.9905654065362444472e7,
   -5.4527124383366679622e7,-4.2556567927661316963e8,2.1504359603506831713e8,
   -4.7813219760450604542e8,2.9534874719608671685e8,-3.8993078578609304285e8,-192972.30908159378754,
   542217.3957487121912,497228.45744635103296,6.1402812082787097079e6,-1.1683573379419080817e7,
   1.2828021335367750141e7,-2.7951613529757662217e7,1.7446969087022672273e7,7.8781516686304468667e6,
   -2.1587079998754752005e7,-1.6849490679541473457e8,8.5150752650034017066e7,
   -1.8931293451696541793e8,1.1693880697471150221e8,-1.5438412539117890239e8,
   -4.5403249816098870184e6,2.5189998147515637443e7,1.8665480652432437982e7,-7.6637038026835126038e7,
   -1.1586704161429928097e8,3.9960332410053616548e8,2.7014191310254451655e8,-6.0885447246837564768e8,
   -1.9373879681383978165e9,-1.585355348385347948e9,2.2052959217402068758e9,1.6031664797902217204e10,
   5.4204783333761088527e9,6.6173293766176602282e9,-1.5311093856770538968e9,-1.8128922983329924627e6,
   7.4442395339094751284e6,6.6727843734876070924e6,8.2664136786099777666e6,-7.7593793748398307909e7,
   1.4400461355975997163e8,-6.9681203513298058244e7,9.5961168062621687891e6,-3.8321435419478807074e8,
   -5.1253336701171174485e8,5.2369821356631440394e8,3.065333482750453467e9,-2.2018018343889167744e9,
   1.1753980947416151061e9,-8.3067108615143376749e8,372144.76223244445904,-2.0532579470341754939e6,
   -1.2537386798820888803e6,3.2309911385994248652e6,1.0231709421835867035e7,-3.0810594799512386565e7,
   -1.5057656208522030344e7,6.4919149685217211347e7,1.364177271170432464e8,8.7953470160201261723e7,
   1.9170372771522122939e8,-1.4658959251269929115e9,-1.3508272536841797085e9,-8.009681406039865469e8,
   2.1361196155481586398e8,-1.9066314546516021369e8,-1.1786640890090446815e8,2.0737534675248092751e8,
   1.6423616384065070167e10,-1.97267925669847962e10,7.6449485819218801377e9,
   -8.0412091087073146488e10,7.510857403420876283e10,1.2159290634234590627e11,
   2.6108738869890329702e10,-9.2897932132034261358e10,5.6050747467484260855e10,
   -1.9745227489305817278e11,1.7965241746853841685e11,-1.8242690625243633653e11,
   -9.8396457234086106407e7,2.4547753218101455906e8,2.4087770740999177633e8,3.517626575635808199e9,
   -6.191938002559500497e9,6.3728475641916155065e9,-1.6488098554103479686e10,
   1.1088103184086401872e10,8.8018479858643064882e9,-9.1676536844650705551e9,
   -7.0872252078862577382e10,7.1490539042216298862e10,-6.8947096498339907615e10,
   7.492046128958219303e10,-7.753594642131302738e10,-8.6983754319659787052e6,1.40256927444037902e8,
   7.1762001420175777546e7,-1.5273264280251153183e9,8.7526242406822606207e8,1.4348079725065685737e9,
   7.8890064396318360194e9,-9.349920902238324147e9,-1.9429963657885575609e10,
   -9.0922604099973744212e9,-1.9849168380025999551e10,1.4270989511086894437e10,
   1.7010411032189988563e9,-2.7489851451327492364e9,-1.6744102199864682181e9,
   -1.7817949603042625525e6,2.0889369450963031612e7,1.1185866110564805257e7,-1.7273719153193572008e8,
   4.9582026957180701077e7,2.1401008618984930712e8,8.5746372638056227673e8,-1.1095265906767746467e9,
   -2.395675343064241287e9,-1.1802192836193605116e9,-2.0422988375025694531e9,5.630864069165689324e9,
   3.7315136912122006803e9,1.8006521469708337534e9,-5.6369631758708288717e8,2.1323545098421114149e8,
   8.7493247426399529921e8,8.9063961054627824228e7,-1.7877229503808826806e10,
   1.1179931157060179839e10,-1.3975335567619231248e10,8.6367721297721445066e10,
   -8.1147472992784131801e10,-1.1976145279174174133e11,-1.2575301402641467984e10,
   9.8273473003494764662e10,-6.8996254183912960642e10,1.6378503182298043458e11,
   -1.8986292255800946939e11,1.9208184450615293338e11,-1.2960672059949865369e8,
   -1.1174053554503666795e8,1.2951288256846997602e8,1.858091411061559432e10,
   -2.5359173243629927458e10,-1.3824461577413794925e9,-9.3460800817285399628e10,
   9.3805566859923421529e10,1.6959728054949700912e11,5.9073316960243211499e10,
   -1.5703376289662421476e10,-2.212225309937900284e9,-1.9750116188437706687e11,
   1.6242770939702071247e11,-1.5105042845190688938e11,3.6155049025646010823e7,
   -2.0088517866437914078e8,-1.3610955792962486094e8,-1.3445033321354164141e9,
   3.8022686218477651876e9,-1.5555048760094114989e9,6.5029264072742232377e9,-3.888821466696554761e9,
   -6.2396373868235471202e9,5.7385355857242091451e8,1.6997967343666454748e10,
   -5.3919780687782454378e10,5.6346008999799632186e9,-4.2360692671688373601e10,
   2.8500259465929803004e10,-5.7523317174590816592e6,-2.3163703515028615355e8,
   -9.2802138705680886226e7,2.7987611189419874401e9,-1.3303743965325819546e9,
   -3.4121484649980382941e8,-1.4072954591910031725e10,1.5410723389485441479e10,
   2.7356554779841228867e10,9.5868720468734264438e9,9.3840405630431729336e9,
   -1.7771085811689197649e10,-1.4731475582242191575e10,1.2539907120410885245e10,
   -1.1874698107253276016e10,542692.80405219101139,-3.7436766267189895728e7,-1.6829596012936111609e7,
   3.5405174885418782973e8,-8.7576228902784927168e7,-1.2707519640721163037e8,
   -1.7489021844424932833e9,2.0176383206432090245e9,3.8198112671195969598e9,1.4882541844281389974e9,
   1.7599809838683808282e9,-5.6351915039881883121e9,-3.9777071456931116607e9,
   -2.9984254115851086272e8,-7.4500417551436527905e8
    }};

    const std::array<double, 225>
    Implementation<BToXsGamma<NLO>>::f22int2coeffs
    {{
   0.0018509779924410270276,9.9144147820230405604,4.2814336341119616213,-21.440562899628940902,
   -98.448296061267621012,-46.679192145877183462,97.249667718770751943,-122.35940515532813674,
   -106.28380898090772938,62.628867732067378838,-93.958782870783115486,346.06681804737775645,
   -316.61893370972456689,103.24256523549438971,-26.162595749979753074,-14872.751073770159495,
   298898.88576444127001,-162830.1797330316567,697803.28815206981859,-792337.16525744979982,
   2.5083612245055595999e6,1.7074631228798882876e7,3.6939885862007903947e6,-1.6825986348151681184e7,
   1.3283356975168157435e7,2.4723252590962497977e8,-1.3531631723534654101e8,1.142491442925697538e7,
   6.2812543773879196617e7,-8.1366715058475577635e7,-5888.493190970520711,118302.25588539728057,
   -64485.782791303194102,276361.84351819033075,-313308.03676714714428,993313.75765731202316,
   6.759915035578605515e6,1.463001340465122145e6,-6.6614284732865889872e6,5.258957486805766306e6,
   9.7886025038904103474e7,-5.3576486087544264879e7,4.524712917299696732e6,2.4868671935665776927e7,
   -3.2215033729358305371e7,-136634.13608775393161,1.54493374785626687e6,-2.0553393370114135538e6,
   1.4057774583566200253e7,-1.2747494322253709123e7,1.9998814927171076001e7,6.740553742381621864e7,
   1.6562823154191508669e8,-3.6139421806241949411e7,1.5270307369496894965e8,3.6882820636352584234e9,
   1.3353708952075290174e8,1.7872645007775517178e8,1.3819835573808907646e7,-9.0033195656080527079e8,
   -55429.440738413906223,961175.73204288507918,-696440.66996352621528,2.183158963716768262e6,
   -4.0792380645614774492e6,8.9384912358411432743e6,4.3329267123166232462e7,5.3636420211996630036e7,
   -4.3584551392659267217e7,5.0866254187796084646e7,1.5865122212858145206e9,4.8846460553399857931e7,
   -1.8657636833695272438e8,-2.0672909866269747796e8,-3.9415868035850367206e8,11010.712553425845205,
   -94581.22864410086648,172308.89580786594657,-1.8782792684011649913e6,1.0084576926864659828e6,
   -1.5369381102462794913e6,-6.6682925363583181389e6,-7.5796867525414633697e6,
   1.4742870605813071696e6,-1.4109924340058980596e7,-1.4181574845209484494e8,1.0334673597116673396e8,
   -1.1071627664571654965e8,-1.327569046401633216e8,4.8202894449283333544e7,-5.7983717946399447864e6,
   1.6644967685148154983e8,-4.1807848097601375993e7,-1.8473212599218757355e8,
   -1.2326574350247716422e8,1.1037248698785939226e9,9.210226536504317261e9,-1.8360945945878849115e9,
   -1.0511544015178108303e10,3.8794853890750482178e9,9.7695698273343704766e10,
   -3.8881321341172616874e10,1.1776164942062392081e10,2.6450713821404160357e10,
   -3.0512083553937758294e10,-2.9956732902954109325e6,5.562211757931405535e7,
   -3.4808843272126387074e7,1.4890289555968915351e8,-1.1708929082220953103e8,5.2573292331471217079e8,
   3.3535294300904809148e9,8.4901075172523107867e8,-3.2872891791801963573e9,2.6600735625746217746e9,
   5.1498522376331526092e10,-2.4732672333541271845e10,3.9700033135784979181e9,
   1.2776789226966079331e10,-1.6405099470519833774e10,-266860.12041648895793,
   -1.5809587778204395304e6,-5.9144662071014241254e6,8.8373359859423361324e7,-6.083265367068793081e7,
   1.6268459772047133959e7,-1.0660362055396017633e8,5.7958428762567494077e8,3.558414525319128418e8,
   4.7452235544103758495e8,3.8089323278243512432e9,-4.7958983398664901751e9,-4.9850001818022340383e8,
   1.3508976462798921633e9,-1.5753063447526921937e9,-53534.914961478101852,1.5354867508193685835e6,
   -379634.19747902808656,7.8809187196750400926e6,-1.9051120650862417946e7,637285.51100015803488,
   3.0770524725143666973e7,4.655859762950385849e7,-3.5678820673028604372e6,8.0971362243109773451e7,
   8.1057768810814212557e8,-5.4200524573915444856e8,1.7815137338009908737e8,4.1302108275977111481e8,
   -2.8430228901040963516e8,6.4831241508396225232e6,-4.5901360012129333729e8,
   -7.1735199046907209674e7,6.8185782461092507654e8,3.094091715689826272e9,1.8395993335514974866e8,
   -1.2184145762422421279e10,4.4098211128841144415e9,1.3410293112682866967e10,
   -6.661452642625989244e9,-1.0848887688920711104e11,3.4308371911305491427e10,
   -2.9575710816393271385e9,-3.1834271811622338628e10,3.5419314652077025856e10,
   -3.9343258043637018866e6,-1.9357741113690389338e8,-1.6154285771448794828e8,
   1.7210438589870628469e8,3.6834011618592181278e9,2.5795680867898870758e9,5.4486725971004258826e9,
   -2.278563144629476472e8,-7.5081629887490499134e9,-1.0596420959502196473e9,
   7.0624223140657605664e10,-2.730777026808419679e10,2.351179123337899567e10,1.417300506684382849e10,
   -1.8759537550232138313e10,1.0935412014588823675e6,-1.9753830980366461888e8,
   -6.4273713668305627033e7,3.1842560480355491748e8,1.8311322414475543812e9,6.5962972439466899918e8,
   -2.8600429415616684316e9,1.7745766311956454083e9,2.9351833223233719199e9,-2.1763714777454040875e9,
   -1.7818718281766462376e10,1.2672421684277214441e9,2.8175901341137028807e9,
   -7.1925275094998248132e9,6.4315717796029008616e9,-173222.72377090499298,-1.4599372396528603285e7,
   -9.7494667564065501013e6,-5.061011101724842655e7,3.6397889601698491743e8,2.1717795109217579543e8,
   6.1849475993382562062e8,-4.9447148849018108669e8,-9.8725157327205511707e8,
   -4.2421183610845372823e8,3.0450575104722378106e9,1.662468029411882364e8,2.4323039971247396768e9,
   6.9111747428617053093e8,-5.4161273267310000476e8,16036.646856039014236,-2.9051336398569834504e6,
   -949857.53291523292127,-5.7215111383426782575e6,4.6980891149328622724e7,2.036441104268170717e7,
   1.9487579224223853122e7,-4.917217211741891997e7,-6.0483831549357655013e7,-8.1279075914450644395e7,
   -2.2554179279373643019e8,2.2196902934300758236e8,8.2093196840044559638e7,-1.9119090441256435969e8,
   1.1492404532552893702e8
    }};

    const std::array<double, 225>
    Implementation<BToXsGamma<NLO>>::f27int1coeffs
    {{
   1.430222424927796887,-51.4954828524410487,-6.6051576122386858943,92.029102045257214035,395.95220675529791832,176.05943473290642169,-449.37389451874300046,582.46336834368025868,488.52240195788447571,
   -247.58987071857043467,483.7487929775183233,-1604.5156987568700341,1378.1036357038927109,-444.33589176886901487,106.77137880685820824,-5.4496836865570700546e8,3.5843973296579379752e8,-3.6748493350256557396e8,
   2.2818648572314095729e8,-3.1418030822368831091e8,3.563649858536344023e8,5.4629695116489302326e8,-7.1985636436099184195e7,9.5613852153615390651e7,-5.143983171863565013e8,9.1153167592525453838e9,-1.6225564461624367565e9,
   -2.9460430577915542203e9,3.2713863377753683265e9,1.1045807924061549656e9,-2.1576663101167379172e8,1.4191525345236111659e8,-1.4549639238553146355e8,9.0344686250209813807e7,-1.2439262576535952805e8,
   1.4109334705880807943e8,2.162924482881903617e8,-2.8500248390056088241e7,3.785591997617721719e7,-2.0366271017547682956e8,3.6089853523424269764e9,-6.4241138705296413455e8,-1.1664144300741704518e9,1.2952252203171383145e9,
   4.3733131608265325377e8,-4.8541242587824425058e9,3.066680534563891325e9,-3.3290730092202566656e9,3.1771656425832593347e9,-3.5328657482698873691e9,2.6711623085999611308e9,-3.4722979231732059721e9,
   1.0597192997108954445e10,1.3977874502900696031e10,-5.486967387224972771e8,1.4133924761625979955e11,2.0148596236641237523e10,-4.4721212246195286228e10,-1.3460394716348490074e9,1.4590852416293162182e9,
   -1.8436052921375571923e9,1.1929431228646846626e9,-1.2520980562232906671e9,8.9994189470386124833e8,-1.2303162522101290708e9,1.0754359175613895664e9,1.1056826764975226815e8,2.6401136061018134217e9,2.8494718085954008621e9,
   -1.267381964960824163e9,6.3786963526174663512e10,1.7542879419302495749e10,-2.3868003374309426617e10,-8.8832917695761634299e9,8.3638557991321967221e7,4.5377537348601449168e8,-2.8647853876914483747e8,
   3.1121685698172931373e8,-3.1877632619847498394e8,3.1035982231841118595e8,-2.6943314179037971119e8,1.6433561795423081966e8,-5.6013518549140688862e8,-1.1372090579200898311e9,-1.6663927855666572298e7,
   -3.9727567859761694524e9,5.0462732062243386777e9,-1.5426344721180838845e8,-5.8439883953073931699e9,-1.052557387735055155e9,-2.1384400236105350694e11,1.4622304546490104432e11,-1.4178485624395038603e11,
   3.1666802645500136152e10,-9.2394750643049630803e10,1.5967788215659933022e11,5.3623499208888352631e11,-4.268071203007096588e11,-5.0685830171315787911e11,-3.9486645409313602877e11,3.3602467567581018309e12,
   3.2371381831981487327e11,-7.1510953809854281992e11,1.2960971279364751053e12,5.8565495215940601342e11,-1.1011334555094499606e11,7.2121420566741249323e10,-7.4385559395643716132e10,4.575087415477440797e10,
   -5.8843882736985581658e10,7.4319927843025854837e10,1.1134238533908952581e11,-1.5594845170216193418e10,1.3525382842514657711e10,-1.1036941213803154059e11,1.878951680762465987e12,-2.6434279712652462247e11,
   -5.5195480024867332348e11,6.4505128594991129803e11,2.2480995220761968779e11,-9.8437912829958546897e9,5.6457975396904004604e9,-6.9969964348546384292e9,1.3924216255458577995e10,-1.2323246487500222474e10,
   2.4207365285716158568e9,-4.4160491091741426471e10,6.4609883961996785346e10,9.4302270119274271484e10,2.5437654418451236458e10,1.7247390917165948085e11,-2.0634932396833594845e11,-1.145190132774056397e11,
   7.8955625076569419468e10,-3.9153624205108096459e9,-2.0707690077341991725e9,1.3652299298475323926e9,-1.3950170279897130415e9,1.815470522106046838e9,-2.9588891474505522383e9,4.2364776075345836424e8,
   -3.2619866311731810847e9,5.9148737975947798904e9,1.0240169610781283191e10,2.9363815256054401722e9,2.9852412271172234849e10,-2.3214330260293939835e10,-1.5463765023580566076e10,2.031767878235915792e10,
   1.6261834743659416043e9,2.388915473206755816e11,-1.751336146604579096e11,1.5327568772146934117e11,-2.1952434116622622256e10,2.4407856757147883614e11,-1.1000323604942852518e11,-6.3924192373297131646e11,
   5.2859767268346714957e11,5.6071195244948761332e11,3.0453237054911370572e11,-3.7525903839422447434e12,-6.5645443026366399285e11,1.3329342654069982773e12,-1.5533126716708941057e12,-5.7228842276474727023e11,
   -1.4519523906006568939e11,8.2544330427039035092e10,-1.0354072684859888513e11,1.7084343944930180821e10,1.8076298480207207434e11,2.2862575949665165815e11,4.3947193898487295057e11,-3.7625744751868545624e11,
   -5.9349848820389866464e11,-5.7838298430085791454e11,2.3572818352128149782e12,3.1522739627330770906e11,5.2408896321185503467e11,6.1007099876516952212e11,5.8307670996736115854e11,4.0733489890777116438e10,
   -4.0912806348807778682e10,2.1336891076447687629e10,1.51367375933492814e10,1.6128447209019433671e11,3.8665479169524832339e10,-1.8605716034979349141e11,1.8925417412992782956e11,1.5931390050237046867e11,
   -4.4817260287081262338e10,-5.4850943313898772515e11,-4.4514282727742672355e11,6.1462968373135208217e11,-3.9984736167079563326e11,-5.0366020614067492077e10,-6.5096847516374668077e9,2.5017903210421222586e9,
   -5.1631816292170266165e9,-4.6684128605023658095e9,3.5310586226140781798e10,2.4083549209593491636e10,5.5631938547760450607e10,-6.0875692309119669687e10,-9.9057001287126511934e10,-7.378789359525384534e10,
   8.1724486968979384088e10,9.2211234226097463099e10,1.6197053456871077345e11,7.343440746337602821e9,5.8981490439029974799e10,6.471140073742371957e8,-6.9363924465320331601e8,3.1999763605321126438e8,
   -1.1352495941691543364e9,5.730313404898461658e9,2.2766152338366495335e9,5.1194527879429432763e9,-6.6850615458026881829e9,-1.2533953603595751803e10,-8.4324624203022541744e9,-9.6209490427814769333e9,
   1.4973730053823531953e10,2.5353923033285711747e10,-1.2146468353158750669e10,4.2334099803814387221e9
    }};

    const std::array<double, 225>
    Implementation<BToXsGamma<NLO>>::f27int2coeffs
    {{
   0.19724041253467676017,26.198459707371890813,9.8019417964937837247,-51.679837587134881575,
   -253.06707379962533923,-116.12371096167627543,239.21090102266641812,-303.24060926420325944,
   -259.33974335852279305,157.40405382309917882,-264.75782989752801171,871.85731269167750417,
   -798.58534616490073209,258.68325216394388837,-66.014915792593664978,-6.1917218163246283936e6,
   -4.2929694701590590509e7,3.5535606608782140966e7,-9.2381207507200819339e7,1.1222214527604050626e8,
   -1.1822621342884956336e8,-2.5240601368806032604e8,1.3381775376700021561e8,
   -8.9330170784072537118e7,1.128161204322267384e8,-2.1424308251208304322e9,1.2060964611965915977e9,
   -3.4014129704295703855e8,-1.3722481477985528521e8,3.0150670321393881678e8,
   -2.4514564711933680357e6,-1.6997179632694559043e7,1.4069323461695152229e7,
   -3.6575471400643169263e7,4.4433959159091451463e7,-4.680755538116108987e7,-9.9936357261853877012e7,
   5.2984878721090173096e7,-3.5365204577346772089e7,4.466518389937368574e7,-8.4824027499463617202e8,
   4.7751516214761437166e8,-1.3466249374537895056e8,-5.4333425496584570115e7,1.1937480001341959755e8,
   -5.5134264900786578215e7,-3.7675959478956310814e8,3.1969714968103598288e8,
   -8.6152521237559809068e8,1.0310860133070221476e9,-1.0411385831208288699e9,-1.8040950644946015e9,
   4.5951790536426063758e8,-1.393995546318670738e9,8.3218869590096784281e8,-2.8054133815577932421e10,
   2.856053219903668584e9,-3.2637902419183495964e9,2.0037910154520784548e9,3.6693546759385478218e9,
   -2.0923981673153993716e7,-1.4502348257453273818e8,1.2104929622617410952e8,
   -3.1041305980617574522e8,3.8672807442464218747e8,-4.0351140078124817494e8,
   -7.8274892908404026468e8,2.110903681271626449e8,-3.7026850317368892604e8,3.583359934234485661e8,
   -1.1970025800945283632e10,1.0352127693985603579e9,2.9413991377494300645e8,1.9817287918747073236e9,
   1.6789905013156770204e9,5.1618994896083133172e6,3.4949425256545942696e7,-2.9755474489551375001e7,
   8.381216582643155492e7,-9.560889756530064515e7,9.515017559616495569e7,1.6874304936184728464e8,
   -9.2313114493756423842e7,1.411877112352425209e8,-7.6079674316650795761e7,1.3205085514653507592e9,
   -9.7341529389015251099e8,8.9185850432559823584e8,5.7700133544009918043e8,-1.4435903359303179487e8,
   -2.4297971575974470546e9,-1.7092227365240495727e10,1.3832652448362846526e10,
   -3.3885097852909322696e10,4.2978176094367117272e10,-4.7045882080760937092e10,
   -1.1191688149943774135e11,6.9519681990776000263e10,-1.4939804880908642205e10,
   5.1194532293252267374e10,-8.4380356681031148726e11,3.9680063692365410715e11,
   -1.7449334758976451306e11,-6.3603965967816862542e10,1.1097031867364374465e11,
   -1.2511126527219075266e9,-8.6514260451746181603e9,7.1890583474160837429e9,
   -1.8700355207169705473e10,2.2472906485029462851e10,-2.3975232082581279322e10,
   -5.0492462610431303006e10,2.6521387180804849291e10,-1.8632299665631940881e10,
   2.2900687425983686642e10,-4.4134620517460136546e11,2.2876545694504997224e11,
   -7.8053693972719048996e10,-2.8106369826948285656e10,6.0787676846578719829e10,
   -1.1184727238012708482e8,-7.3987194814434993413e8,6.5696370395677514023e8,-2.062634013192449442e9,
   2.2776136222866351943e9,-1.9797645851287205598e9,-2.4089985961826360808e9,
   -2.2181976967809536834e8,-5.0306328495118083876e9,8.0939914209340967091e8,
   -3.4511158082304077609e10,3.4668148270643495886e10,-2.4520534180815353121e9,
   -3.8174640343257854888e9,5.9560978948829502582e9,-2.3540876729132625043e7,
   -1.6465069741535153224e8,1.3397644951999073663e8,-3.8016992932589514795e8,5.0581063804537338748e8,
   -4.0429959657349322944e8,-7.6891011164893769683e8,3.3827812757099511319e8,
   -6.6711455702177919532e8,2.6742405721153644402e8,-7.1978693835602961649e9,4.7943389760244744097e9,
   -2.1785428611893750227e9,-1.5136200804306301296e9,1.0060689625955248815e9,2.7143760670091072727e9,
   2.0359653393295351425e10,-1.490423077674483645e10,3.5717577672192535085e10,
   -6.1853824130510497872e10,4.5917611618158873397e10,1.3334994296396093445e11,
   -8.7944864311593259948e10,9.6895514711665934086e9,-4.6080744345217152336e10,
   9.421408570428965835e11,-4.0044691028722482895e11,1.4658026755797355611e11,
   8.0831424706542867827e10,-1.3056786623214592707e11,-1.6497970647791604601e9,
   -1.0182785755289840803e10,1.0009917692800936188e10,-2.4179791595860242644e10,
   1.1301727593431009617e10,-4.0647361523975189806e10,-7.3455274883860045772e10,
   4.3981959983837768417e10,-6.4537996382977136743e9,5.2846242452520608244e10,
   -5.9534908114921549009e11,2.6585428565675946386e11,-1.9412783285404384913e11,
   -2.4934030483151305e10,6.5751613584633033899e10,4.6288608292702077614e8,4.0344447395706140115e9,
   -2.2954361477213641153e9,5.1647920230355754153e9,-1.676048537600857254e10,4.8289698676495714142e9,
   2.6212931944726142404e10,-1.9611577834075486027e10,-1.0790686454661505492e9,
   -2.7951506095730092501e9,1.5712331019916679721e11,-4.6080199406706822882e10,
   1.0309524363005100025e10,2.2538740952452883377e10,-2.4146190853055183094e10,
   -7.3982306089320407709e7,-4.3047854784136843982e8,4.5972621653633638807e8,
   -7.7197031905853483935e8,-4.7760365879522466347e8,-2.3260294030997696666e9,-5.27201129184152302e9,
   4.5472381825368103447e9,3.1189689557021085047e9,4.2844574285691341102e9,-2.5908327207890640011e10,
   4.3549184754389675653e9,-1.5939575584423593031e10,-1.5305620447894595306e9,
   1.4079628439808736233e9,7.3596930169176981228e6,6.2588679702319457907e7,-3.6941502461861927179e7,
   1.3928169854153586804e8,-3.574276102031520974e8,2.5907104459128270016e7,7.1717771496424844843e7,
   9.2820483693983334273e7,5.4757219049159672836e8,2.0002069938382543499e8,2.0534437914839038011e9,
   -1.8441760086436233282e9,4.3276782682320268935e7,7.9816336476471046314e8,-4.4480846275608887503e8
    }};
}
