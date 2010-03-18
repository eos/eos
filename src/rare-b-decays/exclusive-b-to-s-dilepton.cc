/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <src/rare-b-decays/form_factors.hh>
#include <src/utils/concrete_observable.hh>
#include <src/utils/integrate.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>
#include <tr1/functional>
#include <utility>
#include <map>
#include <vector>

#include <gsl/gsl_sf.h>

namespace wf
{
    // Large Recoil

    template <>
    struct Implementation<BToKstarDilepton<LargeRecoil>>
    {
        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter c7;

        Parameter c7prime;

        Parameter c8;

        Parameter c9;

        Parameter c9prime;

        Parameter c10;

        Parameter c10prime;

        Parameter m_b_MSbar;

        Parameter m_c;

        Parameter m_B;

        Parameter m_Kstar;

        double m_l;

        Parameter mu;

        double lambda_t2;

        Complex<double> lambda_u_hat;

        double f_B;

        double f_Kstar_par;

        double f_Kstar_perp;

        double lambda_B_p;

        Parameter a_1_par;

        Parameter a_2_par;

        Parameter a_1_perp;

        Parameter a_2_perp;

        Parameter ckm_A;

        Parameter ckm_lambda;

        double e_q;

        std::tr1::shared_ptr<FormFactors<BToKstar>> form_factors;

        Implementation(const Parameters & p, const ObservableOptions & o) :
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            c7(p["c7"]),
            c8(p["c8"]),
            c9(p["c9"]),
            c10(p["c10"]),
            c7prime(p["c7prime"]),
            c9prime(p["c9prime"]),
            c10prime(p["c10prime"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c(p["mass::c"]),
            m_B(p["mass::B0"]),
            m_Kstar(p["mass::K^*0"]),
            mu(p["mu"]),
            f_B(0.180), // +/- 0.03 GeV, cf. [BHP2008], Table 1, p. 8
            f_Kstar_par(0.225), // +/- 0.005 GeV, cf. [BHP2008], Table 1, p. 8
            f_Kstar_perp(0.185), // +/-0.005 GeV, cf. [BHP2008], Table 1, p. 8
            lambda_B_p(0.458),// +/- 0.115 GeV, cf. [BHP2008], Table 1, p. 8
            a_1_par(p["B->K^*::a_1_par"]),
            a_2_par(p["B->K^*::a_1_par"]),
            a_1_perp(p["B->K^*::a_1_perp"]),
            a_2_perp(p["B->K^*::a_2_perp"]),
            ckm_A(p["CKM::A"]),
            ckm_lambda(p["CKM::lambda"]),
            e_q(-1.0/3.0)
        {
            form_factors = FormFactorFactory<BToKstar>::create(o["form-factors"], p);
            if (! form_factors.get())
                throw std::string("InternalError");

            // TODO: Lepton masses, m_l = m_mu
            m_l = 0.10565836; // (GeV), cf. [PDG2008], p. 13
        }

        double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

        double norm(const double & s) const
        {
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]
            static const double g_fermi = 1.16637e-5; // (Gev^-2 (hbar c)^3), cf. [PDG2008], p.5

            double lambda_t = ckm_A * ckm_lambda * ckm_lambda;

            return std::sqrt(g_fermi * g_fermi * alpha_e * alpha_e / 3.0 / 1024 / std::pow(M_PI, 5.0) / m_B
                    * lambda_t * lambda_t * s_hat(s)
                    * std::sqrt(lambda(m_B * m_B, m_Kstar * m_Kstar, s))); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(double s) const
        {
            return s / m_B / m_B;
        }

        inline double energy(const double & s) const
        {
            return (m_B * m_B + m_Kstar * m_Kstar - s) / (2.0 * m_B);
        }

        inline double m_b_pole() const
        {
            return QCD::mb_pole(m_b_MSbar, mu);
        }

        /* Effective wilson coefficients */
        // cf. [BFS2001], below Eq. (9), p. 4
        double c7eff() const
        {
            return c7() - 1.0/3.0 * c3() - 4.0/9.0 * c4() - 20.0/3.0 * c5() - 80.0/9.0 * c6();
        }

        // cf. [BFS2001], below Eq. (26), p. 8
        double c8eff() const
        {
            return c8() + c3() - 1.0/6.0 * c4() + 20.0 * c5() - 10.0/3.0 * c6();
        }

        /* Charm loop parts */
        // cf. [BFS2001], Eq. (11), p. 4
        Complex<double> h(const double & s, const double & m_q) const
        {
            const double z = 4.0 * m_q * m_q / s;
            const double sqrt1z = std::sqrt(std::abs(z - 1.0));

            double a = 2.0 * std::log(m_q / mu()) - 2.0 / 3.0 - z;
            double b = (2.0 + z) * sqrt1z;
            double rc, ic;
            if (z > 1.0)
            {
                ic = 0.0;
                rc = std::atan(1.0 / sqrt1z);
            }
            else
            {
                ic = -M_PI / 2.0;
                rc = std::log((1.0 + sqrt1z) / std::sqrt(z));
            }

            return -4.0 / 9.0 * (a + b * Complex<double>::Cartesian(rc, ic));
        }

        // cf. [BFS2001], Eq. (11), p. 4 in the limit m_q -> 0
        Complex<double> h0(const double & s) const
        {
            return 4.0 / 9.0 * Complex<double>::Cartesian(2.0 / 3.0 + 2.0 * std::log(2.0 * mu()) - std::log(s), M_PI);
        }

        // cf. [BFS2001], Eq. (10), p. 4
        Complex<double> Y0(const double & s) const
        {
            double Y_c = 4.0 / 3.0 * c1() + c2() + 6.0 * c3() + 60.0 * c5();
            double Y_b = -0.5 * (7.0 * c3() + 4.0 / 3.0 * c4() + 76.0 * c5() + 64.0 / 3.0 * c6());
            double Y_0 = -0.5 * (c3() + 4.0 / 3.0 * c4() + 16.0 * c5() + 64 / 3.0 * c6());
            double Y = 2.0 / 9.0 * (6.0 * c3() + 32.0 * c5() + 32.0 / 3.0 * c6());

            // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
            return Y_c * h(s, m_c) + Y_b * h(s, m_b_pole()) + Y_0 * h0(s) + Y;
        }

        // cf. [BFS2004], ?
        Complex<double> Y0u(const double & s) const
        {
            double a = 4.0 / 3.0 * c1() + c2();

            return a * (h(s, m_c) - h0(s));
        }

        /* Form factors */
        //  cf. [BHP2008], Eq. (E.4), p. 23
        double xi_perp(const double & s) const
        {
            const double factor = m_B / (m_B + m_Kstar);

            return factor * form_factors->v(s_hat(s));
        }

        double xi_par(const double & s) const
        {
            const double factor1 = (m_B + m_Kstar) / (2.0 * energy(s));
            const double factor2 = (1.0 - m_Kstar / m_B);

            return factor1 * form_factors->a_1(s_hat(s)) - factor2 * form_factors->a_2(s_hat(s));
        }

        /* NLO functions */
        // cf. [BFS2001], Eq. (29), p. 8
        Complex<double> B0(const double & s, const double & m_q) const
        {
            double z = 4.0 * m_q * m_q / s;
            double rp, ip;

            if (m_q < 0.0)
            {
                throw InternalError("Implementation<BToKstarDilepton<LargeRecoil>>::B0: m_q < 0!");
            };

            if ((0.0 == m_q) && (0.0 == s))
            {
                throw InternalError("Implementation<BToKstarDilepton<LargeRecoil>>::B0: m_q == 0 & s == 0");
            }

            if (0 == s)
                return Complex<double>::Cartesian(-2.0, 0.0);

            if (z > 1.0)
            {
                rp = -2.0 * std::sqrt(z - 1.0) * std::atan(1.0 / std::sqrt(z - 1.0));
                ip = 0.0;
            }
            else
            {
                rp = std::sqrt(1.0 - z) * std::log((1.0 - std::sqrt(1 - z)) / (1.0 + std::sqrt(1.0 - z)));
                ip = std::sqrt(1.0 - z) * M_PI;
            }

            return Complex<double>::Cartesian(rp, ip);
        }

        // cf. [BFS2001], Eq. (84), p. 30
        double C0(const double & s) const
        {
            static const int points = 40;
            // Integration boundaries of u = (0, u_max]
            static const double dx = (1.0 - 0.0) / points;
            static const double g[2] =
            {
                (1.0 + std::sqrt(3.0 / 5.0)) / 2.0,
                (1.0 - std::sqrt(3.0 / 5.0)) / 2.0
            };

            long double s_hat = s / m_B / m_B;
            long double x1, x2, x;
            long double y[3];
            long double result = 0;

            for (int i = 0; i < points; i++)
            {
                // 0 is lower integration boundary

                x1 = 0.0 + i * dx;
                x2 = 0.0 + (i + 1) * dx;
                for (int j = 0; j < 3; j++)
                {
                    switch (j)
                    {
                        case 0: x = x1 * g[0] + x2 * g[1];
                                break;

                        case 1: x = (x1 + x2) / 2.0;
                                break;

                        case 2: x = x1 * g[1] + x2 * g[0];
                                break;
                    }

                    y[j] = std::log(x * x / (1.0 - s_hat * x * (1.0 - x))) / (1.0 + x * (1.0 - s_hat));
                }

                result += (5.0 * y[0] + 8.0 * y[1] + 5.0 * y[2]) * dx / 18.0;
            };

            return result;
        }

        // cf. [BFS2001], Eqs. (30)-(32), p. 8
        Complex<double> I1(const double & s, const double & u, const double & m_q) const
        {
            if (m_q == 0.0)
                return Complex<double>::Cartesian(1.0, 0.0);

            int status;
            gsl_sf_result res_re, res_im;

            double ubar = 1.0 - u;
            double m_q2 = m_q * m_q;
            double m_B2 = m_B * m_B;

            double a, a2, sign;
            Complex<double> dilogArg, dilog1, dilog2;
            Complex<double> LxpLxm, LypLym;

            if (1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s)) > 0)
            {
                a = (1 - std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s))))
                  / (1 + std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s))));
                LxpLxm = -M_PI* M_PI / 3.0
                    + std::log(a) * (std::log(a) + Complex<double>::Cartesian(0.0, M_PI))
                    + gsl_sf_dilog(-a) + gsl_sf_dilog(-1.0/a);
            }
            else
            {
                a  = sqrt(4.0 * m_q2 / (m_B2 - u * (m_B2 - s)) - 1);
                a2 = a * a;

                if (a2 - 1 > 0)
                    sign = +1.0;
                else
                    sign = -1.0;

                dilogArg = Complex<double>::Cartesian((a2 - 1.0)/(a2 + 1.0), -2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(dilogArg.absolute(), dilogArg.phase(), &res_re, &res_im);
                dilog1 = Complex<double>::Cartesian( res_re.val, res_im.val);

                dilogArg = Complex<double>::Cartesian((a2 - 1.0)/(a2 + 1.0), +2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(dilogArg.absolute(), dilogArg.phase(), &res_re, &res_im);
                dilog2 = Complex<double>::Cartesian(res_re.val, res_im.val);

                LxpLxm = -1.0 / 3.0 * M_PI * M_PI - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                    + dilog1 + dilog2;
            }

            if (1.0 - 4.0 * m_q2 / s > 0)
            {
                a = (1.0 - std::sqrt(1.0 - 4.0 * m_q2 / s))
                  / (1.0 + std::sqrt(1.0 - 4.0 * m_q2 / s));
                LypLym = -1.0 / 3.0 * M_PI * M_PI + std::log(a) * (std::log(a) + Complex<double>::Cartesian(0.0, M_PI))
                    + gsl_sf_dilog(-a) + gsl_sf_dilog(-1./a);
            }
            else
            {
                a  = std::sqrt(4.0 * m_q2 / s - 1.0);
                a2 = a * a;
                if (a2 - 1.0 > 0)
                    sign = +1.0;
                else
                    sign = -1.0;

                dilogArg = Complex<double>::Cartesian((a2 - 1.0) / (a2 + 1.0), -2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(dilogArg.absolute(), dilogArg.phase(), &res_re, &res_im);
                dilog1 = Complex<double>::Cartesian(res_re.val, res_im.val);

                dilogArg = Complex<double>::Cartesian((a2 - 1.0) / (a2 + 1.0), +2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(dilogArg.absolute(), dilogArg.phase(), &res_re, &res_im);
                dilog2 = Complex<double>::Cartesian(res_re.val, res_im.val);

                LypLym = -1.0 / 3.0 * M_PI * M_PI - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                    + dilog1 + dilog2;
            }

            return 1.0 + 2.0 * m_q2 / (ubar * (m_B2 - s)) * (LxpLxm - LypLym);
        }

        // cf. [BFS2001], Eq. (48)
        double phi_K(const double & u, const double & a_1, const double & a_2) const
        {
            double xi = 2.0 * u - 1.0;

            return 6.0 * u * (1 - u) * (1.0 + a_1 * 3.0 * xi + a_2 * (7.5 * xi * xi - 1.5));
        }

        // cf. [BFS2001], Eq. (18), p. 6, modulo the inverse of lambda_B,-
        double T0_perp_m(const double & s, const double & u) const
        {
            static const double e_q = -1.0/3.0;

            double wilson = c3 + 4.0/3.0 * c4 + 16.0 * c5 + 64.0/3.0 * c6;
            return -1.0 * e_q * (4.0 * m_B / m_b_pole()) * wilson;
        }

        // cf. [BFS2001], Eq. (20), p. 6, module the inverse of lambda_B,+
        double Tf_perp_p(const double & s, const double & u) const
        {
            return c7eff() * (2.0 * m_B / (1.0 - u) / energy(s));
        }

        // cf. [BFS2001], Eq. (21), p. 6
        Complex<double> Tf_par_p(const double & s, const double & u) const
        {
            double E = energy(s);

            return (c7eff() + (s / (2.0 * m_B * m_b_pole())) * Y0(s)) * (2.0 * m_B * m_B / (1 - u) / E / E);
        }

        // cf. [BFS2001], Eq. (27), p. 8
        Complex<double> t_perp(const double & s, const double & u, const double & m_q) const
        {
            double ubar = 1.0 - u;
            double E = energy(s);
            double x = ubar * m_B * m_B + u * s;

            Complex<double> result = (2.0 * m_B / ubar / E) * I1(s, u, m_q);
            if (m_q > 0.0)
                result = result + (s / ubar / ubar / E / E) * (B0(x, m_q) - B0(s, m_q));

            return result;
        }

        // cf. [BFS2001], Eq. (28), p. 8
        Complex<double> t_par(const double & s, const double & u, const double & m_q) const
        {
            double ubar = 1.0 - u;
            double E = energy(s);
            double x = ubar * m_B * m_B + u * s;

            Complex<double> result = (2.0 * m_B / ubar / E) * I1(s, u, m_q);
            if (m_q > 0.0)
                result = result + (x / ubar / ubar / E / E) * (B0(x, m_q) - B0(s, m_q));

            return result;
        }

        // cf. [BFS2001], Eq. (23), p. 7, multiplied by phi_K^*,perp
        Complex<double> Tnf_perp_p(const double & s, const double & u) const
        {
            double s_hat = s / m_B / m_B;
            double ubar = 1.0 - u;

            double a = (4.0 / 3.0 / (u + ubar * s_hat)) * c8eff();
            Complex<double> ba = (+2.0 / 3.0) * (-c1 / 6.0 + c2 + 6.0 * c6) * t_perp(s, u, m_c);
            Complex<double> bb = (-1.0 / 3.0)
                * (c3 - c4 / 6.0 + 16.0 * c5 + 10.0/3.0 * c6 - (4.0 * m_b_pole() / m_B) * (c3 - c4/6.0 + 4.0 * c5 - 2.0/3.0 * c6))
                * t_perp(s, u, m_b_pole());
            Complex<double> bc = (-1.0 / 3.0) * (c3 - c4 / 6.0 + 16.0 * c5 - 8.0/3.0 * c6) * t_perp(s, u, 0.0);

            return (a + (m_B / 2.0 / m_b_pole()) * (ba + bb + bc)) * phi_K(u, a_1_perp, a_2_perp);
        }

        // cf. [BFS2001, Eq. (25), p. 7, multiplied by phi_K^*,par
        Complex<double> Tnf_par_p(const double & s, const double & u) const
        {
            Complex<double> a = (+2.0 / 3.0) * (-c1 / 6.0 + c2 + 6.0 * c6) * t_par(s, u, m_c);
            Complex<double> b = (-1.0 / 3.0) * (c3 - c1 / 6.0 + 16.0 * c5 + 10.0/3.0 * c6) * t_par(s, u, m_b_pole());
            Complex<double> c = (-1.0 / 3.0) * (c3 - c4 / 6.0 + 16.0 * c5 - 8.0/3.0 * c6) * t_par(s, u, 0.0);

            return (m_B / m_b_pole()) * (a + b + c) * phi_K(u, a_1_par, a_2_par);
        }

        // cf. [BFS2001], Eq. (26), pp. 7-8, multiplied by phi_K^*,par
        Complex<double> Tnf_par_m(const double & s, const double & u) const
        {
            double s_hat = s / m_B / m_B;
            double ubar = 1.0 - u;
            double x = ubar * m_B * m_B + u * s;

            double a = (e_q * 8.0 / (ubar + u * s_hat)) * c8eff();
            Complex<double> ba = (-c1 / 6.0 + c2 + c4 + 10 * c6) * h(x, m_c);
            Complex<double> bb = (c3 + 5.0/6.0 * c4 + 16.0 * c5 + 22.0/3.0 * c6) * h(x, m_b_pole());
            Complex<double> bc = (c3 + 17.0/6.0 * c4 + 16.0 * c5 + 82.0/3.0 * c6) * h0(x);
            double bd = -8.0 / 27.0 * (-7.5 * c4 + 12.0 * c5 - 32.0 * c6);

            return (6.0 * m_B / m_b_pole()) * (ba + bb + bc + bd) * phi_K(u, a_1_par, a_2_par);
        }

        // cf. [BFS2001], Eq. (36), p. 9
        double L(const double & s) const
        {
            double m_b_pole2 = m_b_pole() * m_b_pole();

            return -1.0 * (m_b_pole2 - s) / s * std::log(1.0 - s / m_b_pole2);
        }

        // cf. [BFS2001], Eq. (54), p. 15
        Complex<double> lambda_B_m_inv(const double & s) const
        {
            if (0.0 == s)
                return Complex<double>::Cartesian(0.0, 0.0);

            static const double Lambda = 0.5; // (GeV), cf. [BFS2001], below Eq. (54), p. 15
            double omega_0 = 2.0 * Lambda / 3.0;
            double x = s / m_B / omega_0;
            double ei = gsl_sf_expint_Ei(x);

            return Complex<double>::Cartesian(-ei, M_PI) * (std::exp(-x) / omega_0);
        }

        // cf. [AAGW2001], Eq. (56), p. 20
        Complex<double> F17(const double & s) const
        {
            // Deviating definition of s_hat == s / m_b_pole^2 !
            double s_hat = s / m_b_pole() / m_b_pole();
            double Lmu = std::log(mu / m_b_pole());
            double Ls = std::log(s_hat);

            // cf. [AAGW2001], Table I, p. 20
            std::vector<std::pair<Complex<double>, Complex<double>>> f17_coefficients =
            {
                std::make_pair(Complex<double>::Cartesian(-0.76730, -0.11418), Complex<double>::Cartesian(-0.00000, -0.00000)),
                std::make_pair(Complex<double>::Cartesian(-0.28480, -0.18278), Complex<double>::Cartesian(-0.00328, +0.02083)),
                std::make_pair(Complex<double>::Cartesian(+0.05611, -0.23357), Complex<double>::Cartesian(+0.01637, +0.02091)),
                std::make_pair(Complex<double>::Cartesian(+0.62438, -0.02744), Complex<double>::Cartesian(+0.07673, +0.00914))
            };

            Complex<double> f17 = Complex<double>::Cartesian(0.0, 0.0);
            double i(0);
            for (auto k = f17_coefficients.begin() ; k != f17_coefficients.end() ; ++k, ++i)
            {
                f17 = f17 + std::pow(s_hat, i) * (k->first + k->second * Ls);
            }

            return -208.0/243.0 * Lmu + f17;
        }

        // cf. [AAGW2001], Eq. (54), p. 19
        Complex<double> F19(const double & s) const
        {
            // Deviating definition of s_hat == s / m_b_pole^2 !
            double s_hat = s / m_b_pole() / m_b_pole();
            double Lmu = std::log(mu / m_b_pole());
            double Ls = std::log(s_hat);
            double Lc = std::log(m_c / m_b_pole());
            double z = (m_c / m_b_pole()) * (m_c / m_b_pole());

            // cf. [AAGW2001], Table I, p. 20
            std::vector<std::pair<Complex<double>, Complex<double>>> f19_coefficients =
            {
                std::make_pair(Complex<double>::Cartesian(-12.7150, +0.09470), Complex<double>::Cartesian(-0.07883, -0.07414)),
                std::make_pair(Complex<double>::Cartesian(-38.7420, -0.67862), Complex<double>::Cartesian(-0.03930, -0.00017)),
                std::make_pair(Complex<double>::Cartesian(-103.830, -2.53880), Complex<double>::Cartesian(-0.04470, +0.00263)),
                std::make_pair(Complex<double>::Cartesian(-313.750, -8.45540), Complex<double>::Cartesian(-0.05113, +0.02275))
            };

            Complex<double> f19 = Complex<double>::Cartesian(0.0, 0.0);
            double i(0);
            for (auto k = f19_coefficients.begin() ; k != f19_coefficients.end() ; ++k, ++i)
            {
                f19 = f19 + std::pow(s_hat, i) * (k->first + k->second * Ls);
            }

            return Lmu * Complex<double>::Cartesian(-1424.0/729.0 + 64.0/27.0 * Lc - 16.0/243.0 * Ls + (16.0/1215.0 - 32.0/(z * 135.0)) * s_hat
                        + (4.0/2835.0 - 8.0/(315.0 * z * z)) * s_hat * s_hat + (16.0/76545.0 - 32.0/(8505.0 * z * z * z)) * s_hat * s_hat * s_hat
                        - 256.0/243.0 * Lmu,
                        16.0/243.0 * M_PI)
                + f19;
        }

        // cf. [AAGW2001], Eq. (56), p. 20
        Complex<double> F27(const double & s) const
        {
            // Deviating definition of s_hat == s / m_b_pole^2 !
            double s_hat = s / m_b_pole() / m_b_pole();
            double Lmu = std::log(mu / m_b_pole());
            double Ls = std::log(s_hat);

            // cf. [AAGW2001], Table II, p. 21
            std::vector<std::pair<Complex<double>, Complex<double>>> f27_coefficients =
            {
                std::make_pair(Complex<double>::Cartesian(+4.60380, +0.68510), Complex<double>::Cartesian(-0.00000, -0.00000)),
                std::make_pair(Complex<double>::Cartesian(+1.70880, +1.09670), Complex<double>::Cartesian(+0.01959, -0.12496)),
                std::make_pair(Complex<double>::Cartesian(-0.33665, +1.40140), Complex<double>::Cartesian(-0.09822, -0.12548)),
                std::make_pair(Complex<double>::Cartesian(-3.74630, +0.16463), Complex<double>::Cartesian(-0.18321, -0.05485))
            };

            Complex<double> f27 = Complex<double>::Cartesian(0.0, 0.0);
            double i(0);
            for (auto k = f27_coefficients.begin() ; k != f27_coefficients.end() ; ++k, ++i)
            {
                f27 = f27 + std::pow(s_hat, i) * (k->first + k->second * Ls);
            }

            return 416.0/81.0 * Lmu + f27;
        }

        // cf. [AAGW2001], Eq. (55), p. 19
        Complex<double> F29(const double & s) const
        {
            // Deviating definition of s_hat == s / m_b_pole^2 !
            double s_hat = s / m_b_pole() / m_b_pole();
            double Lmu = std::log(mu / m_b_pole());
            double Ls = std::log(s_hat);
            double Lc = std::log(m_c / m_b_pole());
            double z = (m_c / m_b_pole()) * (m_c / m_b_pole());

            // cf. [AAGW2001], Table II, p. 21
            std::vector<std::pair<Complex<double>, Complex<double>>> f29_coefficients =
            {
                std::make_pair(Complex<double>::Cartesian(+9.50420, -0.56819), Complex<double>::Cartesian(+0.47298, +0.44483)),
                std::make_pair(Complex<double>::Cartesian(+7.42380, +4.07170), Complex<double>::Cartesian(+0.23581, +0.00104)),
                std::make_pair(Complex<double>::Cartesian(+0.33806, +15.2330), Complex<double>::Cartesian(+0.26821, -0.01577)),
                std::make_pair(Complex<double>::Cartesian(-42.0850, +50.7320), Complex<double>::Cartesian(+0.30680, -0.13652))
            };

            Complex<double> f29 = Complex<double>::Cartesian(0.0, 0.0);
            double i(0);
            for (auto k = f29_coefficients.begin() ; k != f29_coefficients.end() ; ++k, ++i)
            {
                f29 = f29 + std::pow(s_hat, i) * (k->first + k->second * Ls);
            }

            return Lmu * Complex<double>::Cartesian(256.0/243.0 - 128.0/9.0 * Lc + 32.0/81.0 * Ls - (32.0/405.0 - 64.0/(z * 45.0)) * s_hat
                        - (8.0/945.0 - 16.0/(105.0 * z * z)) * s_hat * s_hat - (32.0/25515.0 - 64.0/(2835.0 * z * z * z)) * s_hat * s_hat * s_hat
                        + 512.0/81.0 * Lmu,
                        -32.0/81.0 * M_PI)
                + f29;
        }

        // cf. [BFS2001], Eq. (82), p. 30
        Complex<double> F87(const double & s) const
        {
            double s_hat = s / m_B / m_B;
            double s_hat2 = s_hat * s_hat;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            Complex<double> a = Complex<double>::Cartesian(-32.0 * std::log(mu / m_b_pole()) - 8.0 * s_hat / denom * std::log(s_hat)
                    - 4.0 * (11.0 - 16.0 * s_hat + 8.0 * s_hat2) / denom2,
                    -8.0 * M_PI);
            Complex<double> b = (4.0 / denom / denom2)
                * ((9.0 * s_hat - 5.0 * s_hat2 + 2.0 * s_hat * s_hat2) * B0(s, m_b_pole()) - (4.0 + 2.0 * s_hat) * C0(s));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [BFS2001], Eq. (83), p. 30
        Complex<double> F89(const double & s) const
        {
            double s_hat = s / m_B / m_B;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            double a = 16.0 * std::log(s_hat) / denom + 8.0 * (5.0 - 2.0 * s_hat) / denom2;
            Complex<double> b = (-8.0 * (4.0 - s_hat) / denom / denom2) * ((1.0 + s_hat) * B0(s, m_b_pole()) - 2.0 * C0(s));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [BFS2001], Eq. (15) with a = perp, and [BHP2008], Eq. (C.4)
        Complex<double> tensor_perp(const double & h, const double & s) const
        {
            double nlo_factor = QCD::alpha_s(mu) * QCD::casimir_f / 4.0 / M_PI;
            double scatt_factor = M_PI * M_PI / 3.0 * f_B * f_Kstar_perp / m_B;

            // Here m_b_pole is used, cf. [BFS2001], comment below Eq. (36), p. 9

            /* Form factor corrections */
            Complex<double> ff_0 = (c7eff() + h * c7prime()) + s / (2.0 * m_b_pole() * m_B) * Y0(s);
            // cf. [BFS2001], Eq. (34), p. 9
            double ff_f = (c7eff() + h * c7prime()) * (8.0 * std::log(m_b_pole() / mu) - 4.0 - L(s));
            // cf. [BFS2001], Eq. (37), p. 9
            Complex<double> ff_nf = (-1.0 / QCD::casimir_f) * (
                    (c2 - c1 / 6.0) * F27(s) + c8eff() * F87(s)
                    + (s / (2.0 * m_b_pole() * m_B)) * (c2() * F29(s) + c1() * F19(s) + c8eff() * F89(s)));
            Complex<double> ff = ff_0 + nlo_factor * (ff_f + ff_nf);

            /* Specator scattering, folded with phi_K^*,perp */
            // cf. [BFS2001], Eq. (20), p. 6
            double scatt_f_p = c7eff() * 2.0 * m_B / energy(s) * 3.0 * (1.0 + a_1_perp + a_2_perp);
            // cf. [BFS2001], Eq. (23), p. 7
            Complex<double> scatt_nf_p = integrate(
                    std::tr1::function<Complex<double> (const double &)>(
                        std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Tnf_perp_p, this, s, std::tr1::placeholders::_1)),
                    40, 0.01, 0.99);
            Complex<double> scatt_p = (1.0 / lambda_B_p) * nlo_factor * (scatt_f_p + scatt_nf_p);

            return xi_perp(s) * ff + scatt_factor * scatt_p;
        }

        // cf. [BFS2001], Eq. (15) with a = par, and [BHP2008], Eq. (C.4)
        Complex<double> tensor_par(const double & s) const
        {
            double nlo_factor = QCD::alpha_s(mu) * QCD::casimir_f / 4.0 / M_PI;
            double scatt_factor = M_PI * M_PI / 3.0 * f_B * f_Kstar_perp / m_B * m_Kstar / energy(s);

            // Here m_b_pole is used, cf. [BFS2001], comment below Eq. (36), p. 9

            /* Form factor corrections */
            Complex<double> ff_0 = -1.0 * (c7eff() - c7prime() + m_B / (2.0 * m_b_pole()) * (Y0(s)));
            // cf. [BFS2001], Eq. (35), p. 9
            Complex<double> ff_f = -1.0 * (c7eff() - c7prime()) * (8.0 * std::log(m_b_pole() / mu) - 6.0 + 4.0 * L(s))
                    + (m_B / (2.0 * m_b_pole())) * Y0(s) * (2.0 - 2.0 * L(s));
            // cf. [BFS2001], Eq. (38), p. 9
            Complex<double> ff_nf = (+1.0 / QCD::casimir_f) * (
                    (c2 - c1 / 6.0) * F27(s) + c8eff() * F87(s)
                    + (m_B / (2.0 * m_b_pole())) * (c2() * F29(s) + c1() * F19(s) + c8eff() * F89(s)));
            Complex<double> ff = ff_0 + nlo_factor * (ff_f + ff_nf);

            /* Spectator scattering */
            // cf. [BFS2001], Eq. (18), p. 6
            double scatt_0_m = -e_q * 4.0 * m_B / m_b_pole() * (c3 + 4.0/3.0 * c4 + 16.0 * c5 + 64.0/3.0 * c6);
            // cf. [BFS2001], Eq. (21), p. 6
            Complex<double> scatt_f_p = (c7eff() + (s / (2.0 * m_b_pole() * m_B)) * Y0(s)) * (2.0 * m_B * m_B / energy(s) / energy(s))
                    * 3.0 * (1.0 + a_1_par + a_2_par);
            // cf. [BFS2001], Eq. (25), p. 7
            Complex<double> scatt_nf_p = integrate(
                    std::tr1::function<Complex<double> (const double &)>(std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Tnf_par_p,
                            this, s, std::tr1::placeholders::_1)),
                    40, 0.01, 0.99);
            // cf. [BFS2001], Eq. (26), pp. 7-8
            Complex<double> scatt_nf_m = integrate(
                    std::tr1::function<Complex<double> (const double &)>(std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Tnf_par_m,
                            this, s, std::tr1::placeholders::_1)),
                    40, 0.01, 0.99);
            Complex<double> scatt_p = (1.0 / lambda_B_p) * nlo_factor * (scatt_f_p + scatt_nf_p);
            Complex<double> scatt_m = lambda_B_m_inv(s)  * (scatt_0_m + nlo_factor * scatt_nf_m);

            return xi_par(s) * ff + scatt_factor * (scatt_p + scatt_m);
        }

        /* Amplitudes */
        // cf. [BHP2008], p. 20
        Complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_pole() / m_B;

            double prefactor = -1.0 * m_B * m_B * (1.0 - shat) * (1.0 - shat) / (2.0 * m_Kstar * std::sqrt(shat));
            double wilson = (c9() - c9prime()) + h * (c10() - c10prime());

            return prefactor * (wilson * xi_par(s) - 2.0 * mbhat * tensor_par(s));
        }

        // cf. [BHP2008], p. 20
        Complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_pole() / m_B;

            double prefactor = +std::sqrt(2.0) * m_B * (1.0 - shat);
            double wilson = (c9() + c9prime()) + h * (c10() + c10prime());

            return prefactor * (wilson * xi_perp(s) + (2.0 * mbhat / shat) * tensor_perp(+1.0, s));
        }

        // cf. [BHP2008], p. 20
        Complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_pole() / m_B;

            double prefactor = -std::sqrt(2.0) * m_B * (1.0 - shat);
            double wilson = (c9() - c9prime()) + h * (c10() - c10prime());

            return prefactor * (wilson * xi_perp(s) + (2.0 * mbhat / shat) * tensor_perp(-1.0, s));
        }
    };

    BToKstarDilepton<LargeRecoil>::BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options) :
        PrivateImplementationPattern<BToKstarDilepton<LargeRecoil>>(new Implementation<BToKstarDilepton<LargeRecoil>>(parameters, options))
    {
    }

    BToKstarDilepton<LargeRecoil>::~BToKstarDilepton()
    {
    }

    Complex<double>
    BToKstarDilepton<LargeRecoil>::a_long(const Helicity & h, const double & s) const
    {
        return _imp->a_long(h, s);
    }

    Complex<double>
    BToKstarDilepton<LargeRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        return _imp->a_perp(h, s);
    }

    Complex<double>
    BToKstarDilepton<LargeRecoil>::a_par(const Helicity & h, const double & s) const
    {
        return _imp->a_par(h, s);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_branching_ratio(const double & s) const
    {
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        return differential_decay_width(s) / Gamma;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_decay_width(const double & s) const
    {
        return _imp->norm(s) * _imp->norm(s) * (a_long(left_handed, s).absolute_squared()
            + a_long(right_handed, s).absolute_squared()
            + a_perp(left_handed, s).absolute_squared()
            + a_perp(right_handed, s).absolute_squared()
            + a_par(left_handed, s).absolute_squared()
            + a_par(right_handed, s).absolute_squared());
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 * _imp->norm(s) * _imp->norm(s) / differential_decay_width(s) * (
                (a_par(left_handed, s) * a_perp(left_handed, s).conjugate()).real()
                -(a_par(right_handed, s) * a_perp(right_handed, s).conjugate()).real()
            );
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        double a = a_perp(left_handed, s).absolute_squared() + a_perp(right_handed, s).absolute_squared();
        double b = a_par(left_handed, s).absolute_squared() + a_par(right_handed, s).absolute_squared();

        return (a - b) / (a + b);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        double a = (a_long(left_handed, s) * a_par(left_handed, s).conjugate() + a_long(right_handed, s).conjugate() * a_par(right_handed, s)).absolute();
        double b = std::sqrt(a_long(left_handed, s).absolute_squared() + a_long(right_handed, s).absolute_squared())
                * std::sqrt(a_perp(left_handed, s).absolute_squared() + a_perp(right_handed, s).absolute_squared());

        return a / b;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        double a = (a_long(left_handed, s) * a_perp(left_handed, s).conjugate() - a_long(right_handed, s).conjugate() * a_perp(right_handed, s)).absolute();
        double b = (a_long(left_handed, s).conjugate() * a_par(left_handed, s) + a_long(right_handed, s) * a_par(right_handed, s).conjugate()).absolute();

        return a / b;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return (a_long(left_handed, s).absolute_squared() + a_long(right_handed, s).absolute_squared())
            * _imp->norm(s) * _imp->norm(s) / differential_decay_width(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::tr1::function<double (const double &)> f = std::tr1::bind(std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_branching_ratio),
                this, std::tr1::placeholders::_1);

        return integrate(f, 100, s_min, s_max);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        std::tr1::function<double (const double &)> f = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry), this, std::tr1::placeholders::_1);

        return integrate(f, 100, s_min, s_max);
    }


    // Low Recoil

    template <>
    struct Implementation<BToKstarDilepton<LowRecoil>>
    {
        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter c7;

        Parameter c7prime;

        Parameter c9;

        Parameter c9prime;

        Parameter c10;

        Parameter c10prime;

        Parameter m_b_MSbar;

        Parameter m_c;

        Parameter m_B;

        Parameter m_Kstar;

        Parameter mu;

        Parameter ckm_A;

        Parameter ckm_lambda;

        std::tr1::shared_ptr<FormFactors<BToKstar>> form_factors;

        Implementation(const Parameters & p, const ObservableOptions & o) :
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            c7(p["c7"]),
            c9(p["c9"]),
            c10(p["c10"]),
            c7prime(p["c7prime"]),
            c9prime(p["c9prime"]),
            c10prime(p["c10prime"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c(p["mass::c"]),
            m_B(p["mass::B0"]),
            m_Kstar(p["mass::K^*0"]),
            mu(p["mu"]),
            ckm_A(p["CKM::A"]),
            ckm_lambda(p["CKM::lambda"])
        {
            form_factors = FormFactorFactory<BToKstar>::create(o["form-factors"], p);
            if (! form_factors.get())
                throw std::string("InternalError");
        }

        // Helpers
        Complex<double> c7eff(double s) const
        {
            // TODO: Neglecting contributions ~alpha_s / 4.0 / M_PI. These do involve spectator scattering,
            // cf. [BFS2001] Eq. (29), p. 8, and Eqs. (82)-(84), p. 30
            double c7eff0 = c7 - 4.0/9.0 * c3 - 4.0/3.0 * c4 + 1.0/9.0 * c5 + 1.0/3.0 * c6;
            return Complex<double>::Cartesian(c7eff0, 0); // cf. [GP2004] Eq. (56), p. 10
        }

        Complex<double> c9eff(double s) const
        {
            // For r_i, g_i cf. [GP2004], Eqs. (27)-(29), p. 6
            // All occurences of sqrt(r_c) in Eq. need to be replaced by r_c. Cf. also the footnote
            // on p. 6.
            double r_b = std::sqrt(4.0 * m_b_MSbar * m_b_MSbar / s - 1.0);
            double r_c = std::sqrt(1.0 - 4.0 * s / m_c / m_c);
            Complex<double> g_0 = Complex<double>::Cartesian(std::log(s / mu / mu), -M_PI) * (1.0 / 6.0) - 5.0 / 18.0;
            double g_m_b = std::log(m_b_MSbar * m_b_MSbar / mu / mu) / 6.0 - 5.0 / 18.0 - 2.0 * m_b_MSbar * m_b_MSbar / 3.0 / s
                + r_b / 3.0 * (1.0 + 2.0 * m_b_MSbar * m_b_MSbar / s) * std::atan(1.0 / r_b);
            Complex<double> g_m_c = std::log(m_c * m_c / mu / mu) / 6.0 - 5.0 / 18.0 - 2.0 * m_c * m_c / 3.0 / s
                + r_c / 6.0 * (1.0 + 2.0 * m_c * m_c / s) * Complex<double>::Cartesian(std::log((1.0 + r_c) / (1.0 - r_c)), -M_PI);

            Complex<double> c9eff0 = c9() - (c1 + c2 / 3.0) * (8.0 * g_0 - 4.0 / 3.0) - c3() * (20.0 / 3.0 * g_0 - 16.0 / 3.0 * g_m_b + 2.0 / 27.0)
                + c4() * (4.0 / 3.0 * g_0 + 16.0 / 3.0 * g_m_b + 14.0 / 9.0) - c5() * (8.0 * g_0 - 4.0 * g_m_b - 14.0 / 27.0)
                - c6() * (8.0 / 3.0 * g_0 - 4.0 / 3.0 * g_m_b + 2.0 / 9.0);

            // TODO: Neglecting contributions ~alpha_s / 4.0 / M_PI. These do involve spectator scattering,
            // cf. [BFS2001] Eq. (29), p. 8, and Eqs. (82)-(84), p. 30
            return c9eff0; // cf. [GP2004] Eq. (55), p. 10
        }

        double kappa_1() const
        {
            static const double c0v // cf. [GP] Eq. (48)
                = 1.0 - QCD::alpha_s(mu) * QCD::casimir_f / (4.0 * M_PI) * (3.0 * std::log(mu / m_b_MSbar) + 4.0);
            static const double d0v // cf. [GP] Eq. (A30)
                = QCD::alpha_s(mu) * QCD::casimir_f / (2.0 * M_PI) * (std::log(mu / m_b_MSbar) + 1.0);

            // TODO: [GP] uses m_b(\mu). Which m_b? Using m_b_MSbar for the time
            // being.
            return (1.0 + 2.0 * d0v / c0v) * m_b_MSbar / m_B; // cf. [GP] Eq. (A24)
        }

        double norm(const double & s) const
        {
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]
            static const double g_fermi = 1.16637e-5; // (Gev^-2 (hbar c)^3), cf. [PDG2008], p.5
            double lambda_t = ckm_A * ckm_lambda * ckm_lambda;

            return std::sqrt(g_fermi * g_fermi * alpha_e * alpha_e / 3.0 / 1024 / std::pow(M_PI, 5.0) / m_B
                    * lambda_t * lambda_t * s_hat(s)
                    * std::sqrt(lambda(m_B * m_B, m_Kstar * m_Kstar, s))); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(const double & s) const
        {
            return s / m_B / m_B;
        }

        // Amplitudes
        Complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstar_hat = m_Kstar / m_B;
            Complex<double> wilson = (c9eff(s) - c9prime()) + h * (c10() - c10prime()) + kappa_1() * (c7eff(s) + c7prime()) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, -0.5 * m_B * m_B * m_B / m_Kstar / std::sqrt(s));
            double formfactor = lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s)) * form_factors->a_1(s_hat(s)) - (1 - s_hat(s)) * form_factors->a_2(s_hat(s));

            return prefactor * wilson * formfactor; // cf. [BHvD2010], Eq. (??)
        }

        Complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstar_hat = m_Kstar / m_B;
            Complex<double> wilson = (c9eff(s) + c9prime()) + h * (c10() + c10prime()) + kappa_1() * (c7eff(s) - c7prime()) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, std::sqrt(2 * lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s))) * m_B);

            return prefactor * wilson * form_factors->v(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }

        Complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            Complex<double> wilson = (c9eff(s) - c9prime()) + h * (c10() - c10prime()) + kappa_1() * (c7eff(s) - c7prime()) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, -std::sqrt(2) * m_B);

            return prefactor * wilson * form_factors->a_1(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }
    };

    BToKstarDilepton<LowRecoil>::BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options) :
        PrivateImplementationPattern<BToKstarDilepton<LowRecoil>>(new Implementation<BToKstarDilepton<LowRecoil>>(parameters, options))
    {
    }

    BToKstarDilepton<LowRecoil>::~BToKstarDilepton()
    {
    }

    Complex<double>
    BToKstarDilepton<LowRecoil>::a_long(const Helicity & h, const double & s) const
    {
        return _imp->a_long(h, s);
    }

    Complex<double>
    BToKstarDilepton<LowRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        return _imp->a_perp(h, s);
    }

    Complex<double>
    BToKstarDilepton<LowRecoil>::a_par(const Helicity & h, const double & s) const
    {
        return _imp->a_par(h, s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_branching_ratio(const double & s) const
    {
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        return differential_decay_width(s) * (_imp->norm(s) * _imp->norm(s) / Gamma);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_decay_width(const double & s) const
    {
        return (a_long(left_handed, s).absolute_squared()
            + a_long(right_handed, s).absolute_squared()
            + a_perp(left_handed, s).absolute_squared()
            + a_perp(right_handed, s).absolute_squared()
            + a_par(left_handed, s).absolute_squared()
            + a_par(right_handed, s).absolute_squared());
    }

    double
    BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 / differential_decay_width(s) * (
                (a_par(left_handed, s) * a_perp(left_handed, s).conjugate()).real()
                -(a_par(right_handed, s) * a_perp(right_handed, s).conjugate()).real()
            );
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        double a = a_perp(left_handed, s).absolute_squared() + a_perp(right_handed, s).absolute_squared();
        double b = a_par(left_handed, s).absolute_squared() + a_par(right_handed, s).absolute_squared();

        return (a - b) / (a + b);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        double a = (a_long(left_handed, s) * a_par(left_handed, s).conjugate() + a_long(right_handed, s).conjugate() * a_par(right_handed, s)).absolute();
        double b = std::sqrt((a_long(left_handed, s).absolute_squared() + a_long(right_handed, s).absolute_squared())
                * (a_perp(left_handed, s).absolute_squared() + a_perp(right_handed, s).absolute_squared()));

        return a / b;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        double a = (a_long(left_handed, s) * a_perp(left_handed, s).conjugate() - a_long(right_handed, s).conjugate() * a_perp(right_handed, s)).absolute();
        double b = (a_long(left_handed, s).conjugate() * a_par(left_handed, s) + a_long(right_handed, s) * a_par(right_handed, s).conjugate()).absolute();

        return a / b;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return (a_long(left_handed, s).absolute_squared() + a_long(right_handed, s).absolute_squared())
            / differential_decay_width(s);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::tr1::function<double (const double &)> f = std::tr1::bind(std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_branching_ratio),
                this, std::tr1::placeholders::_1);

        return integrate(f, 100, s_min, s_max);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        std::tr1::function<double (const double &)> f = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry), this, std::tr1::placeholders::_1);

        return integrate(f, 100, s_min, s_max);
    }
}
