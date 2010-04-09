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
    using std::norm;

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

        complex<double> lambda_u_hat;

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

        inline double mu_f() const
        {
            static const double Lambda_QCD = 0.5; // (GeV)
            return std::sqrt(mu * Lambda_QCD);
        }

        inline double m_b_pole() const
        {
            // Actually use the PS mass at mu_f = 2 GeV
            return 4.6;
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
        complex<double> h(const double & s, const double & m_q) const
        {
            if (m_q < 1e-4)
                return h0(s);

            const double z = 4.0 * m_q * m_q / s;
            if (z < 1e-10)
                return complex<double>(-4.0/9.0 * (2.0 * std::log(m_q / mu) + 1.0), 0.0);

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

            return -4.0 / 9.0 * (a + b * complex<double>(rc, ic));
        }

        // cf. [BFS2001], Eq. (11), p. 4 in the limit m_q -> 0
        complex<double> h0(const double & s) const
        {
            return 4.0 / 9.0 * complex<double>(2.0 / 3.0 + 2.0 * std::log(2.0 * mu()) - std::log(s), M_PI);
        }

        // cf. [BFS2001], Eq. (10), p. 4
        complex<double> Y0(const double & s) const
        {
            double Y_c = 4.0 / 3.0 * c1() + c2() + 6.0 * c3() + 60.0 * c5();
            double Y_b = -0.5 * (7.0 * c3() + 4.0 / 3.0 * c4() + 76.0 * c5() + 64.0 / 3.0 * c6());
            double Y_0 = -0.5 * (c3() + 4.0 / 3.0 * c4() + 16.0 * c5() + 64 / 3.0 * c6());
            double Y = 2.0 / 9.0 * (6.0 * c3() + 32.0 * c5() + 32.0 / 3.0 * c6());

            // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
            return Y_c * h(s, m_c) + Y_b * h(s, m_b_pole()) + Y_0 * h0(s) + Y;
        }

        // cf. [BFS2004], ?
        complex<double> Y0u(const double & s) const
        {
            double a = 4.0 / 3.0 * c1() + c2();

            return a * (h(s, m_c) - h0(s));
        }

        /* Form factors */
        //  cf. [BHP2008], Eq. (E.4), p. 23
        double xi_perp(const double & s) const
        {
            const double factor = m_B / (m_B + m_Kstar);
            double result = factor * form_factors->v(s_hat(s));

            return result;
        }

        double xi_par(const double & s) const
        {
            const double factor1 = (m_B + m_Kstar) / (2.0 * energy(s));
            const double factor2 = (1.0 - m_Kstar / m_B);
            double result = factor1 * form_factors->a_1(s_hat(s)) - factor2 * form_factors->a_2(s_hat(s));

            return result;
        }

        /* NLO functions */
        // cf. [BFS2001], Eq. (29), p. 8
        complex<double> B0(const double & s, const double & m_q) const
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
                return complex<double>(-2.0, 0.0);

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

            return complex<double>(rp, ip);
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
        complex<double> I1(const double & s, const double & u, const double & m_q) const
        {
            if (m_q == 0.0)
                return complex<double>(1.0, 0.0);

            int status;
            gsl_sf_result res_re, res_im;

            double ubar = 1.0 - u;
            double m_q2 = m_q * m_q;
            double m_B2 = m_B * m_B;

            double a, a2, sign;
            complex<double> dilogArg, dilog1, dilog2;
            complex<double> LxpLxm, LypLym;

            if (1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s)) > 0)
            {
                a = (1 - std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s))))
                  / (1 + std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s))));
                LxpLxm = -M_PI* M_PI / 3.0
                    + std::log(a) * (std::log(a) + complex<double>(0.0, M_PI))
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

                dilogArg = complex<double>((a2 - 1.0)/(a2 + 1.0), -2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                dilog1 = complex<double>( res_re.val, res_im.val);

                dilogArg = complex<double>((a2 - 1.0)/(a2 + 1.0), +2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                dilog2 = complex<double>(res_re.val, res_im.val);

                LxpLxm = -1.0 / 3.0 * M_PI * M_PI - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                    + dilog1 + dilog2;
            }

            if (1.0 - 4.0 * m_q2 / s > 0)
            {
                a = (1.0 - std::sqrt(1.0 - 4.0 * m_q2 / s))
                  / (1.0 + std::sqrt(1.0 - 4.0 * m_q2 / s));
                LypLym = -1.0 / 3.0 * M_PI * M_PI + std::log(a) * (std::log(a) + complex<double>(0.0, M_PI))
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

                dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), -2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                dilog1 = complex<double>(res_re.val, res_im.val);

                dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), +2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                dilog2 = complex<double>(res_re.val, res_im.val);

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
            double wilson = c3 + 4.0/3.0 * c4 + 16.0 * c5 + 64.0/3.0 * c6;
            return -1.0 * e_q * (4.0 * m_B / m_b_pole()) * wilson;
        }

        // cf. [BFS2001], Eq. (20), p. 6, module the inverse of lambda_B,+
        double Tf_perp_p(const double & h, const double & s, const double & u) const
        {
            return (c7eff() + h * c7prime()) * (2.0 * m_B / (1.0 - u) / energy(s));
        }

        // cf. [BFS2001], Eq. (21), p. 6
        complex<double> Tf_par_p(const double & s, const double & u) const
        {
            double E = energy(s);

            return ((c7eff() - c7prime) + (s / (2.0 * m_B * m_b_pole())) * Y0(s)) * (2.0 * m_B * m_B / (1 - u) / E / E);
        }

        // cf. [BFS2001], Eq. (27), p. 8
        complex<double> t_perp(const double & s, const double & u, const double & m_q) const
        {
            if (0.0 == s)
                return t_perp_0(u, m_q);

            double ubar = 1.0 - u;
            double E = energy(s);
            double x = ubar * m_B * m_B + u * s;

            complex<double> result = (2.0 * m_B / ubar / E) * I1(s, u, m_q);
            if (m_q > 0.0)
                result = result + (s / ubar / ubar / E / E) * (B0(x, m_q) - B0(s, m_q));

            return result;
        }

        complex<double> t_perp_0(const double & u, const double & m_q) const
        {
            double ubar = 1.0 - u;
            double m_q2 = m_q * m_q, m_B2 = m_B * m_B;
            double a, a2, sign;
            complex<double> dilogArg, dilog1, dilog2;
            complex<double> LxpLxm;
            int status;
            gsl_sf_result res_re, res_im;

            if (m_q > 0)
            { // m != 0
                if (1.0 - 4.0 * m_q2 / (m_B2 - u * m_B2) > 0)
                {
                    a = (1.0 - std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * m_B2)))
                        / (1.0 + std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * m_B2)));
                    LxpLxm = -M_PI * M_PI / 3.0 + std::log(a) * (std::log(a) + complex<double>(0.0, M_PI))
                        + gsl_sf_dilog(-a) + gsl_sf_dilog(-1.0/a);
                }
                else
                {
                    a = std::sqrt(4.0 * m_q2 / (m_B2 - u * m_B2) - 1.0);
                    a2 = a * a;

                    if (a2 - 1 > 0)
                        sign = +1.0;
                    else
                        sign = -1.0;

                    dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), -2.0 * a / (a2 + 1.0));
                    status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                    dilog1 = complex<double>(res_re.val, res_im.val);
                    dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), +2.0 * a / (a2 + 1.0));
                    status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                    dilog2 = complex<double>(res_re.val, res_im.val);

                    LxpLxm = -M_PI * M_PI / 3.0 - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                        + dilog1 + dilog2;
                }

                return 4.0 / ubar * (1.0 + 2.0 * m_q2 / ubar / m_B2 * LxpLxm);
            }
            else
            {
                return complex<double>(4.0 / ubar, 0.0);
            }
        }

        // cf. [BFS2001], Eq. (28), p. 8
        complex<double> t_par(const double & s, const double & u, const double & m_q) const
        {
            double ubar = 1.0 - u;
            double E = energy(s);
            double x = ubar * m_B * m_B + u * s;

            complex<double> result = (2.0 * m_B / ubar / E) * I1(s, u, m_q);
            if (m_q > 0.0)
                result = result + (x / ubar / ubar / E / E) * (B0(x, m_q) - B0(s, m_q));

            return result;
        }

        // cf. [BFS2001], Eq. (23), p. 7, multiplied by phi_K^*,perp
        complex<double> Tnf_perp_p(const double & s, const double & u) const
        {
            double s_hat = s / m_B / m_B;
            double ubar = 1.0 - u;

            double a = (4.0 / 3.0 / (u + ubar * s_hat)) * c8eff();
            complex<double> ba = (+2.0 / 3.0) * (-c1 / 6.0 + c2 + 6.0 * c6) * t_perp(s, u, m_c);
            complex<double> bb = (-1.0 / 3.0)
                * (c3 - c4 / 6.0 + 16.0 * c5 + 10.0/3.0 * c6 - (4.0 * m_b_pole() / m_B) * (c3 - c4/6.0 + 4.0 * c5 - 2.0/3.0 * c6))
                * t_perp(s, u, m_b_pole());
            complex<double> bc = (-1.0 / 3.0) * (c3 - c4 / 6.0 + 16.0 * c5 - 8.0/3.0 * c6) * t_perp(s, u, 0.0);

            return (a + (m_B / 2.0 / m_b_pole()) * (ba + bb + bc)) * phi_K(u, a_1_perp, a_2_perp);
        }

        // cf. [BFS2001, Eq. (25), p. 7, multiplied by phi_K^*,par
        complex<double> Tnf_par_p(const double & s, const double & u) const
        {
            complex<double> a = (+2.0 / 3.0) * (-c1 / 6.0 + c2 + 6.0 * c6) * t_par(s, u, m_c);
            complex<double> b = (-1.0 / 3.0) * (c3 - c1 / 6.0 + 16.0 * c5 + 10.0/3.0 * c6) * t_par(s, u, m_b_pole());
            complex<double> c = (-1.0 / 3.0) * (c3 - c4 / 6.0 + 16.0 * c5 - 8.0/3.0 * c6) * t_par(s, u, 0.0);

            return (m_B / m_b_pole()) * (a + b + c) * phi_K(u, a_1_par, a_2_par);
        }

        // cf. [BFS2001], Eq. (26), pp. 7-8, multiplied by phi_K^*,par
        complex<double> Tnf_par_m(const double & s, const double & u) const
        {
            double s_hat = s / m_B / m_B;
            double ubar = 1.0 - u;
            double x = ubar * m_B * m_B + u * s;

            double a = (e_q * 8.0 / (ubar + u * s_hat)) * c8eff();
            complex<double> ba = (-c1 / 6.0 + c2 + c4 + 10 * c6) * h(x, m_c);
            complex<double> bb = (c3 + 5.0/6.0 * c4 + 16.0 * c5 + 22.0/3.0 * c6) * h(x, m_b_pole());
            complex<double> bc = (c3 + 17.0/6.0 * c4 + 16.0 * c5 + 82.0/3.0 * c6) * h0(x);
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
        complex<double> lambda_B_m_inv(const double & s) const
        {
            if (0.0 == s)
                return complex<double>(0.0, 0.0);

            double omega_0 = lambda_B_p;
            double x = s / m_B / omega_0;
            double ei = gsl_sf_expint_Ei(x);

            complex<double> result = complex<double>(-ei, M_PI) * (std::exp(-x) / omega_0);

            return result;
        }

        // cf. [AAGW2001], Eq. (56), p. 20
        complex<double> F17(const double & s) const
        {
            // Deviating definition of s_hat == s / m_b_pole^2 !
            double s_hat = s / m_b_pole() / m_b_pole();
            double Lmu = std::log(mu / m_b_pole());
            double Ls = std::log(s_hat);

            // cf. [AAGW2001], Table I, p. 20
            std::vector<std::pair<complex<double>, complex<double>>> f17_coefficients =
            {
                std::make_pair(complex<double>(-0.76730, -0.11418), complex<double>(-0.00000, -0.00000)),
                std::make_pair(complex<double>(-0.28480, -0.18278), complex<double>(-0.00328, +0.02083)),
                std::make_pair(complex<double>(+0.05611, -0.23357), complex<double>(+0.01637, +0.02091)),
                std::make_pair(complex<double>(+0.62438, -0.02744), complex<double>(+0.07673, +0.00914))
            };

            complex<double> f17 = complex<double>(0.0, 0.0);
            double i(0);
            for (auto k = f17_coefficients.begin() ; k != f17_coefficients.end() ; ++k, ++i)
            {
                f17 = f17 + std::pow(s_hat, i) * (k->first + k->second * Ls);
            }

            return -208.0/243.0 * Lmu + f17;
        }

        // cf. [AAGW2001], Eq. (54), p. 19
        complex<double> F19(const double & s) const
        {
            // Deviating definition of s_hat == s / m_b_pole^2 !
            double s_hat = s / m_b_pole() / m_b_pole();
            double Lmu = std::log(mu / m_b_pole());
            double Ls = std::log(s_hat);
            double Lc = std::log(m_c / m_b_pole());
            double z = (m_c / m_b_pole()) * (m_c / m_b_pole());

            // cf. [AAGW2001], Table I, p. 20
            std::vector<std::pair<complex<double>, complex<double>>> f19_coefficients =
            {
                std::make_pair(complex<double>(-12.7150, +0.09470), complex<double>(-0.07883, -0.07414)),
                std::make_pair(complex<double>(-38.7420, -0.67862), complex<double>(-0.03930, -0.00017)),
                std::make_pair(complex<double>(-103.830, -2.53880), complex<double>(-0.04470, +0.00263)),
                std::make_pair(complex<double>(-313.750, -8.45540), complex<double>(-0.05113, +0.02275))
            };

            complex<double> f19 = complex<double>(0.0, 0.0);
            double i(0);
            for (auto k = f19_coefficients.begin() ; k != f19_coefficients.end() ; ++k, ++i)
            {
                f19 = f19 + std::pow(s_hat, i) * (k->first + k->second * Ls);
            }

            complex<double> result =
                        Lmu * complex<double>(-1424.0/729.0 + 64.0/27.0 * Lc - 16.0/243.0 * Ls + (16.0/1215.0 - 32.0/(z * 135.0)) * s_hat
                        + (4.0/2835.0 - 8.0/(315.0 * z * z)) * s_hat * s_hat + (16.0/76545.0 - 32.0/(8505.0 * z * z * z)) * s_hat * s_hat * s_hat
                        - 256.0/243.0 * Lmu,
                        16.0/243.0 * M_PI)
                + f19;

            return result;
        }

        // cf. [AAGW2001], Eq. (56), p. 20
        complex<double> F27(const double & s) const
        {
            // Deviating definition of s_hat == s / m_b_pole^2 !
            double s_hat = s / m_b_pole() / m_b_pole();
            double Lmu = std::log(mu / m_b_pole());
            double Ls = std::log(s_hat);

            // cf. [AAGW2001], Table II, p. 21
            std::vector<std::pair<complex<double>, complex<double>>> f27_coefficients =
            {
                std::make_pair(complex<double>(+4.60380, +0.68510), complex<double>(-0.00000, -0.00000)),
                std::make_pair(complex<double>(+1.70880, +1.09670), complex<double>(+0.01959, -0.12496)),
                std::make_pair(complex<double>(-0.33665, +1.40140), complex<double>(-0.09822, -0.12548)),
                std::make_pair(complex<double>(-3.74630, +0.16463), complex<double>(-0.18321, -0.05485))
            };

            complex<double> f27 = complex<double>(0.0, 0.0);
            double i(0);
            for (auto k = f27_coefficients.begin() ; k != f27_coefficients.end() ; ++k, ++i)
            {
                f27 = f27 + std::pow(s_hat, i) * (k->first + k->second * Ls);
            }

            complex<double> result = 416.0/81.0 * Lmu + f27;

            return result;
        }

        // cf. [AAGW2001], Eq. (55), p. 19
        complex<double> F29(const double & s) const
        {
            // Deviating definition of s_hat == s / m_b_pole^2 !
            double s_hat = s / m_b_pole() / m_b_pole();
            double Lmu = std::log(mu / m_b_pole());
            double Ls = std::log(s_hat);
            double Lc = std::log(m_c / m_b_pole());
            double z = (m_c / m_b_pole()) * (m_c / m_b_pole());

            // cf. [AAGW2001], Table II, p. 21
            std::vector<std::pair<complex<double>, complex<double>>> f29_coefficients =
            {
                std::make_pair(complex<double>(+9.50420, -0.56819), complex<double>(+0.47298, +0.44483)),
                std::make_pair(complex<double>(+7.42380, +4.07170), complex<double>(+0.23581, +0.00104)),
                std::make_pair(complex<double>(+0.33806, +15.2330), complex<double>(+0.26821, -0.01577)),
                std::make_pair(complex<double>(-42.0850, +50.7320), complex<double>(+0.30680, -0.13652))
            };

            complex<double> f29 = complex<double>(0.0, 0.0);
            double i(0);
            for (auto k = f29_coefficients.begin() ; k != f29_coefficients.end() ; ++k, ++i)
            {
                f29 = f29 + std::pow(s_hat, i) * (k->first + k->second * Ls);
            }

            complex<double> result =
                Lmu * complex<double>(256.0/243.0 - 128.0/9.0 * Lc + 32.0/81.0 * Ls - (32.0/405.0 - 64.0/(z * 45.0)) * s_hat
                        - (8.0/945.0 - 16.0/(105.0 * z * z)) * s_hat * s_hat - (32.0/25515.0 - 64.0/(2835.0 * z * z * z)) * s_hat * s_hat * s_hat
                        + 512.0/81.0 * Lmu,
                        -32.0/81.0 * M_PI)
                + f29;

            return result;
        }

        // cf. [BFS2001], Eq. (82), p. 30
        complex<double> F87(const double & s) const
        {
            double s_hat = s / m_B / m_B;
            double s_hat2 = s_hat * s_hat;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            complex<double> a = complex<double>(-32.0 * std::log(mu / m_b_pole()) - 8.0 * s_hat / denom * std::log(s_hat)
                    - 4.0 * (11.0 - 16.0 * s_hat + 8.0 * s_hat2) / denom2,
                    -8.0 * M_PI);
            complex<double> b = (4.0 / denom / denom2)
                * ((9.0 * s_hat - 5.0 * s_hat2 + 2.0 * s_hat * s_hat2) * B0(s, m_b_pole()) - (4.0 + 2.0 * s_hat) * C0(s));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [BFS2001], Eq. (83), p. 30
        complex<double> F89(const double & s) const
        {
            double s_hat = s / m_B / m_B;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            double a = 16.0 * std::log(s_hat) / denom + 8.0 * (5.0 - 2.0 * s_hat) / denom2;
            complex<double> b = (-8.0 * (4.0 - s_hat) / denom / denom2) * ((1.0 + s_hat) * B0(s, m_b_pole()) - 2.0 * C0(s));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [BFS2004], Eq. (51) the integrand of the first term only,
        // or [FM2002], Eq. (15) times a factor of 3.0, respectively
        double Twa_perp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;

            return phi_K(u, a_1_perp, a_2_perp) / (ubar + u * s_hat);
        }

        // cf. [FM2002], Eq. (17), the integrand only
        double Xperp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;
            double denom = ubar + u * s_hat;

            return phi_K(u, a_1_perp, a_2_perp) * (1 / denom + 1 / denom / denom) / 3.0;
        }

        // cf. [FM2002], Eq. (22), p. 9
        complex<double> FV(const double & s) const
        {
            return 3.0 / 4.0 * (
                    (-c1 / 6.0 + c2 + c4 + 10.0 * c6) * h(s, m_c)
                    + (c3 + 5.0/6.0 * c4 + 16.0 * c5 + 22.0/3.0 * c6) * h(s, m_b_pole())
                    + (c3 + 17.0/6.0 * c4 + 16.0 * c5 + 82.0/3.0 * c6) * h0(s)
                    - 8.0/27.0 * (-7.5 * c4 + 12 * c5 - 32 * c6));
        }

        // cf. [BFS2004], Eq. (52), the integrand of the second term only
        complex<double> Thsa_1_perp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;
            double x = (ubar + u * s_hat) * m_B * m_B;

            return phi_K(u, a_1_perp, a_2_perp) / (ubar + u * s_hat) * FV(x);
        }

        // cf. [BFS2004], Eq. (52), the integrand of the third term only.
        // the v integration has been executed analytically
        complex<double> Thsa_2_perp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;
            double x = (ubar + u * s_hat) * m_B * m_B;

            return 3.0 * u * u * (1.0 + a_1_par * (4.0 * u - 3.0) + a_2_par * (15.0 * u * u - 20.0 * u + 6.0))
                * FV(x);
        }

        // cf. [BFS2001], Eq. (15) with a = perp, and [BHP2008], Eq. (C.4)
        complex<double> tensor_perp(const double & h, const double & s) const
        {
            // cf. [BFS2004], paragraph below Eq. (42)
            double ff_nlo_factor = QCD::alpha_s(mu) * QCD::casimir_f / 4.0 / M_PI;
            double scatt_nlo_factor = QCD::alpha_s(mu_f()) * QCD::casimir_f / 4.0 / M_PI;

            double scatt_factor = M_PI * M_PI / 3.0 * f_B * f_Kstar_perp / m_B;
            double shat = s_hat(s);

            // Here m_b_pole is used, cf. [BFS2001], comment below Eq. (36), p. 9

            /* Form factor corrections */
            complex<double> ff_0 = (c7eff() + h * c7prime()) + s / (2.0 * m_b_pole() * m_B) * Y0(s);
            // cf. [BFS2001], Eq. (34), p. 9
            double ff_f = (c7eff() + h * c7prime()) * (8.0 * std::log(m_b_pole() / mu) - 4.0 * (1.0 - mu_f() / m_b_pole()) - L(s));
            // cf. [BFS2001], Eq. (37), p. 9
            complex<double> ff_nf = (-1.0 / QCD::casimir_f) * (
                    (c2 - c1 / 6.0) * F27(s) + c8eff() * F87(s)
                    + (s / (2.0 * m_b_pole() * m_B)) * (c2() * F29(s) + c1() * F19(s) + c8eff() * F89(s)));
            complex<double> ff = ff_0 + ff_nlo_factor * (ff_f + ff_nf);

            /* Specator scattering, folded with phi_K^*,perp */
            // cf. [BFS2001], Eq. (20), p. 6
            double scatt_f_p = (c7eff() + h * c7prime) * 2.0 * m_B / energy(s) * 3.0 * (1.0 + a_1_perp + a_2_perp);
            // cf. [BFS2001], Eq. (23), p. 7
            complex<double> scatt_nf_p = integrate(
                    std::tr1::function<complex<double> (const double &)>(
                        std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Tnf_perp_p, this, s, std::tr1::placeholders::_1)),
                    40, 0.001, 0.999);
            complex<double> scatt_p = (1.0 / lambda_B_p) * scatt_nlo_factor * (scatt_f_p + scatt_nf_p);

            /* Weak annihilation */
            // cf. [BFS2004], Eq. (51), p. 26
            double wa = (e_q * 2.0 * M_PI * M_PI * f_B / (3.0 * m_b_pole() * m_B)) * (
                    - f_Kstar_perp * (c3 + 4.0/3.0 * c4 + 4.0 * c5 + 16.0/3.0 * c6) * integrate(
                        std::tr1::function<double (const double &)>(
                            std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Twa_perp, this, shat, std::tr1::placeholders::_1)),
                        40, 0.0, 1.0)
                    + f_Kstar_par * m_Kstar / ((1 - shat) * lambda_B_p));
            /* Hard spectator scattering */
            complex<double> hsa = e_q * scatt_nlo_factor * M_PI * M_PI * f_B / (3.0 * m_b_pole() * m_B) * (
                    12.0 * c8eff() * m_b_pole() / m_B * f_Kstar_perp * integrate(
                        std::tr1::function<double (const double &)>(
                            std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Xperp, this, shat, std::tr1::placeholders::_1)),
                        40, 0.0, 1.0)
                    + 8.0 * f_Kstar_perp * integrate(
                        std::tr1::function<complex<double> (const double &)>(
                            std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Thsa_1_perp, this, shat, std::tr1::placeholders::_1)),
                        40, 0.0, 1.0)
                    - 4.0 * m_Kstar * f_Kstar_par / ((1.0 - shat) * lambda_B_p) * integrate(
                        std::tr1::function<complex<double> (const double &)>(
                            std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Thsa_2_perp, this, shat, std::tr1::placeholders::_1)),
                        40, 0.0, 1.0));

            complex<double> result = xi_perp(s) * ff + scatt_factor * scatt_p + wa + hsa;

            return result;
        }

        // cf. [BFS2001], Eq. (15) with a = par, and [BHP2008], Eq. (C.4)
        complex<double> tensor_par(const double & s) const
        {
            // cf. [BFS2004], paragraph below Eq. (42)
            double ff_nlo_factor = QCD::alpha_s(mu) * QCD::casimir_f / 4.0 / M_PI;
            double scatt_nlo_factor = QCD::alpha_s(mu_f()) * QCD::casimir_f / 4.0 / M_PI;

            double scatt_factor = M_PI * M_PI / 3.0 * f_B * f_Kstar_perp / m_B * m_Kstar / energy(s);

            // Here m_b_pole is used, cf. [BFS2001], comment below Eq. (36), p. 9

            /* Form factor corrections */
            complex<double> ff_0 = -1.0 * (c7eff() - c7prime() + m_B / (2.0 * m_b_pole()) * (Y0(s)));
            // cf. [BFS2004], Eq. (44), p. 24
            complex<double> ff_f = -1.0 * (c7eff() - c7prime()) * (8.0 * std::log(m_b_pole() / mu) + 2.0 * L(s) - 4.0 * (1.0 - mu_f() / m_b_pole()))
                    + (m_B / (2.0 * m_b_pole())) * Y0(s) * (2.0 - 2.0 * L(s));
            // cf. [BFS2001], Eq. (38), p. 9
            complex<double> ff_nf = (+1.0 / QCD::casimir_f) * (
                    (c2 - c1 / 6.0) * F27(s) + c8eff() * F87(s)
                    + (m_B / (2.0 * m_b_pole())) * (c2() * F29(s) + c1() * F19(s) + c8eff() * F89(s)));
            complex<double> ff = ff_0 + ff_nlo_factor * (ff_f + ff_nf);

            /* Spectator scattering */
            // cf. [BFS2001], Eq. (18), p. 6
            double scatt_0_m = -e_q * 4.0 * m_B / m_b_pole() * (c3 + 4.0/3.0 * c4 + 16.0 * c5 + 64.0/3.0 * c6);
            // cf. [BFS2001], Eq. (21), p. 6
            complex<double> scatt_f_p = (c7eff() - c7prime + (s / (2.0 * m_b_pole() * m_B)) * Y0(s)) * (2.0 * m_B * m_B / energy(s) / energy(s))
                    * 3.0 * (1.0 + a_1_par + a_2_par);
            scatt_f_p = complex<double>(((c7eff() - c7prime) * 4.0 * m_B / energy(s)) * 3.0 * (1.0 + a_1_par + a_2_par), 0.0);
            // cf. [BFS2001], Eq. (25), p. 7
            complex<double> scatt_nf_p = integrate(
                    std::tr1::function<complex<double> (const double &)>(std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Tnf_par_p,
                            this, s, std::tr1::placeholders::_1)),
                    40, 0.001, 0.999);
            // cf. [BFS2001], Eq. (26), pp. 7-8
            complex<double> scatt_nf_m = integrate(
                    std::tr1::function<complex<double> (const double &)>(std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Tnf_par_m,
                            this, s, std::tr1::placeholders::_1)),
                    40, 0.001, 0.999);
            complex<double> scatt_p = (1.0 / lambda_B_p) * scatt_nlo_factor * (scatt_f_p + scatt_nf_p);
            complex<double> scatt_m = lambda_B_m_inv(s)  * (scatt_0_m + scatt_nlo_factor * scatt_nf_m);

            complex<double> result = xi_par(s) * ff + scatt_factor * (scatt_p + scatt_m);

            return result;
        }

        /* Amplitudes */
        // cf. [BHP2008], p. 20
        complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_pole() / m_B;
            double E = 0.5 * (m_B - s / m_B);
            double mKhat = m_Kstar / m_B;
            double lambdahat = lambda(1.0, mKhat * mKhat, s);

            double wilson = (c9() - c9prime()) + h * (c10() - c10prime());
            double prefactor = -1.0 / (2.0 * m_Kstar * std::sqrt(s));

            double a = wilson * ((m_B * m_B - m_Kstar * m_Kstar - s) * 2.0 * energy(s) * xi_perp(s)
                -lambda(m_B * m_B, m_Kstar * m_Kstar, s) * m_B / (m_B * m_B - m_Kstar * m_Kstar) * (xi_perp(s) - xi_par(s)));
            complex<double> b = 2 * m_b_pole() * (((m_B * m_B + 3 * m_Kstar * m_Kstar - s) * 2.0 * energy(s) / m_B
                        - lambda(m_B * m_B, m_Kstar * m_Kstar, s) / (m_B * m_B - m_Kstar * m_Kstar)) * tensor_perp(-1.0, s)
                    - lambda(m_B * m_B, m_Kstar * m_Kstar, s) / (m_B * m_B - m_Kstar * m_Kstar) * tensor_par(s));
            return prefactor * (a + b);
        }

        // cf. [BHP2008], p. 20
        complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_pole() / m_B;
            double mKhat = m_Kstar / m_B;

            double prefactor = +std::sqrt(2.0) * m_B * std::sqrt(lambda(1.0, mKhat * mKhat, shat));
            double wilson = (c9() + c9prime()) + h * (c10() + c10prime());

            return prefactor * (wilson * xi_perp(s) + (2.0 * mbhat / shat) * tensor_perp(+1.0, s));
        }

        // cf. [BHP2008], p. 20
        complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_pole() / m_B;
            double mKhat = m_Kstar / m_B;

            double prefactor = -std::sqrt(2.0) * m_B * (1.0 - shat);
            double wilson = (c9() - c9prime()) + h * (c10() - c10prime());

            return prefactor * (wilson * xi_perp(s) + (2.0 * mbhat / shat) * (1.0 - mKhat * mKhat) * tensor_perp(-1.0, s));
        }
    };

    BToKstarDilepton<LargeRecoil>::BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options) :
        PrivateImplementationPattern<BToKstarDilepton<LargeRecoil>>(new Implementation<BToKstarDilepton<LargeRecoil>>(parameters, options))
    {
    }

    BToKstarDilepton<LargeRecoil>::~BToKstarDilepton()
    {
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_long(const Helicity & h, const double & s) const
    {
        return _imp->a_long(h, s);
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        return _imp->a_perp(h, s);
    }

    complex<double>
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
        return _imp->norm(s) * _imp->norm(s) * (norm(a_long(left_handed, s))
            + norm(a_long(right_handed, s))
            + norm(a_perp(left_handed, s))
            + norm(a_perp(right_handed, s))
            + norm(a_par(left_handed, s))
            + norm(a_par(right_handed, s)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 * _imp->norm(s) * _imp->norm(s) / differential_decay_width(s)
            * (real(a_par(left_handed, s) * conj(a_perp(left_handed, s))) - real(a_par(right_handed, s) * conj(a_perp(right_handed, s))));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_unnormalized_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 * _imp->norm(s) * _imp->norm(s)
            * (real(a_par(left_handed, s) * conj(a_perp(left_handed, s))) - real(a_par(right_handed, s) * conj(a_perp(right_handed, s))));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        double a = norm(a_perp(left_handed, s)) + norm(a_perp(right_handed, s));
        double b = norm(a_par(left_handed, s)) + norm(a_par(right_handed, s));

        return (a - b) / (a + b);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        double a = abs(a_long(left_handed, s) * conj(a_par(left_handed, s)) + conj(a_long(right_handed, s)) * a_par(right_handed, s));
        double b = std::sqrt(norm(a_long(left_handed, s)) + norm(a_long(right_handed, s)))
                * std::sqrt(norm(a_perp(left_handed, s)) + norm(a_perp(right_handed, s)));

        return a / b;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        double a = abs(a_long(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_long(right_handed, s)) * a_perp(right_handed, s));
        double b = abs(conj(a_long(left_handed, s)) * a_par(left_handed, s) + a_long(right_handed, s) * conj(a_par(right_handed, s)));

        return a / b;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return (norm(a_long(left_handed, s)) + norm(a_long(right_handed, s))) * _imp->norm(s) * _imp->norm(s) / differential_decay_width(s);
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
        std::tr1::function<double (const double &)> num = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_unnormalized_forward_backward_asymmetry), this, std::tr1::placeholders::_1);

        std::tr1::function<double (const double &)> denom = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_decay_width), this, std::tr1::placeholders::_1);

        return integrate(num, 100, s_min, s_max) / integrate(denom, 100, s_min, s_max);
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

        Parameter c8;

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
            ckm_A(p["CKM::A"]),
            ckm_lambda(p["CKM::lambda"])
        {
            form_factors = FormFactorFactory<BToKstar>::create(o["form-factors"], p);
            if (! form_factors.get())
                throw std::string("InternalError");
        }

        /* Charm loop parts */
        // cf. [BFS2001], Eq. (11), p. 4
        complex<double> h(const double & s, const double & m_q) const
        {
            if (m_q < 1e-4)
                return h0(s);

            const double z = 4.0 * m_q * m_q / s;
            if (z < 1e-10)
                return complex<double>(-4.0/9.0 * (2.0 * std::log(m_q / mu) + 1.0), 0.0);

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

            return -4.0 / 9.0 * (a + b * complex<double>(rc, ic));
        }

        // cf. [BFS2001], Eq. (11), p. 4 in the limit m_q -> 0
        complex<double> h0(const double & s) const
        {
            return 4.0 / 9.0 * complex<double>(2.0 / 3.0 + 2.0 * std::log(2.0 * mu()) - std::log(s), M_PI);
        }

        // cf. [BFS2001], Eq. (29), p. 8
        complex<double> B0(const double & s, const double & m_q) const
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
                return complex<double>(-2.0, 0.0);

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

            return complex<double>(rp, ip);
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

        // cf. [GP2004], between Eqs. (41) and (42), p. 8
        complex<double> A(const double & s) const
        {
            return -104.0/243.0 * 2.0 * std::log(m_b_MSbar / mu) + complex<double>(0.736, 0.836);
        }

        // cf. [GP2004], between Eqs. (41) and (42), p. 8
        complex<double> B(const double & s) const
        {
            const double l = std::log(std::pow(m_b_MSbar / mu, 2.0));
            const double z = 4 * std::pow(m_b_MSbar, 2.0) / s;

            complex<double> result = complex<double>(z - 34.0, -17 * M_PI) * l + 8.0 * l * l - 17.0 * std::log(z / 4.0) * l
                - (2.0 + z) * std::sqrt(z - 1.0) * std::atan(1.0 / (std::sqrt(z - 1.0))) * l;
            result = 8.0 / 243.0 * result;
            result = result + complex<double>(-1.332, 3.058);

            return result;
        }

        // cf. [GP2004], between Eqs. (41) and (42), p. 8
        complex<double> C(const double & s) const
        {
            const double zeta3 = 1.20206;

            return complex<double>(-16.0/81.0 * std::log(s / std::pow(mu, 2.0)) + 428.0/243.0 - 64.0/27.0 * zeta3, 16.0 / 81.0 * M_PI);
        }

        // cf. [BFS2001], Eq. (82), p. 30
        complex<double> F87(const double & s) const
        {
            double s_hat = s / m_B / m_B;
            double s_hat2 = s_hat * s_hat;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            complex<double> a = complex<double>(-32.0 * std::log(mu / m_b_MSbar) - 8.0 * s_hat / denom * std::log(s_hat)
                    - 4.0 * (11.0 - 16.0 * s_hat + 8.0 * s_hat2) / denom2,
                    -8.0 * M_PI);
            complex<double> b = (4.0 / denom / denom2)
                * ((9.0 * s_hat - 5.0 * s_hat2 + 2.0 * s_hat * s_hat2) * B0(s, m_b_MSbar) - (4.0 + 2.0 * s_hat) * C0(s));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [BFS2001], Eq. (83), p. 30
        complex<double> F89(const double & s) const
        {
            double s_hat = s / m_B / m_B;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            double a = 16.0 * std::log(s_hat) / denom + 8.0 * (5.0 - 2.0 * s_hat) / denom2;
            complex<double> b = (-8.0 * (4.0 - s_hat) / denom / denom2) * ((1.0 + s_hat) * B0(s, m_b_MSbar) - 2.0 * C0(s));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [GP2004], Eq. (56)
        complex<double> c7eff(double s) const
        {
            // TODO: Neglecting contributions ~alpha_s / 4.0 / M_PI. These do involve spectator scattering,
            // cf. [BFS2001] Eq. (29), p. 8, and Eqs. (82)-(84), p. 30
            double lo = - 1.0/3.0 * c3 - 4.0/9.0 * c4 + 20.0/3.0 * c5 + 80.0/9.0 * c6;
            complex<double> nlo = (c1() - 6.0 * c2()) * A(s) - c8() * F87(s);

            return c7() + lo + (QCD::alpha_s(mu) / (4.0 * M_PI)) * nlo;
        }


        // cf. [GP2004], Eq. (55), p. 10
        complex<double> c9eff(const double & s) const
        {
            double c_0 = 4.0 / 3.0 * c1() + c2() + 5.5 * c3() - 2.0 / 3.0 * c4() + 52.0 * c5() + 32.0 / 3.0 * c6();
            double c_b = -0.5 * (7.0 * c3() + 4.0 / 3.0 * c4() + 76.0 * c5() + 64.0 / 3.0 * c6());
            double c = 2.0 / 9.0 * (6.0 * c3() + 32.0 * c5() + 32.0 / 3.0 * c6());

            // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
            complex<double> lo = c_b * h(s, m_b_MSbar) + c_0 * h0(s) + c;
            complex<double> nlo = c1() * (B(s) + 4.0 * C(s)) - 3.0 * c2() * (2.0 * B(s) - C(s)) - c8() * F89(s);

            return c9() + lo + (QCD::alpha_s(mu) / (4.0 * M_PI)) * nlo;
        }

        double kappa_1() const
        {
            // cf. [BHvD2010], Eq. (?), p. ?
            return (1.0 - 2.0 * QCD::alpha_s(mu) / (3.0 * M_PI) * std::log(mu / m_b_MSbar));
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
        complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstar_hat = m_Kstar / m_B;
            complex<double> wilson = (c9eff(s) - c9prime()) + h * (c10() - c10prime()) + kappa_1() * (c7eff(s) + c7prime()) * (2 * m_b_MSbar * m_B / s);
            complex<double> prefactor = complex<double>(0.0, -0.5 * m_B * m_B * m_B / m_Kstar / std::sqrt(s));
            double formfactor = lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s)) * form_factors->a_1(s_hat(s)) - (1 - s_hat(s)) * form_factors->a_2(s_hat(s));

            return prefactor * wilson * formfactor; // cf. [BHvD2010], Eq. (??)
        }

        complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstar_hat = m_Kstar / m_B;
            complex<double> wilson = (c9eff(s) + c9prime()) + h * (c10() + c10prime()) + kappa_1() * (c7eff(s) - c7prime()) * (2 * m_b_MSbar * m_B / s);
            complex<double> prefactor = complex<double>(0.0, std::sqrt(2 * lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s))) * m_B);

            return prefactor * wilson * form_factors->v(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }

        complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            complex<double> wilson = (c9eff(s) - c9prime()) + h * (c10() - c10prime()) + kappa_1() * (c7eff(s) - c7prime()) * (2 * m_b_MSbar * m_B / s);
            complex<double> prefactor = complex<double>(0.0, -std::sqrt(2) * m_B);

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

    complex<double>
    BToKstarDilepton<LowRecoil>::a_long(const Helicity & h, const double & s) const
    {
        return _imp->a_long(h, s);
    }

    complex<double>
    BToKstarDilepton<LowRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        return _imp->a_perp(h, s);
    }

    complex<double>
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
        return (norm(a_long(left_handed, s))
            + norm(a_long(right_handed, s))
            + norm(a_perp(left_handed, s))
            + norm(a_perp(right_handed, s))
            + norm(a_par(left_handed, s))
            + norm(a_par(right_handed, s)));
    }

    double
    BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 / differential_decay_width(s) * (real(a_par(left_handed, s) * conj(a_perp(left_handed, s))) - real(a_par(right_handed, s) * conj(a_perp(right_handed, s))));
    }

    double
    BToKstarDilepton<LowRecoil>::differential_unnormalized_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 * (real(a_par(left_handed, s) * conj(a_perp(left_handed, s))) - real(a_par(right_handed, s) * conj(a_perp(right_handed, s))));
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        double a = norm(a_perp(left_handed, s)) + norm(a_perp(right_handed, s));
        double b = norm(a_par(left_handed, s)) + norm(a_par(right_handed, s));

        return (a - b) / (a + b);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        double a = abs(a_long(left_handed, s) * conj(a_par(left_handed, s)) + conj(a_long(right_handed, s)) * a_par(right_handed, s));
        double b = std::sqrt((norm(a_long(left_handed, s)) + norm(a_long(right_handed, s))))
                * std::sqrt((norm(a_perp(left_handed, s)) + norm(a_perp(right_handed, s))));

        return a / b;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        double a = abs(a_long(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_long(right_handed, s)) * a_perp(right_handed, s));
        double b = abs(conj(a_long(left_handed, s)) * a_par(left_handed, s) + a_long(right_handed, s) * conj(a_par(right_handed, s)));

        return a / b;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return (norm(a_long(left_handed, s)) + norm(a_long(right_handed, s)))
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
        std::tr1::function<double (const double &)> num = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_unnormalized_forward_backward_asymmetry), this, std::tr1::placeholders::_1);

        std::tr1::function<double (const double &)> denom = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_decay_width), this, std::tr1::placeholders::_1);

        return integrate(num, 1000, s_min, s_max) / integrate(denom, 1000, s_min, s_max);
    }
}
