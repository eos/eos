/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Viktor Kuschke
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_CHARM_LOOPS_IMPL_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_CHARM_LOOPS_IMPL_HH 1

#include <eos/maths/complex.hh>
#include <eos/maths/interpolation.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>

#include <vector>

namespace eos
{
    namespace agv_2019a
    {
        // Often used constants

        static const double lnhalf  = - 0.693147180559945309;
        static const complex<double> lnmhalf = - 0.693147180559945309 + M_PI * 1.0i;
        static const complex<double> lnm1 = M_PI * 1.0i;
        static const complex<double> lnm2 = 0.693147180559945309 + M_PI * 1.0i;
        static const double ln2     = 0.693147180559945309;
        static const double ln4     = 2.0 * ln2;
        static const double ln64    = 6.0 * ln2;
        static const double ln256   = 8.0 * ln2;
        static const double ln2squ  = power_of<2>(ln2);
        static const double ln2cube = power_of<3>(ln2);
        static const constexpr double pisqu   = power_of<2>(M_PI);
        static const constexpr double li2half = 0.5822405264650124; // dilog(0.5)
        static const constexpr double li3half = 0.5372131936080402; // trilog(0.5)
        static const constexpr double zeta3   = 1.2020569031595943;
        static const double wx3 = 2.0 + std::sqrt(3.0);
        static const double wx4 = 2.0 - std::sqrt(3.0);
        static const double wx4squ = wx4 * wx4;
        static const double wx3inv = wx4;
        static const double wx4inv = wx3;

        struct CharmLoopsParameters
        {
            double muhat;

            complex<double> s_eps; // dilepton invariant mass divided by bottom quark mass m_b^2 including the epsilon prescription: q^2 / (m_b^2 - i * eps)
            complex<double> z_eps; // charm quark mass divided by bottom quark mass, squared, including the epsilon prescription: (m_c^2 - i * epsilon) / (m_b^2 - i * eps)

            // only quark mass dependend variables, cf. [AGV:2019A] p. 15 eq. (3.20)
            complex<double> xa;
            complex<double> xb;
            complex<double> xc;
            complex<double> xd;
            complex<double> xe;

            // q^2 and quark mass dependend variables, cf. [AGV:2019A] p. 15 eq. (3.20)
            complex<double> ya;
            complex<double> yb;
            complex<double> yc;
            complex<double> yd;
            complex<double> ye;

            /*!
            * Input parameters
            *
            * @param muhat              renormalization scale divided by bottom quark mass m_b
            * @param s                  dilepton invariant mass divided by squared bottom quark mass m_b
            * @param z                  squared charm quark mass divided by squared bottom quark mass
            * @param feynepsilonhat     epsilon prescription divided by bottom quark mass m_b squared
            *
            */

            CharmLoopsParameters(const double & muhat, const complex<double> & s, const double & z, const double & feynepsilonhat);
        };

        // Helper functions

        // Signum
        inline double my_sign(const double & x)
        {
            if (x < 0.0)
            {
                return - 1.0;
            }
            if (x > 0.0)
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        }

        // Heaviside Theta
        inline double my_HT(const double & x)
        {
            if (x == 0.0)
            {
                throw InternalError("Ill-defined Theta(0.0)");
            }
            if (x > 0.0)
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        }

        // Functions depending on the Heaviside step function

        // Triangle function T(a,b;x) from [FTW:2016A] p. 7 eq. (3.3)
        inline double T(const complex<double> & a, const complex<double> & b, const complex<double> & x)
        {
            const complex<double> amb = a - b;
            const complex<double> xconj = std::conj(x);
            const complex<double> aconj = std::conj(a);
            const double denom = (xconj * amb).imag();

            if (denom == 0.0)
            {
                return 0.0;
            }
            else
            {
                const double arg1 = (xconj * a).imag() / denom;
                const double arg2 = 1.0 - (xconj * a).imag() / denom;
                const double arg3 = - 1.0 - (aconj * b).imag() / denom; // note the minus sign in front of (aconj * b).imag() / denom: in [FTW:2016A] it is missing which seems to be a typo

                if ( arg1 < 0.0 || arg2 < 0.0 || arg3 < 0.0)
                {
                    return 0.0;
                }
                if (arg1 == 0.0 || arg2 == 0.0 || arg3 == 0.0)
                {
                    throw InternalError("Ill-defined Theta(0.0)");
                }
                else
                {
                    return 1.0;
                }
            }
        }

        inline complex<double> p(const complex<double> & x1, const complex<double> & x2)
        {
            const complex<double> x1conj = std::conj(x1);
            const double denom = (x1 - x2 + x2 * x1conj).imag();

            if (denom == 0.0)
            {
                throw InternalError("0 in denominator");
            }
            else
            {
                return (x1 - 1.0) * x2.imag() / denom;
            }
        }

        // r function from [FTW:2016A] p. 11 eq. (4.5)
        inline double r(const complex<double> & a, const complex<double> & b)
        {
            const complex<double> aconjb = std::conj(a) * b;
            const double denom = aconjb.imag();
            const double num = power_of<2>(std::abs(a)) * b.imag() - power_of<2>(std::abs(b)) * a.imag();

            if (denom == 0.0)
            {
                throw InternalError("0 in denominator");
            }
            else
            {
                return num / denom;
            }
        }

        // H1 and H2 function from [FTW:2016A] p. 11 eq. (4.5) and p. 12 eq. (4.9), respectively

        inline double H1(const complex<double> & a, const complex<double> & b)
        {
            const complex<double> aconjb = std::conj(a) * b;

            if (aconjb.imag() == 0.0)
            {
                return 0.0;
            }
            else
            {
                double min;

                if (1.0 <= power_of<2>(std::abs(a)) * b.imag() / aconjb.imag())
                {
                    min = 1.0;
                }
                else
                {
                   min = power_of<2>(std::abs(a)) * b.imag() / aconjb.imag();
                }

                const double r_val = agv_2019a::r(a,b);
                const double arg = min - r_val;

                if (r_val < 0.0 || arg < 0.0)
                {
                    return 0.0;
                }
                if (r_val == 0.0 || arg == 0.0)
                {
                    throw InternalError("Ill-defined Theta(0.0)");
                }
                else
                {
                    return 1.0;
                }
            }
        }

        inline double H2(const complex<double> & a, const complex<double> & b)
        {
            const complex<double> aconjb = std::conj(a) * b;

            if (aconjb.imag() == 0.0)
            {
                return 0.0;
            }
            else
            {
                const double r_val = agv_2019a::r(a,b);
                const double imaimb = a.imag() * b.imag();

                if (r_val < 0.0 || r_val > 1.0 || imaimb > 0.0)
                {
                    return 0.0;
                }
                if (r_val == 0.0 || r_val == 1.0 || imaimb == 0.0)
                {
                    throw InternalError("Ill-defined Theta(0.0)");
                }
                else
                {
                    return 1.0;
                }
            }
        }

        // LO functions
        complex<double> f190(const CharmLoopsParameters & clp);
        complex<double> f290(const CharmLoopsParameters & clp);

        // Counterterms
        complex<double> f17ctQs(const CharmLoopsParameters & );
        complex<double> f17ctQc(const CharmLoopsParameters & );
        complex<double> f17ctQb(const CharmLoopsParameters & clp);
        complex<double> f19ctQs(const CharmLoopsParameters & clp);
        complex<double> f19ctQc(const CharmLoopsParameters & clp);
        complex<double> f19ctQb(const CharmLoopsParameters & clp);
        complex<double> f27ctQs(const CharmLoopsParameters & );
        complex<double> f27ctQc(const CharmLoopsParameters & );
        complex<double> f27ctQb(const CharmLoopsParameters & clp);
        complex<double> f29ctQs(const CharmLoopsParameters & clp);
        complex<double> f29ctQc(const CharmLoopsParameters & clp);
        complex<double> f29ctQb(const CharmLoopsParameters & clp);

        // Two-loop functions
        complex<double> f17a(const CharmLoopsParameters & clp);
        complex<double> f19a(const CharmLoopsParameters & clp);
        complex<double> f27a(const CharmLoopsParameters & clp);
        complex<double> f29a(const CharmLoopsParameters & clp);

        complex<double> f17b(const CharmLoopsParameters & clp);
        complex<double> f19b(const CharmLoopsParameters & clp);
        complex<double> f27b(const CharmLoopsParameters & clp);
        complex<double> f29b(const CharmLoopsParameters & clp);

        complex<double> f17c(const CharmLoopsParameters & clp);
        complex<double> f19c(const CharmLoopsParameters & clp);
        complex<double> f27c(const CharmLoopsParameters & clp);
        complex<double> f29c(const CharmLoopsParameters & clp);

        complex<double> weight4_wx3_wx4(const CharmLoopsParameters & clp, const complex<double> & wx);
        complex<double> weight4_w4_w5_w7(const CharmLoopsParameters & clp, const complex<double> & w);
        complex<double> GPLweight4Parts(const CharmLoopsParameters & clp);

        complex<double> f17d(const CharmLoopsParameters & clp);
        complex<double> f19d(const CharmLoopsParameters & clp);
        complex<double> f27d(const CharmLoopsParameters & clp);
        complex<double> f29d(const CharmLoopsParameters & clp);

        complex<double> f17e(const CharmLoopsParameters & );
        complex<double> f19e(const CharmLoopsParameters & clp);
        complex<double> f27e(const CharmLoopsParameters & );
        complex<double> f29e(const CharmLoopsParameters & clp);

        complex<double> F17_Qc(const CharmLoopsParameters & clp);
        complex<double> F19_Qc(const CharmLoopsParameters & clp);
        complex<double> F27_Qc(const CharmLoopsParameters & clp);
        complex<double> F29_Qc(const CharmLoopsParameters & clp);
        complex<double> F17_Qsb(const CharmLoopsParameters & clp);
        complex<double> F19_Qsb(const CharmLoopsParameters & clp);
        complex<double> F27_Qsb(const CharmLoopsParameters & clp);
        complex<double> F29_Qsb(const CharmLoopsParameters & clp);

        // Helper functions
        complex<double> F17_Qc(const complex<double> & s, const double & mu, const double & m_c, const double & m_b);
        complex<double> F19_Qc(const complex<double> & s, const double & mu, const double & m_c, const double & m_b);
        complex<double> F27_Qc(const complex<double> & s, const double & mu, const double & m_c, const double & m_b);
        complex<double> F29_Qc(const complex<double> & s, const double & mu, const double & m_c, const double & m_b);
        complex<double> F17_Qsb(const complex<double> & s, const double & mu, const double & m_c, const double & m_b);
        complex<double> F19_Qsb(const complex<double> & s, const double & mu, const double & m_c, const double & m_b);
        complex<double> F27_Qsb(const complex<double> & s, const double & mu, const double & m_c, const double & m_b);
        complex<double> F29_Qsb(const complex<double> & s, const double & mu, const double & m_c, const double & m_b);

        complex<double> delta_c7_Qc(const complex<double> & s, const double & mu, const double & alpha_s, const double & m_c, const double & m_b, const WilsonCoefficients<BToS> & wc, bool use_nlo = true);
        complex<double> delta_c7(const complex<double> & s, const double & mu, const double & alpha_s, const double & m_c, const double & m_b, const WilsonCoefficients<BToS> & wc, bool use_nlo = true);
        complex<double> delta_c9_Qc(const complex<double> & s, const double & mu, const double & alpha_s, const double & m_c, const double & m_b, const WilsonCoefficients<BToS> & wc, bool use_nlo = true);
        complex<double> delta_c9(const complex<double> & s, const double & mu, const double & alpha_s, const double & m_c, const double & m_b, const WilsonCoefficients<BToS> & wc, bool use_nlo = true);
    }
}

#endif
