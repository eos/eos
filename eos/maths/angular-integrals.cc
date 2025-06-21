/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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

#include <eos/maths/angular-integrals.hh>
#include <eos/maths/power-of.hh>

#include <algorithm>
#include <cmath>

namespace eos
{
    double wigner_3j(const int & j1, const int & m1, const int & j2, const int & m2, const int & j3, const int & m3)
    {
        if (m1 + m2 + m3 != 0 || j1 < 0 || j2 < 0 || j3 < 0) return 0.0;

        double ret = 0.0;
        const double pref = (((j1 - j2 - m3) % 2) ? -1.0 : 1.0)
                    * std::sqrt(std::tgamma(j1 + j2 - j3 + 1)) * std::sqrt(std::tgamma(j1 - j2 + j3 + 1))
                    * std::sqrt(std::tgamma(-j1 + j2 + j3 + 1)) / std::sqrt(std::tgamma(1 + j1 + j2 + j3 + 1))
                    * std::sqrt(std::tgamma(j1 + m1 + 1)) * std::sqrt(std::tgamma(j1 - m1 + 1)) * std::sqrt(std::tgamma(j2 + m2 + 1))
                    * std::sqrt(std::tgamma(j2 - m2 + 1)) * std::sqrt(std::tgamma(j3 + m3 + 1)) * std::sqrt(std::tgamma(j3 - m3 + 1));
        const int lowerlim = std::max({ 0, j2 - j3 - m1, j1 - j3 + m2 });
        const int upperlim = std::min({ j1 + j2 - j3, j1 - m1, j2 + m2});

        for(int i = lowerlim; i <= upperlim; i++)
        {
            ret += ((i % 2) ? -1.0 : 1.0) / (std::tgamma(i + 1) * std::tgamma(j1 + j2 - j3 - i + 1) * std::tgamma(j1 - m1 - i + 1) * std::tgamma(j2 + m2 - i + 1) * std::tgamma(j3 - j2 + m1 + i + 1) * std::tgamma(j3 - j1 - m2 + i + 1));
        }
        return ret * pref;
    }

    // Integral over 1 associated Legendre polynomial
    double legendre_integral(const int & l, const int & m)
    {
        const int absm = (m >= 0) ? m : -m;
        if ( l == 0 && absm == 0 ) return 2.0;
        if ( l > 0  && l < absm ) return 0.0;
        if ( l > 0 && absm >= 0 && (l % 2) != (absm % 2)) return 0.0;

        const double reflection = (m >= 0) ? 1.0 : std::tgamma(l + m + 1) / std::tgamma(l - m + 1) * ((absm % 2) ? -1.0 : 1.0);
        const double sign = (l % 2) ? -1.0 : 1.0;
        return reflection * sign * absm * std::pow(2.0, absm - 1) * std::tgamma(l / 2.0) * std::tgamma((l + absm + 1) / 2.0) / std::tgamma((l + 3) / 2.0) / std::tgamma((l - absm) / 2 + 1);
    }

    // Integral over 2 associated Legendre polynomials
    double two_legendre_integral(const int & l1, const int & m1, const int & l2, const int & m2)
    {
        if ((m1 > l1 && m1 > 0) || (m2 > l2 && m2 > 0) || (-m1 > l1 && m1 < 0) || (-m2 > l2 && m2 < 0)) return 0.0;
        if ((m1 == m2) && (l1 != l2)) return 0.0;
        if ((m1 == m2) && (l1 == l2)) return 2.0 * std::tgamma(l1 + m1 + 1) / std::tgamma(l1 - m1 + 1) / (2 * l1 + 1);
        if ((l2 == 0) && (m2 == 0)) return legendre_integral(l1, m1);
        if ((l1 == 0) && (m1 == 0)) return legendre_integral(l2, m2);

        double ret = 0.0;
        const int m12 = m1 + m2;
        const double pref = ((m12 % 2) ? -1.0 : 1.0) * std::sqrt(std::tgamma(l1 + m1 + 1) * std::tgamma(l2 + m2 + 1) / (std::tgamma(l1 - m1 + 1) * std::tgamma(l2 - m2 + 1)));

        const int lowerlim = std::max({ m12, l1 - l2, l2 - l1 });
        for (int l12 = lowerlim; l12 <= l1 + l2; l12++)
        {
            ret += (2 * l12 + 1) * legendre_integral(l12, m12) * wigner_3j(l1, 0, l2, 0, l12, 0) * wigner_3j(l1, m1, l2, m2, l12, -m12) * std::sqrt(std::tgamma(l12 - m12 + 1) / std::tgamma(l12 + m12 + 1));
        }
        return pref * ret;
    }

    // Integral over 3 associated Legendre polynomials
    double three_legendre_integral(const int & l1, const int & m1, const int & l2, const int & m2, const int & l3, const int & m3)
    {
        if ((m1 > l1 && m1 > 0) || (m2 > l2 && m2 > 0) || (m3 > l3 && m3 > 0) || (-m1 > l1 && m1 < 0) || (-m2 > l2 && m2 < 0) || (-m3 > l3 && m3 < 0)) return 0.0;
        if ((l3 == 0) && (m3 == 0)) return two_legendre_integral(l1,m1,l2,m2);
        if ((l2 == 0) && (m2 == 0)) return two_legendre_integral(l1,m1,l3,m3);
        if ((l1 == 0) && (m1 == 0)) return two_legendre_integral(l2,m2,l3,m3);

        double ret = 0.0;
        const int m12 = m1 + m2;
        const double pref = ((m12 % 2) ? -1.0 : 1.0) * std::sqrt(std::tgamma(l1 + m1 + 1) * std::tgamma(l2 + m2 + 1) / (std::tgamma(l1 - m1 + 1) * std::tgamma(l2 - m2 + 1)));

        const int lowerlim = std::max({ m12, l1 - l2, l2 - l1 });
        for (int l12 = lowerlim; l12 <= l1 + l2; l12++)
        {
            ret += (2 * l12 + 1) * two_legendre_integral(l12, m12, l3, m3) * wigner_3j(l1, 0, l2, 0, l12, 0) * wigner_3j(l1, m1, l2, m2, l12, -m12)
                * std::sqrt(std::tgamma(l12 - m12 + 1) / std::tgamma(l12 + m12 + 1));
        }
        return pref * ret;
    }
}
