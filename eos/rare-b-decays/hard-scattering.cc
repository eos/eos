/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#include <eos/rare-b-decays/hard-scattering.hh>

#include <gsl/gsl_sf_dilog.h>

namespace eos
{
    complex<double>
    HardScattering::I1(const double & s, const double & u, const double & m_q, const double & m_B)
    {
        if (m_q == 0.0)
            return complex<double>(1.0, 0.0);

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
            gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
            dilog1 = complex<double>( res_re.val, res_im.val);

            dilogArg = complex<double>((a2 - 1.0)/(a2 + 1.0), +2.0 * a / (a2 + 1.0));
            gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
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
            gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
            dilog1 = complex<double>(res_re.val, res_im.val);

            dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), +2.0 * a / (a2 + 1.0));
            gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
            dilog2 = complex<double>(res_re.val, res_im.val);

            LypLym = -1.0 / 3.0 * M_PI * M_PI - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                + dilog1 + dilog2;
        }

        return 1.0 + 2.0 * m_q2 / (ubar * (m_B2 - s)) * (LxpLxm - LypLym);
    }
}
