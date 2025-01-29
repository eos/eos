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

#ifndef EOS_GUARD_EOS_MODELS_TOP_LOOPS_HH
#define EOS_GUARD_EOS_MODELS_TOP_LOOPS_HH 1

namespace eos
{
    /*!
     * Helper class to bundle top-loop functions in the SM which are relevant to @f$b\to s@f$ transitions.
     */
    struct TopLoops
    {
            ///@{
            /*!
             * One-loop functions in the MSbar scheme at the renormalization scale @f$\mu_t@f$.
             *
             * @param x_t        @f$= m_t(\mu_t)^2 / m_W^2@f$
             */
            static double A0(const double & x_t);

            static double B0(const double & x_t);

            static double C0(const double & x_t);

            static double D0(const double & x_t);

            static double E0(const double & x_t);

            static double F0(const double & x_t);
            ///@}

            ///@{
            /*!
             * Two-loop functions in the MSbar scheme at the renormalization scale @f$\mu_t@f$.
             *
             * @param x_t        @f$= m_t(\mu_t)^2 / m_W^2@f$
             * @param log_t      @f$= \ln(\mu_t / m_t(\mu_t))@f$
             */
            static double A1(const double & x_t, const double & log_t);

            static double B1(const double & x_t, const double & log_t);

            static double C1(const double & x_t, const double & log_t);

            static double D1(const double & x_t, const double & log_t);

            static double E1(const double & x_t, const double & log_t);

            static double F1(const double & x_t, const double & log_t);

            static double G1(const double & x_t, const double & log_t);
            ///@}
    };
} // namespace eos

#endif
