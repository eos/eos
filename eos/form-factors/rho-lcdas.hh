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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_RHO_LCDAS_HH
#define EOS_GUARD_EOS_FORM_FACTORS_RHO_LCDAS_HH 1

#include <eos/form-factors/vec-lcdas.hh>
#include <eos/maths/gegenbauer-polynomial.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    class RhoLCDAs :
        public VectorLCDAs,
        public PrivateImplementationPattern<RhoLCDAs>
    {
        public:
            RhoLCDAs(const Parameters &, const Options &);
            ~RhoLCDAs();

            static VectorLCDAs * make(const Parameters &, const Options &);

            /* Gegenbauer polynomials */
            const GegenbauerPolynomial gp_2_3o2;
            const GegenbauerPolynomial gp_4_3o2;

            /* Twist 2 LCDA (even) para Gegenbauer coefficients */
            double a1para(const double & /*mu*/) const override { return 0.0; }
            double a2para(const double & mu) const override;
            double a3para(const double & /*mu*/) const override { return 0.0; }
            double a4para(const double & mu) const override;
            double fpara() const override;

            /* Twist 2 LCDA (even) perp Gegenbauer coefficients */
            double a1perp(const double & /*mu*/) const override { return 0.0; }
            double a2perp(const double & mu) const override;
            double a3perp(const double & /*mu*/) const override { return 0.0; }
            double a4perp(const double & mu) const override;
            double fperp(const double & mu) const override;

            /* Twist 2 LCDAs */
            double phipara(const double & u, const double & mu) const override;
            double phiperp(const double & u, const double & mu) const override;

            /* Internal diagnostics */
            Diagnostics diagnostics() const;
    };
}

#endif
