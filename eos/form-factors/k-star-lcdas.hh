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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_K_STAR_LCDAS_HH
#define EOS_GUARD_EOS_FORM_FACTORS_K_STAR_LCDAS_HH 1

#include <eos/form-factors/vec-lcdas.hh>
#include <eos/maths/gegenbauer-polynomial.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    class AntiKStarLCDAs :
        public VectorLCDAs,
        public PrivateImplementationPattern<AntiKStarLCDAs>
    {
        public:
            AntiKStarLCDAs(const Parameters &, const Options &);
            ~AntiKStarLCDAs();

            static VectorLCDAs * make(const Parameters &, const Options &);

            /* Twist 2 LCDA (even) para Gegenbauer coefficients */
            double a1para(const double & mu) const override;
            double a2para(const double & mu) const override;
            double a3para(const double & mu) const override;
            double a4para(const double & mu) const override;
            double fpara() const override;

            /* Twist 2 LCDA (even) perp Gegenbauer coefficients */
            double a1perp(const double & mu) const override;
            double a2perp(const double & mu) const override;
            double a3perp(const double & mu) const override;
            double a4perp(const double & mu) const override;
            double fperp(const double & mu) const override;

            /* Twist 3 LCDA parameters */
            double kappa3para(const double & mu) const override;
            double omega3para(const double & mu) const override;
            double lambda3para(const double & mu) const override;
            double zeta3para(const double & mu) const override;
            double lambda3paratilde(const double & mu) const override;
            double omega3paratilde(const double & mu) const override;
            double kappa3perp(const double & mu) const override;
            double omega3perp(const double & mu) const override;
            double lambda3perp(const double & mu) const override;

            /* Twist 4 LCDA parameters */
            double zeta4para(const double & mu) const override;
            double omega4paratilde(const double & mu) const override;
            double zeta4perp(const double & mu) const override;
            double zeta4perptilde(const double & mu) const override;
            double kappa4para(const double & mu) const override;
            double kappa4perp(const double & mu) const override;

            /* Twist 2 LCDAs */
            double phipara(const double & u, const double & mu) const override;
            double phiperp(const double & u, const double & mu) const override;

            /* Twist 3 two particle LCDAs */
            virtual double phi3para(const double & u, const double & mu) const override;
            virtual double phi3perp(const double & u, const double & mu) const override;
            virtual double psi3para(const double & u, const double & mu) const override;
            virtual double psi3perp(const double & u, const double & mu) const override;

            /* Twist 3 three particle LCDAs */
            double Phi3para(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Phi3paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Phi3perp(const double & u1, const double & u2, const double & u3, const double & mu) const override;

            /* Twist 4 three particle chiral even LCDAs */
            double Psi4para(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Psi4paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Phi4para(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Phi4paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Xi4para(const double & u1, const double & u2, const double & u3, const double & mu) const override;

            /* Twist 4 three particle chiral odd LCDAs */
            double Psi4perp(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Psi4perptilde(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perp1(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perp2(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perp3(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perp4(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perptilde1(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perptilde2(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perptilde3(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perptilde4(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Xi4perp(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }

            /* Internal diagnostics */
            Diagnostics diagnostics() const;
    };

    class KStarLCDAs :
        public VectorLCDAs,
        public PrivateImplementationPattern<KStarLCDAs>
    {
        public:
            KStarLCDAs(const Parameters &, const Options &);
            ~KStarLCDAs();

            static VectorLCDAs * make(const Parameters &, const Options &);

            /* Twist 2 LCDA (even) para Gegenbauer coefficients */
            double a1para(const double & mu) const override;
            double a2para(const double & mu) const override;
            double a3para(const double & mu) const override;
            double a4para(const double & mu) const override;
            double fpara() const override;

            /* Twist 2 LCDA (even) perp Gegenbauer coefficients */
            double a1perp(const double & mu) const override;
            double a2perp(const double & mu) const override;
            double a3perp(const double & mu) const override;
            double a4perp(const double & mu) const override;
            double fperp(const double & mu) const override;

            /* Twist 3 LCDA parameters */
            double kappa3para(const double & mu) const override;
            double omega3para(const double & mu) const override;
            double lambda3para(const double & mu) const override;
            double zeta3para(const double & mu) const override;
            double lambda3paratilde(const double & mu) const override;
            double omega3paratilde(const double & mu) const override;
            double kappa3perp(const double & mu) const override;
            double omega3perp(const double & mu) const override;
            double lambda3perp(const double & mu) const override;

            /* Twist 4 LCDA parameters */
            double zeta4para(const double & mu) const override;
            double omega4paratilde(const double & mu) const override;
            double zeta4perp(const double & mu) const override;
            double zeta4perptilde(const double & mu) const override;
            double kappa4para(const double & mu) const override;
            double kappa4perp(const double & mu) const override;

            /* Twist 2 LCDAs */
            double phipara(const double & u, const double & mu) const override;
            double phiperp(const double & u, const double & mu) const override;

            /* Twist 3 two particle LCDAs */
            virtual double phi3para(const double & u, const double & mu) const override;
            virtual double phi3perp(const double & u, const double & mu) const override;
            virtual double psi3para(const double & u, const double & mu) const override;
            virtual double psi3perp(const double & u, const double & mu) const override;

            /* Twist 3 three particle LCDAs */
            double Phi3para(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Phi3paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Phi3perp(const double & u1, const double & u2, const double & u3, const double & mu) const override;

            /* Twist 4 three particle chiral even LCDAs */
            double Psi4para(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Psi4paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Phi4para(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Phi4paratilde(const double & u1, const double & u2, const double & u3, const double & mu) const override;
            double Xi4para(const double & u1, const double & u2, const double & u3, const double & mu) const override;

            /* Twist 4 three particle chiral odd LCDAs */
            double Psi4perp(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Psi4perptilde(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perp1(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perp2(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perp3(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perp4(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perptilde1(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perptilde2(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perptilde3(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Phi4perptilde4(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }
            double Xi4perp(const double & /*u1*/, const double & /*u2*/, const double & /*u3*/, const double & /*mu*/) const override { return 0.0; }

            /* Internal diagnostics */
            Diagnostics diagnostics() const;
    };
}

#endif
