/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <src/rare-b-decays/exclusive-b-to-s-gamma.hh>
#include <src/utils/options.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <cmath>

namespace eos
{
    template <>
    struct Implementation<BToKstarGamma>
    {
        Parameter abs_c7;

        Parameter abs_c7prime;

        Parameter ckm_A;

        Parameter ckm_lambda;

        Parameter ckm_etabar;

        Parameter ckm_rhobar;

        Implementation(const Parameters & p) :
            abs_c7(p["Abs{c7}"]),
            abs_c7prime(p["Abs{c7'}"]),
            ckm_A(p["CKM::A"]),
            ckm_lambda(p["CKM::lambda"]),
            ckm_etabar(p["CKM::etabar"]),
            ckm_rhobar(p["CKM::rhobar"])
        {
        }

        double beta() const
        {
            double A = ckm_A(), lambda = ckm_lambda(), etabar = ckm_etabar(), rhobar = ckm_rhobar();
            double A2 = A * A, lambda2 = lambda * lambda, lambda4 = lambda2 * lambda2;

            complex<double> a = complex<double>(rhobar, etabar);
            complex<double> b = 1.0 - A2 * lambda4 * conj(a);
            double c = std::sqrt((1.0 - A2 * lambda4) / (1.0 - lambda2));
            double d = std::pow(1.0 - A2 * lambda4 * rhobar, 2.0) + std::pow(A2 * lambda4 * etabar, 2.0);

            return -1.0 * arg(1.0 - a * b * (c / d));
        }

        double s_kstar_gamma() const
        {
            double r = abs_c7prime() / abs_c7();

            return -2.0 * r / (1.0 + r * r) * std::sin(2.0 * beta());
        }
    };

    BToKstarGamma::BToKstarGamma(const Parameters & parameters, const Options &) :
        PrivateImplementationPattern<BToKstarGamma>(new Implementation<BToKstarGamma>(parameters))
    {
    }

    BToKstarGamma::~BToKstarGamma()
    {
    }

    double
    BToKstarGamma::s_kstar_gamma() const
    {
        return _imp->s_kstar_gamma();
    }
}
