/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2020 Danny van Dyk
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

#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    using std::sin;
    using std::cos;
    using std::sqrt;

    // B -> P
    NonlocalFormFactor<nff::PToP>::~NonlocalFormFactor()
    {
    }

    complex<double>
    NonlocalFormFactor<nff::PToP>::jpsi_residues_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the J/psi residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalFormFactor<nff::PToP>::psi2s_residues_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the psi(2S) residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalFormFactor<nff::PToP>::moments_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the LCSR moments has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    // B -> V
    NonlocalFormFactor<nff::PToV>::~NonlocalFormFactor()
    {
    }

    complex<double>
    NonlocalFormFactor<nff::PToV>::jpsi_residues_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the J/psi residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalFormFactor<nff::PToV>::psi2s_residues_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the psi(2S) residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalFormFactor<nff::PToV>::moments_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the LCSR moments has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    namespace nff_utils
    {

        complex<double> z(const double & q2, complex<double> s_plus, complex<double> s_0)
        {
            return (pow(s_plus - q2, 0.5) - pow(s_plus - s_0, 0.5)) / (pow(s_plus - q2, 0.5) + pow(s_plus - s_0, 0.5));
        }

        // Blaschke factor capturing the two poles for J/psi and psi(2S).
        complex<double> blaschke_cc(const complex<double> z, const complex<double> z_Jpsi, const complex<double> z_psi2S)
        {
            return (z - z_Jpsi)/(1.0 - z * std::conj(z_Jpsi)) * (z - z_psi2S)/(1.0 - z * std::conj(z_psi2S));
        }

        //Expansion in z monomials (they form a basis on the unit circle)
        complex<double> P(complex<double> z, const complex<double> alpha[4])
        {
            return 1.0 / sqrt(2*M_PI) * (alpha[0] + z * (alpha[1] + z * (alpha[2] + z * alpha[3])));
        }

        //Expansion in polynomials orthogonal on the arc of the unit circle (zXY, zXY*)
        complex<double> PGvDV2020(complex<double> z, const complex<double> zXY, const complex<double> alpha[4])
        {

            const double alphaXY  = std::abs(std::arg(zXY));
            const double alphaXY2 = power_of<2>(alphaXY);
            const double alphaXY3 = power_of<3>(alphaXY);

            double denom = 2 * alphaXY2 + cos(2 * alphaXY) - 1;

            const complex<double> P0z = 1.0 / sqrt(2 * alphaXY);
            const complex<double> P1z = (z - sin(alphaXY) / alphaXY) * sqrt(alphaXY / denom);
            const complex<double> P2z = (z * z + z * sin(alphaXY) * (sin(2 * alphaXY) - 2 * alphaXY) / denom +
                                        2 * sin(alphaXY) * (sin(alphaXY) - alphaXY * cos(alphaXY)) / denom) *
                                        sqrt(2 * denom / (-9 * alphaXY + 8 * alphaXY3 + 8 * alphaXY * cos(2 * alphaXY) +
                                        alphaXY * cos(4 * alphaXY) + 4 * sin(2 * alphaXY) - 2 * sin(4 * alphaXY)));

            denom = -9 * alphaXY + 8 * alphaXY3 + 8 * alphaXY * cos(2 * alphaXY) + alphaXY * cos(4 * alphaXY) + 4 * sin(2 * alphaXY) - 2 * sin(4 * alphaXY);

            const complex<double> P3z_norm = 12.0 * sqrt(
                (2 * alphaXY - sin(2 * alphaXY)) * (2 * (-1 + alphaXY2 + cos(2 * alphaXY)) + alphaXY * sin(2 * alphaXY)) /
                ( 419 + 32 * alphaXY2 * (-65 + 36 * alphaXY2) + 16 * (-53 + 108 * alphaXY2) * cos(2 * alphaXY) +
                    4 * (151 + 72 * alphaXY2) * cos(4 * alphaXY) + cos(8 * alphaXY) + 16 * ((-11 + 4 * alphaXY2) * cos(6*alphaXY) +
                    96 * alphaXY * (5 * cos(alphaXY) + cos(3 * alphaXY)) * power_of<3>(sin(alphaXY))))
            );
            const complex<double> P3z_0 = - sin(alphaXY) * (15 + 8 * alphaXY2 + 16 * (-1 + alphaXY2) * cos(2 * alphaXY) + cos(4 * alphaXY) -
                24 * alphaXY * sin(2 * alphaXY)) / (3 * denom);
            const complex<double> P3z_1 = sin(alphaXY) * (-12 * alphaXY * cos(alphaXY) + 9 * sin(alphaXY) + sin(3 * alphaXY)) /
                (12 * (-1 + alphaXY2 + cos(2 * alphaXY)) + 6 * alphaXY * sin(2 * alphaXY));
            const complex<double> P3z_2 = sin(alphaXY) * (9 - 24 * alphaXY2 - 16 * cos(2 * alphaXY) + 7 * cos(4 * alphaXY) +
                4 * alphaXY * (4 * sin(2 * alphaXY) + sin(4 * alphaXY))) / (3 * denom);
            const complex<double> P3z = P3z_norm * (z * z * z + P3z_2 * z * z + P3z_1 * z + P3z_0);

            return alpha[0] * P0z + alpha[1] * P1z + alpha[2] * P2z + alpha[3] * P3z;
        }

    }

}
