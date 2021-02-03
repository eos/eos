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
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

using std::sin;
using std::cos;
using std::sqrt;


namespace eos
{
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
        complex<double> P(complex<double> z, const complex<double> & alpha_0, const complex<double> & alpha_1, const complex<double> & alpha_2)
        {
            return 1.0 / sqrt(2*M_PI) * (alpha_0 + alpha_1*z + alpha_2*z*z);
        }

        //Expansion in polynomials orthogonal on the arc of the unit circle (zXY, zXY*)
        complex<double> PGvDV2020(complex<double> z, const complex<double> zXY,
            const complex<double> & alpha_0, const complex<double> & alpha_1, const complex<double> & alpha_2)
        {

            const double alphaXY = std::abs(std::arg(zXY));

            const double denom = 2*pow(alphaXY, 2) + cos(2*alphaXY) - 1;

            const complex<double> P0z = 1.0/sqrt(2*alphaXY);
            const complex<double> P1z = (z - sin(alphaXY)/alphaXY)*sqrt(alphaXY/denom);
            const complex<double> P2z = ( z*z + z*sin(alphaXY)*(sin(2*alphaXY)-2*alphaXY)/denom +
                                        2*sin(alphaXY)*(sin(alphaXY)-alphaXY*cos(alphaXY))/denom) *
                                        sqrt( 2*denom/(-9*alphaXY + 8*pow(alphaXY,3) + 8*alphaXY*cos(2*alphaXY) +
                                        alphaXY*cos(4*alphaXY) + 4*sin(2*alphaXY) - 2*sin(4*alphaXY)) );

            return alpha_0*P0z + alpha_1*P1z + alpha_2*P2z;
        }

    }

}
