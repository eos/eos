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

#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>
#include <eos/maths/power-of.hh>
#include <eos/maths/szego-polynomial.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <array>

namespace eos
{
    using std::sin;
    using std::cos;
    using std::sqrt;

    // B -> P
    NonlocalFormFactor<PToP>::~NonlocalFormFactor()
    {
    }

    complex<double>
    NonlocalFormFactor<PToP>::jpsi_residues_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the J/psi residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalFormFactor<PToP>::psi2s_residues_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the psi(2S) residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalFormFactor<PToP>::moments_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the LCSR moments has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    // B -> V
    NonlocalFormFactor<PToV>::~NonlocalFormFactor()
    {
    }

    complex<double>
    NonlocalFormFactor<PToV>::jpsi_residues_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the J/psi residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalFormFactor<PToV>::psi2s_residues_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the psi(2S) residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalFormFactor<PToV>::moments_not_implemented() const
    {
        throw InternalError("A NonlocalFormFactor without implementation of the LCSR moments has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    namespace nff_utils
    {

        complex<double> z(const complex<double> & q2, complex<double> s_plus, complex<double> s_0)
        {
            return (pow(s_plus - q2, 0.5) - pow(s_plus - s_0, 0.5)) / (pow(s_plus - q2, 0.5) + pow(s_plus - s_0, 0.5));
        }

        complex<double> z(const double & q2, complex<double> s_plus, complex<double> s_0)
        {
            return z(complex<double>(q2, 0.0), s_plus, s_0);
        }

        // Blaschke factor capturing the two poles for J/psi and psi(2S).
        complex<double> blaschke_cc(const complex<double> & z, const complex<double> & z_Jpsi, const complex<double> & z_psi2S)
        {
            return (z - z_Jpsi)/(1.0 - z * std::conj(z_Jpsi)) * (z - z_psi2S)/(1.0 - z * std::conj(z_psi2S));
        }
    }

    std::shared_ptr<SzegoPolynomial<5u>>
    PolynomialsFactory::create(const std::string & opt_q)
    {
        switch (opt_q[0])
        {
        case 's':
            // These values are computed using t_0 = 4 GeV^2, m_Bs = 5.366 GeV and m_phi = 1.020 GeV
            return std::make_shared<SzegoPolynomial<5u>>(SzegoPolynomial<5u>::FlatMeasure(2.18309));
            break;
        default: // opt_q = u, d
            // These values are computed using t_0 = 4 GeV^2, m_B = 5.279 GeV and m_K* = 0.896 GeV
            return std::make_shared<SzegoPolynomial<5u>>(SzegoPolynomial<5u>::FlatMeasure(2.27631));
            break;
        }
    }
}
