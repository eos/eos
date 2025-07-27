/*
 * Copyright (c) 2019 Stephan Kürten
 * Copyright (c) 2019 Danny van Dyk
 * Copyright (c) 2021-2023 Méril Reboud
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

#ifndef EOS_GUARD_EOS_UTILS_KMATRIX_HH
#define EOS_GUARD_EOS_UTILS_KMATRIX_HH 1

#include <eos/maths/complex.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/parameters.hh>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <array>
#include <memory>
#include <vector>

namespace eos
{
    template <unsigned nchannels_, unsigned nresonances_> class KMatrix
    {
        public:
            struct Channel;
            struct Resonance;

            std::array<std::shared_ptr<KMatrix::Channel>, nchannels_>     _channels;
            std::array<std::shared_ptr<KMatrix::Resonance>, nresonances_> _resonances;
            // The non-resonant contribution is described with a set of constants
            std::array<std::array<Parameter, nchannels_>, nchannels_>     _bkgcst;

            const std::string & _prefix;

            // Khat contains the normalized K-matrix entries
            // as described in [CBHKSS:1995A].
            gsl_matrix_complex * _Khat;
            // That is defined as T = (1-i*rho*K)^(-1)*K
            gsl_matrix_complex * _That;

            // Two temporary matrices and a permutation are needed to evaluate That
            gsl_matrix_complex * _tmp_1;
            gsl_matrix_complex * _tmp_2;
            gsl_matrix_complex * _tmp_3;
            gsl_permutation *    _perm;

            // Constructor
            KMatrix(std::array<std::shared_ptr<KMatrix::Channel>, nchannels_> channels, std::array<std::shared_ptr<KMatrix::Resonance>, nresonances_> resonances,
                    std::array<std::array<Parameter, nchannels_>, nchannels_> bkgcst, const std::string & prefix);

            // Destructor
            ~KMatrix();

            // Adapt s to avoid resonnances masses
            complex<double> adapt_s(const complex<double> s) const;

            // Return the rowindex^th row of the T matrix on the first Riemann sheet, i.e. defined as T = n * (1 - i * K * rho * n * n)^(-1) * K * n
            // rowindex corresponds to the initial channel
            // If second_sheet is set to true, the value is calculated on the second Riemann sheet
            std::array<complex<double>, nchannels_> tmatrix_row(unsigned rowindex, const complex<double> s, bool second_sheet = false) const;


            // Return the K matrix partial and total widths of a resonance.
            // Note that these widths do not necessarily correspond to the experimental ones.
            double partial_width(unsigned resonance, unsigned channel) const;

            double width(unsigned resonance) const;

            // Spectral function of a resonance.
            double spectral_function(unsigned resonance, const double & s) const;
    };

    template <unsigned nchannels_, unsigned nresonances_> struct KMatrix<nchannels_, nresonances_>::Channel
    {
            std::string _name;
            // Masses of the two final state particles
            Parameter   _m1;
            Parameter   _m2;

            // Properties of the centrifugal barrier factors, cf (50.26)
            unsigned  _l_orbital;
            Parameter _q0;

            std::array<Parameter, nresonances_> _g0s;

            Channel(std::string name, Parameter m1, Parameter m2, unsigned l_orbital, Parameter q0, std::array<Parameter, nresonances_> g0s) :
                _name(name),
                _m1(m1),
                _m2(m2),
                _l_orbital(l_orbital),
                _q0(q0),
                _g0s(g0s)
            {
                if (m1.evaluate() < 0 || m2.evaluate() < 0)
                {
                    throw InternalError("Masses of channel '" + _name + "' must be positive");
                }
                if (q0.evaluate() < 0)
                {
                    throw InternalError("Scale parameter of the channel'" + _name + "' must be positive");
                }
            }

            // Phase space factor
            virtual complex<double> rho(const complex<double> & s) = 0;

            // Analytic continuation of the phase space factor
            virtual complex<double> chew_mandelstam(const complex<double> & s) = 0;
    };

    template <unsigned nchannels_, unsigned nresonances_> struct KMatrix<nchannels_, nresonances_>::Resonance
    {
            std::string _name;

            // Mass of the resonance
            Parameter _m;

            Resonance(std::string name, Parameter m) :
                _name(name),
                _m(m)
            {
                if (m.evaluate() < 0)
                {
                    throw InternalError("Mass of resonance '" + _name + "' must be positive");
                }
            }
    };

    namespace kmatrix_utils
    {
        // We follow the notation of the PDG's resonance review, i.e. z = (q / q0)
        complex<double> blatt_weisskopf_factor(const unsigned & l, const complex<double> & z);
    } // namespace kmatrix_utils
} // namespace eos

#endif
