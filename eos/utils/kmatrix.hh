/*
 * Copyright (c) 2019 Stephan Kürten
 * Copyright (c) 2019 Danny van Dyk
 * Copyright (c) 2021 Méril Reboud
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


#include <eos/utils/complex.hh>
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
    template <unsigned nchannels_, unsigned nresonances_>
    class KMatrix
    {
        public:

            struct Channel;
            struct Resonance;

            std::vector<std::shared_ptr<KMatrix::Channel>> _channels;
            std::vector<std::shared_ptr<KMatrix::Resonance>> _resonances;
            std::vector<std::vector<Parameter>> _bkgcst;

            const std::string & _prefix;

            // Khat contains the normalized K-matrix entries
            // as described in [CBHKSS:1995A].
            gsl_matrix_complex * _Khat;
            // That is defined as T = (1-i*rho*K)^(-1)*K
            gsl_matrix_complex * _That;

            // Two temporary matrices and a permutation are needed to evaluate That
            gsl_matrix_complex * _tmp_1;
            gsl_matrix_complex * _tmp_2;
            gsl_permutation * _perm;

            // Constructor
            KMatrix(std::initializer_list<std::shared_ptr<KMatrix::Channel>> channels,
                std::initializer_list<std::shared_ptr<KMatrix::Resonance>> resonances,
                std::vector<std::vector<Parameter>> bkgcst,
                const std::string & prefix);

            // Destuctor
            ~KMatrix();

            // Adapt s to avoid resonnances masses
            double adapt_s(const double s) const;

            // Return rowindex^th row of the T matrix defined as T = (1-i*rho*K)^(-1)*K
            // rowindex corresponds to the initial channel
            std::array<complex<double>, nchannels_> tmatrix_row(unsigned rowindex, const double s) const;

            // Return the K matrix partial and total widths of a resonance.
            // Note that these widths do not necessarily correspond to the experimental ones.
            double partial_width(unsigned resonance, unsigned channel) const;

            double width(unsigned resonance) const;
    };


    template <unsigned nchannels_, unsigned nresonances_>
    struct KMatrix<nchannels_, nresonances_>::Channel
    {
        std::string _name;
        //Masses of the two final state particles
        double _m1;
        double _m2;

        unsigned _N_orbital;

        std::vector<Parameter> _g0s;

        Channel(std::string name, double m1, double m2, unsigned N_orbital, std::vector<Parameter> g0s) :
            _name(name),
            _m1(m1),
            _m2(m2),
            _N_orbital(N_orbital),
            _g0s(g0s)
        {
            if (m1 < 0 || m2 < 0)
                throw InternalError("Masses of channel '" + _name + "' must be positive");
        };


        // Phase space factor
        virtual double beta(const double & s) = 0;

        // Analytic continuation of the phase space factor
        virtual complex<double> rho(const double & s) = 0;
    };


    template <unsigned nchannels_, unsigned nresonances_>
    struct KMatrix<nchannels_, nresonances_>::Resonance
    {
        std::string _name;

        // Mass of the resonance
        Parameter _m;

        Resonance(std::string name, Parameter m) :
            _name(name),
            _m(m)
        {
            if (m < 0)
                throw InternalError("Mass of resonance '" + _name + "' must be positive");
        };
    };
}

#endif
