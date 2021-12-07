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

#ifndef EOS_GUARD_EOS_UTILS_KMATRIX_IMPL_HH
#define EOS_GUARD_EOS_UTILS_KMATRIX_IMPL_HH 1

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/kmatrix.hh>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <limits>

namespace eos
{
    template <unsigned nchannels_, unsigned nresonances_>
    KMatrix<nchannels_, nresonances_>::KMatrix(std::initializer_list<std::shared_ptr<KMatrix::Channel>> channels,
                                               std::initializer_list<std::shared_ptr<KMatrix::Resonance>> resonances,
                                               std::vector<std::vector<Parameter>> bkgcst,
                                               const std::string & prefix) :
        _channels(channels),
        _resonances(resonances),
        _bkgcst(bkgcst),
        _prefix(prefix),
        _Khat(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
        _That(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
        _tmp_1(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
        _tmp_2(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
        _perm(gsl_permutation_calloc(nchannels_))
    {
        // Perform size checks
        if (channels.size() != nchannels_)
            throw InternalError("The size of the channels array does not match nchannels_.");
        if (bkgcst.size() != nchannels_)
            throw InternalError("The size of the array of background constants does not match nchannels_.");
        if (bkgcst.size() != 0 and bkgcst[0].size() != nchannels_)
            throw InternalError("The array of background constants is not square.");
        if (resonances.size() != nresonances_)
            throw InternalError("The size of the resonances array does not match nresonances_.");

        // Perform pointer checks
        if (_Khat == nullptr)
            throw std::bad_alloc();
        if (_That == nullptr)
            throw std::bad_alloc();
        if (_tmp_1 == nullptr)
            throw std::bad_alloc();
        if (_tmp_2 == nullptr)
            throw std::bad_alloc();
        if (_perm == nullptr)
            throw std::bad_alloc();
    }

    template <unsigned nchannels_, unsigned nresonances_>
    KMatrix<nchannels_, nresonances_>::~KMatrix()
    {
        if (_perm)
            gsl_permutation_free(_perm);
        _perm = nullptr;

        if (_tmp_1)
            gsl_matrix_complex_free(_tmp_1);
        _tmp_1 = nullptr;

        if (_tmp_2)
            gsl_matrix_complex_free(_tmp_2);
        _tmp_2 = nullptr;

        if (_That)
            gsl_matrix_complex_free(_That);
        _That = nullptr;

        if (_Khat)
            gsl_matrix_complex_free(_Khat);
        _Khat = nullptr;
    }


    // Adapt s to avoid resonances masses
    template <unsigned nchannels_, unsigned nresonances_>
    double
    KMatrix<nchannels_, nresonances_>::adapt_s(const double s) const
    {
        // Disallowed range around resonance masses
        const double minimal_distance = 1.0e-7;
        const auto & resonances = this->_resonances;

        double adapted_s = s;

        for (size_t a = 0 ; a < nresonances_ ; ++a)
        {
            const double mres_a = resonances[a]->_m;

            for (size_t b = 0 ; b < a ; ++b)
            {
                const double mres_b = resonances[b]->_m;

                if (abs(mres_a * mres_a - mres_b * mres_b) < minimal_distance)
                    throw InternalError("The resonances masses are degenerate.");
            }

            if (abs(mres_a * mres_a - s) < minimal_distance) [[unlikely]]
            {
                if (s > mres_a * mres_a)
                    adapted_s = mres_a * mres_a + minimal_distance;
                else
                    adapted_s = mres_a * mres_a - minimal_distance;
            }
        }

        return adapted_s;
    }


    // Return the row corresponding to the index rowindex of the T matrix defined as T = (1-i*rho*K)^(-1)*K
    template <unsigned nchannels_, unsigned nresonances_>
    std::array<complex<double>, nchannels_>
    KMatrix<nchannels_, nresonances_>::tmatrix_row(unsigned rowindex, double _s) const
    {
        std::array<complex<double>, nchannels_> tmatrixrow;
        const auto & channels = this->_channels;
        const auto & resonances = this->_resonances;
        const auto & bkgcst = this->_bkgcst;
        // Adapt s to avoid pole in the K matrix
        const double s = adapt_s(_s);

        ///////////////////
        //1. Fill tmp1 = rho matrix
        ///////////////////
        gsl_matrix_complex_set_zero(_tmp_1);
        gsl_matrix_complex_set_zero(_tmp_2);

        for (size_t i = 0 ; i < nchannels_ ; ++i)
        {
            complex<double> rhoentry = channels[i]->rho(s);
            gsl_matrix_complex_set(_tmp_1,  i,  i, gsl_complex_rect(rhoentry.real(), rhoentry.imag()));
        }

        ///////////////////
        // 2. Fill Khat
        ///////////////////
        for (size_t i = 0 ; i < nchannels_ ; ++i)
        {
            for (size_t j = 0 ; j < nchannels_ ; ++j)
            {
                complex<double> entry = complex<double>(bkgcst[i][j].evaluate(), 0.0);

                for (size_t a = 0 ; a < nresonances_ ; ++a)
                {
                    const double mres = resonances[a]->_m;

                    Parameter g0rci = channels[i]->_g0s[a];
                    Parameter g0rcj = channels[j]->_g0s[a];

                    entry += g0rci * g0rcj / (mres * mres - s);

                }

                gsl_matrix_complex_set(_Khat, i, j, gsl_complex_rect(entry.real(), entry.imag()));
            }
        }

        ///////////////////
        // 3. Compute That
        ///////////////////
        static const gsl_complex one  = gsl_complex_rect(1.0, 0.0);
        static const gsl_complex minusi  = gsl_complex_rect(0.0, -1.0);
        static const gsl_complex zero = gsl_complex_rect(0.0, 0.0);

        // 3a. Set 1.0 to the diagonal elements of tmp_2
        for (unsigned i = 0 ; i < nchannels_ ; ++i)
        {
            gsl_matrix_complex_set(_tmp_2, i, i, one);
        }

        // 3b. multiply -i*rho with Khat from the right and add tmp_2
        // -> tmp_2 = -i * tmp_1 * Khat + one * tmp_2 ; no transpose anywhere
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, minusi, _tmp_1, _Khat, one, _tmp_2);

        // 3c. invert tmp_2 and assign the result to tmp_1
        int signum = 0;
        gsl_permutation_init(_perm);
        gsl_linalg_complex_LU_decomp(_tmp_2, _perm, &signum);
        gsl_linalg_complex_LU_invert(_tmp_2, _perm, _tmp_1);

        // 3d. calculate That = Khat / (1 - i rho * Khat) + zero * That
        //                    = Khat * tmp_1 + zero * That
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, _Khat, _tmp_1, zero, _That);


        ///////////////////
        // 4. Extract T matrix row
        ///////////////////
        for (size_t i = 0 ; i < nchannels_ ; ++i)
        {
            gsl_complex value = gsl_matrix_complex_get(_That, rowindex, i);
            complex<double> cvalue(GSL_REAL(value), GSL_IMAG(value));
            tmatrixrow[i] = cvalue;
        }

        return tmatrixrow;
    }

    template <unsigned nchannels_, unsigned nresonances_>
    double
    KMatrix<nchannels_, nresonances_>::partial_width(unsigned resonance, unsigned channel) const
    {
        double mres = this->_resonances[resonance]->_m;
        double rho  = real(this->_channels[channel]->rho(mres*mres));

        return power_of<2>((double)this->_channels[channel]->_g0s[resonance]) / mres * rho;
    }

    template <unsigned nchannels_, unsigned nresonances_>
    double
    KMatrix<nchannels_, nresonances_>::width(unsigned resonance) const
    {
        double result = 0.0;

        for (size_t i = 0; i < nchannels_; i++)
        {
            result += KMatrix<nchannels_, nresonances_>::partial_width(resonance, i);
        }

        return result;
    }
}
#endif
