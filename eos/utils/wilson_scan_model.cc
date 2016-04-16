/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2015 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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

#include <eos/utils/complex.hh>
#include <eos/utils/log.hh>
#include <eos/utils/matrix.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/wilson_scan_model.hh>

#include <cmath>

namespace eos
{
    using std::sqrt;

    namespace wcimplementation
    {
        complex<double> polar(const Parameter & abs, const Parameter & arg) { return std::polar(abs(), arg()); }
        complex<double> cartesian(const Parameter & re, const Parameter & im) { return complex<double>(re(), im()); }
        complex<double> polar_negative(const Parameter & abs, const Parameter & arg) { return std::polar(abs(), arg() + M_PI); }
        complex<double> cartesian_negative(const Parameter & re, const Parameter & im) { return complex<double>(-re(), -im()); }
        complex<double> zero() { return complex<double>(0.0, 0.0); }
    }

    /* b->s Wilson coefficients */
    WilsonScanComponent<components::DeltaBS1>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _alpha_s_Z__deltabs1(p["QCD::alpha_s(MZ)"], u),
        _mu_b__deltabs1(p["QCD::mu_b"], u),
        _m_Z__deltabs1(p["mass::Z"], u),
        _mu__deltabs1(p["mu"], u),
        /* b->s */
        _c1(p["c1"], u),
        _c2(p["c2"], u),
        _c3(p["c3"], u),
        _c4(p["c4"], u),
        _c5(p["c5"], u),
        _c6(p["c6"], u),
        _re_c7(p["Re{c7}"], u),
        _im_c7(p["Im{c7}"], u),
        _c8(p["c8"], u),
        _re_c9(p["Re{c9}"], u),
        _im_c9(p["Im{c9}"], u),
        _re_c10(p["Re{c10}"], u),
        _im_c10(p["Im{c10}"], u),
        _re_c7prime(p["Re{c7'}"], u),
        _im_c7prime(p["Im{c7'}"], u),
        _c8prime(p["c8'"], u),
        _re_c9prime(p["Re{c9'}"], u),
        _im_c9prime(p["Im{c9'}"], u),
        _re_c10prime(p["Re{c10'}"], u),
        _im_c10prime(p["Im{c10'}"], u),
        _re_cS(p["Re{cS}"], u),
        _im_cS(p["Im{cS}"], u),
        _re_cSprime(p["Re{cS'}"], u),
        _im_cSprime(p["Im{cS'}"], u),
        _re_cP(p["Re{cP}"], u),
        _im_cP(p["Im{cP}"], u),
        _re_cPprime(p["Re{cP'}"], u),
        _im_cPprime(p["Im{cP'}"], u),
        _re_cT(p["Re{cT}"], u),
        _im_cT(p["Im{cT}"], u),
        _re_cT5(p["Re{cT5}"], u),
        _im_cT5(p["Im{cT5}"], u),
        _c7(std::bind(&wcimplementation::cartesian, _re_c7, _im_c7)),
        _c9(std::bind(&wcimplementation::cartesian, _re_c9, _im_c9)),
        _c10(std::bind(&wcimplementation::cartesian, _re_c10, _im_c10)),
        _c7prime(std::bind(&wcimplementation::cartesian, _re_c7prime, _im_c7prime)),
        _c9prime(std::bind(&wcimplementation::cartesian, _re_c9prime, _im_c9prime)),
        _c10prime(std::bind(&wcimplementation::cartesian, _re_c10prime, _im_c10prime)),
        _cS(std::bind(&wcimplementation::cartesian, _re_cS, _im_cS)),
        _cSprime(std::bind(&wcimplementation::cartesian, _re_cSprime, _im_cSprime)),
        _cP(std::bind(&wcimplementation::cartesian, _re_cP, _im_cP)),
        _cPprime(std::bind(&wcimplementation::cartesian, _re_cPprime, _im_cPprime)),
        _cT(std::bind(&wcimplementation::cartesian, _re_cT, _im_cT)),
        _cT5(std::bind(&wcimplementation::cartesian, _re_cT5, _im_cT5))
    {
    }

    WilsonCoefficients<BToS>
    WilsonScanComponent<components::DeltaBS1>::wilson_coefficients_b_to_s(const bool & cp_conjugate) const
    {
        double alpha_s = 0.0;
        if (_mu__deltabs1 < _mu_b__deltabs1)
        {
            alpha_s = QCD::alpha_s(_mu_b__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
            alpha_s = QCD::alpha_s(_mu__deltabs1, alpha_s, _mu_b__deltabs1, QCD::beta_function_nf_4);
        }
        else
        {
            alpha_s = QCD::alpha_s(_mu__deltabs1, _alpha_s_Z__deltabs1, _m_Z__deltabs1, QCD::beta_function_nf_5);
        }

        complex<double> a_s = alpha_s / 4.0 / M_PI;
        WilsonCoefficients<BToS> result;
        result._sm_like_coefficients = std::array<std::complex<double>, 15>
        {
            _c1(), _c2(), _c3(), _c4(), _c5(), _c6(),
            0.0, 0.0, 0.0, 0.0, 0.0,
            a_s * _c7(), a_s * _c8(), a_s * _c9(), a_s * _c10()
        };
        result._primed_coefficients = std::array<std::complex<double>, 15>
        {
            /* we only consider c7', c8', c9' and c10' */
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            a_s * _c7prime(), a_s * _c8prime(), a_s * _c9prime(), a_s * _c10prime()
        };
        result._scalar_tensor_coefficients = std::array<std::complex<double>, 6>
        {
            _cS(), _cSprime(), _cP(), _cPprime(), _cT(), _cT5()
        };
        result._alpha_s = alpha_s;

        if (cp_conjugate)
        {
            for (auto & c : result._sm_like_coefficients)
            {
                c = conj(c);
            }

            for (auto & c : result._primed_coefficients)
            {
                c = conj(c);
            }

            for (auto & c : result._scalar_tensor_coefficients)
            {
                c = conj(c);
            }
        }

        return result;
    }

    /* b->u Wilson coefficients */
    WilsonScanComponent<components::DeltaBU1>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _re_csl(p["b->ulnu::Re{cSL}"], u),
        _im_csl(p["b->ulnu::Im{cSL}"], u),
        _re_csr(p["b->ulnu::Re{cSR}"], u),
        _im_csr(p["b->ulnu::Im{cSR}"], u),
        _re_cvl(p["b->ulnu::Re{cVL}"], u),
        _im_cvl(p["b->ulnu::Im{cVL}"], u),
        _re_cvr(p["b->ulnu::Re{cVR}"], u),
        _im_cvr(p["b->ulnu::Im{cVR}"], u),
        _re_ct(p["b->ulnu::Re{cT}"], u),
        _im_ct(p["b->ulnu::Im{cT}"], u),
        _csl(std::bind(&wcimplementation::cartesian, _re_csl, _im_csl)),
        _csr(std::bind(&wcimplementation::cartesian, _re_csr, _im_csr)),
        _cvl(std::bind(&wcimplementation::cartesian, _re_cvl, _im_cvl)),
        _cvr(std::bind(&wcimplementation::cartesian, _re_cvr, _im_cvr)),
        _ct(std::bind(&wcimplementation::cartesian, _re_ct, _im_ct))
    {
    }

    WilsonCoefficients<BToU>
    WilsonScanComponent<components::DeltaBU1>::wilson_coefficients_b_to_u(const bool & cp_conjugate) const
    {
        WilsonCoefficients<BToU> result
        {
            {{
                _cvl(), _cvr(), _csl(), _csr(), _ct()
            }},
        };

        if (cp_conjugate)
        {
            for (auto c = result._coefficients.begin(), c_end = result._coefficients.end() ; c != c_end ; ++c)
            {
                *c = conj(*c);
            }
        }

        return result;
    }

    ConstrainedWilsonScanComponent::ConstrainedWilsonScanComponent(const Parameters & p, const Options & o, ParameterUser & u) :
        WilsonScanComponent<components::DeltaBS1>(p, o, u)
    {
        _cT = std::bind(&wcimplementation::zero);
        _cT5 = std::bind(&wcimplementation::zero);
        _cP = std::bind(&wcimplementation::cartesian_negative, _re_cS, _im_cS);
        _cPprime = std::bind(&wcimplementation::cartesian, _re_cSprime, _im_cSprime);

        u.drop(_re_cP.id()); u.drop(_im_cP.id());
        u.drop(_re_cPprime.id()); u.drop(_im_cPprime.id());
        u.drop(_re_cT.id()); u.drop(_im_cT.id());
        u.drop(_re_cT5.id()); u.drop(_im_cT5.id());
    }

    WilsonScanModel::WilsonScanModel(const Parameters & parameters, const Options & options) :
        SMComponent<components::CKM>(parameters, *this),
        SMComponent<components::QCD>(parameters, *this),
        WilsonScanComponent<components::DeltaBS1>(parameters, options, *this),
        WilsonScanComponent<components::DeltaBU1>(parameters, options, *this)
    {
    }

    WilsonScanModel::~WilsonScanModel()
    {
    }

    std::shared_ptr<Model>
    WilsonScanModel::make(const Parameters & parameters, const Options & options)
    {
        return std::shared_ptr<Model>(new WilsonScanModel(parameters, options));
    }

    ConstrainedWilsonScanModel::ConstrainedWilsonScanModel(const Parameters & parameters, const Options & options) :
        SMComponent<components::CKM>(parameters, *this),
        SMComponent<components::QCD>(parameters, *this),
        ConstrainedWilsonScanComponent(parameters, options, *this),
        WilsonScanComponent<components::DeltaBU1>(parameters, options, *this)
    {
    }

    ConstrainedWilsonScanModel::~ConstrainedWilsonScanModel()
    {
    }

    std::shared_ptr<Model>
    ConstrainedWilsonScanModel::make(const Parameters & parameters, const Options & options)
    {
        return std::shared_ptr<Model>(new ConstrainedWilsonScanModel(parameters, options));
    }
}
