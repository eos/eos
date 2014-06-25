/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013 Danny van Dyk
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
    }

    WilsonScanComponent::WilsonScanComponent(const Parameters & p, const Options & o, ParameterUser & u) :
        _alpha_s_Z__deltab1(p["QCD::alpha_s(MZ)"], u),
        _mu_b__deltab1(p["QCD::mu_b"], u),
        _m_Z__deltab1(p["mass::Z"], u),
        _mu__deltab1(p["mu"], u),
        _c1(p["c1"], u),
        _c2(p["c2"], u),
        _c3(p["c3"], u),
        _c4(p["c4"], u),
        _c5(p["c5"], u),
        _c6(p["c6"], u),
        _abs_c7(p["Abs{c7}"]),
        _arg_c7(p["Arg{c7}"]),
        _re_c7(p["Re{c7}"]),
        _im_c7(p["Im{c7}"]),
        _c8(p["c8"], u),
        _abs_c9(p["Abs{c9}"]),
        _arg_c9(p["Arg{c9}"]),
        _re_c9(p["Re{c9}"]),
        _im_c9(p["Im{c9}"]),
        _abs_c10(p["Abs{c10}"]),
        _arg_c10(p["Arg{c10}"]),
        _re_c10(p["Re{c10}"]),
        _im_c10(p["Im{c10}"]),
        _abs_c7prime(p["Abs{c7'}"]),
        _arg_c7prime(p["Arg{c7'}"]),
        _re_c7prime(p["Re{c7'}"]),
        _im_c7prime(p["Im{c7'}"]),
        _c8prime(p["c8'"], u),
        _abs_c9prime(p["Abs{c9'}"]),
        _arg_c9prime(p["Arg{c9'}"]),
        _re_c9prime(p["Re{c9'}"]),
        _im_c9prime(p["Im{c9'}"]),
        _abs_c10prime(p["Abs{c10'}"]),
        _arg_c10prime(p["Arg{c10'}"]),
        _re_c10prime(p["Re{c10'}"]),
        _im_c10prime(p["Im{c10'}"]),
        _abs_cS(p["Abs{cS}"]),
        _arg_cS(p["Arg{cS}"]),
        _re_cS(p["Re{cS}"]),
        _im_cS(p["Im{cS}"]),
        _abs_cSprime(p["Abs{cS'}"]),
        _arg_cSprime(p["Arg{cS'}"]),
        _re_cSprime(p["Re{cS'}"]),
        _im_cSprime(p["Im{cS'}"]),
        _abs_cP(p["Abs{cP}"]),
        _arg_cP(p["Arg{cP}"]),
        _re_cP(p["Re{cP}"]),
        _im_cP(p["Im{cP}"]),
        _abs_cPprime(p["Abs{cP'}"]),
        _arg_cPprime(p["Arg{cP'}"]),
        _re_cPprime(p["Re{cP'}"]),
        _im_cPprime(p["Im{cP'}"]),
        _abs_cT(p["Abs{cT}"]),
        _arg_cT(p["Arg{cT}"]),
        _re_cT(p["Re{cT}"]),
        _im_cT(p["Im{cT}"]),
        _abs_cT5(p["Abs{cT5}"]),
        _arg_cT5(p["Arg{cT5}"]),
        _re_cT5(p["Re{cT5}"]),
        _im_cT5(p["Im{cT5}"])
    {
        if ("polar" == o.get("scan-mode", "polar"))
        {
            _c7 = std::bind(&wcimplementation::polar, _abs_c7, _arg_c7);
            u.uses(_abs_c7.id()); u.uses(_arg_c7.id());
            _c9 = std::bind(&wcimplementation::polar, _abs_c9, _arg_c9);
            u.uses(_abs_c9.id()); u.uses(_arg_c9.id());
            _c10 = std::bind(&wcimplementation::polar, _abs_c10, _arg_c10);
            u.uses(_abs_c10.id()); u.uses(_arg_c10.id());
            _c7prime = std::bind(&wcimplementation::polar, _abs_c7prime, _arg_c7prime);
            u.uses(_abs_c7prime.id()); u.uses(_arg_c7prime.id());
            _c9prime = std::bind(&wcimplementation::polar, _abs_c9prime, _arg_c9prime);
            u.uses(_abs_c9prime.id()); u.uses(_arg_c9prime.id());
            _c10prime = std::bind(&wcimplementation::polar, _abs_c10prime, _arg_c10prime);
            u.uses(_abs_c10prime.id()); u.uses(_arg_c10prime.id());
            _cS = std::bind(&wcimplementation::polar, _abs_cS, _arg_cS);
            u.uses(_abs_cS.id()); u.uses(_arg_cS.id());
            _cSprime = std::bind(&wcimplementation::polar, _abs_cSprime, _arg_cSprime);
            u.uses(_abs_cSprime.id()); u.uses(_arg_cSprime.id());
            _cP = std::bind(&wcimplementation::polar, _abs_cP, _arg_cP);
            u.uses(_abs_cP.id()); u.uses(_arg_cP.id());
            _cPprime = std::bind(&wcimplementation::polar, _abs_cPprime, _arg_cPprime);
            u.uses(_abs_cPprime.id()); u.uses(_arg_cPprime.id());
            _cT = std::bind(&wcimplementation::polar, _abs_cT, _arg_cT);
            u.uses(_abs_cT.id()); u.uses(_arg_cT.id());
            _cT5 = std::bind(&wcimplementation::polar, _abs_cT5, _arg_cT5);
            u.uses(_abs_cT5.id()); u.uses(_arg_cT5.id());
        }
        else if ("cartesian" == o.get("scan-mode", "polar"))
        {
            _c7 = std::bind(&wcimplementation::cartesian, _re_c7, _im_c7);
            u.uses(_re_c7.id()); u.uses(_im_c7.id());
            _c9 = std::bind(&wcimplementation::cartesian, _re_c9, _im_c9);
            u.uses(_re_c9.id()); u.uses(_im_c9.id());
            _c10 = std::bind(&wcimplementation::cartesian, _re_c10, _im_c10);
            u.uses(_re_c10.id()); u.uses(_im_c10.id());
            _c7prime = std::bind(&wcimplementation::cartesian, _re_c7prime, _im_c7prime);
            u.uses(_re_c7prime.id()); u.uses(_im_c7prime.id());
            _c9prime = std::bind(&wcimplementation::cartesian, _re_c9prime, _im_c9prime);
            u.uses(_re_c9prime.id()); u.uses(_im_c9prime.id());
            _c10prime = std::bind(&wcimplementation::cartesian, _re_c10prime, _im_c10prime);
            u.uses(_re_c10prime.id()); u.uses(_im_c10prime.id());
            _cS = std::bind(&wcimplementation::cartesian, _re_cS, _im_cS);
            u.uses(_re_cS.id()); u.uses(_im_cS.id());
            _cSprime = std::bind(&wcimplementation::cartesian, _re_cSprime, _im_cSprime);
            u.uses(_re_cSprime.id()); u.uses(_im_cSprime.id());
            _cP = std::bind(&wcimplementation::cartesian, _re_cP, _im_cP);
            u.uses(_re_cP.id()); u.uses(_im_cP.id());
            _cPprime = std::bind(&wcimplementation::cartesian, _re_cPprime, _im_cPprime);
            u.uses(_re_cPprime.id()); u.uses(_im_cPprime.id());
            _cT = std::bind(&wcimplementation::cartesian, _re_cT, _im_cT);
            u.uses(_re_cT.id()); u.uses(_im_cT.id());
            _cT5 = std::bind(&wcimplementation::cartesian, _re_cT5, _im_cT5);
            u.uses(_re_cT5.id()); u.uses(_im_cT5.id());
        }
        else
        {
            throw InternalError("scan-mode = '" + stringify(o.get("scan-mode", "polar")) + "' is not a valid scan mode for WilsonScanModel");
        }
    }

    /* b->s Wilson coefficients */
    WilsonCoefficients<BToS>
    WilsonScanComponent::wilson_coefficients_b_to_s(const bool & cp_conjugate) const
    {
        double alpha_s = 0.0;
        if (_mu__deltab1 < _mu_b__deltab1)
        {
            alpha_s = QCD::alpha_s(_mu_b__deltab1, _alpha_s_Z__deltab1, _m_Z__deltab1, QCD::beta_function_nf_5);
            alpha_s = QCD::alpha_s(_mu__deltab1, alpha_s, _mu_b__deltab1, QCD::beta_function_nf_4);
        }
        else
        {
            alpha_s = QCD::alpha_s(_mu__deltab1, _alpha_s_Z__deltab1, _m_Z__deltab1, QCD::beta_function_nf_5);
        }

        complex<double> a_s = alpha_s / 4.0 / M_PI;
        WilsonCoefficients<BToS> result
        {
            {{
                _c1(), _c2(), _c3(), _c4(), _c5(), _c6(),
                0.0, 0.0, 0.0, 0.0, 0.0,
                a_s * _c7(), a_s * _c8(), a_s * _c9(), a_s * _c10()
            }},
            {{
                /* we only consider c7', c8', c9' and c10' */
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
                a_s * _c7prime(), a_s * _c8prime(), a_s * _c9prime(), a_s * _c10prime()
            }},
            {{
                _cS(), _cSprime(), _cP(), _cPprime(), _cT(), _cT5()
            }},
            alpha_s
        };

        if (cp_conjugate)
        {
            for (auto & c : result._sm_like_coefficients)
                c = conj(c);
            for (auto & c : result._primed_coefficients)
                c = conj(c);
            for (auto & c : result._scalar_tensor_coefficients)
                c = conj(c);
        }

        return result;
    }

    WilsonScanModel::WilsonScanModel(const Parameters & parameters, const Options & options) :
        SMComponent<components::CKM>(parameters, *this),
        SMComponent<components::QCD>(parameters, *this),
        WilsonScanComponent(parameters, options, *this)
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
}
