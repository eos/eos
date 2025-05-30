/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2024 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014, 2018 Christoph Bobeth
 * Copyright (c) 2018 Ahmet Kokulu
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

#include <eos/maths/complex.hh>
#include <eos/maths/matrix.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/wet.hh>
#include <eos/utils/log.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

#include <cmath>

namespace eos
{
    using std::sqrt;

    namespace wcimplementation
    {
        complex<double>
        polar(const Parameter & abs, const Parameter & arg)
        {
            return std::polar(abs(), arg());
        }

        complex<double>
        cartesian(const Parameter & re, const Parameter & im)
        {
            return complex<double>(re(), im());
        }

        complex<double>
        polar_negative(const Parameter & abs, const Parameter & arg)
        {
            return std::polar(abs(), arg() + M_PI);
        }

        complex<double>
        cartesian_negative(const Parameter & re, const Parameter & im)
        {
            return complex<double>(-re(), -im());
        }

        complex<double>
        zero()
        {
            return complex<double>(0.0, 0.0);
        }
    } // namespace wcimplementation

    /* Charged-current semileptonic sectors (Delta S = 1) */

    /* [ubar s] [lbar nu] Wilson coefficients */
    WilsonScanComponent<components::WET::USLNu>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _e_re_csl(p["usenue::Re{cSL}"], u),
        _e_im_csl(p["usenue::Im{cSL}"], u),
        _e_re_csr(p["usenue::Re{cSR}"], u),
        _e_im_csr(p["usenue::Im{cSR}"], u),
        _e_re_cvl(p["usenue::Re{cVL}"], u),
        _e_im_cvl(p["usenue::Im{cVL}"], u),
        _e_re_cvr(p["usenue::Re{cVR}"], u),
        _e_im_cvr(p["usenue::Im{cVR}"], u),
        _e_re_ct(p["usenue::Re{cT}"], u),
        _e_im_ct(p["usenue::Im{cT}"], u),

        _mu_re_csl(p["usmunumu::Re{cSL}"], u),
        _mu_im_csl(p["usmunumu::Im{cSL}"], u),
        _mu_re_csr(p["usmunumu::Re{cSR}"], u),
        _mu_im_csr(p["usmunumu::Im{cSR}"], u),
        _mu_re_cvl(p["usmunumu::Re{cVL}"], u),
        _mu_im_cvl(p["usmunumu::Im{cVL}"], u),
        _mu_re_cvr(p["usmunumu::Re{cVR}"], u),
        _mu_im_cvr(p["usmunumu::Im{cVR}"], u),
        _mu_re_ct(p["usmunumu::Re{cT}"], u),
        _mu_im_ct(p["usmunumu::Im{cT}"], u),

        _tau_re_csl(p["ustaunutau::Re{cSL}"], u),
        _tau_im_csl(p["ustaunutau::Im{cSL}"], u),
        _tau_re_csr(p["ustaunutau::Re{cSR}"], u),
        _tau_im_csr(p["ustaunutau::Im{cSR}"], u),
        _tau_re_cvl(p["ustaunutau::Re{cVL}"], u),
        _tau_im_cvl(p["ustaunutau::Im{cVL}"], u),
        _tau_re_cvr(p["ustaunutau::Re{cVR}"], u),
        _tau_im_cvr(p["ustaunutau::Im{cVR}"], u),
        _tau_re_ct(p["ustaunutau::Re{cT}"], u),
        _tau_im_ct(p["ustaunutau::Im{cT}"], u),

        _e_csl(std::bind(&wcimplementation::cartesian, _e_re_csl, _e_im_csl)),
        _e_csr(std::bind(&wcimplementation::cartesian, _e_re_csr, _e_im_csr)),
        _e_cvl(std::bind(&wcimplementation::cartesian, _e_re_cvl, _e_im_cvl)),
        _e_cvr(std::bind(&wcimplementation::cartesian, _e_re_cvr, _e_im_cvr)),
        _e_ct(std::bind(&wcimplementation::cartesian, _e_re_ct, _e_im_ct)),

        _mu_csl(std::bind(&wcimplementation::cartesian, _mu_re_csl, _mu_im_csl)),
        _mu_csr(std::bind(&wcimplementation::cartesian, _mu_re_csr, _mu_im_csr)),
        _mu_cvl(std::bind(&wcimplementation::cartesian, _mu_re_cvl, _mu_im_cvl)),
        _mu_cvr(std::bind(&wcimplementation::cartesian, _mu_re_cvr, _mu_im_cvr)),
        _mu_ct(std::bind(&wcimplementation::cartesian, _mu_re_ct, _mu_im_ct)),

        _tau_csl(std::bind(&wcimplementation::cartesian, _tau_re_csl, _tau_im_csl)),
        _tau_csr(std::bind(&wcimplementation::cartesian, _tau_re_csr, _tau_im_csr)),
        _tau_cvl(std::bind(&wcimplementation::cartesian, _tau_re_cvl, _tau_im_cvl)),
        _tau_cvr(std::bind(&wcimplementation::cartesian, _tau_re_cvr, _tau_im_cvr)),
        _tau_ct(std::bind(&wcimplementation::cartesian, _tau_re_ct, _tau_im_ct))
    {
    }

    WilsonCoefficients<ChargedCurrent>
    WilsonScanComponent<components::WET::USLNu>::wet_uslnu(LeptonFlavor lepton_flavor, const bool & cp_conjugate) const
    {
        std::function<complex<double>()> cvl;
        std::function<complex<double>()> cvr;
        std::function<complex<double>()> csl;
        std::function<complex<double>()> csr;
        std::function<complex<double>()> ct;

        if (LeptonFlavor::electron == lepton_flavor)
        {
            cvl = _e_cvl;
            cvr = _e_cvr;
            csl = _e_csl;
            csr = _e_csr;
            ct  = _e_ct;
        }
        else if (LeptonFlavor::muon == lepton_flavor)
        {
            cvl = _mu_cvl;
            cvr = _mu_cvr;
            csl = _mu_csl;
            csr = _mu_csr;
            ct  = _mu_ct;
        }
        else if (LeptonFlavor::tauon == lepton_flavor)
        {
            cvl = _tau_cvl;
            cvr = _tau_cvr;
            csl = _tau_csl;
            csr = _tau_csr;
            ct  = _tau_ct;
        }
        else
        {
            throw InternalError("WilsonScan implements 'e', 'mu' and 'tau' lepton flavors");
        }

        WilsonCoefficients<ChargedCurrent> result{
            { { cvl(), cvr(), csl(), csr(), ct() } },
        };

        if (cp_conjugate)
        {
            for (auto & _coefficient : result._coefficients)
            {
                _coefficient = conj(_coefficient);
            }
        }

        return result;
    }

    /* Charged-current semileptonic sectors (Delta C = 1) */

    // [dbar c] [nubar l] Wilson coefficients
    WilsonScanComponent<components::WET::DCNuL>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _e_re_csl(p["dcnuee::Re{cSL}"], u),
        _e_im_csl(p["dcnuee::Im{cSL}"], u),
        _e_re_csr(p["dcnuee::Re{cSR}"], u),
        _e_im_csr(p["dcnuee::Im{cSR}"], u),
        _e_re_cvl(p["dcnuee::Re{cVL}"], u),
        _e_im_cvl(p["dcnuee::Im{cVL}"], u),
        _e_re_cvr(p["dcnuee::Re{cVR}"], u),
        _e_im_cvr(p["dcnuee::Im{cVR}"], u),
        _e_re_ct(p["dcnuee::Re{cT}"], u),
        _e_im_ct(p["dcnuee::Im{cT}"], u),

        _mu_re_csl(p["dcnumumu::Re{cSL}"], u),
        _mu_im_csl(p["dcnumumu::Im{cSL}"], u),
        _mu_re_csr(p["dcnumumu::Re{cSR}"], u),
        _mu_im_csr(p["dcnumumu::Im{cSR}"], u),
        _mu_re_cvl(p["dcnumumu::Re{cVL}"], u),
        _mu_im_cvl(p["dcnumumu::Im{cVL}"], u),
        _mu_re_cvr(p["dcnumumu::Re{cVR}"], u),
        _mu_im_cvr(p["dcnumumu::Im{cVR}"], u),
        _mu_re_ct(p["dcnumumu::Re{cT}"], u),
        _mu_im_ct(p["dcnumumu::Im{cT}"], u),

        _tau_re_csl(p["dcnutautau::Re{cSL}"], u),
        _tau_im_csl(p["dcnutautau::Im{cSL}"], u),
        _tau_re_csr(p["dcnutautau::Re{cSR}"], u),
        _tau_im_csr(p["dcnutautau::Im{cSR}"], u),
        _tau_re_cvl(p["dcnutautau::Re{cVL}"], u),
        _tau_im_cvl(p["dcnutautau::Im{cVL}"], u),
        _tau_re_cvr(p["dcnutautau::Re{cVR}"], u),
        _tau_im_cvr(p["dcnutautau::Im{cVR}"], u),
        _tau_re_ct(p["dcnutautau::Re{cT}"], u),
        _tau_im_ct(p["dcnutautau::Im{cT}"], u),

        _e_csl(std::bind(&wcimplementation::cartesian, _e_re_csl, _e_im_csl)),
        _e_csr(std::bind(&wcimplementation::cartesian, _e_re_csr, _e_im_csr)),
        _e_cvl(std::bind(&wcimplementation::cartesian, _e_re_cvl, _e_im_cvl)),
        _e_cvr(std::bind(&wcimplementation::cartesian, _e_re_cvr, _e_im_cvr)),
        _e_ct(std::bind(&wcimplementation::cartesian, _e_re_ct, _e_im_ct)),

        _mu_csl(std::bind(&wcimplementation::cartesian, _mu_re_csl, _mu_im_csl)),
        _mu_csr(std::bind(&wcimplementation::cartesian, _mu_re_csr, _mu_im_csr)),
        _mu_cvl(std::bind(&wcimplementation::cartesian, _mu_re_cvl, _mu_im_cvl)),
        _mu_cvr(std::bind(&wcimplementation::cartesian, _mu_re_cvr, _mu_im_cvr)),
        _mu_ct(std::bind(&wcimplementation::cartesian, _mu_re_ct, _mu_im_ct)),

        _tau_csl(std::bind(&wcimplementation::cartesian, _tau_re_csl, _tau_im_csl)),
        _tau_csr(std::bind(&wcimplementation::cartesian, _tau_re_csr, _tau_im_csr)),
        _tau_cvl(std::bind(&wcimplementation::cartesian, _tau_re_cvl, _tau_im_cvl)),
        _tau_cvr(std::bind(&wcimplementation::cartesian, _tau_re_cvr, _tau_im_cvr)),
        _tau_ct(std::bind(&wcimplementation::cartesian, _tau_re_ct, _tau_im_ct))
    {
    }

    WilsonCoefficients<bern::ClassII>
    WilsonScanComponent<components::WET::DCNuL>::wet_dcnul(LeptonFlavor lepton_flavor, const bool & cp_conjugate) const
    {
        std::function<complex<double>()> cvl;
        std::function<complex<double>()> cvr;
        std::function<complex<double>()> csl;
        std::function<complex<double>()> csr;
        std::function<complex<double>()> ct;

        if (LeptonFlavor::electron == lepton_flavor)
        {
            cvl = _e_cvl;
            cvr = _e_cvr;
            csl = _e_csl;
            csr = _e_csr;
            ct  = _e_ct;
        }
        else if (LeptonFlavor::muon == lepton_flavor)
        {
            cvl = _mu_cvl;
            cvr = _mu_cvr;
            csl = _mu_csl;
            csr = _mu_csr;
            ct  = _mu_ct;
        }
        else if (LeptonFlavor::tauon == lepton_flavor)
        {
            cvl = _tau_cvl;
            cvr = _tau_cvr;
            csl = _tau_csl;
            csr = _tau_csr;
            ct  = _tau_ct;
        }
        else
        {
            throw InternalError("WilsonScan implements 'e', 'mu' and 'tau' lepton flavors");
        }

        WilsonCoefficients<ChargedCurrent> result{
            { { cvl(), cvr(), csl(), csr(), ct() } },
        };

        if (cp_conjugate)
        {
            for (auto & _coefficient : result._coefficients)
            {
                _coefficient = conj(_coefficient);
            }
        }

        return result;
    }

    /* [sbar c] [nubar l] Wilson coefficients */
    WilsonScanComponent<components::WET::SCNuL>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _e_re_csl(p["scnuee::Re{cSL}"], u),
        _e_im_csl(p["scnuee::Im{cSL}"], u),
        _e_re_csr(p["scnuee::Re{cSR}"], u),
        _e_im_csr(p["scnuee::Im{cSR}"], u),
        _e_re_cvl(p["scnuee::Re{cVL}"], u),
        _e_im_cvl(p["scnuee::Im{cVL}"], u),
        _e_re_cvr(p["scnuee::Re{cVR}"], u),
        _e_im_cvr(p["scnuee::Im{cVR}"], u),
        _e_re_ct(p["scnuee::Re{cT}"], u),
        _e_im_ct(p["scnuee::Im{cT}"], u),

        _mu_re_csl(p["scnumumu::Re{cSL}"], u),
        _mu_im_csl(p["scnumumu::Im{cSL}"], u),
        _mu_re_csr(p["scnumumu::Re{cSR}"], u),
        _mu_im_csr(p["scnumumu::Im{cSR}"], u),
        _mu_re_cvl(p["scnumumu::Re{cVL}"], u),
        _mu_im_cvl(p["scnumumu::Im{cVL}"], u),
        _mu_re_cvr(p["scnumumu::Re{cVR}"], u),
        _mu_im_cvr(p["scnumumu::Im{cVR}"], u),
        _mu_re_ct(p["scnumumu::Re{cT}"], u),
        _mu_im_ct(p["scnumumu::Im{cT}"], u),

        _tau_re_csl(p["scnutautau::Re{cSL}"], u),
        _tau_im_csl(p["scnutautau::Im{cSL}"], u),
        _tau_re_csr(p["scnutautau::Re{cSR}"], u),
        _tau_im_csr(p["scnutautau::Im{cSR}"], u),
        _tau_re_cvl(p["scnutautau::Re{cVL}"], u),
        _tau_im_cvl(p["scnutautau::Im{cVL}"], u),
        _tau_re_cvr(p["scnutautau::Re{cVR}"], u),
        _tau_im_cvr(p["scnutautau::Im{cVR}"], u),
        _tau_re_ct(p["scnutautau::Re{cT}"], u),
        _tau_im_ct(p["scnutautau::Im{cT}"], u),

        _e_csl(std::bind(&wcimplementation::cartesian, _e_re_csl, _e_im_csl)),
        _e_csr(std::bind(&wcimplementation::cartesian, _e_re_csr, _e_im_csr)),
        _e_cvl(std::bind(&wcimplementation::cartesian, _e_re_cvl, _e_im_cvl)),
        _e_cvr(std::bind(&wcimplementation::cartesian, _e_re_cvr, _e_im_cvr)),
        _e_ct(std::bind(&wcimplementation::cartesian, _e_re_ct, _e_im_ct)),

        _mu_csl(std::bind(&wcimplementation::cartesian, _mu_re_csl, _mu_im_csl)),
        _mu_csr(std::bind(&wcimplementation::cartesian, _mu_re_csr, _mu_im_csr)),
        _mu_cvl(std::bind(&wcimplementation::cartesian, _mu_re_cvl, _mu_im_cvl)),
        _mu_cvr(std::bind(&wcimplementation::cartesian, _mu_re_cvr, _mu_im_cvr)),
        _mu_ct(std::bind(&wcimplementation::cartesian, _mu_re_ct, _mu_im_ct)),

        _tau_csl(std::bind(&wcimplementation::cartesian, _tau_re_csl, _tau_im_csl)),
        _tau_csr(std::bind(&wcimplementation::cartesian, _tau_re_csr, _tau_im_csr)),
        _tau_cvl(std::bind(&wcimplementation::cartesian, _tau_re_cvl, _tau_im_cvl)),
        _tau_cvr(std::bind(&wcimplementation::cartesian, _tau_re_cvr, _tau_im_cvr)),
        _tau_ct(std::bind(&wcimplementation::cartesian, _tau_re_ct, _tau_im_ct))
    {
    }

    WilsonCoefficients<bern::ClassII>
    WilsonScanComponent<components::WET::SCNuL>::wet_scnul(LeptonFlavor lepton_flavor, const bool & cp_conjugate) const
    {
        std::function<complex<double>()> cvl;
        std::function<complex<double>()> cvr;
        std::function<complex<double>()> csl;
        std::function<complex<double>()> csr;
        std::function<complex<double>()> ct;

        if (LeptonFlavor::electron == lepton_flavor)
        {
            cvl = _e_cvl;
            cvr = _e_cvr;
            csl = _e_csl;
            csr = _e_csr;
            ct  = _e_ct;
        }
        else if (LeptonFlavor::muon == lepton_flavor)
        {
            cvl = _mu_cvl;
            cvr = _mu_cvr;
            csl = _mu_csl;
            csr = _mu_csr;
            ct  = _mu_ct;
        }
        else if (LeptonFlavor::tauon == lepton_flavor)
        {
            cvl = _tau_cvl;
            cvr = _tau_cvr;
            csl = _tau_csl;
            csr = _tau_csr;
            ct  = _tau_ct;
        }
        else
        {
            throw InternalError("WilsonScan implements 'e', 'mu' and 'tau' lepton flavors");
        }

        WilsonCoefficients<ChargedCurrent> result{
            { { cvl(), cvr(), csl(), csr(), ct() } },
        };

        if (cp_conjugate)
        {
            for (auto & _coefficient : result._coefficients)
            {
                _coefficient = conj(_coefficient);
            }
        }

        return result;
    }

    /* Charged-current semileptonic sectors (Delta B = 1) */

    /* [ubar b] [lbar nu] Wilson coefficients */
    WilsonScanComponent<components::WET::UBLNu>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _e_re_csl(p["ubenue::Re{cSL}"], u),
        _e_im_csl(p["ubenue::Im{cSL}"], u),
        _e_re_csr(p["ubenue::Re{cSR}"], u),
        _e_im_csr(p["ubenue::Im{cSR}"], u),
        _e_re_cvl(p["ubenue::Re{cVL}"], u),
        _e_im_cvl(p["ubenue::Im{cVL}"], u),
        _e_re_cvr(p["ubenue::Re{cVR}"], u),
        _e_im_cvr(p["ubenue::Im{cVR}"], u),
        _e_re_ct(p["ubenue::Re{cT}"], u),
        _e_im_ct(p["ubenue::Im{cT}"], u),

        _mu_re_csl(p["ubmunumu::Re{cSL}"], u),
        _mu_im_csl(p["ubmunumu::Im{cSL}"], u),
        _mu_re_csr(p["ubmunumu::Re{cSR}"], u),
        _mu_im_csr(p["ubmunumu::Im{cSR}"], u),
        _mu_re_cvl(p["ubmunumu::Re{cVL}"], u),
        _mu_im_cvl(p["ubmunumu::Im{cVL}"], u),
        _mu_re_cvr(p["ubmunumu::Re{cVR}"], u),
        _mu_im_cvr(p["ubmunumu::Im{cVR}"], u),
        _mu_re_ct(p["ubmunumu::Re{cT}"], u),
        _mu_im_ct(p["ubmunumu::Im{cT}"], u),

        _tau_re_csl(p["ubtaunutau::Re{cSL}"], u),
        _tau_im_csl(p["ubtaunutau::Im{cSL}"], u),
        _tau_re_csr(p["ubtaunutau::Re{cSR}"], u),
        _tau_im_csr(p["ubtaunutau::Im{cSR}"], u),
        _tau_re_cvl(p["ubtaunutau::Re{cVL}"], u),
        _tau_im_cvl(p["ubtaunutau::Im{cVL}"], u),
        _tau_re_cvr(p["ubtaunutau::Re{cVR}"], u),
        _tau_im_cvr(p["ubtaunutau::Im{cVR}"], u),
        _tau_re_ct(p["ubtaunutau::Re{cT}"], u),
        _tau_im_ct(p["ubtaunutau::Im{cT}"], u),

        _e_csl(std::bind(&wcimplementation::cartesian, _e_re_csl, _e_im_csl)),
        _e_csr(std::bind(&wcimplementation::cartesian, _e_re_csr, _e_im_csr)),
        _e_cvl(std::bind(&wcimplementation::cartesian, _e_re_cvl, _e_im_cvl)),
        _e_cvr(std::bind(&wcimplementation::cartesian, _e_re_cvr, _e_im_cvr)),
        _e_ct(std::bind(&wcimplementation::cartesian, _e_re_ct, _e_im_ct)),

        _mu_csl(std::bind(&wcimplementation::cartesian, _mu_re_csl, _mu_im_csl)),
        _mu_csr(std::bind(&wcimplementation::cartesian, _mu_re_csr, _mu_im_csr)),
        _mu_cvl(std::bind(&wcimplementation::cartesian, _mu_re_cvl, _mu_im_cvl)),
        _mu_cvr(std::bind(&wcimplementation::cartesian, _mu_re_cvr, _mu_im_cvr)),
        _mu_ct(std::bind(&wcimplementation::cartesian, _mu_re_ct, _mu_im_ct)),

        _tau_csl(std::bind(&wcimplementation::cartesian, _tau_re_csl, _tau_im_csl)),
        _tau_csr(std::bind(&wcimplementation::cartesian, _tau_re_csr, _tau_im_csr)),
        _tau_cvl(std::bind(&wcimplementation::cartesian, _tau_re_cvl, _tau_im_cvl)),
        _tau_cvr(std::bind(&wcimplementation::cartesian, _tau_re_cvr, _tau_im_cvr)),
        _tau_ct(std::bind(&wcimplementation::cartesian, _tau_re_ct, _tau_im_ct))
    {
    }

    WilsonCoefficients<ChargedCurrent>
    WilsonScanComponent<components::WET::UBLNu>::wet_ublnu(LeptonFlavor lepton_flavor, const bool & cp_conjugate) const
    {
        std::function<complex<double>()> cvl;
        std::function<complex<double>()> cvr;
        std::function<complex<double>()> csl;
        std::function<complex<double>()> csr;
        std::function<complex<double>()> ct;

        if (LeptonFlavor::electron == lepton_flavor)
        {
            cvl = _e_cvl;
            cvr = _e_cvr;
            csl = _e_csl;
            csr = _e_csr;
            ct  = _e_ct;
        }
        else if (LeptonFlavor::muon == lepton_flavor)
        {
            cvl = _mu_cvl;
            cvr = _mu_cvr;
            csl = _mu_csl;
            csr = _mu_csr;
            ct  = _mu_ct;
        }
        else if (LeptonFlavor::tauon == lepton_flavor)
        {
            cvl = _tau_cvl;
            cvr = _tau_cvr;
            csl = _tau_csl;
            csr = _tau_csr;
            ct  = _tau_ct;
        }
        else
        {
            throw InternalError("WilsonScan implements 'e', 'mu' and 'tau' lepton flavors");
        }

        WilsonCoefficients<ChargedCurrent> result{
            { { cvl(), cvr(), csl(), csr(), ct() } },
        };

        if (cp_conjugate)
        {
            for (auto & _coefficient : result._coefficients)
            {
                _coefficient = conj(_coefficient);
            }
        }

        return result;
    }

    /* [cbar b] [lbar nu] Wilson coefficients */
    WilsonScanComponent<components::WET::CBLNu>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _e_re_csl(p["cbenue::Re{cSL}"], u),
        _e_im_csl(p["cbenue::Im{cSL}"], u),
        _e_re_csr(p["cbenue::Re{cSR}"], u),
        _e_im_csr(p["cbenue::Im{cSR}"], u),
        _e_re_cvl(p["cbenue::Re{cVL}"], u),
        _e_im_cvl(p["cbenue::Im{cVL}"], u),
        _e_re_cvr(p["cbenue::Re{cVR}"], u),
        _e_im_cvr(p["cbenue::Im{cVR}"], u),
        _e_re_ct(p["cbenue::Re{cT}"], u),
        _e_im_ct(p["cbenue::Im{cT}"], u),

        _mu_re_csl(p["cbmunumu::Re{cSL}"], u),
        _mu_im_csl(p["cbmunumu::Im{cSL}"], u),
        _mu_re_csr(p["cbmunumu::Re{cSR}"], u),
        _mu_im_csr(p["cbmunumu::Im{cSR}"], u),
        _mu_re_cvl(p["cbmunumu::Re{cVL}"], u),
        _mu_im_cvl(p["cbmunumu::Im{cVL}"], u),
        _mu_re_cvr(p["cbmunumu::Re{cVR}"], u),
        _mu_im_cvr(p["cbmunumu::Im{cVR}"], u),
        _mu_re_ct(p["cbmunumu::Re{cT}"], u),
        _mu_im_ct(p["cbmunumu::Im{cT}"], u),

        _tau_re_csl(p["cbtaunutau::Re{cSL}"], u),
        _tau_im_csl(p["cbtaunutau::Im{cSL}"], u),
        _tau_re_csr(p["cbtaunutau::Re{cSR}"], u),
        _tau_im_csr(p["cbtaunutau::Im{cSR}"], u),
        _tau_re_cvl(p["cbtaunutau::Re{cVL}"], u),
        _tau_im_cvl(p["cbtaunutau::Im{cVL}"], u),
        _tau_re_cvr(p["cbtaunutau::Re{cVR}"], u),
        _tau_im_cvr(p["cbtaunutau::Im{cVR}"], u),
        _tau_re_ct(p["cbtaunutau::Re{cT}"], u),
        _tau_im_ct(p["cbtaunutau::Im{cT}"], u),

        _e_csl(std::bind(&wcimplementation::cartesian, _e_re_csl, _e_im_csl)),
        _e_csr(std::bind(&wcimplementation::cartesian, _e_re_csr, _e_im_csr)),
        _e_cvl(std::bind(&wcimplementation::cartesian, _e_re_cvl, _e_im_cvl)),
        _e_cvr(std::bind(&wcimplementation::cartesian, _e_re_cvr, _e_im_cvr)),
        _e_ct(std::bind(&wcimplementation::cartesian, _e_re_ct, _e_im_ct)),

        _mu_csl(std::bind(&wcimplementation::cartesian, _mu_re_csl, _mu_im_csl)),
        _mu_csr(std::bind(&wcimplementation::cartesian, _mu_re_csr, _mu_im_csr)),
        _mu_cvl(std::bind(&wcimplementation::cartesian, _mu_re_cvl, _mu_im_cvl)),
        _mu_cvr(std::bind(&wcimplementation::cartesian, _mu_re_cvr, _mu_im_cvr)),
        _mu_ct(std::bind(&wcimplementation::cartesian, _mu_re_ct, _mu_im_ct)),

        _tau_csl(std::bind(&wcimplementation::cartesian, _tau_re_csl, _tau_im_csl)),
        _tau_csr(std::bind(&wcimplementation::cartesian, _tau_re_csr, _tau_im_csr)),
        _tau_cvl(std::bind(&wcimplementation::cartesian, _tau_re_cvl, _tau_im_cvl)),
        _tau_cvr(std::bind(&wcimplementation::cartesian, _tau_re_cvr, _tau_im_cvr)),
        _tau_ct(std::bind(&wcimplementation::cartesian, _tau_re_ct, _tau_im_ct))
    {
    }

    WilsonCoefficients<ChargedCurrent>
    WilsonScanComponent<components::WET::CBLNu>::wet_cblnu(LeptonFlavor lepton_flavor, const bool & cp_conjugate) const
    {
        std::function<complex<double>()> cvl;
        std::function<complex<double>()> cvr;
        std::function<complex<double>()> csl;
        std::function<complex<double>()> csr;
        std::function<complex<double>()> ct;

        if (LeptonFlavor::electron == lepton_flavor)
        {
            cvl = _e_cvl;
            cvr = _e_cvr;
            csl = _e_csl;
            csr = _e_csr;
            ct  = _e_ct;
        }
        else if (LeptonFlavor::muon == lepton_flavor)
        {
            cvl = _mu_cvl;
            cvr = _mu_cvr;
            csl = _mu_csl;
            csr = _mu_csr;
            ct  = _mu_ct;
        }
        else if (LeptonFlavor::tauon == lepton_flavor)
        {
            cvl = _tau_cvl;
            cvr = _tau_cvr;
            csl = _tau_csl;
            csr = _tau_csr;
            ct  = _tau_ct;
        }
        else
        {
            throw InternalError("WilsonScan implements 'e', 'mu' and 'tau' lepton flavors");
        }

        WilsonCoefficients<ChargedCurrent> result{
            { { cvl(), cvr(), csl(), csr(), ct() } },
        };

        if (cp_conjugate)
        {
            for (auto & _coefficient : result._coefficients)
            {
                _coefficient = conj(_coefficient);
            }
        }

        return result;
    }

    /* Neutral-current semileptonic sectors (Delta B = 1) */

    /* [sbar b] [nubar nu] Wilson coefficients */
    WilsonScanComponent<components::WET::SBNuNu>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _re_cvl(p["sbnunu::Re{cVL}"], u),
        _im_cvl(p["sbnunu::Im{cVL}"], u),
        _re_cvr(p["sbnunu::Re{cVR}"], u),
        _im_cvr(p["sbnunu::Im{cVR}"], u),
        _re_csl(p["sbnunu::Re{cSL}"], u),
        _im_csl(p["sbnunu::Im{cSL}"], u),
        _re_csr(p["sbnunu::Re{cSR}"], u),
        _im_csr(p["sbnunu::Im{cSR}"], u),
        _re_ctl(p["sbnunu::Re{cTL}"], u),
        _im_ctl(p["sbnunu::Im{cTL}"], u)
    {
    }

    WilsonCoefficients<wc::SBNuNu>
    WilsonScanComponent<components::WET::SBNuNu>::wet_sbnunu(const bool & cp_conjugate) const
    {
        WilsonCoefficients<wc::SBNuNu> result;

        result._coefficients = std::array<complex<double>, 5>{
            { complex<double>(_re_cvl, _im_cvl),
             complex<double>(_re_cvr, _im_cvr),
             complex<double>(_re_csl, _im_csl),
             complex<double>(_re_csr, _im_csr),
             complex<double>(_re_ctl, _im_ctl) }
        };

        if (cp_conjugate)
        {
            for (auto & _coefficient : result._coefficients)
            {
                _coefficient = conj(_coefficient);
            }
        }

        return result;
    }

    /* Hadronic sectors (Delta B = 1) */

    std::array<std::tuple<UsedParameter, UsedParameter>, 20>
    make_wet_parameters_classIII(const Parameters & p, ParameterUser & u, const std::string & prefix)
    {
        auto result =
                std::array<std::tuple<UsedParameter, UsedParameter>, 20>{ std::make_tuple(UsedParameter(p[prefix + "::Re{c1}"], u), UsedParameter(p[prefix + "::Im{c1}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c2}"], u), UsedParameter(p[prefix + "::Im{c2}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c3}"], u), UsedParameter(p[prefix + "::Im{c3}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c4}"], u), UsedParameter(p[prefix + "::Im{c4}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c5}"], u), UsedParameter(p[prefix + "::Im{c5}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c6}"], u), UsedParameter(p[prefix + "::Im{c6}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c7}"], u), UsedParameter(p[prefix + "::Im{c7}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c8}"], u), UsedParameter(p[prefix + "::Im{c8}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c9}"], u), UsedParameter(p[prefix + "::Im{c9}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c10}"], u), UsedParameter(p[prefix + "::Im{c10}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c1'}"], u), UsedParameter(p[prefix + "::Im{c1'}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c2'}"], u), UsedParameter(p[prefix + "::Im{c2'}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c3'}"], u), UsedParameter(p[prefix + "::Im{c3'}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c4'}"], u), UsedParameter(p[prefix + "::Im{c4'}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c5'}"], u), UsedParameter(p[prefix + "::Im{c5'}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c6'}"], u), UsedParameter(p[prefix + "::Im{c6'}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c7'}"], u), UsedParameter(p[prefix + "::Im{c7'}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c8'}"], u), UsedParameter(p[prefix + "::Im{c8'}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c9'}"], u), UsedParameter(p[prefix + "::Im{c9'}"], u)),
                                                                          std::make_tuple(UsedParameter(p[prefix + "::Re{c10'}"], u), UsedParameter(p[prefix + "::Im{c10'}"], u)) };

        return result;
    }

    /* [dbar b] [cbar u] Wilson coefficients */
    WilsonScanComponent<components::WET::DBCU>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _dbcu_parameters(make_wet_parameters_classIII(p, u, "dbcu"))
    {
    }

    WilsonCoefficients<wc::DBCU>
    WilsonScanComponent<components::WET::DBCU>::wet_dbcu(const bool & cp_conjugate) const
    {
        const double sign = cp_conjugate ? -1.0 : 1.0;

        WilsonCoefficients<wc::DBCU> result;
        for (unsigned i = 0; i < 10; ++i)
        {
            result._unprimed[i] = complex<double>(std::get<0>(_dbcu_parameters[i]), sign * std::get<1>(_dbcu_parameters[i]));
            result._primed[i]   = complex<double>(std::get<0>(_dbcu_parameters[i + 10]), sign * std::get<1>(_dbcu_parameters[i + 10]));
        }

        return result;
    }

    /* [sbar b] [cbar u] Wilson coefficients */
    WilsonScanComponent<components::WET::SBCU>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _sbcu_parameters(make_wet_parameters_classIII(p, u, "sbcu"))
    {
    }

    WilsonCoefficients<wc::SBCU>
    WilsonScanComponent<components::WET::SBCU>::wet_sbcu(const bool & cp_conjugate) const
    {
        const double sign = cp_conjugate ? -1.0 : 1.0;

        WilsonCoefficients<wc::SBCU> result;
        for (unsigned i = 0; i < 10; ++i)
        {
            result._unprimed[i] = complex<double>(std::get<0>(_sbcu_parameters[i]), sign * std::get<1>(_sbcu_parameters[i]));
            result._primed[i]   = complex<double>(std::get<0>(_sbcu_parameters[i + 10]), sign * std::get<1>(_sbcu_parameters[i + 10]));
        }

        return result;
    }

    /* Hadronic sectors (Delta B = 2) */

    /* [sbar b] [sbar b] Wilson coefficients */
    WilsonScanComponent<components::WET::SBSB>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _re_sbsb_c1__deltab2(p["sbsb::Re{c1}"], u),
        _im_sbsb_c1__deltab2(p["sbsb::Im{c1}"], u),
        _re_sbsb_c2__deltab2(p["sbsb::Re{c2}"], u),
        _im_sbsb_c2__deltab2(p["sbsb::Im{c2}"], u),
        _re_sbsb_c3__deltab2(p["sbsb::Re{c3}"], u),
        _im_sbsb_c3__deltab2(p["sbsb::Im{c3}"], u),
        _re_sbsb_c4__deltab2(p["sbsb::Re{c4}"], u),
        _im_sbsb_c4__deltab2(p["sbsb::Im{c4}"], u),
        _re_sbsb_c5__deltab2(p["sbsb::Re{c5}"], u),
        _im_sbsb_c5__deltab2(p["sbsb::Im{c5}"], u),
        _re_sbsb_c1p__deltab2(p["sbsb::Re{c1'}"], u),
        _im_sbsb_c1p__deltab2(p["sbsb::Im{c1'}"], u),
        _re_sbsb_c2p__deltab2(p["sbsb::Re{c2'}"], u),
        _im_sbsb_c2p__deltab2(p["sbsb::Im{c2'}"], u),
        _re_sbsb_c3p__deltab2(p["sbsb::Re{c3'}"], u),
        _im_sbsb_c3p__deltab2(p["sbsb::Im{c3'}"], u)
    {
    }

    WilsonCoefficients<wc::SBSB>
    WilsonScanComponent<components::WET::SBSB>::wet_sbsb() const
    {
        WilsonCoefficients<wc::SBSB> result;

        result._coefficients = std::array<complex<double>, 8>{
            { complex<double>(_re_sbsb_c1__deltab2(), _im_sbsb_c1__deltab2()),
             complex<double>(_re_sbsb_c2__deltab2(), _im_sbsb_c2__deltab2()),
             complex<double>(_re_sbsb_c3__deltab2(), _im_sbsb_c3__deltab2()),
             complex<double>(_re_sbsb_c4__deltab2(), _im_sbsb_c4__deltab2()),
             complex<double>(_re_sbsb_c5__deltab2(), _im_sbsb_c5__deltab2()),
             complex<double>(_re_sbsb_c1p__deltab2(), _im_sbsb_c1p__deltab2()),
             complex<double>(_re_sbsb_c2p__deltab2(), _im_sbsb_c2p__deltab2()),
             complex<double>(_re_sbsb_c3p__deltab2(), _im_sbsb_c3p__deltab2()) }
        };

        return result;
    }

    /* Old-style WET sectors */

    /* b->s Wilson coefficients */
    WilsonScanComponent<components::DeltaBS1>::WilsonScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _alpha_s_Z__deltabs1(p["QCD::alpha_s(MZ)"], u),
        _mu_b__deltabs1(p["QCD::mu_b"], u),
        _m_Z__deltabs1(p["mass::Z"], u),
        _mu__deltabs1(p["sb::mu"], u),
        /* b->s */
        _c1(p["b->s::c1"], u),
        _c2(p["b->s::c2"], u),
        _c3(p["b->s::c3"], u),
        _c4(p["b->s::c4"], u),
        _c5(p["b->s::c5"], u),
        _c6(p["b->s::c6"], u),
        _re_c7(p["b->s::Re{c7}"], u),
        _im_c7(p["b->s::Im{c7}"], u),
        _re_c7prime(p["b->s::Re{c7'}"], u),
        _im_c7prime(p["b->s::Im{c7'}"], u),
        _c8(p["b->s::c8"], u),
        _c8prime(p["b->s::c8'"], u),
        /* b->see */
        _e_re_c9(p["b->see::Re{c9}"], u),
        _e_im_c9(p["b->see::Im{c9}"], u),
        _e_re_c10(p["b->see::Re{c10}"], u),
        _e_im_c10(p["b->see::Im{c10}"], u),
        _e_re_c9prime(p["b->see::Re{c9'}"], u),
        _e_im_c9prime(p["b->see::Im{c9'}"], u),
        _e_re_c10prime(p["b->see::Re{c10'}"], u),
        _e_im_c10prime(p["b->see::Im{c10'}"], u),
        _e_re_cS(p["b->see::Re{cS}"], u),
        _e_im_cS(p["b->see::Im{cS}"], u),
        _e_re_cSprime(p["b->see::Re{cS'}"], u),
        _e_im_cSprime(p["b->see::Im{cS'}"], u),
        _e_re_cP(p["b->see::Re{cP}"], u),
        _e_im_cP(p["b->see::Im{cP}"], u),
        _e_re_cPprime(p["b->see::Re{cP'}"], u),
        _e_im_cPprime(p["b->see::Im{cP'}"], u),
        _e_re_cT(p["b->see::Re{cT}"], u),
        _e_im_cT(p["b->see::Im{cT}"], u),
        _e_re_cT5(p["b->see::Re{cT5}"], u),
        _e_im_cT5(p["b->see::Im{cT5}"], u),
        /* b->smumu */
        _mu_re_c9(p["b->smumu::Re{c9}"], u),
        _mu_im_c9(p["b->smumu::Im{c9}"], u),
        _mu_re_c10(p["b->smumu::Re{c10}"], u),
        _mu_im_c10(p["b->smumu::Im{c10}"], u),
        _mu_re_c9prime(p["b->smumu::Re{c9'}"], u),
        _mu_im_c9prime(p["b->smumu::Im{c9'}"], u),
        _mu_re_c10prime(p["b->smumu::Re{c10'}"], u),
        _mu_im_c10prime(p["b->smumu::Im{c10'}"], u),
        _mu_re_cS(p["b->smumu::Re{cS}"], u),
        _mu_im_cS(p["b->smumu::Im{cS}"], u),
        _mu_re_cSprime(p["b->smumu::Re{cS'}"], u),
        _mu_im_cSprime(p["b->smumu::Im{cS'}"], u),
        _mu_re_cP(p["b->smumu::Re{cP}"], u),
        _mu_im_cP(p["b->smumu::Im{cP}"], u),
        _mu_re_cPprime(p["b->smumu::Re{cP'}"], u),
        _mu_im_cPprime(p["b->smumu::Im{cP'}"], u),
        _mu_re_cT(p["b->smumu::Re{cT}"], u),
        _mu_im_cT(p["b->smumu::Im{cT}"], u),
        _mu_re_cT5(p["b->smumu::Re{cT5}"], u),
        _mu_im_cT5(p["b->smumu::Im{cT5}"], u),


        /* functions for b->sgamma */
        _c7(std::bind(&wcimplementation::cartesian, _re_c7, _im_c7)),
        _c7prime(std::bind(&wcimplementation::cartesian, _re_c7prime, _im_c7prime)),

        /* functions for b->see */
        _e_c9(std::bind(&wcimplementation::cartesian, _e_re_c9, _e_im_c9)),
        _e_c10(std::bind(&wcimplementation::cartesian, _e_re_c10, _e_im_c10)),
        _e_c9prime(std::bind(&wcimplementation::cartesian, _e_re_c9prime, _e_im_c9prime)),
        _e_c10prime(std::bind(&wcimplementation::cartesian, _e_re_c10prime, _e_im_c10prime)),
        _e_cS(std::bind(&wcimplementation::cartesian, _e_re_cS, _e_im_cS)),
        _e_cSprime(std::bind(&wcimplementation::cartesian, _e_re_cSprime, _e_im_cSprime)),
        _e_cP(std::bind(&wcimplementation::cartesian, _e_re_cP, _e_im_cP)),
        _e_cPprime(std::bind(&wcimplementation::cartesian, _e_re_cPprime, _e_im_cPprime)),
        _e_cT(std::bind(&wcimplementation::cartesian, _e_re_cT, _e_im_cT)),
        _e_cT5(std::bind(&wcimplementation::cartesian, _e_re_cT5, _e_im_cT5)),

        /* functions for b->smumu */
        _mu_c9(std::bind(&wcimplementation::cartesian, _mu_re_c9, _mu_im_c9)),
        _mu_c10(std::bind(&wcimplementation::cartesian, _mu_re_c10, _mu_im_c10)),
        _mu_c9prime(std::bind(&wcimplementation::cartesian, _mu_re_c9prime, _mu_im_c9prime)),
        _mu_c10prime(std::bind(&wcimplementation::cartesian, _mu_re_c10prime, _mu_im_c10prime)),
        _mu_cS(std::bind(&wcimplementation::cartesian, _mu_re_cS, _mu_im_cS)),
        _mu_cSprime(std::bind(&wcimplementation::cartesian, _mu_re_cSprime, _mu_im_cSprime)),
        _mu_cP(std::bind(&wcimplementation::cartesian, _mu_re_cP, _mu_im_cP)),
        _mu_cPprime(std::bind(&wcimplementation::cartesian, _mu_re_cPprime, _mu_im_cPprime)),
        _mu_cT(std::bind(&wcimplementation::cartesian, _mu_re_cT, _mu_im_cT)),
        _mu_cT5(std::bind(&wcimplementation::cartesian, _mu_re_cT5, _mu_im_cT5))
    {
    }

    WilsonCoefficients<BToS>
    WilsonScanComponent<components::DeltaBS1>::wilson_coefficients_b_to_s(const double & /*mu*/, const LeptonFlavor & lepton_flavor, const bool & cp_conjugate) const
    {
        std::function<complex<double>()> c9, c9prime;
        std::function<complex<double>()> c10, c10prime;
        std::function<complex<double>()> cS, cSprime;
        std::function<complex<double>()> cP, cPprime;
        std::function<complex<double>()> cT, cT5;

        if (LeptonFlavor::electron == lepton_flavor)
        {
            c9       = _e_c9;
            c9prime  = _e_c9prime;
            c10      = _e_c10;
            c10prime = _e_c10prime;
            cS       = _e_cS;
            cSprime  = _e_cSprime;
            cP       = _e_cP;
            cPprime  = _e_cPprime;
            cT       = _e_cT;
            cT5      = _e_cT5;
        }
        else if (LeptonFlavor::muon == lepton_flavor)
        {
            c9       = _mu_c9;
            c9prime  = _mu_c9prime;
            c10      = _mu_c10;
            c10prime = _mu_c10prime;
            cS       = _mu_cS;
            cSprime  = _mu_cSprime;
            cP       = _mu_cP;
            cPprime  = _mu_cPprime;
            cT       = _mu_cT;
            cT5      = _mu_cT5;
        }
        else
        {
            throw InternalError("WilsonScan presently only implements 'e' and 'mu' lepton flavors");
        }

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

        complex<double>          a_s = alpha_s / 4.0 / M_PI;
        WilsonCoefficients<BToS> result;
        result._sm_like_coefficients = std::array<std::complex<double>, 15>{
            { _c1(), _c2(), _c3(), _c4(), _c5(), _c6(), 0.0, 0.0, 0.0, 0.0, 0.0, a_s * _c7(), a_s * _c8(), a_s * c9(), a_s * c10() }
        };
        result._primed_coefficients = std::array<std::complex<double>, 15>{
            { /* we only consider c7', c8', c9' and c10' */
              0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0,
             0.0, a_s * _c7prime(),
             a_s * _c8prime(),
             a_s * c9prime(),
             a_s * c10prime() }
        };
        result._scalar_tensor_coefficients = std::array<std::complex<double>, 6>{
            { cS(), cSprime(), cP(), cPprime(), cT(), cT5() }
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

    ConstrainedWilsonScanComponent::ConstrainedWilsonScanComponent(const Parameters & p, const Options & o, ParameterUser & u) :
        WilsonScanComponent<components::DeltaBS1>(p, o, u)
    {
        /* b->smee */
        _e_cT      = std::bind(&wcimplementation::zero);
        _e_cT5     = std::bind(&wcimplementation::zero);
        _e_cP      = std::bind(&wcimplementation::cartesian_negative, _e_re_cS, _e_im_cS);
        _e_cPprime = std::bind(&wcimplementation::cartesian, _e_re_cSprime, _e_im_cSprime);

        u.drop(_e_re_cP.id());
        u.drop(_e_im_cP.id());
        u.drop(_e_re_cPprime.id());
        u.drop(_e_im_cPprime.id());
        u.drop(_e_re_cT.id());
        u.drop(_e_im_cT.id());
        u.drop(_e_re_cT5.id());
        u.drop(_e_im_cT5.id());

        /* b->smumu */
        _mu_cT      = std::bind(&wcimplementation::zero);
        _mu_cT5     = std::bind(&wcimplementation::zero);
        _mu_cP      = std::bind(&wcimplementation::cartesian_negative, _mu_re_cS, _mu_im_cS);
        _mu_cPprime = std::bind(&wcimplementation::cartesian, _mu_re_cSprime, _mu_im_cSprime);

        u.drop(_mu_re_cP.id());
        u.drop(_mu_im_cP.id());
        u.drop(_mu_re_cPprime.id());
        u.drop(_mu_im_cPprime.id());
        u.drop(_mu_re_cT.id());
        u.drop(_mu_im_cT.id());
        u.drop(_mu_re_cT5.id());
        u.drop(_mu_im_cT5.id());
    }

    WilsonScanModel::WilsonScanModel(const Parameters & parameters, const Options & options) :
        CKMScanComponent(parameters, options, *this),
        SMComponent<components::QCD>(parameters, *this),
        // Charged-current semileptonic sectors (Delta S = 1)
        WilsonScanComponent<components::WET::USLNu>(parameters, options, *this),
        // Charged-current semileptonic sectors (Delta C = 1)
        WilsonScanComponent<components::WET::DCNuL>(parameters, options, *this),
        WilsonScanComponent<components::WET::SCNuL>(parameters, options, *this),
        // Charged-current semileptonic sectors (Delta B = 1)
        WilsonScanComponent<components::WET::UBLNu>(parameters, options, *this),
        WilsonScanComponent<components::WET::CBLNu>(parameters, options, *this),
        // Neutral-current semileptonic sectors (Delta B = 1)
        WilsonScanComponent<components::WET::SBNuNu>(parameters, options, *this),
        // Hadronic sectors (Delta B = 1)
        WilsonScanComponent<components::WET::DBCU>(parameters, options, *this),
        WilsonScanComponent<components::WET::SBCU>(parameters, options, *this),
        // Hadronic sectors (Delta B = 2)
        WilsonScanComponent<components::WET::SBSB>(parameters, options, *this),
        // Old-style WET sectors
        WilsonScanComponent<components::DeltaBS1>(parameters, options, *this)
    {
    }

    WilsonScanModel::~WilsonScanModel() {}

    std::shared_ptr<Model>
    WilsonScanModel::make(const Parameters & parameters, const Options & options)
    {
        return std::shared_ptr<Model>(new WilsonScanModel(parameters, options));
    }

    ConstrainedWilsonScanModel::ConstrainedWilsonScanModel(const Parameters & parameters, const Options & options) :
        CKMScanComponent(parameters, options, *this),
        SMComponent<components::QCD>(parameters, *this),
        // Charged-current semileptonic sectors (Delta S = 1)
        WilsonScanComponent<components::WET::USLNu>(parameters, options, *this),
        // Charged-current semileptonic sectors (Delta C = 1)
        WilsonScanComponent<components::WET::DCNuL>(parameters, options, *this),
        WilsonScanComponent<components::WET::SCNuL>(parameters, options, *this),
        // Charged-current semileptonic sectors (Delta B = 1)
        WilsonScanComponent<components::WET::UBLNu>(parameters, options, *this),
        WilsonScanComponent<components::WET::CBLNu>(parameters, options, *this),
        // Neutral-current semileptonic sectors (Delta B = 1)
        WilsonScanComponent<components::WET::SBNuNu>(parameters, options, *this),
        // Hadronic sectors (Delta B = 1)
        WilsonScanComponent<components::WET::DBCU>(parameters, options, *this),
        WilsonScanComponent<components::WET::SBCU>(parameters, options, *this),
        // Hadronic sectors (Delta B = 2)
        WilsonScanComponent<components::WET::SBSB>(parameters, options, *this),
        // Old-style WET sectors
        ConstrainedWilsonScanComponent(parameters, options, *this)
    {
    }

    ConstrainedWilsonScanModel::~ConstrainedWilsonScanModel() {}

    std::shared_ptr<Model>
    ConstrainedWilsonScanModel::make(const Parameters & parameters, const Options & options)
    {
        return std::shared_ptr<Model>(new ConstrainedWilsonScanModel(parameters, options));
    }
} // namespace eos
