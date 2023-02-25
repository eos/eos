/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014, 2015 Danny van Dyk
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

#include <eos/models/ckm.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/log.hh>
#include <eos/maths/matrix.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

#include <cmath>

namespace eos
{
    using std::sqrt;

    CKMScanComponent::CKMScanComponent(const Parameters & p, const Options &, ParameterUser & u) :
        _v_ud_abs__ckm(p["CKM::abs(V_ud)"], u),
        _v_ud_arg__ckm(p["CKM::arg(V_ud)"], u),
        _v_us_abs__ckm(p["CKM::abs(V_us)"], u),
        _v_us_arg__ckm(p["CKM::arg(V_us)"], u),
        _v_ub_abs__ckm(p["CKM::abs(V_ub)"], u),
        _v_ub_arg__ckm(p["CKM::arg(V_ub)"], u),
        _v_cd_abs__ckm(p["CKM::abs(V_cd)"], u),
        _v_cd_arg__ckm(p["CKM::arg(V_cd)"], u),
        _v_cs_abs__ckm(p["CKM::abs(V_cs)"], u),
        _v_cs_arg__ckm(p["CKM::arg(V_cs)"], u),
        _v_cb_abs__ckm(p["CKM::abs(V_cb)"], u),
        _v_cb_arg__ckm(p["CKM::arg(V_cb)"], u),
        _v_td_abs__ckm(p["CKM::abs(V_td)"], u),
        _v_td_arg__ckm(p["CKM::arg(V_td)"], u),
        _v_ts_abs__ckm(p["CKM::abs(V_ts)"], u),
        _v_ts_arg__ckm(p["CKM::arg(V_ts)"], u),
        _v_tb_abs__ckm(p["CKM::abs(V_tb)"], u),
        _v_tb_arg__ckm(p["CKM::arg(V_tb)"], u)
    {
    }

    /* CKM matrix elements */
    complex<double>
    CKMScanComponent::ckm_ud() const
    {
        return std::polar(_v_ud_abs__ckm(), _v_ud_arg__ckm());
    }

    complex<double>
    CKMScanComponent::ckm_us() const
    {
        return std::polar(_v_us_abs__ckm(), _v_us_arg__ckm());
    }

    complex<double>
    CKMScanComponent::ckm_ub() const
    {
        return std::polar(_v_ub_abs__ckm(), _v_ub_arg__ckm());
    }

    complex<double>
    CKMScanComponent::ckm_cd() const
    {
        return std::polar(_v_cd_abs__ckm(), _v_cd_arg__ckm());
    }

    complex<double>
    CKMScanComponent::ckm_cs() const
    {
        return std::polar(_v_cs_abs__ckm(), _v_cs_arg__ckm());
    }

    complex<double>
    CKMScanComponent::ckm_cb() const
    {
        return std::polar(_v_cb_abs__ckm(), _v_cb_arg__ckm());
    }

    complex<double>
    CKMScanComponent::ckm_td() const
    {
        return std::polar(_v_td_abs__ckm(), _v_td_arg__ckm());
    }

    complex<double>
    CKMScanComponent::ckm_ts() const
    {
        return std::polar(_v_ts_abs__ckm(), _v_ts_arg__ckm());
    }

    complex<double>
    CKMScanComponent::ckm_tb() const
    {
        return std::polar(_v_tb_abs__ckm(), _v_tb_arg__ckm());
    }

    CKMScanModel::CKMScanModel(const Parameters & parameters, const Options & options) :
        CKMScanComponent(parameters, options, *this),
        SMComponent<components::QCD>(parameters, *this),
        SMComponent<components::WET::SBSB>(parameters, *this),
        SMComponent<components::DeltaBS1>(parameters, *this),
        SMComponent<components::WET::CBLNu>(parameters, *this),
        SMComponent<components::WET::UBLNu>(parameters, *this),
        SMComponent<components::WET::SBNuNu>(parameters, *this),
        SMComponent<components::WET::SBCU>(parameters, *this),
        SMComponent<components::WET::DBCU>(parameters, *this)
    {
    }

    CKMScanModel::~CKMScanModel()
    {
    }

    std::shared_ptr<Model>
    CKMScanModel::make(const Parameters & parameters, const Options & options)
    {
        return std::shared_ptr<Model>(new CKMScanModel(parameters, options));
    }
}
