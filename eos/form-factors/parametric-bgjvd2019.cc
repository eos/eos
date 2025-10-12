/* vim: set sw=4 sts=4 tw=120 et foldmethod=syntax : */

/*
 * Copyright (c) 2022-2025 Danny van Dyk
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

#include <eos/form-factors/parametric-bgjvd2019-impl.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    using namespace std::literals::string_literals;

    std::string
    HQETFormFactorBase::_sslp_prefix(const std::string & prefix)
    {
        if ("B(*)->D(*)" == prefix)
            return prefix;

        if (_opt_sslp_limit.value())
            return "B(*)->D(*)";

        return prefix;
    }

    HQETFormFactorBase::HQETFormFactorBase(const Parameters & p, const Options & o, const std::string & prefix) :
        _model(Model::make("SM", p, o)),
        _mBar(p[prefix + "::mBar@HQET"], *this),
        _a(p[prefix + "::a@HQET"], *this),
        _opt_lp_model(o, "model-lp"_ok, { "power-series", "exponential" }, "power-series"),
        _opt_lp_zorder(o, option_specifications, "z-order-lp"_ok),
        _enable_lp_z3(1.0 ? _opt_lp_zorder.value() >= 3 : 0.0),
        _enable_lp_z4(1.0 ? _opt_lp_zorder.value() >= 4 : 0.0),
        _enable_lp_z5(1.0 ? _opt_lp_zorder.value() >= 5 : 0.0),
        _opt_slp_zorder(o, option_specifications, "z-order-slp"_ok),
        _enable_slp_z2(1.0 ? _opt_slp_zorder.value() >= 2 : 0.0),
        _opt_sslp_zorder(o, option_specifications, "z-order-sslp"_ok),
        _enable_sslp_z1(1.0 ? _opt_sslp_zorder.value() >= 1 : 0.0),
        _enable_sslp_z2(1.0 ? _opt_sslp_zorder.value() >= 2 : 0.0),
        _opt_sslp_limit(o, option_specifications, "SU3F-limit-sslp"_ok),
        _xipone(p[prefix + "::xi'(1)@HQET"], *this),
        _xippone(p[prefix + "::xi''(1)@HQET"], *this),
        _xipppone(p[prefix + "::xi'''(1)@HQET"], *this),
        _xippppone(p[prefix + "::xi''''(1)@HQET"], *this),
        _xipppppone(p[prefix + "::xi'''''(1)@HQET"], *this),
        _chi2one(p[prefix + "::chi_2(1)@HQET"], *this),
        _chi2pone(p[prefix + "::chi_2'(1)@HQET"], *this),
        _chi2ppone(p[prefix + "::chi_2''(1)@HQET"], *this),
        _chi3pone(p[prefix + "::chi_3'(1)@HQET"], *this),
        _chi3ppone(p[prefix + "::chi_3''(1)@HQET"], *this),
        _etaone(p[prefix + "::eta(1)@HQET"], *this),
        _etapone(p[prefix + "::eta'(1)@HQET"], *this),
        _etappone(p[prefix + "::eta''(1)@HQET"], *this),
        _l1one(p[_sslp_prefix(prefix) + "::l_1(1)@HQET"], *this),
        _l1pone(p[_sslp_prefix(prefix) + "::l_1'(1)@HQET"], *this),
        _l1ppone(p[_sslp_prefix(prefix) + "::l_1''(1)@HQET"], *this),
        _l2one(p[_sslp_prefix(prefix) + "::l_2(1)@HQET"], *this),
        _l2pone(p[_sslp_prefix(prefix) + "::l_2'(1)@HQET"], *this),
        _l2ppone(p[_sslp_prefix(prefix) + "::l_2''(1)@HQET"], *this),
        _l3one(p[_sslp_prefix(prefix) + "::l_3(1)@HQET"], *this),
        _l3pone(p[_sslp_prefix(prefix) + "::l_3'(1)@HQET"], *this),
        _l3ppone(p[_sslp_prefix(prefix) + "::l_3''(1)@HQET"], *this),
        _l4one(p[_sslp_prefix(prefix) + "::l_4(1)@HQET"], *this),
        _l4pone(p[_sslp_prefix(prefix) + "::l_4'(1)@HQET"], *this),
        _l4ppone(p[_sslp_prefix(prefix) + "::l_4''(1)@HQET"], *this),
        _l5one(p[_sslp_prefix(prefix) + "::l_5(1)@HQET"], *this),
        _l5pone(p[_sslp_prefix(prefix) + "::l_5'(1)@HQET"], *this),
        _l5ppone(p[_sslp_prefix(prefix) + "::l_5''(1)@HQET"], *this),
        _l6one(p[_sslp_prefix(prefix) + "::l_6(1)@HQET"], *this),
        _l6pone(p[_sslp_prefix(prefix) + "::l_6'(1)@HQET"], *this),
        _l6ppone(p[_sslp_prefix(prefix) + "::l_6''(1)@HQET"], *this)
    {
        using std::placeholders::_1;

        if (_opt_lp_model.value() == "exponential")
        {
            _xi = [=, this](const double & q2) -> double
            {
                return _xi_exponential(q2);
            };
        }
        else
        {
            _xi = [=, this](const double & q2) -> double
            {
                return _xi_power_series(q2);
            };
        }
    }

    HQETFormFactorBase::~HQETFormFactorBase() = default;

    const std::set<ReferenceName>
    HQETFormFactorBase::references
    {
        "BLPR:2017A"_rn,
        "JS:2018A"_rn,
        "BJvD:2019A"_rn,
        "BGJvD:2019A"_rn
    };

    const std::vector<OptionSpecification>
    HQETFormFactorBase::option_specifications
    {
        { "z-order-lp"_ok,      { "2"s, "3"s, "4"s, "5"s }, "3"s     },
        { "z-order-slp"_ok,     { "1"s, "2"s },             "2"s     },
        { "z-order-sslp"_ok,    { "0"s, "1"s, "2"s },       "1"s     },
        { "SU3F-limit-sslp"_ok, { "true"s, "false"s },      "false"s }
    };

    // uses a power series ansatz
    double
    HQETFormFactorBase::_xi_power_series(const double & q2) const
    {
        const double a = _a(), a2 = a * a, a3 = a * a2, a4 = a2 * a2, a5 = a3 * a2;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_z(q2) - z_0);
        const double z2   =  z *  z;
        const double z3   = z2 *  z * _enable_lp_z3;
        const double z4   = z2 * z2 * _enable_lp_z4;
        const double z5   = z3 * z2 * _enable_lp_z5;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                          + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2
                          + (2.0 +       a) * power_of<4>(1.0 + a) / (2.0 * a3) * z3
                          + (5.0 + 3.0 * a) * power_of<5>(1.0 + a) / (8.0 * a4) * z4
                          + (3.0 + 2.0 * a) * power_of<6>(1.0 + a) / (8.0 * a5) * z5;

        const double wm12 =   4.0                  * power_of<4>(1.0 + a) / a2         * z2
                          + ( 6.0 +  2.0 * a     ) * power_of<5>(1.0 + a) / a3         * z3
                          + (25.0 + 14.0 * a + a2) * power_of<6>(1.0 + a) / (4.0 * a4) * z4
                          + (11.0 +  8.0 * a + a2) * power_of<7>(1.0 + a) / (2.0 * a5) * z5;

        const double wm13 =   8.0                  * power_of<6>(1.0 + a) / a3         * z3
                          + (18.0 +  6.0 * a     ) * power_of<7>(1.0 + a) / a4         * z4
                          + (51.0 + 30.0 * a + a2) * power_of<8>(1.0 + a) / (2.0 * a5) * z5;

        const double wm14 =  16.0             * power_of<8>(1.0 + a) / a4 * z4
                          + (48.0 + 16.0 * a) * power_of<9>(1.0 + a) / a5 * z5;

        const double wm15 = 32.0 * power_of<5>(1.0 + a) / a5 * z5;

        return 1.0
            + _xipone             * wm11
            + _xippone    / 2.0   * wm12
            + _xipppone   / 6.0   * wm13
            + _xippppone  / 24.0  * wm14
            + _xipppppone / 120.0 * wm15;
    }

    // uses an exponential ansatz and expands in (w-1) first, then in z*
    double
    HQETFormFactorBase::_xi_exponential(const double & q2) const
    {
        const double a = _a(), a2 = a * a, a3 = a * a2, a4 = a2 * a2, a5 = a3 * a2;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_z(q2) - z_0);
        const double z2   =  z *  z;
        const double z3   = z2 *  z * _enable_lp_z3;
        const double z4   = z2 * z2 * _enable_lp_z4;
        const double z5   = z3 * z2 * _enable_lp_z5;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                          + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2
                          + (2.0 +       a) * power_of<4>(1.0 + a) / (2.0 * a3) * z3
                          + (5.0 + 3.0 * a) * power_of<5>(1.0 + a) / (8.0 * a4) * z4
                          + (3.0 + 2.0 * a) * power_of<6>(1.0 + a) / (8.0 * a5) * z5;

        const double wm12 =   4.0                  * power_of<4>(1.0 + a) / a2         * z2
                          + ( 6.0 +  2.0 * a     ) * power_of<5>(1.0 + a) / a3         * z3
                          + (25.0 + 14.0 * a + a2) * power_of<6>(1.0 + a) / (4.0 * a4) * z4
                          + (11.0 +  8.0 * a + a2) * power_of<7>(1.0 + a) / (2.0 * a5) * z5;

        const double wm13 =   8.0                  * power_of<6>(1.0 + a) / a3         * z3
                          + (18.0 +  6.0 * a     ) * power_of<7>(1.0 + a) / a4         * z4
                          + (51.0 + 30.0 * a + a2) * power_of<8>(1.0 + a) / (2.0 * a5) * z5;

        const double wm14 =  16.0             * power_of<8>(1.0 + a) / a4 * z4
                          + (48.0 + 16.0 * a) * power_of<9>(1.0 + a) / a5 * z5;

        const double wm15 = 32.0 * power_of<5>(1.0 + a) / a5 * z5;

        return (1.0
            + _xipone              * wm11
            - _xipone              * wm12
            + _xipone * 2.0 /  3.0 * wm13
            - _xipone       /  3.0 * wm14
            + _xipone * 2.0 / 15.0 * wm15)
            * (1.0 + _xippone      * wm11);
    }

    double
    HQETFormFactorBase::_chi2(const double & q2) const
    {
        const double a = _a(), a2 = a * a;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_z(q2) - z_0);
        const double z2   =  z * z * _enable_slp_z2;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                            + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2;

        const double wm12 =   4.0           * power_of<4>(1.0 + a) / a2         * z2;

        return _chi2one + _chi2pone * wm11 + _chi2ppone / 2.0 * wm12;
    }

    double
    HQETFormFactorBase::_chi3(const double & q2) const
    {
        const double a = _a(), a2 = a * a;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_z(q2) - z_0);
        const double z2   =  z * z * _enable_slp_z2;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                            + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2;

        const double wm12 =   4.0           * power_of<4>(1.0 + a) / a2         * z2;

        return 0.0 + _chi3pone * wm11 + _chi3ppone / 2.0 * wm12;
    }

    double
    HQETFormFactorBase::_eta(const double & q2) const
    {
        const double a = _a(), a2 = a * a;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_z(q2) - z_0);
        const double z2   =  z * z * _enable_slp_z2;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                            + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2;

        const double wm12 =   4.0           * power_of<4>(1.0 + a) / a2         * z2;

        return _etaone + _etapone * wm11 + _etappone / 2.0 * wm12;
    }

    /* Power corrections */
    double
    HQETFormFactorBase::_l1(const double & w) const
    {
        const double a = _a(), a2 = a * a;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
        const double z2   =  z * z * _enable_sslp_z2;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                            + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2;

        const double wm12 =   4.0           * power_of<4>(1.0 + a) / a2         * z2;

        return _l1one + _l1pone * wm11 + _l1ppone / 2.0 * wm12;
    }

    double
    HQETFormFactorBase::_l2(const double & w) const
    {
        const double a = _a(), a2 = a * a;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
        const double z2   =  z * z * _enable_sslp_z2;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                            + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2;

        const double wm12 =   4.0           * power_of<4>(1.0 + a) / a2         * z2;

        return _l2one + _l2pone * wm11 + _l2ppone / 2.0 * wm12;
    }

    double
    HQETFormFactorBase::_l3(const double & w) const
    {
        const double a = _a(), a2 = a * a;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
        const double z2   =  z * z * _enable_sslp_z2;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                            + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2;

        const double wm12 =   4.0           * power_of<4>(1.0 + a) / a2         * z2;

        return _l3one + _l3pone * wm11 + _l3ppone / 2.0 * wm12;
    }

    double
    HQETFormFactorBase::_l4(const double & w) const
    {
        const double a = _a(), a2 = a * a;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
        const double z2   =  z * z * _enable_sslp_z2;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                            + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2;

        const double wm12 =   4.0           * power_of<4>(1.0 + a) / a2         * z2;

        return _l4one + _l4pone * wm11 + _l4ppone / 2.0 * wm12;
    }

    double
    HQETFormFactorBase::_l5(const double & w) const
    {
        const double a = _a(), a2 = a * a;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
        const double z2   =  z * z * _enable_sslp_z2;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                            + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2;

        const double wm12 =   4.0           * power_of<4>(1.0 + a) / a2         * z2;

        return _l5one + _l5pone * wm11 + _l5ppone / 2.0 * wm12;
    }

    double
    HQETFormFactorBase::_l6(const double & w) const
    {
        const double a = _a(), a2 = a * a;

        // expansion in z around z_0
        const double  z_0 = (1.0 - a) / (1.0 + a);
        const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
        const double z2   =  z * z * _enable_sslp_z2;

        const double wm11 =  2.0            * power_of<2>(1.0 + a) / a          * z
                            + (3.0 +       a) * power_of<3>(1.0 + a) / (2.0 * a2) * z2;

        const double wm12 =   4.0           * power_of<4>(1.0 + a) / a2         * z2;

        return _l6one + _l6pone * wm11 + _l6ppone / 2.0 * wm12;
    }

    /* Wilson Coefficients */

    double
    HQETFormFactorBase::_CS(const double & w, const double & z) const
    {
        const double z2  = z * z;
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = 2.0 * z * (w - wz) * _Omega(w, z);
        result -= (w - 1.0) * (z + 1.0) * (z + 1.0) * _r(w);
        result += (z2 - 1.0) * lnz;

        return result / (3.0 * z * (w - wz));
    }

    double
    HQETFormFactorBase::_CP(const double & w, const double & z) const
    {
        const double z2  = z * z;
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = 2.0 * z * (w - wz) * _Omega(w, z);
        result -= (w + 1.0) * (z - 1.0) * (z - 1.0) * _r(w);
        result += (z2 - 1.0) * lnz;

        return result / (3.0 * z * (w - wz));
    }

    double
    HQETFormFactorBase::_CV1(const double & w, const double & z) const
    {
        const double z2  = z * z;
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = 2.0 * (w + 1.0) * ((3.0 * w - 1.0) * z - z2 - 1.0) * _r(w);
        result += (12.0 * z * (wz - w) - (z2 - 1.0) * lnz);
        result += 4.0 * z * (w - wz) * _Omega(w, z);

        return result / (6.0 * z * (w - wz));
    }

    double
    HQETFormFactorBase::_CV2(const double & w, const double & z) const
    {
        const double z2  = z * z, z3 = z2 * z;
        const double w2  = w * w;
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = ((4.0 * w2 + 2.0 * w) * z2 - (2.0 * w2 + 5.0 * w - 1.0) * z - (1.0 + w) * z3 + 2.0) * _r(w);
        result += z * (2.0 * (z - 1.0) * (wz - w) + (z2 - (4.0 * w - 2.0) * z + (-2.0 * w + 3)) * lnz);

        return -1.0 * result / (6.0 * z2 * power_of<2>(w - wz));
    }

    double
    HQETFormFactorBase::_CV3(const double & w, const double & z) const
    {
        const double z2  = z * z, z3 = z2 * z;
        const double w2  = w * w;
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = (-2.0 * z3 + (2.0 * w2 + 5.0 * w - 1.0) * z2 - (4.0 * w2 + 2.0 * w) * z + w + 1.0) * _r(w);
        result += 2.0 * z * (z - 1.0) * (wz - w) + ((-2.0 * w + 3.0) * z2 + (-4.0 * w + 2.0) * z + 1.0) * lnz;

        return +1.0 * result / (6.0 * z * power_of<2>(w - wz));
    }

    double
    HQETFormFactorBase::_CA1(const double & w, const double & z) const
    {
        const double z2  = z * z;
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = 2.0 * (w - 1.0) * ((3.0 * w + 1.0) * z - z2 - 1.0) * _r(w);
        result += (12.0 * z * (wz - w) - (z2 - 1.0) * lnz);
        result += 4.0 * z * (w - wz) * _Omega(w, z);

        return result / (6.0 * z * (w - wz));
    }

    double
    HQETFormFactorBase::_CA2(const double & w, const double & z) const
    {
        const double z2  = z * z, z3 = z2 * z;
        const double w2  = w * w;
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = ((4.0 * w2 - 2.0 * w) * z2 + (2.0 * w2 - 5.0 * w - 1.0) * z + (1.0 - w) * z3 + 2.0) * _r(w);
        result += z * (2.0 * (z + 1.0) * (wz - w) + (z2 - (4.0 * w + 2.0) * z + (2.0 * w + 3)) * lnz);

        return -1.0 * result / (6.0 * z2 * power_of<2>(w - wz));
    }

    double
    HQETFormFactorBase::_CA3(const double & w, const double & z) const
    {
        const double z2  = z * z, z3 = z2 * z;
        const double w2  = w * w;
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = (2.0 * z3 + (2.0 * w2 - 5.0 * w - 1.0) * z2 + (4.0 * w2 - 2.0 * w) * z - w + 1.0) * _r(w);
        result += 2.0 * z * (z + 1.0) * (wz - w) - ((2.0 * w + 3.0) * z2 - (4.0 * w + 2.0) * z + 1.0) * lnz;

        return +1.0 * result / (6.0 * z * power_of<2>(w - wz));
    }

    double
    HQETFormFactorBase::_CT1(const double & w, const double & z) const
    {
        const double z2  = z * z;
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = (w - 1.0) * ((4.0 * w + 2.0) * z - z2 - 1.0) * _r(w);
        result += 6.0 * z * (wz - w) - (z2 - 1.0) * lnz;
        result += 2.0 * z * (w - wz) * _Omega(w, z);

        return +1.0 / (3.0 * z * (w - wz)) * result;
    }

    double
    HQETFormFactorBase::_CT2(const double & w, const double & z) const
    {
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = (1.0 - w * z) * _r(w) + z * lnz;

        return +2.0 / (3.0 * z * (w - wz)) * result;
    }

    double
    HQETFormFactorBase::_CT3(const double & w, const double & z) const
    {
        const double wz  = _wz(z);
        const double lnz = std::log(z);

        double result = (w - z) * _r(w) + lnz;

        return +2.0 / (3.0 * (w - wz)) * result;
    }

    // P->P
    template class HQETFormFactors<BToD,   PToP>;
    template class HQETFormFactors<BsToDs, PToP>;

    // P->V
    template class HQETFormFactors<BToDstar,   PToV>;
    template class HQETFormFactors<BsToDsstar, PToV>;

    // P->V
    template class HQETFormFactors<BstarToD,  VToP>;

    // V->V
    template class HQETFormFactors<BstarToDstar,  VToV>;
}
