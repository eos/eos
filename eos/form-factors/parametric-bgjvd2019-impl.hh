/* vim: set sw=4 sts=4 tw=120 et foldmethod=syntax : */

/*
 * Copyright (c) 2018-2024 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGJvD2019_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGJvD2019_IMPL_HH 1

#include <eos/form-factors/parametric-bgjvd2019.hh>
#include <eos/utils/kinematic.hh>
#include <eos/models/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/maths/polylog.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/log.hh>

#include <cmath>
#include <limits>

namespace eos
{
    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>::_w(const double & q2) const
    {
        const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
        const double m_P = this->_m_P(), m_P2 = power_of<2>(m_P);

        return (m_B2 + m_P2 - q2) / (2.0 * m_B * m_P);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>::_q2(const double & w) const
    {
        const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
        const double m_P = this->_m_P(), m_P2 = power_of<2>(m_P);

        return m_B2 + m_P2 - 2.0 * m_B * m_P * w;
    }

    /* HQET form factors h_i */

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>::_h_p(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;

        double result = 1.0 + as * (_CV1(w, z) + (w + 1.0) / 2.0 * (_CV2(w, z) + _CV3(w, z)));
        result += eps_c * (L1);
        result += eps_b * (L1);
        result += eps_c * eps_c * _l1(w);

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>::_h_m(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L4 = 2.0 * eta - 1.0;

        double result = (0.0 + as * (w + 1.0) / 2.0 * (_CV2(w, z) - _CV3(w, z)));
        result += eps_c * L4;
        result -= eps_b * L4;
        result += eps_c * eps_c * _l4(w);

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>:: _h_S(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
        const double L4 = 2.0 * eta - 1.0;

        double result = (1.0 + as * _CS(w, z));
        result += eps_c * (L1 - (w - 1.0) / (w + 1.0) * L4);
        result += eps_b * (L1 - (w - 1.0) / (w + 1.0) * L4);
        result += eps_c * eps_c * (_l1(w) - (w - 1.0) / (w + 1.0) * _l4(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>::_h_T(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
        const double L4 = 2.0 * eta - 1.0;

        double result = 1.0 + as * (_CT1(w, z) - _CT2(w, z) + _CT3(w, z));
        result += eps_c * (L1 - L4);
        result += eps_b * (L1 - L4);
        result += eps_c * eps_c * (_l1(w) - _l4(w));

        return result * xi;
    }


    template <typename Process_>
    HQETFormFactors<Process_, PToP>::HQETFormFactors(const Parameters & p, const Options & o) :
        HQETFormFactorBase(p, o, Process_::hqe_prefix),
        _m_B(p[Process_::name_B], *static_cast<ParameterUser *>(this)),
        _m_P(p[Process_::name_P], *static_cast<ParameterUser *>(this))
    {
        static const Log::OneTimeMessage message_HQET_FFs_PToP
        (
            std::string("HQETFormFactors<") + Process_::label + ",PToP>",
            ll_warning,
            "This form factor parametrization is not a general one and requires careful attention. "
            "By default, it returns zeros for all form factors."
        );
    }

    template <typename Process_>
    HQETFormFactors<Process_, PToP>::~HQETFormFactors() = default;

    template <typename Process_>
    FormFactors<PToP> *
    HQETFormFactors<Process_, PToP>::make(const Parameters & parameters, const Options & options)
    {
        return new HQETFormFactors<Process_, PToP>(parameters, options);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>::f_p(const double & q2) const
    {
        const double r = _m_P / _m_B;

        // cf. [FKKM2008], eq. (22)
        return 1.0 / (2.0 * sqrt(r)) * ((1.0 + r) * _h_p(q2) - (1.0 - r) * _h_m(q2));
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>::f_m(const double & q2) const
    {
        const double r = _m_P / _m_B;

        // cf. [FKKM2008], eq. (22)
        return 1.0 / (2.0 * sqrt(r)) * ((1.0 + r) * _h_m(q2) - (1.0 - r) * _h_p(q2));
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>::f_0(const double & q2) const
    {
        // We do not use the relation between f_0 and the (scale-dependent) h_S.
        return f_p(q2) + q2 / (_m_B * _m_B - _m_P * _m_P) * f_m(q2);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>::f_t(const double & q2) const
    {
        const double r = _m_P / _m_B;

        // cf. [BJvD2019], eq. (A7)
        return (1.0 + r) / (2.0 * sqrt(r)) * _h_T(q2);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToP>::f_plus_T(const double & q2) const
    {
        return f_t(q2) * q2 / _m_B / (_m_B + _m_P);
    }

    template <typename Process_>
    Diagnostics
    HQETFormFactors<Process_, PToP>::diagnostics() const
    {
        Diagnostics results;

        // Inputs
        {
            const double m_b = _m_b_pole();
            const double m_c = _m_c_pole();
            const double z   = m_c / m_b;
            const double wz  = _wz(z);

            results.add(Diagnostics::Entry{ z,  "z = m_c_pole / m_b_pole" });
            results.add(Diagnostics::Entry{ wz, "w_z"                     });
        }

        // Switches
        {
            results.add(Diagnostics::Entry{ _enable_lp_z3,  "enable LP  z^3 terms" });
            results.add(Diagnostics::Entry{ _enable_lp_z4,  "enable LP  z^4 terms" });
            results.add(Diagnostics::Entry{ _enable_lp_z5,  "enable LP  z^5 terms" });
            results.add(Diagnostics::Entry{ _enable_slp_z2, "enable SLP z^2 terms" });
        }

        // z
        {
            results.add(Diagnostics::Entry{ _z(_q2(1.10)), "z(w = 1.10)" });
            results.add(Diagnostics::Entry{ _z(_q2(1.05)), "z(w = 1.05)" });
            results.add(Diagnostics::Entry{ _z(_q2(1.00)), "z(w = 1.00)" });
        }

        // xi
        {
            results.add(Diagnostics::Entry{ _xi(_q2(2.10)), "xi(w = 2.10)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.60)), "xi(w = 1.60)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.10)), "xi(w = 1.10)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.05)), "xi(w = 1.05)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.00)), "xi(w = 1.00)" });
        }

        // chi2
        {
            results.add(Diagnostics::Entry{ _chi2(_q2(2.10)), "chi2(w = 2.10)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.60)), "chi2(w = 1.60)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.10)), "chi2(w = 1.10)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.05)), "chi2(w = 1.05)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.00)), "chi2(w = 1.00)" });
        }

        // chi3
        {
            results.add(Diagnostics::Entry{ _chi3(_q2(2.10)), "chi3(w = 2.10)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.60)), "chi3(w = 1.60)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.10)), "chi3(w = 1.10)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.05)), "chi3(w = 1.05)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.00)), "chi3(w = 1.00)" });
        }

        // eta
        {
            results.add(Diagnostics::Entry{ _eta(_q2(2.10)), "eta(w = 2.10)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.60)), "eta(w = 1.60)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.10)), "eta(w = 1.10)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.05)), "eta(w = 1.05)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.00)), "eta(w = 1.00)" });
        }

        // r(w)
        {
            results.add(Diagnostics::Entry{ _r(1.1),     "r(w = 1.1)"     });
            results.add(Diagnostics::Entry{ _r(1.0007),  "r(w = 1.0007)"  });
            results.add(Diagnostics::Entry{ _r(1.0001),  "r(w = 1.0001)"  });
            results.add(Diagnostics::Entry{ _r(1.00005), "r(w = 1.00005)" });
            results.add(Diagnostics::Entry{ _r(1.0),     "r(w = 1.0)"     });
        }

        // Omega(w, z = 0.25)
        {
            results.add(Diagnostics::Entry{ _Omega(1.1,     0.25), "Omega(w = 1.1,     z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0007,  0.25), "Omega(w = 1.0007,  z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0001,  0.25), "Omega(w = 1.0001,  z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.00005, 0.25), "Omega(w = 1.00005, z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0,     0.25), "Omega(w = 1.0,     z = 0.25)" });
        }

        // Omega(w, z = 0.20)
        {
            results.add(Diagnostics::Entry{ _Omega(1.1,     0.20), "Omega(w = 1.1,     z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0007,  0.20), "Omega(w = 1.0007,  z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0001,  0.20), "Omega(w = 1.0001,  z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.00005, 0.20), "Omega(w = 1.00005, z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0,     0.20), "Omega(w = 1.0,     z = 0.20)" });
        }

        // WCs at w = 1.2, z = 0.20
        {
            results.add(Diagnostics::Entry{ _CS( 1.2, 0.20), "C_{S  }(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CP( 1.2, 0.20), "C_{P  }(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV1(1.2, 0.20), "C_{V_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV2(1.2, 0.20), "C_{V_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV3(1.2, 0.20), "C_{V_3}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA1(1.2, 0.20), "C_{A_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA2(1.2, 0.20), "C_{A_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA3(1.2, 0.20), "C_{A_3}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT1(1.2, 0.20), "C_{T_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT2(1.2, 0.20), "C_{T_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT3(1.2, 0.20), "C_{T_3}(w = 1.2, z = 0.20)" });
        }

        // WCs at w = 1.0, z = 0.25
        {
            results.add(Diagnostics::Entry{ _CS( 1.0, 0.25), "C_{S  }(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CP( 1.0, 0.25), "C_{P  }(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV1(1.0, 0.25), "C_{V_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV2(1.0, 0.25), "C_{V_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV3(1.0, 0.25), "C_{V_3}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA1(1.0, 0.25), "C_{A_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA2(1.0, 0.25), "C_{A_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA3(1.0, 0.25), "C_{A_3}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT1(1.0, 0.25), "C_{T_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT2(1.0, 0.25), "C_{T_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT3(1.0, 0.25), "C_{T_3}(w = 1.0, z = 0.25)" });
        }

        // HQET definition of the form factors
        {
            results.add(Diagnostics::Entry{ _h_p(_q2(1.4)), "h_+(w = 1.4)" });
            results.add(Diagnostics::Entry{ _h_m(_q2(1.4)), "h_-(w = 1.4)" });
            results.add(Diagnostics::Entry{ _h_T(_q2(1.4)), "h_T(w = 1.4)" });

            results.add(Diagnostics::Entry{ _h_p(_q2(1.2)), "h_+(w = 1.2)" });
            results.add(Diagnostics::Entry{ _h_m(_q2(1.2)), "h_-(w = 1.2)" });
            results.add(Diagnostics::Entry{ _h_T(_q2(1.2)), "h_T(w = 1.2)" });

            results.add(Diagnostics::Entry{ _h_p(_q2(1.0)), "h_+(w = 1.0)" });
            results.add(Diagnostics::Entry{ _h_m(_q2(1.0)), "h_-(w = 1.0)" });
            results.add(Diagnostics::Entry{ _h_T(_q2(1.0)), "h_T(w = 1.0)" });
        }

        return results;
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    HQETFormFactors<Process_, PToP>::begin_options()
    {
        return option_specifications.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    HQETFormFactors<Process_, PToP>::end_options()
    {
        return option_specifications.cend();
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::_w(const double & q2) const
    {
        const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
        const double m_V = this->_m_V(), m_V2 = power_of<2>(m_V);

        return (m_B2 + m_V2 - q2) / (2.0 * m_B * m_V);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::_q2(const double & w) const
    {
        const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
        const double m_V = this->_m_V(), m_V2 = power_of<2>(m_V);

        return m_B2 + m_V2 - 2.0 * m_B * m_V * w;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::_h_a1(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
        const double L2 = -4.0 * chi3;
        const double L4 = 2.0 * eta - 1.0;
        const double L5 = -1.0;

        double result = (1.0 + as * _CA1(w, z));
        result += eps_c * (L2 - L5 * (w - 1.0) / (w + 1.0));
        result += eps_b * (L1 - L4 * (w - 1.0) / (w + 1.0));
        result += eps_c * eps_c * (_l2(w) - (w - 1.0) / (w + 1.0) * _l5(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::_h_a2(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);

        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L3 = 4.0 * chi2;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (0.0 + as * _CA2(w, z));
        result += eps_c * (L3 + L6);
        result += eps_c * eps_c * (_l3(w) + _l6(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::_h_a3(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
        const double L2 = -4.0 * chi3;
        const double L3 = 4.0 * chi2;
        const double L4 = 2.0 * eta - 1.0;
        const double L5 = -1.0;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (1.0 + as * (_CA1(w, z) +_CA3(w, z)));
        result += eps_c * (L2 - L3 + L6 - L5);
        result += eps_b * (L1 - L4);
        result += eps_c * eps_c * (_l2(w) - _l3(w) + _l6(w) - _l5(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::_h_v(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
        const double L2 = -4.0 * chi3;
        const double L4 = 2.0 * eta - 1.0;
        const double L5 = -1.0;

        double result = (1.0 + as * _CV1(w, z));
        result += eps_c * (L2 - L5);
        result += eps_b * (L1 - L4);
        result += eps_c * eps_c * (_l2(w) - _l5(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::_h_t1(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
        const double L2 = -4.0 * chi3;

        double result = (1.0 + as * (_CT1(w, z) + (w - 1.0) / 2.0 * (_CT2(w, z) - _CT3(w, z))));
        result += eps_c * L2;
        result += eps_b * L1;
        result += eps_c * eps_c * _l2(w);

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::_h_t2(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L4 = 2.0 * eta - 1.0;
        const double L5 = -1.0;

        double result = (0.0 + as * (w + 1.0) / 2.0 * (_CT2(w, z) + _CT3(w, z)));
        result += eps_c * L5;
        result -= eps_b * L4;
        result += eps_c * eps_c * _l5(w);

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::_h_t3(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);

        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L3 = 4.0 * chi2;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (0.0 + as * _CT2(w, z));
        result += eps_c * (L6 - L3);
        result += eps_c * eps_c * (_l6(w) - _l3(w));

        return result * xi;
    }

    template <typename Process_>
    HQETFormFactors<Process_, PToV>::HQETFormFactors(const Parameters & p, const Options & o) :
        HQETFormFactorBase(p, o, Process_::hqe_prefix),
        _m_B(p[Process_::name_B], *static_cast<ParameterUser *>(this)),
        _m_V(p[Process_::name_V], *static_cast<ParameterUser *>(this))
    {
        static const Log::OneTimeMessage message_HQET_FFs_PToV
        (
            std::string("HQETFormFactors<") + Process_::label + ",PToV>",
            ll_warning,
            "This form factor parametrization is not a general one and requires careful attention. "
            "By default, it returns zeros for all form factors."
        );
    }

    template <typename Process_>
    HQETFormFactors<Process_, PToV>::~HQETFormFactors() = default;

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    HQETFormFactors<Process_, PToV>::begin_options()
    {
        return option_specifications.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    HQETFormFactors<Process_, PToV>::end_options()
    {
        return option_specifications.cend();
    }


    template <typename Process_>
    FormFactors<PToV> *
    HQETFormFactors<Process_, PToV>::make(const Parameters & parameters, const Options & options)
    {
        return new HQETFormFactors(parameters, options);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::v(const double & q2) const
    {
        const double r = _m_V / _m_B;

        // cf. [FKKM2008], eq. (22)
        return (1.0 + r) / 2.0 / sqrt(r) * _h_v(q2);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::a_0(const double & q2) const
    {
        const double r = _m_V / _m_B;
        const double w = _w(q2);

        return 1.0 / (2.0 * sqrt(r)) * ((1.0 + w) * _h_a1(q2) + (r * w - 1.0) * _h_a2(q2) + (r - w) * _h_a3(q2));
        // cf. [FKKM2008], eq. (22)
        //const double a_30 = (1.0 + r * r - 2.0 * r * w) / (4.0 * r * sqrt(r)) * (r * _h_a2(q2) - _h_a3(q2));
        //return a_3(q2) - a_30;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::a_1(const double & q2) const
    {
        const double r = _m_V / _m_B;
        const double w = _w(q2);

        // cf. [FKKM2008], eq. (22)
        return sqrt(r) * (1.0 + w) / (1.0 + r) * _h_a1(q2);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::a_2(const double & q2) const
    {
        const double r = _m_V / _m_B;

        // cf. [FKKM2008], eq. (22)
        return (1.0 + r) / (2.0 * sqrt(r)) * (r * _h_a2(q2) + _h_a3(q2));
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::a_3(const double & q2) const
    {
        const double r = _m_V / _m_B;

        // cf. [FKKM2008], below eq. (6)
        return ((1.0 + r) * a_1(q2) - (1.0 - r) * a_2(q2)) / (2.0 * r);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::a_12(const double & q2) const
    {
        const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
        const double m_V = this->_m_V(), m_V2 = power_of<2>(m_V);
        const double lambda = eos::lambda(m_B2, m_V2, q2);

        double result = (m_B + m_V) * (m_B + m_V) * (m_B2 - m_V2 - q2) * a_1(q2) - lambda * a_2(q2);
        result /= 16.0 * m_B * m_V2 * (m_B + m_V);

        return result;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::t_1(const double & q2) const
    {
        const double r = _m_V / _m_B;

        return -1.0 / (2.0 * sqrt(r)) * ((1.0 - r) * _h_t2(q2) - (1.0 + r) * _h_t1(q2));
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::t_2(const double & q2) const
    {
        const double r = _m_V / _m_B;
        const double w = _w(q2);

        return +1.0 / (2.0 * sqrt(r)) * (2.0 * r * (w + 1.0) / (1.0 + r) * _h_t1(q2) - 2.0 * r * (w - 1.0) / (1.0 - r) * _h_t2(q2));
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::t_3(const double & q2) const
    {
        const double r = _m_V / _m_B;

        return +1.0 / (2.0 * sqrt(r)) * ((1.0 - r) * _h_t1(q2) - (1.0 + r) * _h_t2(q2) + (1.0 - r * r) * _h_t3(q2));
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::t_23(const double & q2) const
    {
        const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
        const double m_V = this->_m_V(), m_V2 = power_of<2>(m_V);
        const double lambda = eos::lambda(m_B2, m_V2, q2);

        return ((m_B2 - m_V2) * (m_B2 + 3.0 * m_V2 - q2) * t_2(q2) - lambda * t_3(q2)) / (8.0 * m_B * m_V2 * (m_B - m_V));
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::f_perp(const double &) const
    {
        return 0.;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::f_para(const double &) const
    {
        return 0.;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::f_long(const double &) const
    {
        return 0.;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::f_perp_T(const double &) const
    {
        return 0.;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::f_para_T(const double &) const
    {
        return 0.;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, PToV>::f_long_T(const double &) const
    {
        return 0.;
    }

    template <typename Process_>
    Diagnostics
    HQETFormFactors<Process_, PToV>::diagnostics() const
    {
        Diagnostics results;

        // Inputs
        {
            const double m_b = _m_b_pole();
            const double m_c = _m_c_pole();
            const double z   = m_c / m_b;
            const double wz  = _wz(z);

            results.add(Diagnostics::Entry{ z,  "z = m_c_pole / m_b_pole" });
            results.add(Diagnostics::Entry{ wz, "w_z"                     });
        }

        // Switches
        {
            results.add(Diagnostics::Entry{ _enable_lp_z3,  "enable LP  z^3 terms" });
            results.add(Diagnostics::Entry{ _enable_lp_z4,  "enable LP  z^4 terms" });
            results.add(Diagnostics::Entry{ _enable_lp_z5,  "enable LP  z^5 terms" });
            results.add(Diagnostics::Entry{ _enable_slp_z2, "enable SLP z^2 terms" });
        }

        // z
        {
            results.add(Diagnostics::Entry{ _z(_q2(1.10)), "z(w = 1.10)" });
            results.add(Diagnostics::Entry{ _z(_q2(1.05)), "z(w = 1.05)" });
            results.add(Diagnostics::Entry{ _z(_q2(1.00)), "z(w = 1.00)" });
        }

        // xi
        {
            results.add(Diagnostics::Entry{ _xi(_q2(2.10)), "xi(w = 2.10)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.60)), "xi(w = 1.60)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.10)), "xi(w = 1.10)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.05)), "xi(w = 1.05)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.00)), "xi(w = 1.00)" });
        }

        // chi2
        {
            results.add(Diagnostics::Entry{ _chi2(_q2(2.10)), "chi2(w = 2.10)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.60)), "chi2(w = 1.60)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.10)), "chi2(w = 1.10)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.05)), "chi2(w = 1.05)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.00)), "chi2(w = 1.00)" });
        }

        // chi3
        {
            results.add(Diagnostics::Entry{ _chi3(_q2(2.10)), "chi3(w = 2.10)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.60)), "chi3(w = 1.60)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.10)), "chi3(w = 1.10)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.05)), "chi3(w = 1.05)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.00)), "chi3(w = 1.00)" });
        }

        // eta
        {
            results.add(Diagnostics::Entry{ _eta(_q2(2.10)), "eta(w = 2.10)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.60)), "eta(w = 1.60)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.10)), "eta(w = 1.10)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.05)), "eta(w = 1.05)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.00)), "eta(w = 1.00)" });
        }

        // r(w)
        {
            results.add(Diagnostics::Entry{ _r(1.1),     "r(w = 1.1)"     });
            results.add(Diagnostics::Entry{ _r(1.0007),  "r(w = 1.0007)"  });
            results.add(Diagnostics::Entry{ _r(1.0001),  "r(w = 1.0001)"  });
            results.add(Diagnostics::Entry{ _r(1.00005), "r(w = 1.00005)" });
            results.add(Diagnostics::Entry{ _r(1.0),     "r(w = 1.0)"     });
        }

        // Omega(w, z = 0.25)
        {
            results.add(Diagnostics::Entry{ _Omega(1.1,     0.25), "Omega(w = 1.1,     z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0007,  0.25), "Omega(w = 1.0007,  z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0001,  0.25), "Omega(w = 1.0001,  z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.00005, 0.25), "Omega(w = 1.00005, z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0,     0.25), "Omega(w = 1.0,     z = 0.25)" });
        }

        // Omega(w, z = 0.20)
        {
            results.add(Diagnostics::Entry{ _Omega(1.1,     0.20), "Omega(w = 1.1,     z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0007,  0.20), "Omega(w = 1.0007,  z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0001,  0.20), "Omega(w = 1.0001,  z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.00005, 0.20), "Omega(w = 1.00005, z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0,     0.20), "Omega(w = 1.0,     z = 0.20)" });
        }

        // WCs at w = 1.2, z = 0.20
        {
            results.add(Diagnostics::Entry{ _CS( 1.2, 0.20), "C_{S  }(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CP( 1.2, 0.20), "C_{P  }(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV1(1.2, 0.20), "C_{V_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV2(1.2, 0.20), "C_{V_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV3(1.2, 0.20), "C_{V_3}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA1(1.2, 0.20), "C_{A_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA2(1.2, 0.20), "C_{A_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA3(1.2, 0.20), "C_{A_3}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT1(1.2, 0.20), "C_{T_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT2(1.2, 0.20), "C_{T_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT3(1.2, 0.20), "C_{T_3}(w = 1.2, z = 0.20)" });
        }

        // WCs at w = 1.0, z = 0.25
        {
            results.add(Diagnostics::Entry{ _CS( 1.0, 0.25), "C_{S  }(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CP( 1.0, 0.25), "C_{P  }(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV1(1.0, 0.25), "C_{V_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV2(1.0, 0.25), "C_{V_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV3(1.0, 0.25), "C_{V_3}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA1(1.0, 0.25), "C_{A_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA2(1.0, 0.25), "C_{A_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA3(1.0, 0.25), "C_{A_3}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT1(1.0, 0.25), "C_{T_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT2(1.0, 0.25), "C_{T_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT3(1.0, 0.25), "C_{T_3}(w = 1.0, z = 0.25)" });
        }

        // HQET definition of the form factors
        {
            results.add(Diagnostics::Entry{ _h_a1(_q2(1.4)), "h_A1(w = 1.4)" });
            results.add(Diagnostics::Entry{ _h_a2(_q2(1.4)), "h_A2(w = 1.4)" });
            results.add(Diagnostics::Entry{ _h_a3(_q2(1.4)), "h_A3(w = 1.4)" });
            results.add(Diagnostics::Entry{ _h_v (_q2(1.4)), "h_V (w = 1.4)" });
            results.add(Diagnostics::Entry{ _h_t1(_q2(1.4)), "h_T1(w = 1.4)" });
            results.add(Diagnostics::Entry{ _h_t2(_q2(1.4)), "h_T2(w = 1.4)" });
            results.add(Diagnostics::Entry{ _h_t3(_q2(1.4)), "h_T3(w = 1.4)" });

            results.add(Diagnostics::Entry{ _h_a1(_q2(1.2)), "h_A1(w = 1.2)" });
            results.add(Diagnostics::Entry{ _h_a2(_q2(1.2)), "h_A2(w = 1.2)" });
            results.add(Diagnostics::Entry{ _h_a3(_q2(1.2)), "h_A3(w = 1.2)" });
            results.add(Diagnostics::Entry{ _h_v (_q2(1.2)), "h_V (w = 1.2)" });
            results.add(Diagnostics::Entry{ _h_t1(_q2(1.2)), "h_T1(w = 1.2)" });
            results.add(Diagnostics::Entry{ _h_t2(_q2(1.2)), "h_T2(w = 1.2)" });
            results.add(Diagnostics::Entry{ _h_t3(_q2(1.2)), "h_T3(w = 1.2)" });

            results.add(Diagnostics::Entry{ _h_a1(_q2(1.0)), "h_A1(w = 1.0)" });
            results.add(Diagnostics::Entry{ _h_a2(_q2(1.0)), "h_A2(w = 1.0)" });
            results.add(Diagnostics::Entry{ _h_a3(_q2(1.0)), "h_A3(w = 1.0)" });
            results.add(Diagnostics::Entry{ _h_v (_q2(1.0)), "h_V (w = 1.0)" });
            results.add(Diagnostics::Entry{ _h_t1(_q2(1.0)), "h_T1(w = 1.0)" });
            results.add(Diagnostics::Entry{ _h_t2(_q2(1.0)), "h_T2(w = 1.0)" });
            results.add(Diagnostics::Entry{ _h_t3(_q2(1.0)), "h_T3(w = 1.0)" });
        }

        return results;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToP>::_w(const double & q2) const
    {
        const double m_Bst = this->_m_Bst(), m_Bst2 = power_of<2>(m_Bst);
        const double m_P   = this->_m_P(),   m_P2   = power_of<2>(m_P);

        return (m_Bst2 + m_P2 - q2) / (2.0 * m_Bst * m_P);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToP>::_q2(const double & w) const
    {
        const double m_Bst = this->_m_Bst(), m_Bst2 = power_of<2>(m_Bst);
        const double m_P   = this->_m_P(),   m_P2   = power_of<2>(m_P);

        return m_Bst2 + m_P2 - 2.0 * m_Bst * m_P * w;
    }

    template <typename Process_>
    HQETFormFactors<Process_, VToP>::HQETFormFactors(const Parameters & p, const Options & o) :
        HQETFormFactorBase(p, o, Process_::hqe_prefix),
        _m_Bst(p[Process_::name_Bst], *static_cast<ParameterUser *>(this)),
        _m_P(p[Process_::name_P], *static_cast<ParameterUser *>(this))
    {
        static const Log::OneTimeMessage message_HQET_FFs_VToP
        (
            std::string("HQETFormFactors<") + Process_::label + ",VToP>",
            ll_warning,
            "This form factor parametrization is not a general one and requires careful attention. "
            "By default, it returns zeros for all form factors."
        );
    }

    template <typename Process_>
    HQETFormFactors<Process_, VToP>::~HQETFormFactors() = default;

    template <typename Process_>
    FormFactors<VToP> *
    HQETFormFactors<Process_, VToP>::make(const Parameters & parameters, const Options & options)
    {
        return new HQETFormFactors<Process_, VToP>(parameters, options);
    }

    /* HQET form factors h_i */

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToP>::h_abar_1(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
        const double L2 = -4.0 * chi3;
        const double L4 = 2.0 * eta - 1.0;
        const double L5 = -1.0;

        double result = (1.0 + as * _CA1(w, z));
        result += eps_c * (L1 - L4 * (w - 1.0) / (w + 1.0));
        result += eps_b * (L2 - L5 * (w - 1.0) / (w + 1.0));
        result += eps_c * eps_c * (_l1(w) - _l4(w) * (w - 1.0) / (w + 1.0));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToP>::h_abar_2(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L3 = 4.0 * chi2;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (0.0 - as * _CA3(w, z));
        result += eps_b * (L3 + L6);

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToP>::h_abar_3(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
        const double L2 = -4.0 * chi3;
        const double L3 = 4.0 * chi2;
        const double L4 = 2.0 * eta - 1.0;
        const double L5 = -1.0;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (1.0 + as * (_CA1(w, z) - _CA2(w, z)));
        result += eps_b * (L2 - L3 + L6 - L5);
        result += eps_c * (L1 - L4);
        result += eps_c * eps_c * (_l1(w) - _l4(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToP>::h_vbar(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
        const double L2 = -4.0 * chi3;
        const double L4 = 2.0 * eta - 1.0;
        const double L5 = -1.0;

        double result = (1.0 + as * _CV1(w, z));
        result += eps_b * (L2 - L5);
        result += eps_c * (L1 - L4);
        result += eps_c * eps_c * (_l1(w) - _l4(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToP>::h_tbar_1(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
        const double L2 = -4.0 * chi3;

        double result = (1.0 + as * (_CT1(w, z) - (w - 1.0) / 2.0 * (_CT2(w, z) - _CT3(w, z))));
        result += eps_b * (L2);
        result += eps_c * (L1);
        result += eps_c * eps_c * (_l1(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToP>::h_tbar_2(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L4 = 2.0 * eta - 1.0;
        const double L5 = -1.0;

        double result = (0.0 - as * (w + 1.0) / 2.0 * (_CT2(w, z) + _CT3(w, z)));
        result += eps_b * L5;
        result += eps_c * (-L4);
        result += eps_c * eps_c * (-_l4(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToP>::h_tbar_3(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L3 = 4.0 * chi2;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (0.0 - as * _CT3(w, z));
        result += eps_b * (L6 - L3);

        return result * xi;
    }

    template <typename Process_>
    Diagnostics
    HQETFormFactors<Process_, VToP>::diagnostics() const
    {
        Diagnostics results;

        // Inputs
        {
            const double m_b = _m_b_pole();
            const double m_c = _m_c_pole();
            const double z   = m_c / m_b;
            const double wz  = _wz(z);

            results.add(Diagnostics::Entry{ z,  "z = m_c / m_b" });
            results.add(Diagnostics::Entry{ wz, "w_z"           });
        }

        // Switches
        {
            results.add(Diagnostics::Entry{ _enable_lp_z3,  "enable LP  z^3 terms" });
            results.add(Diagnostics::Entry{ _enable_lp_z4,  "enable LP  z^4 terms" });
            results.add(Diagnostics::Entry{ _enable_lp_z5,  "enable LP  z^5 terms" });
            results.add(Diagnostics::Entry{ _enable_slp_z2, "enable SLP z^2 terms" });
        }

        // z
        {
            results.add(Diagnostics::Entry{ _z(_q2(1.10)), "z(w = 1.10)" });
            results.add(Diagnostics::Entry{ _z(_q2(1.05)), "z(w = 1.05)" });
            results.add(Diagnostics::Entry{ _z(_q2(1.00)), "z(w = 1.00)" });
        }

        // xi
        {
            results.add(Diagnostics::Entry{ _xi(_q2(2.10)), "xi(w = 2.10)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.60)), "xi(w = 1.60)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.10)), "xi(w = 1.10)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.05)), "xi(w = 1.05)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.00)), "xi(w = 1.00)" });
        }

        // chi2
        {
            results.add(Diagnostics::Entry{ _chi2(_q2(2.10)), "chi2(w = 2.10)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.60)), "chi2(w = 1.60)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.10)), "chi2(w = 1.10)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.05)), "chi2(w = 1.05)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.00)), "chi2(w = 1.00)" });
        }

        // chi3
        {
            results.add(Diagnostics::Entry{ _chi3(_q2(2.10)), "chi3(w = 2.10)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.60)), "chi3(w = 1.60)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.10)), "chi3(w = 1.10)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.05)), "chi3(w = 1.05)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.00)), "chi3(w = 1.00)" });
        }

        // eta
        {
            results.add(Diagnostics::Entry{ _eta(_q2(2.10)), "eta(w = 2.10)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.60)), "eta(w = 1.60)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.10)), "eta(w = 1.10)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.05)), "eta(w = 1.05)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.00)), "eta(w = 1.00)" });
        }

        // r(w)
        {
            results.add(Diagnostics::Entry{ _r(1.1),     "r(w = 1.1)"     });
            results.add(Diagnostics::Entry{ _r(1.0007),  "r(w = 1.0007)"  });
            results.add(Diagnostics::Entry{ _r(1.0001),  "r(w = 1.0001)"  });
            results.add(Diagnostics::Entry{ _r(1.00005), "r(w = 1.00005)" });
            results.add(Diagnostics::Entry{ _r(1.0),     "r(w = 1.0)"     });
        }

        // Omega(w, z = 0.25)
        {
            results.add(Diagnostics::Entry{ _Omega(1.1,     0.25), "Omega(w = 1.1,     z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0007,  0.25), "Omega(w = 1.0007,  z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0001,  0.25), "Omega(w = 1.0001,  z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.00005, 0.25), "Omega(w = 1.00005, z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0,     0.25), "Omega(w = 1.0,     z = 0.25)" });
        }

        // Omega(w, z = 0.20)
        {
            results.add(Diagnostics::Entry{ _Omega(1.1,     0.20), "Omega(w = 1.1,     z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0007,  0.20), "Omega(w = 1.0007,  z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0001,  0.20), "Omega(w = 1.0001,  z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.00005, 0.20), "Omega(w = 1.00005, z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0,     0.20), "Omega(w = 1.0,     z = 0.20)" });
        }

        // WCs at w = 1.2, z = 0.20
        {
            results.add(Diagnostics::Entry{ _CS( 1.2, 0.20), "C_{S  }(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CP( 1.2, 0.20), "C_{P  }(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV1(1.2, 0.20), "C_{V_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV2(1.2, 0.20), "C_{V_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV3(1.2, 0.20), "C_{V_3}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA1(1.2, 0.20), "C_{A_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA2(1.2, 0.20), "C_{A_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA3(1.2, 0.20), "C_{A_3}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT1(1.2, 0.20), "C_{T_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT2(1.2, 0.20), "C_{T_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT3(1.2, 0.20), "C_{T_3}(w = 1.2, z = 0.20)" });
        }

        // WCs at w = 1.0, z = 0.25
        {
            results.add(Diagnostics::Entry{ _CS( 1.0, 0.25), "C_{S  }(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CP( 1.0, 0.25), "C_{P  }(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV1(1.0, 0.25), "C_{V_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV2(1.0, 0.25), "C_{V_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV3(1.0, 0.25), "C_{V_3}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA1(1.0, 0.25), "C_{A_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA2(1.0, 0.25), "C_{A_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA3(1.0, 0.25), "C_{A_3}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT1(1.0, 0.25), "C_{T_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT2(1.0, 0.25), "C_{T_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT3(1.0, 0.25), "C_{T_3}(w = 1.0, z = 0.25)" });
        }

        // HQET definition of the form factors
        {
            results.add(Diagnostics::Entry{ h_abar_1(_q2(1.4)), "h_Abar1(w = 1.4)" });
            results.add(Diagnostics::Entry{ h_abar_2(_q2(1.4)), "h_Abar2(w = 1.4)" });
            results.add(Diagnostics::Entry{ h_abar_3(_q2(1.4)), "h_Abar3(w = 1.4)" });
            results.add(Diagnostics::Entry{ h_vbar (_q2(1.4)),  "h_Vbar (w = 1.4)" });
            results.add(Diagnostics::Entry{ h_tbar_1(_q2(1.4)), "h_Tbar1(w = 1.4)" });
            results.add(Diagnostics::Entry{ h_tbar_2(_q2(1.4)), "h_Tbar2(w = 1.4)" });
            results.add(Diagnostics::Entry{ h_tbar_3(_q2(1.4)), "h_Tbar3(w = 1.4)" });

            results.add(Diagnostics::Entry{ h_abar_1(_q2(1.2)), "h_Abar1(w = 1.2)" });
            results.add(Diagnostics::Entry{ h_abar_2(_q2(1.2)), "h_Abar2(w = 1.2)" });
            results.add(Diagnostics::Entry{ h_abar_3(_q2(1.2)), "h_Abar3(w = 1.2)" });
            results.add(Diagnostics::Entry{ h_vbar (_q2(1.2)),  "h_Vbar (w = 1.2)" });
            results.add(Diagnostics::Entry{ h_tbar_1(_q2(1.2)), "h_Tbar1(w = 1.2)" });
            results.add(Diagnostics::Entry{ h_tbar_2(_q2(1.2)), "h_Tbar2(w = 1.2)" });
            results.add(Diagnostics::Entry{ h_tbar_3(_q2(1.2)), "h_Tbar3(w = 1.2)" });

            results.add(Diagnostics::Entry{ h_abar_1(_q2(1.0)), "h_Abar1(w = 1.0)" });
            results.add(Diagnostics::Entry{ h_abar_2(_q2(1.0)), "h_Abar2(w = 1.0)" });
            results.add(Diagnostics::Entry{ h_abar_3(_q2(1.0)), "h_Abar3(w = 1.0)" });
            results.add(Diagnostics::Entry{ h_vbar (_q2(1.0)),  "h_Vbar (w = 1.0)" });
            results.add(Diagnostics::Entry{ h_tbar_1(_q2(1.0)), "h_Tbar1(w = 1.0)" });
            results.add(Diagnostics::Entry{ h_tbar_2(_q2(1.0)), "h_Tbar2(w = 1.0)" });
            results.add(Diagnostics::Entry{ h_tbar_3(_q2(1.0)), "h_Tbar3(w = 1.0)" });
        }

        return results;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::_w(const double & q2) const
    {
        static constexpr double mV1 = Process_::m_V1, mV12 = power_of<2>(Process_::m_V1);
        static constexpr double mV2 = Process_::m_V2, mV22 = power_of<2>(Process_::m_V2);

        return (mV12 + mV22 - q2) / (2.0 * mV1 * mV2);
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::_q2(const double & w) const
    {
        static constexpr double mV1 = Process_::m_V1, mV12 = power_of<2>(Process_::m_V1);
        static constexpr double mV2 = Process_::m_V2, mV22 = power_of<2>(Process_::m_V2);

        return mV12 + mV22 - 2.0 * mV1 * mV2 * w;
    }

    template <typename Process_>
    HQETFormFactors<Process_, VToV>::HQETFormFactors(const Parameters & p, const Options & o) :
        HQETFormFactorBase(p, o, Process_::hqe_prefix)
    {
        static const Log::OneTimeMessage message_HQET_FFs_VToV
        (
            std::string("HQETFormFactors<") + Process_::label + ",VToV>",
            ll_warning,
            "This form factor parametrization is not a general one and requires careful attention. "
            "By default, it returns zeros for all form factors."
        );
    }

    template <typename Process_>
    HQETFormFactors<Process_, VToV>::~HQETFormFactors() = default;

    template <typename Process_>
    FormFactors<VToV> *
    HQETFormFactors<Process_, VToV>::make(const Parameters & parameters, const Options & options)
    {
        return new HQETFormFactors(parameters, options);
    }

    /* HQET form factors h_i */

    // vector current
    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::h_1(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L2 = -4.0 * chi3;

        double result = 1.0 + as * (_CV1(w, z) + (w + 1.0) / 2.0 * (_CV2(w, z) + _CV3(w, z)));
        result += eps_c * L2;
        result += eps_b * L2;
        result += eps_c * eps_c * _l2(w);

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::h_2(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L5 = -1.0;

        double result = as * (w + 1.0) / 2.0 * (_CV2(w, z) - _CV3(w, z));
        result += eps_c * L5;
        result -= eps_b * L5;
        result += eps_c * eps_c * _l5(w);

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::h_3(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L2 = -4.0 * chi3;
        const double L3 = 4.0 * chi2;
        const double L5 = -1.0;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (1.0 + as * _CV1(w, z));
        result += eps_c * (L2 + L5 + (w - 1.0) * L3 - (w + 1.0) * L6);
        result += eps_b * (L2 - L5);
        result += eps_c * eps_c * (_l2(w) + _l5(w) + (w - 1.0) * _l3(w) - (w + 1.0) * _l6(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::h_4(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L2 = -4.0 * chi3;
        const double L3 = 4.0 * chi2;
        const double L5 = -1.0;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (1.0 + as * _CV1(w, z));
        result += eps_b * (L2 + L5 + (w - 1.0) * L3 - (w + 1.0) * L6);
        result += eps_c * (L2 - L5);
        result += eps_c * eps_c * (_l2(w) - _l5(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::h_5(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);

        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L3 = 4.0 * chi2;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (0.0 - as * _CV2(w, z));
        result += eps_c * (L3 - L6);
        result += eps_c * eps_c * (_l3(w) - _l6(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::h_6(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L3 = 4.0 * chi2;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (0.0 - as * _CV3(w, z));
        result += eps_b * (L3 - L6);

        return result * xi;
    }

    // axial current
    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::h_7(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double chi3 = _chi3(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L2 = -4.0 * chi3;

        double result = 1.0 + as * (_CA1(w, z) + (w - 1.0) / 2.0 * (_CA2(w, z) - _CA3(w, z)));
        result += eps_b * L2;
        result += eps_c * L2;
        result += eps_c * eps_c * _l2(w);

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::h_8(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L5 = -1.0;

        double result = as * (w + 1.0) / 2.0 * (_CA2(w, z) + _CA3(w, z));
        result += eps_c * L5;
        result -= eps_b * L5;
        result += eps_c * eps_c * _l5(w);

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::h_9(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);

        const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L3 = 4.0 * chi2;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (0.0 - as * _CA2(w, z));
        result += eps_c * (L3 - L6);
        result += eps_c * eps_c * (_l3(w) - _l6(w));

        return result * xi;
    }

    template <typename Process_>
    double
    HQETFormFactors<Process_, VToV>::h_10(const double & q2) const
    {
        const double m_b_pole = _m_b_pole();
        const double m_c_pole = _m_c_pole();

        const double w = this->_w(q2);
        const double z = m_c_pole / m_b_pole;

        const double as = _alpha_s() / M_PI;

        const double xi  = _xi(q2);
        const double eta = _eta(q2);
        const double chi2 = _chi2(q2);

        const double eps_b = _LambdaBar() / (2.0 * m_b_pole);

        // chi_1 is absorbed into def. of xi for LP and LV
        const double L3 = 4.0 * chi2;
        const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

        double result = (0.0 + as * _CA3(w, z));
        result += eps_b * (L3 - L6);

        return result * xi;
    }

    template <typename Process_>
    Diagnostics
    HQETFormFactors<Process_, VToV>::diagnostics() const
    {
        Diagnostics results;

        // Inputs
        {
            const double m_b = _m_b_pole();
            const double m_c = _m_c_pole();
            const double z   = m_c / m_b;
            const double wz  = _wz(z);

            results.add(Diagnostics::Entry{ z,  "z = m_c / m_b" });
            results.add(Diagnostics::Entry{ wz, "w_z"           });
        }

        // Switches
        {
            results.add(Diagnostics::Entry{ _enable_lp_z3,  "enable LP  z^3 terms" });
            results.add(Diagnostics::Entry{ _enable_lp_z4,  "enable LP  z^4 terms" });
            results.add(Diagnostics::Entry{ _enable_lp_z5,  "enable LP  z^5 terms" });
            results.add(Diagnostics::Entry{ _enable_slp_z2, "enable SLP z^2 terms" });
        }

        // z
        {
            results.add(Diagnostics::Entry{ _z(_q2(1.10)), "z(w = 1.10)" });
            results.add(Diagnostics::Entry{ _z(_q2(1.05)), "z(w = 1.05)" });
            results.add(Diagnostics::Entry{ _z(_q2(1.00)), "z(w = 1.00)" });
        }

        // xi
        {
            results.add(Diagnostics::Entry{ _xi(_q2(2.10)), "xi(w = 2.10)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.60)), "xi(w = 1.60)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.10)), "xi(w = 1.10)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.05)), "xi(w = 1.05)" });
            results.add(Diagnostics::Entry{ _xi(_q2(1.00)), "xi(w = 1.00)" });
        }

        // chi2
        {
            results.add(Diagnostics::Entry{ _chi2(_q2(2.10)), "chi2(w = 2.10)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.60)), "chi2(w = 1.60)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.10)), "chi2(w = 1.10)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.05)), "chi2(w = 1.05)" });
            results.add(Diagnostics::Entry{ _chi2(_q2(1.00)), "chi2(w = 1.00)" });
        }

        // chi3
        {
            results.add(Diagnostics::Entry{ _chi3(_q2(2.10)), "chi3(w = 2.10)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.60)), "chi3(w = 1.60)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.10)), "chi3(w = 1.10)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.05)), "chi3(w = 1.05)" });
            results.add(Diagnostics::Entry{ _chi3(_q2(1.00)), "chi3(w = 1.00)" });
        }

        // eta
        {
            results.add(Diagnostics::Entry{ _eta(_q2(2.10)), "eta(w = 2.10)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.60)), "eta(w = 1.60)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.10)), "eta(w = 1.10)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.05)), "eta(w = 1.05)" });
            results.add(Diagnostics::Entry{ _eta(_q2(1.00)), "eta(w = 1.00)" });
        }

        // r(w)
        {
            results.add(Diagnostics::Entry{ _r(1.1),     "r(w = 1.1)"     });
            results.add(Diagnostics::Entry{ _r(1.0007),  "r(w = 1.0007)"  });
            results.add(Diagnostics::Entry{ _r(1.0001),  "r(w = 1.0001)"  });
            results.add(Diagnostics::Entry{ _r(1.00005), "r(w = 1.00005)" });
            results.add(Diagnostics::Entry{ _r(1.0),     "r(w = 1.0)"     });
        }

        // Omega(w, z = 0.25)
        {
            results.add(Diagnostics::Entry{ _Omega(1.1,     0.25), "Omega(w = 1.1,     z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0007,  0.25), "Omega(w = 1.0007,  z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0001,  0.25), "Omega(w = 1.0001,  z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.00005, 0.25), "Omega(w = 1.00005, z = 0.25)" });
            results.add(Diagnostics::Entry{ _Omega(1.0,     0.25), "Omega(w = 1.0,     z = 0.25)" });
        }

        // Omega(w, z = 0.20)
        {
            results.add(Diagnostics::Entry{ _Omega(1.1,     0.20), "Omega(w = 1.1,     z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0007,  0.20), "Omega(w = 1.0007,  z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0001,  0.20), "Omega(w = 1.0001,  z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.00005, 0.20), "Omega(w = 1.00005, z = 0.20)" });
            results.add(Diagnostics::Entry{ _Omega(1.0,     0.20), "Omega(w = 1.0,     z = 0.20)" });
        }

        // WCs at w = 1.2, z = 0.20
        {
            results.add(Diagnostics::Entry{ _CS( 1.2, 0.20), "C_{S  }(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CP( 1.2, 0.20), "C_{P  }(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV1(1.2, 0.20), "C_{V_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV2(1.2, 0.20), "C_{V_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CV3(1.2, 0.20), "C_{V_3}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA1(1.2, 0.20), "C_{A_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA2(1.2, 0.20), "C_{A_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CA3(1.2, 0.20), "C_{A_3}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT1(1.2, 0.20), "C_{T_1}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT2(1.2, 0.20), "C_{T_2}(w = 1.2, z = 0.20)" });
            results.add(Diagnostics::Entry{ _CT3(1.2, 0.20), "C_{T_3}(w = 1.2, z = 0.20)" });
        }

        // WCs at w = 1.0, z = 0.25
        {
            results.add(Diagnostics::Entry{ _CS( 1.0, 0.25), "C_{S  }(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CP( 1.0, 0.25), "C_{P  }(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV1(1.0, 0.25), "C_{V_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV2(1.0, 0.25), "C_{V_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CV3(1.0, 0.25), "C_{V_3}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA1(1.0, 0.25), "C_{A_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA2(1.0, 0.25), "C_{A_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CA3(1.0, 0.25), "C_{A_3}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT1(1.0, 0.25), "C_{T_1}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT2(1.0, 0.25), "C_{T_2}(w = 1.0, z = 0.25)" });
            results.add(Diagnostics::Entry{ _CT3(1.0, 0.25), "C_{T_3}(w = 1.0, z = 0.25)" });
        }

        // HQET definition of the form factors
        {
            results.add(Diagnostics::Entry{ h_1(_q2(1.4)),  "h_1 (w = 1.4)" });
            results.add(Diagnostics::Entry{ h_2(_q2(1.4)),  "h_2 (w = 1.4)" });
            results.add(Diagnostics::Entry{ h_3(_q2(1.4)),  "h_3 (w = 1.4)" });
            results.add(Diagnostics::Entry{ h_4(_q2(1.4)),  "h_4 (w = 1.4)" });
            results.add(Diagnostics::Entry{ h_5(_q2(1.4)),  "h_5 (w = 1.4)" });
            results.add(Diagnostics::Entry{ h_6(_q2(1.4)),  "h_6 (w = 1.4)" });
            results.add(Diagnostics::Entry{ h_7(_q2(1.4)),  "h_7 (w = 1.4)" });
            results.add(Diagnostics::Entry{ h_8(_q2(1.4)),  "h_8 (w = 1.4)" });
            results.add(Diagnostics::Entry{ h_9(_q2(1.4)),  "h_9 (w = 1.4)" });
            results.add(Diagnostics::Entry{ h_10(_q2(1.4)), "h_10(w = 1.4)" });

            results.add(Diagnostics::Entry{ h_1(_q2(1.2)),  "h_1 (w = 1.2)" });
            results.add(Diagnostics::Entry{ h_2(_q2(1.2)),  "h_2 (w = 1.2)" });
            results.add(Diagnostics::Entry{ h_3(_q2(1.2)),  "h_3 (w = 1.2)" });
            results.add(Diagnostics::Entry{ h_4(_q2(1.2)),  "h_4 (w = 1.2)" });
            results.add(Diagnostics::Entry{ h_5(_q2(1.2)),  "h_5 (w = 1.2)" });
            results.add(Diagnostics::Entry{ h_6(_q2(1.2)),  "h_6 (w = 1.2)" });
            results.add(Diagnostics::Entry{ h_7(_q2(1.2)),  "h_7 (w = 1.2)" });
            results.add(Diagnostics::Entry{ h_8(_q2(1.2)),  "h_8 (w = 1.2)" });
            results.add(Diagnostics::Entry{ h_9(_q2(1.2)),  "h_9 (w = 1.2)" });
            results.add(Diagnostics::Entry{ h_10(_q2(1.2)), "h_10(w = 1.2)" });

            results.add(Diagnostics::Entry{ h_1(_q2(1.0)),  "h_1 (w = 1.0)" });
            results.add(Diagnostics::Entry{ h_2(_q2(1.0)),  "h_2 (w = 1.0)" });
            results.add(Diagnostics::Entry{ h_3(_q2(1.0)),  "h_3 (w = 1.0)" });
            results.add(Diagnostics::Entry{ h_4(_q2(1.0)),  "h_4 (w = 1.0)" });
            results.add(Diagnostics::Entry{ h_5(_q2(1.0)),  "h_5 (w = 1.0)" });
            results.add(Diagnostics::Entry{ h_6(_q2(1.0)),  "h_6 (w = 1.0)" });
            results.add(Diagnostics::Entry{ h_7(_q2(1.0)),  "h_7 (w = 1.0)" });
            results.add(Diagnostics::Entry{ h_8(_q2(1.0)),  "h_8 (w = 1.0)" });
            results.add(Diagnostics::Entry{ h_9(_q2(1.0)),  "h_9 (w = 1.0)" });
            results.add(Diagnostics::Entry{ h_10(_q2(1.0)), "h_10(w = 1.0)" });
        }

        return results;
    }
}

#endif
