/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Danny van Dyk
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

#include <eos/form-factors/parametric-ksvd2025.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>
#include <eos/maths/integrate.hh>

#include <functional>
#include <numeric>

namespace eos
{
    /* Vacuum -> pi pi */
    KSvD2025FormFactors<VacuumToKPi>::KSvD2025FormFactors(const Parameters & p, const Options & /*o*/) :
        _b_fp{{
            UsedParameter(p[_coeff_name("+", "2")], *this),
            UsedParameter(p[_coeff_name("+", "3")], *this),
            UsedParameter(p[_coeff_name("+", "4")], *this),
            UsedParameter(p[_coeff_name("+", "5")], *this),
            UsedParameter(p[_coeff_name("+", "6")], *this),
            UsedParameter(p[_coeff_name("+", "7")], *this),
            UsedParameter(p[_coeff_name("+", "8")], *this),
            UsedParameter(p[_coeff_name("+", "9")], *this)
        }},
        _M_fp{{
            UsedParameter(p["0->Kpi::M_(+,0)@KSvD2025"], *this),
            UsedParameter(p["0->Kpi::M_(+,1)@KSvD2025"], *this)
        }},
        _G_fp{{
            UsedParameter(p["0->Kpi::M_(+,0)@KSvD2025"], *this),
            UsedParameter(p["0->Kpi::M_(+,1)@KSvD2025"], *this)
        }},
        _m_K(p["mass::K_d"], *this),
        _m_pi(p["mass::pi^-"], *this),
        _t_0(p["0->Kpi::t_0@KSvD2025"], *this),
        _hbar(p["QM::hbar"], *this)
    {
    }

    KSvD2025FormFactors<VacuumToKPi>::~KSvD2025FormFactors() = default;

    FormFactors<VacuumToPP> *
    KSvD2025FormFactors<VacuumToKPi>::make(const Parameters & p, const Options & o)
    {
        return new KSvD2025FormFactors<VacuumToKPi>(p, o);
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::z(const complex<double> & q2) const
    {
        return _z(q2, this->_t_0());
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::dzdq2(const complex<double> & q2) const
    {
        const double t_p = this->_t_p();
        const double t_0 = this->_t_0();
        return -sqrt(t_p - t_0) / (sqrt(t_p - q2) * power_of<2>(sqrt(t_p - q2) + sqrt(t_p - t_0)));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::dzdq2_II(const complex<double> & q2) const
    {
        // dz/dq2 on the second Riemann sheet
        const double t_p = this->_t_p();
        const double t_0 = this->_t_0();
        return sqrt(t_p - t_0) / (sqrt(t_p - q2) * power_of<2>(sqrt(t_p - q2) - sqrt(t_p - t_0)));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::series_m(const complex<double> & z, const std::array<double, 10u> & c) const
    {
        std::array<complex<double>, 10> zvalues;
        complex<double> current_zv = 1.0;
        for (complex<double> & zv: zvalues)
        {
            zv = current_zv;
            current_zv *= z;
        }

        return std::inner_product(c.cbegin(), c.cend(), zvalues.cbegin(), complex<double>(0.0, 0.0));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::w_p(const complex<double> & z) const
    {
        return power_of<2>(1.0 + z) * pow((1.0 - z), 5.0 / 2.0);
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::phitilde_p(const complex<double> & z, const double & chi) const
    {
        // The weight function ``(1.0 + z)^2 * (1.0 - z)^(+5/2)`` has been cancelled against the outer function
        // to remove unphysical singularities and correct the asymptotic behaviour
        const double t_p     = this->_t_p();
        const double t_0     = this->_t_0();
        const double tfactor = 1.0 - t_0 / t_p;
        const double Q2      = 1.0;

        // cf. [BL:1998A], eq. (5.2), p. 11
        return 1.0
            / sqrt(12.0 * M_PI * t_p * chi)
            * pow(tfactor, 5.0 / 4.0) * pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5)
            * pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0)
            / power_of<2>(1.0 - z);
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::phitildeprime_p(const complex<double> & z, const double & chi) const
    {
        // Derivative of the modified outer function phitilde with respect to z
        const double t_p      = this->_t_p();
        const double t_0      = this->_t_0();
        const double tfactor  = 1.0 - t_0 / t_p;
        const double Q2       = 1.0;
        const double Q2factor = 1.0 + Q2 / t_p;

        return +1.0 * pow(tfactor, 5.0 / 4.0)
                    * (-11 * sqrt(Q2factor) * power_of<2>(z - 1.0) - tfactor * (1.0 + z) * (11.0 * z - 3.0) + sqrt(tfactor) * (z - 1.0) * (-1.0 + 9.0 * sqrt(Q2factor) + 11.0 * (1 + sqrt(Q2factor)) * z ))
                    / (4 * sqrt(3 * M_PI * t_p * chi) * power_of<3>(z - 1.0) * power_of<4>(sqrt(Q2factor) * (z - 1.0) - sqrt(tfactor) * (1.0 + z)) * pow(1.0 - z + sqrt(tfactor) * (1.0 + z), 3.0/2.0));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::w_z(const complex<double> & z) const
    {
        // TODO
        return 1.0;
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::phitilde_z(const complex<double> & z, const double & chi) const
    {
        // TODO
        return 1.0;
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::phitildeprime_z(const complex<double> & z, const double & chi) const
    {
        // TODO
        return 1.0;
    }

    double
    KSvD2025FormFactors<VacuumToKPi>::_b0_fp(const double & chi_p) const
    {
        // determine the coefficient b_0 of f_+(q^2) by imposing that Im f_+(q^2) ~ sqrt(q^2 - t_+)^3
        double b0 = 0.0;

        return b0;
    }

    double
    KSvD2025FormFactors<VacuumToKPi>::_b0_fz(const double & chi_p, const double & chi_z) const
    {
        // determine the coefficient b_0 of f_0(q^2) by imposing that f_+(0) = f_0(0)
        double b0 = 0.0;

        return b0;
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_p(const double & q2) const
    {
        static const double eps = 1.0e-12;
        return f_p(complex<double>(q2, eps));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_p(const complex<double> & q2) const
    {
        // TODO
        return 1.0;
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_0(const double & /*q2*/) const
    {
        return 0.0; // vanishes in our approximation
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_0(const complex<double> & /*q2*/) const
    {
        return 0.0; // vanishes in our approximation
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_t(const double & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_t(const complex<double> & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    double KSvD2025FormFactors<VacuumToKPi>::dispersive_integrand(const double & alpha) const
    {
        // TODO
        return 0.0;
    }

    double KSvD2025FormFactors<VacuumToKPi>::saturation() const
    {
        std::function<double (const double &)> f = [this](const double & alpha) -> double { return this->dispersive_integrand(alpha); };
        return integrate<GSL::QAGS>(f, -M_PI, M_PI) / (2.0 * M_PI);
    }

    const std::set<ReferenceName>
    KSvD2025FormFactors<VacuumToKPi>::references
    {
    };

    const std::vector<OptionSpecification>
    KSvD2025FormFactors<VacuumToKPi>::option_specifications
    {
    };

    std::vector<OptionSpecification>::const_iterator
    KSvD2025FormFactors<VacuumToKPi>::begin_options()
    {
        return option_specifications.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    KSvD2025FormFactors<VacuumToKPi>::end_options()
    {
        return option_specifications.cend();
    }

    template class KSvD2025FormFactors<VacuumToKPi>;
}
