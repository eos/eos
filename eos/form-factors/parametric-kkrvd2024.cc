/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020-2024 Danny van Dyk
 * Copyright (c) 2024 Matthew J. Kirk
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

#include <eos/form-factors/parametric-kkrvd2024.hh>
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
    KKRvD2024FormFactors<VacuumToPiPi>::KKRvD2024FormFactors(const Parameters & p, const Options & /*o*/) :
        _b_fp_I1{{
            UsedParameter(p[_par_name("+", "1", "2")], *this),
            UsedParameter(p[_par_name("+", "1", "3")], *this),
            UsedParameter(p[_par_name("+", "1", "4")], *this),
            UsedParameter(p[_par_name("+", "1", "5")], *this),
            UsedParameter(p[_par_name("+", "1", "6")], *this),
            UsedParameter(p[_par_name("+", "1", "7")], *this),
            UsedParameter(p[_par_name("+", "1", "8")], *this),
            UsedParameter(p[_par_name("+", "1", "9")], *this)
        }},
        _M_fp_I1(p["0->pipi::M_(+,1)@KKRvD2024"], *this),
        _G_fp_I1(p["0->pipi::Gamma_(+,1)@KKRvD2024"], *this),
        _m_pi(p["mass::pi^+"], *this),
        _t_0(p["0->pipi::t_0@KKRvD2024"], *this),
        _hbar(p["QM::hbar"], *this)
    {
    }

    KKRvD2024FormFactors<VacuumToPiPi>::~KKRvD2024FormFactors() = default;

    FormFactors<VacuumToPP> *
    KKRvD2024FormFactors<VacuumToPiPi>::make(const Parameters & p, const Options & o)
    {
        return new KKRvD2024FormFactors<VacuumToPiPi>(p, o);
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::z(const complex<double> & q2) const
    {
        return _z(q2, this->_t_0());
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::dzdq2(const complex<double> & q2) const
    {
        const double t_p = this->_t_p();
        const double t_0 = this->_t_0();
        return -sqrt(t_p - t_0) / (sqrt(t_p - q2) * power_of<2>(sqrt(t_p - q2) + sqrt(t_p - t_0)));
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::dzdq2_II(const complex<double> & q2) const
    {
        // dz/dq2 on the second Riemann sheet
        const double t_p = this->_t_p();
        const double t_0 = this->_t_0();
        return sqrt(t_p - t_0) / (sqrt(t_p - q2) * power_of<2>(sqrt(t_p - q2) - sqrt(t_p - t_0)));
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::w(const complex<double> & z) const
    {
        return power_of<2>(1.0 + z) * pow((1.0 - z), 5.0/2.0);
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::phitilde_p(const complex<double> & z, const double & chi) const
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
    KKRvD2024FormFactors<VacuumToPiPi>::phitildeprime_p(const complex<double> & z, const double & chi) const
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
    KKRvD2024FormFactors<VacuumToPiPi>::series_m(const complex<double> & z, const std::array<double, 10u> & c) const
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

    double
    KKRvD2024FormFactors<VacuumToPiPi>::_b0_fp_I1(const double & chi, const complex<double> & zr) const
    {
        const complex<double> z0 = this->z(0.0);

        const complex<double> phitilde_z0      = this->phitilde_p(z0, chi);
        const complex<double> phitilde_m1      = this->phitilde_p(-1.0, chi);
        const complex<double> phitildeprime_m1 = this->phitildeprime_p(-1.0, chi);

        const complex<double> x_z0      = 1.0 / (phitilde_z0 * std::norm(z0 - zr));
        const complex<double> x_m1      = 1.0 / (phitilde_m1 * std::norm(1.0 + zr));
        const complex<double> xprime_m1 = (2.0 * (1 + std::real(zr)) * phitilde_m1 - std::norm(1.0 + zr) * phitildeprime_m1 ) / power_of<2>(std::norm(1.0 + zr) * phitilde_m1);

        std::array<double, 10> b;
        b[0] = 0.0;
        b[1] = 0.0;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);

        complex<double> sum = 0.0;
        for (auto i = 0u; i < b.size(); i++)
        {
            sum += b[i] * ( (pow(-1, i) * i + pow(z0, i-1)) * x_m1 - (pow(-1, i) + pow(z0, i-1)) * xprime_m1 );
        };

        const double b0 = std::real((x_m1 - xprime_m1 - x_z0 * z0 * sum) / (x_z0 * (x_m1 - (1.0 + z0) * xprime_m1)));
        return b0;
    }

    double
    KKRvD2024FormFactors<VacuumToPiPi>::_b1_fp_I1(const double & chi, const complex<double> & zr) const
    {
        const complex<double> z0 = this->z(0.0);

        const complex<double> phitilde_z0      = this->phitilde_p(z0, chi);
        const complex<double> phitilde_m1      = this->phitilde_p(-1.0, chi);
        const complex<double> phitildeprime_m1 = this->phitildeprime_p(-1.0, chi);

        const complex<double> x_z0      = 1.0 / (phitilde_z0 * std::norm(z0 - zr));
        const complex<double> x_m1      = 1.0 / (phitilde_m1 * std::norm(1.0 + zr));
        const complex<double> xprime_m1 = (2.0 * (1 + std::real(zr)) * phitilde_m1 - std::norm(1.0 + zr) * phitildeprime_m1 ) / power_of<2>(std::norm(1.0 + zr) * phitilde_m1);

        std::array<double, 10> b;
        b[0] = 0.0;
        b[1] = 0.0;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);

        complex<double> sum = 0.0;
        for (auto i = 0u; i < b.size(); i++)
        {
            sum += b[i] * (pow(-1, i) * i * x_m1 + (pow(-1, i+1) + pow(z0, i)) * xprime_m1);
        };

        const double b1 = std::real((sum * x_z0 -  xprime_m1) / ( x_z0 * (x_m1 - (1.0 + z0) * xprime_m1)));
        return b1;
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::f_p(const double & q2) const
    {
        static const double eps = 1.0e-12;
        return f_p(complex<double>(q2, eps));
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::f_p(const complex<double> & q2) const
    {
        const auto z        = this->z(q2);
        const auto chi      = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A
        const auto phitilde = this->phitilde_p(z, chi);

        // Super-threshold pole location
        const auto zr = this->_zr(this->_M_fp_I1(), this->_G_fp_I1());

        // prepare expansion coefficients
        std::array<double, 10> b;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);
        // Fix b[0] and b[1] to enforce F(q2=0) = 1 and F'(q2=t+) = 0
        b[0] = _b0_fp_I1(chi, zr);
        b[1] = _b1_fp_I1(chi, zr);

        const auto series = this->series_m(z, b);

        return series / (z - zr) / (z - std::conj(zr)) /  phitilde;
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::f_t(const double & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::f_0(const double & /*q2*/) const
    {
        return 0.0; // vanishes in our approximation
    }

    double
    KKRvD2024FormFactors<VacuumToPiPi>::b_0() const
    {
        const auto chi = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.Av
        const auto zr  = this->_zr(this->_M_fp_I1(), this->_G_fp_I1());

        return _b0_fp_I1(chi, zr);
    }

    double
    KKRvD2024FormFactors<VacuumToPiPi>::b_1() const
    {
        const auto chi = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A
        const auto zr  = this->_zr(this->_M_fp_I1(), this->_G_fp_I1());

        return _b1_fp_I1(chi, zr);
    }

    double
    KKRvD2024FormFactors<VacuumToPiPi>::dFdq2_q2eq0() const
    {
        const double z0 = std::real(this->z(0.0));
        const double chi = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A

        const double phitilde_z0      = std::real(this->phitilde_p(z0, chi));
        const double phitildeprime_z0 = std::real(this->phitildeprime_p(z0, chi));

        // Super-threshold pole location
        const auto zr =  this->_zr(this->_M_fp_I1(), this->_G_fp_I1());

        // prepare expansion coefficients
        std::array<double, 10> b;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);
        // Fix b[0] and b[1] to enforce F(q2=0) = 1 and F'(q2=t+) = 0
        b[0] = _b0_fp_I1(chi, zr);
        b[1] = _b1_fp_I1(chi, zr);

        const double sum1 = std::real(series_m(z0, b));
        double sum2 = 0.0;
        for (auto i = 0u; i < b.size(); i++)
        {
            sum2 += b[i] * i * std::pow(z0, i-1);
        };

        const double xprime_z0 = ( 2.0 * (std::real(zr) - z0) * phitilde_z0 - std::norm(z0 - zr) * phitildeprime_z0 ) / power_of<2>(std::norm(z0 - zr) * phitilde_z0);

        const double dFdz_z0 = sum1 * xprime_z0 + sum2 / (phitilde_z0 * std::norm(z0 - zr));
        return dFdz_z0 * std::real(this->dzdq2(0.0));
    }

    double
    KKRvD2024FormFactors<VacuumToPiPi>::r_pi_squared() const
    {
        return 6.0 * this->dFdq2_q2eq0() * power_of<2>(this->_hbarc());
    }

    complex<double> KKRvD2024FormFactors<VacuumToPiPi>::residue_rho() const
    {
        // Super-threshold pole location
        const auto zr       = this->_zr(this->_M_fp_I1(), this->_G_fp_I1());
        const auto chi      = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A
        const auto phitilde = this->phitilde_p(zr, chi);

        // prepare expansion coefficients
        std::array<double, 10> b;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);
        // Fix b[0] and b[1] to enforce F(q2=0) = 1 and F'(q2=t+) = 0
        b[0] = _b0_fp_I1(chi, zr);
        b[1] = _b1_fp_I1(chi, zr);

        const auto series = this->series_m(zr, b);

        return series / (zr - std::conj(zr)) /  phitilde;
    }

    double KKRvD2024FormFactors<VacuumToPiPi>::re_residue_rho() const
    {
        return std::real(this->residue_rho());
    }

    double KKRvD2024FormFactors<VacuumToPiPi>::im_residue_rho() const
    {
        return std::imag(this->residue_rho());
    }

    complex<double> KKRvD2024FormFactors<VacuumToPiPi>::residue_rho_q2() const
    {
        const auto s_rho = power_of<2>(complex<double>(this->_M_fp_I1(), -this->_G_fp_I1()/2));
        return this->residue_rho() / this->dzdq2_II(s_rho);
    }

    double KKRvD2024FormFactors<VacuumToPiPi>::re_residue_rho_q2() const
    {
        return std::real(this->residue_rho_q2());
    }

    double KKRvD2024FormFactors<VacuumToPiPi>::im_residue_rho_q2() const
    {
        return std::imag(this->residue_rho_q2());
    }

    double KKRvD2024FormFactors<VacuumToPiPi>::dispersive_integrand(const double & alpha) const
    {
        const complex<double> z = std::polar(1.0, alpha);
        const complex<double> w = this->w(z);

        const auto chi = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A

        // Super-threshold pole location
        const auto zr =  this->_zr(this->_M_fp_I1(), this->_G_fp_I1());

        // prepare expansion coefficients
        std::array<double, 10> b;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);
        // Fix b[0] and b[1] to enforce F(q2=0) = 1 and F'(q2=t+) = 0
        b[0] = _b0_fp_I1(chi, zr);
        b[1] = _b1_fp_I1(chi, zr);

        const auto series = this->series_m(z, b);

        return std::norm(w * series / (z - zr) / (z - std::conj(zr)));
    }

    double KKRvD2024FormFactors<VacuumToPiPi>::saturation() const
    {
        std::function<double (const double &)> f = [this](const double & alpha) -> double { return this->dispersive_integrand(alpha); };
        return integrate<GSL::QAGS>(f, -M_PI, M_PI) / (2.0 * M_PI);
    }

    const std::set<ReferenceName>
    KKRvD2024FormFactors<VacuumToPiPi>::references
    {
        "BL:1998A"_rn
    };

    const std::vector<OptionSpecification>
    KKRvD2024FormFactors<VacuumToPiPi>::option_specifications
    {
    };

    std::vector<OptionSpecification>::const_iterator
    KKRvD2024FormFactors<VacuumToPiPi>::begin_options()
    {
        return option_specifications.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    KKRvD2024FormFactors<VacuumToPiPi>::end_options()
    {
        return option_specifications.cend();
    }

    template class KKRvD2024FormFactors<VacuumToPiPi>;
}
