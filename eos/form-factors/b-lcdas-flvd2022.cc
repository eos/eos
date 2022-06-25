/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Philip Lüghausen
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

#include <eos/form-factors/b-lcdas-flvd2022.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>

namespace eos
{
    namespace b_lcdas
    {
        FLvD2022::FLvD2022(const Parameters & p, const Options & o) :
            opt_q(o, "q", { "u", "s" }, "u"),
            opt_gminus(o, "gminus", { "zero", "WW-limit" }, "WW-limit"),
            switch_gminus(1.0),
            mu_0(p[parameter("mu_0")], *this),
            omega_0(p[parameter("omega_0")], *this),
            a({
                UsedParameter(p[parameter("a^phi+_0")], *this),
                UsedParameter(p[parameter("a^phi+_1")], *this),
                UsedParameter(p[parameter("a^phi+_2")], *this),
                UsedParameter(p[parameter("a^phi+_3")], *this),
                UsedParameter(p[parameter("a^phi+_4")], *this),
                UsedParameter(p[parameter("a^phi+_5")], *this),
                UsedParameter(p[parameter("a^phi+_6")], *this),
                UsedParameter(p[parameter("a^phi+_7")], *this),
                UsedParameter(p[parameter("a^phi+_8")], *this)
            })
        {
            if (opt_gminus.value() == "zero")
            {
                switch_gminus = 0.0;
            }
        }

        std::string
        FLvD2022::parameter(const char * _name) const
        {
            qnp::Name name(_name);

            if (opt_q.value() == "s")
                return QualifiedName(qnp::Prefix("B_s"), name, qnp::Suffix("FLvD2022")).str();

            return QualifiedName(qnp::Prefix("B_u"), name, qnp::Suffix("FLvD2022")).str();
        }

        BMesonLCDAs *
        FLvD2022::make(const Parameters & parameters, const Options & options)
        {
            return new FLvD2022(parameters, options);
        }

        std::tuple<BMesonLCDAs::CoefficientIterator, BMesonLCDAs::CoefficientIterator>
        FLvD2022::coefficient_range(const double & mu) const
        {
            // (WIP) only implement mu == mu_0 for the time being
            if (std::abs(mu_0 - mu) > 1.0e-12)
            {
                throw InternalError("FLvD2022: The case mu != mu_0 is not yet implemented");
            }

            // copy values to array of doubles
            static thread_local std::array<double, 9> values;
            for (size_t i = 0; i < values.size(); i++)
            {
                values[i] = a[i]; // evaluates UsedParameter
            }

            return {values.begin(), values.end()};
        }


        /* Leading twist two-particle LCDAs */

        double
        FLvD2022::phi_plus(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_minus(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }


        double
        FLvD2022::phi_bar(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar_d1(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }


        /* Next-to-leading twist two-particle LCDAs */

        double
        FLvD2022::g_minusWW(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_minusWW_d1(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_minusWW_d2(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_plus(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_plus_d1(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_plus_d2(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_bar(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_bar_d1(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_bar_d2(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::g_bar_d3(const double & omega) const
        {
            throw InternalError("Function not yet implemented");
        }

        /* Leading twist three-particle LCDAs */

        double
        FLvD2022::phi_3(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_4(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar_3(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar_4(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar2_3(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar2_4(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar_bar_3(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::phi_bar_bar_4(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::psi_bar_4(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::psi_bar_bar_4(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }


        double
        FLvD2022::chi_bar_4(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::chi_bar_bar_4(const double & omega_1, const double & omega_2) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::inverse_lambda_plus() const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::psi_A(const double & omega, const double & xi) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::psi_V(const double & omega, const double & xi) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::X_A(const double & omega, const double & xi) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::Y_A(const double & omega, const double & xi) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::Xbar_A(const double & omega, const double & xi) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022::Ybar_A(const double & omega, const double & xi) const
        {
            throw InternalError("Function not yet implemented");
        }

        Diagnostics
        FLvD2022::diagnostics() const
        {
            Diagnostics results;
            // add diagnostic results here
            return results;
        }
    }

    template <>
    struct WrappedForwardIteratorTraits<BMesonLCDAs::CoefficientIteratorTag>
    {
        using UnderlyingIterator = std::array<double, 9>::const_iterator;
    };
    template class WrappedForwardIterator<BMesonLCDAs::CoefficientIteratorTag, const double &>;
}
