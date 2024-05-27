/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Danny van Dyk
 * Copyright (c) 2024 Stefan Meiser
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

#include <eos/b-decays/lifetime.hh>
#include <eos/maths/matrix.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <complex>

namespace eos
{
    using std::real;

    template <>
    struct Implementation<Lifetime>
    {
        SpecifiedOption opt_model;
        std::shared_ptr<Model> model;

        UsedParameter hbar;
        UsedParameter g_fermi;

        QuarkFlavorOption opt_q;
        UsedParameter m_B;
        UsedParameter f_B;

        UsedParameter mu_dbcu;
        UsedParameter mu_sbcu;

        double switch_pauli_interference_dbcu;
        double switch_pauli_interference_sbcu;
        double switch_weak_exchange_dbcu;
        double switch_weak_exchange_sbcu;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_model(o, options, "model"),
            model(Model::make(opt_model.value(), p, o)),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            opt_q(o, options, "q"),
            m_B(p["mass::B_" + opt_q.str()], u),
            f_B(p["decay-constant::B_" + opt_q.str()], u),
            mu_dbcu(p["dbcu::mu"], u),
            mu_sbcu(p["sbcu::mu"], u)
        {
            Context ctx("When constructing a B meson lifetime observable");

            switch (opt_q.value())
            {
                case QuarkFlavor::up:
                    switch_pauli_interference_dbcu = 1.0;
                    switch_pauli_interference_sbcu = 1.0;
                    switch_weak_exchange_dbcu      = 0.0;
                    switch_weak_exchange_sbcu      = 0.0;
                    break;
                case QuarkFlavor::down:
                    switch_pauli_interference_dbcu = 0.0;
                    switch_pauli_interference_sbcu = 0.0;
                    switch_weak_exchange_dbcu      = 1.0;
                    switch_weak_exchange_sbcu      = 0.0;
                    break;
                case QuarkFlavor::strange:
                    switch_pauli_interference_dbcu = 0.0;
                    switch_pauli_interference_sbcu = 0.0;
                    switch_weak_exchange_dbcu      = 0.0;
                    switch_weak_exchange_sbcu      = 1.0;
                    break;
                default:
                    throw InternalError("Invalid quark flavor: " + stringify(opt_q.value()));
            }
            u.uses(*model);
        }

        inline double convert_GeV_to_inv_ps() const { return 1.0e-12 / hbar(); }

        // cf. [LMPR:2022A], eqs. (2.22) to (2.27)
        // provided by A. Rusov as a Mathematica expression
        std::array<std::array<complex<double>, 20u>, 20u>
        A_pauli_interference(const double & /*mu*/, const double & sqrtrho) const
        {
            const double rho = sqrtrho * sqrtrho;

            // Matrix elements of operators [bbar Gamma q] [qbar Gamma b]
            // are only known from HQET sum rules. For the time being use,
            // constant values for these matrix elements, mostly following
            // [LMPR:2022A]. For the matrix elements vanishing in vacuum
            // insertion approximation, we use 10% as the naive bag factor.
            const double me  =  f_B * f_B * m_B * m_B;
            const double me1 =  1.0 * me;
            const double me2 =  1.0 * me;
            const double me3 =  0.1 * me;
            const double me4 =  0.1 * me;
            const double me5 = -1.0 * me;
            const double me6 = -1.0 * me;
            const double me7 =  0.1 * me;
            const double me8 =  0.1 * me;

            std::array<std::array<complex<double>, 20u>, 20u> result
            {{
                {me1 + 6.0 * me3, 3.0 * me1, -0.5 * ((me1 + 6.0 * me3) * sqrtrho), (-3.0 * me1 * sqrtrho) / 2.0, -0.25 * ((me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) * sqrtrho), (-3.0 * (me5 - 2.0 * me6) * sqrtrho) / 4.0, (-me5 + 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) / 4.0, (-3.0 * (me5 - 2.0 * me6)) / 4.0, 3.0 * (me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)), 9.0 * (me5 - 2.0 * me6), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {3.0 * me1, me1 + 6.0 * me3, (-3.0 * me1 * sqrtrho) / 2.0, -0.5 * ((me1 + 6.0 * me3) * sqrtrho), (-3.0 * (me5 - 2.0 * me6) * sqrtrho) / 4.0, -0.25 * ((me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) * sqrtrho), (-3.0 * (me5 - 2.0 * me6)) / 4.0, (-me5 + 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) / 4.0, 9.0 * (me5 - 2.0 * me6), 3.0 * (me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {-0.5 * ((me1 + 6.0 * me3) * sqrtrho), (-3.0 * me1 * sqrtrho) / 2.0, (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho) + 6.0 * (me3 + 2.0 * me3 * rho - 2.0 * me4 * (2.0 + rho))) / 6.0, me1 * (0.5 + rho) - me2 * (2.0 + rho), (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho + 3.0 * (me7 + 2.0 * me8 - me7 * rho + 4.0 * me8 * rho))) / 12.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho)) / 4.0, -0.5 * ((me6 + 6.0 * me8) * sqrtrho), (-3.0 * me6 * sqrtrho) / 2.0, -2.0 * (me5 - me6 + 6.0 * me7 - 6.0 * me8) * sqrtrho, -6.0 * (me5 - me6) * sqrtrho, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {(-3.0 * me1 * sqrtrho) / 2.0, -0.5 * ((me1 + 6.0 * me3) * sqrtrho), me1 * (0.5 + rho) - me2 * (2.0 + rho), (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho) + 6.0 * (me3 + 2.0 * me3 * rho - 2.0 * me4 * (2.0 + rho))) / 6.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho)) / 4.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho + 3.0 * (me7 + 2.0 * me8 - me7 * rho + 4.0 * me8 * rho))) / 12.0, (-3.0 * me6 * sqrtrho) / 2.0, -0.5 * ((me6 + 6.0 * me8) * sqrtrho), -6.0 * (me5 - me6) * sqrtrho, -2.0 * (me5 - me6 + 6.0 * me7 - 6.0 * me8) * sqrtrho, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {-0.25 * ((me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) * sqrtrho), (-3.0 * (me5 - 2.0 * me6) * sqrtrho) / 4.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho + 3.0 * (me7 + 2.0 * me8 - me7 * rho + 4.0 * me8 * rho))) / 12.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho)) / 4.0, (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho) + 6.0 * (me3 + 2.0 * me3 * rho - 2.0 * me4 * (2.0 + rho))) / 24.0, (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho)) / 8.0, ((me1 - 2.0 * (me2 - 3.0 * me3 + 6.0 * me4)) * sqrtrho) / 8.0, (3.0 * (me1 - 2.0 * me2) * sqrtrho) / 8.0, -0.5 * ((me1 + 2.0 * (me2 + 3.0 * me3 + 6.0 * me4)) * sqrtrho), (-3.0 * (me1 + 2.0 * me2) * sqrtrho) / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {(-3.0 * (me5 - 2.0 * me6) * sqrtrho) / 4.0, -0.25 * ((me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) * sqrtrho), (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho)) / 4.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho + 3.0 * (me7 + 2.0 * me8 - me7 * rho + 4.0 * me8 * rho))) / 12.0, (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho)) / 8.0, (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho) + 6.0 * (me3 + 2.0 * me3 * rho - 2.0 * me4 * (2.0 + rho))) / 24.0, (3.0 * (me1 - 2.0 * me2) * sqrtrho) / 8.0, ((me1 - 2.0 * (me2 - 3.0 * me3 + 6.0 * me4)) * sqrtrho) / 8.0, (-3.0 * (me1 + 2.0 * me2) * sqrtrho) / 2.0, -0.5 * ((me1 + 2.0 * (me2 + 3.0 * me3 + 6.0 * me4)) * sqrtrho), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {(-me5 + 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) / 4.0, (-3.0 * (me5 - 2.0 * me6)) / 4.0, -0.5 * ((me6 + 6.0 * me8) * sqrtrho), (-3.0 * me6 * sqrtrho) / 2.0, ((me1 - 2.0 * (me2 - 3.0 * me3 + 6.0 * me4)) * sqrtrho) / 8.0, (3.0 * (me1 - 2.0 * me2) * sqrtrho) / 8.0, (me1 * (2.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho - 3.0 * me3 * (2.0 + rho) + 6.0 * (me4 + 2.0 * me4 * rho))) / 24.0, (me1 * (2.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho)) / 8.0, (me1 * (-4.0 + rho) - 2.0 * (me2 - 3.0 * me3 * (-4.0 + rho) + 2.0 * me2 * rho + 6.0 * (me4 + 2.0 * me4 * rho))) / 6.0, (me1 * (-4.0 + rho)) / 2.0 - me2 * (1.0 + 2.0 * rho), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {(-3.0 * (me5 - 2.0 * me6)) / 4.0, (-me5 + 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) / 4.0, (-3.0 * me6 * sqrtrho) / 2.0, -0.5 * ((me6 + 6.0 * me8) * sqrtrho), (3.0 * (me1 - 2.0 * me2) * sqrtrho) / 8.0, ((me1 - 2.0 * (me2 - 3.0 * me3 + 6.0 * me4)) * sqrtrho) / 8.0, (me1 * (2.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho)) / 8.0, (me1 * (2.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho - 3.0 * me3 * (2.0 + rho) + 6.0 * (me4 + 2.0 * me4 * rho))) / 24.0, (me1 * (-4.0 + rho)) / 2.0 - me2 * (1.0 + 2.0 * rho), (me1 * (-4.0 + rho) - 2.0 * (me2 - 3.0 * me3 * (-4.0 + rho) + 2.0 * me2 * rho + 6.0 * (me4 + 2.0 * me4 * rho))) / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {3.0 * (me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)), 9.0 * (me5 - 2.0 * me6), -2.0 * (me5 - me6 + 6.0 * me7 - 6.0 * me8) * sqrtrho, -6.0 * (me5 - me6) * sqrtrho, -0.5 * ((me1 + 2.0 * (me2 + 3.0 * me3 + 6.0 * me4)) * sqrtrho), (-3.0 * (me1 + 2.0 * me2) * sqrtrho) / 2.0, (me1 * (-4.0 + rho) - 2.0 * (me2 - 3.0 * me3 * (-4.0 + rho) + 2.0 * me2 * rho + 6.0 * (me4 + 2.0 * me4 * rho))) / 6.0, (me1 * (-4.0 + rho)) / 2.0 - me2 * (1.0 + 2.0 * rho), (2.0 * (me1 * (14.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho - 3.0 * me3 * (14.0 + rho) + 6.0 * (me4 + 2.0 * me4 * rho)))) / 3.0, 2.0 * me1 * (14.0 + rho) - 4.0 * (me2 + 2.0 * me2 * rho), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {9.0 * (me5 - 2.0 * me6), 3.0 * (me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)), -6.0 * (me5 - me6) * sqrtrho, -2.0 * (me5 - me6 + 6.0 * me7 - 6.0 * me8) * sqrtrho, (-3.0 * (me1 + 2.0 * me2) * sqrtrho) / 2.0, -0.5 * ((me1 + 2.0 * (me2 + 3.0 * me3 + 6.0 * me4)) * sqrtrho), (me1 * (-4.0 + rho)) / 2.0 - me2 * (1.0 + 2.0 * rho), (me1 * (-4.0 + rho) - 2.0 * (me2 - 3.0 * me3 * (-4.0 + rho) + 2.0 * me2 * rho + 6.0 * (me4 + 2.0 * me4 * rho))) / 6.0, 2.0 * me1 * (14.0 + rho) - 4.0 * (me2 + 2.0 * me2 * rho), (2.0 * (me1 * (14.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho - 3.0 * me3 * (14.0 + rho) + 6.0 * (me4 + 2.0 * me4 * rho)))) / 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, me1 + 6.0 * me3, 3.0 * me1, -0.5 * ((me1 + 6.0 * me3) * sqrtrho), (-3.0 * me1 * sqrtrho) / 2.0, -0.25 * ((me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) * sqrtrho), (-3.0 * (me5 - 2.0 * me6) * sqrtrho) / 4.0, (-me5 + 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) / 4.0, (-3.0 * (me5 - 2.0 * me6)) / 4.0, 3.0 * (me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)), 9.0 * (me5 - 2.0 * me6)},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0 * me1, me1 + 6.0 * me3, (-3.0 * me1 * sqrtrho) / 2.0, -0.5 * ((me1 + 6.0 * me3) * sqrtrho), (-3.0 * (me5 - 2.0 * me6) * sqrtrho) / 4.0, -0.25 * ((me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) * sqrtrho), (-3.0 * (me5 - 2.0 * me6)) / 4.0, (-me5 + 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) / 4.0, 9.0 * (me5 - 2.0 * me6), 3.0 * (me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8))},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5 * ((me1 + 6.0 * me3) * sqrtrho), (-3.0 * me1 * sqrtrho) / 2.0, (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho) + 6.0 * (me3 + 2.0 * me3 * rho - 2.0 * me4 * (2.0 + rho))) / 6.0, me1 * (0.5 + rho) - me2 * (2.0 + rho), (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho + 3.0 * (me7 + 2.0 * me8 - me7 * rho + 4.0 * me8 * rho))) / 12.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho)) / 4.0, -0.5 * ((me6 + 6.0 * me8) * sqrtrho), (-3.0 * me6 * sqrtrho) / 2.0, -2.0 * (me5 - me6 + 6.0 * me7 - 6.0 * me8) * sqrtrho, -6.0 * (me5 - me6) * sqrtrho},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (-3.0 * me1 * sqrtrho) / 2.0, -0.5 * ((me1 + 6.0 * me3) * sqrtrho), me1 * (0.5 + rho) - me2 * (2.0 + rho), (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho) + 6.0 * (me3 + 2.0 * me3 * rho - 2.0 * me4 * (2.0 + rho))) / 6.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho)) / 4.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho + 3.0 * (me7 + 2.0 * me8 - me7 * rho + 4.0 * me8 * rho))) / 12.0, (-3.0 * me6 * sqrtrho) / 2.0, -0.5 * ((me6 + 6.0 * me8) * sqrtrho), -6.0 * (me5 - me6) * sqrtrho, -2.0 * (me5 - me6 + 6.0 * me7 - 6.0 * me8) * sqrtrho},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25 * ((me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) * sqrtrho), (-3.0 * (me5 - 2.0 * me6) * sqrtrho) / 4.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho + 3.0 * (me7 + 2.0 * me8 - me7 * rho + 4.0 * me8 * rho))) / 12.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho)) / 4.0, (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho) + 6.0 * (me3 + 2.0 * me3 * rho - 2.0 * me4 * (2.0 + rho))) / 24.0, (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho)) / 8.0, ((me1 - 2.0 * (me2 - 3.0 * me3 + 6.0 * me4)) * sqrtrho) / 8.0, (3.0 * (me1 - 2.0 * me2) * sqrtrho) / 8.0, -0.5 * ((me1 + 2.0 * (me2 + 3.0 * me3 + 6.0 * me4)) * sqrtrho), (-3.0 * (me1 + 2.0 * me2) * sqrtrho) / 2.0},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (-3.0 * (me5 - 2.0 * me6) * sqrtrho) / 4.0, -0.25 * ((me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) * sqrtrho), (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho)) / 4.0, (me5 * (-1.0 + rho) - 2.0 * (me6 + 2.0 * me6 * rho + 3.0 * (me7 + 2.0 * me8 - me7 * rho + 4.0 * me8 * rho))) / 12.0, (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho)) / 8.0, (me1 + 2.0 * me1 * rho - 2.0 * me2 * (2.0 + rho) + 6.0 * (me3 + 2.0 * me3 * rho - 2.0 * me4 * (2.0 + rho))) / 24.0, (3.0 * (me1 - 2.0 * me2) * sqrtrho) / 8.0, ((me1 - 2.0 * (me2 - 3.0 * me3 + 6.0 * me4)) * sqrtrho) / 8.0, (-3.0 * (me1 + 2.0 * me2) * sqrtrho) / 2.0, -0.5 * ((me1 + 2.0 * (me2 + 3.0 * me3 + 6.0 * me4)) * sqrtrho)},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (-me5 + 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) / 4.0, (-3.0 * (me5 - 2.0 * me6)) / 4.0, -0.5 * ((me6 + 6.0 * me8) * sqrtrho), (-3.0 * me6 * sqrtrho) / 2.0, ((me1 - 2.0 * (me2 - 3.0 * me3 + 6.0 * me4)) * sqrtrho) / 8.0, (3.0 * (me1 - 2.0 * me2) * sqrtrho) / 8.0, (me1 * (2.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho - 3.0 * me3 * (2.0 + rho) + 6.0 * (me4 + 2.0 * me4 * rho))) / 24.0, (me1 * (2.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho)) / 8.0, (me1 * (-4.0 + rho) - 2.0 * (me2 - 3.0 * me3 * (-4.0 + rho) + 2.0 * me2 * rho + 6.0 * (me4 + 2.0 * me4 * rho))) / 6.0, (me1 * (-4.0 + rho)) / 2.0 - me2 * (1.0 + 2.0 * rho)},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (-3.0 * (me5 - 2.0 * me6)) / 4.0, (-me5 + 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)) / 4.0, (-3.0 * me6 * sqrtrho) / 2.0, -0.5 * ((me6 + 6.0 * me8) * sqrtrho), (3.0 * (me1 - 2.0 * me2) * sqrtrho) / 8.0, ((me1 - 2.0 * (me2 - 3.0 * me3 + 6.0 * me4)) * sqrtrho) / 8.0, (me1 * (2.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho)) / 8.0, (me1 * (2.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho - 3.0 * me3 * (2.0 + rho) + 6.0 * (me4 + 2.0 * me4 * rho))) / 24.0, (me1 * (-4.0 + rho)) / 2.0 - me2 * (1.0 + 2.0 * rho), (me1 * (-4.0 + rho) - 2.0 * (me2 - 3.0 * me3 * (-4.0 + rho) + 2.0 * me2 * rho + 6.0 * (me4 + 2.0 * me4 * rho))) / 6.0},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0 * (me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)), 9.0 * (me5 - 2.0 * me6), -2.0 * (me5 - me6 + 6.0 * me7 - 6.0 * me8) * sqrtrho, -6.0 * (me5 - me6) * sqrtrho, -0.5 * ((me1 + 2.0 * (me2 + 3.0 * me3 + 6.0 * me4)) * sqrtrho), (-3.0 * (me1 + 2.0 * me2) * sqrtrho) / 2.0, (me1 * (-4.0 + rho) - 2.0 * (me2 - 3.0 * me3 * (-4.0 + rho) + 2.0 * me2 * rho + 6.0 * (me4 + 2.0 * me4 * rho))) / 6.0, (me1 * (-4.0 + rho)) / 2.0 - me2 * (1.0 + 2.0 * rho), (2.0 * (me1 * (14.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho - 3.0 * me3 * (14.0 + rho) + 6.0 * (me4 + 2.0 * me4 * rho)))) / 3.0, 2.0 * me1 * (14.0 + rho) - 4.0 * (me2 + 2.0 * me2 * rho)},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0 * (me5 - 2.0 * me6), 3.0 * (me5 - 2.0 * (me6 - 3.0 * me7 + 6.0 * me8)), -6.0 * (me5 - me6) * sqrtrho, -2.0 * (me5 - me6 + 6.0 * me7 - 6.0 * me8) * sqrtrho, (-3.0 * (me1 + 2.0 * me2) * sqrtrho) / 2.0, -0.5 * ((me1 + 2.0 * (me2 + 3.0 * me3 + 6.0 * me4)) * sqrtrho), (me1 * (-4.0 + rho)) / 2.0 - me2 * (1.0 + 2.0 * rho), (me1 * (-4.0 + rho) - 2.0 * (me2 - 3.0 * me3 * (-4.0 + rho) + 2.0 * me2 * rho + 6.0 * (me4 + 2.0 * me4 * rho))) / 6.0, 2.0 * me1 * (14.0 + rho) - 4.0 * (me2 + 2.0 * me2 * rho), (2.0 * (me1 * (14.0 + rho) - 2.0 * (me2 + 2.0 * me2 * rho - 3.0 * me3 * (14.0 + rho) + 6.0 * (me4 + 2.0 * me4 * rho)))) / 3.0}
            }};

            return result;
        }

        // cf. [LMPR:2022A], eqs. (2.28) to (2.33)
        // provided by A. Rusov as a Mathematica expression
        std::array<std::array<complex<double>, 20u>, 20u>
        A_weak_exchange(const double & /*mu*/, const double & sqrtrho) const
        {
            const double rho = sqrtrho * sqrtrho;

            // Matrix elements of operators [bbar Gamma q] [qbar Gamma b]
            // are only known from HQET sum rules. For the time being use,
            // constant values for these matrix elements, mostly following
            // [LMPR:2022A]. For the matrix elements vanishing in vacuum
            // insertion approximation, we use 10% as the naive bag factor.
            const double me  =  f_B * f_B * m_B * m_B;
            const double me1 =  1.0 * me;
            const double me2 =  1.0 * me;
            const double me3 =  0.1 * me;
            const double me4 =  0.1 * me;
            const double me5 = -1.0 * me;
            const double me6 = -1.0 * me;
            const double me7 =  0.1 * me;
            const double me8 =  0.1 * me;

            std::array<std::array<complex<double>, 20u>, 20u> result
            {{
                {(-(me1 * (2.0 + rho)) + 2.0 * (me2 + 2.0 * me2 * rho - 3 * me3 * (2 + rho) + 6 * (me4 + 2 * me4 * rho))) / 6., me2 + 2 * me2 * rho - (me1 * (2 + rho)) / 2., -((me2 + 6 * me4) * sqrtrho), -3 * me2 * sqrtrho, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho - 3 * me7 * (2 + rho) + 6 * (me8 + 2 * me8 * rho))) / 12., (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho)) / 4., ((me5 - 2 * (me6 - 3 * me7 + 6 * me8)) * sqrtrho) / 4., (3 * (me5 - 2 * me6) * sqrtrho) / 4., -((me5 + 2 * (me6 + 3 * me7 + 6 * me8)) * sqrtrho), -3 * (me5 + 2 * me6) * sqrtrho},
                {me2 + 2 * me2 * rho - (me1 * (2 + rho)) / 2., (-3 * me1 * (2 + rho)) / 2. + 3 * me2 * (1 + 2 * rho), -3 * me2 * sqrtrho, -9 * me2 * sqrtrho, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho)) / 4., (3 * (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho))) / 4., (3 * (me5 - 2 * me6) * sqrtrho) / 4., (9 * (me5 - 2 * me6) * sqrtrho) / 4., -3 * (me5 + 2 * me6) * sqrtrho, -9 * (me5 + 2 * me6) * sqrtrho},
                {-((me2 + 6 * me4) * sqrtrho), -3 * me2 * sqrtrho, 2 * (me2 + 6 * me4), 6 * me2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ((me6 + 6 * me8) * sqrtrho) / 2., (3 * me6 * sqrtrho) / 2., (me6 + 6 * me8) / 2., (3 * me6) / 2., 6 * (me6 + 6 * me8), 18 * me6},
                {-3 * me2 * sqrtrho, -9 * me2 * sqrtrho, 6 * me2, 18 * me2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3 * me6 * sqrtrho) / 2., (9 * me6 * sqrtrho) / 2., (3 * me6) / 2., (9 * me6) / 2., 18 * me6, 54 * me6},
                {0, 0, 0, 0, (-(me1 * (2 + rho)) + 2 * (me2 + 2 * me2 * rho - 3 * me3 * (2 + rho) + 6 * (me4 + 2 * me4 * rho))) / 24., (-(me1 * (2 + rho)) + 2 * me2 * (1 + 2 * rho)) / 8., -0.125 * ((me1 - 2 * (me2 - 3 * me3 + 6 * me4)) * sqrtrho), (-3 * (me1 - 2 * me2) * sqrtrho) / 8., ((me1 + 2 * (me2 + 3 * me3 + 6 * me4)) * sqrtrho) / 2., (3 * (me1 + 2 * me2) * sqrtrho) / 2., (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho - 3 * me7 * (2 + rho) + 6 * (me8 + 2 * me8 * rho))) / 12., (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho)) / 4., ((me6 + 6 * me8) * sqrtrho) / 2., (3 * me6 * sqrtrho) / 2., 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, (-(me1 * (2 + rho)) + 2 * me2 * (1 + 2 * rho)) / 8., (-3 * (me1 * (2 + rho) - 2 * (me2 + 2 * me2 * rho))) / 8., (-3 * (me1 - 2 * me2) * sqrtrho) / 8., (-9 * (me1 - 2 * me2) * sqrtrho) / 8., (3 * (me1 + 2 * me2) * sqrtrho) / 2., (9 * (me1 + 2 * me2) * sqrtrho) / 2., (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho)) / 4., (3 * (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho))) / 4., (3 * me6 * sqrtrho) / 2., (9 * me6 * sqrtrho) / 2., 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, -0.125 * ((me1 - 2 * (me2 - 3 * me3 + 6 * me4)) * sqrtrho), (-3 * (me1 - 2 * me2) * sqrtrho) / 8., (2 * me2 * (2 + rho) - me1 * (1 + 2 * rho) - 6 * (me3 + 2 * me3 * rho - 2 * me4 * (2 + rho))) / 24., (2 * me2 * (2 + rho) - me1 * (1 + 2 * rho)) / 8., (me1 - 2 * me2 * (-4 + rho) + 2 * me1 * rho + 6 * (me3 + 8 * me4 + 2 * me3 * rho - 2 * me4 * rho)) / 6., -(me2 * (-4 + rho)) + me1 * (0.5 + rho), ((me5 - 2 * (me6 - 3 * me7 + 6 * me8)) * sqrtrho) / 4., (3 * (me5 - 2 * me6) * sqrtrho) / 4., (me6 + 6 * me8) / 2., (3 * me6) / 2., 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, (-3 * (me1 - 2 * me2) * sqrtrho) / 8., (-9 * (me1 - 2 * me2) * sqrtrho) / 8., (2 * me2 * (2 + rho) - me1 * (1 + 2 * rho)) / 8., (-3 * (me1 + 2 * me1 * rho - 2 * me2 * (2 + rho))) / 8., -(me2 * (-4 + rho)) + me1 * (0.5 + rho), (3 * (me1 - 2 * me2 * (-4 + rho) + 2 * me1 * rho)) / 2., (3 * (me5 - 2 * me6) * sqrtrho) / 4., (9 * (me5 - 2 * me6) * sqrtrho) / 4., (3 * me6) / 2., (9 * me6) / 2., 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, ((me1 + 2 * (me2 + 3 * me3 + 6 * me4)) * sqrtrho) / 2., (3 * (me1 + 2 * me2) * sqrtrho) / 2., (me1 - 2 * me2 * (-4 + rho) + 2 * me1 * rho + 6 * (me3 + 8 * me4 + 2 * me3 * rho - 2 * me4 * rho)) / 6., -(me2 * (-4 + rho)) + me1 * (0.5 + rho), (-2 * (me1 + 2 * me1 * rho - 2 * me2 * (14 + rho) - 12 * me4 * (14 + rho) + 6 * me3 * (1 + 2 * rho))) / 3., 4 * me2 * (14 + rho) - 2 * me1 * (1 + 2 * rho), -((me5 + 2 * (me6 + 3 * me7 + 6 * me8)) * sqrtrho), -3 * (me5 + 2 * me6) * sqrtrho, 6 * (me6 + 6 * me8), 18 * me6, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, (3 * (me1 + 2 * me2) * sqrtrho) / 2., (9 * (me1 + 2 * me2) * sqrtrho) / 2., -(me2 * (-4 + rho)) + me1 * (0.5 + rho), (3 * (me1 - 2 * me2 * (-4 + rho) + 2 * me1 * rho)) / 2., 4 * me2 * (14 + rho) - 2 * me1 * (1 + 2 * rho), -6 * (me1 + 2 * me1 * rho - 2 * me2 * (14 + rho)), -3 * (me5 + 2 * me6) * sqrtrho, -9 * (me5 + 2 * me6) * sqrtrho, 18 * me6, 54 * me6, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho - 3 * me7 * (2 + rho) + 6 * (me8 + 2 * me8 * rho))) / 12., (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho)) / 4., ((me5 - 2 * (me6 - 3 * me7 + 6 * me8)) * sqrtrho) / 4., (3 * (me5 - 2 * me6) * sqrtrho) / 4., -((me5 + 2 * (me6 + 3 * me7 + 6 * me8)) * sqrtrho), -3 * (me5 + 2 * me6) * sqrtrho, (-(me1 * (2 + rho)) + 2 * (me2 + 2 * me2 * rho - 3 * me3 * (2 + rho) + 6 * (me4 + 2 * me4 * rho))) / 6., me2 + 2 * me2 * rho - (me1 * (2 + rho)) / 2., -((me2 + 6 * me4) * sqrtrho), -3 * me2 * sqrtrho, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho)) / 4., (3 * (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho))) / 4., (3 * (me5 - 2 * me6) * sqrtrho) / 4., (9 * (me5 - 2 * me6) * sqrtrho) / 4., -3 * (me5 + 2 * me6) * sqrtrho, -9 * (me5 + 2 * me6) * sqrtrho, me2 + 2 * me2 * rho - (me1 * (2 + rho)) / 2., (-3 * me1 * (2 + rho)) / 2. + 3 * me2 * (1 + 2 * rho), -3 * me2 * sqrtrho, -9 * me2 * sqrtrho, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, ((me6 + 6 * me8) * sqrtrho) / 2., (3 * me6 * sqrtrho) / 2., (me6 + 6 * me8) / 2., (3 * me6) / 2., 6 * (me6 + 6 * me8), 18 * me6, -((me2 + 6 * me4) * sqrtrho), -3 * me2 * sqrtrho, 2 * (me2 + 6 * me4), 6 * me2, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, (3 * me6 * sqrtrho) / 2., (9 * me6 * sqrtrho) / 2., (3 * me6) / 2., (9 * me6) / 2., 18 * me6, 54 * me6, -3 * me2 * sqrtrho, -9 * me2 * sqrtrho, 6 * me2, 18 * me2, 0, 0, 0, 0, 0, 0},
                {(me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho - 3 * me7 * (2 + rho) + 6 * (me8 + 2 * me8 * rho))) / 12., (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho)) / 4., ((me6 + 6 * me8) * sqrtrho) / 2., (3 * me6 * sqrtrho) / 2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-(me1 * (2 + rho)) + 2 * (me2 + 2 * me2 * rho - 3 * me3 * (2 + rho) + 6 * (me4 + 2 * me4 * rho))) / 24., (-(me1 * (2 + rho)) + 2 * me2 * (1 + 2 * rho)) / 8., -0.125 * ((me1 - 2 * (me2 - 3 * me3 + 6 * me4)) * sqrtrho), (-3 * (me1 - 2 * me2) * sqrtrho) / 8., ((me1 + 2 * (me2 + 3 * me3 + 6 * me4)) * sqrtrho) / 2., (3 * (me1 + 2 * me2) * sqrtrho) / 2.},
                {(me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho)) / 4., (3 * (me5 * (2 + rho) - 2 * (me6 + 2 * me6 * rho))) / 4., (3 * me6 * sqrtrho) / 2., (9 * me6 * sqrtrho) / 2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-(me1 * (2 + rho)) + 2 * me2 * (1 + 2 * rho)) / 8., (-3 * (me1 * (2 + rho) - 2 * (me2 + 2 * me2 * rho))) / 8., (-3 * (me1 - 2 * me2) * sqrtrho) / 8., (-9 * (me1 - 2 * me2) * sqrtrho) / 8., (3 * (me1 + 2 * me2) * sqrtrho) / 2., (9 * (me1 + 2 * me2) * sqrtrho) / 2.},
                {((me5 - 2 * (me6 - 3 * me7 + 6 * me8)) * sqrtrho) / 4., (3 * (me5 - 2 * me6) * sqrtrho) / 4., (me6 + 6 * me8) / 2., (3 * me6) / 2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.125 * ((me1 - 2 * (me2 - 3 * me3 + 6 * me4)) * sqrtrho), (-3 * (me1 - 2 * me2) * sqrtrho) / 8., (2 * me2 * (2 + rho) - me1 * (1 + 2 * rho) - 6 * (me3 + 2 * me3 * rho - 2 * me4 * (2 + rho))) / 24., (2 * me2 * (2 + rho) - me1 * (1 + 2 * rho)) / 8., (me1 - 2 * me2 * (-4 + rho) + 2 * me1 * rho + 6 * (me3 + 8 * me4 + 2 * me3 * rho - 2 * me4 * rho)) / 6., -(me2 * (-4 + rho)) + me1 * (0.5 + rho)},
                {(3 * (me5 - 2 * me6) * sqrtrho) / 4., (9 * (me5 - 2 * me6) * sqrtrho) / 4., (3 * me6) / 2., (9 * me6) / 2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-3 * (me1 - 2 * me2) * sqrtrho) / 8., (-9 * (me1 - 2 * me2) * sqrtrho) / 8., (2 * me2 * (2 + rho) - me1 * (1 + 2 * rho)) / 8., (-3 * (me1 + 2 * me1 * rho - 2 * me2 * (2 + rho))) / 8., -(me2 * (-4 + rho)) + me1 * (0.5 + rho), (3 * (me1 - 2 * me2 * (-4 + rho) + 2 * me1 * rho)) / 2.},
                {-((me5 + 2 * (me6 + 3 * me7 + 6 * me8)) * sqrtrho), -3 * (me5 + 2 * me6) * sqrtrho, 6 * (me6 + 6 * me8), 18 * me6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ((me1 + 2 * (me2 + 3 * me3 + 6 * me4)) * sqrtrho) / 2., (3 * (me1 + 2 * me2) * sqrtrho) / 2., (me1 - 2 * me2 * (-4 + rho) + 2 * me1 * rho + 6 * (me3 + 8 * me4 + 2 * me3 * rho - 2 * me4 * rho)) / 6., -(me2 * (-4 + rho)) + me1 * (0.5 + rho), (-2 * (me1 + 2 * me1 * rho - 2 * me2 * (14 + rho) - 12 * me4 * (14 + rho) + 6 * me3 * (1 + 2 * rho))) / 3., 4 * me2 * (14 + rho) - 2 * me1 * (1 + 2 * rho)},
                {-3 * (me5 + 2 * me6) * sqrtrho, -9 * (me5 + 2 * me6) * sqrtrho, 18 * me6, 54 * me6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3 * (me1 + 2 * me2) * sqrtrho) / 2., (9 * (me1 + 2 * me2) * sqrtrho) / 2., -(me2 * (-4 + rho)) + me1 * (0.5 + rho), (3 * (me1 - 2 * me2 * (-4 + rho) + 2 * me1 * rho)) / 2., 4 * me2 * (14 + rho) - 2 * me1 * (1 + 2 * rho), -6 * (me1 + 2 * me1 * rho - 2 * me2 * (14 + rho))}
            }};

            return result;
        }

        double decay_width_dbcu_dim6_lo() const
        {
            const double mu      = mu_dbcu();
            const double m_b     = model->m_b_msbar(mu);
            const double sqrtrho = model->m_c_msbar(mu) / m_b;
            const double rho     = sqrtrho * sqrtrho;
            const auto   wc      = model->wet_dbcu(mu);
            const auto   A_pi_cd = complex<double>(switch_pauli_interference_dbcu) * A_pauli_interference(mu, sqrtrho);
            const auto   A_we_cu = complex<double>(switch_weak_exchange_dbcu)      * A_weak_exchange(mu, sqrtrho);

            // transforming the Wilson coefficients from the EOS basis
            // to the basis used in [LMPR:2022A], eqs. (2.1) to (2.6).
            const std::array<std::complex<double>, 20u> C
            {
                wc.c2() / 2.0 + 8.0 * wc.c4(),
                wc.c1() - wc.c2() / 6.0 + 16.0 * wc.c3() - (8.0 * wc.c4()) / 3.0,
                -4.0 * wc.c10() - wc.c6() / 4.0,
                1.0 / 12.0 * (16.0 * wc.c10() - 6.0 * wc.c5() + wc.c6() - 96.0 * wc.c9()),
                -wc.c2() - 4.0 * wc.c4(),
                1.0 / 3.0 * (-6.0 * wc.c1() + wc.c2() + 4.0 * (wc.c4() - 6.0 * wc.c3())),
                32.0 * wc.c10() - wc.c6() / 4.0 - 3.0 * wc.c8(),
                -((32.0 * wc.c10()) / 3.0) - wc.c5() / 2.0 + wc.c6() / 12.0 - 6.0 * wc.c7() + wc.c8() + 64.0 * wc.c9(),
                -8.0 * wc.c10() - wc.c6() / 16.0 + wc.c8() / 4.0,
                1.0 / 48.0 * (128.0 * wc.c10() - 6.0 * wc.c5() + wc.c6() + 24.0 * wc.c7() - 4.0 * wc.c8() - 768.0 * wc.c9()),
                wc.c2p() / 2.0 + 8.0 * wc.c4p(), wc.c1p() - wc.c2p() / 6.0 + 16.0 * wc.c3p() - (8.0 * wc.c4p()) / 3.0,
                -4.0 * wc.c10p() - wc.c6p() / 4.0, 1.0 / 12.0 * (16.0 * wc.c10p() - 6.0 * wc.c5p() + wc.c6p() - 96.0 * wc.c9p()),
                -wc.c2p() - 4.0 * wc.c4p(),
                1.0 / 3.0 * (-6.0 * wc.c1p() + wc.c2p() + 4.0 * (wc.c4p() - 6.0 * wc.c3p())), 32.0 * wc.c10p() - wc.c6p() / 4.0 - 3.0 * wc.c8p(),
                -((32.0 * wc.c10p()) / 3.0) - wc.c5p() / 2.0 + wc.c6p() / 12.0 - 6.0 * wc.c7p() + wc.c8p() + 64.0 * wc.c9p(),
                -8.0 * wc.c10p() - wc.c6p() / 16.0 + wc.c8p() / 4.0, 1.0 / 48.0 * (128.0 * wc.c10p() - 6.0 * wc.c5p() + wc.c6p() + 24.0 * wc.c7p() - 4.0 * wc.c8p() - 768.0 * wc.c9p())
            };
            std::array<std::complex<double>, 20u> Cconj;
            complex<double> (*conj)(const std::complex<double> &) = &std::conj<double>;
            std::transform(C.cbegin(), C.cend(), Cconj.begin(), conj);

            const double ckm = abs(model->ckm_ud() * model->ckm_cb());

            const auto A = A_pi_cd + A_we_cu;

            return power_of<2>(g_fermi * m_b * ckm * (1.0 - rho)) / (12.0 * m_B * M_PI * hbar) * real(dot(Cconj, (A * C))) * 1.0e-12;
        }

        double decay_width_sbcu_dim6_lo() const
        {
            const double mu      = mu_sbcu();
            const double m_b     = model->m_b_msbar(mu);
            const double sqrtrho = model->m_c_msbar(mu) / m_b;
            const double rho     = sqrtrho * sqrtrho;
            const auto   wc      = model->wet_sbcu(mu);
            const auto   A_pi_cs = complex<double>(switch_pauli_interference_sbcu) * A_pauli_interference(mu, sqrtrho);
            const auto   A_we_cu = complex<double>(switch_weak_exchange_sbcu)      * A_weak_exchange(mu, sqrtrho);

            // transforming the Wilson coefficients from the EOS basis
            // to the basis used in [LMPR:2022A], eqs. (2.1) to (2.6).
            const std::array<std::complex<double>, 20u> C
            {
                wc.c2() / 2.0 + 8.0 * wc.c4(),
                wc.c1() - wc.c2() / 6.0 + 16.0 * wc.c3() - (8.0 * wc.c4()) / 3.0,
                -4.0 * wc.c10() - wc.c6() / 4.0,
                1.0 / 12.0 * (16.0 * wc.c10() - 6.0 * wc.c5() + wc.c6() - 96.0 * wc.c9()),
                -wc.c2() - 4.0 * wc.c4(),
                1.0 / 3.0 * (-6.0 * wc.c1() + wc.c2() + 4.0 * (wc.c4() - 6.0 * wc.c3())),
                32.0 * wc.c10() - wc.c6() / 4.0 - 3.0 * wc.c8(),
                -((32.0 * wc.c10()) / 3.0) - wc.c5() / 2.0 + wc.c6() / 12.0 - 6.0 * wc.c7() + wc.c8() + 64.0 * wc.c9(),
                -8.0 * wc.c10() - wc.c6() / 16.0 + wc.c8() / 4.0,
                1.0 / 48.0 * (128.0 * wc.c10() - 6.0 * wc.c5() + wc.c6() + 24.0 * wc.c7() - 4.0 * wc.c8() - 768.0 * wc.c9()),
                wc.c2p() / 2.0 + 8.0 * wc.c4p(), wc.c1p() - wc.c2p() / 6.0 + 16.0 * wc.c3p() - (8.0 * wc.c4p()) / 3.0,
                -4.0 * wc.c10p() - wc.c6p() / 4.0, 1.0 / 12.0 * (16.0 * wc.c10p() - 6.0 * wc.c5p() + wc.c6p() - 96.0 * wc.c9p()),
                -wc.c2p() - 4.0 * wc.c4p(),
                1.0 / 3.0 * (-6.0 * wc.c1p() + wc.c2p() + 4.0 * (wc.c4p() - 6.0 * wc.c3p())), 32.0 * wc.c10p() - wc.c6p() / 4.0 - 3.0 * wc.c8p(),
                -((32.0 * wc.c10p()) / 3.0) - wc.c5p() / 2.0 + wc.c6p() / 12.0 - 6.0 * wc.c7p() + wc.c8p() + 64.0 * wc.c9p(),
                -8.0 * wc.c10p() - wc.c6p() / 16.0 + wc.c8p() / 4.0, 1.0 / 48.0 * (128.0 * wc.c10p() - 6.0 * wc.c5p() + wc.c6p() + 24.0 * wc.c7p() - 4.0 * wc.c8p() - 768.0 * wc.c9p())
            };
            std::array<std::complex<double>, 20u> Cconj;
            complex<double> (*conj)(const std::complex<double> &) = &std::conj<double>;
            std::transform(C.cbegin(), C.cend(), Cconj.begin(), conj);

            const double ckm = abs(model->ckm_us() * model->ckm_cb());

            const auto A = A_pi_cs + A_we_cu;

            return power_of<2>(g_fermi * m_b * ckm * (1.0 - rho)) / (12.0 * m_B * M_PI * hbar) * real(dot(Cconj, (A * C))) * 1.0e-12;
        }

        double decay_width_dbcu_dim3_lo() const
        {
            const double mu  = mu_dbcu();
            const double m_b = model->m_b_msbar(mu);
            const double z   = model->m_c_msbar(mu) / m_b, log_z = std::log(z),
                         z2  = z * z, z4 = z2 * z2, z6 = z4 * z2, z8 = z4 * z4;
            const double ckm = abs(model->ckm_ud() * model->ckm_cb());
            const auto   wc  = model->wet_dbcu(mu);

            std::array<std::array<complex<double>, 5u>, 5u> G
            {{
                {
                    -6.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    -60.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    -6.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -36.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    336.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z)
                },
                {
                    -60.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    -816.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    -60.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -144.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    768.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z)
                },
                {
                    -6.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -60.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -3.0 / 2.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    0,
                    -60.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z)
                },
                {
                    -36.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -144.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    0,
                    -36.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    576.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z)
                },
                {
                    336.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    768.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -60.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    576.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    -12480.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z)
                }
            }};

            const std::array<std::complex<double>, 5u> C_sing
            {
                wc.c1(),  wc.c3(),  wc.c5(),  wc.c7(),  wc.c9()
            };

            std::array<std::complex<double>, 5u> Cconj_sing;
            complex<double> (*conj_sing)(const std::complex<double> &) = &std::conj<double>;
            std::transform(C_sing.cbegin(), C_sing.cend(), Cconj_sing.begin(), conj_sing);

            const std::array<std::complex<double>, 5u> C_sing_p
            {
                wc.c1p(), wc.c3p(), wc.c5p(), wc.c7p(), wc.c9p()
            };

            std::array<std::complex<double>, 5u> Cconj_sing_p;
            complex<double> (*conj_sing_p)(const std::complex<double> &) = &std::conj<double>;
            std::transform(C_sing_p.cbegin(), C_sing_p.cend(), Cconj_sing_p.begin(), conj_sing_p);

            const std::array<std::complex<double>, 5u> C_oct
            {
                wc.c2(),  wc.c4(),  wc.c6(),  wc.c8(),  wc.c10()
            };

            std::array<std::complex<double>, 5u> Cconj_oct;
            complex<double> (*conj_oct)(const std::complex<double> &) = &std::conj<double>;
            std::transform(C_oct.cbegin(), C_oct.cend(), Cconj_oct.begin(), conj_oct);

            const std::array<std::complex<double>, 5u> C_oct_p
            {
                wc.c2p(), wc.c4p(), wc.c6p(), wc.c8p(), wc.c10p()
            };

            std::array<std::complex<double>, 5u> Cconj_oct_p;
            complex<double> (*conj_oct_p)(const std::complex<double> &) = &std::conj<double>;
            std::transform(C_oct_p.cbegin(), C_oct_p.cend(), Cconj_oct_p.begin(), conj_oct_p);

            const double Gamma_0      = g_fermi * g_fermi * ckm * ckm * power_of<5>(m_b) / (192.0 * power_of<3>(M_PI));
            const double Gamma_sing   = real(dot(Cconj_sing, (G * C_sing)));
            const double Gamma_sing_p = real(dot(Cconj_sing_p, (G * C_sing_p)));
            const double Gamma_oct    = 2.0 / 9.0 * real(dot(Cconj_oct, (G * C_oct)));
            const double Gamma_oct_p  = 2.0 / 9.0 * real(dot(Cconj_oct_p, (G * C_oct_p)));

            return Gamma_0 * (Gamma_sing + Gamma_sing_p + Gamma_oct + Gamma_oct_p) * convert_GeV_to_inv_ps();
        }

        double decay_width_sbcu_dim3_lo() const
        {
            const double mu  = mu_sbcu();
            const double m_b = model->m_b_msbar(mu);
            const double z   = model->m_c_msbar(mu) / m_b, log_z = std::log(z),
                         z2  = z * z, z4 = z2 * z2, z6 = z4 * z2, z8 = z4 * z4;
            const double ckm = abs(model->ckm_us() * model->ckm_cb());
            const auto   wc  = model->wet_sbcu(mu);

            std::array<std::array<complex<double>, 5u>, 5u> G
            {{
                {
                    -6.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    -60.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    -6.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -36.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    336.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z)
                },
                {
                    -60.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    -816.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    -60.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -144.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    768.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z)
                },
                {
                    -6.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -60.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -3.0 / 2.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    0,
                    -60.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z)
                },
                {
                    -36.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -144.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    0,
                    -36.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    576.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z)
                },
                {
                    336.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    768.0 * z * (-1.0 - 9.0 * z2 + 9.0 * z4 + z6 - 12.0 * (z2 + z4) * log_z),
                    -60.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    576.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z),
                    -12480.0 * (-1.0 + 8.0 * z2 - 8.0 * z6 + z8 + 24.0 * z4 * log_z)
                }
            }};

            // Singlet Wilson coefficients
            const std::array<std::complex<double>, 5u> C_singlet
            {
                wc.c1(),  wc.c3(),  wc.c5(),  wc.c7(),  wc.c9()
            };
            std::array<std::complex<double>, 5u> Cconj_singlet;
            std::transform(C_singlet.cbegin(), C_singlet.cend(), Cconj_singlet.begin(), static_cast<complex<double> (*)(const std::complex<double> &)>(&std::conj<double>));

            // Singlet Wilson coefficients (helicity-flipped)
            const std::array<std::complex<double>, 5u> C_singlet_p
            {
                wc.c1p(), wc.c3p(), wc.c5p(), wc.c7p(), wc.c9p()
            };
            std::array<std::complex<double>, 5u> Cconj_singlet_p;
            std::transform(C_singlet_p.cbegin(), C_singlet_p.cend(), Cconj_singlet_p.begin(), static_cast<complex<double> (*)(const std::complex<double> &)>(&std::conj<double>));

            // Octet Wilson coefficients
            const std::array<std::complex<double>, 5u> C_octet
            {
                wc.c2(),  wc.c4(),  wc.c6(),  wc.c8(),  wc.c10()
            };
            std::array<std::complex<double>, 5u> Cconj_octet;
            std::transform(C_octet.cbegin(), C_octet.cend(), Cconj_octet.begin(), static_cast<complex<double> (*)(const std::complex<double> &)>(&std::conj<double>));

            // Octet Wilson coefficients (helicity-flipped)
            const std::array<std::complex<double>, 5u> C_octet_p
            {
                wc.c2p(), wc.c4p(), wc.c6p(), wc.c8p(), wc.c10p()
            };
            std::array<std::complex<double>, 5u> Cconj_octet_p;
            std::transform(C_octet_p.cbegin(), C_octet_p.cend(), Cconj_octet_p.begin(), static_cast<complex<double> (*)(const std::complex<double> &)>(&std::conj<double>));

            const double Gamma_0      = g_fermi * g_fermi * ckm * ckm * power_of<5>(m_b) / (192.0 * power_of<3>(M_PI));
            const double Gamma_sing   = real(dot(Cconj_singlet,   (G * C_singlet)));
            const double Gamma_sing_p = real(dot(Cconj_singlet_p, (G * C_singlet_p)));
            const double Gamma_oct    = 2.0 / 9.0 * real(dot(Cconj_octet,   (G * C_octet)));
            const double Gamma_oct_p  = 2.0 / 9.0 * real(dot(Cconj_octet_p, (G * C_octet_p)));

            return Gamma_0 * (Gamma_sing + Gamma_sing_p + Gamma_oct + Gamma_oct_p) * convert_GeV_to_inv_ps();
        }
    };

    const std::vector<OptionSpecification>
    Implementation<Lifetime>::options
    {
        Model::option_specification(),
        { "q", { "u", "d", "s" }, "" }
    };

    Lifetime::Lifetime(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<Lifetime>(new Implementation<Lifetime>(parameters, options, *this))
    {
    }

    Lifetime::~Lifetime()
    {
    }

    double
    Lifetime::decay_width_dbcu_dim3_lo() const
    {
        return _imp->decay_width_dbcu_dim3_lo();
    }

    double
    Lifetime::decay_width_sbcu_dim3_lo() const
    {
        return _imp->decay_width_sbcu_dim3_lo();
    }

    /*
    * The partial decay width computed with ``decay_width_dbcu`` and ``decay_width_sbcu`` can be negative.
    * This issue arises, since
    *  a) the computation is only a partial one (i.e., further terms can arise);
    *  b) the computation is done to leading-order in alpha_s.
    * For the purpose of determining bounds on the size of the WET Wilson coefficients,
    * we would like to compare this partial decay width with the measured total decay
    * width. Returning the absolute value suffices for this comparison.
    */
    double
    Lifetime::decay_width_dbcu_dim6_lo() const
    {
        return std::fabs(_imp->decay_width_dbcu_dim6_lo());
    }

    double
    Lifetime::decay_width_sbcu_dim6_lo() const
    {
        return std::fabs(_imp->decay_width_sbcu_dim6_lo());
    }

    const std::set<ReferenceName>
    Lifetime::references
    {
        "LMPR:2022A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    Lifetime::begin_options()
    {
        return Implementation<Lifetime>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    Lifetime::end_options()
    {
        return Implementation<Lifetime>::options.cend();
    }
}
