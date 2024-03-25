/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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


#include <eos/nonleptonic-amplitudes/topological-amplitudes.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/options-impl.hh>

#include <map>

namespace eos
{
    using std::sqrt;

    NonleptonicAmplitudes<PToPP> *
    TopologicalRepresentation<PToPP>::make(const Parameters & p, const Options & o)
    {
        return new TopologicalRepresentation<PToPP>(p, o);
    }


    const std::vector<OptionSpecification>
    TopologicalRepresentation<PToPP>::options
    {
        Model::option_specification(),
        { "q", { "u", "d", "s" }, "d" },
        { "P1", { "pi^0", "pi^+", "K_d", "K_u", "eta_q", "eta'_q", "eta_s", "eta'_s" }, "pi^0" },
        { "P2", { "pi^0", "pi^+", "K_d", "K_u", "eta_q", "eta'_q", "eta_s", "eta'_s" }, "pi^0" },
    };

    complex<double>
    TopologicalRepresentation<PToPP>::tree_amplitude() const
    {
        complex<double> T_tda = 0.0;
        double T = this->T(), C = this->C(),
            A = this->A(), E = this->E(),
            TES = this->TES(), TAS = this->TAS(),
            TS = this->TS(), TPA = this->TPA(),
            TP = this->TP(), TSS = this->TSS();

        for (unsigned i = 0; i < 3; i++)
        {
            for (unsigned j = 0; j < 3; j++)
            {
                for (unsigned k = 0; k < 3; k++)
                {
                    for (unsigned l = 0; l < 3; l++)
                    {
                        T_tda += T * B[i] * P1[i][j] * Hbar[j][l][k] * P2[k][l];
                        // T_tda += C * ...;
                    }
                }
            }
        }

        return T_tda;
    }

    complex<double>
    TopologicalRepresentation<PToPP>::penguin_amplitude() const
    {
        complex<double> P_tda = 0.0;
        double P = this->P(), PT = this->PT(),
            S = this->S(), PC = this->PC(),
            PTA = this->PTA(), PA = this->PA(),
            PTE = this->PTE(), PAS = this->PAS(),
            PSS = this->PSS(), PES = this->PES();

        for (unsigned i = 0; i < 3; i++)
        {
            for (unsigned j = 0; j < 3; j++)
            {
                for (unsigned k = 0; k < 3; k++)
                {
                    for (unsigned l = 0; l < 3; l++)
                    {
                        P_tda += P * B[i] * P1[i][j] * P2[j][k] * H1tilde[k];
                        // P_tda += PT * ...;
                    }
                }
            }
        }

        return P_tda;
    }

    complex<double>
    TopologicalRepresentation<PToPP>::amplitude() const
    {
        return complex<double>(0.0, 1.0) * Gfermi() / sqrt(2.0) * (
            this->tree_amplitude() + this->penguin_amplitude()
        );
    }
}
