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
        { "cp-conjugate", { "true", "false" },  "false" },
        { "B_bar", { "true", "false"},  "false" },
        { "q", { "u", "d", "s" }, "" },
        { "P1", { "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_S", "K_u", "Kbar_u", "eta", "eta_prime" }, "" },
        { "P2", { "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_S", "K_u", "Kbar_u", "eta", "eta_prime" }, "" },
    };

    complex<double>
    TopologicalRepresentation<PToPP>::tree_amplitude(su3f::rank2 p1, su3f::rank2 p2) const
    {
        complex<double> T_tda = 0.0;
        const complex<double> T   = complex<double>(this->re_T(),   this->im_T()),
                              C   = complex<double>(this->re_C(),   this->im_C()),
                              A   = complex<double>(this->re_A(),   this->im_A()),
                              E   = complex<double>(this->re_E(),   this->im_E()),
                              TES = complex<double>(this->re_TES(), this->im_TES()),
                              TAS = complex<double>(this->re_TAS(), this->im_TAS()),
                              TS  = complex<double>(this->re_TS(),  this->im_TS()),
                              TPA = complex<double>(this->re_TPA(), this->im_TPA()),
                              TP  = complex<double>(this->re_TP(),  this->im_TP()),
                              TSS = complex<double>(this->re_TSS(), this->im_TSS());

        for (unsigned i = 0; i < 3; i++)
        {
            for (unsigned j = 0; j < 3; j++)
            {
                for (unsigned k = 0; k < 3; k++)
                {
                    for (unsigned l = 0; l < 3; l++)
                    {
                        T_tda += T * B[i] * p1[i][j] * Hbar[j][l][k] * p2[k][l];
                        T_tda += C * B[i] * p1[i][j] * Hbar[l][j][k] * p2[k][l];
                        T_tda += A * B[i] * Hbar[i][l][j] * p1[j][k] * p2[k][l];
                        T_tda += E * B[i] * Hbar[l][i][j] * p1[j][k] * p2[k][l];
                        T_tda += TES * B[i] * Hbar[i][j][l] * p1[l][j] * p2[k][k];
                        T_tda += TAS * B[i] * Hbar[j][i][l] * p1[l][j] * p2[k][k];
                        T_tda += TS * B[i] * p1[i][j] * Hbar[l][j][l] * p2[k][k];
                        T_tda += TPA * B[i] * Hbar[l][i][l] * p1[j][k] * p2[k][j];
                        T_tda += TP * B[i] * p1[i][j] * p2[j][k] * Hbar[l][k][l];
                        T_tda += TSS * B[i] * Hbar[l][i][l] * p1[j][j] * p2[k][k];
                    }
                }
            }
        }

        return T_tda;
    }

    complex<double>
    TopologicalRepresentation<PToPP>::penguin_amplitude(su3f::rank2 p1, su3f::rank2 p2) const
    {
        complex<double> P_tda = 0.0;
        const complex<double> P   = complex<double>(this->re_P(),   this->im_P()),
                              PT  = complex<double>(this->re_PT(),  this->im_PT()),
                              S   = complex<double>(this->re_S(),   this->im_S()),
                              PC  = complex<double>(this->re_PC(),  this->im_PC()),
                              PTA = complex<double>(this->re_PTA(), this->im_PTA()),
                              PA  = complex<double>(this->re_PA(),  this->im_PA()),
                              PTE = complex<double>(this->re_PTE(), this->im_PTE()),
                              PAS = complex<double>(this->re_PAS(), this->im_PAS()),
                              PSS = complex<double>(this->re_PSS(), this->im_PSS()),
                              PES = complex<double>(this->re_PES(), this->im_PES());

        for (unsigned i = 0; i < 3; i++)
        {
            for (unsigned j = 0; j < 3; j++)
            {
                for (unsigned k = 0; k < 3; k++)
                {
                    P_tda += P * B[i] * p1[i][j] * p2[j][k] * H1tilde[k];
                    P_tda += S * B[i] * p1[i][j] * H1tilde[j] * p2[k][k];
                    P_tda += PA * B[i] * H1tilde[i] * p1[j][k] * p2[k][j];
                    P_tda += PSS * B[i] * H1tilde[i] * p1[j][j] * p2[k][k];

                    for (unsigned l = 0; l < 3; l++)
                    {
                        P_tda += PT * B[i] * p1[i][j] * H3tilde[j][l][k] * p2[k][l];
                        P_tda += PC * B[i] * p1[i][j] * H3tilde[l][j][k] * p2[k][l];
                        P_tda += PTA * B[i] * H3tilde[i][l][j] * p1[j][k] * p2[k][l];
                        P_tda += PTE * B[i] * H3tilde[j][i][k] * p1[k][l] * p2[l][j];
                        P_tda += PAS * B[i] * H3tilde[j][i][l] * p1[l][j] * p2[k][k];
                        P_tda += PES * B[i] * H3tilde[i][j][l] * p1[l][j] * p2[k][k];
                    }
                }
            }
        }

        return P_tda;
    }

    complex<double>
    TopologicalRepresentation<PToPP>::ordered_amplitude() const
    {
        this->update();

        if (opt_B_bar.value())
        {
            su3f::transpose(P1);
            su3f::transpose(P2);
        }

        return complex<double>(0.0, 1.0) * Gfermi() / sqrt(2.0) * (
            this->tree_amplitude(P1, P2) + this->penguin_amplitude(P1, P2)
        );
    }

    complex<double>
    TopologicalRepresentation<PToPP>::inverse_amplitude() const
    {
        this->update();

        if (opt_B_bar.value())
        {
            su3f::transpose(P1);
            su3f::transpose(P2);
        }

        return complex<double>(0.0, 1.0) * Gfermi() / sqrt(2.0) * (
            this->tree_amplitude(P2, P1) + this->penguin_amplitude(P2, P1)
        );
    }
}
