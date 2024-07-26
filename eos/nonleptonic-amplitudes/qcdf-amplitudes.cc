/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Marta Burgos
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


#include <eos/nonleptonic-amplitudes/qcdf-amplitudes.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/options-impl.hh>

#include <map>

namespace eos
{
    using std::sqrt;

    NonleptonicAmplitudes<PToPP> *
    QCDFRepresentation<PToPP>::make(const Parameters & p, const Options & o)
    {
        return new QCDFRepresentation<PToPP>(p, o);
    }


    const std::vector<OptionSpecification>
    QCDFRepresentation<PToPP>::options
    {
        Model::option_specification(),
        { "cp-conjugate", { "true", "false" },  "false" },
        { "q", { "u", "d", "s" } },
        { "P1", { "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_u", "Kbar_u", "eta", "eta_prime" } },
        { "P2", { "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_u", "Kbar_u", "eta", "eta_prime" } },
    };

    complex<double>
    QCDFRepresentation<PToPP>::alpha_amplitude(su3f::rank2 p1, su3f::rank2 p2) const
    {
        complex<double> A_alpha_qcdf = 0.0;

        complex<double>               alpha1     = complex<double>(this->re_alpha1(),    this->im_alpha1()),
                                      alpha2     = complex<double>(this->re_alpha2(),    this->im_alpha2()),
                                      b1      = complex<double>(this->re_b1(),     this->im_b1()),
                                      b2      = complex<double>(this->re_b2(),     this->im_b2()),
                                      bS1     = complex<double>(this->re_bS1(),    this->im_bS1()),
                                      bS2     = complex<double>(this->re_bS2(),    this->im_bS2()),

                                      alpha3_u   = complex<double>(this->re_alpha3_u(),   this->im_alpha3_u()),
                                      alpha3_c   = complex<double>(this->re_alpha3_c(),   this->im_alpha3_c()),
                                      alpha4_u   = complex<double>(this->re_alpha4_u(),   this->im_alpha4_u()),
                                      alpha4_c   = complex<double>(this->re_alpha4_c(),   this->im_alpha4_c()),
                                      b4_u    = complex<double>(this->re_b4_u(),    this->im_b4_u()),
                                      b4_c    = complex<double>(this->re_b4_c(),    this->im_b4_c()),
                                      bS4_u   = complex<double>(this->re_bS4_u(),   this->im_bS4_u()),
                                      bS4_c   = complex<double>(this->re_bS4_c(),   this->im_bS4_c()),

                                      alpha3EW_c = complex<double>(this->re_alpha3EW_c(), this->im_alpha3EW_c()),
                                      alpha4EW_c = complex<double>(this->re_alpha4EW_c(), this->im_alpha4EW_c()),
                                      b3EW_c  = complex<double>(this->re_b3EW_c(),  this->im_b3EW_c()),
                                      bS3EW_c = complex<double>(this->re_bS3EW_c(), this->im_bS3EW_c()),
                                      b4EW_c     = complex<double>(this->re_b4EW_c(),     this->im_b4EW_c()),
                                      bS4EW_c    = complex<double>(this->re_bS4EW_c(),    this->im_bS4EW_c());

                T[0] = alpha1;
                T[1] = alpha2;
                T[2] = b2;
                T[3] = b1;
                T[4] = bS2;
                T[5] = bS1;
 
                P1_u[0] = alpha4_u;
                P1_u[1] = alpha3_u;
                P1_u[3] = b4_u;
                P1_u[5] = bS4_u;

                P1_c[0] = alpha4_c;
                P1_c[1] = alpha3_c;
                P1_c[3] = b4_c;
                P1_c[5] = bS4_c;

                P2_c[0] = alpha4EW_c;
                P2_c[1] = alpha3EW_c;
                P2_c[2] = b3EW_c;
                P2_c[3] = b4EW_c;
                P2_c[4] = bS3EW_c;
                P2_c[5] = bS4EW_c;

                for (unsigned i = 0; i < 3; i++)
                { 
                    for (unsigned j = 0; j < 3; j++)
                    { 
                        for (unsigned a = 0; a < 6; a++)
                        { 
                            C_u[a][i][j] = T[a] * U[i][j] + P1_u[a] * I[i][j];
                            C_c[a][i][j] = 3.0/2.0 * P2_c[a] * U[i][j] + P1_c[a] * I[i][j] ;
                        }
                    }
                }
        
        for (unsigned i = 0; i < 3; i++)
        {
            for (unsigned j = 0; j < 3; j++)
            {
                for (unsigned k = 0; k < 3; k++)
                {
                    for (unsigned l = 0; l < 3; l++)
                    {

                        A_alpha_qcdf += B[i] * p1[i][j] * C_u[0][j][k] * p2[k][l] * Lambda_u[l];
                        A_alpha_qcdf += B[i] * p1[i][j] * Lambda_u[j] * C_u[1][l][k] * p2[k][l];

                        A_alpha_qcdf += B[i] * p1[i][j] * C_c[0][j][k] * p2[k][l] * Lambda_c[l];
                        A_alpha_qcdf += B[i] * p1[i][j] * Lambda_c[j] * C_c[1][l][k] * p2[k][l];

                    }
                }
            }
        }

        return complex<double>(0.0, 1.0) * Gfermi() / sqrt(2.0) * A_alpha_qcdf;
    }

    complex<double>
    QCDFRepresentation<PToPP>::b_amplitude(su3f::rank2 p1, su3f::rank2 p2) const
    {
        complex<double> A_b_qcdf = 0.0;

        complex<double> alpha1     = complex<double>(this->re_alpha1(),    this->im_alpha1()),
                                      alpha2     = complex<double>(this->re_alpha2(),    this->im_alpha2()),
                                      b1      = complex<double>(this->re_b1(),     this->im_b1()),
                                      b2      = complex<double>(this->re_b2(),     this->im_b2()),
                                      bS1     = complex<double>(this->re_bS1(),    this->im_bS1()),
                                      bS2     = complex<double>(this->re_bS2(),    this->im_bS2()),

                                      alpha3_u   = complex<double>(this->re_alpha3_u(),   this->im_alpha3_u()),
                                      alpha3_c   = complex<double>(this->re_alpha3_c(),   this->im_alpha3_c()),
                                      alpha4_u   = complex<double>(this->re_alpha4_u(),   this->im_alpha4_u()),
                                      alpha4_c   = complex<double>(this->re_alpha4_c(),   this->im_alpha4_c()),
                                      b4_u    = complex<double>(this->re_b4_u(),    this->im_b4_u()),
                                      b4_c    = complex<double>(this->re_b4_c(),    this->im_b4_c()),
                                      bS4_u   = complex<double>(this->re_bS4_u(),   this->im_bS4_u()),
                                      bS4_c   = complex<double>(this->re_bS4_c(),   this->im_bS4_c()),

                                      alpha3EW_c = complex<double>(this->re_alpha3EW_c(), this->im_alpha3EW_c()),
                                      alpha4EW_c = complex<double>(this->re_alpha4EW_c(), this->im_alpha4EW_c()),
                                      b3EW_c  = complex<double>(this->re_b3EW_c(),  this->im_b3EW_c()),
                                      bS3EW_c = complex<double>(this->re_bS3EW_c(), this->im_bS3EW_c()),
                                      b4EW_c     = complex<double>(this->re_b4EW_c(),     this->im_b4EW_c()),
                                      bS4EW_c    = complex<double>(this->re_bS4EW_c(),    this->im_bS4EW_c());

                T[0] = alpha1;
                T[1] = alpha2;
                T[2] = b2;
                T[3] = b1;
                T[4] = bS2;
                T[5] = bS1;
 
                P1_u[0] = alpha4_u;
                P1_u[1] = alpha3_u;
                P1_u[3] = b4_u;
                P1_u[5] = bS4_u;

                P1_c[0] = alpha4_c;
                P1_c[1] = alpha3_c;
                P1_c[3] = b4_c;
                P1_c[5] = bS4_c;

                P2_c[0] = alpha4EW_c;
                P2_c[1] = alpha3EW_c;
                P2_c[2] = b3EW_c;
                P2_c[3] = b4EW_c;
                P2_c[4] = bS3EW_c;
                P2_c[5] = bS4EW_c;

                for (unsigned i = 0; i < 3; i++)
                { 
                    for (unsigned j = 0; j < 3; j++)
                    { 
                        for (unsigned a = 0; a < 6; a++)
                        { 
                            C_u[a][i][j] = T[a] * U[i][j] + P1_u[a] * I[i][j];
                            C_c[a][i][j] = 3.0/2.0 * P2_c[a] * U[i][j] + P1_c[a] * I[i][j] ;
                        }
                    }
                }
        
        for (unsigned i = 0; i < 3; i++)
        {
            for (unsigned j = 0; j < 3; j++)
            {
                for (unsigned k = 0; k < 3; k++)
                {
                    for (unsigned l = 0; l < 3; l++)
                    {

                        A_b_qcdf += B[i] * C_u[2][i][k] * p1[k][l] * p2[l][j] * Lambda_u[j];
                        A_b_qcdf += B[i] * Lambda_u[i] * C_u[3][l][k] * p1[k][j] * p2[j][l];
                        A_b_qcdf += B[i] * C_u[4][i][k] * p1[k][j] * Lambda_u[j] * p2[l][l];
                        A_b_qcdf += B[i] * Lambda_u[i] * C_u[5][j][k] * p1[k][j] * p2[l][l];

                        A_b_qcdf += B[i] * C_c[2][i][k] * p1[k][l] * p2[l][j] * Lambda_c[j];
                        A_b_qcdf += B[i] * Lambda_c[i] * C_c[3][l][k] * p1[k][j] * p2[j][l];
                        A_b_qcdf += B[i] * C_c[4][i][k] * p1[k][j] * Lambda_c[j] * p2[l][l];
                        A_b_qcdf += B[i] * Lambda_c[i] * C_c[5][j][k] * p1[k][j] * p2[l][l];

                    }
                }
            }
        }

        return complex<double>(0.0, 1.0) * Gfermi() / sqrt(2.0) * A_b_qcdf;
    }

    complex<double>
    QCDFRepresentation<PToPP>::ordered_amplitude() const
    {
        this->update();

        return power_of<2>(mB()) * FP1() * fP2() * this->alpha_amplitude(P1, P2) + fB() * fP1() * fP2() * this->b_amplitude(P1, P2);
        
    }

    complex<double>
    QCDFRepresentation<PToPP>::inverse_amplitude() const
    {
        this->update();

        return power_of<2>(mB()) * FP2() * fP1() * this->alpha_amplitude(P2, P1) + fB() * fP1() * fP2() * this->b_amplitude(P2, P1);
        
    }
}
